"""Lower the dunes that a storm reaches but does not overtop

This module provides functions for adding collision-regime dune erosion to Barrier3D. Barrier3D
lowers a dune only where a storm overtops it (Rhigh above the crest). Sallenger (2000) puts a storm
in the collision regime when the water reaches the dune toe but not the crest
(D_low < R_high < D_high), where real dunes lose sand to wave impact on the dune face. The
DuneCollisionManager adds that loss storm by storm, inside Barrier3D's own storm loop, using:

1. the wave-impact rule of Larson et al. (2004), per storm: dV = 4 Cs (R_high - z_toe)^2 t / T, the
   storm's peak R_high held over its duration t, T its peak period, and
2. the crest lowering over the dune's cross-shore width: dH = dV / W_dune, both dune rows lowered in
   proportion, no lower than the toe.

Calibration on the 1996-2009 lidar at Hatteras Island, NC:
scripts/real-world-cases/Hatteras/dune-collision/README.md.

References
----------

.. [1] Sallenger, A. H. (2000). Storm impact scale for barrier islands. Journal of Coastal
       Research, 16(3), 890-895.
.. [2] Larson, M., Erikson, L., & Hanson, H. (2004). An analytical model to predict dune erosion
       due to wave impact. Coastal Engineering, 51(8-9), 675-696.

Notes
-----
Units are Barrier3D's: elevations in dam above MHW, dune heights in dam above the berm. The toe is
per dune cell and fixed for the run.

Barrier3D has no hook between storms, so the manager replaces two methods on the Barrier3d instance,
each still calling the original: DuneGrowth (once a year, before the storms), after which it saves
the pre-storm crest; and DiffuseDunes (after each storm's overwash), before which it lowers the
collision cells (update). Later storms in the same year see the lowered dunes.

Outside the conditions tested at Hatteras (1996-2025 storms), the manager raises
DuneCollisionOutsideTestedRange: toes below MHW, or Cs, storm duration, peak period or water above
the toe beyond the limits set as DuneCollisionManager's class constants.
Tests: scripts/real-world-cases/Hatteras/dune-collision/3-module-tests/.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
"""

import numpy as np

__all__ = [
    "DuneCollisionManager",
    "DuneCollisionOutsideTestedRange",
    "collision_cells",
    "larson_volume",
    "lower_dunes",
]


class DuneCollisionOutsideTestedRange(ValueError):
    """The run has left the conditions the dune collision module was tested in; its results would not hold."""


def collision_cells(crest_prestorm, crest, toe_elevation, r_high):
    """
    The dune cells a storm leaves in the collision regime (Sallenger, 2000): not overwashed by
    Barrier3D, the water above the toe, and the crest still above the water.

    Parameters
    ----------
    crest_prestorm: array of float
        Crest elevation before the year's storms, as Barrier3D finds the overwashed cells
        [dam MHW]
    crest: array of float
        Crest elevation now, after this storm's overwash lowering [dam MHW]
    toe_elevation: array of float
        Dune toe elevation [dam MHW]
    r_high: float
        The storm's peak runup elevation [dam MHW]

    Returns
    -------
    array of bool
        True where the dune is in collision
    """

    return (
        (crest_prestorm >= r_high)  # not overwashed by Barrier3D
        & (r_high > toe_elevation)
        & (crest > r_high)
    )


def larson_volume(r_high, toe_elevation, period, duration, cs):
    """
    Following Larson et al. (2004), the dune volume eroded by wave impact over one storm,
    with the storm's peak runup held over its duration.

    Parameters
    ----------
    r_high: float
        The storm's peak runup elevation [dam MHW]
    toe_elevation: array of float
        Dune toe elevation [dam MHW]
    period: float
        Peak wave period [s]
    duration: float
        Storm duration [hr]
    cs: float
        Larson et al. (2004) impact coefficient [-]

    Returns
    -------
    array of float
        Eroded volume per metre of shore [m^3/m]
    """

    excess_m = (r_high - toe_elevation) * 10  # dam to m
    return 4 * cs * excess_m**2 * duration * 3600.0 / period


def lower_dunes(dune_rows, volume, toe_elevation, berm_elevation, dune_width, dune_restart):
    """
    Lower each dune cell's crest by the eroded volume over the dune width, no lower than its
    toe, both dune rows in proportion to their height. Cells with no volume are returned
    exactly as given.

    Parameters
    ----------
    dune_rows: array of float
        Dune heights above the berm, (cells, dune rows) [dam]
    volume: array of float
        Eroded volume per metre of shore [m^3/m]
    toe_elevation: array of float
        Dune toe elevation [dam MHW]
    berm_elevation: float
        Berm elevation [dam MHW]
    dune_width: float
        Cross-shore width the eroded volume is spread over [dam]
    dune_restart: float
        Barrier3D's smallest dune height [dam]

    Returns
    -------
    array of float
        New dune heights above the berm, (cells, dune rows) [dam]
    """

    top = dune_rows.max(axis=1)
    crest = top + berm_elevation
    new_crest = np.maximum(crest - volume / (dune_width * 10) / 10, toe_elevation)
    new_top = np.maximum(new_crest - berm_elevation, dune_restart)
    scale = new_top / np.where(top > 0, top, 1.0)

    new_rows = np.minimum(np.maximum(dune_rows * scale[:, None], dune_restart), dune_rows)

    # no volume, no change: the crest-to-height round trip above is not exact in floating point
    return np.where((np.asarray(volume) > 0)[:, None], new_rows, dune_rows)


class DuneCollisionManager:
    """Lower the dunes that each storm reaches but does not overtop.

    Examples
    --------
    # >>> from cascade.dune_collision import DuneCollisionManager
    # >>> collision = DuneCollisionManager(barrier3d, toe_elevation_dam, cs=1.9e-5)
    # >>> barrier3d.update()  # collision erosion runs within each storm

    The manager installs itself on ``barrier3d``; each later ``barrier3d.update()`` includes it.
    """

    # The tested range; outside it the manager raises DuneCollisionOutsideTestedRange
    CS_MAX = 2.5e-4  # the highest Cs run at Hatteras [-]
    MAX_STORM_DURATION_HR = 200.0  # the storm series run reach 193 h [hr]
    STORM_PERIOD_RANGE_S = (5.0, 20.0)  # the storm series run span 7-19 s [s]
    MAX_WATER_ABOVE_TOE_DAM = 0.5  # 5 m; the storms run put water up to 4.6 m above a toe [dam]

    def __init__(self, barrier3d, toe_elevation, cs, dune_width=2.0):
        """The DuneCollisionManager module

        Parameters
        ----------
        barrier3d: Barrier3d, required
            The domain to add collision erosion to
        toe_elevation: array of float, required
            Dune toe elevation per alongshore dune cell [dam MHW]
        cs: float, required
            Larson et al. (2004) impact coefficient [-]
        dune_width: float, optional
            Cross-shore width the eroded volume is spread over [dam]; Barrier3D's two dune rows
        """

        n = barrier3d.BarrierLength
        toe = np.asarray(toe_elevation, dtype=float)
        if toe.shape != (n,):
            raise ValueError(
                f"toe_elevation needs one value per dune cell ({n}), got shape {toe.shape}"
            )
        if not 0 <= cs <= self.CS_MAX:
            raise DuneCollisionOutsideTestedRange(
                f"cs {cs:g} is outside the tested range [0, {self.CS_MAX:g}]; the module was calibrated "
                "and run at Hatteras with Cs from 1.0e-5 to 2.5e-4 "
                "(scripts/real-world-cases/Hatteras/dune-collision/3-module-tests/)."
            )
        if np.nanmin(toe) < 0:
            raise DuneCollisionOutsideTestedRange(
                f"a dune toe lies {-np.nanmin(toe) * 10:.2f} m below MHW; the module was tested only "
                "with toes at or above MHW (the toe sits on the dry beach)."
            )

        self._barrier3d = barrier3d
        self._toe_elevation = toe
        self._cs = cs
        self._dune_width = dune_width
        self._crest_prestorm = None  # dam MHW, set each year after dune growth
        self._storm = 0  # the year's storms handled so far
        # Per model year: mean crest lowering [dam] and eroded volume [dam^3/dam of shore]
        self._crest_lowering_TS = []
        self._volume_TS = []

        # Stand in for Barrier3D's two methods on this instance (module Notes)
        barrier3d.DuneGrowth = self.dune_growth
        barrier3d.DiffuseDunes = self.diffuse_dunes

    def dune_growth(self, DuneDomain, t):
        """Barrier3D's DuneGrowth; then keep the pre-storm crest and restart the storm count."""

        b3d = self._barrier3d
        DuneDomain, Dmax, Qdg = type(b3d).DuneGrowth(b3d, DuneDomain, t)

        # as Barrier3D's DuneDomainCrest, from which it finds the overwashed cells
        crest = DuneDomain[t].max(axis=1)
        self._crest_prestorm = np.maximum(crest, b3d._DuneRestart) + b3d._BermEl
        self._storm = 0
        self._crest_lowering_TS.append(0.0)
        self._volume_TS.append(0.0)

        return DuneDomain, Dmax, Qdg

    def diffuse_dunes(self, DuneDomain, t):
        """Collision lowering for the storm just run (update); then Barrier3D's DiffuseDunes."""

        self.update(DuneDomain, t)

        return type(self._barrier3d).DiffuseDunes(self._barrier3d, DuneDomain, t)

    def update(self, DuneDomain, t):
        """Lower the dunes the year's next storm left in the collision regime.

        Parameters
        ----------
        DuneDomain: array of float
            Barrier3D's dune domain, (years, cells, dune rows) [dam above the berm]; changed in place
        t: int
            Barrier3D's time index
        """

        b3d = self._barrier3d
        storms = np.flatnonzero(b3d._StormSeries[:, 0] == t)
        if self._storm >= len(storms):
            return
        r_high, period, duration = b3d._StormSeries[storms[self._storm], [1, 3, 4]]
        self._storm += 1

        dune_rows = DuneDomain[t]  # (cells, dune rows)
        crest = dune_rows.max(axis=1) + b3d._BermEl  # dam MHW
        hit = collision_cells(self._crest_prestorm, crest, self._toe_elevation, r_high)
        if not hit.any():
            return
        self._check_storm(r_high, self._toe_elevation[hit], period, duration)

        volume = larson_volume(r_high, self._toe_elevation[hit], period, duration, self._cs)
        new_rows = lower_dunes(
            dune_rows[hit],
            volume,
            self._toe_elevation[hit],
            b3d._BermEl,
            self._dune_width,
            b3d._DuneRestart,
        )

        n = len(dune_rows)
        self._crest_lowering_TS[-1] += float(
            (dune_rows[hit].max(axis=1) - new_rows.max(axis=1)).sum() / n
        )
        self._volume_TS[-1] += float((dune_rows[hit] - new_rows).sum() / n)
        dune_rows[hit] = new_rows

    def _check_storm(self, r_high, toe, period, duration):
        """Raise DuneCollisionOutsideTestedRange if a storm that lowers dunes is outside the tested range."""
        if duration > self.MAX_STORM_DURATION_HR:
            raise DuneCollisionOutsideTestedRange(
                f"a storm of {duration:g} h lowers dunes; the module was tested on storms up to "
                f"{self.MAX_STORM_DURATION_HR:g} h, and it holds the peak water level over the whole "
                "storm, so longer storms are increasingly overstated."
            )
        lo, hi = self.STORM_PERIOD_RANGE_S
        if not lo <= period <= hi:
            raise DuneCollisionOutsideTestedRange(
                f"a storm with a {period:g} s peak period lowers dunes; the module was tested on "
                f"periods of {lo:g}-{hi:g} s."
            )
        above = float(np.max(r_high - toe))
        if above > self.MAX_WATER_ABOVE_TOE_DAM:
            raise DuneCollisionOutsideTestedRange(
                f"a storm puts the water {above * 10:.2f} m above a dune toe; the module was tested "
                f"up to {self.MAX_WATER_ABOVE_TOE_DAM * 10:g} m (Larson's volume grows with its square)."
            )

    @property
    def toe_elevation(self):
        return self._toe_elevation

    @property
    def cs(self):
        return self._cs

    @property
    def dune_width(self):
        return self._dune_width

    @property
    def crest_lowering_TS(self):
        return self._crest_lowering_TS

    @property
    def volume_TS(self):
        return self._volume_TS
