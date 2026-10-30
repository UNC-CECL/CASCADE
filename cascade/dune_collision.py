"""Lower the dunes that a storm reaches but does not overtop

This module provides functions for adding collision-regime dune erosion to Barrier3D. Barrier3D
lowers a dune only where a storm overtops it (Rhigh above the crest). Sallenger (2000) puts a storm
in the collision regime when the water reaches the dune toe but not the crest
(D_low < R_high < D_high), where real dunes lose sand to wave impact on the dune face. The
DuneCollisionManager adds that loss storm by storm, inside Barrier3D's own storm loop, using:

1. the wave-impact rule of Larson et al. (2004), per storm: dV = 4 Cs (R_high - z_toe)^2 t / T, the
   storm's peak R_high held over its duration t, T its peak period, and
2. optionally, the beach in front of the dune: dV scaled by max(0, 1 - W / W0), W the beach width, as
   the dunes behind wide beaches lost no sand to collision in the Hatteras lidar; W fixed, or
   (beach_follows_shoreline) narrowing as the shoreline retreats toward a dune held in place, and
3. the crest lowering over the dune's cross-shore width: dH = dV / W_dune, both dune rows lowered in
   proportion, no lower than the toe.

Calibration on the 1996-2009 lidar at Hatteras Island, NC:
scripts/real-world-cases/Hatteras/dune-storm-response/README.md.

References
----------

.. [1] Sallenger, A. H. (2000). Storm impact scale for barrier islands. Journal of Coastal
       Research, 16(3), 890-895.
.. [2] Larson, M., Erikson, L., & Hanson, H. (2004). An analytical model to predict dune erosion
       due to wave impact. Coastal Engineering, 51(8-9), 675-696.

Notes
-----
Barrier3D's units: elevations in dam relative to MHW, dune heights in dam above the berm. The toe and
the beach width are taken per alongshore dune cell, in dam, and held for the run.

Barrier3D's source is not changed, and Barrier3D has no step between storms that CASCADE can call.
So the manager stands in for two of Barrier3D's methods on the instance it is given, each running
Barrier3D's own method as before: DuneGrowth, called once a year just before the storms, after
which the manager keeps the pre-storm crest Barrier3D uses to find the overwashed cells; and
DiffuseDunes, called once after each storm's overwash lowering, before which the manager lowers the
cells that storm left in the collision regime (update). A cell is in collision when Barrier3D did
not overwash it (pre-storm crest at or above R_high), the water reaches its toe, and its current
crest still stands above R_high. Later storms in the same year see the lowered dunes.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
"""

import numpy as np

__all__ = [
    "DuneCollisionManager",
    "beach_factor",
    "collision_cells",
    "larson_volume",
    "lower_dunes",
]


def beach_factor(beach_width, beach_width_limit):
    """
    The share of the collision volume a dune loses behind its beach: 1 on no beach, falling
    linearly to 0 at ``beach_width_limit`` and beyond.

    Parameters
    ----------
    beach_width: array of float
        Beach width per alongshore dune cell [dam]
    beach_width_limit: float
        The beach width W0 at and beyond which a dune loses nothing to collision [dam]

    Returns
    -------
    array of float
        max(0, 1 - W / W0) per dune cell [-]
    """

    return np.clip(1 - beach_width / beach_width_limit, 0, None)


def collision_cells(crest_prestorm, crest, toe_elevation, r_high, beach):
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
    beach: array of float
        Beach factor per cell (beach_factor) [-]; cells at 0 are left out

    Returns
    -------
    array of bool
        True where the dune is in collision
    """

    return (
        (crest_prestorm >= r_high)  # not overwashed by Barrier3D
        & (r_high > toe_elevation)
        & (crest > r_high)
        & (beach > 0)
    )


def larson_volume(r_high, toe_elevation, period, duration, cs, beach):
    """
    Following Larson et al. (2004), the dune volume eroded by wave impact over one storm,
    with the storm's peak runup held over its duration and scaled by the beach in front of
    the dune.

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
    beach: array of float
        Beach factor per cell (beach_factor) [-]

    Returns
    -------
    array of float
        Eroded volume per metre of shore [m^3/m]
    """

    excess_m = (r_high - toe_elevation) * 10  # dam to m
    return 4 * cs * excess_m**2 * duration * 3600.0 / period * beach


def lower_dunes(dune_rows, volume, toe_elevation, berm_elevation, dune_width, dune_restart):
    """
    Lower each dune cell's crest by the eroded volume over the dune width, no lower than its
    toe, both dune rows in proportion to their height.

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
    scale = np.where(top > 0, new_top / np.where(top > 0, top, 1.0), 1.0)

    return np.minimum(np.maximum(dune_rows * scale[:, None], dune_restart), dune_rows)


class DuneCollisionManager:
    """Lower the dunes that each storm reaches but does not overtop.

    Examples
    --------
    # >>> from cascade.dune_collision import DuneCollisionManager
    # >>> collision = DuneCollisionManager(barrier3d, toe_elevation_dam, cs=1.24e-4)
    # >>> barrier3d.update()  # collision erosion runs within each storm

    The manager installs itself on ``barrier3d``; each later ``barrier3d.update()`` includes it.
    """

    def __init__(
        self,
        barrier3d,
        toe_elevation,
        cs,
        dune_width=2.0,
        beach_width=None,
        beach_width_limit=4.5,
        beach_follows_shoreline=False,
    ):
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
        beach_width: array of float, optional
            Beach width per alongshore dune cell [dam]; None leaves the beach out
        beach_width_limit: float, optional
            The beach width W0 at and beyond which a dune loses nothing to collision [dam]
        beach_follows_shoreline: bool, optional
            If True, the dune is held in place and the beach takes the shoreline's change: each
            year's beach width is beach_width less the shoreline's landward move since the start
            (x_s_TS), no less than 0, so it narrows as the shoreline retreats and widens with
            nourishment
        """

        n = barrier3d.BarrierLength
        toe = np.asarray(toe_elevation, dtype=float)
        if toe.shape != (n,):
            raise ValueError(
                f"toe_elevation needs one value per dune cell ({n}), got shape {toe.shape}"
            )
        if beach_width is None:
            beach = np.ones(n)
        else:
            beach_width = np.asarray(beach_width, dtype=float)
            if beach_width.shape != (n,):
                raise ValueError(
                    f"beach_width needs one value per dune cell ({n}), "
                    f"got shape {beach_width.shape}"
                )
            beach = beach_factor(beach_width, beach_width_limit)

        self._barrier3d = barrier3d
        self._toe_elevation = toe
        self._beach_width = beach_width
        self._beach_width_limit = beach_width_limit
        self._beach = beach
        self._beach_follows_shoreline = beach_follows_shoreline and beach_width is not None
        self._shoreline_start = barrier3d.x_s_TS[-1]  # dam, landward positive
        self._beach_width_TS = []  # per model year, mean beach width [dam]
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
        if self._beach_follows_shoreline:
            retreat = b3d.x_s_TS[-1] - self._shoreline_start
            width = np.maximum(self._beach_width - retreat, 0)
            self._beach = beach_factor(width, self._beach_width_limit)
            self._beach_width_TS.append(float(width.mean()))
        self._crest_lowering_TS.append(0.0)
        self._volume_TS.append(0.0)

        return DuneDomain, Dmax, Qdg

    def diffuse_dunes(self, DuneDomain, t):
        """Collision lowering for the storm just run (update); then Barrier3D's DiffuseDunes."""

        b3d = self._barrier3d
        self.update(b3d, DuneDomain, t)

        return type(b3d).DiffuseDunes(b3d, DuneDomain, t)

    def update(self, barrier3d, DuneDomain, t):
        """Lower the dunes the year's next storm left in the collision regime.

        Parameters
        ----------
        barrier3d: Barrier3d
            The domain, mid storm loop (its storm series is read)
        DuneDomain: array of float
            Barrier3D's dune domain, (years, cells, dune rows) [dam above the berm]; changed in place
        t: int
            Barrier3D's time index
        """

        storms = np.flatnonzero(barrier3d._StormSeries[:, 0] == t)
        if self._storm >= len(storms):
            return
        r_high, period, duration = barrier3d._StormSeries[storms[self._storm], [1, 3, 4]]
        self._storm += 1

        dune_rows = DuneDomain[t]  # (cells, dune rows)
        crest = dune_rows.max(axis=1) + barrier3d._BermEl  # dam MHW
        hit = collision_cells(
            self._crest_prestorm, crest, self._toe_elevation, r_high, self._beach
        )
        if not hit.any():
            return

        volume = larson_volume(
            r_high, self._toe_elevation[hit], period, duration, self._cs, self._beach[hit]
        )
        new_rows = lower_dunes(
            dune_rows[hit],
            volume,
            self._toe_elevation[hit],
            barrier3d._BermEl,
            self._dune_width,
            barrier3d._DuneRestart,
        )

        n = len(dune_rows)
        self._crest_lowering_TS[-1] += float(
            (dune_rows[hit].max(axis=1) - new_rows.max(axis=1)).sum() / n
        )
        self._volume_TS[-1] += float((dune_rows[hit] - new_rows).sum() / n)
        dune_rows[hit] = new_rows

    @property
    def toe_elevation(self):
        return self._toe_elevation

    @property
    def beach_width(self):
        return self._beach_width

    @property
    def beach_width_limit(self):
        return self._beach_width_limit

    @property
    def beach_follows_shoreline(self):
        return self._beach_follows_shoreline

    @property
    def beach_width_TS(self):
        return self._beach_width_TS

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
