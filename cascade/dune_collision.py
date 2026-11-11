"""Lower the dunes that a storm reaches but does not overtop

This module provides functions for adding collision-regime dune erosion to Barrier3D. Barrier3D
lowers a dune only where a storm overtops it (Rhigh above the crest). Sallenger (2000) puts a storm
in the collision regime when the water reaches the dune toe but not the crest
(D_low < R_high < D_high), where real dunes lose sand to wave impact on the dune face. The
DuneCollisionManager adds that loss storm by storm, inside Barrier3D's own storm loop, using:

1. the wave-impact rule of Larson et al. (2004), per storm: dV = 4 Cs (R_high - z_toe)^2 t / T, the
   storm's peak R_high held over its duration t, T its peak period, and
2. the crest lowering over the dune's cross-shore width: dH = dV / W_dune, both dune rows lowered in
   proportion, no lower than the toe (the default; ``rows`` can instead take it all from the front
   row, lower both rows by the same height, or, "lateral", let the face retreat through the front
   row: it keeps its height until the sand it holds above the toe is gone, then drops to the toe).

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


ROW_MODES = ("proportional", "front", "equal")
# The manager also takes "lateral" (erode_laterally): the front cell keeps its height while the face retreats
# through it, which needs the eroded volume kept from storm to storm, so it is not a lower_dunes mode
MANAGER_ROW_MODES = ROW_MODES + ("lateral",)
M3_PER_M_PER_DAM2 = 100.0  # a 1 dam x 1 dam dune cross-section is 100 m^3 per metre of shore


def erode_laterally(dune_rows, deficit, volume, toe_elevation, berm_elevation, dune_restart, row_width):
    """
    Erode each dune cell's rows from the front, as a dune face retreating landward: a row keeps
    its height while the eroded volume builds up against the sand it holds above the toe, and
    drops to the toe (no lower than DuneRestart) once that sand is gone. What is left of the
    storm's volume goes on to the next row. The volume each row has lost so far (the deficit)
    is carried from storm to storm.

    Parameters
    ----------
    dune_rows: array of float
        Dune heights above the berm, (cells, dune rows), front (seaward) row first [dam]
    deficit: array of float
        Volume already eroded from each row, (cells, dune rows) [m^3/m]
    volume: array of float
        This storm's eroded volume per metre of shore [m^3/m]
    toe_elevation: array of float
        Dune toe elevation [dam MHW]
    berm_elevation: float
        Berm elevation [dam MHW]
    dune_restart: float
        Barrier3D's smallest dune height [dam]
    row_width: float
        Cross-shore width of one dune row [dam]

    Returns
    -------
    tuple of array
        New dune heights above the berm and the new deficit, both (cells, dune rows)
    """

    rows = np.array(dune_rows, dtype=float)
    lost = np.array(deficit, dtype=float)
    left = np.array(volume, dtype=float)
    toe = np.asarray(toe_elevation, dtype=float)
    floor = np.maximum(toe - berm_elevation, dune_restart)  # dam above the berm
    eroding = left > 0  # only a storm that erodes a cell removes a row (growth or diffusion can lower a row below its deficit)
    for q in range(rows.shape[1]):
        holds = np.maximum(rows[:, q] + berm_elevation - toe, 0.0) * row_width * M3_PER_M_PER_DAM2
        take = np.minimum(left, np.maximum(holds - lost[:, q], 0.0))
        lost[:, q] += take
        left -= take
        gone = eroding & (holds > 0) & (lost[:, q] >= holds * (1 - 1e-12))
        rows[gone, q] = np.minimum(floor[gone], rows[gone, q])
        lost[gone, q] = 0.0
    return rows, lost


def refill_deficit(grown_rows, previous_rows, deficit, row_width):
    """
    Dune growth refills the eroded face first: each row's growth this year pays back its
    deficit, and only growth beyond that raises the row.

    Parameters
    ----------
    grown_rows: array of float
        Dune heights after Barrier3D's growth, (cells, dune rows) [dam]
    previous_rows: array of float
        The same before it [dam]
    deficit: array of float
        Volume eroded from each row and not yet refilled, (cells, dune rows) [m^3/m]
    row_width: float
        Cross-shore width of one dune row [dam]

    Returns
    -------
    tuple of array
        The rows after the refill and the deficit left, both (cells, dune rows)
    """
    grown = np.asarray(grown_rows, dtype=float)
    growth = np.maximum(grown - np.asarray(previous_rows, dtype=float), 0.0) * row_width * M3_PER_M_PER_DAM2
    repay = np.minimum(np.asarray(deficit, dtype=float), growth)
    return grown - repay / (row_width * M3_PER_M_PER_DAM2), np.asarray(deficit, dtype=float) - repay


def lower_dunes(dune_rows, volume, toe_elevation, berm_elevation, dune_width, dune_restart,
                rows="proportional"):
    """
    Lower each dune cell by the eroded volume, no row lower than its toe. Cells with no volume
    are returned exactly as given. How the volume is taken from the dune rows (``rows``):

    - "proportional": the crest drops by the volume over the dune width, and both rows are
      lowered in proportion to their height (the module's behaviour before the option existed)
    - "front": all of it from the front (seaward, first) row, over that one row's width; the
      rows behind are unchanged
    - "equal": every row drops by the same height, the volume over the dune width

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
    rows: str, optional
        How the volume is taken from the dune rows: "proportional", "front" or "equal"

    Returns
    -------
    array of float
        New dune heights above the berm, (cells, dune rows) [dam]
    """

    if rows not in ROW_MODES:
        raise ValueError(f"rows must be one of {ROW_MODES}, got {rows!r}")
    volume = np.asarray(volume, dtype=float)
    floor = np.maximum(np.asarray(toe_elevation, dtype=float) - berm_elevation, dune_restart)  # dam above the berm

    if rows == "proportional":
        top = dune_rows.max(axis=1)
        crest = top + berm_elevation
        new_crest = np.maximum(crest - volume / (dune_width * 10) / 10, toe_elevation)
        new_top = np.maximum(new_crest - berm_elevation, dune_restart)
        scale = new_top / np.where(top > 0, top, 1.0)
        new_rows = np.minimum(np.maximum(dune_rows * scale[:, None], dune_restart), dune_rows)
    elif rows == "front":
        row_width = dune_width / dune_rows.shape[1]  # dam
        new_rows = dune_rows.copy()
        lowered = np.maximum(dune_rows[:, 0] - volume / (row_width * 10) / 10, floor)
        new_rows[:, 0] = np.minimum(lowered, dune_rows[:, 0])
    else:  # "equal"
        drop = volume / (dune_width * 10) / 10  # dam
        new_rows = np.minimum(np.maximum(dune_rows - drop[:, None], floor[:, None]), dune_rows)

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
    # A row counts as collapsed (dune_row_loss_TS) only if it stood at least this far above its toe before the
    # storm: below it, a cell already flattened and topped up by Barrier3D's alongshore diffusion would be counted
    # again at every storm [dam]
    COLLAPSE_MIN_HEIGHT_DAM = 0.05

    def __init__(self, barrier3d, toe_elevation, cs, dune_width=2.0, rows="proportional"):
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
        rows: str, optional
            Which dune rows lose the volume (lower_dunes): "proportional" (both, in proportion
            to their height), "front" (the seaward row only) or "equal" (both by the same height);
            or "lateral" (erode_laterally): the face retreats through the front row, which keeps
            its height until its sand above the toe is gone, then through the next. The volume
            lost is kept per row, refilled first by dune growth (refill_deficit) and moved with
            the rows when Barrier3D migrates the dune line
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

        if rows not in MANAGER_ROW_MODES:
            raise ValueError(f"rows must be one of {MANAGER_ROW_MODES}, got {rows!r}")

        self._barrier3d = barrier3d
        self._rows = rows
        # lateral: the volume lost by each row and not yet refilled [m^3/m], (cells, dune rows)
        self._deficit = np.zeros((n, barrier3d._DuneWidth)) if rows == "lateral" else None
        self._toe_elevation = toe
        self._cs = cs
        self._dune_width = dune_width
        self._crest_prestorm = None  # dam MHW, set each year after dune growth
        self._storm = 0  # the year's storms handled so far
        # Per model year: mean crest lowering [dam] and eroded volume [dam^3/dam of shore]
        self._crest_lowering_TS = []
        self._volume_TS = []
        # Per model year and dune cell, as Barrier3D's DuneDomain (index: Barrier3D's time index): the volume the
        # collision rule eroded [dam^3/dam], the storms that reached the cell in the collision regime, the rows
        # eroded through to the toe (dune cell collapse), and, lateral, each row's deficit [dam^3/dam]
        nt, nrows = barrier3d._TMAX, barrier3d._DuneWidth
        self._collision_volume_TS = np.zeros((nt, n))
        self._collision_storms_TS = np.zeros((nt, n), dtype=int)
        self._dune_row_loss_TS = np.zeros((nt, n, nrows), dtype=int)
        self._collision_deficit_TS = np.zeros((nt, n, nrows)) if rows == "lateral" else None

        # Stand in for Barrier3D's two methods on this instance (module Notes)
        barrier3d.DuneGrowth = self.dune_growth
        barrier3d.DiffuseDunes = self.diffuse_dunes
        if rows == "lateral":  # the deficit moves with the rows
            barrier3d.migrate_dunes = self.migrate_dunes

    def _row_width(self):
        return self._dune_width / self._barrier3d._DuneWidth  # dam

    def dune_growth(self, DuneDomain, t):
        """Barrier3D's DuneGrowth (lateral: refilling the eroded face first); then keep the
        pre-storm crest and restart the storm count."""

        b3d = self._barrier3d
        lateral = getattr(self, "_rows", "proportional") == "lateral"
        before = DuneDomain[t - 1].copy() if lateral else None
        DuneDomain, Dmax, Qdg = type(b3d).DuneGrowth(b3d, DuneDomain, t)
        if lateral and self._deficit.any():
            DuneDomain[t], self._deficit = refill_deficit(DuneDomain[t], before, self._deficit, self._row_width())
        if lateral and getattr(self, "_collision_deficit_TS", None) is not None and t < len(self._collision_deficit_TS):
            self._collision_deficit_TS[t] = self._deficit / M3_PER_M_PER_DAM2  # after the refill; storms update it

        # as Barrier3D's DuneDomainCrest, from which it finds the overwashed cells
        crest = DuneDomain[t].max(axis=1)
        self._crest_prestorm = np.maximum(crest, b3d._DuneRestart) + b3d._BermEl
        self._storm = 0
        self._crest_lowering_TS.append(0.0)
        self._volume_TS.append(0.0)

        return DuneDomain, Dmax, Qdg

    def migrate_dunes(self, shoreline_change_aggregate, cellular_shoreline_change, time_index):
        """Barrier3D's migrate_dunes, with each row's deficit moved as the rows are (lateral):
        a new front row on progradation starts with none, and the front row's goes with it
        when the shoreline erodes."""

        b3d = self._barrier3d
        before = b3d._DuneDomain[time_index].copy()
        out = type(b3d).migrate_dunes(b3d, shoreline_change_aggregate, cellular_shoreline_change, time_index)
        if np.array_equal(before, b3d._DuneDomain[time_index]):
            return out
        for _ in range(int(cellular_shoreline_change)):
            if shoreline_change_aggregate > 0:
                self._deficit = np.roll(self._deficit, 1, axis=1)
                self._deficit[:, 0] = 0.0
            elif shoreline_change_aggregate < 0:
                self._deficit = np.roll(self._deficit, -1, axis=1)
                self._deficit[:, -1] = 0.0
        return out

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
        n = len(dune_rows)
        mode = getattr(self, "_rows", "proportional")  # managers saved before the option: as before
        if mode == "lateral":
            new_rows, self._deficit[hit] = erode_laterally(
                dune_rows[hit],
                self._deficit[hit],
                volume,
                self._toe_elevation[hit],
                b3d._BermEl,
                b3d._DuneRestart,
                self._row_width(),
            )
            eroded = float(np.sum(volume) / M3_PER_M_PER_DAM2 / n)  # the face's loss, as dam^3/dam
        else:
            new_rows = lower_dunes(
                dune_rows[hit],
                volume,
                self._toe_elevation[hit],
                b3d._BermEl,
                self._dune_width,
                b3d._DuneRestart,
                rows=mode,
            )
            eroded = float((dune_rows[hit] - new_rows).sum() / n)

        self._crest_lowering_TS[-1] += float(
            (dune_rows[hit].max(axis=1) - new_rows.max(axis=1)).sum() / n
        )
        self._volume_TS[-1] += eroded
        self._track(t, hit, volume, dune_rows[hit], new_rows)
        dune_rows[hit] = new_rows

    def _track(self, t, hit, volume, old_rows, new_rows):
        """The per-cell series: this storm's volume, its collision cells and the rows it eroded
        through to the toe (lateral: and the deficit after it)."""
        if getattr(self, "_collision_volume_TS", None) is None or t >= len(self._collision_volume_TS):
            return  # a manager saved before the series existed, or past the run's end
        b3d = self._barrier3d
        self._collision_volume_TS[t, hit] += np.asarray(volume) / M3_PER_M_PER_DAM2
        self._collision_storms_TS[t, hit] += 1
        floor = np.maximum(self._toe_elevation[hit] - b3d._BermEl, b3d._DuneRestart)[:, None]
        collapsed = (old_rows >= floor + self.COLLAPSE_MIN_HEIGHT_DAM) & (new_rows <= floor + 1e-9)
        self._dune_row_loss_TS[t, hit] += collapsed.astype(int)
        if self._collision_deficit_TS is not None:
            self._collision_deficit_TS[t] = self._deficit / M3_PER_M_PER_DAM2

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
    def rows(self):
        return getattr(self, "_rows", "proportional")

    @property
    def deficit(self):
        """lateral: the volume each row has lost and not yet refilled [m^3/m]; None otherwise"""
        return getattr(self, "_deficit", None)

    @property
    def collision_volume_TS(self):
        """Volume the collision rule eroded per model year and dune cell, (years, cells) [dam^3/dam]"""
        return getattr(self, "_collision_volume_TS", None)

    @property
    def collision_storms_TS(self):
        """Storms that reached each cell in the collision regime, per model year, (years, cells)"""
        return getattr(self, "_collision_storms_TS", None)

    @property
    def dune_row_loss_TS(self):
        """Dune rows eroded through to the toe (dune cell collapse) from at least COLLAPSE_MIN_HEIGHT_DAM
        above it, per model year, cell and row, (years, cells, dune rows)"""
        return getattr(self, "_dune_row_loss_TS", None)

    @property
    def collision_deficit_TS(self):
        """lateral: each row's deficit at the end of each model year, (years, cells, dune rows)
        [dam^3/dam]; None otherwise"""
        return getattr(self, "_collision_deficit_TS", None)

    @property
    def crest_lowering_TS(self):
        return self._crest_lowering_TS

    @property
    def volume_TS(self):
        return self._volume_TS
