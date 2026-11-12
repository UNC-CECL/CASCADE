"""
Village dune lines that follow the shoreline until the beach is first nourished.

CASCADE's BeachDuneManager holds a village's dune line fixed (Barrier3D's dune_migration_on False) until the beach
width reaches its threshold, which Cascade sets to 0 m for every village: the beach is assumed to be held by
nourishment, so shoreline change only narrows it. Hatteras had no nourishment before 2014 (Rodanthe; Buxton from
2017), and the 1996-2009 lidar shows the village dune lines moving landward with the shoreline as on the open coast
(dune-collision/1-observations/8-dune-change-by-survey/HAT_dune_line_position.py, 2026-11-12). So, by default since
2026-11-12, a village manager's threshold is out of reach until its domain is first nourished (the dune line then
migrates every year, as outside the villages), and CASCADE's own 0 m from the year after: the fixed line belongs to a
nourished beach.

The managers are switched to this subclass in place after Cascade builds them (build_cascade), so a saved model state
still pickles (the class lives in an importable module).
"""
from __future__ import annotations

import numpy as np
from cascade.beach_dune_manager import BeachDuneManager

UNTIL_NOURISHED_THRESHOLD_M = 1e9      # a beach this wide never exists: migration stays on
CASCADE_THRESHOLD_M = 0.0              # cascade.py's own, for every village


class MigrateUntilNourishedManager(BeachDuneManager):
    """A BeachDuneManager whose dune line migrates until its domain's first nourishment."""

    def update(self, barrier3d, nourish_now, rebuild_dune_now, nourishment_interval):
        t = barrier3d.time_index
        nourished = bool(np.nansum(np.asarray(self._nourishment_TS[:max(t - 1, 0)], float)) > 0)
        self._beach_width_threshold = CASCADE_THRESHOLD_M if nourished else UNTIL_NOURISHED_THRESHOLD_M
        return super().update(barrier3d, nourish_now, rebuild_dune_now, nourishment_interval)


def migrate_until_nourished(cascade):
    """Switches every beach/dune manager of a built Cascade to MigrateUntilNourishedManager; returns how many."""
    n = 0
    for manager in cascade.nourishments:
        if manager is not None and type(manager) is BeachDuneManager:
            manager.__class__ = MigrateUntilNourishedManager
            n += 1
    return n
