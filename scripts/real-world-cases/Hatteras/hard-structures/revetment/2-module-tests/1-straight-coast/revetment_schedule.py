"""
The simple revetment: a fixed dune line, a dune-height floor, and a shoreline held at the wall.

    from revetment_schedule import RevetmentSchedule
    wall = RevetmentSchedule(domains=[6, 7, 8], build_year=1, crest_height_m=1.2)
    wall.before_update(cascade); cascade.update(); wall.after_update(cascade)

Applied by the runner around each cascade.update(), like the nourishment schedule;
CASCADE, BRIE and Barrier3D are not edited. The wall is repaired, so it never
fails: once the shoreline reaches it, the shoreline is put back at the wall each
year through CASCADE's nourishment path, and the sand that takes is tallied.
Reasons in ../../README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-12
"""

from __future__ import annotations

from cascade.beach_dune_manager import shoreface_nourishment


# One wall over a block of domains; holds state and a per-year record
class RevetmentSchedule:
    def __init__(self, domains, build_year, crest_height_m, start_gap_m=None):
        self.domains = list(domains)
        self.build_year = build_year            # Barrier3D time_index of the first held year
        self.crest_height_m = crest_height_m    # above the berm
        self.start_gap_m = start_gap_m          # None = Barrier3D's implicit beach width
        self.gap_m = {}
        self.records = []

    # Before the year: build if due, then hold the dune line and floor the crest
    def before_update(self, cascade):
        for d in self.domains:
            b3d = cascade.barrier3d[d]
            if b3d.time_index < self.build_year:
                continue
            if d not in self.gap_m:
                self._build(b3d, d)
            b3d.dune_migration_on = False
            b3d._DuneRestart = self.crest_height_m / 10

    # After the year: move the gap with the shoreline, hold the shoreline at the wall, record
    def after_update(self, cascade):
        pinned_any = False
        for d in self.domains:
            b3d = cascade.barrier3d[d]
            if d not in self.gap_m:
                continue
            self.gap_m[d] -= (b3d.x_s_TS[-1] - b3d.x_s_TS[-2]) * 10
            pinned_m, sand_m3 = 0.0, 0.0
            if self.gap_m[d] < 0:
                pinned_m, sand_m3 = self._pin(b3d, -self.gap_m[d])
                self.gap_m[d] = 0.0
                pinned_any = True
            t = b3d.time_index - 1
            crest = b3d.DuneDomain[t, :, :].max(axis=1)
            self.records.append(dict(
                year=t, domain=d, gap_m=self.gap_m[d], pinned_m=pinned_m,
                sand_added_m3_per_m=sand_m3, crest_min_m=float(crest.min()) * 10,
                crest_mean_m=float(crest.mean()) * 10, x_s_m=b3d.x_s * 10,
            ))
        if pinned_any:
            self._sync_brie(cascade)

    # Place the wall: raise low dune cells to the crest, set the starting gap
    def _build(self, b3d, d):
        t = b3d.time_index - 1
        dunes = b3d.DuneDomain
        dunes[t][dunes[t] < self.crest_height_m / 10] = self.crest_height_m / 10
        b3d.DuneDomain = dunes
        if self.start_gap_m is None:
            self.gap_m[d] = int(b3d._BermEl / b3d._beta) * 10
        else:
            self.gap_m[d] = self.start_gap_m

    # Put the shoreline back at the wall with the sand a nourishment would need
    def _pin(self, b3d, overshoot_m):
        height = 2 * b3d.h_b_TS[-1] + b3d.DShoreface        # active profile, dam
        volume = overshoot_m / 10 * height / 2              # dam^3/dam
        b3d.x_s, b3d.s_sf_TS[-1], _ = shoreface_nourishment(
            b3d.x_s, b3d.x_t, volume, b3d.h_b_TS[-1], b3d.DShoreface, 0)
        b3d.x_s_TS[-1] = b3d.x_s
        return overshoot_m, volume * 100                    # m, m^3/m

    # Tell BRIE about the moved shorelines, as CASCADE does after a nourishment
    def _sync_brie(self, cascade):
        b3ds = cascade.barrier3d
        cascade._brie_coupler.update_brie_for_human_modifications(
            [b.x_t_TS[-1] for b in b3ds], [b.x_s_TS[-1] for b in b3ds],
            [b.x_b_TS[-1] for b in b3ds], [b.h_b_TS[-1] for b in b3ds],
            [b.s_sf_TS[-1] for b in b3ds])
