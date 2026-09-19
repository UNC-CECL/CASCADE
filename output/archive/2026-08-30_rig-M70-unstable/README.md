# This directory's name and its contents disagree — read before using it

Written 2026-09-02. **Nothing was renamed or deleted**, because the contents
cannot be attributed with confidence and guessing would be worse than the
current mess. Do not cite a number out of here without resolving the below.

## What is wrong

The directory says `rig_M70_f0.6_UNSTABLE_sweep_leftover`. Every file inside is
named `HAT_1967_2018_edge_calibrated_groin_*`, which is the name of a
**different run** — the one in the sibling directory
`HAT_1967_2018_edge_calibrated_groin/`.

It is not a copy of that sibling. The `.npz`, the `_shoreline_matrix.npy` and
the `_groin_diagnostics.csv` here all differ byte-for-byte from the sibling's,
so this is real output from a real run, wearing another run's filenames.

## Three timestamps, not one run

    14:21  HAT_1967_2018_edge_calibrated_groin_PLOT_*.png/.gif   (6 files)
    15:42  HAT_1967_2018_edge_calibrated_groin.npz               (4 data files)
           ...._shoreline_matrix.npy, ...._groin_diagnostics.csv,
           ...._historical_BN_log.csv
    21:16  PLOT_*.png                                            (5 files)

All 2026-08-30. **The prefixed figures predate the data files beside them by
81 minutes**, so they cannot have been drawn from that `.npz`. The unprefixed
`PLOT_*.png` set at 21:16 is a third thing again, and differs from the prefixed
set of the same names.

So the directory holds at least two, probably three, points in time. A run
directory in this project is flat and holds exactly one run's output — that is
what `guard_run_dir` enforces — and this one predates that rule.

## What would settle it

`_groin_diagnostics.csv` records per-year trapping but **not** M or f, so the
data files cannot be attributed from their own contents. The candidates are the
groin sweep cells under `output/groin_sweep/1984_2004_edgeBE/M70_*`, and the
rig runs described in `hard-structures/groin/GROIN_PLAN.md`. Comparing this
`_shoreline_matrix.npy` against those would identify it.

## Why it was kept

`UNSTABLE` in the name suggests it was retained as a negative result. Nothing in
`scripts/` or `hard-structures/` references this directory — zero matches for
either `UNSTABLE_sweep_leftover` or `rig_M70` — so nothing breaks if it is
resolved or removed. The decision is whose ever remembers what the M = 70 run
was for.
