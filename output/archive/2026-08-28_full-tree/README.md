# Archived output — everything under `output/` as it stood on 2026-08-28

**Old results. Do not use anything in here for analysis.** This is the complete
`output/` tree, moved aside intact on 2026-08-28 so a fresh set of runs could
start against re-measured forcing without old and new results sharing a
directory. Nothing was deleted: every `.npz` model state, every figure and every
CSV is exactly as it was.

Live results are in the sibling directories — `raw_runs/`, `comparisons/`,
`groin_sweep/`, `sensitivity_analysis/`. If a number came from here, it predates
the 2026-08-28 forcing rebuild and the `5b26cb5` code changes.

```
comparisons/           56 MB   relocation + offset-mode figures
driver/               2.2 MB   run logs and driver manifests
groin_sweep/           28 MB   both periods, plus validation and figures
raw_runs/              12 GB   the runs themselves, and run_index.csv
sensitivity_analysis/  16 MB   three wave-sensitivity sweeps
```

## Why all of it went, not just the stale period

### `raw_runs/1984_2004` — 6 runs, definitively stale

`run_index.csv` records these at `topo_dune_version = v2`. That topography was
**deleted on 2026-08-27** when `1984-start` was cleared and re-extracted as a
new `v1` from the same DEM against a new pick set.

They also consumed the pre-2026-08-28 setbacks. `RoadSetback_1984_dunestart.csv`
was re-measured that morning against the live `v1` interior, and **15 of 83
road-bearing domains moved**, up to 25 m — including GIS 80, one of the three
roadways the relocation logic acts on, and GIS 10-13, where a 30→10 m shift
halves `int(setback/10)` and lands the road on a different row.

Two independent reasons; either alone would be enough.

### `raw_runs/2004_2024` — 30 runs, a subtler case

Their **inputs are still valid**. `2009_v5` in the index is only the pre-rename
name for today's `2004-start/v1` (renamed 2026-08-26, see
`1-barrier3d-domains/LINEAGE.md`), and `RoadSetback_2004_dunestart.csv` came
back **byte-identical** from the 2026-08-28 rebuild, verified by diff.

What moved is the **code**. They ran at `613d332b`, before `5b26cb5` changed:

| file | lines |
|---|---:|
| `cascade/groin.py` | 24 |
| `scripts/hatteras_ms/HAT_hindcast_1984_2024.py` | 1715 |
| `scripts/hatteras_ms/HAT_hindcast_config.py` | 511 |
| `scripts/hatteras_site_config.py` | 784 |
| `scripts/hatteras_ms/HAT_run_all.py` | 529 |

Every row in `run_index.csv`, both periods, also carries `git_dirty = True`, so
the recorded commit does not pin them even to that tree. Same inputs, different
model — which is exactly the old-vs-new ambiguity this move exists to remove.

### The derived products

`comparisons/`, `groin_sweep/` and `sensitivity_analysis/` are all built from
the runs above. Left in place they would have described data no longer in the
live tree, which is worse than not having them.

### The two older archives came too

`raw_runs/1984_2004/superseded_pre1996mosaic_20260827/` and
`superseded_presetbackfix_20260827/` are nested inside this one now. Their own
READMEs still explain what each was superseded by; they were already frozen and
nothing reads them.

## What replaced it

`output/` was rebuilt as an empty skeleton — `raw_runs/{1984_2004,2004_2024}/`,
`comparisons/`, `groin_sweep/`, `sensitivity_analysis/`, `driver/logs/` — with a
fresh `raw_runs/run_index.csv` holding **only the 44-column header**, byte-identical
to the archived file's header so nothing downstream has to change.

## One consequence in git

`.gitignore` ignores `output/*` wholesale and re-includes only
`output/raw_runs/`. This archive is not on that path, so **everything in here is
now ignored**. 331 files that were tracked or staged before the move show as
deletions in `git status`.

Of those, 113 were committed and remain recoverable from git history. The other
**218 were staged additions that had never been committed** — the
`superseded_pre1996mosaic_20260827` import — so for those this directory on disk
is the only copy. That is fine as long as this directory is not deleted, and it
is the reason "archive intact" was chosen over stripping the bulk.

If you ever want the small text products tracked again, move the run tree back
under `output/raw_runs/<start>_<end>/` — the re-include pattern is
`!output/raw_runs/[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]/`, so only
period-shaped directory names are matched.
