# 6-result — what the model does with the reconstruction

**v2 against v3 under the same hindcast** (Hannah's advisor, 2026-09-09: a
side-by-side under the modules' automatic behaviour, full management,
calibrated). Two run pairs, each run on both versions by
`HAT_run_version_pair.py` (scripts `6-result/`) on the same code the same day,
into `output/raw_runs/version-pair/<version>/1984_2004/calibBE/`:

| pair | run | what |
|---|---|---|
| emergent | `HAT_1984_2004_calibBE_road_bdm_groin` | full management, calibBE, groin on; the roadway and beach-dune modules decide on their own |
| prescribed | `HAT_1984_2004_calibBE_road_reloc_bdm_groin` | the same with the recorded 1989 (GIS 84-87) and 1999 (GIS 9-14) relocations imposed: the control |

`HAT_compare_versions.py` reads the two runs' saved state and writes, per pair,
two figures in `../figures/6-result/island/` — `HAT_compare_v2_v3_<pair>_relocations.png`
(the setback each version starts with; every relocation by year against the
recorded events; relocations per year) and `HAT_compare_v2_v3_<pair>_geometry.png`
(interior width per domain at 1984 and 2004; island-mean width through time and
the difference between versions on its own axis; cumulative overwash; mean
interior elevation at 2004 and, below it, the difference) — and
`version_compare_<pair>.csv` (one row per domain, both versions and the
difference), plus `HAT_compare_versions.txt` with the skill, the relocation
counts and timing, and the island medians, and a `CAPTIONS.md` beside the
figures. Nothing is re-run or re-scored. The two were one seven-panel figure
until 2026-09-10, when they were split and sized to a 190 mm column.

Why the emergent pair was re-run rather than read from the calibration tree:
the v2 run there (2026-09-07) and the v3 run in `behindroad-copy` (09-08) sat
on different commits, with the pipeline and the live 1984 setback CSV changed
between them, so their differences would not have been the topography's alone.

**Through time, side by side.** `HAT_version_pair_gif.py` (same scripts folder)
renders the relocation comparison's animations with the two panels being v2
and v3 under one scenario, into
`output/comparisons/relocation_1984_2004/v2_vs_v3/<scenario>/<place>/`, one
folder per place (the island, the two event blocks, Pea Island where rows are
added, Avon to Tri-Village where they are removed) with `topography.gif` and
`dune-and-road.gif` in each. That is where the inserted and removed cells can
be watched doing something.

**The relocation comparison, side by side.** `HAT_version_pair_report.py` (same
scripts folder) reads the v2 and v3 relocation-comparison sets
(`output/comparisons/relocation_1984_2004/<version>/calibBE_groin/tables/`)
and the v3 footprint audit, and writes
`output/comparisons/relocation_1984_2004/v2_vs_v3/report.txt` with `tables/`
beside it: every section of a per-version report with a v2 column, a v3
column and the difference, headed by the identity of the four runs. Nothing
is re-run or re-scored; it must be re-run after `HAT_relocation_comparison.py`
is re-run on either version.
