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
`../figures/6-result/island/HAT_compare_v2_v3_<pair>.png` (setback at 1984;
every relocation by year against the recorded events; relocations per year;
interior width at 1984 and 2004; island-mean width and cumulative overwash
through time; mean interior elevation at 2004 and the v3 − v2 difference) and
`version_compare_<pair>.csv` (one row per domain, both versions and the
difference), plus `HAT_compare_versions.txt` with the skill, the relocation
counts and timing, and the island medians. Nothing is re-run or re-scored.

Why the emergent pair was re-run rather than read from the calibration tree:
the v2 run there (2026-09-07) and the v3 run in `behindroad-copy` (09-08) sat
on different commits, with the pipeline and the live 1984 setback CSV changed
between them, so their differences would not have been the topography's alone.
