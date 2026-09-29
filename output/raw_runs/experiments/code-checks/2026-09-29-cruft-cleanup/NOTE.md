# 2026-09-29: cruft removed from the cascade/ package. Hatteras results unchanged

Hannah, 2026-09-29: "remove the leftover cruft too" (`HATTERAS_CASCADE_CHANGES.md` section 6).

| item | done |
|---|---|
| `cascade/res_manager.py` | deleted: an empty file that nothing referenced and upstream does not have |
| `print`s in `roadway_manager.check_sandbag_need` | brought into line with upstream `main`. Removed the two local-only prints ("Road close enough for sandbags", "Road far enough away") and two commented-out debug prints; the "Sandbags would be added" message is shortened to upstream's form. The prints upstream also has ("Roadway relocated", "Elevation is low enough for sandbags", "Sandbags would be added…") are **kept**: removing them would add differences from upstream, not remove them. |
| formatting drift | `black` (the formatter in upstream's pre-commit) run on `cascade.py`, `cascade_groin.py`, `roadway_manager.py`, `beach_dune_manager.py`, `brie_coupler.py`, `chom_coupler.py` and `tools/plotters.py`. Black only rewrites layout, and it confirms the rewritten code parses to the same program before saving. `groin.py` was left out because another session is editing it. |

Not touched: `check_sandbag_need`'s *logic* differs from upstream (threshold 0.08 vs 0.101, all dune rows vs row 0). That is behaviour, not cruft, and sandbags are off in the Hatteras runs.

## Checks

- Every module compiles. Both `Cascade` classes run 3 years on their defaults.
- The 2010-2024 edgeBE full-management run, re-run here with the matrix's settings, is **bit-identical** to `raw_runs/matrix/2010_2024/edgeBE/HAT_2010_2024_edgeBE_offsetmetres_road_bdm_nourish_nogroin`: shoreline matrix, `shoreline_change_rate.csv`, `road_management.csv` and `nourishment_log.csv` all identical.

Log: `output/logs/scratch/codecheck_cruft_cleanup.log`.
