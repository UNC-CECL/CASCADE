# Source/sink (background erosion) — the zone rework

The source/sink field is being rebuilt step by step for the DEM-to-DEM periods (calibration 1996-2009, test 2009-2025), around zones on the net-change target. The scripts are in `scripts/input_prep/7-source-sink/`; their README tracks the steps. The field the model actually runs with is in `scripts/site_layer/hatteras_site_config.py`, not here.

## What is here

| path | what it is |
|---|---|
| `2-calibrate/1996_2009__2009_2025/zone_residuals_edgeBE.csv` | step 1: per period and GIS domain, the observed and modelled net change, and the residual (observed minus model) in m and m/yr over the run. The model is the edgeBE full-management run, no groin. Written by `2-calibrate/be_zone_residuals.py` |
| `3-figures/1996_2009__2009_2025/1-field/1-rate/` | step 1 as rates (m/yr over the run, what the zones are picked on): `zone_residuals_rate_edgeBE.png` (the residual per domain as bars, calibration over test), the same panels alone (`..._calibration_1996_2009.png`, `..._test_2009_2025.png`), and `zone_residuals_rate_overlap_edgeBE.png` (both periods in one panel, with a strip of the domains that are a source in both, a sink in both, or of opposite sign) |
| `3-figures/1996_2009__2009_2025/1-field/2-net-change/` | step 1 as net position change from the starting shoreline (m, what the target and its scores are in): the same four figures, `zone_residuals_net_change_*` |
| `2-calibrate/1996_2009__2009_2025/zone_candidates_edgeBE.csv` | step 2: one row per candidate zone, by set (shared, calibration, test): kind, GIS range, width, each period's mean and peak residual, and whether it still qualifies at ±1 m/yr. Written by `2-calibrate/be_zone_candidates.py` |
| `3-figures/1996_2009__2009_2025/2-method/zone_candidates_edgeBE.png` | step 2: the candidate zones over the residuals, shared, calibration only and test only |

Captions are in each figure folder's `supporting/CAPTIONS.md`.

The folder numbers match the script stages: tables in `2-calibrate/`, figures in `3-figures/`. Each period pair has its own folder, named `<p1start>_<p1end>__<p2start>_<p2end>`; resolve the paths through `scripts/site_layer/hat_source_sink.py`, never by typing them.

## Removed 2026-10-12

The earlier zone calibration's products: the 1984/2004 and 1996/2010 pairs' tables, figures and export, the `prebe/` config backups and `archive/superseded_20260914/`. They are in git history before 2026-10-12 (the figures as PDF only). The model never read them. The config still carries that calibration as the `calibBE` preset (1984 and 2004), and `be_apply_fit_to_config.py`, `plot_be_*.py` and `export_be_calibration.py` read the removed files, so they no longer run until the rework replaces them.
