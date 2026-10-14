# matrix - the production runs

Only current windows live here. Superseded windows go to `D:\archive\CASCADE_offload\output\raw_runs\matrix\` (see the bottom).

**Re-run from 2026-10-13 with the new dune settings.** Done: `1996_2009/` and `2009_2025/`, zeroBE and edgeBE,
ends re-solved (1996 +0.9769 / +10.8017; 2009 +29.6576 / +20.7639 m/yr), and `1996_2025/` (full window,
zeroBE + edgeBE, ends +134.7 / +67.50 on its LRR target; GIS 1 noise-limited, best probe residual -0.026 m/yr).
The management ladder is re-run for all three windows and its figures redrawn. Not re-run: domainBE (waits on the
groin). New: in every managed 2009_2025 run the GIS 14 road drowns in model year 9 (20% of road-side cells below
0 m MHW); it did not under the 10-12 settings. Earlier the same day every run was deleted, not moved
(D: was not connected).
They used the 3.0 m MHW rebuild height and per-cell dune ceilings, both replaced that day
(design height 3.94 m MHW, one island-wide Dmaxel of 5.5 m NAVD88). Their rows are in
`../retired_runs.csv`; 38 of their small text files had never been committed, so not all are in
git history. `figures/` was kept and is drawn from the deleted runs. The windows and rules below
stand for the re-run.

## Current windows (the 2026-10-05 calibration/test plan)

| window | role | presets | target |
|---|---|---|---|
| `1996_2009/` | **calibration**: 1996–2008, ending 1 Jan 2009 | zeroBE, edgeBE (nogroin) | CoastSat net change from the 1996 offset window (1996-04-12..1997-04-12) to 2008 ±6 months (2007-07-01..2009-06-30), LOWESS-7 |
| `2009_2025/` | **test**: 2009–2025, ending 1 Jan 2026 (17 yr), fills on | zeroBE, edgeBE (nogroin) | CoastSat net change from the 2009 offset window (2009-02-17..2010-02-17) to 2025 ±6 months (2024-07-01..2026-06-30; CoastSat ends 2026-01-13), LOWESS-7 |
| `1996_2025/` | **full window**: one 30-yr run, through Dec 2025 | zeroBE, edgeBE (nogroin) | CoastSat LRR 1996-2025 (end solve) and net change from the 1996 offset window to the test end window; figures in `output/comparisons/full_window/` |

Rules since 2026-10-12 (the same in both periods): the island offset is the CoastSat mean over the DEM midpoint ±6 months (offset shoreline/v3); the target's start window is that same window; the end window is the last modelled year ±6 months. All storms are kept (the pre-DEM storm question is open with the advisor). The plan, scores and open decisions are in `scripts/hatteras_ms/DEM_TO_DEM_CALIBRATION.md`. Each run draws its own net-change figure under `<run>/figures/`.

**Management ladder (2026-10-12, edgeBE only):** each period also holds the edgeBE runs with less management: natural (`noroad_nobdm`), roads only (`road_nobdm`) and, for test and the full window, roads + beach/dune without fills (`road_bdm_nonourish`). All use the window's full-management end rates. The full-window rungs come from `HAT_full_window_1996_2025.py final --scenario`. Figures: `figures/<window>/management_ladder_{net_change,rate}_edgeBE_<window>.png` (`matrix_management_ladder.py`), scores in `figures/supporting/ladder_scores.csv`.

**Clean re-run 2026-10-12.** Every earlier matrix run (v2 ±1-yr offsets) was deleted at Hannah's request; their tracked text files are in git history. domainBE (BE set 1, blocking groin) is not re-run: set 1 and the groin pin were fitted on the old targets and wait for the groin decision.

## Moved off C: on 2026-10-07

`1996_2010/`, `2010_2024/` (the 14-yr windows, superseded 2026-10-03), `1996_2015/`, `2010_2026/` (the candidate windows, superseded by the calibration/test plan 2026-10-05) and `figures/` (management-ladder figures on 1996_2015 / 2010_2026). All of them are on `D:\archive\CASCADE_offload\output\raw_runs\matrix\`, same paths, and their tracked text files are in git history before the commit that removed them.
