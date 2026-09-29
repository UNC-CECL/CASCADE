# Which figure do I open?

*Written 2026-09-29 by `scripts/repo_tools/hat_write_figure_index.py`, which checks every path below against the disk. Re-run it after adding or renaming a figure.*

Three things to fix before reading any of these:

1. **Which window.** Five window folders sit as peers across two chains and one is context only — [`data/hatteras_init/5-scr/WINDOWS.md`](data/hatteras_init/5-scr/WINDOWS.md).
2. **Which estimator.** A rate turned into a distance is named by the window it was *fitted* on: **total** = same window, **projected** = carried onto another, **observed** = no rate — [`data/hatteras_init/5-scr/3-rates/README.md`](data/hatteras_init/5-scr/3-rates/README.md).
3. **Which units.** `output/comparisons/model_vs_observed/` is in **m/yr**. Every other comparison tree is in **metres**. The same comparison exists in both, and they are not interchangeable.

Every figure carries its quantity, window and method in its title, and a full caption in `supporting/CAPTIONS.md` beside it.

## What did the shoreline do?

CoastSat satellite waterline, observations only. No model anywhere in these.

| for | open | note |
|---|---|---|
| The rate, m/yr | `data/hatteras_init/5-scr/3-rates/coastsat/lrr/<w>/lrr_<w>.png` <br>*windows:* 1996_2010, 2010_2024, 1996_2024, 1984_2004, 2004_2024 | OLS slope through every satellite position in the window. **This is the model's scoring target.** |
| That rate as a distance, m | `data/hatteras_init/5-scr/3-rates/coastsat/total_change/<w>/total_change_<w>.png` <br>*windows:* 1996_2010, 2010_2024, 1996_2024 | the window's own rate x its own years, against the observed change. TOTAL change — nothing extrapolated. |
| The long-term rate applied to a half, m | `data/hatteras_init/5-scr/3-rates/coastsat/projected/<w>/projected_<w>.png` <br>*windows:* 1996_2010, 2010_2024 | the 1996–2024 rate x 14 yr, on a window it was NOT fitted on. PROJECTED. 1996_2010 and 2010_2024 only. |
| Two snapshots differenced, m and m/yr | `data/hatteras_init/5-scr/3-rates/coastsat/endpoint/<w>/coastsat_endpoint_<w>.png` <br>*windows:* 1996_2010, 2010_2024, 1996_2024, 1984_2004, 2004_2024 | mean position ±6 months about each dune-line date. No rate fit. |
| The rate in 5-year bins | `data/hatteras_init/5-scr/3-rates/coastsat/5yr_bins/<w>/lrr_5yr_bins_<w>.png` <br>*windows:* 1996_2010, 2010_2024, 1996_2024 | is the trend steady inside the window? |
| How much the alongshore smoothing changes it | `data/hatteras_init/5-scr/3-rates/coastsat/total_change/<w>/smoothed/` <br>*windows:* 1996_2010, 2010_2024, 1996_2024 | and `projected/<w>/smoothed/`. Read the bias, not r — a smoother inflates r on both sides. |
| The full-record rate, Cape Point to Virginia, m/yr | `data/hatteras_init/5-scr/3-rates/coastsat/lrr/1984_2025_obx/lrr_obx_1984_2025.png` | 1984–2025, per transect, 153 km of coast, well past the model domain; with an overview and four regional maps (`..._map_*.png`) beside it. A hand-off to the Murray lab (2026-09-27), NOT a model window or target. |

## What did the dune line do?

The digitized dune line, observations only.

| for | open | note |
|---|---|---|
| Net change between the two lines | `data/hatteras_init/5-scr/3-rates/duneline/endpoint/<w>/duneline_endpoint_<w>.png` <br>*windows:* 1996_2010, 2010_2024, 1996_2024, 1984_2004, 2004_2024 | end line minus start line. Measured, never fitted. |
| Where each line actually sat | `data/hatteras_init/5-scr/4-comparisons/duneline_positions/` | maps, imagery zooms, distance to NC-12, beach width. |

## Did the dune line move with the shoreline?

Both observations on one panel, the gap between them shaded as beach-width change. **Metres.**

| for | open | note |
|---|---|---|
| Shoreline as two snapshots | `data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/coastsat_endpoint_vs_duneline_endpoint/<w>/coastsat_endpoint_vs_duneline_<w>_alongshore.png` <br>*windows:* 1996_2010, 2010_2024, 1996_2024, 1984_2004, 2004_2024 | observed vs observed; `..._scatter.png` beside it. |
| Shoreline as its OWN window's trend | `data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/coastsat_total_change_vs_duneline_endpoint/<w>/coastsat_total_change_vs_duneline_<w>_two_panel.png` <br>*windows:* 1996_2010, 2010_2024, 1996_2024 | also `_shaded_gap` and `_overlay`; the dune-interval value is a column in `domain_comparison.csv`. |
| Shoreline as the LONG-TERM trend carried onto a half | `data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/coastsat_projected_vs_duneline_endpoint/<w>/coastsat_projected_vs_duneline_<w>_two_panel.png` <br>*windows:* 1996_2010, 2010_2024 | the 1996–2024 LRR × 14 yr against the dune line measured over that half. 1996_2010 and 2010_2024 only. |
| One long-term prediction vs two dune-line outcomes | `data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/coastsat_projected_vs_duneline_endpoint/all_windows_stacked/coastsat_projected_vs_duneline_1996_2010_2024_halves_overlay.png` | the two halves stacked. The shoreline side is IDENTICAL in both panels, so every difference between them is the dune line's. Read beside the `coastsat_total_change_...` sheet of the same name. |
| The same two halves, each on its OWN rate | `data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/coastsat_total_change_vs_duneline_endpoint/all_windows_stacked/coastsat_total_change_vs_duneline_1996_2010_2024_halves_overlay.png` | the counterpart of the row above. The pair separates what the long-term trend PREDICTS from what it was FITTED on. |
| Does smoothing change any of it? | `data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/smoothed_loess7/` | both sheets again with BOTH curves LOESS-smoothed at 7 domains (3.5 km). Read the beach width, not r — a symmetric smoother inflates r on both sides. |
| The whole period above its two halves | `data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/coastsat_total_change_vs_duneline_endpoint/all_windows_stacked/coastsat_total_change_vs_duneline_1996_2010_2024_stacked.png` | three panels, shoreline and dune line as lines with the gap shaded. A different question from the two-panel sheets above. |
| Did they change pace together? | `data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/coastsat_total_change_vs_duneline_endpoint/change_between_periods/coastsat_total_change_vs_duneline_change_between_periods.png` | second half minus first half, on both sides. Only for the total product: the projected one uses the same rate in both halves, so its difference is zero by construction. |

## How does the model compare to the observations?

**These are in m/yr, not metres** — the one tree that is. Sliced by which observation, not by estimator; the cross-reference below says which estimator each folder uses.

| for | open | note |
|---|---|---|
| vs the CoastSat shoreline | `output/comparisons/model_vs_observed/vs_shoreline/domain_means/model_vs_shoreline_means_<w>.png` <br>*windows:* 1996_2010, 2010_2024, 1984_2004, 2004_2024 | `_grid.png` puts all four windows on one sheet. |
| vs the shoreline, as graded | `output/comparisons/model_vs_observed/vs_shoreline/smoothed/model_vs_shoreline_smoothed_<w>.png` <br>*windows:* 1996_2010, 2010_2024, 1984_2004, 2004_2024 | the LOESS form the runner actually scores. |
| vs the dune line | `output/comparisons/model_vs_observed/vs_duneline/endpoint_net_change/model_vs_duneline_netchange_<w>.png` <br>*windows:* 1996_2010, 2010_2024, 1984_2004, 2004_2024 | and `net_change_smoothed/` beside it. |
| vs both at once | `output/comparisons/model_vs_observed/vs_shoreline_and_duneline/change_rate/model_vs_shoreline_and_duneline_rate_<w>.png` <br>*windows:* 1996_2010, 2010_2024, 1984_2004, 2004_2024 | each solve drawn in its own target's estimator; the same in metres (x window years) under `net_change/`. |
| Does the answer survive changing the estimator or the solve? | `output/comparisons/model_vs_observed/sensitivity/` | `ends-swapped`, `dune-raw-solve`, `mixed-estimator`. The arm is in every filename. |

## Which target should the model be graded on?

CoastSat against the dune line, with the runs, as **net change in metres** over each 14-yr window.

| for | open | note |
|---|---|---|
| Start here | `output/comparisons/target_comparison/projected/paired/target_and_own_run_projected_<w>.png` <br>*windows:* 1996_2010, 2010_2024 | each target with its own run. **`projected/` is the target in use**: the 1996–2024 LRR x 14 yr. |
| The same on each window's own rate | `output/comparisons/target_comparison/total_change/paired/target_and_own_run_total_change_<w>.png` <br>*windows:* 1996_2010, 2010_2024 | kept for the record — what the runner grades against. |
| Neither target fitted anywhere | `output/comparisons/target_comparison/projected/ends_unsolved/` | zeroBE: no source/sink term in any domain, so all 90 are the model's own response. |
| Does the grading window matter? | `output/comparisons/target_comparison/smoothing_scale/projected_vs_model_<w>.png` <br>*windows:* 1996_2010, 2010_2024 | **No** — and r was never the number to read. See its PROVENANCE.md. |
| Both targets AND the model, smoothed | `output/comparisons/target_comparison/smoothed_loess7_with_cascade/` | the two smoothed sheets with the zeroBE run over them in dark green. Nothing in that run was fitted to either target, so all 90 domains are the model's own response. |
| The numbers | `output/comparisons/target_comparison/projected/tables/skill.csv` | bias, RMSE and r per window, model set and target. |
| The same, as rates (m/yr) | `output/comparisons/target_comparison/projected/change_rate/paired/target_and_own_run_projected_rate_<w>.png` <br>*windows:* 1996_2010, 2010_2024 | every figure above divided back by 14 yr: the CoastSat LRR, the dune line's measured rate (no scaling to the window), the model's endpoint rate. `_rate` in every stem; same layout, own `tables/skill.csv` in m/yr. |
| ...on each window's own rate | `output/comparisons/target_comparison/total_change/change_rate/paired/target_and_own_run_total_change_rate_<w>.png` <br>*windows:* 1996_2010, 2010_2024 | the m/yr twin of `total_change/`. |

## Cross-reference: what `model_vs_observed/` actually plots

That tree is sliced by which **observation** the model is held against, while `3-rates/` and `4-comparisons/` are sliced by **estimator**. Same underlying observations, different question, so the folder names do not line up. This is the translation:

| folder | observation | model estimator | reading |
|---|---|---|---|
| `vs_shoreline/domain_means` | CoastSat | OLS rate (`lrr_m_yr`) | raw domain means |
| `vs_shoreline/smoothed` | CoastSat | OLS rate (`lrr_m_yr`) | spliced LOESS — the form the runner grades |
| `vs_duneline/endpoint_net_change` | dune line | endpoint rate (`change_rate_m_yr`) | raw domain means |
| `vs_duneline/net_change_smoothed` | dune line | endpoint rate (`change_rate_m_yr`) | spliced LOESS |
| `vs_shoreline_and_duneline/change_rate` | both | endpoint rate, both solves | raw domain means |
| `vs_shoreline_and_duneline/net_change` | both | net change in metres (rate x window years), both solves | raw domain means |
| `sensitivity/mixed-estimator` | dune line | **OLS** rate against an **endpoint** observation | deliberately mismatched, as a sensitivity |

All of it in **m/yr** except `vs_shoreline_and_duneline/net_change/`, the combined panel in metres. To see the full comparison as a distance in metres, use `output/comparisons/target_comparison/` (its m/yr twin is under each target mode's `change_rate/`).

## Where the rest lives

| | |
|---|---|
| Finished figures for the paper | `output/figures/<subject>/` |
| How the repo is laid out | [`ORGANIZATION.md`](ORGANIZATION.md) |
| Figure house style | `scripts/site_layer/hat_figure_style.py`, `figure_making/STYLE.md`, `output/figures/style/` |
| Model runs | `output/raw_runs/`, indexed in `run_index.csv` |
