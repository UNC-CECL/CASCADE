# 2026-10-03: comparison figures on the 14-yr windows

**Do not use for analysis.** These are the comparison outputs drawn on the 1996-2010 and
2010-2024 windows (the matrix of 09-28/09-29). They were replaced when the model windows
moved to 1996-2015 and 2010-2026 (decided 2026-10-02) and the matrix was re-run on 10-03
(2017 Buxton fill, 1.62 M cy Rodanthe).

What was replaced, and by what:

| here | was | replaced by |
|---|---|---|
| `comparisons/hindcast_calibrated/*.png` | `output/comparisons/hindcast_calibrated/` | the same names, 1996-2015 + 2010-2026, model smoothed like the target |
| `comparisons/model_vs_observed/vs_shoreline/` | `output/comparisons/model_vs_observed/vs_shoreline/` | the same folder on the new windows |
| `comparisons/model_vs_observed/{tables,runs_used.csv,y_bounds.txt}` | the same, COPIED | they also score the dune-line figures that were left in place |
| `comparisons/matrix_vs_observed/` | `rate_and_position_change/` and `start_and_end_positions/` `{1996_2010,2010_2024}`, `scores.csv`, `y_bounds.txt` | the `1996_2015`/`2010_2026` folders |
| `comparisons/target_comparison/total_change/` | `output/comparisons/target_comparison/total_change/` | `total_change/runs_vs_coastsat/` (CoastSat only) |
| `comparisons/scenario_grid/` | the png and a stale 09-11 CAPTIONS.md | the png and a new CAPTIONS.md |
| `raw_runs_matrix_figures/` | `output/raw_runs/matrix/figures/{1996_2010,2010_2024}`, `supporting/ladder_scores.csv` | the `1996_2015`/`2010_2026` ladders |

Not archived, because nothing replaces them: the dune-line folders of `model_vs_observed/`
(`vs_duneline/`, `vs_shoreline_and_duneline/`, `sensitivity/`), and `target_comparison/`
`projected/`, `smoothing_scale/` and `smoothed_lowess7_with_cascade/`. There is no dune line
for 2015 or 2026, and the 1996-2024 LRR solve exists on the 14-yr windows only. They stay in
place as the 14-yr record. The scripts reproduce them with `--windows 14yr` /
`select_windows("14yr")`.

Methods that changed between these figures and their replacements: the model line is now
smoothed like the observation wherever the observation is smoothed (Hannah, 2026-10-03), and
the target_comparison CoastSat target on the new windows is each window's own LRR, not the
1996-2024 LRR x 14.
