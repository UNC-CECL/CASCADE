# matrix_vs_observed - every scenario in the run matrix against the observations

The run matrix (`output/raw_runs/matrix/`), every run scored and
drawn against the observed rate and position change. Nothing here runs the
model.

| script | question | writes to |
|---|---|---|
| `matrix_vs_observed.py` | each matrix run against the observations, plus scores | each run's `figures/vs_observed/`, `output/comparisons/matrix_vs_observed/` |
| `matrix_management_ladder.py` | how the fit changes as management is added, one rung per panel | `output/raw_runs/matrix/figures/` |

`matrix_management_ladder.py` imports its loaders (`matrix_runs`, `rates`,
`score`, ...) from `matrix_vs_observed.py`.

## The scripts in detail

Moved here from `scripts/analyze_output/README.md` on 2026-09-30, when the
scripts were grouped by question.

## matrix_vs_observed.py

Every option A matrix run against the observations, two figures per run
(Hannah, 2026-09-27). The runner draws only rate figures into each run's
`figures/`, so neither of these existed for the matrix. Other scripts
(`matrix_management_ladder.py`) import its loaders.

**Rate and position change** ("Where are these position plots?"), in the
layout the wave experiments use (`HAT_metres_2_wave_sensitivity_plot`):

- (a) *rate:* the model's OLS rate (`lrr_m_yr`) against the CoastSat LRR
  target, 7-domain LOWESS with the southern 10 raw, as scored;
- (b) *position:* the model's position change over the window, endpoint rate
  x 14 yr, against the observed CoastSat change: the mean position over the
  last calendar year minus that over the first, smoothed at 7 domains (10
  until 2026-09-28) (`5-scr/3-rates/coastsat/total_change/<w>/smoothed`).

**Start and end positions** ("the starting island position with the end
modeled position and the observed end position", both CoastSat and the dune
line). Drawn *relative to the start line* (Hannah chose this): the island's own
position swings ~6 km along the reach with its planform while the changes are
tens of metres, so on an absolute axis the four lines coincide. The start
position is the zero line; each end is metres seaward (+) or landward (-) of
it, per domain:

- *model end:* the run's last annual shoreline minus its first
  (`shoreline_matrix.npy`, sign flipped to seaward +);
- *CoastSat end:* the observed change per domain, unsmoothed (the
  total_change domain means, window 0), with its 7-domain LOWESS as a faint
  line;
- *dune-line end:* the runner's own end-year target, the dune-line change
  between the start and end vintages (`hindcast.build_shoreline_target`;
  1997 -> 2009 for 1996-2010, 2009 -> 2023 for 2010-2024). Its survey interval
  is not the calendar window (11.6 and 14.1 yr); it is reported in the
  caption, not rescaled.

**One y axis per panel type, across every run** (Hannah, 2026-09-27): rate
±7.5 m/yr (as the runner's real-domains rate figure,
`rerender_run_figures.py --ylim-real`); position change and start/end
positions each symmetric about zero, set from the widest run or observation
and rounded up to 10 m. `y_bounds.txt` says which.

Writes, in subfolders by figure, then window, then preset:

```
each matrix run's figures/vs_observed/
    rate_and_position_change.png, start_and_end_positions.png
output/comparisons/matrix_vs_observed/
    rate_and_position_change/<window>/
        scenarios_rate_and_position_<preset>_<window>.png   every scenario
        <preset>/rate_and_position_<preset>_<scenario>[_reloc]_<window>.png
    start_and_end_positions/<window>/<preset>/
        start_and_end_positions_<preset>_<scenario>[_reloc]_<window>.png
    scores.csv, y_bounds.txt, README.md
```

<details><summary>Function notes (the original docstrings)</summary>

**`observed_change()`**

```
CoastSat position change per domain, seaward +; window 0 is the
unsmoothed domain mean.
```

**`model_end_change()`**

```
Last annual shoreline minus the first, real domains, seaward + (x_s
increases landward, so the sign flips).
```

</details>

## matrix_management_ladder.py

The matrix runs of one window and preset in order of increasing management,
so the effect of each setting can be seen building up. Hannah, 2026-09-29:
"all of the runs in sequential order of increasing management so we can see
the effect of each setting"; then "remove the step panels ... show all of the
curves in each progressive step, also do a lrr version and position change
version".

**The ladder:** natural -> road management only -> + beach and dune
management -> + nourishment fills. Each rung is the matrix run that adds one
layer to the rung above. A rung the window does not have is skipped rather
than drawn twice: 1996-2010 has no fills (the driver skips `full_no_fill`
there as identical to `full_management`). What each rung adds is read from the
runs themselves (scenario and the fills the run applied), not typed. Left off
the ladder: beach/dune-only, a side branch, and the historical-relocation runs
(Hannah, 2026-09-29: "remove the historical relocations panel"; in 1996-2010
it matched full management to 0.0002 m/yr, since relocation moves the road,
not the shoreline; 2010-2024 has no relocation event).

**One panel per rung, cumulative:** panel k draws rungs 1..k, each in its own
shade of the house blue ramp (lighter = less managed), the newest rung
heaviest, over the observation. A rung keeps its colour in every panel, so a
curve can be followed down. The colours are fixed per layer in every window
(Hannah, 2026-09-29: "keep colours consistent"), so 1996-2010 (no fills) never
uses the darkest.

**Two versions:**

- *rate:* the OLS rate (`lrr_m_yr`) against the CoastSat LRR scoring target
  (7-domain LOWESS, raw means GIS 1-10), as the runner scores it; ±7.5 m/yr,
  as in the runner's own figure.
- *position:* the position change over the window (endpoint rate x 14 yr)
  against the observed CoastSat change (mean position in the last calendar
  year minus the first, smoothed at 7 domains); ±130 m, as in
  `output/comparisons/matrix_vs_observed/`.

Interior (GIS 2-89) RMSE and bias of the newest rung are in each row title;
every rung's are in `supporting/ladder_scores.csv`.

Writes `output/raw_runs/matrix/figures/<window>/management_ladder_{rate,position}_<preset>_<window>.png`,
the PDF and CAPTIONS.md under that folder's `supporting/`, and
`output/raw_runs/matrix/figures/supporting/ladder_scores.csv`. Each run also
deletes the earlier step-panel versions of the figure.

<details><summary>Function notes (the original docstrings)</summary>

**`ladder()`**

```
The rungs for one window and preset, least managed first, each adding
one or more layers to the one above; a run that adds nothing is dropped.
```

</details>
