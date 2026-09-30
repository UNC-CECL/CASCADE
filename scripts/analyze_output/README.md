# analyze_output - questions that span more than one run

Nothing here runs the model. Each script reads finished runs and writes to
`output/comparisons/`, resolved as `hat_figure_style.COMPARISONS_ROOT`
(2026-09-18); a figure finished for the manuscript also goes to
`output/figures/` (numbered layout, map in its README.md), from `scripts/figure_making/`.

```
compare_runs/
    rate_windows.py      LIVE. Observed vs modelled shoreline-change rate,
                             every window, the three end-solve model sets
                             -> comparisons/model_vs_observed/
    compare_runs.py      a general run-vs-run / run-vs-CoastSat tool.
                             RUNS_TO_COMPARE is EMPTY: every example in it is
                             commented out. Fill it before running
                             -> comparisons/<COMPARISON_NAME>/
overwash/
    compare_overwash_figures.py   a run's overwash (Qow) against the observed
    compare_overwash_observed.py  record; both resolve their run through
                                  run_registry (HAT_1984_2004_calibBE_road_bdm_groin,
                                  read from raw_runs/archive/2026-09-24-pre-metres/
                                  since 2026-09-30) -> comparisons/overwash/
    superseded_20260918/          plot_overwash.py (dead absolute path) and an
                                  early copy of Roya's Pea Island script
smoothing_vs_cascade/
    smoothing_vs_cascade.py   what the LOWESS smoothing does to the
                             comparison. Names HAT_1984_2004_SQ_BE_Hs2p0, which
                             no longer exists anywhere under raw_runs/: point it
                             at a current run before running
```

The observed overwash record and its own figures are
`scripts/input_prep/8-overwash-analysis/`; this folder only compares a RUN
against it.

The smoothing figures and the source/sink-zone figures were deleted from
`output/comparisons/` on 2026-09-17 (unregenerable, uncited); running either
script recreates its folder, so give it a current run first.

The run tree is addressed through `cascade_pipeline.run_registry`, never by
building a path: a run's NAME describes its scenario and its PATH describes its
forcing, and the two are easy to get wrong by hand.

## The scripts in detail

Each script's header says what it does and how to run it; the reasoning, the
choices behind it and its history are here, one section per script. Moved out
of the scripts on 2026-09-30, when they were brought in line with
`scripts/STYLE.md`.

Shared history: every script here drew in matplotlib's defaults until
2026-09-17, when each started calling `apply_style()` from
`site_layer/hat_figure_style.py` at import. Their paths used to be absolute
literals typed on one machine; they are all found by searching upward for the
repo root now (ORGANIZATION.md rule 5), so they follow the checkout and survive
a file changing depth.

### compare_runs/adoption_before_after.py

Hannah, 2026-09-28: "show me the comparison when the matrix finishes". It
compares the matrix before and after the 2026-09-28 adoption:

| side | runs | Barrier3D and forcing |
|---|---|---|
| before | `output/raw_runs/archive/2026-09-28-pre-ceiling/matrix/` | `fix/route-overwash-axis-swap` (49fd069), Dmaxel default (3.4 m NAVD88), storms `v3_72`, the pre-adoption LOWESS-7 ends |
| after | `output/raw_runs/matrix/` | `hatteras/adopted` (overwash fixes + per-cell dune ceilings), storms `v3_trim24`, the ends re-solved on it (`end-domain-boundaries/2026-09-28-ends-resolved-adopted`) |

For every matrix run (both windows, both presets, every scenario) it scores:

- **shoreline:** interior (GIS 2-89) RMSE and bias of the model LRR against the
  CoastSat LOWESS-7 target (`run_registry.skill_vs_target`, as the runner
  scores), and the spatial correlation r;
- **overwash:** against the imagery (`8-overwash-analysis`), each run dated by
  its own storm file: POD, POFD, PSS, timing r, space r;
- **dunes:** for 1996-2010 runs, the 2010 dune crest minus the 2009 lidar.

Beside the scores it writes per-domain tables: `cells_<side>.csv` (every image
x domain, observed and model overwash) and `crest_<side>.csv` (end-of-run crest
per GIS domain, m MHW). `crest_lidar_2009.csv` is the crest from the 2010-start
dune file. The figures are `adoption_before_after_figures.py`.

**Each side is scored under the Barrier3D it ran on**, so the storm sharing
uses that version's DuneGaps and DuneGrowth. `score --side before` must run
with `PYTHONPATH=<Barrier3D at 49fd069 + the ceiling feature, off>` (the
worktree `../Barrier3D-dune-ceiling`); `score --side after` uses the editable
install. The script refuses to score a side under the wrong one.

### compare_runs/adoption_before_after_figures.py

Hannah, 2026-09-28: "make figures of the before and after comparison". It
reads the tables `adoption_before_after.py` writes (`scores_`, `cells_`,
`crest_<side>.csv`, `crest_lidar_2009.csv`) and the runs' own
`shoreline_change_rate.csv`. edgeBE only: zeroBE is within 0.1 of it on every
score (see the README beside the tables). The relocation arms are left out;
in both windows they score the same as their non-relocation twins.

| figure | shows |
|---|---|
| `adoption_scorecard.png` | every score, before -> after, per scenario |
| `adoption_shoreline_alongshore` | model LRR vs CoastSat LOWESS-7, managed + natural |
| `adoption_overwash_map_<scenario>` | image x domain: hit / miss / false alarm |
| `adoption_overwash_by_image` | domains overwashed per image, grouped bars |
| `adoption_dune_crest_2010` | the 1996-2010 runs' 2010 crest vs the 2009 lidar |

Output: `output/comparisons/adoption_2026-09-28/figures/`. On each run it also
deletes the retired `adoption_overwash_by_domain` figure if one is left over.

### compare_runs/matrix_management_ladder.py

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

### compare_runs/matrix_vs_observed.py

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

### compare_runs/offset_source_comparison.py

Asked by Hannah on 2026-09-28: a simple comparison of the model output
started from the dune-line offset against the model output started from the
shoreline offset, for 1996-2010 and 2010-2024. The main question is how much
the island's *orientation* in the offset affects the outcome.

**Runs** (no new runs): the full-management pair of each period from
`output/raw_runs/experiments/island-offset/2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a/runs/{duneline,shoreline}_full_management/<period>/zeroBE/<run>/`.
Option A waves, zeroBE ends, relocations and groins off; the only thing that
differs within a pair is the island offset.

**Offsets:** `2-brie-offset/<start>/{duneline,shoreline}/<version>/*_unpadded.csv`,
at the version each run's metadata records. BRIE adds the offset to x_s, so
larger = more *landward*; everything uses the *seaward* position s = -offset,
with the mean removed (a uniform shift does not change what BRIE does, and the
builds are not on a common datum).

**Quantities,** per GIS domain (500 m), shoreline start minus dune-line start:

- orientation, theta = atan(ds/dx), degrees;
- turning, d(theta)/dx, degrees per km; positive where the line bends
  landward (an embayment), negative at a bulge;
- model change, each run's LRR x 14 yr (m, seaward positive).

Statistics are on the interior, GIS 2-89.

**Target** (second version of the change figure, Hannah 2026-09-29):
projected shoreline change, the CoastSat LRR fitted on 1996-2024, LOWESS over
7 domains (southern 10 raw), x 14 yr: one profile, the same in both periods.
The model stays unsmoothed.

**Label placement** (Hannah, 2026-09-28): every label stays out of the data.
The villages and shoals are named in a strip above the highest line, and the
groin and piers are drawn below that strip and named in the legend. Draw
order, bottom to top: village and shoal bands, grid, zero line, groin and
piers, the two runs, labels; the grid has no ticks inside the label strip.

Output, `output/comparisons/offset_source/`:

| file | shows |
|---|---|
| `offset_source_model_change_full_management.png` | the two runs, both periods |
| `offset_source_model_change_vs_projected_full_management.png` | the same with the projected target on top (2026-09-29) |
| `offset_source_difference_full_management.png` | profiles |
| `offset_source_orientation_vs_model_full_management.png` | scatter |
| `tables/summary.csv`, `per_domain.csv`, `vs_projected.csv` | the numbers |

<details><summary>Function notes (the original docstrings)</summary>

**`projected_target()`**

```
Projected shoreline change (m): the 1996-2024 CoastSat LRR target, built
as the runner builds it at LOWESS_DOMAINS, x YEARS.
```

**`vs_target()`**

```
Each run against the projected target, interior GIS 2-89: bias and RMS
residual are the numbers to read; explained and r beside them.
```

**`fig_change_only()`**

```
Panels (a, b) of the profile figure on their own: the two runs' change.
With `obs`, the second version: the projected target drawn on top of them.
Every label stays out of the data (Hannah, 2026-09-28): the villages and
shoals are named in a strip above the highest line, and the groin and
piers are drawn below that strip and named in the legend, not on the lines.
```

</details>

### compare_runs/smoothed_lowess7_with_cascade.py

The two smoothed halves-overlay sheets with the CASCADE hindcast drawn over
them in dark green, so the model's alongshore behaviour can be read against
both candidate targets at the scale the model resolves. Built 2026-09-22
(Hannah, by interview).

It is the observations-only pair in
`5-scr/4-comparisons/shoreline_vs_duneline/smoothed_lowess7/` plus one line.
Everything about the two observed curves (the 7-domain LOWESS, the GIS 1-10
raw splice, both sides smoothed at transect resolution, the faint raw domain
means behind them) is imported from that script rather than re-implemented,
so the sheets differ in exactly one thing: the green line.

**Why this lives in `output/` and not beside its twin** (Hannah, 2026-09-22):
`data/hatteras_init/5-scr/4-comparisons/` is observations only; a figure
carrying model output is a product, and ORGANIZATION.md rule 1 puts products
in `output/`. `target_comparison/` already holds exactly this kind of figure
(both candidate targets with the hindcast over them), so this is its
smoothed, two-panel sibling.

**The run:** zeroBE, full management, groin off, one per period:
`HAT_1996_2010_zeroBE_offsetmetres_road_bdm_nogroin` and
`HAT_2010_2024_zeroBE_offsetmetres_road_bdm_nourish_nogroin`. zeroBE and not
the headline edgeBE matrix run (Hannah's choice): edgeBE has its two end
domains solved against the CoastSat target, so at GIS 1 and 90 the model would
be partly fitted to one of the two things it is being compared against.
zeroBE carries no source/sink term in any domain, so all 90 are the model's
own response and neither target was fitted anywhere in it. The script checks
the run name for `zeroBE` before drawing; a caption saying "nothing was
fitted" over a solved run would be false.

**The model line:** net change in metres over the window = the run's own
endpoint rate (`change_rate_m_yr`) x 14 yr, the same conversion
`target_comparison.py` uses. It is smoothed at the same 7-domain window as the
two observed curves, because a raw model line against two smoothed targets
would make the gap between them partly an artefact of the treatment. One
asymmetry remains and is stated rather than hidden: the observed sides are
smoothed at transect resolution (~10 CoastSat and 5 dune transects per domain)
while the model exists only per domain, so its LOWESS runs over 90 points
rather than ~900. At a 3.5 km window the fitted curve barely notices the
difference in density, but it is not literally the same operation.

Output, `output/comparisons/target_comparison/smoothed_lowess7_with_cascade/`:
`lowess7_projected_vs_duneline_with_cascade_1996_2010_2024.png`,
`lowess7_total_change_vs_duneline_with_cascade_1996_2010_2024.png`,
`domain_values.csv`, `runs_used.csv`, `PROVENANCE.md`, `README.md`, `supporting/`.

<details><summary>Function notes (the original docstrings)</summary>

**`model_series()`**

```
The run's net change in metres, raw and smoothed, plus its provenance.

Smoothed at the same width as the observed curves; see the module
docstring for the resolution asymmetry that remains.
```

</details>

### compare_runs/target_comparison.py

Which observation should CASCADE be graded against? The two candidate
targets, the CoastSat shoreline and the digitized dune line, side by side
with the model, as net change in position over each model window, 1996-2010
and 2010-2024. Built 2026-09-19 (Hannah, by interview). The loaders are
`rate_windows.py`'s (imported), so the observations and runs are exactly the
ones `model_vs_observed` draws as rates; `smoothing_scale.py` and
`smoothed_lowess7_with_cascade.py` import this script in turn.

**Everything over the model period (14 yr per window):**

- *CoastSat target:* the window's LRR (`3-rates/coastsat/lrr/<window>`, the
  runner's scoring series) x 14 yr, the projected change;
- *dune-line target:* the measured dune-line net change
  (`3-rates/duneline/endpoint/<window>`) projected to the model years: its
  rate over the survey interval (11.6 yr for 1997-10 to 2009-05, 14.1 yr for
  2009-05 to 2023-07) x 14 yr;
- *model:* the run's own net change over its 14 years, unchanged (endpoint
  rate x 14 = last annual shoreline minus first).

Seaward positive, metres. The gap between the two targets is the beach-width
change they imply (CoastSat minus dune line): solid grey where the beach
widened, hatched where it narrowed.

**Raw and smoothed:** the lines are the raw domain means. `tables/skill.csv`
scores the model against both targets both raw and with the scoring target's
LOWESS treatment (raw means D1-10, 7-domain LOWESS beyond), the form the runs
are graded in.

**The CoastSat target's two modes** (2026-09-19, Hannah), named in the
vocabulary settled by interview on 2026-09-21: a rate turned into a distance
is named by the window it was *fitted* on, never by the arithmetic.

- `projected/`: the 1996-2024 LRR x 14 yr in both windows. The rate is carried
  onto windows it was not fitted on, so it is a projection. **The target in
  use**, paired with runs whose ends were solved against it.
- `total_change/`: each window's own LRR x 14 yr, as the runner grades. The
  rate is evaluated over the window it was fitted on, so nothing is
  extrapolated. Kept for the record.

These were `coastsat_full_period_lrr/` and `coastsat_subperiod_lrr/` until
2026-09-21; those names described the fit window but not what was done with
it. `--coastsat-target full` / `subperiod` still work as aliases. The fit
window is in every method string, so a figure pulled out of its folder still
says where the rate came from (Hannah, 2026-09-21); `total` has no single fit
window across the two panels, so it is filled in per window.

**Three model sets, one subfolder each** (Hannah has not chosen the target):

- `ends_solved_on_coastsat/`: the edgeBE run with its end domains solved
  against the CoastSat target. In projected mode this is the converged step of
  the end-domain solve against the 1996-2024 LRR, read from its
  `loop_log.csv`.
- `ends_solved_on_duneline/`: the 09-18 dune edge-solve run (mean3), end
  domains solved against the dune line.
- `ends_unsolved/` (2026-09-21, Hannah, for her advisor): the zeroBE arm of
  the same matrix cell, with no source/sink term in any domain, the two ends
  included, so all 90 domains are the model's own response and neither target
  was fitted. Also `unsolved_run_and_targets_<window>.png` (and `_smoothed`):
  that one run against both targets, the paired form with a single line. It is
  the only figure here whose GIS 1 and 90 mean anything.

Each holds `target_comparison_..._1996_2010_2024.png`: 1996-2010 above
2010-2024, both targets and the model on one y axis. The stem carries the
target mode and the model set (Hannah, 2026-09-21), because six folders used
to write the same basename.

`paired/` (2026-09-19, Hannah, style B of three rendered candidates): one
figure per window, `target_and_own_run_<window>.png`, each target with its own
run only: CoastSat with the CoastSat-solved run, the dune line with the
dune-solved run. The target is the house fill (blue seaward, red landward),
its run the black line, the misfit the gap between them. The end values each
run carries are in the panel titles; the runs have no other source/sink term
(GIS 2-89 are zero, checked from the index: `be_nonzero_domains == 2`).
`paired_smoothed/` (2026-09-19, Hannah): the fill is the target as graded
(raw domain means over GIS 1-10, the LOWESS beyond), with the raw domain means
as dots. The smoothed versions are drawn in projected mode only.

**Current end-domain solve:** `end-domain-boundaries/2026-09-29-ends-solved-on-lrr-1996-2024-split12`.
It was re-solved under option A on 2026-09-27 (Hannah: redo target_comparison
for the new wave climate and offset), on the adopted model on 2026-09-28
(Barrier3D `hatteras/adopted`, storms `v3_trim24`), after the dune-cap fix on
2026-09-29, and after the storms changed to `v3_split12_trim24` later that day.
The earlier solves, newest first:
`2026-09-29-ends-solved-on-lrr-1996-2024-dunecap`,
`2026-09-28-ends-solved-on-lrr-1996-2024-adopted`,
`2026-09-27-ends-solved-on-lrr-1996-2024-option-a`, and the ÷10 solve
`2026-09-19-end-domains-solved-on-lrr-1996-2024`.

**Fixed y range:** ±100 m on every figure here. It was 80 m from 2026-09-19 and
100 m from 2026-09-22 (Hannah), to match the metre figures in `3-rates` and
`4-comparisons` so the three trees can be laid side by side. Anything beyond
it is named in the caption and marked with a triangle at the edge. The legend
has two columns (2026-09-29): at three, the long target labels ran off the
page. Observed lines are pale and thick and model lines dark and thin
(Hannah, 2026-09-19: the same-hue solid/dashed pair was hard to tell apart).

**As a rate:** `--units rate` (2026-09-29, Hannah: "do the same net change
subfolder for target_comparison", choosing to keep the metres figures and add
m/yr beside them). Every metres column is divided by the 14 model years, which
undoes the x 14 exactly: the CoastSat LRR, the dune line's measured rate over
its own survey interval (no scaling to the window), and the model's endpoint
rate. Same figures, same skill table (bias/RMSE in m/yr), written to
`<projected|total_change>/change_rate/` with `_rate` after the target mode in
every stem; y axis ±8 m/yr, as model_vs_observed's rates.

Output, `output/comparisons/target_comparison/`: `README.md`,
`runs_used.csv`, `tables/domain_values_<w>.csv`, `tables/skill.csv`, and
`<model set>/target_comparison_1996_2010_2024.png` (PDF and CAPTIONS.md under
`supporting/`).

<details><summary>Function notes (the original docstrings)</summary>

**`_targets_line()`**

```
"CoastSat: ... · dune line: ... measured, scaled to 14 yr" (2026-09-22).

The header used to name only the CoastSat side, so the red line's dates
and its real interval -- which is not 14 yr -- were in the caption alone.
```

**`cs_label()`**

```
The CoastSat target named by the window its rate was FITTED on
(the 2026-09-21 vocabulary): PROJECTED when the 1996-2024 rate is carried
onto a 14-yr half, TOTAL when each window uses its own.
```

**`full_period_runs()`**

```
window -> (run_name, tag): the converged step of the 1996-2024 CoastSat
edge solve, from its loop_log.csv.
```

**`over_note()`**

```
' Beyond ±100 m, off the axis: ...' for the (frame, columns, label)
triples one figure draws, naming each out-of-range domain; '' if none.

The canvas marker is `mark_offaxis` at each draw site; this is the words
that go with it.
```

**`to_rate()`**

```
Every metres column over the model years: the x 14 undone. Column
names are kept so the drawing code reads either frame.
```

**`window_values()`**

```
Per domain, metres over the model period: both targets raw and LOWESS,
and each model set's net change.
```

**`_pad_title()`**

```
Lift the centred title clear of the fill bars.

draw_fills puts its bars at 1.025 in axes fractions with the year above
them, so a title at the default pad lands on "2022 fill". _title() has
already set the bold letter at the left; re-setting only the centred
string keeps it and moves both (pad is per-axes in matplotlib).
```

**`paired_figure()`**

```
Each target with the run solved on it, ONE FIGURE PER WINDOW, one panel
per target (Hannah, 2026-09-19, style B of three rendered candidates):
the target as the house fill (blue seaward, red landward), its run as the
black line, the misfit the gap between them. Both windows on one y axis.

smoothed=True (2026-09-19, Hannah): the fill is the target AS GRADED (raw
domain means over GIS 1-10, the rw.TARGET_WINDOW-domain LOWESS beyond, the form the runs
and the edge solve are scored against), the raw domain means as dots over
it; written to paired_smoothed/.
```

**`unsolved_figure()`**

```
The UNSOLVED run against both targets, one figure per window (Hannah,
2026-09-21, for her advisor). The paired form, except that there is only
one run: the same zeroBE line is drawn in both panels, because no part of
it was fitted to either target. (a) the CoastSat target as the fill, (b)
the dune-line target, the run in black over each.

With no end solve the two end domains are the model's own response too,
so this is the only figure here whose GIS 1 and GIS 90 mean anything.
```

**`load_model_sets()`**

```
key -> ({window: model frame}, index rows) for the three model sets.
Factored out of main 2026-09-21 so smoothing_scale.py reads exactly
the same runs; it depends on CS_MODE, which must be set first.
```

</details>

### compare_runs/smoothing_scale.py

The modelled net change in shoreline position against the change projected
from the CoastSat LRR, with the *rate* smoothed at several widths before it
is projected. Built 2026-09-21 (Hannah, by interview, after the same sweep on
the observations alone in `3-rates/coastsat/total_change/<window>/smoothed/`).

**The figure is the point** (Hannah, 2026-09-21):
`projected_vs_model_<window>.png` is one alongshore panel per model period,
all 90 domains, with the smoothing widths laid over each other on a
light-to-dark blue ramp (`SMOOTH_RAMP`) and the model in black. The black line
is the run untouched (nothing about the model responds to the smoothing), so
it is the one fixed thing in the panel, and the spread of the blue family
around it is the whole result.

The runs are graded against a target that is not the raw rate: raw domain
means over GIS 1-10, a LOWESS of the transect rates beyond
(`cascade_pipeline.coastsat_lowess`, via `hindcast.build_target_table`). The
darkest curve is that grading window; the palest is no smoothing at all. Over
GIS 1-10 all the curves coincide, because the splice keeps the raw domain
means there whatever the window. The widths are 0 (raw), 3, 5, 7 and 10
domains (0, 1.5, 2.5, 3.5, 5.0 km): the four the projected-LRR product uses,
plus 7, the grading window since 2026-09-28 (10 was the grading window until
then and stays in the sweep for comparison).

**The table behind the figure:** `tables/skill_by_window.csv` scores every
combination in two forms:

- *as_graded:* the target smoothed, the model raw. What the runner actually
  does, and what the figures draw; `skill.csv`'s `coastsat_lowess`
  generalised to every width.
- *scale_matched:* target and model both smoothed at the same width, the only
  form in which the two sides are treated alike. Kept for the record, not
  drawn.

**The null, and why it is not optional here:** interior r is 0.05-0.33. That
is the regime where a symmetric smoother inflates correlation hardest: it
strips high-frequency variance that the two sides do not share, so r rises and
RMSE falls at every window whether or not the model is any good. Without a
baseline the sweep draws a curve that looks like "the model improves at
coarser scales" and means nothing. So every r is reported beside
`r_null_p95`: the 95th percentile of r over `N_NULL` phase-randomised
surrogates of the same model series (same mean, same variance, same
alongshore autocorrelation, no relation to the target), each put through the
identical smoothing and splice. An r above that band is skill the smoother
cannot manufacture; an r inside it is not. For the as_graded form the raw
model is the same series at every window, so its null is redrawn per window
only because the target changed. The bias is close to smoothing-invariant and
needs no null; it is the one number in the table that a wider window cannot
flatter.

**One asymmetry, on the record:** the target's LOWESS is fitted at transect
resolution (~906 points) and then averaged to domains. The model exists only
at 90 domains, so smoothing it means LOWESS over those 90 values at the same
physical width (frac = window / n). Same width, coarser resolution. That is
why as_graded is reported too: it involves no model-side smoothing at all.

**Structure diagnostic** (`tables/target_structure.csv`): how much alongshore
structure each width removes from every LRR field on disk (1984-2004,
2004-2024, 1996-2024, 1996-2010, 2010-2024), reported as the SD removed in
m/yr, not as a share of variance. The domain-mean variance is dominated by the
long-wavelength swings, so a wiggle plainly visible on a figure reads as a few
per cent of it and the percentage badly undersells the effect (Hannah caught
this 2026-09-21, comparing against the 1984-2004 panels of
`input_prep/6-scr-smooth/lowess_method_comparison.py`, where the same LOWESS
removes roughly twice as much). The full reading of it is written into
`PROVENANCE.md` by the script.

**Scope:** the full-period CoastSat target (`target_comparison/projected/`,
the target in use, the 1996-2024 LRR x 14 yr in both windows) and the dune
line beside it, against all three model sets. `ends_unsolved` is the headline:
the zeroBE arm carries no source/sink term in any domain, so neither target
was fitted anywhere in it and all 90 domains are held out. The two solved sets
are swept too, which answers whether the edge solve still buys anything once
the grading is done at 5 km. Interior GIS 2-89 throughout, as the run index
scores. The script sets `target_comparison` to its canonical `projected` mode
(not an alias, since that module compares the mode by equality).

Output, `output/comparisons/target_comparison/smoothing_scale/`:
`projected_vs_model_<window>.png` (PDF and caption under `supporting/`),
`tables/skill_by_window.csv` (every width x model set x target x form: n,
bias, RMSE, r, r_null_p95), `tables/domain_values_<window>.csv`,
`tables/target_structure.csv`, `runs_used.csv`, `PROVENANCE.md`.

**Fixed 2026-09-30:** the caption string held mis-encoded characters ("GIS 1â10"); re-decoded to –, —, × and ±.

<details><summary>Function notes (the original docstrings)</summary>

**`_along()`**

```
Per-transect along-coast distance in metres, the convention every
target build here uses: each domain's transects spread evenly across its
500 m band, ordered within the domain by `order_col`.
```

**`duneline_transects()`**

```
(domain ids, along-coast m, rate m/yr) for the two-survey dune rate,
read and ordered exactly as rw.load_dune_endpoint_target does.
```

**`field_structure()`**

```
How much alongshore structure a LOWESS of each width takes out of one
LRR field, and whether there is independent error for it to average.

Reported as the SD removed IN m/yr, not as a share of variance: the
domain-mean variance is dominated by the long-wavelength swings, so a
wiggle that is plainly visible on the figure reads as a few per cent of it
and the percentage badly undersells the effect (Hannah caught this
2026-09-21, comparing against the 1984-2004 panels of
input_prep/6-scr-smooth/lowess_method_comparison.py).

Returns a dict, or None when the window's LRR has not been built.
```

**`target_structure()`**

```
field_structure over every LRR window on disk, so the grading window's
effect on the TARGET can be read against the other windows -- in
particular the 1984-2004 field the method-comparison figure draws, where
the same LOWESS removes roughly twice as much.
```

**`smooth_domain_series()`**

```
The model's analogue of the target's pass: lowess over the 90 per-domain
values at the same physical width, with the same GIS 1..skip splice. The
model has no transect resolution to smooth at -- see the asymmetry note in
the module docstring.
```

**`phase_randomise()`**

```
A surrogate with y's mean, variance and alongshore autocorrelation but
randomised phases, so it carries no relation to the target. The amplitude
spectrum is kept and only the phases are redrawn.
```

**`figure()`**

```
The alongshore picture (Hannah, 2026-09-21): the model's own net change
against the projected change from the LRR, ONE PANEL per model period with
every smoothing width laid over it on a light-to-dark ramp.

The model line is the run untouched, so it is the one fixed thing in the
panel: the spread of the blue family around it is the whole result.
```

</details>

### compare_runs/rate_windows.py

The hindcast windows, the model against an observation: the CoastSat
waterline and the digitised dune line, each with the run whose end domains
were solved on it drawn over it. edgeBE, full management, groin off.
`target_comparison.py`, `smoothing_scale.py` and
`smoothed_lowess7_with_cascade.py` import its loaders and constants.

**One folder, one script** (Hannah, 2026-09-17: "more organized and clear, and
potentially condensed, especially with the naming"). Until 2026-09-17 this was
two scripts and two output trees, `observed_vs_modeled_windows/` (CoastSat,
09-15) and `duneline_vs_modeled_windows/` (the dune line, 09-16), with two
vocabularies for one axis, five stems that all said "model", and the dune tree
carrying four copies of its own layout under `sensitivity/`. Now there is one
tree under `output/comparisons/model_vs_observed/` (named `rate_windows/`
until 2026-09-18), organised by observation and then reading, one naming rule
for every file, and each sensitivity drawn once.

**Why it is here and not in 5-scr:** the observed-only figure lives with the
observations (`data/hatteras_init/5-scr/3-rates/coastsat/lrr/`, drawn by
`coastsat_lrr_windows.py`). Once a run's curve is on the panel the figure
spans runs, and every cross-run figure is filed under `output/comparisons/`.
This script imports the observed drawing from the 5-scr producer (found
through `scr_paths`) rather than copying it, so the two cannot drift apart.
Those imports used to be two hand-built paths naming folders the 2026-09-22
reorganisation removed; a stale path does not fail where it is written, only
when the module is executed.

**The runs, since 2026-09-27:** the option A metres matrix
(Hs 2.0 / Tp 7.5 / asym 0.6 / high-angle 0.5, ends solved for it), 1996 and
2010 only (Hannah: "make the model vs observed figures for the new matrix
runs"). 1984-2004 and 2004-2024 have no metres run, so their panels say "model
not yet run for this window". The ÷10-offset runs drawn until then
(`HAT_<w>_edgeBE_road_bdm[_nourish]_nogroin`, 1984 on version-pair/v2) are in
`raw_runs/archive/2026-09-24-pre-metres/`, and their figures in
`output/archive/2026-09-27_model-vs-observed-div10/`. For the record, the
2026-09-15 set was: 1984-2004 the version-pair/v2 arm (topography v2),
1996-2010 calibration arm (offsets v1, the re-digitized 1997 line), 2004-2024
calibration arm with nourishment, and 2010-2024 calibration arm (run
2026-09-16 once the 2009 dune line gave a 2010 offset). Run folders are
resolved through `cascade_pipeline.run_registry` with the arm named
explicitly; the registry raises, listing the arms a run is under, rather than
guessing. The run index is keyed on (run_name, kind, tag) since 2026-09-16,
so the legacy arm name is translated.

**Three model sets, named for where their ends were solved.** The matrix runs
have their two end domains solved against CoastSat, so against the CoastSat
target the model meets the observation at GIS 1 and 90 by construction. A
dune-line target deserves the same (Hannah, 2026-09-17: "a fair comparison"),
so the dune-line figures draw runs whose ends were solved on the dune line:

- `coastsat`: the matrix;
- `dune-mean3`: the dune solve, mean of the end domain and its two inward
  neighbours (the main dune-line set);
- `dune-raw`: the dune solve, the end domain's own value (sensitivity).

The dune solve in use is
`end-domain-boundaries/2026-09-29-ends-solved-on-duneline-split12`. Its
history, newest first: re-solved after the storms changed to
`v3_split12_trim24` (2026-09-29, later); after the dune-cap fix
(`2026-09-29-ends-solved-on-duneline-dunecap`); on the adopted model
(`2026-09-28-ends-solved-on-duneline-adopted`: Barrier3D `hatteras/adopted`,
storms `v3_trim24`); under option A
(`2026-09-27-ends-solved-on-duneline-option-a`, pre-adoption); and the ÷10
solve on the re-digitized lines
(`2026-09-18-end-domains-solved-on-redigitized-duneline`, 1984-2004 carried
over from the 09-16 solve, whose lines did not change; each row of
`solved.csv` names its own run tag, so a carried-over row points back into the
09-16 experiment).

*Corrected 2026-09-30:* the source carried a note that the dune-solved sets
were ÷10 runs and not drawn. That stopped being true when the dune line was
re-solved under option A on 2026-09-27 (`DUNE_SOLVE_CURRENT = True`), and
every set is drawn. The captions named the 09-16 experiment until the same
date; they now name the solve in use.

**The observations:**

- *coastsat:* the CoastSat transect LRR, an OLS slope through ~250 satellite
  dates per transect. Two readings: `means`, the per-domain means as the 5-scr
  figure draws them (sign-coloured line, fill, ±1 std); and `lowess`, the
  scoring target as the fill (a 7-domain LOWESS, 10 until 2026-09-28, of the
  transect rates north of D10 and the raw means over D1-10, as
  `cascade_pipeline.hindcast.build_target_table` makes it for the runner),
  with the means as dots over it. The width follows the runner's
  `TARGET_WINDOW` (Hannah, 2026-09-28: "redo the target comparison with LOWESS
  7"); the CoastSat-solved end rates were re-solved against the LOWESS-7 value
  at GIS 90 the same day (`end-domain-boundaries/2026-09-28-ends-resolved-lowess7`).
- *duneline:* the digitised dune line, the feature CASCADE's shoreline
  actually is (a dune line behind a fixed berm). Vintages from
  `hat_topo_version.DUNE_LINE_FOR_YEAR` (1996 reads the 1997 line, 2010 the
  2009, 2024 the 2023); stations from
  `2-brie-offset/raw_offsets/<vintage>_duneline_offset_raw.csv`, read as the
  hindcast's end-year target loader reads them; seaward positive. Survey dates
  from `coastsat_vs_duneline.KNOWN_SURVEY_DATES`; 2023 has no known flight
  date and is centred on 2023-07-01, flagged in every caption that uses it.
  Read from the stored product `5-scr/3-rates/duneline/endpoint/<window>/`
  (2026-09-18), not computed here, so this figure and the stored numbers
  cannot disagree. Two readings: `endpoint` (two surveys differenced, per
  domain, over the interval) and `endpoint-lowess` (the same per transect,
  then the scoring target's treatment; LOWESS frac 0.111 vs CoastSat's
  0.110). A third reading, `lrr`, was retired 2026-09-18 with
  `3-rates/duneline_lrr` (Hannah: "these should not be lrr, they would just be
  endpoint, we are tracking net change").
- *both:* the two observations as lines on one panel, no fill, both as net
  change between the same two dune-line dates (2026-09-18, Hannah): CoastSat
  blue (`3-rates/coastsat/endpoint`, the mean position within ±6 months of
  each date, differenced) and the dune line red, each given the scoring
  target's LOWESS treatment; a model line per solve, all as the endpoint rate.
  The CoastSat LRR, the model's actual scoring target, is in `vs_shoreline/`.
  Drawn twice since 2026-09-29 (Hannah: "subfolder showing these plots as
  change rate and also net position change"): `both` in m/yr
  (`change_rate/`), and `both-netchange` in metres (`net_change/`), every line
  x the window's calendar span (14 or 20 yr). The model's is its own
  last-minus-first displacement exactly; the two observations are measured
  over the dune-line survey interval (11.6 / 14.1 yr) and scaled to the
  window, as `target_comparison.py` scales them; the interval stays in the
  caption.

**The model line:** `lrr_m_yr`, the OLS slope over the run's annual
shorelines (the estimator CoastSat and the run index use), or
`change_rate_m_yr`, the endpoint rate, last annual shoreline minus first (the
like-for-like estimator for a two-survey reading). Each reading is paired with
its own estimator; the one pairing that mixes them (endpoint observation, OLS
model) is kept as a sensitivity. The model line is black, not the site
config's model orange (Hannah, 2026-09-15): the observed line already carries
two hues and a fill, and a third hue on top read as noise. Black sits over
both fills and survives greyscale.

**Naming:** `model_vs_<feature>_<reading>_<start>_<end>.png`, and `_grid` for
the 2 x 2 by period (1984-start left, 1996-start right, earlier window above).
The feature is shoreline (CoastSat) or duneline; the reading means, smoothed or
netchange (2026-09-18, Hannah: the old `coastsat_` / `duneline_endpoint_` /
`both_` stems never said a model was being compared). The internal variant
keys (`coastsat/means`, ...) are unchanged; `OUTPUT_FOLDER` maps them to the
folders below. A sensitivity arm is also in the file stem (Hannah,
2026-09-21): the arm used to live only in the folder path, so the same
basename existed under the main level and each sensitivity and the copies
were indistinguishable once moved. Every figure folder keeps its PDFs and
CAPTIONS.md under `supporting/`. The legend has one entry per row: two abreast
ran past the page edge at 190 mm (2026-09-16).

Output, `output/comparisons/model_vs_observed/`:

```
vs_shoreline/domain_means/          model_vs_shoreline_means_<w>.png      ends solved on CoastSat
vs_shoreline/smoothed/              model_vs_shoreline_smoothed_<w>.png
vs_duneline/endpoint_net_change/    model_vs_duneline_netchange_<w>.png   ends solved on the dune line (mean3)
vs_duneline/net_change_smoothed/    model_vs_duneline_netchange_smoothed_<w>.png
vs_shoreline_and_duneline/          both targets, both solves
    change_rate/    model_vs_shoreline_and_duneline_rate_<w>.png       m/yr
    net_change/     model_vs_shoreline_and_duneline_netchange_<w>.png  metres
tables/             domain_rates_<w>.csv  every reading, every model set, the residual against each
                    skill.csv             bias and RMSE, GIS 2-89, per window x model set x estimator x target
runs_used.csv       one row per window per model set: run, arm, folder, timestamp, commit,
                    topography, offsets, and the dune line's vintages, dates and interval
y_bounds.txt        the shared y range and the rule behind it
sensitivity/
    ends-swapped/       vs_shoreline/* on the dune-solved runs, vs_duneline/* on the
                        CoastSat-solved runs: each target against the other solve
    dune-raw-solve/     vs_duneline/* on the raw-reading dune solve
    mixed-estimator/    model_ols_vs_duneline_netchange_<w>.png: the net-change
                        observation against the model's OLS rate
```

<details><summary>Function notes (the original docstrings)</summary>

**`_import_by_path()`**

```
A script in the input-prep tree, not a package; its drawing or its
readers are what is reused.
```

**`dune_solved_runs()`**

```
window -> (run_name, arm) from the dune solve's solved.csv, where the
arm is the experiment tag the registry translates.
```

**`load_model()`**

```
Both per-domain estimators for one run, and the provenance row for
runs_used.csv. (None, row) where no run exists.

preset names the source/sink preset the run was filed under; it defaults
to PRESET (edgeBE), the only one this module's own figures draw. It is a
parameter so a caller can read the zeroBE arm of the same matrix cell
(target_comparison's ends_unsolved set, 2026-09-21).
```

**`load_coastsat_target()`**

```
The per-domain scoring target: raw means D1-10, the TARGET_WINDOW LOWESS
beyond, from build_target_table on the window's transect_lrr_full.csv.
```

**`load_dune_endpoint()`**

```
Two surveys differenced per domain, seaward positive, from the stored
product 3-rates/duneline/endpoint/<window>/ (duneline_endpoint.py builds
it; 2026-09-18). Returns a frame (domain_number / mean_lrr / std_lrr, the
columns obs.draw_panel expects, holding the RATE in m/yr) and the
vintages, dates and interval it was built from.
```

**`load_coastsat_endpoint()`**

```
CoastSat NET CHANGE at the dune-line dates (3-rates/coastsat/endpoint,
the rate over the survey interval): (domain frame, target frame). The
target goes through the SAME builder as the LRR target, only reading the
endpoint rate column, so the two differ in the estimator alone.
```

**`load_dune_endpoint_target()`**

```
The two-survey rate per transect (from the stored product), then the
scoring target's treatment. Returns domain_number / target_lrr_m_yr /
source.
```

**`Observation()`**

```
Everything observed for one window: the CoastSat means and target,
the dune line in its three readings.
```

**`shared_bounds()`**

```
One half-range for every panel: the 5-scr rule (largest |rate| plus
1 m, rounded up) over every observed reading and both estimators of
every model set drawn.
```

**`shared_bounds_net()`**

```
The metres half-range for the net-change panels: the largest |rate x
span| over both endpoint observations and every model set's endpoint
rate, every window, plus NET_PAD_M, rounded up to NET_STEP_M.
```

**`_draw_observed()`**

```
The observation: the 5-scr panel as it draws itself (means, with the
std lines), a line with its fill (endpoint), or with a target frame the
target as the fill and the per-domain values as dots over it.
```

**`_draw_both()`**

```
Axes, bands and structures from the observed panel drawn empty, then
the two targets as lines, x scale (1 for rates, the span for net change).
The model lines go on afterwards.
```

**`_panel()`**

```
One window on one axes: the observation in the variant's reading and
the model line(s) over it. Returns whether any model line was drawn.
```

**`add_legend()`**

```
One entry per row: the model label carries the estimator and the
solve, and two abreast ran past the page edge at 190 mm (2026-09-16).
```

**`_stem_tag()`**

```
The sensitivity arm as a filename token, "" for the main level.

Hannah, 2026-09-21: the arm lived only in the folder path, so
`model_vs_shoreline_means_1996_2010.png` existed under the main level AND
under each sensitivity, and the three were indistinguishable once moved.
`sensitivity/ends-swapped` -> `ends-swapped`.
```

**`write_tables()`**

```
domain_rates_<w>.csv: every reading of the observation, both
estimators of every model set, the residual against each. skill.csv:
bias and RMSE over GIS 2-89 per window x model set x estimator x target.
```

**`fair_rows()`**

```
Each target scored against the runs solved on it, with its own
estimator: the three columns of the README's skill table.
```

</details>

### compare_runs/compare_runs.py

A general run-vs-run and run-vs-CoastSat tool: loads up to four finished
CASCADE runs from their saved shoreline change rate CSVs, overlays them on one
figure, and adds the LOWESS-smoothed CoastSat LRR for reference.
**`RUNS_TO_COMPARE` ships empty**: fill it (and `COMPARISON_NAME`) before
running.

Each run must already have been run by the hindcast runner, which writes
`<run_dir>/tables/shoreline_change_rate.csv` (before the 2026-09-10 layout
change, `<run_dir>/{run_name}_shoreline_change_rate.csv`; `run_layout`
resolves either). Columns: `gis_domain | change_rate_m_yr | lrr_m_yr | lrr_r2`.
The run folder is resolved from (run_name, period, preset, arm) by
`cascade_pipeline.run_registry.find_run_dir`, never joined by hand; it raises
naming the arms a run *is* under, which reads as "it is over there" rather than
"it never existed". The domain ids are read from the file's own `gis_domain`
column: this used to slice padded indices out of a 120-row CSV and pair them
positionally with a hand-built 1-90, so a run written in the other alongshore
order would have shifted every rate against its label with nothing raised.

**Which rate column is read is a method choice, not a spelling.** The file
carries `lrr_m_yr` (an OLS slope through the annual shoreline positions) and
`change_rate_m_yr` (the endpoint rate). This script reads `lrr_m_yr` because
the CoastSat target it is plotted against is an OLS rate too; putting an
endpoint rate up against an OLS one would move the difference between two
estimators into the residual panel, where it reads as model error.

**`RUNS_TO_COMPARE` fields:**

| field | meaning |
|---|---|
| `run_name` | folder name and filename prefix for this run; must match the name the run was made under |
| `period` | the run's period directory, e.g. `"1996_2010"` |
| `preset` | source/sink preset the run was made under (`calibBE`, `edgeBE`, `zeroBE`), the directory below the period |
| `arm` | optional forcing arm. Omit for the calibration arm. Naming no arm never silently searches: a run forced off the calibration wave climate must be asked for by name |
| `label` | legend label |
| `start_year` | which CoastSat period is drawn solid for this run, and which panel it appears in in the two-period figure |
| `sort_key` | optional; orders the run along the colour gradient (lower = lighter), typically the swept parameter. Without it, list order is used |
| `color` | optional explicit colour that overrides the gradient for this run only |
| `run_dir` | optional escape hatch: full path to a run outside the raw_runs tree (another drive, a collaborator's folder, an archive). Overrides the resolver for this run only. Prefer period/preset: a hand-typed path is how a figure ends up drawing a run other than the one it names |

Colours are sampled from `RUN_COLORMAP` (`YlOrRd`) at each run's rank along
`sort_key`, within `RUN_COLORMAP_RANGE` (0.35-0.95, avoiding near-white and
near-black), so reordering or adding runs never means re-picking hex codes.
`plasma` or `inferno` keep more contrast for more than five runs.

**The four runs the list used to hold are gone.** `HAT_1984_2004_SStest`,
`HAT_1984_2004_FinalSS`, `HAT_2004_2024_SStest` and `HAT_2004_2024_FinalSS`
were addressed by absolute paths into `output/raw_runs/source&sink_tests/`
and the flat `raw_runs/<name>` level; neither exists any more, and none of the
four names appears in `run_index.csv` or in `superseded_20260828/` (now
`output/archive/2026-08-28_full-tree/`). Their four figures sat at
`output/comparisons/source_sink_zones/` (2026-06-19) until 2026-09-17, when
they were deleted as unregenerable; nothing cited them (checked 2026-09-02 and
2026-09-17).

**CoastSat:** each `COASTSAT_DATASETS` entry points to a
`transect_lrr_full.csv` (one row per ~50 m transect). LOWESS is applied at
transect resolution, then averaged to domains. The period whose start matches
most runs is drawn solid (active) and the other faded (reference);
`ACTIVE_PERIOD_START` overrides that. Domains 1-10
(`LOWESS_SKIP_SOUTHERN_DOMAINS`) show raw transect scatter instead of LOWESS,
because Oregon Inlet boundary effects dominate there and smoothing hides the
sharp gradient; this matches the hindcast script. `LOWESS_WINDOW_DOMAINS`
lists one or two widths (7 since 2026-09-28; [7, 10] before), each with its
own style.

**Colour palette:** CoastSat is the cool blue family, light to dark = less to
more processed (transect scatter and domain means very light blue #9ECAE1,
7-domain LOWESS #6BAED6, 10-domain LOWESS #08519C). CASCADE runs are the warm
orange-red gradient. Blue against orange-red is the most colour-blind-safe
pairing (deuteranopia and protanopia); within each family lighter to darker
encodes less to more smoothing or a lower to higher parameter; the period is
carried by line style, not colour.

**Output:** `output/comparisons/<COMPARISON_NAME>/`:
`<name>_diagnostic.png` (quick multi-run check), `<name>_annotated.png` (the
publication figure with geographic annotations), the two-period figure
(1984-start runs left, 2004-start right, one y axis and one legend; skipped
for a period with no runs), and `<name>_residuals.png` (each run minus the
active CoastSat LOWESS). Per-run fit statistics are shown in the legend.

**History of fixes recorded in the code:**

- The 2004 CoastSat dataset pointed to `2004_2024_specific_dates`, a folder
  that did not exist, so the load failed silently and every 2004 run was
  compared against the 1984 CoastSat data in every figure (residuals reaching
  +9 m/yr instead of ±1-2, no transect scatter in the 2004 panel). The script
  now stops if a run's CoastSat period failed to load.
- The pier label heights were 85 and 70 (meant as percentages) in
  axes-fraction coordinates, which placed labels 85x the axes height above the
  plot; the tight bounding box grew to 300+ inches on save and crashed with an
  out-of-memory error. They are 0.85 and 0.70 now.
- Legends and captions placed below the canvas (negative y) made
  `bbox_inches="tight"` expand to include them, which on the wide two-period
  figure produced a ~490-inch canvas and another out-of-memory crash. Every
  figure now reserves its margins first and places legend and caption in
  figure-fraction coordinates on the canvas, saving with the figure's own
  bounds. The two-period figure's bottom margin is 0.20 (xlabel ~0.04, legend
  ~0.10, caption ~0.02); it was 0.30 before the legend was deduplicated.
  `loc="best"` used to put the diagnostic legend on top of the data.
- The two-period legend is deduplicated by visual style: both panels draw
  their CoastSat in the same colours, so one entry per style is enough and the
  period is dropped from the label (each panel's title names its period).
  Deduplicating by label text used to list every CoastSat entry twice.

**Dead code** (left in place; it runs but nothing uses it): `TRANSECT_DATASETS`
is not referenced anywhere (`COASTSAT_DATASETS` already carries the transect
data), `_gis_to_pad` is never called, and `sys`/`Path` are imported twice.

<details><summary>Function notes (the original docstrings)</summary>

**`assign_run_colors()`**

```
Auto-generate a light->dark color gradient for a list of run configs.

Runs are ranked by `sort_key` if every run provides one; otherwise by
their position in the list (so simply listing runs low-to-high parameter
value works without setting sort_key explicitly). Colors are sampled
evenly across RUN_COLORMAP_RANGE of RUN_COLORMAP, lightest first.

Any run with an explicit, non-None `color` field keeps that color
untouched and is excluded from the rank-based assignment - useful for
pinning one run (e.g. a baseline) to a fixed color while the rest of the
sweep auto-generates.

Parameters
----------
runs : list of dict - entries from RUNS_TO_COMPARE (or equivalent)

Returns
-------
list of dict - same entries, each with a resolved 'color' key (hex string)
```

**`load_run_rates()`**

```
Load the shoreline change rate CSV produced by HAT_hindcast_1984_2024_old version.py.

Parameters
----------
run_name : str       — used to resolve the default folder AND, inside it,
                        the rate CSV, whose name run_layout knows in both
                        the current and the pre-2026-09-10 layout. Which
                        is used is unchanged regardless of how the folder
                        was resolved.
period   : str, optional — "1984_2004" or "2004_2024".
preset   : str, optional — source/sink preset the run was made under.
arm      : str, optional — forcing arm; None means the calibration arm.
run_dir  : str, optional — ESCAPE HATCH. Full path to the folder holding
                        that CSV, for a run outside the raw_runs tree.
                        If given it OVERRIDES the resolver entirely.
                        Otherwise period and preset are both required and
                        find_run_dir locates the run, raising with what
                        IS on disk when it is absent.

Returns
-------
gis_ids   : int array — GIS domain IDs, READ FROM the CSV's gis_domain
                        column rather than assumed. Normally 1-90; a
                        short array means the run wrote fewer rows, which
                        is warned about rather than padded over.
rates_myr : float array, same shape — the run's LRR in m/yr, from
                        RUN_RATE_COL, aligned to gis_ids by construction.
run_dir   : str                                   — full path to run folder (resolved)
```

**`load_transect_data()`**

```
Load individual transect LRR values from transect_lrr_full.csv and derive
along-coast distance by spreading each domain's transects evenly across its
500 m band (mirrors 6-scr-smooth/lowess_method_comparison.py: load_transect_csv).

Returns
-------
domain_ids    : int array   — CASCADE domain ID for each transect
lrr_values    : float array — LRR (m/yr) for each transect
along_coast_m : float array — cumulative along-coast distance (m)
All three arrays share the same length (one entry per transect).
Returns (None, None, None) on load failure.
```

**`lowess_smooth_transect_to_domains()`**

```
Apply LOWESS at transect resolution using physical along-coast distance (m) as x,
then aggregate smoothed values to CASCADE domain resolution by averaging within
each domain.  Mirrors smooth_transect_df() + aggregate_to_domains() from
6-scr-smooth/lowess_method_comparison.py.

window_domains is converted to km (× DOMAIN_SPACING_M) so the physical window
is consistent regardless of transect density.

Returns
-------
gis_x    : int array   — domain IDs that have at least one transect
smoothed : float array — domain-averaged smoothed LRR (m/yr), same length as gis_x
frac     : float       — LOWESS frac used (for logging)
```

**`splice_lowess_with_raw_south()`**

```
Trim a LOWESS curve so it starts north of the southernmost `skip_n`
domains, leaving that southern zone to show raw transect scatter only.

Ported from HAT_hindcast_1984_2024_old version.py's function of the same name for
visual consistency between the two scripts - domains 1-skip_n are
boundary-affected (Oregon Inlet dynamics) and LOWESS smoothing there can
obscure the sharp gradient rather than clarify it, so the LOWESS line is
simply not drawn there; the raw scatter (already restricted to this same
zone via RAW_LRR_SOUTHERN_ONLY) carries the signal instead.

Parameters
----------
win_gis_x    : int array   - GIS domain IDs from the LOWESS result
win_smoothed : float array - LOWESS-smoothed LRR (m/yr)
skip_n       : int         - domains 1..skip_n are excluded from the
                              returned line. Defaults to
                              LOWESS_SKIP_SOUTHERN_DOMAINS.

Returns
-------
plot_x, plot_y : arrays - the LOWESS curve restricted to domains > skip_n
```

**`add_geographic_annotations()`**

```
Draw the standard Hatteras geographic annotation layer onto ax.
X-axis must be in GIS domain IDs (1–90).
```

**`plot_diagnostic()`**

```
Quick diagnostic plot — all runs + CoastSat + geographic annotations on
one panel. Uses the same drawing helper as plot_annotated/plot_two_period
so all comparison figures stay visually consistent.

Previously this plot had NO geographic annotations and used
loc="best" for the legend, which with 4+ runs landed on top of the
data (see uploaded screenshot) - both fixed below.
```

**`_draw_comparison_panel()`**

```
Draw geographic annotations + CoastSat curves + model run lines onto a
single axis. Shared by plot_annotated (one panel) and plot_two_period
(two panels, called once per period) so both stay visually identical.

Parameters
----------
show_reference_period : bool
    True  -> the inactive CoastSat period is still drawn, faded, with a
             "(ref)" label (useful in plot_annotated's single-panel view,
             where showing the other period for context is informative).
    False -> only the active period (cs["period_start"] == active_period)
             is drawn at all; the inactive period is skipped entirely.
             This is what plot_two_period uses, since each panel is
             already dedicated to one period - showing the other period
             there just duplicates what the OTHER panel is for.

Returns
-------
model_handles, cs_handles : lists of Line2D proxies for the legend
```

**`plot_two_period()`**

```
Two-panel figure: left = all runs with start_year=1984 (vs. 1984-2004
CoastSat), right = all runs with start_year=2004 (vs. 2004-2024 CoastSat).
Shares one y-axis range and one combined legend below both panels.

If every run shares the same start_year, the empty panel is skipped and
a single-panel figure is produced instead (so this is always safe to call
regardless of what's in RUNS_TO_COMPARE).
```

**`plot_residuals()`**

```
Residual panel — each model run minus the active CoastSat LOWESS curve.
Helps identify where each run over- or under-predicts relative to observations.
Only produced if PLOT_RESIDUALS = True.
```

</details>

### overwash/compare_overwash_observed.py and overwash/compare_overwash_figures.py

A run's overwash (Qow) against the observed record from the imagery. Both
read one run's `.npz` and the observation workbook (`Overwash_Matrix` sheet,
resolved by `site_layer/hat_overwash.py`) and write to
`output/comparisons/overwash/`. `compare_overwash_figures.py` is the fuller
of the two: the same two figures in a revised palette, plus the spatial
figures. The two are largely duplicates, and a future cleanup could keep one.

| figure | switch | shows |
|---|---|---|
| 1, stacked | `PLOT_STACKED` | model Qow heatmap above (every model year, imagery years highlighted), observed overwash below (binary, imagery years; hatched where there is no image); shared domain axis and island-section bar |
| 2, contingency | `PLOT_CONTINGENCY` | each image x domain classified hit / miss / false alarm / correct rejection; summary statistics printed |
| 3 and 4, spatial | `PLOT_SPATIAL` (figures script only) | overwash frequency alongshore, model (mean with IQR shading) against observed (bars): normalised on one panel, and as two panels |

**The run:** `HAT_1984_2004_calibBE_road_bdm_groin`, 1984-2004, calibBE,
resolved through `cascade_pipeline.run_registry`. *Fixed 2026-09-30:* the run
moved to `raw_runs/archive/2026-09-24-pre-metres/` when the matrix was
archived on 2026-09-24, and both scripts crashed on import until they were
pointed at it there (`RUN_KIND`, `RUN_TAG`). Whether to compare a current
matrix run instead is an open choice.

**Threshold:** a model cell counts as overwash if Qow > `QOW_THRESHOLD`
(dam³/yr). Start at 0 (any non-zero flux) and raise it if there are too many
false alarms.

**Imagery-year remap** (figures script): some imagery post-dates the storm it
captures, so `YEAR_REMAP` maps an image year to the model year that best
represents the event visible in it. The May 2004 imagery shows Hurricane
Isabel's (September 2003) deposits, so 2004 is remapped to 2003. 1996 imagery
is flagged as poor quality (`POOR_QUALITY_YEARS`).

**Palette** (figures script): muted coastal earth and ocean tones, with
consistent warmth across all elements; the section bar uses desaturated
versions of the warm and cool families in the contingency legend. Section
bar: warm linen #CABB9E for villages, dusty maritime blue #8DAFC2 between
them. Contingency: hit deep maritime green #2E7D5A, miss deep maritime blue
#2A5F8F, false alarm warm brick #A85840, correct rejection warm cream #F2EFE9,
no imagery warm greige #BEB9B4. Observed overwash uses the same brick as a
false alarm. The look was modelled on AGU / Nature Geoscience figures.

**To use on another run:** set `RUN_NAME`, `RUN_PERIOD`, `RUN_PRESET` (and
`RUN_KIND`/`RUN_TAG` if it is not in the matrix), then adjust `START_YEAR`,
`END_YEAR`, the domain constants and `SECTIONS` to match it, and choose
`QOW_THRESHOLD`. If the run is not where you say, `find_run_dir` raises naming
where it *is* on disk.

**History:** the paths used to be absolute literals under a folder spelling
that no longer exists (`input_preperation`), so the scripts could not read
their observations and crashed on the first save; they are anchored on the
repo root now. The workbook left `scripts/input_prep/8-overwash-analysis/`
on 2026-09-10 and is resolved by `site_layer/hat_overwash.py` since
2026-09-18. The `.npz` loader substitutes stand-in classes for anything the
pickled model references but cannot import, and is duplicated from
`plot_overwash.py` so each script stands alone.

<details><summary>Function notes (the original docstrings)</summary>

**`plot_contingency()`**

```
For each imagery year × domain where obs data exists, classify as:
  Hit (H)             : Qow > T  AND  obs = 1   (green)
  Miss (M)            : Qow <= T AND  obs = 1   (blue)
  False Alarm (FA)    : Qow > T  AND  obs = 0   (orange)
  Correct Rejection(C): Qow <= T AND  obs = 0   (light grey)
  Not Assessed (NA)   : obs = NaN               (medium grey)
Rows: imagery years only. Columns: all 90 domains.
```

</details>

<details><summary>Function notes (the original docstrings)</summary>

**`compute_spatial_data()`**

```
Returns per-domain summary arrays used by both spatial profile figures.

qow_mean  : mean annual Qow per domain across all model years (dam³/yr)
qow_p25   : 25th percentile of annual Qow per domain
qow_p75   : 75th percentile of annual Qow per domain
obs_freq  : fraction of assessed imagery years with overwash per domain
            (NaN where no assessed imagery exists for that domain)
bar_colors: colour per domain based on section type (village vs inter-village)
```

</details>

### smoothing_vs_cascade/smoothing_vs_cascade.py

Applies LOWESS smoothing to the CoastSat shoreline change rates of both
periods (1984-2004, 2004-2024) and draws a CASCADE run against them. LOWESS is
applied to each period's series independently; it keeps the large-scale
alongshore pattern and removes per-domain noise.

| figure | shows |
|---|---|
| `overview_smoothed.png` | both periods, raw CoastSat (faded) with the LOWESS overlay |
| `smoothed_only_comparison.png` | both periods, LOWESS lines only (for presentations) |
| `combined_periods.png` | both periods' smoothed lines on one panel |
| `smoothing_sensitivity_<period>.png` | one panel per LOWESS window |
| `window_comparison.png` | raw CoastSat with every window in `COMPARE_WINDOWS_DOMAINS` overlaid, for choosing a window |
| `cascade_vs_lowess.png` | the CASCADE run against raw and LOWESS CoastSat, per period |
| `cascade_by_window.png` | the CASCADE run against one LOWESS window per panel |

Also `coastsat_smoothed_table.csv`. Output:
`output/comparisons/smoothing_vs_cascade/1984_2004/`.

**Bandwidth:** `LOWESS_FRAC` is the fraction of the data in each local fit:
0.10 is ~9 domains (more local), 0.167 ~15 domains (the default), 0.20 ~18
domains (smoother, loses finer structure). Window sizes in domains convert to
frac as n / 90, over the whole island regardless of how many domains a CSV
has data for: 5 domains ≈ 0.056 (2.5 km), 10 ≈ 0.111 (5 km), 15 ≈ 0.167
(7.5 km).

**`CASCADE_RUNS`:** one `dict(label, period, csv)` per run, where `period`
matches a CoastSat period label ("1984–2004") and `csv` is
`run_rate_csv("<run_name>")`, which finds the rate CSV in either run-folder
layout (`base=` for a run outside `output/raw_runs`). Add one colour to
`C_CASCADE` per run. An empty list skips the model figures.

**Cannot draw the model today** (open, needs your choice):

- `CASCADE_RUNS` names `HAT_1984_2004_SQ_BE_Hs2p0`, which no longer exists
  anywhere under `raw_runs/`, so the model is skipped and only the CoastSat
  figures are drawn.
- Even with a current run, `load_cascade_rate` expects columns
  `gis_domain_id` and `model_rate_m_per_yr`, but run tables have carried
  `gis_domain`, `lrr_m_yr` and `change_rate_m_yr` since the 2026-09-10 layout
  change. Which rate to read (OLS or endpoint) is a method choice, so this was
  not changed blind.

**Fixed 2026-09-30:** the script drew every figure and then crashed printing
its summary on a Windows console (arrows that cp1252 cannot encode); it now
switches its output to UTF-8 first.

**Styling note:** after `apply_style()` it sets its own rcParams (Arial,
dotted grid, no top/right spines) and its own period colours, so its figures
are not fully in the house style. The smoothing windows run cool to warm, teal
(5 domains) to amber (10) to crimson (15), so more smoothing reads as a warmer
colour, and `cascade_by_window.png` uses the same colour per window as
`window_comparison.png` so the figures can be read against each other. Legends
use `loc="best"` because the free corner differs between periods.

**CSV quirk:** some exports prepend the file name to the first column
("domain_lrr_summary.csvdomain_number"); the loader strips it.

**History:** every path used to be an absolute literal; the output one had
lost its drive and wrote figures to `C:\scripts\`, and the input ones still
spelled the folder `input_preperation`. They are anchored on the repo root
now, and the CoastSat CSVs are resolved through `hat_observed_rates.py`
(2026-09-18).

<details><summary>Function notes (the original docstrings)</summary>

**`run_rate_csv()`**

```
The shoreline change rate CSV inside one run folder.

RESOLVED, NOT JOINED. That file is tables/shoreline_change_rate.csv in the
new run layout and {run_name}_shoreline_change_rate.csv in the old one;
run_layout.resolve returns whichever is on disk, so a half-migrated tree
reads either way. `base` is the folder the run folder sits in -- pass it
for a run outside output/raw_runs.
```

**`load_cascade_rate()`**

```
Load a CASCADE shoreline change rate CSV produced by HAT_hindcast_1984_2024_old version.py.

Expects columns:
    gis_domain_id       – integer 1–90 for real domains, NaN for buffer rows
    model_rate_m_per_yr – shoreline change rate in m/yr

Returns a DataFrame with columns [domain, model_rate] filtered to
DOMAIN_MIN–DOMAIN_MAX, or None if the file is missing.
```

**`domains_to_frac()`**

```
Convert a window size in CASCADE domains to a LOWESS frac value.

Always divides by the total island domain range (DOMAIN_MAX - DOMAIN_MIN + 1)
so that fracs are consistent regardless of how many domains have valid data
in a given CSV.  This ensures:
  5  domains → frac ≈ 0.056
  10 domains → frac ≈ 0.111
  15 domains → frac ≈ 0.167  (matches the LOWESS_FRAC default)
```

**`add_annotations()`**

```
Add all geographic reference annotations to an axis.

Layer order (bottom → top):
  1. Wimble Shoals influence zone  (hatched amber fill, bottom label)
  2. Community shaded spans        (steel-blue fill, top labels)
  3. Village center lines          (dashed gray,  y=0.88)
  4. Pier lines                    (dash-dot blue, y=0.76, rotated)
  5. Groin lines                   (dotted red,    y=0.76, rotated)

All label y-positions use blended axes-fraction coordinates so they
stay fixed relative to the panel height regardless of data range.
```

**`annotation_legend_handles()`**

```
Return proxy artists explaining the annotation layer types.
Append these to a plot's legend handle list so readers can decode
all reference marks without scanning every label individually.
```

**`plot_smoothed_only()`**

```
2-panel: LOWESS smoothed lines only, no raw data.
Cleanest version for presentations or dissertation figures.
```

**`plot_combined_periods()`**

```
Single panel: both periods of smoothed CoastSat on one axis.
Good for directly comparing the two periods.
```

**`plot_smoothing_sensitivity()`**

```
3-panel showing the effect of each LOWESS window on CoastSat data.
One smoothed line per panel so individual window behavior is clear.
Fracs are derived from COMPARE_WINDOWS_DOMAINS / 90 domains, matching
exactly what plot_window_comparison uses.
```

**`plot_window_comparison()`**

```
2-panel figure (one per period) showing the raw CoastSat LRR (faded)
overlaid by LOWESS-smoothed lines for each window size in
COMPARE_WINDOWS_DOMAINS.  All smoothed lines share one panel so spatial
patterns and differences between window choices are directly visible.

Window sizes are specified in CASCADE domains and converted to LOWESS fracs
using the actual number of valid domains in each dataset.

Parameters
----------
cs_1984, cs_2004 : pd.DataFrame or None
    Loaded CoastSat data for each period.
out_path : str
    Full path for the comparison PNG.
window_domains : list of int
    Number of CASCADE domains for each smoothing window to compare.
    Defaults to COMPARE_WINDOWS_DOMAINS from CONFIG.
```

**`plot_cascade_vs_lowess()`**

```
One panel per CoastSat period (1984–2004, 2004–2024).

Each panel shows:
  • Raw CoastSat LRR (faded, period color)
  • Three LOWESS-smoothed CoastSat curves (green / orange / purple)
  • CASCADE modeled change rate(s) for that period (thick black line)

CASCADE runs are matched to panels by their 'period' key in CASCADE_RUNS.
If no CASCADE run exists for a period, that panel still shows the CoastSat
smoothed curves alone.

Parameters
----------
cs_1984, cs_2004 : pd.DataFrame or None
    CoastSat data for each period.
cascade_runs : list of dict
    Loaded CASCADE rate DataFrames, each with an extra 'label' and
    'period' key (same structure as CASCADE_RUNS but with 'df' added).
out_path : str
    Full path for the comparison PNG.
window_domains : list of int
    LOWESS window sizes in CASCADE domains (from COMPARE_WINDOWS_DOMAINS).
```

**`plot_cascade_by_window()`**

```
3-panel figure (one per LOWESS window) for a single CoastSat period.

Each panel shows:
  • Raw CoastSat LRR (faded, period color)
  • ONE LOWESS-smoothed CoastSat curve (bold, period color)
  • CASCADE modeled rate (thick black)

This isolates the model-vs-observation comparison for each smoothing
choice so you can judge fit quality without the visual clutter of
seeing all three windows simultaneously.

Panels share the same y-axis limits so differences in smoothing
level — not axis scaling — drive the visual comparison.

Parameters
----------
cs_df : pd.DataFrame or None
    CoastSat data for the target period.
cascade_runs : list of dict
    Loaded CASCADE run dicts (with 'df', 'label', 'period' keys).
    Only runs whose 'period' matches period_label are plotted.
period_label : str
    e.g. "1984–2004".  Used for title and CASCADE run matching.
out_path : str
    Full path for the comparison PNG.
window_domains : list of int
    LOWESS window sizes in CASCADE domains.
```

</details>

