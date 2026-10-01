# hindcast_vs_observed - the hindcast against the observations

The hindcast (edgeBE and zeroBE, full management, groin off, each window run
with its end-solved boundaries) drawn against what was observed: the CoastSat
shoreline and the dune line, as rates and as net change. Nothing here runs the
model.

| script | question | writes to `output/comparisons/` |
|---|---|---|
| `rate_windows.py` | model vs each observation, every window | `model_vs_observed/` |
| `target_comparison.py` | which observation should the model be graded against? | `target_comparison/` |
| `smoothing_scale.py` | does the LOWESS width used for grading matter? | `target_comparison/smoothing_scale/` |
| `smoothed_lowess7_with_cascade.py` | the smoothed halves-overlay sheets with the model drawn over | `target_comparison/smoothed_lowess7_with_cascade/` |

`rate_windows.py` holds the shared loaders; the other three import it (and
`target_comparison.py`), which is why the four live together. Run
`rate_windows.py` first only if you also want its figures; the imports do not
need it to have run.

The script name `rate_windows.py` predates the output folder's rename to
`model_vs_observed/` (2026-09-18) and was kept.

## The scripts in detail

Moved here from `scripts/analyze_output/README.md` on 2026-09-30, when the
scripts were grouped by question.

## rate_windows.py

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

## target_comparison.py

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

## smoothing_scale.py

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

## smoothed_lowess7_with_cascade.py

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
