# target_comparison — CoastSat or the dune line as the CASCADE target?

> **Lost?** [`FIGURES.md`](../../../FIGURES.md) is the one-page index of which figure answers which question.

The two candidate targets and the hindcast, as **net change in position (m)
over each 14-yr model window**, 1996–2010 and 2010–2024. Built 2026-09-19
(Hannah, by interview) by `scripts/analyze_output/compare_runs/HAT_target_comparison.py`,
which reuses `HAT_rate_windows.py`'s loaders, so the observations and runs are
the ones `model_vs_observed/` draws as rates.

| line | what | over 14 yr how |
|---|---|---|
| blue | CoastSat target — in `projected/`, the 1996–2024 LRR (a rate carried onto windows it was not fitted on); in `total_change/`, the window's own LRR (the runner's scoring series) | LRR x 14 |
| red | dune-line target: the MEASURED net change between the digitized lines (1997-10→2009-05, 11.6 yr; 2009-05→2023-07, 14.1 yr, end date assumed) | its rate x 14, scaled to the model years |
| black | CASCADE: the run's own net change over the window, unchanged | endpoint rate x 14 |

The space between the two targets is the beach-width change they imply
(solid grey widened, hatched narrowed). Lines are raw domain means;
`tables/skill.csv` also scores against the LOESS-smoothed targets, the form
the runs are graded in.

**Two versions of the CoastSat target (2026-09-19, Hannah; renamed
2026-09-21 by interview).** A rate turned into a distance is named by the
window it was **fitted** on, never by the arithmetic — so the folders say what
was done to the rate, not just which window it came from:

| was | now | why |
|---|---|---|
| `coastsat_full_period_lrr/` | `projected/` | the 1996–2024 rate is carried onto two 14-yr windows it was **not** fitted on |
| `coastsat_subperiod_lrr/` | `total_change/` | each window's own rate over its own years — nothing extrapolated |

The observations-only twin of `projected/` is
`data/hatteras_init/5-scr/3-rates/coastsat/projected/`, and of `total_change/`
is `.../coastsat/total_change/`.

- `projected/` — **the target in use**: the 1996–2024 LRR x 14 yr
  in BOTH windows, paired with runs whose ends were solved against it
  (`raw_runs/experiments/2026-09-19-edgesolve-lrr1996_2024/`: GIS 1 / 90 =
  +28.5 / +24.5 in 1996–2010, +37.1 / +25.4 in 2010–2024).
  `projected/paired_smoothed/` is the same pairing with both
  targets drawn as GRADED (raw domain means over GIS 1–10, 10-domain LOESS
  beyond) as the fill and the raw domain means as dots (Hannah, 2026-09-19).
- `total_change/` — kept for the record: each window's own LRR (what
  the runner grades against), paired with the matrix runs.

**`smoothed_loess7_with_cascade/` (2026-09-22, Hannah, by interview).** The
two-panel smoothed sheets from
`5-scr/4-comparisons/shoreline_vs_duneline/smoothed_loess7/` with the zeroBE
run over them in dark green — 1996–2010 above 2010–2024, all three curves
LOESS-smoothed at 7 domains. The model-vs-dune-line bias is the same in both
sheets (+6.6 m in 1996–2010, −14.5 m in 2010–2024); what changes is the
distance to the SHORELINE target, −14.8 m against the long-term projection in
2010–2024 but −28.6 m against that window's own rate. See its README.

**`smoothing_scale/` (2026-09-21, Hannah, by interview).** Does the grading
window matter? The 10-domain LOESS the target is built with had never been
examined, and the `coastsat` → `coastsat_loess` rows in `skill.csv` show it is
worth 3–4 m of RMSE. `HAT_smoothing_scale.py` sweeps it over raw / 1.5 / 2.5 /
5.0 km against all three model sets, in two forms (target smoothed and model
raw, as the runner grades; and both smoothed), each r carried beside the 95th
percentile of 1000 phase-randomised surrogates of the same model series.

**The answer is that it does not, and that r was never the number to read.**
RMSE falls at every window for every model set, but so does the null: over all
96 rows, no r clears its own null band, and the gap *widens* with the window.
The apparent gain from smoothing is the smoother, not the model. What the
sweep leaves standing is the **bias**, which barely moves with the window
(unsolved run vs CoastSat: −11.0 → −10.6 m in 1996–2010, −14.8 → −14.4 m in
2010–2024) and is the one number a wider window cannot flatter. See
`smoothing_scale/PROVENANCE.md` for the tables and the reading rule.

The alongshore companion to this, on the observations alone, is
`data/hatteras_init/5-scr/3-rates/coastsat/total_change/1996_2024/smoothed/`.

**Every figure title carries quantity, window and method** (Hannah,
2026-09-21), and the method names the source and the window the rate was
*fitted* on, so a figure pulled out of its folder still says where its target
came from:

> Projected shoreline change, 1996–2010 (CoastSat LRR **1996–2024** × 14 yr)
> Total shoreline change, 1996–2010 (CoastSat LRR **1996–2010** × 14 yr)

One year apart in the parenthetical, and that is the entire difference between
the two targets.

**Every figure here uses one y axis, ±80 m** (Hannah, 2026-09-19); each caption names the few domain values that run off it (Cape Point in the sub-period 2010–2024 figures, and the dune line at GIS 12 in 1996–2010).

The dune-line target is the same in both (sub-period: 1997→2009 and
2009→2023, each scaled to 14 yr). Each version has the layout below.
`python ... HAT_target_comparison.py --coastsat-target projected|total`
(`full` and `subperiod` still work, as the pre-2026-09-21 names).

```
ends_solved_on_coastsat/   the matrix edgeBE run (end domains solved on CoastSat)
ends_solved_on_duneline/   the 09-18 dune edge-solve run, mean3 (solved on the dune line)
ends_unsolved/             the zeroBE arm of the same matrix cell: NO source/sink term in
                           ANY domain, the two ends included (2026-09-21, Hannah, for her
                           advisor). Nothing was fitted to either target, so all 90 domains
                           are the model's own response and the ends are readable.
    target_comparison_1996_2010_2024.png   1996–2010 above 2010–2024; PDF, CAPTIONS under supporting/
    (ends_unsolved only) unsolved_run_and_targets_<window>.png and _smoothed:
                           the paired form with ONE line — the same run in both panels,
                           (a) against the CoastSat target, (b) against the dune line
paired/                    START HERE: each target with ITS OWN run only (CoastSat with the
                           CoastSat-solved run, the dune line with the dune-solved run), both
                           one figure per window, all on the same ±80 m axis; per window
                           (a) the CoastSat target as the house fill with its run in black,
                           (b) the dune-line target the same way; end values in the figure title.
                           Chosen 2026-09-19 over dots+lines and residuals (style B)
    target_and_own_run_1996_2010.png
    target_and_own_run_2010_2024.png
tables/
    domain_values_<w>.csv  per domain: both targets raw and LOESS (m), the dune line's
                           measured change and interval, both model sets (m)
    skill.csv              model minus target, GIS 2–89: bias, RMSE (m), r; per window,
                           model set, target (coastsat, coastsat_loess, duneline, duneline_loess)
runs_used.csv              the runs, with their index provenance
```

Hannah had not chosen the target when this was built; the two subfolders are
the two answers side by side.

**No other source/sink correction.** The two solved sets are edgeBE: the only
source/sink terms are the two end domains, GIS 1 and GIS 90 (checked from the
run index, `be_nonzero_domains == 2`; the paired figure refuses to draw
otherwise). GIS 2–89 carry none, so the interior of each model line is the
model's own response. `ends_unsolved/` is zeroBE — `be_nonzero_domains == 0`,
checked the same way — so there GIS 1 and GIS 90 are the model's response too.

**`ends_unsolved/` (2026-09-21, Hannah, for her advisor).** The same matrix
cell as `ends_solved_on_coastsat/` (full management, groin off, `1984-start`
v2 / `2004-start` v1, offset v1), run with no boundary term at all:
`matrix/1996_2010/zeroBE/HAT_1996_2010_zeroBE_road_bdm_nogroin` and
`matrix/2010_2024/zeroBE/HAT_2010_2024_zeroBE_road_bdm_nourish_nogroin`.
Nothing in it was fitted to either candidate target, so the CoastSat target —
the 1996–2024 LRR, the same curve in both windows — and the dune line are both
held out. Interior GIS 2–89, model minus target:

| window | vs CoastSat | vs the dune line |
|---|---|---|
| 1996–2010 | bias −11.0 m, RMSE 22.7 m, r 0.24 | bias +10.3 m, RMSE 27.7 m, r 0.25 |
| 2010–2024 | bias −14.8 m, RMSE 26.8 m, r 0.10 | bias −14.3 m, RMSE 27.6 m, r 0.05 |

Removing the end solve costs little in the interior: against CoastSat the RMSE
goes 21.4 → 22.7 m in 1996–2010 and 24.9 → 26.8 m in 2010–2024. In 1996–2010
the unsolved run sits between the two solved runs against BOTH targets, so the
edge solve there is doing little the interior notices.

| window | ends solved on CoastSat (GIS 1 / 90, m/yr) | ends solved on the dune line |
|---|---|---|
| 1996–2010 | +32.2 / +10.0 | −16.0 / −3.3 |
| 2010–2024 | +72.6 / +31.3 | +20.1 / +19.2 |
