# target_comparison — CoastSat or the dune line as the CASCADE target?

> **Redrawn 2026-09-27 on option A** (Hannah: "Now that we have the new wave climate/offset, dont we need to redo th analyses in here").
>
> **What every model set now is:**
> - metres island offset, Hs 2.0 / Tp 7.5 / asymmetry 0.6 / high-angle 0.5;
> - the CoastSat-solved and unsolved sets are the option A matrix runs;
> - the ends were re-solved under option A against the dune line (`raw_runs/experiments/end-domain-boundaries/2026-09-27-ends-solved-on-duneline-option-a/`) and against the 1996–2024 LRR (`.../2026-09-27-ends-solved-on-lrr-1996-2024-option-a/`).
>
> **Where the old version is:** the /10-era tree is in `output/archive/2026-09-27_target-comparison-div10/`. Every number below comes from the redraw.
>
> **What changed from the /10 version:**
> - In 1996–2010 the model's correlation with the CoastSat target now clears its null at every smoothing width, for every model set (`smoothing_scale/`). In 2010–2024 it still clears it nowhere.
> - The dune-line solve still barely moves the interior (RMSE within 1 m of the CoastSat-solved run).

> **Lost?** [`FIGURES.md`](../../../FIGURES.md) is the one-page index of which figure answers which question.

The two candidate targets and the hindcast, as **net change in position (m)
over each 14-yr model window**, 1996–2010 and 2010–2024. Built 2026-09-19
(Hannah, by interview) by `scripts/analyze_output/compare_runs/target_comparison.py`,
which reuses `rate_windows.py`'s loaders, so the observations and runs are
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
  (`raw_runs/experiments/end-domain-boundaries/2026-09-27-ends-solved-on-lrr-1996-2024-option-a/`: GIS 1 / 90 =
  +4.5 / +27.5 in 1996–2010, +4.8 / +20.5 in 2010–2024; /10 era +28.5 / +24.5 and +37.1 / +25.4).
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
sheets (+10.7 m in 1996–2010, −10.7 m in 2010–2024; /10 era +6.6 and −14.5);
what changes is the distance to the SHORELINE target, −11.0 m against the
long-term projection in 2010–2024 but −24.7 m against that window's own rate
(/10 era −14.8 and −28.6). Numbers in its `PROVENANCE.md`.

**`smoothing_scale/` (2026-09-21, Hannah, by interview).** Does the grading
window matter? The 10-domain LOESS the target is built with had never been
examined, and the `coastsat` → `coastsat_loess` rows in `skill.csv` show it is
worth 3–4 m of RMSE. `smoothing_scale.py` sweeps it over raw / 1.5 / 2.5 /
5.0 km against all three model sets, in two forms (target smoothed and model
raw, as the runner grades; and both smoothed), each r carried beside the 95th
percentile of 1000 phase-randomised surrogates of the same model series.

**The answer is that the window does not matter.** RMSE falls at every window
for every model set, but so does the null band. What the sweep leaves standing is
the **bias**, which barely moves with the window (unsolved run vs CoastSat:
−7.2 → −6.9 m in 1996–2010, −10.7 → −10.3 m in 2010–2024) and is the one
number a wider window cannot flatter.

**Option A changes the r reading for 1996–2010** (2026-09-27):
- Every 1996–2010 row now clears its null, at every width, for every model set: r 0.35–0.49 against a null of about 0.18–0.29.
- 2010–2024 still clears nowhere.
- In the /10 version no r cleared its null in any of the 96 rows, which is why "r was never the number to read". That held for the /10 runs, not for option A in 1996–2010.

See `smoothing_scale/PROVENANCE.md` for the tables and the reading rule.

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

**Every figure here uses one y axis**, ±80 m until 2026-09-27 and ±100 m since the option A redraw (Hannah, 2026-09-19: one axis). Each caption names the few domain values that run off it; Cape Point in the 2010–2024 total-change figures is one of them.

The dune-line target is the same in both (sub-period: 1997→2009 and
2009→2023, each scaled to 14 yr). Each version has the layout below.
`python ... target_comparison.py --coastsat-target projected|total`
(`full` and `subperiod` still work, as the pre-2026-09-21 names).

```
ends_solved_on_coastsat/   total_change/: the option A matrix edgeBE run; projected/: the
                           option A 1996-2024 LRR solve (end domains solved on CoastSat)
ends_solved_on_duneline/   the option A dune-line solve, mean3 (2026-09-27)
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
cell as `ends_solved_on_coastsat/` (full management, groin off, option A),
run with no boundary term at all:
`matrix/1996_2010/zeroBE/HAT_1996_2010_zeroBE_offsetmetres_road_bdm_nogroin` and
`matrix/2010_2024/zeroBE/HAT_2010_2024_zeroBE_offsetmetres_road_bdm_nourish_nogroin`.
Nothing in it was fitted to either candidate target, so the CoastSat target —
the 1996–2024 LRR, the same curve in both windows — and the dune line are both
held out. Interior GIS 2–89, model minus target (projected/, raw domain means):

| window | vs CoastSat | vs the dune line |
|---|---|---|
| 1996–2010 | bias −7.2 m, RMSE 19.6 m, r 0.35 | bias +14.0 m, RMSE 27.7 m, r 0.39 |
| 2010–2024 | bias −10.7 m, RMSE 24.1 m, r 0.18 | bias −10.2 m, RMSE 23.9 m, r 0.22 |

/10 era: 1996–2010 −11.0 / 22.7 / 0.24 and +10.3 / 27.7 / 0.25; 2010–2024 −14.8 /
26.8 / 0.10 and −14.3 / 27.6 / 0.05.

Removing the end solve still costs little in the interior. Against CoastSat, the RMSE goes from 18.3 m (ends solved on CoastSat) to 19.6 m in 1996–2010, and from 23.5 m to 24.1 m in 2010–2024.

In 2010–2024, all three model sets sit the same distance from both targets: bias about −10 m and RMSE about 24 m against either one.

End rates (m/yr, GIS 1 / GIS 90):

| window | ends solved on the 1996–2024 LRR (projected/) | ends solved on the window's own LRR (total_change/, the option A matrix) | ends solved on the dune line (mean3) |
|---|---|---|---|
| 1996–2010 | +4.5 / +27.5 | +4.84 / +17.55 | −3.0 / +7.6 |
| 2010–2024 | +4.8 / +20.5 | +18.8 / +24.54 | +3.4 / +15.0 |

/10 era: CoastSat (window's own) +32.2 / +10.0 and +72.6 / +31.3; dune line −16.0 / −3.3 and +20.1 / +19.2.

## How this relates to the other study (added 2026-09-28)

`output/comparisons/target_comparison/` and `output/raw_runs/experiments/island-offset/2026-09-28-offset-source-duneline-vs-shoreline-option-a/` ask different questions of the same option A model.

- **target_comparison** holds the model fixed (dune-line offset) and changes the ruler: which observation should CASCADE be graded against, CoastSat or the dune line? Every run is graded against both targets. There are three sets of end rates (solved on CoastSat, solved on the dune line, unsolved). Scores are bias, RMSE and r in metres over 14 yr.
- **the offset-source study** changes the model's starting island: is the BRIE offset built from the dune line or from the CoastSat shoreline? Every run uses unsolved (zeroBE) ends, and each run is graded on its own feature: the dune-line offset against the dune line, the shoreline offset against CoastSat. Scores are the share of variation explained, then bias and r. Observations are smoothed at 7 domains and the model is left unsmoothed.

| | target_comparison | offset-source study |
|---|---|---|
| varies | target, and which target the ends were solved on | offset source (model input) |
| offset | dune line in every run | dune line vs shoreline (1996/ and 2010/shoreline/v1) |
| ends | CoastSat-solved, dune-solved, unsolved | unsolved (zeroBE) only |
| graded against | both targets | each run's own feature |

**Where they overlap.** The offset study's dune-line, full-management run is effectively the same model as target_comparison's `ends_unsolved` set (option A, zeroBE, dune-line offset). Graded against the dune line in metres, the two agree up to one scaling choice:

| model minus dune line | offset-source study | target_comparison (smoothed_loess7_with_cascade) |
|---|---|---|
| 1996–2010 | +9.6 m | +10.7 m |
| 2010–2024 | −10.5 m | −10.7 m |

- **Why 1996–2010 differs.** The 1997→2009 dune lines span 11.6 yr. target_comparison scales that change to 14 yr; the offset study compares it unscaled. 11.6 → 14 is a factor of 1.2, which is the 1 m gap. 2009→2023 spans 14.1 yr, so 2010–2024 nearly agrees.
- **The shoreline gradings are not directly comparable.** The offset study compares the model's LRR × 14 with CoastSat LRR × 14; target_comparison uses the model's end-minus-start change.

**Consistent findings:**
- **1996–2010:** the model sits roughly on CoastSat but about 10 m seaward of the dune line. The dune line retreats island-wide while the model holds, and that retreat may be an imagery artefact.
- **2010–2024:** the model erodes past CoastSat (about −24 m against the window's own rate, the 2021 step).
- **Only the offset study adds:** a shoreline-built offset fits slightly better at the domain scale (+3–4 points raw, level when smoothed). So the offset source is a small lever and not the cause of the 2010–2024 failure.
