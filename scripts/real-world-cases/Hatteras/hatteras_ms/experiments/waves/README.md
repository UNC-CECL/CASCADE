# Wave climate

Kept on 2026-10-23: HAT_wave_grid.py. The other scripts this README describes were retired (see ../RETIRED_STUDIES.md; recover from commit a7c4358c).

Which wave climate (Hs, asymmetry, high-angle fraction) fits each window, and the figures that recommend one.

Each study drives the unchanged runner (`../../HAT_hindcast.py`); see `../README.md` for how
the experiments are run and filed, and how they import each other across folders.

| Script | What it asks |
|---|---|
| `HAT_wave_grid.py` | Which wave climate fits each window best, on the current dunes, storms, buffers and groin? |
| `HAT_wave_grid_smoothed_score.py` | Four-parameter wave grid, scored on the smoothed model. |
| `HAT_wave_grid_smoothed_score_plot.py` | Figures for the smoothed-score wave grid: the best settings per window and shared. |
| `HAT_wave_grid_fixed_ends.py` | The four-parameter wave grid with the end domains held fixed. |
| `HAT_wave_shortlist_ends_solved.py` | Do different wave settings win once the end domains are solved per setting? |
| `HAT_wave_recommendation_figures.py` | Figures behind the recommended hindcast wave climate. |
| `HAT_metres_2_wave_sensitivity.py` | Wave-climate sensitivity with the offset in metres, natural scenario, both windows. |
| `HAT_metres_2_wave_sensitivity_plot.py` | The figures for the natural-scenario wave sensitivity. |

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### HAT_wave_grid_smoothed_score.py

Four-parameter wave grid, scored on the smoothed model.

From the script's original header:

```text
Four-parameter wave grid, scored on the smoothed model (2026-09-25).

Designed with Hannah on 2026-09-25, after the 2026-09-24 step-2 study
(one parameter at a time plus two 2-D grids) left no search over all four
wave parameters together, and none at all under full management:

    score     share of the alongshore variation explained, 1 - SSE/SST, with
              the MODEL SMOOTHED LIKE THE COASTSAT TARGET (LOWESS over 10
              domains, the southern 10 raw: common.smooth_like_target),
              interior GIS 2-89, against the window's CoastSat LRR target.
              Bias, RMSE, correlation and the raw (unsmoothed) score beside it.
    coarse    Hs 0.75 1 1.5 2  x  Tp 7 8 10  x  asym 0.5 0.7 0.9
              x  high-angle 0.3 0.45 0.55  = 108 per period x scenario
              (Hs 0.65 and Tp 12 left out: they drowned the barrier)
    scope     natural and full management, 1996-2010 first, then 2010-2024
    refine    per period x scenario, a 3x3x3x3 grid at half the coarse step
              around the best coarse-or-refine setting (the midpoints to the
              neighbouring coarse values), launched automatically
    cross     the top 5 per period x scenario that lack a run in the other
              window are run there, so the shared pick has candidates
    shared    one setting for both windows: the lowest mean of smoothed RMSE
              / that window's flat-line RMSE, among settings run in both

Every run here is made on the Barrier3D route_overwash fix (checked at
start). The step-2 runs are not reused as grid cells (most predate the fix);
they are rescored on the smoothed output into tables/step2_rescored_smoothed.csv
for comparison.

WHERE: output/raw_runs/experiments/wave-climate/2026-09-25-wave-grid-smoothed-score/
    README.md, tables/, figures/, logs/<phase>_<scenario>/<period>/<settings>.log
    runs/<phase>_<scenario>/<period>/zeroBE/<run_name>/   (on disk only)
    phase = coarse, refine, cross

USAGE (from the project root):
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/waves/HAT_wave_grid_smoothed_score.py run all
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/waves/HAT_wave_grid_smoothed_score.py run coarse --periods 1996
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/waves/HAT_wave_grid_smoothed_score.py score
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/waves/HAT_wave_grid_smoothed_score.py rescore-step2
```

<details><summary>Function notes (the original docstrings)</summary>

**`refine_values()`**

```text
The best value and the midpoints to its coarse neighbours (half the
coarse step), clipped at the ends of the coarse range.
```

**`drop_already_run()`**

```text
Skip a cell whose settings already have a result (any phase) in that
period and scenario: a refine grid overlaps the coarse one at its centre.
```

</details>

### HAT_wave_grid_smoothed_score_plot.py

Figures for the smoothed-score wave grid: the best settings per window and shared.

From the script's original header:

```text
Figures for wave-climate/2026-09-25-wave-grid-smoothed-score: the best settings on the
smoothed score, for each window and shared, one figure per scenario.

Two sources, drawn by the same code:
    grid    this study's runs (tables/all_runs.csv), once the sweep has run
    step2   the 2026-09-24 step-2 runs rescored on the smoothed output
            (tables/step2_rescored_smoothed.csv): drawn first, 2026-09-25,
            while the grid was still running (Hannah asked for the figures)

Writes, under output/raw_runs/experiments/wave-climate/2026-09-25-wave-grid-smoothed-score/figures/:
    best/<source>/per_period/best_by_period_<scenario>_<source>.png
    best/<source>/shared/best_shared_<scenario>_<source>.png
each with supporting/ (PDF, CAPTIONS.md, the data CSV).

Reads tables only; it never rebuilds the run index, so it is safe to run
while the sweep is going.

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/waves/HAT_wave_grid_smoothed_score_plot.py [--source step2|grid|both]
```

Notes that were in the code:

```text
Which score ranks the runs (2026-09-27, Hannah: pick on the UNSMOOTHED model,
as the runner's own RMSE and every matrix run are scored). "raw" draws the
per-domain model as a dark line with no smoothed curve; figures go to
figures/best/<source>_raw/.
```

<details><summary>Function notes (the original docstrings)</summary>

**`picks()`**

```text
Best per window x scenario (smoothed share explained), and the shared
setting per scenario (lowest mean smoothed RMSE / flat line, run in both).
```

**`fig_top()`**

```text
The top TOP_N runs by smoothed score in each window x scenario, drawn as
the smoothed profiles that were scored (added 2026-09-25: Hannah asked why
the best figures looked the same as the raw-scored ones).
```

</details>

### HAT_wave_grid_fixed_ends.py

The four-parameter wave grid with the end domains held fixed.

From the script's original header:

```text
The four-parameter wave grid with the end domains FIXED (2026-09-27).

Hannah, 2026-09-26/27: sweep the wave climate with the end source/sink
terms held fixed, the ends first re-solved under the metres offset. The
fixed pair per window comes from
end-domain-boundaries/2026-09-27-ends-resolved-metres-offset/tables/ends.json
(solved at the step-2 baseline waves, full management, against CoastSat) and
is imposed in every run -- natural and managed alike -- through the edgeBE
preset and HAT_BE_OVERRIDE.

Everything else is HAT_wave_grid_smoothed_score's, unchanged (the same coarse
grid, refine, cross-runs and smoothed score), so this study and the zeroBE
grid (wave-climate/2026-09-25-wave-grid-smoothed-score) compare one to one.
This module only points that driver at its own folder, preset and ends.

WHERE: output/raw_runs/experiments/wave-climate/2026-09-27-wave-grid-fixed-ends/

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/waves/HAT_wave_grid_fixed_ends.py run all --jobs 8
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/waves/HAT_wave_grid_fixed_ends.py score
```

<details><summary>Function notes (the original docstrings)</summary>

**`build_targeted()`**

```text
The targeted list for one window (2026-09-27): the top n on the raw score
in that window from this study's runs on the superseded ends (read from
archive_tag if already archived, else from tables/all_runs.csv) and from the
zeroBE grid, the other window's top 5 here, and the adopted candidates.
```

**`run_targeted()`**

```text
The 2010-2024 rerun, TARGETED (2026-09-27, Hannah: "is there not a faster
version"): after the 2010 ends were re-solved at Hs 2 (GIS 1 +137.6 -> +18.8),
only the settings that can change a decision are re-run, as phase "target":
the top 15 on the raw score in 2010-2024 from the superseded fixed-ends runs
(archive/2026-09-27-fixed-ends-2010-ends-solved-at-hs1) and from the zeroBE
grid, the 1996-2010 top 5 (for the one-setting and one-parameter picks),
and the adopted candidates. The list is tables/targeted_2010_cells.csv.
Coarse cells already run on the new ends are kept.
```

</details>

### HAT_wave_shortlist_ends_solved.py

Do different wave settings win once the end domains are solved per setting?

From the script's original header:

```text
Wave shortlist with the end domains solved per setting (2026-09-26).

Hannah, 2026-09-26: "what about when you solve for the ends, are there
different wave parameters that perform the best?" The 2026-09-25 wave grid
ran zeroBE (nothing imposed at GIS 1 or 90), and the stored edgeBE values
(HATTERAS_BE_EDGE_ONLY) were solved at Hs 2.5 on the old /10 offset, so they
do not carry over. What the ends must carry depends on the waves, so each
wave setting gets its own solve. Chosen with Hannah:

    shortlist  the top 10 zeroBE settings (smoothed score) per window x
               scenario from wave-climate/2026-09-25-wave-grid-smoothed-score
    scenarios  natural and full management, both windows (40 chains)
    ends       solved against each window's CoastSat LRR, as the matrix end
               values were: GIS 1 against the raw domain mean, GIS 90 against
               the LOWESS-10 value (the target table's own splice)
    solve      step 0 is the setting's zeroBE grid run (ends 0, 0). Each end
               is stepped on its own (they are 89 domains apart): step 1 from
               the 2026-09-11 response (about 0.09 m/yr of residual per m/yr
               imposed at GIS 1, 0.13 at GIS 90), then the secant through the
               last two probes. Lockstep: every chain runs its next probe
               before any is solved again. Converged at |residual| <= 0.02
               m/yr at both ends, or stop after MAX_STEPS probes
    score      the converged (or last) run, as the grid: share of the
               alongshore variation explained by the model smoothed like the
               CoastSat target, interior GIS 2-89; bias, r, the raw score, the
               imposed ends and their residuals beside it
    fixed      metres offset (dune line v1), edgeBE preset with HAT_BE_OVERRIDE,
               no groin, no relocations, Barrier3D route_overwash fix

WHERE: output/raw_runs/experiments/wave-climate/2026-09-26-wave-shortlist-ends-solved/
    README.md, tables/{shortlist,solve_log,all_runs}.csv, figures/,
    logs/<scenario>/step<k>/<period>_<settings>.log
    runs/<scenario>_step<k>/<period>/edgeBE/<run_name>/     (on disk only)

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/waves/HAT_wave_shortlist_ends_solved.py run [--jobs 8]
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/waves/HAT_wave_shortlist_ends_solved.py score
```

Notes that were in the code:

```text
The first pass (plain secant, 5 probes) converged 6 of 40 chains. GIS 90
was fine; GIS 1 in 2010-2024 is not smooth -- imposed 0-25 m/yr gives
residuals near zero, anything above ~25 gives +5 to +15 whatever the value
-- so the secant took steps to -466 and +335 m/yr, and four probes drowned
the barrier. The resume keeps every probe already run and steps each end on
its own:
bracketed  (a probe on each side of the target): interpolate between the
closest pair, held at least 10% inside it so it always shrinks
otherwise  secant through the two latest probes, capped at +-STEP_CAP
converged  an end within TOL keeps its value
A drowned probe carries no residual and is left out of the history.
```

<details><summary>Function notes (the original docstrings)</summary>

**`find_run()`**

```text
The run a probe made: its folder under the step's tag, matched on the
wave settings in its metadata (the run name leaves defaults out).
```

**`history()`**

```text
[(ends, residuals)] for one chain: its zeroBE run, then every probe
that produced a run (from solve_log.csv).
```

</details>

### HAT_wave_recommendation_figures.py

Figures behind the recommended hindcast wave climate.

From the script's original header:

```text
Figures behind the recommended hindcast wave climate (2026-09-27).

Hannah, 2026-09-27: "given all of the different tests ... what do you suggest
as the best wave parameters to use for the model hindcast. Provide figures and
a clear and detailed explanation". Draws, from the studies already run (no new
runs), into output/raw_runs/experiments/wave-climate/2026-09-27-wave-recommendation/figures/:

  1_score_by_parameter.png      best raw score reachable at each value of each
                                parameter (the complete zeroBE grid, coarse +
                                refine), both windows, both scenarios
  2_agreement_across_tests.png  the best 1996-2010 setting of each search, by
                                parameter
  3_recommended_vs_coastsat.png the recommended setting against CoastSat,
                                rate and position change, both windows
  4_high_angle_roughness.png    domain-to-domain roughness and score against
                                the high-angle fraction
  5_window_2010.png             the recommended setting scored on 2010-2024
                                and on 2010-2020 (before the 2021 CoastSat step)

Scores are the RAW share of the alongshore variation explained (the per-domain
model against the CoastSat LOWESS-10 target, interior GIS 2-89), as the runner
and the matrix runs are scored.
```

Notes that were in the code:

```text
coarse phase only: the one full factorial. Refine values (Hs 1.25, Tp 7.5,
...) were run next to a few settings only, so their "best" is not
comparable (it drew a false dip at Tp 7.5).
```

```text
Target, observed change and the header scores at common.SMOOTH_DOMAINS
(7 since 2026-09-28); the table's scores were made at 10, so re-scored here.
```

```text
The step-2 high-angle sweep: the one clean series through 0.5 (0.1-0.55,
Hs 1.0, Tp 8, asym 0.8, ends zero), so every point differs in the
high-angle fraction only.
```

<details><summary>Function notes (the original docstrings)</summary>

**`observed_change_smoothed()`**

```text
Observed end-minus-start change, smoothed as the target is
(common.smooth_like_target: LOWESS at common.SMOOTH_DOMAINS, the southern
10 raw). The 5-scr table carries 0/3/5/10 only, so 7 is built from its raw
(window 0) column (Hannah, 2026-09-28: the group's range is 7).
```

**`draw_shoals()`**

```text
Shoal zones as faint hatched boxes, as coastsat_lrr_windows.draw_shoals
(5-scr/3-rates) draws them, named at the bottom when `label`.
```

**`draw_fills()`**

```text
A bar over each enabled model-input fill in the window, the year on it,
just inside the top of the panel (the title sits above the frame).
```

**`fig3()`**

```text
`runs` {(period, scenario): run folder} and `ends` {period: (GIS 1, GIS 90)}
draw another end solve (2026-09-28: the position-change solve); by default
the fixed-ends sweep's option-A runs and ENDS.
```

**`fig6()`**

```text
2010-2024: the same setting as 1996-2010 against the one allowed change,
Hs 2.0 -> 2.5, each on the ends solved for it (2026-09-27, Hannah).
```

</details>

### HAT_metres_2_wave_sensitivity.py

Wave-climate sensitivity with the offset in metres, natural scenario, both windows.

From the script's original header:

```text
HAT_metres_2_wave_sensitivity.py -- wave-climate sensitivity, natural scenario, offset in metres
THE QUESTION (Hannah, 2026-09-24, by interview). With the island offset in
metres, how does the model's alongshore shoreline change respond to each wave
parameter when nothing human acts on the island, in both canonical windows?

THE DESIGN (every choice Hannah's)
    scenario     natural: no road management, no beach/dune management, no
                 fills, no relocations; no groin; zeroBE; offset in metres
                 (dune line, CURRENT v1)
    periods      1996-2010 and 2010-2024, each against its own CoastSat LRR
    baseline     Hs 1.0 m, Tp 8 s, asymmetry 0.8, high-angle fraction 0.45
                 (re-centring on asymmetry 0.7 pending, 2026-09-25)
                 (the best metres point of experiments/2026-09-24-island-
                 offset-scale-wave-tuning, found under full management)
    stage 1      one parameter at a time around the baseline:
                   wave_height  Hs    0.65 0.75 1.0 1.25 1.5 2.0 2.5 3.0
                   high_angle   ahf   0.1 0.2 0.3 0.4 0.45 0.5 0.55
                   asymmetry    asym  0.3 0.4 0.5 0.6 0.7 0.8 0.9
                   wave_period  Tp    6 7 8 10 12
                 plus the baseline under full_management, per period
    stage 2      a 5 x 5 grid per period over the two parameters whose range
                 moves the share of alongshore variation explained the most
                 (mean over both periods; a drowned run counts as the worst
                 score of its period). Each axis: the 5 stage-1 values centred
                 on that parameter's best value, shifted inward at the ends.
                 The other two stay at the baseline. Cells already run are
                 reused, not re-run.

LAYOUT (output/raw_runs/experiments/wave-climate/2026-09-24-metres-2-wave-sensitivity/)
    README.md, tables/, figures/, logs/<group>/<period>/<settings>.log
    runs/<group>/<period>/zeroBE/<run_name>/     group = baseline, wave_height,
                                            high_angle, asymmetry, wave_period,
                                            baseline_full_management,
                                            grid_<p1>_x_<p2>
    The run folders stay on disk only: this scenario's run name is ~85
    characters and sits twice in each path, past Windows' 260, so git cannot
    index them (the same choice as the offset-scale study).

USAGE
    python HAT_metres_2_wave_sensitivity.py run stage1 [--jobs 6] [--dry-run]
    python HAT_metres_2_wave_sensitivity.py run stage1 --scenario full_management
    python HAT_metres_2_wave_sensitivity.py score
    python HAT_metres_2_wave_sensitivity.py score-window --start 2010 --end 2020
    python HAT_metres_2_wave_sensitivity.py run stage2 [--jobs 6] [--dry-run]
    python HAT_metres_2_wave_sensitivity.py run combo --hs 1 --tp 8 --asym 0.7 --ahf 0.4
    python HAT_metres_2_wave_sensitivity.py run grid --pair wave_height wave_period \
        --values1 1 1.25 1.5 2 2.5 --values2 6 7 8 10 12 --periods 1996
```

Notes that were in the code:

```text
Shared with the offset-scale study: the target, the alongshore scores, the
drowning reader, the number spelling and the keep-awake.
```

```text
PENDING (2026-09-25): Hannah chose asymmetry 0.7 for the Buxton dip (GIS 6-7)
and asked to re-centre the one-at-a-time sweeps on it; the re-run was
stopped to test a combination first (runs/combos/). Until it runs, the
baseline stays 0.8: every figure and the scoring select on it, and 0.7-
centred sweeps do not exist yet. To re-centre: set 0.7 here, then
`run stage1` and `run stage1 --scenario full_management`.
```

```text
Stage 2 (the two grids) was chosen and run around the first baseline,
asymmetry 0.8, and stays there: selecting or drawing it on BASELINE would
find no runs.
```

```text
THE FULL-MANAGEMENT SWEEP (Hannah, 2026-09-24, after stage 1 showed
management halves the 2010-2024 bias): the same stage-1 values under
full_management, filed as full_management_<parameter>/ beside the natural
folders. Its baseline is baseline_full_management/, already run.
```

```text
A clean finish, or a drowned barrier, is a result: re-running either
would only hit the runner's existing-folder guard and overwrite the
log with that refusal (found when the sweep was paused 2026-09-24).
...and a run that finished and only failed to replace the shared
run_index.csv (a Windows lock between parallel runs, 2026-09-24):
its outputs are complete and scored.
```

```text
THE SELECTION RULE, and why it changed (2026-09-24). The first rule took
each parameter's range of variation explained over BOTH periods and counted
a run without a score as the worst score of its period. It chose Hs x Tp
for the wrong reasons: 2010-2024 natural runs explain -800% to -2800% at
every setting (a ~-4.5 m/yr bias), so its ranges measure how badly a
setting fails; the drowned Tp 12 became each period's worst score; and two
crashes were scored as worst although the rule named drownings only. That
selection is kept in tables/stage2_selection_first_rule.csv and its 12
finished runs in grid_wave_height_x_wave_period/. Hannah then chose: the
range over runs that SURVIVED in 1996-2010 only, best value from the same
runs; the grid still covers both periods.
```

```text
The crash_check* logs are the NUMBA_BOUNDSCHECK / JIT-off diagnostics
of 2026-09-24, not sweep cells: scored as cells they appeared as
unscored duplicates of the baseline (and put a failure mark on it).
```

```text
which Barrier3D: recorded by the runner since 2026-09-24; a run
without the field predates it and ran on the unfixed model
```

```text
Added 2026-09-24 (Hannah): score the 2010-2024 runs on 2010-2020 too. The
CoastSat 2010-2024 target is lifted from ~0 to +1.07 m/yr by the island-wide
+17 m step into 2021, which no model run can make. The model's rate over
the window is the same OLS estimator the runner uses (shoreline.compute_lrr)
on the first (end - start + 1) annual states; the observed one is
5-scr/3-rates/coastsat/lrr/<start>_<end>/, LOWESS at the runner's window (7 since 2026-09-28).
```

```text
One hand-picked setting of all four (added 2026-09-25, Hannah: "test a
combo" before re-running stage 1). Filed as runs/combos/ and
runs/full_management_combos/; scored with everything else.
```

<details><summary>Function notes (the original docstrings)</summary>

**`stop_reason()`**

```text
Why a run has no score: drowned, crashed, or the error it raised.

A run that dies with no Python traceback and no drowning is the silent
access violation in Barrier3D's jitted route_overwash (an out-of-bounds
read once a domain has prograded; memory note of 2026-09-11, where it was
found with the groin). Seen here with no groin, in 2010-2024 natural runs.
```

**`grid_cells()`**

```text
A named grid: every combination of two parameters' values, the other
two at the baseline, in the given periods; cells already run anywhere in
the study (natural scenario, same settings) are reused, not re-run.

Added 2026-09-24 to finish the Hs x Tp grid for 1996-2010 (Hannah), which
the first stage-2 rule chose and which was stopped when the rule was
replaced; its best cell (Hs 1.25, Tp 10) turned out to be the best
natural 1996-2010 run of the study.
```

</details>

### HAT_metres_2_wave_sensitivity_plot.py

The figures for the natural-scenario wave sensitivity.

From the script's original header:

```text
HAT_metres_2_wave_sensitivity_plot.py -- the figures for the natural-scenario wave sensitivity
Reads tables/all_runs.csv and tables/stage2_selection.csv (written by
HAT_metres_2_wave_sensitivity.py), each scored run's tables/shoreline_change_rate.csv,
the CoastSat LRR target and the observed CoastSat position change. Writes,
under output/raw_runs/experiments/wave-climate/2026-09-24-metres-2-wave-sensitivity/figures/:

  stage1/scores_by_<parameter>_both_periods.png
      share of alongshore variation explained, bias and RMSE against the
      parameter, both periods on one axis (cross-period consistency)
  alongshore/<period>/{natural,full_management}/rate_and_position_change_by_<parameter>[_full_management]_<period>.png
      the modelled rate and position change along the island for each value
      against CoastSat
  management/natural_vs_full_management_baseline_both_periods.png
      the baseline under the natural scenario and under full management
  stage2/grid_<p1>_x_<p2>_both_periods.png
  best/best_settings_by_period.png, best/best_settings_both_periods.png
      the best settings found for each window, and one setting for both
      the grid as lines: the score against the first parameter, one line per
      value of the second, a panel per period
```

Notes that were in the code:

```text
The two grids: the one the stage-2 rule chose (both windows), and the Hs x
Tp grid the first rule chose, stopped, then finished for 1996-2010 only.
```

```text
Split and redrawn 2026-09-25 (Hannah): one figure per scenario and per rule
(best/per_period/, best/shared/), drawn for the screen, a colour per
scenario, and the model smoothed like the target as a faint dashed line.
```

```text
Hannah's combo test (run combo --hs 1 --tp 8 --asym 0.7 --ahf 0.4) closed a
2x2 whose other corners were already run: the old baseline, and each of the
two single changes that give the Buxton dip its observed depth.
Only the two Hannah asked to see (2026-09-25): the old baseline and her
combination. The single-change corners are in tables/combo_asym0.7_ahf0.4_2x2.csv
and the explorer's "Asym x high-angle" set.
```

<details><summary>Function notes (the original docstrings)</summary>

**`no_score_marker()`**

```text
A drowned barrier is a model result; a crash is not (Barrier3D's
route_overwash access violation), so they get different marks.
```

**`best_both_periods()`**

```text
{(scenario, period): row} for the one setting, run in both periods, with
the lowest mean of RMSE / that window's flat-line RMSE. Share explained is
not averaged: 2010-2024's values run to -2800% and would decide alone.
```

</details>
