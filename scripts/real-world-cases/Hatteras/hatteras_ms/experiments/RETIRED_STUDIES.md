# Retired experiment studies

Removed on 2026-10-23 when the hindcast moved to CASCADE-obx. Each asked a question that is now
answered and built into the model; the runs, figures and NOTE.md write-ups stay on D: under
`output/raw_runs/experiments/<theme>/`. The scripts were written for the old model, and several no
longer ran (the /10 `asrun` offset mode and `run_index.csv` are gone).

To recover a script, it is in commit `a7c4358c` under `scripts/real-world-cases/Hatteras/hatteras_ms/experiments/`
(before the theme folders):

    git show a7c4358c:scripts/real-world-cases/Hatteras/hatteras_ms/experiments/<script>.py

The shared helpers the matrix analysis and figures used from `HAT_metres_1_offset_units.py`
(`coastsat_target`, `interior`, `smooth_like_target`) moved to `site_layer/hat_scoring.py`. The
relocation results write-up moved to `hatteras_ms/docs/RELOCATION_COMPARISON_RESULTS.md`.

## The retired scripts

Also retired the same day: `analyze_output/compare_runs/adoption_2026-09-28/` (`adoption_before_after.py`, `adoption_before_after_figures.py`, its README), the one-off comparison of the matrix before and after the 2026-09-28 model adoption. Its runs were archived and their .npz files deleted; its results stay in `output/comparisons/adoption_2026-09-28/`. Recover with `git show a7c4358c:scripts/real-world-cases/Hatteras/analyze_output/compare_runs/adoption_2026-09-28/<file>`.

**code_checks/**: `HAT_adopt_dune_ceiling_check.py`, `HAT_barrier3d_gap_momentum_fix.py`, `HAT_metres_3_overwash_fix.py`, `HAT_metres_3_overwash_fix_plot_explained.py`

**dunes/**: `HAT_cape_growth_rate.py`, `HAT_cape_growth_rate_figure.py`, `HAT_dune_ceiling_callback.py`, `HAT_dune_ceiling_domain_value.py`, `HAT_dune_ceiling_per_domain.py`, `HAT_dune_ceiling_range.py`, `HAT_dune_ceiling_rebuild.py`, `HAT_dune_crest_dem_vs_model.py`, `HAT_dune_design_height.py`, `HAT_dune_findings_figures.py`, `HAT_dune_findings_simple_figures.py`, `HAT_dune_growth_alongshore_sensitivity.py`, `HAT_dune_growth_alongshore_survey_fit.py`, `HAT_dune_growth_beach_width_check.py`, `HAT_dune_growth_rate.py`, `HAT_dune_recovery_diagnosis.py`, `HAT_dune_recovery_rate.py`, `HAT_dune_recovery_slow_growth.py`, `HAT_measured_ceiling_slow_growth.py`, `HAT_shoreline_rescore_new_ends.py`, `HAT_zerobe_test_growth.py`

**ends/**: `HAT_buffer_orientation.py`, `HAT_peaisland_extension.py`, `HAT_position_change_ends_figure.py`, `HAT_resolve_ends_metres.py`, `HAT_resolve_ends_on_position_change.py`

**management/**: `HAT_buxton_2017_fill.py`, `HAT_buxton_2017_fill_plot.py`, `HAT_buxton_2017_shift_north.py`, `HAT_buxton_2017_shift_north_plot.py`, `HAT_fill_footprint_coastsat.py`, `HAT_fill_footprint_coastsat_plot.py`, `HAT_nourishment_volume_check.py`

**offset/**: `HAT_metres_1_offset_units.py`, `HAT_metres_1_offset_units_plot.py`, `HAT_offset_source_0922_figures.py`, `HAT_offset_source_comparison.py`, `HAT_offset_source_comparison_div10.py`, `HAT_offset_source_comparison_option_a.py`, `HAT_offset_source_shoreline_v2.py`

**overwash/**: `HAT_excess_overwash_diagnosis.py`, `HAT_false_alarm_locations.py`, `HAT_isabel_pattern_scores.py`, `HAT_overwash_scoring_tolerance.py`, `HAT_soundside_overwash.py`

**relocation/**: `HAT_relocation_comparison.py`, `HAT_relocation_dune_position_check.py`, `HAT_relocation_period_compare.py`

**storms/**: `HAT_storm_event_splitting.py`, `HAT_storm_height_test.py`, `HAT_storm_length_selection.py`, `HAT_storm_max_duration.py`, `HAT_storm_notebook_rule.py`, `HAT_storm_rlow_duration.py`, `HAT_storm_rlow_duration_plot.py`, `HAT_storm_runup_by_slope.py`, `HAT_trim_length_adopted.py`

**waves/**: `HAT_metres_2_wave_sensitivity.py`, `HAT_metres_2_wave_sensitivity_plot.py`, `HAT_wave_grid_fixed_ends.py`, `HAT_wave_grid_smoothed_score.py`, `HAT_wave_grid_smoothed_score_plot.py`, `HAT_wave_recommendation_figures.py`, `HAT_wave_shortlist_ends_solved.py`

**windows/**: `HAT_candidate_windows.py`, `HAT_candidate_windows_plot.py`, `HAT_full_window_1996_2025.py`, `HAT_window_1996_2010.py`

## The theme READMEs, as they were

The READMEs of the folders removed whole, kept here because they were never committed.

---

## Code checks

Checks that a Barrier3D fix or adoption changes what it should and nothing else.

Each study drives the unchanged runner (`../../HAT_hindcast.py`); see `../README.md` for how
the experiments are run and filed, and how they import each other across folders.

| Script | What it asks |
|---|---|
| `HAT_metres_3_overwash_fix.py` | How much does the Barrier3D route_overwash fix move the results? |
| `HAT_metres_3_overwash_fix_plot_explained.py` | A picture of the Barrier3D route_overwash bug and what fixing it changed. |
| `HAT_barrier3d_gap_momentum_fix.py` | How much do the three Barrier3D overwash fixes move the results? |
| `HAT_adopt_dune_ceiling_check.py` | Does the committed per-cell dune ceiling reproduce the experiment that proposed it? |

### The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

#### HAT_metres_3_overwash_fix.py

How much does the Barrier3D route_overwash fix move the results?

From the script's original header:

```text
HAT_metres_3_overwash_fix.py -- how much does the Barrier3D route_overwash fix move the results?
THE BUG (found 2026-09-24): Barrier3D/barrier3d/barrier3d.py, route_overwash,
the subaerial test indexed Elevation[TS, i, d+1:d+10] (row i, columns d+1..d+9)
where Elevation[TS, d+1:d+10, i] (the nine cells landward of the flow) is
meant: the wrong cells whenever i < rows, out of bounds whenever i >= rows.
The fix is one line, commit 49fd069 on Barrier3D branch
fix/route-overwash-axis-swap (local). Barrier3D is installed editable, so the
branch checked out IS the model every run uses.

THE CHECK (Hannah: "patch it on a branch and measure the impact")
    patched            8 runs on the fix branch, each the twin of a run
                       already made unpatched
    patched_boundscheck the natural 1996 and managed 2010 baselines again with
                       NUMBA_BOUNDSCHECK=1: does the fixed model read out of
                       bounds anywhere else?
    unpatched          the two /10 twins re-run on master with today's code:
                       their archived matrix runs were made on older CASCADE
                       code, so they are not a clean control. The six metres
                       twins were made today with today's code and are.
    compare            runner scores, share of variation explained, and the
                       per-domain LRR difference, patched minus unpatched

    Every launch records the Barrier3D branch and commit it ran on, and
    refuses to run a patched member off the fix branch or an unpatched one on it.

WHERE: output/raw_runs/experiments/code-checks/2026-09-24-metres-3-barrier3d-overwash-fix/
           runs/<variant>_<member>/<period>/<preset>/<run_name>/   runs (on disk only)
           logs/ (+ launches.jsonl), tables/comparison.csv, NOTE.md

USAGE
    python HAT_metres_3_overwash_fix.py run patched            (fix branch checked out)
    python HAT_metres_3_overwash_fix.py run unpatched          (master checked out)
    python HAT_metres_3_overwash_fix.py compare
```

Notes that were in the code:

```text
member -> (period, scenario, wave settings or None for the model defaults,
extra env, the unpatched twin's run folder or None)
```

#### HAT_metres_3_overwash_fix_plot_explained.py

A picture of the Barrier3D route_overwash bug and what fixing it changed.

From the script's original header:

```text
HAT_metres_3_overwash_fix_plot_explained.py -- a picture of the Barrier3D route_overwash bug
Hannah asked for a figure that shows the bug plainly. Five panels:

  (a) the question the model means to ask: from the cell carrying overwash sand,
      look at the nine cells LANDWARD of it, down its own column
  (b) what the code actually looks at: row and column swapped, so a strip
      ACROSS the island somewhere else
  (c) on a narrow island (fewer rows than columns) that strip lies off the
      grid: the read is outside the island's memory
  (d) every compared run's score with the bug and fixed
  (e) what fixing it changes along the island in the run it moved most,
      natural 2010-2024 (experiments/code-checks/2026-09-24-metres-3-barrier3d-overwash-fix)

The grids in (a)-(c) are schematic (a small domain, not to scale); the
indexing is the real one from barrier3d.py line 1092.

Output: output/raw_runs/experiments/code-checks/2026-09-24-metres-3-barrier3d-overwash-fix/figures/
        route_overwash_bug_explained.png
```

<details><summary>Function notes (the original docstrings)</summary>

**`title()`**

```text
Letter and title as one left-aligned line: the house _title centres the
title, which collides with the letter on these narrow panels.
```

**`panel_offgrid()`**

```text
A narrow island: fewer rows (8) than the column number of the cell (12),
so the swapped read, row 12, lies below the end of the grid.
```

</details>

#### HAT_barrier3d_gap_momentum_fix.py

How much do the three Barrier3D overwash fixes move the results?

From the script's original header:

```text
HAT_barrier3d_gap_momentum_fix.py -- how much do the three overwash fixes move the results?
THE FIXES (found 2026-09-27 with the storm replay,
scripts/figure_making/model/storm_replay.py; committed 2026-09-28 on Barrier3D
branch fix/overwash-gaps-momentum, 2.0.2.dev1, LOCAL ONLY - not pushed):
    990c3bd  DuneGaps dropped the last overtopped cell of the last gap and any
             single-cell gap
    015f11e  gap discharge was set on start:stop although stop is inclusive
    e929e65  the inundation momentum constant C was reset to 0 before routing
             (upstream b11b880, the 2024 Numba refactor)

THE BRANCH IS NOT CHECKED OUT in the main Barrier3D repository, which stays on
fix/route-overwash-axis-swap (49fd069), the code every matrix run used. It is
a git worktree at ../Barrier3D-overwashfix, and a fixed run reaches it through
PYTHONPATH, which puts the worktree ahead of the editable install. The hindcast
runner, notebook and config are NOT modified: the run records the Barrier3D it
actually imported (run_registry.barrier3d_provenance), and this driver refuses
a run whose log shows any other.

THE CHECK
    fixed      four runs on the fix branch, each the twin of a matrix run made
               2026-09-27 on the same CASCADE code (edgeBE, option A waves):
               natural and full_management, 1996-2010 and 2010-2024
    controls   the matrix runs themselves. A re-run of the natural 1996 run on
               today's unchanged code reproduced its shoreline matrix to 0.0 m
               (2026-09-28), so they are clean controls.
    compare    net shoreline change and LRR skill per domain, and overwash:
               domain-years with overwash, and the observed-imagery hit rate
               (8-overwash-analysis/4-vs-model), fixed against control

WHERE: output/raw_runs/experiments/code-checks/2026-09-28-barrier3d-overwash-gap-momentum-fix/
           runs/fixed_<member>/<period>/edgeBE/<run_name>/     the runs
           logs/, tables/comparison.csv, figures/, NOTE.md

USAGE
    python HAT_barrier3d_gap_momentum_fix.py run [--workers 4]
    python HAT_barrier3d_gap_momentum_fix.py compare
```

#### HAT_adopt_dune_ceiling_check.py

Does the committed per-cell dune ceiling reproduce the experiment that proposed it?

From the script's original header:

```text
HAT_adopt_dune_ceiling_check.py -- does the committed per-cell ceiling reproduce the experiment?
Step 4 of adopting the per-cell dune ceiling (Hannah, 2026-09-28).

    reproduce  the four trim24 per-cell cases of
               storms-overwash-and-dunes/2026-09-28-dune-ceiling-per-domain
               (1996/2010 x full_management/natural) run on Barrier3D branch
               feature/per-cell-dune-ceiling (worktree ../Barrier3D-dune-ceiling,
               PYTHONPATH) with DuneCeilingFromStart switched on, against the
               experiment's runs, which set the same ceilings by wrapping the
               model in-process. The shoreline matrix and every domain's dune
               domain must be identical.
    off        with the setting off, the branch must equal the code every run
               uses now (49fd069): the natural 1996 run against the matrix run's
               shoreline matrix is not usable (the matrix was re-run on other
               ends), so the unit-level check is in the Barrier3D tests.

Nothing in the main code changes: the setting reaches the model through
cascade.brie_coupler.set_yaml in the run's own process, and the storm file is
swapped as in the experiments.

WHERE: output/raw_runs/experiments/code-checks/2026-09-28-per-cell-dune-ceiling-reproduces/
```

---

## Management

The 2017 Buxton fill, nourishment volumes and the fill footprint against CoastSat.

Each study drives the unchanged runner (`../../HAT_hindcast.py`); see `../README.md` for how
the experiments are run and filed, and how they import each other across folders.

| Script | What it asks |
|---|---|
| `HAT_buxton_2017_fill.py` | Does adding the 2017 Buxton fill improve the 2010-2026 hindcast? |
| `HAT_buxton_2017_fill_plot.py` | Before and after the 2017 Buxton fill: model rates with and without it against CoastSat, 2010-2026. |
| `HAT_buxton_2017_shift_north.py` | Does moving the 2017 Buxton fill two domains north (GIS 6-16 -> 8-18) improve the 2009-2025 test run? |
| `HAT_buxton_2017_shift_north_plot.py` | Net change 2009-2025 with the 2017 Buxton fill shifted two domains north, in the run-figure layout. |
| `HAT_nourishment_volume_check.py` | Are the beach nourishment volumes applied as intended, for every fill in both windows? |
| `HAT_fill_footprint_coastsat.py` | Does placing each fill on its CoastSat-observed footprint improve the 2009-2025 test run? |
| `HAT_fill_footprint_coastsat_plot.py` | Net change 2009-2025 with the fills on their reported and on their CoastSat-observed footprints. |

### The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

#### HAT_buxton_2017_fill.py

Does adding the 2017 Buxton fill improve the 2010-2026 hindcast?

```text
Runs 2010-2026 edgeBE full_management twice: with the fill list as committed
(with2017, which must reproduce the matrix run) and with Buxton 2017 dropped
from HATTERAS_NOURISHMENT_PROJECTS in the run's own process (without2017).
Scored against the CoastSat 2010-2026 LOWESS-7 target over the interior and
over the fill footprint, GIS 6-15.

NOTHING IN THE MAIN CODE CHANGES.

    python HAT_buxton_2017_fill.py run
    python HAT_buxton_2017_fill.py score

WHERE: output/raw_runs/experiments/management/2026-10-03-buxton-2017-fill/
```

#### HAT_buxton_2017_fill_plot.py

The before-and-after figure for HAT_buxton_2017_fill.py: model rates with and
without the fill against CoastSat over GIS 1-30, and the residual per domain.
`--smoothed` passes the model rates through the target's own smoothing first
(7-domain LOWESS, raw at GIS 1-10).

```text
    python HAT_buxton_2017_fill_plot.py [--smoothed]
```

#### HAT_nourishment_volume_check.py

Are the beach nourishment volumes applied as intended, for every fill in both windows?

```text
Follows each fill from Hatteras_BN_data.xlsx to the model shoreline: sheet
volume and length against the site config, cubic yards to m^3/m, the volume
each BeachDuneManager recorded against the schedule, the shoreline step
against 2V / (2 h_b + D_sf), and the total placed. Reads saved
full_management runs (1996-2015 from its own folder, 2010-2026 from the
Buxton experiment); runs nothing. A step larger than the formula in a
community domain is overwash the filter returned to the shoreface in the
same model step, not a fill error.

    python HAT_nourishment_volume_check.py

WHERE: output/raw_runs/experiments/management/2026-10-03-nourishment-volume-check/
```

---

## Island offset

How the island offset is set (units, source) and how the shoreline it gives compares. `HAT_metres_1_offset_units.py` is also the shared helper module many other scripts import as `common`.

Each study drives the unchanged runner (`../../HAT_hindcast.py`); see `../README.md` for how
the experiments are run and filed, and how they import each other across folders.

| Script | What it asks |
|---|---|
| `HAT_metres_1_offset_units.py` | Island-offset scale against wave-climate tuning: can a full-scale offset be re-tuned to match? |
| `HAT_metres_1_offset_units_plot.py` | The figures for the offset-scale study. |
| `HAT_offset_source_comparison.py` | How much does the output change when the island offset comes from the dune line rather than the CoastSat shoreline? |
| `HAT_offset_source_comparison_div10.py` | Dune line against shoreline as the island offset, rebuilt on the old ÷10 offset. |
| `HAT_offset_source_comparison_option_a.py` | Does the 09-25 offset-source result hold at the option A waves? |
| `HAT_offset_source_shoreline_v2.py` | Shoreline offset v1 against v2, beside the dune line, on the adopted setup. |
| `HAT_offset_source_0922_figures.py` | The 2026-09-22 shoreline-offset study, drawn in the island-offset house form. |

### The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

#### HAT_metres_1_offset_units.py

Island-offset scale against wave-climate tuning: can a full-scale offset be re-tuned to match?

From the script's original header:

```text
HAT_metres_1_offset_units.py -- island-offset scale against wave-climate tuning
THE QUESTION (Hannah, 2026-09-24). BRIE's shoreline is in metres, and so is
the island-offset file, but every calibrated run hands BRIE the offset divided
by ten (`offset_mode: asrun`, a units error). Put back at full scale, can the
model be re-tuned through its wave climate to match the /10 runs' skill?

TWO SWEEPS, ONE STUDY
    wave_height  Hs 0.5-3.0 m at the default wave angles, every offset
                 scale, dune-line and shoreline sources
    wave_angle   wave asymmetry 0.5-0.8 x high-angle fraction 0.1-0.4 at
                 Hs 1.0 m, every offset scale, dune-line source
    fixed        1996-2010, zeroBE, full_management, no groin, relocations
                 off, base geometry, Tp 8 s

OFFSET SCALE -- the folder label, and the value the runner reads
    div10              HAT_OFFSET_MODE=asrun      offset / 10 (every calibrated run)
    metres             HAT_OFFSET_MODE=metres     the measurement as is
    metres-detrended   HAT_OFFSET_MODE=detrended  metres, linear trend removed

LAYOUT (output/raw_runs/experiments/<STUDY>/)
    README.md                     question, design, how to read labels, results
    tables/                       written by `score`, every setting a column;
                                  observed_target.csv: the target's mean, sd
                                  and flat-line RMSE
    figures/<sweep or combined>/  written by HAT_metres_1_offset_units_plot.py
    logs/<sweep>/<scale>_<source>/Hs<h>_asymmetry<a>_highangle<f>.log
    logs/drivers/                 this script's own console logs
    runs_<sweep>/<scale>_<source>/1996_2010/zeroBE/<run_name>/

    The run folder's tag is <STUDY>/runs_<sweep>/<scale>_<source>, the
    runner's three-level maximum. The run NAME is the runner's and leaves out
    whatever is at its default (no offset token for div10, no asym token at
    0.7, no ahf token at 0.1); the folder, the log name and the tables always
    carry every setting.

USAGE
    python HAT_metres_1_offset_units.py run wave_height [--scales ...] [--hs ...]
    python HAT_metres_1_offset_units.py run wave_angle
    python HAT_metres_1_offset_units.py score
    (--jobs N, --dry-run, --overwrite on `run`)
```

Notes that were in the code:

```text
Run after the first grid, as their own commands (see README):
--scales metres metres-detrended --hs 0.5 0.75
--scales div10 --hs 0.5
```

```text
metres was still improving at 0.4, so extended (2026-09-24):
--scales metres --high-fraction 0.45 0.5
```

```text
A clean log means the run exists; the runner would refuse it anyway,
and that refusal used to read as a failure.
```

```text
THE TARGET AND THE ALONGSHORE SCORES

The runner scores RMSE and bias; these add how much of the observed
alongshore variation a run explains. The target is rebuilt exactly as the
runner builds it (section 8 of HAT_hindcast.py), and every run's
RMSE is recomputed from it and checked against the runner's, so the two
sets of scores are provably on the same target.
```

```text
`end` scores a window inside a period (2010-2020 inside 2010-2024,
2026-09-24): the CoastSat table must exist under lrr/<start>_<end>/.
```

<details><summary>Function notes (the original docstrings)</summary>

**`run_env()`**

```text
The environment one run reads, built the way HAT_run_matrix builds it:
every HAT_* variable named here, none inherited from the shell.
```

**`coastsat_target()`**

```text
The CoastSat LRR target, GIS 1-90, as the runner builds it for a start
year (section 8 of the runner: LOWESS at 10 domains, the southern 10 raw).
Shared with scripts/real-world-cases/Hatteras/hatteras_ms/experiments/waves/HAT_metres_2_wave_sensitivity.py.
```

**`smooth_like_target()`**

```text
A model series (per GIS domain) smoothed as the CoastSat target is:
LOWESS over a SMOOTH_DOMAINS window (7 since 2026-09-28, 10 before; frac
SMOOTH_DOMAINS/90 on the 90 domains, matching the
target's 0.110), the southern 10 domains left raw as the target leaves
them. Added 2026-09-25 for the smoothed score (Hannah); shared by the
step-2 figures and wave-climate/2026-09-25-wave-grid-smoothed-score.
```

**`alongshore_scores()`**

```text
How much of the observed alongshore variation a run explains.

variance_explained          1 - sum((m - o)^2) / sum((o - mean o)^2).
                            1 is perfect; 0 is no better than a flat
                            line at the observed mean; negative is worse.
                            Bias counts against it.
pattern_variance_explained  the same with each series' own mean removed
                            first: the pattern alone, bias forgiven.
r_alongshore                correlation of the two alongshore series
sd_ratio                    model sd / observed sd: below 1 the model
                            varies less along the island than observed
```

</details>

#### HAT_metres_1_offset_units_plot.py

The figures for the offset-scale study.

From the script's original header:

```text
HAT_metres_1_offset_units_plot.py -- the figures for the offset-scale study
Reads tables/all_runs.csv (written by `HAT_metres_1_offset_units.py score`),
each scored run's tables/shoreline_change_rate.csv, and the CoastSat target the
runner scores against. Writes, under figures/ in the study folder:

  wave_height/rmse_bias_vs_wave_height_by_offset_scale_1996_2010.png
      the offset each scale hands BRIE; interior RMSE and bias against Hs
  wave_angle/rmse_bias_vs_high_angle_fraction_by_asymmetry_1996_2010.png
      RMSE and bias against the high-angle fraction, one line per
      asymmetry, one column per offset scale
  combined/best_rmse_by_offset_scale_1996_2010.png
      the RMSE range each scale reaches in each sweep, best run labelled
  combined/bias_vs_rmse_all_runs_1996_2010.png
      every scored run of both sweeps, bias against RMSE
  combined/alongshore_rate_best_runs_vs_coastsat_1996_2010.png
      the best run of each scale against the CoastSat target, GIS 1-90
  combined/variance_explained_by_offset_scale_1996_2010.png
      the share of the observed alongshore variation each scale explains
  combined/spread_vs_correlation_all_runs_1996_2010.png
      each run's alongshore spread against its correlation with CoastSat
  alongshore_sensitivity/<scale>/rate_and_position_change_by_<parameter>_<scale>_1996_2010.png
      one parameter moved (Hs, high-angle fraction, asymmetry), the others at
      their defaults: (a) rate vs the CoastSat LRR target, (b) position change
      vs the observed CoastSat change 1996 -> 2010

Scores are the runner's (interior GIS 2-89, LRR, CoastSat LOWESS 7-domain since 2026-09-28, 10 before
target), plus the alongshore-variation scores `score` adds from the same
target (study.coastsat_target, checked there against the runner's RMSE).
```

Notes that were in the code:

```text
drowned cells, at the foot of the panel, so a missing point says why
drowned cells, at the foot of the panel in the scale's colour, one
row per scale, so a missing point says which run and why
```

```text
label every half metre plus 0.75; the other Hs run values (0.6,
0.65) are unlabelled minor ticks, too close to label
```

```text
the runner's RMSE, reproduced from the drawn curves: a check that
this figure shows what was scored
```

```text
Lines of equal pattern skill (bias removed): skill = 2 r s - s^2 with s
the sd ratio, so r = (skill + s^2) / (2 s).
```

```text
Observed position change over the window, from 5-scr: mean CoastSat position
over calendar 2010 minus calendar 1996, seaward positive, LOWESS-smoothed at
10 domains to match the rate target's window.
```

```text
One parameter moved, the other two at the calibration defaults
(Hs 2.5 is not a default here: the wave-angle sweep runs at Hs 1.0).
```

<details><summary>Function notes (the original docstrings)</summary>

**`flat_line_rmse()`**

```text
RMSE of predicting the observed interior mean at every domain: the
score a model with no alongshore pattern at all would get.
```

</details>

#### HAT_offset_source_comparison.py

How much does the output change when the island offset comes from the dune line rather than the CoastSat shoreline?

From the script's original header:

```text
Dune line vs shoreline as the island offset (orientation), 1996-2010 (2026-09-25).

Asked by Hannah on 2026-09-25: how much does setting the island's planform
from the dune line rather than the CoastSat shoreline change the output?

    offsets   duneline (1996/duneline/v1) and shoreline (1996/shoreline/v1),
              both metres, the Hermite wrap-around written in the file
    waves     Hs 1.0 m, Tp 8 s, asymmetry 0.8 -- the best managed 1996-2010
              setting found (tuned with the DUNE-LINE offset) -- and, so the
              shoreline offset gets a fair chance, high-angle fraction
              0.3 0.4 0.45 0.5 0.55 (the lever that mattered most)
    scope     natural and full management, 1996-2010
    score     as wave-climate/2026-09-25-wave-grid-smoothed-score: share of the alongshore
              variation explained by the model SMOOTHED like the CoastSat
              target, interior GIS 2-89; raw score, bias and r beside it
Both offsets are run fresh here (20 runs) so every run in the comparison is
on the same code and the same Barrier3D (the route_overwash fix).

WHERE: output/raw_runs/experiments/island-offset/2026-09-25-metres-offset-duneline-vs-shoreline-waves-hs1-tp8/
    README.md, tables/all_runs.csv, figures/, logs/<source>_<scenario>/<settings>.log
    runs/<source>_<scenario>/1996_2010/zeroBE/<run_name>/   (on disk only)

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/offset/HAT_offset_source_comparison.py run
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/offset/HAT_offset_source_comparison.py score
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/offset/HAT_offset_source_comparison.py plot
```

Notes that were in the code:

```text
Figure labels, overridden by a study that reuses this driver
(HAT_offset_source_comparison_div10.py).
```

```text
a superseded build records the source alone
(hatteras_site_config.island_offset_version)
```

```text
Text sized for reading the figures side by side (Hannah, 2026-09-29:
"make all the text larger"); the canvas stays 16 in wide.
```

```text
One y-axis for every island-offset study (09-22, 09-25, 09-28 option A),
set from the largest of them: the change figures, and the difference one
```

```text
legend_title (the end correction) goes to the caption, not the canvas:
it read as a heading for the legend entries (2026-09-29)
```

```text
Every observed and modelled line first, so the three comparison figures
share ONE y-axis and read against each other (Hannah, 2026-09-29).
```

```text
No headroom above the data: the one value past the label line is the
observed +96 m at GIS 1-2 (2010-2024), where no label sits; headroom for
it pushed the top to 160 m and flattened every line.
```

```text
its own quantity, so its own axis, with room above the bars for
the village and shoal labels
```

```text
THE FIGURES ARE FULL MANAGEMENT x BOTH PERIODS (Hannah, 2026-09-28: "showing
only full management and both periods per figure", as the option A study
draws them). 1996's run is the study's own; 2010's comes from `run-2010`.
```

```text
Observations for the house-form figures: smoothed at 7 domains, the research
group's range (Hannah, 2026-09-28), the southern 10 domains left raw.
```

<details><summary>Function notes (the original docstrings)</summary>

**`house_figures()`**

```text
The four island-offset figures, in ONE form for every study that asks the
dune-line-or-shoreline question (Hannah, 2026-09-28: "ensure the figures
among these experiments are consistent ... so it is easier to compare").

panels   [{"label": "Natural 1996–2010", "start": 1996, "end": 2010,
           "rates": {"duneline": df, "shoreline": df}}, ...], one row each;
           df is a run's tables/shoreline_change_rate.csv indexed by
           gis_domain, or None where that arm has no run
note     the study's own sentence for every caption (waves, offset, ends)
suffix   the stem ending, e.g. "1996_2010" or "full_management"
shoreline_window  the mean-shoreline windows the shoreline offset was built
         on, for the captions; the default is the v1 (calendar) builds

Net change in metres; observations smoothed with LOWESS over LOWESS_DOMAINS
(southern SKIP_SOUTHERN raw), the model unsmoothed; the model's ENABLED
fills marked above each panel; no scores on the figures.
```

**`cmd_plot()`**

```text
This study's figures through house_figures: (a) full management
1996-2010, (b) full management 2010-2024, at the headline setting.
```

**`coastsat_target_lowess()`**

```text
The CoastSat LRR target built as the runner builds it, at LOWESS_DOMAINS,
the rate fitted on `window` ("1996_2010", or "1996_2024" for the long-term rate).
```

**`smooth_lowess()`**

```text
A per-domain series smoothed as the target is: LOWESS over LOWESS_DOMAINS,
the southern SKIP_SOUTHERN left raw.
```

**`limits()`**

```text
The FIXED y-limits, the same in every island-offset study so the
figures compare across studies (Hannah, 2026-09-29), widened to `step`
only if a study's data leaves them, with a warning.
```

</details>

#### HAT_offset_source_comparison_div10.py

Dune line against shoreline as the island offset, rebuilt on the old ÷10 offset.

From the script's original header:

```text
Dune line vs shoreline as the island offset, on the OLD ÷10 offset (2026-09-28).

Asked by Hannah on 2026-09-28: rebuild the offset-source experiment on the ÷10
island offset from before the metres fix, clearly labelled, "because I am
trying to track the changes as I was experimenting".

    offsets   the builds from before the metres fix, as the runs of the time
              read them: offset_mode "asrun" (offset / 10, the historical unit
              error) on superseded_20260924_pre-metres/v1 of each source.
              asrun divides the buffers too, so only the build a ÷10 run was
              made from reproduces it (1996/shoreline/v1/PROVENANCE.md)
    waves     the ÷10-era defaults, the calibration the model ran at before
              2026-09-24: Hs 2.5 m, Tp 8 s, asymmetry 0.7, high-angle 0.1
              (Hannah's choice). Offset AND waves differ from the 09-25 study,
              so the two cannot separate their effects
    ends      zeroBE; relocations and groins off, as on 09-25
    scope     natural and full management, 1996-2010; 4 runs
    Barrier3D the current one, with the route_overwash fix. The ÷10-era runs
              used the unfixed router, which segfaulted on some storms
    figures   the 09-25 driver's house form: metres of net change, observations
              LOWESS 7, no scores on the figures

WHERE: output/raw_runs/experiments/island-offset/2026-09-28-div10-offset-duneline-vs-shoreline-rebuild/

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/offset/HAT_offset_source_comparison_div10.py run --jobs 4
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/offset/HAT_offset_source_comparison_div10.py score
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/offset/HAT_offset_source_comparison_div10.py plot
```

Notes that were in the code:

```text
2010: the dune line has a pre-metres build; the 2010 shoreline offset was
first built 2026-09-28, so asrun reads its v1. On the 1996 builds the two
differ only in the 15+15 buffer domains, not the real ones.
```

#### HAT_offset_source_comparison_option_a.py

Does the 09-25 offset-source result hold at the option A waves?

From the script's original header:

```text
Dune line vs shoreline as the island offset, re-run at the option A waves (2026-09-28).

Asked by Hannah on 2026-09-28: does the 09-25 result (shoreline offset matches
the dune line under management, beats it in the natural run) hold at the
waves adopted on 09-27?

    offsets   duneline (1996/duneline/v1) and shoreline (1996/shoreline/v1),
              both metres. There is no 2010 shoreline offset, so 1996-2010 only
    waves     option A: Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5 --
              one setting, no sweep (Hannah's choice)
    ends      zeroBE, as on 09-25: option A's edgeBE ends were solved on the
              dune-line offset and would favour it
    scope     natural and full management; 4 runs, relocations and groins off
    score     RAW share of the alongshore variation explained, interior GIS 2-89
              (the score option A was chosen on); smoothed, bias and r beside it

Everything but the settings and the headline score is the 09-25 driver
(HAT_offset_source_comparison), pointed at this study's folder.

WHERE: output/raw_runs/experiments/island-offset/2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a/

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/offset/HAT_offset_source_comparison_option_a.py run --jobs 4
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/offset/HAT_offset_source_comparison_option_a.py score
```

Notes that were in the code:

```text
EACH OFFSET GRADED ON ITS OWN FEATURE (Hannah, 2026-09-28). A run started
from the dune line is scored against the dune line's change, a run started
from the shoreline against the shoreline's; full management, both periods.
dune line   net change in metres: the observed dune-line endpoint change
(1997->2009, 2009->2023 as surveyed, 11.6 and 14.1 yr) against
the model's endpoint change over its 14 calendar years. The
interval mismatch is reported, not corrected
([[cascade-period-is-the-calendar-year]]).
shoreline   the CoastSat LRR target against the model's LRR, as before.
The two scores are on different targets and different estimators (Hannah's
choice), so they rank each offset against its own feature; they are not a
head-to-head.
```

```text
The research group's alongshore smoothing range is 7 domains (Hannah,
2026-09-28), for the CoastSat target and the dune line alike. The runner and
the 09-25 helpers still smooth at 10; this study sets 7 here and leaves them.
```

```text
smoothed as the CoastSat target is: LOWESS over 7 domains, the
southern 10 raw (Hannah, 2026-09-28); the model stays raw
```

```text
PROJECTED shoreline change: the 1996-2024 LRR x 14 yr, the same
observed profile in both periods, against the same model runs
```

<details><summary>Function notes (the original docstrings)</summary>

**`coastsat_target7()`**

```text
The CoastSat LRR target built as the runner builds it, at 7 domains.
`window` fits the LRR on another span ("1996_2024" for the long-term rate);
default the period's own.
```

**`smooth7()`**

```text
A per-domain series smoothed as the target is: LOWESS over 7 domains,
the southern 10 left raw (common.smooth_like_target, at 7).
```

**`cmd_plot_own()`**

```text
This study's figures through the shared house_figures (09-25 driver):
(a) full management 1996-2010, (b) full management 2010-2024.
```

</details>

#### HAT_offset_source_shoreline_v2.py

Shoreline offset v1 against v2, beside the dune line, on the adopted setup.

From the script's original header:

```text
Shoreline offset v1 vs v2, beside the dune line, on the adopted setup (2026-09-29).

Asked by Hannah on 2026-09-29: re-run the shoreline arm on shoreline offset
v2, the CoastSat mean over +/-1 yr of the start DEM's lidar flights
(1995-10-12..1997-10-12 for 1996, 2008-08-17..2010-08-17 for 2010), which
became CURRENT that day. Every shoreline run before then read v1 (the
calendar means, 1995-1997 and 2009-2011).

The option A study (2026-09-28) ran before three adopted changes -- the
per-cell dune ceilings and the beach/dune cap fix (09-28, later that day)
and the split12 storm files (09-29) -- so its v1 runs cannot be set beside
a v2 run made now. Hannah chose a clean three-way study instead: all three
offsets on today's setup, so v1 -> v2 is the only difference between the
shoreline arms.

    arms      duneline     2-brie-offset/<year>/duneline/CURRENT (v1)
              shoreline_v1 2-brie-offset/<year>/shoreline/v1, pinned by
                           HAT_OFFSET_VERSION_<year>_SHORELINE=v1
              shoreline_v2 2-brie-offset/<year>/shoreline/v2 (CURRENT)
    waves     option A (Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5)
    ends      zeroBE, as in the option A study (edgeBE ends were solved on
              the dune-line offset and would favour it)
    scope     natural and full management, 1996-2010 and 2010-2024;
              relocations and groins off; 12 runs
    setup     everything else is the code default: dune ceilings, storm
              series (v3_split12_trim24), beach/dune cap

Every run's metadata is checked for the offset it actually read.

SCORES, as in the option A study
    vs CoastSat   each period's own CoastSat LRR, LOWESS 7 domains: raw and
                  smoothed share of the alongshore variation explained,
                  interior GIS 2-89, bias and r (the option A headline)
    own feature   dune-line arm against dune-line net change (m); shoreline
                  arms against total (own LRR x 14 yr) and projected
                  (1996-2024 LRR x 14 yr) shoreline change

WHERE: output/raw_runs/experiments/island-offset/2026-09-29-shoreline-offset-v1-vs-v2-adopted-setup/

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/offset/HAT_offset_source_shoreline_v2.py run --jobs 4
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/offset/HAT_offset_source_shoreline_v2.py score
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/offset/HAT_offset_source_shoreline_v2.py plot
```

<details><summary>Function notes (the original docstrings)</summary>

**`cmd_plot()`**

```text
The house figures (full management, both periods), once per shoreline
version, and the v2-minus-v1 difference.
```

**`diff_figure()`**

```text
Model total change on shoreline v2 minus v1, per domain: natural and
full management, both periods, on the difference figures' fixed axis.
```

</details>

#### HAT_offset_source_0922_figures.py

The 2026-09-22 shoreline-offset study, drawn in the island-offset house form.

From the script's original header:

```text
The 2026-09-22 shoreline-offset study, drawn in the island-offset house form (2026-09-28).

Hannah, 2026-09-28: make the figures across the island-offset experiments
consistent with the most recent version, so they are easier to compare. The
09-22 study had only the runner's per-run figures. This draws its runs with
the shared house_figures (HAT_offset_source_comparison.py) and runs nothing.

    shoreline  the study's own run,
               island-offset/2026-09-22-div10-offset-shoreline-trial-original/
               1996_2010/edgeBE/HAT_1996_2010_edgeBE_road_bdm_nogroin
    dune line  its matrix control, since archived:
               archive/2026-09-24-pre-metres/matrix/1996_2010/edgeBE/
               HAT_1996_2010_edgeBE_road_bdm_nogroin
    both       ÷10 offset (asrun) on the pre-metres builds, Hs 2.5 / Tp 8 /
               asymmetry 0.7 / high-angle 0.1, edgeBE +32.2 / +10.0 m/yr,
               full management, no groin, the unfixed Barrier3D. The control
               ran 2026-09-18 (7e492c3), the shoreline arm 2026-09-22 (c8ed4a1)

The road_reloc_bdm arm is left out: no relocation falls inside 1996-2010, so
it is the same run (NOTE.md), and it has no archived control.

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/offset/HAT_offset_source_0922_figures.py
```

Notes that were in the code:

```text
(b) 2010-2024 was never run in this study: the panel carries the
observation alone so the layout matches the others (full management x
both periods, Hannah 2026-09-28); no run is added to a record
```

---

## Overwash scoring

How modelled overwash is scored against the imagery, and why the model overwashes more than observed.

Each study drives the unchanged runner (`../../HAT_hindcast.py`); see `../README.md` for how
the experiments are run and filed, and how they import each other across folders.

| Script | What it asks |
|---|---|
| `HAT_excess_overwash_diagnosis.py` | Why does the model overwash more than the imagery shows? |
| `HAT_overwash_scoring_tolerance.py` | How much do the overwash scores depend on how strictly modelled overwash is matched to the imagery? |
| `HAT_soundside_overwash.py` | How much of the observed overwash could not have come from the ocean side, and so is likely sound-side? |
| `HAT_false_alarm_locations.py` | Where does slower dune growth add its false alarms: on dunes the surveys show were low, or elsewhere? |
| `HAT_isabel_pattern_scores.py` | How well does each dune setting reproduce where Isabel overwashed, beside how well it follows the surveyed dunes? |

### The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

#### HAT_excess_overwash_diagnosis.py

Why does the model overwash more than the imagery shows?

From the script's original header:

```text
HAT_excess_overwash_diagnosis.py -- why does the model overwash more than the imagery shows?
THE QUESTION (Hannah, 2026-09-28). Every storm series, the committed one
included, overwashes far more domains than the observed record in most image
windows (storms-overwash-and-dunes/2026-09-28-storm-length-selection, stage 1).

THREE EXPLANATIONS, EACH TESTED DIRECTLY (read-only: existing runs and inputs)
    1  volume   the extra overwash is real but too small to see in imagery
    2  dunes    the model's dune row is lower than the real foredune, either
                from the extraction (a clipped search window leaves the true
                crest in interior row 0, behind a lower dune row) or from the
                dunes changing during the run
    3  storms   Rhigh (Stockdon R2 on WIS Hs, slope 0.06) clears the dunes by
                too much

    `hidden`  per domain, the start-of-run gap between the dune row's crest
              and the highest of the dune row + first N interior rows: the
              foredune height the dune row does not carry
    `cells`   per image x domain cell (managed runs, drop72 and trim24): the
              storms credited with its overwash, their Rhigh, the pre-storm
              crest (after growth, as Barrier3D tests it), the margin, the
              volume, and the domain's hidden crest
    `summary` how false alarms and hits differ on each of those

WHERE: output/raw_runs/experiments/storms-overwash-and-dunes/2026-09-28-excess-overwash-diagnosis/
```

<details><summary>Function notes (the original docstrings)</summary>

**`hidden()`**

```text
Start-of-run dune-row crest vs the highest ground in the dune row plus
the first N_ROWS interior rows, per domain (m MHW, alongshore medians).
```

**`storm_table()`**

```text
Per (pad, storm): Rhigh, pre-storm crest (after growth), margin, and the
overwash share credited to it (S.storm_shares).
```

</details>

#### HAT_overwash_scoring_tolerance.py

How much do the overwash scores depend on how strictly modelled overwash is matched to the imagery?

Re-scores existing runs, no new runs. The imagery is about one image a year at
irregular dates (1998-10 to 2004-05 is one window); the model records overwash
per model year. Rules, each adding to the last: strict (current: a model year
is dated by its largest storm); any storm in the model year; one image either
side; one domain either side. Tolerance is symmetric (a false alarm needs no
observed overwash nearby). "Wrong place" is the same model output shifted
alongshore by 10-80 domains, averaged: a score that survives the shift is
timing, not location.

Results 2026-10-13 (`output/raw_runs/experiments/storms-overwash/2026-10-13-overwash-scoring-tolerance/scores.csv`), skill / skill with the overwash in the wrong place:

| 1996-2009 | strict | any storm | +1 image | +1 domain |
|---|---|---|---|---|
| single 5.5 m | 0.594 / 0.557 | 0.589 / 0.552 | 0.714 / 0.661 | 0.810 / 0.766 |
| per-cell | 0.575 / 0.455 | 0.573 / 0.412 | 0.710 / 0.537 | 0.782 / 0.668 |
| per-domain median | 0.460 / 0.206 | 0.454 / 0.188 | 0.561 / 0.242 | 0.801 / 0.432 |

| 2009-2025 | strict | any storm | +1 image | +1 domain |
|---|---|---|---|---|
| single 5.5 m | 0.120 / 0.108 | 0.054 / 0.052 | 0.113 / 0.104 | 0.258 / 0.251 |
| per-cell | 0.195 / 0.170 | 0.020 / 0.009 | 0.147 / 0.101 | 0.436 / 0.383 |
| per-domain median | 0.111 / 0.102 | 0.001 / 0.013 | 0.088 / 0.075 | 0.242 / 0.253 |

Most of the skill score is timing: moving the single-ceiling run's overwash to
the wrong domains costs 0.04 of its 0.59. Looser matching raises every score
and the wrong-place score with it, so it adds no information about location.
"Any storm in the model year" lowers the 2009-2025 scores, because a model
year's overwash then counts in both of two images a few months apart (2009,
2014, 2023).

The imagery records NEW overwash since the previous image (Hannah, 2026-10-13),
so the window-per-image comparison is the right frame and old fans are never
re-counted; dating a model year by its largest storm stays.

Timing and location, scored apart (`timing_location_scores.csv` in the same
folder). Timing r: per image, share of domains with new overwash, observed
against modelled. Location r: per domain, share of images with new overwash.

| | 1996-2009 timing r | location r | modelled / observed share | 2009-2025 timing r | location r | modelled / observed share |
|---|---|---|---|---|---|---|
| single 5.5 m (adopted) | 0.98 | 0.29 | 0.125 / 0.095 | 0.24 | -0.06 | 0.098 / 0.227 |
| single 12.1 m (DEM max) | 0.23 | 0.12 | 0.015 | 0.37 | 0.04 | 0.065 |
| per-cell | 0.85 | 0.64 | 0.230 | 0.04 | 0.10 | 0.263 |
| per-domain median | 0.86 | 0.48 | 0.113 | -0.07 | 0.01 | 0.156 |

1996-2009 has 8 images and its timing r rests mostly on the one after Isabel.

#### HAT_soundside_overwash.py

How much of the observed overwash could not have come from the ocean side, and so is likely sound-side?

The overwash record has no direction, and Barrier3D overwashes only from the
ocean (storm Rhigh against the dune crest; bay water never enters). For each
observed new overwash (image, domain): the highest ocean run-up in the image's
window (storm summary, m NAVD88) against the domain's LOWEST surveyed crest at
the period start (1996 or 2009 dune file). Run-up below it makes ocean
overwash implausible. NOAA sound-gauge peaks (Oregon Inlet 8652587, USCG
Hatteras 8654467) for named storms mark high sound water (>= 1 m above MHW).
Scored on the adopted single 5.5 m run.

Results 2026-10-13 (`output/raw_runs/experiments/storms-overwash/2026-10-13-soundside-overwash/`):

| | 1996-2009 | 2009-2025 |
|---|---|---|
| observed new-overwash domain-images | 68 | 157 |
| ocean run-up below the lowest surveyed crest | 16 (24%) | 70 (45%) |
| ... of those, with high sound water in the window | 8 | 36 |
| model hit rate where the ocean could reach the dune | 0.87 (45/52) | 0.34 (30/87) |
| model hit rate where it could not | 0/16 | 0/70 |

Irene's image (2011-08-27; Oregon Inlet 2.02 m) carries 30 of the 70. The
2024-03-14 image carries 26 with no high sound water (the 2023-24 winter
peaked at 0.69-0.79 m in the sound, 0.84 m at Duck): there the likelier cause
is dunes lowered in the 15 years since the 2009 survey. So "implausible" is an
upper bound on sound-side overwash: about 12-24% of calibration events and
23-45% of test events. Separately, the largest 2010-2025 ocean run-up is 4.66 m
NAVD88, below the adopted 5.5 m ceiling: once a model dune reaches the ceiling
no test-period storm can overtop it from the ocean, which is why the model also
misses most of the ocean-plausible test events.

---

## Road relocation

When and where the road was relocated, and how the model reads it. `RELOCATION_COMPARISON_RESULTS.md` is what the set concluded.

Each study drives the unchanged runner (`../../HAT_hindcast.py`); see `../README.md` for how
the experiments are run and filed, and how they import each other across folders.

| Script | What it asks |
|---|---|
| `HAT_relocation_comparison.py` | Does CASCADE relocate NC-12 where and when history did? |
| `HAT_relocation_period_compare.py` | One relocation event, two hindcast windows: does the start year change whether CASCADE reproduces it? |
| `HAT_relocation_dune_position_check.py` | Does the dune line ever reach NC-12? Original, observed and modelled positions, per domain. |

### The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

#### HAT_relocation_comparison.py

Does CASCADE relocate NC-12 where and when history did?

From the script's original header:

```text
Does CASCADE relocate NC-12 where and when history did?

THE QUESTION
    `roadway_manager` relocates the road on its own when the dune line
    overruns it -- `road_relocation_checks` fires the moment the setback goes
    negative. Separately, the pipeline can PRESCRIBE the two historical
    relocations (1989, GIS 84-87; 1999, GIS 9-14) as measured displacements.
    This script runs the two against each other:

        arm A  relocations OFF -- the module decides on its own
        arm B  relocations ON  -- the measured displacements are applied

    and asks whether A reproduces B's timing and footprint unaided.

WHAT IS DELIBERATELY NOT DONE
    Arm A runs the module EXACTLY as built. `_road_relocation_setback` is left
    at its initialised value -- the run's starting setback -- rather than being
    set to a design distance or to the measured post-relocation alignment.
    That matters because six of the ten historical domains (10-13, 85-87)
    start at a setback of 0.0 m, so in those domains a relocation puts the road
    back where it already was and the trigger re-fires the next year the dune
    line moves landward. That ratcheting is a property of the module, and
    reporting it is the point; designing around it would hide it.

WHAT THE COMPARISON CAN AND CANNOT SAY
    CASCADE's trigger is purely geometric: dune line overruns road. There is no
    storm damage, no cost, and no maintenance decision in it. So a match here
    means "the modelled physics would have overrun NC-12 near that year", NOT
    "NCDOT would have moved the road then". The second question is outside what
    this module represents, and no configuration of it gets there.

BACKGROUND EROSION
    Both arms run under whichever source/sink preset the runs were driven with.
    Note that `edgeBE` carries rates on GIS 1 and 90 ONLY, so at every domain
    under test here edgeBE and zeroBE are the same forcing: the dune-line
    retreat that fires the trigger comes entirely from Barrier3D/BRIE dynamics.
    Only `calibBE` puts a background-erosion term on the relocation domains,
    which makes a calibBE re-run the natural sensitivity test once those
    source/sink terms are updated.

TIME INDEXING
    `RoadwayManager` writes every time series at `time_index - 1`, and the
    loop applies a historical event before the update for `start_year +
    time_step`. So index i in `_road_setback_TS` is calendar year
    START_YEAR + i. The prescribed arm is used to CHECK that rather than
    assume it: arm B's setbacks must jump at exactly 1989 and 1999.

THE PERIOD (2026-09-15)
    `--period <start year>` selects a hindcast window from HATTERAS_PERIODS
    (1984 -> 1984-2004, the default; 1996 -> 1996-2010). Only the relocation
    events INSIDE the window are scored: a 1996 start scores the 1999 event
    alone, because the 1989 event is already in its derived setback file. The
    output root follows the period (relocation_<start>_<end>/), the event
    animations are drawn only for events in the window, and the independent
    position cross-check stays at the 2004 measurement -- the period's END for
    1984-2004 and year 8 of 14 for 1996-2010 -- because it is the only
    surveyed road position; the tables say which year it is.

USAGE
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/relocation/HAT_relocation_comparison.py
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/relocation/HAT_relocation_comparison.py --period 1996 --preset edgeBE
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/relocation/HAT_relocation_comparison.py --arm-a DIR --arm-b DIR
```

Notes that were in the code:

```text
Anchored by SEARCHING UPWARD for the project root rather than by
counting parent directories (2026-09-13). A counted depth is correct
only while the file stays where it was written, and these moved into
subfolders of hatteras_ms. Six files here already did it this way.
```

```text
THE PERIOD. Module globals because every scorer, the arm names and the
output root read them; set ONCE by set_period() from --period before any of
that runs. 1984-2004 is the default so the six 2026-09-01 sets and every
regenerate line in RELOCATION_COMPARISON_RESULTS.md still mean what they did.
```

```text
The ONE surveyed road position: RoadOffset_2004_domains.csv, measured on
the 2008 NC-12 line against 2004-start row 0. It is the end of 1984-2004
and the middle of 1996-2010, and the same number in both, which is what
makes the position check comparable across periods.
```

```text
The scenario both arms run. `full_management` is the status-quo hindcast,
and it is the only scenario where a relocation arm is meaningful and the
village management is also present.
```

```text
One directory per preset underneath it. The background-erosion preset is
the axis this comparison is repeated over -- edgeBE and zeroBE differ only
at GIS 1 and 90, while calibBE is the only one carrying a source/sink term
on the relocation domains themselves -- so the artifacts have to be kept
apart or the second run silently overwrites the first.
relocation/<start>_<end>/ (one relocation tree since 2026-09-17; before that
relocation_<start>_<end>/ per window at the top level); re-pointed by
set_period().
```

```text
Tolerance windows for the hit/miss matrix. Two are reported rather than one
because the answer is sensitive to it and a single number would hide that.
```

```text
These animations are READ, not watched. The question they answer -- in which
year does this domain's road step landward, and does the model do it when
history did -- needs the viewer to hold one frame long enough to find the
domain, read the year clock, and compare the two panels. At the shared
default of 3 fps that is 333 ms per frame and the eye cannot do it; a 20-year
run is over in seven seconds.

1 fps gives a second per model year and a ~21 s loop. Slower than a general
shoreline animation wants, which is why this overrides rather than changing
GifConfig's default: the per-run shoreline GIFs show a smooth trend where 3
fps reads fine, and it is only the discrete, dated relocation events that
need dwell time.
```

```text
Alongshore windows for the animation. The two event windows are padded a
few domains beyond the event footprint so the relocating stretch is seen
against road that is NOT relocating -- an unpadded window shows every
domain stepping at once and reads as a global effect.
```

```text
Windows for the topographic raster. The full island is included despite
being 4100 cells wide and a few hundred tall -- rendered wide and short,
that IS the shape of Hatteras, and it is the view that reads as the island
rather than as a chart. The event windows carry the cross-shore detail.
ONE FOLDER PER PLACE inside a set (2026-09-09, Hannah): the two animations
of a window sit together under a readable name, and the eight CSVs under
tables/, so a set reads as report + tables + places instead of fifteen files.
```

```text
Recorded on every topographic frame. The measured island planform spans
6.3 km of cross-shore offset across the real domains, but Cascade's
shoreline_offset reaches BRIE in decameters where BRIE reads metres, so
the run carries a tenth of it. Stated on the figure rather than silently
corrected: the animation shows the island the model actually ran.
```

```text
np.load changes nothing on disk, but Cascade.save() os.chdir()s into the
run directory, so anything downstream that assumes cwd is the repo root
must not rely on it. Absolute paths are used throughout this file.
```

```text
Only meaningful where the trigger never fired. On a domain that
DID relocate the setback is reset by _apply_relocation, so its
minimum is an artifact of the reset rather than a near miss.
```

```text
The modelled position in the check year, per arm, or None if the
manager had stopped by then. For 1984-2004 this is the last year;
for 1996-2010 it is mid-window, and the column says which year.
```

```text
Bound the search by the last MANAGED year, dated from the elevation
series -- see last_managed_index. Unwritten trailing years would
otherwise contribute a spurious jump back to zero.
```

```text
Runs are filed [<forcing arm>/]<period>/<preset>/. Both relocation arms
share a preset by construction -- arm_names derives both from the one
token -- so the directory is resolved once and used for both. Resolved
rather than joined: the join had no slot for the forcing-arm component,
and "arm" here means the relocation switch, not that one.
```

```text
A STALE report is worse than a missing one -- see _report_header. Delete
any previous report BEFORE doing the work, so a run that dies partway
leaves this folder with no report rather than with the previous run's.
```

<details><summary>Function notes (the original docstrings)</summary>

**`set_period()`**

```text
Points the module at one hindcast window. Raises on a year that is not
a HATTERAS_PERIODS key, so a typo cannot score an empty window.
```

**`arm_names()`**

```text
The two run directory names this comparison reads, for one preset.

Built with the same token rule the hindcast derives RUN_NAME with in its
section 7.5: the arms are identical in every token except `reloc`. Both
are derived from ONE preset argument on purpose -- an arm A read against
an arm B from a different preset would produce a clean-looking comparison
of two different forcings, and nothing downstream would notice, because
every check in this file is a difference between the arms.

Args:
    preset: A source/sink preset key or deprecated alias.

Returns:
    (arm_a_name, arm_b_name), relocations off and on.
```

**`windows_in_period()`**

```text
The windows worth drawing for this period: the island, plus each event
window whose event year is one the period actually scores. A 1996 start
has no 1989 event to show; drawing that block would animate two identical
arms and read as a null result.
```

**`load_shoreline_matrix()`**

```text
Loads the plan-view shoreline matrix a run saved beside its model.

Args:
    run_dir: The run directory.

Returns:
    A 2-D [n_years, total_domains] array in metres, raw x_s_TS
    convention, or None if the run did not save one.
```

**`back_barrier_matrix()`**

```text
Back-barrier shoreline per domain per year, metres, raw convention.

The shoreline matrix the run saves carries only x_s. Drawing the island
with any width needs x_b as well, on the same sign convention so the two
can share an axis.

Args:
    cascade: A finished Cascade instance.

Returns:
    A [n_years, n_domains] array in metres, or None if x_b_TS is absent.
```

**`load_cascade()`**

```text
Loads the pickled Cascade a run wrote with `cascade.save(run_dir)`.

Args:
    run_dir: Directory holding exactly one .npz model state.

Returns:
    The Cascade instance.

Raises:
    FileNotFoundError: If the directory holds no .npz, which means the run
        was driven with HAT_SAVE_MODEL_STATE=false and cannot be compared.
```

**`road_series()`**

```text
Pulls the per-domain roadway time series out of a finished run.

Reads the managers CASCADE already holds rather than re-deriving anything.
Only domains the run actually managed are returned: a domain outside
`roadway_management_module` has a RoadwayManager object that was never
called, and its all-zero series would read as "never relocated" rather
than "never asked".

Args:
    cascade: A Cascade instance after its run.
    geometry: DomainGeometry describing the padded array.
    first_gis: First GIS domain carrying road.
    last_gis: Last GIS domain carrying road.

Returns:
    A {gis: dict} mapping, each dict holding setback, relocated, elevation
    (numpy arrays indexed by years-since-START_YEAR) plus the drowned and
    relocation_blocked flags.
```

**`historical_targets()`**

```text
Maps each domain a relocation event moves to that event's year.

Read from HATTERAS_ROAD_EVENTS rather than retyped, so dropping a domain
from an event (as GIS 15 was) drops it from the scoring too.

Args:
    start_year: First calendar year of the period.
    end_year: Last calendar year of the period.

Returns:
    A {gis: event_year} dict for relocation events inside the period.
```

**`first_relocation_year()`**

```text
Calendar year of the first modelled relocation, or None.

Args:
    relocated_ts: `_road_relocated_TS`, indexed by years since start_year.
    start_year: First calendar year of the run.

Returns:
    The calendar year, or None if the domain never relocated.
```

**`relocation_margin()`**

```text
How close a road came to firing the relocation trigger, and never did.

THE TRIGGER, EXACTLY. `roadway_manager.road_relocation_checks` does

    road_setback = road_setback + dune_migrated      # dune_migrated < 0 landward
    if road_setback < 0:  relocate

and the caller supplies

    dune_migration = barrier3d.ShorelineChangeTS[t-1] * 10       # m

`ShorelineChangeTS` counts WHOLE dam cells, so the setback only ever moves
in 10 m steps and the test is STRICT. Both facts matter for a margin:

  * a setback sitting at exactly 0.0 has NOT fired. The dune line has
    reached the road and stopped there. It needs one more full cell.
  * so the extra landward migration needed is `min_setback + 10 m`, not
    `min_setback`. Reporting the setback alone would say GIS 84 needed 0 m
    more, which is wrong -- it needed one more cell.

Args:
    entry: One road_series() value, with its "setback" array.
    start_year: First calendar year of the run.

Returns:
    dict with the closest approach, the year it happened, the cells still
    between dune and road there, and the extra migration that would fire.
```

**`near_miss_table()`**

```text
Every managed domain that never relocated, ranked by how close it came.

Covers the CONTROL domains as well as the historical ones on purpose. The
comparison's headline false-positive count is 0/45, which reads as "the
module is appropriately conservative" -- but a control domain sitting one
cell from firing is a different statement about robustness than one sitting
thirty cells away, and only this table separates them.

Args:
    series: road_series() output for the free-running arm.
    targets: {gis: historical_event_year}.
    start_year: First calendar year of the run.

Returns:
    DataFrame of non-relocating domains, closest first.
```

**`score_first_year()`**

```text
Builds the first-relocation-year table -- the primary result.

Robust to the ratcheting described in the module docstring: however many
times a domain relocates afterwards, the FIRST firing is the model's
answer to "when did the dune line reach the road".

Args:
    series: road_series() output for the free-running arm.
    targets: {gis: historical_event_year}.
    start_year: First calendar year of the run.

Returns:
    A DataFrame, one row per historical domain.
```

**`score_confusion()`**

```text
Hit/miss matrix at one tolerance window.

A HIT is a historical domain that relocated within +/-tolerance years of
its event. A FALSE POSITIVE is a managed road domain history never
relocated that relocated anyway, at any time. The false-positive count is
the half of this that a per-domain error table cannot show: a module that
relocates everywhere scores perfectly on the historical domains while
saying nothing.

Args:
    series: road_series() output for the free-running arm.
    targets: {gis: historical_event_year}.
    start_year: First calendar year of the run.
    tolerance: Half-width of the match window, in years.

Returns:
    A dict of counts and rates.
```

**`last_managed_index()`**

```text
Index of the last year the manager actually ran for this domain.

Dated from `_road_ele_TS`, NOT from the setback series. A setback of
exactly 0.0 m is legitimate here -- it means the road sits on the dune
line, which is where six of the ten historical domains start -- so a
nonzero test on the setback would report those domains as never managed.
Road ELEVATION has no such ambiguity: the module stops managing the moment
it drops below 0 m MHW, so its last non-zero entry dates the last managed
year. This is the same signal `summarise_road_management` dates from.

Args:
    entry: One road_series() value, carrying "elevation".

Returns:
    The index, or -1 if the manager never ran.
```

**`score_trajectories()`**

```text
Per-domain setback trajectory comparison between the two arms.

Compared over the window BOTH arms were still managing the domain. If one
arm's road drowns, its series stops being written; extending the
comparison past that point would score an unwritten zero against a live
setback and report a difference that is really the end of the record.

Args:
    series_a: road_series() for the free-running arm.
    series_b: road_series() for the prescribed arm.
    targets: {gis: historical_event_year}.
    start_year: First calendar year of the run.
    check_2004: {gis: measured_setback_m}, the independent cross-check.

Returns:
    A (summary_df, long_df) tuple. long_df is one row per domain-year and
    is what the GIF and any trajectory plot read; it carries a `managed`
    flag per arm so a plot can stop the line where the record stops.
```

**`check_determinism()`**

```text
Confirms the two arms are identical before the first prescribed event.

Both arms are the same configuration up to the first event, and Barrier3D
runs on a seeded RNG with a prescribed storm file, so every series must
agree exactly through `first_event_year - 1`. A divergence there is a bug,
not a result, and everything downstream would be uninterpretable -- so
this is checked rather than assumed.

Args:
    series_a: road_series() for the free-running arm.
    series_b: road_series() for the prescribed arm.
    first_event_year: Calendar year of the earliest prescribed event.
    start_year: First calendar year of the run.

Returns:
    A (ok, offenders) tuple; offenders lists the GIS domains that diverged.
```

**`check_event_indexing()`**

```text
Confirms index i really is calendar year start_year + i.

Arm B's setback is displaced by a measured amount at the event year and by
dune migration alone in every other year, so the largest single-year jump
in the prescribed arm must land on the event year. If it does not, the
time indexing in this file is wrong and every year reported here is off.

Args:
    series_b: road_series() for the prescribed arm.
    targets: {gis: historical_event_year}.
    start_year: First calendar year of the run.

Returns:
    A DataFrame with the located jump year beside the expected one.
```

**`_Tee()`**

```text
Mirrors everything printed to the console into a buffer.

The report is the run's own console output rather than a second rendering
of the same numbers, which is the point: a separately-composed summary can
disagree with the CSVs beside it, and this cannot.
```

**`default_out_dir()`**

```text
OUTPUT_ROOT/<topo version>/<preset>[_groin] (2026-09-09, Hannah: the
folder must say which dune-topo version a comparison was made on).

The version is READ from arm A's run metadata, never typed, so a set cannot
land unlabelled; the two arms must agree on it. The groin token follows
the arm names: a `_groin` arm (not `nogroin`) gets a `_groin` folder, the
convention the six 2026-09-01 sets used.
```

**`_report_header()`**

```text
Provenance block written above the captured output.

WHY THIS EXISTS. On 2026-08-25 a re-run of this comparison rewrote every
CSV and GIF in output/comparisons/relocation_1984_2004/<preset>/ (deleted 2026-10-23) and left
the report.txt from 2026-08-22 sitting beside them -- the script had lost
its report-writing step, so nothing overwrote it. For three days that
folder held a report describing DIFFERENT runs from the CSVs next to it,
with nothing on its face to say so. It even named the pre-restructure flat
run paths, which by then did not exist.

So every report now carries the identity of the two runs it was built
from. A report whose arms do not match the runs on disk is visible at a
glance instead of having to be inferred from file mtimes.
```

</details>

#### HAT_relocation_period_compare.py

One relocation event, two hindcast windows: does the start year change whether CASCADE reproduces it?

From the script's original header:

```text
One relocation event, two hindcast windows: does the start year change
whether CASCADE reproduces it?

THE QUESTION
    The 1999 NC-12 relocation (GIS 9-14) sits inside BOTH the 1984-2004 and
    the 1996-2010 hindcast windows. Each window has its own emergent-vs-
    prescribed comparison (HAT_relocation_comparison.py, one set per preset),
    scored on its own terms. This script reads those two sets side by side
    for the domains ONE event moved and asks what the start year changed:

      * how much dune retreat each window accumulates at the road before
        the event year -- 15 model years from a 1984 start, 3 from 1996;
      * whether the free-running arm fires at all, and when, in each;
      * how each window's modelled 2004 road position compares with the one
        surveyed position, which is the END of one window and the MIDDLE of
        the other.

    It re-scores nothing. Every number is read from the per-period tables,
    so a disagreement between this report and a per-period report is a
    stale set, not a second opinion.

WHAT THE TWO WINDOWS SHARE, AND WHAT THEY DO NOT
    Same topography product (1984-start, one dune-topo version, read from the
    sets and required to agree), same road line (1978 for both starts), same
    code, same relocation target. They differ in the start year, the storm
    series, the offset survey (1984 line vs the 1997 line), and the setback
    file: a 1996 start reads the 1984 setbacks with the 1989 event already
    applied, which does not touch GIS 9-14. So at the event domains the two
    windows start the road in the SAME place and differ only in how many
    years of modelled retreat precede 1999.

USAGE
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/relocation/HAT_relocation_period_compare.py
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/relocation/HAT_relocation_period_compare.py --presets zeroBE edgeBE --version v2
```

Notes that were in the code:

```text
A domain the free arm relocated BEFORE the event has had its
setback reset to the relocation target, so its pre-event
setback no longer measures retreat. Reported as NaN, and the
first-year table says when it fired.
```

<details><summary>Function notes (the original docstrings)</summary>

**`domain_table()`**

```text
One row per (period, domain): retreat before the event, the free arm's
answer, and the position check. Everything read, nothing re-scored.
```

**`event_recall()`**

```text
Hits among the event domains, from the FIRST modelled relocation year.
The per-period confusion.csv counts every event in the window and any
relocation year; this restricts to one event and uses the first firing,
which is the model's answer to 'when did the dune reach the road'.
```

**`trajectory_figure()`**

```text
One panel per event domain. Each window is a colour (the vintage pair:
the earlier start in the RdBu red, the later in the blue); the free arm
is solid, the prescribed arm dashed; the event year is a vertical rule;
the surveyed 2004 position is a black marker. A line stops where that
arm stopped managing the road.
```

</details>

#### HAT_relocation_dune_position_check.py

Does the dune line ever reach NC-12? Original, observed and modelled positions, per domain.

From the script's original header:

```text
Does the dune line ever reach NC-12? Original vs observed vs modelled.

WHY THIS FIGURE EXISTS

    The 1984-2004 relocation comparison reports that almost no domain relocates
    unaided once the road setbacks are measured against `1984-start` row 0. The
    obvious objection is that the roads are still only 3-6 cells behind the
    dune, so something must be wrong. This figure is the check: it puts the
    three cross-shore positions on one axis, per domain, so the claim can be
    read off rather than argued.

        original   the 1984 dune line -- the run's own year-0 position, and the
                   datum every other quantity here is measured from
        observed   the surveyed 2004 dune line, from
                   2-brie-offset/raw_offsets/2004_duneline_offset_raw.csv minus
                   the 1984 file. This is the SAME target the run scores its
                   misfit against (cascade_pipeline.hindcast.
                   build_shoreline_target), not a second opinion.
        modelled   where the run put the dune line in 2004

    and draws NC-12 as the 20 m band it occupies, at the setback the model was
    initialised with.

WHAT "CORRECT BEHAVIOUR" LOOKS LIKE HERE

    `road_relocation_checks` relocates when the setback goes STRICTLY negative,
    and the setback is driven by

        dune_migration = barrier3d.ShorelineChangeTS[t-1] * 10       # m

    -- whole 10 m cells, because ShorelineChangeTS counts cells. So the road is
    overrun only when the dune line travels PAST the near edge of the road band.
    A domain whose observed and modelled 2004 dune lines both stop short of that
    edge SHOULD not relocate, and a model that agrees with the survey about
    where the dune line got to is behaving correctly even though it misses the
    historical relocation. That is the distinction this figure is drawn to make
    visible: a miss caused by the TRIGGER (geometric overrun) is not the same as
    a miss caused by the PHYSICS (dune line in the wrong place).

SIGN CONVENTION
    +x is LANDWARD throughout, matching x_s_TS and the raw offset files. The
    1984 dune line is 0 by construction on every domain.

USAGE
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/relocation/HAT_relocation_dune_position_check.py
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/relocation/HAT_relocation_dune_position_check.py --preset calibBE
```

Notes that were in the code:

```text
HOUSE STYLE: one typeface and one palette across every figure in this
project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
found by searching upward (ORGANIZATION.md rule 5). This file drew in
matplotlib's defaults until 2026-09-17 -- it never called apply_style().
```

```text
Anchored by SEARCHING UPWARD for the project root rather than by
counting parent directories (2026-09-13). A counted depth is correct
only while the file stays where it was written, and these moved into
subfolders of hatteras_ms. Six files here already did it this way.
```

```text
Palette. Deliberately not a rainbow: the three positions are one family
(where is the dune line) and the road is the thing they are compared against,
so the road is the only warm colour on the figure.
```

```text
Observed: difference two ABSOLUTE surveyed distances. The padded offset
files each subtract their own year's minimum, so differencing THOSE is
not a shoreline change -- this is the same call build_shoreline_target
makes for the run's own misfit line.
```

```text
`margin_m` mirrors HAT_relocation_comparison.relocation_margin():
the trigger is `setback < 0` STRICT and the setback moves in whole
10 m cells, so the extra migration needed is the closest approach
plus one cell. Meaningless where the road did relocate -- the reset
in _apply_relocation puts an artificial minimum in the series.
```

```text
---- panel A: every managed domain -------------------------------------
The road band is drawn as a bar from its near (seaward) edge landward,
because that near edge is the thing the dune line has to cross.
```

```text
Both labels go ABOVE the road bar. An earlier version put the relocation
count below the lowest marker, where it collided with the tick labels on
exactly the two domains that relocate.

For a domain that never fired, the useful number is not "it did not
relocate" but HOW CLOSE it came: the extra landward dune migration that
would have fired the trigger. The setback moves in whole 10 m cells and
the test is `< 0` strict, so a road sitting at setback 0 still needs one
more full cell -- hence min_setback + 10, not min_setback.
```

```text
under the dune-topo version the run was made on (2026-09-09), read from
the run's metadata; the layout is OUT_ROOT/<version>/dune_position_check/
```

```text
---- the decomposition this figure was drawn to produce ----------------
Splitting the misses by CAUSE is the whole point. An island-wide misfit
near zero hides it: the model tracks the survey on average and still
under-predicts badly on exactly the domains under test.
```

---

## Storms

How the storm series is built: maximum duration, trim length, event splitting, height and run-up.

Each study drives the unchanged runner (`../../HAT_hindcast.py`); see `../README.md` for how
the experiments are run and filed, and how they import each other across folders.

| Script | What it asks |
|---|---|
| `HAT_storm_max_duration.py` | Why does a storm series with longer events drown the barrier? |
| `HAT_storm_length_selection.py` | Which storm duration rule should the hindcast use? |
| `HAT_storm_event_splitting.py` | Should back-to-back storms be separate events? |
| `HAT_storm_height_test.py` | Are the storm water levels too low for Hatteras? |
| `HAT_storm_rlow_duration.py` | Do Rlow and duration, defined as Barrier3D expects them, improve the model? |
| `HAT_storm_rlow_duration_plot.py` | How the equivalent duration is set, drawn for Isabel 2003 and the March 2018 nor'easter. |
| `HAT_storm_runup_by_slope.py` | Does storm run-up from each domain's own beach slope improve the overwash skill? |
| `HAT_storm_notebook_rule.py` | Does running the storms by the notebook's own rule (no split, no trim) change the model's skill? |
| `HAT_trim_length_adopted.py` | Does the storm trim length still matter once the dunes are realistic? |

### The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

#### HAT_storm_max_duration.py

Why does a storm series with longer events drown the barrier?

From the script's original header:

```text
HAT_storm_max_duration.py -- why does a storm series with longer events drown the barrier?
THE QUESTION (Hannah, 2026-09-28). The storm series keeps events of 8-72 h.
72 was chosen because runs on 96, 120 and 240 h series "did not work": the
barrier drowned and the simulation ended (PROVENANCE.md in 3-storms/). But the
builder DROPS a longer event rather than shortening it, so the 72 h series is
missing 29 events 1996-2024, among them Isabel 2003 (Rhigh 5.22 m MHW), the
March 2018 nor'easter, Florence 2018 and Dennis 1999. Which domain drowns on a
longer series, when, and by what mechanism?

NOTHING IN THE MAIN CODE CHANGES
    The variants are built by the builder's own functions, read out of
    historical_storm_creation_v3_HAT.py by `ast` (the script runs at import,
    so it cannot be imported), with only max_storm_dur and the save location
    changed. They are written HERE, never into hindcast_storms/. The 72 h
    variant is rebuilt too and must equal the committed series exactly.

    A run is the unchanged hindcast runner, executed by this file's
    `_launch` action, which points HATTERAS_PERIODS[start]["storm_file"] at the
    variant in its own process first. Barrier3D is the current one
    (49fd069), as in every matrix run.

VARIANTS  (name -> max hours; "trim" keeps a longer event, cut to the 72 h
around its peak, which the builder cannot do)
    72     the committed series (check only; the matrix runs are its controls)
    96, 120, 240, nocap
    72trim

WHERE: output/raw_runs/experiments/storms-overwash-and-dunes/2026-09-28-storm-max-duration/
           storms/<window>/<window>_storms_v3_<variant>.npy (+ _summary.csv)
           runs/<variant>_<scenario>/<period>/edgeBE/<run_name>/
           logs/, tables/, figures/, NOTE.md

USAGE
    python HAT_storm_max_duration.py build
    python HAT_storm_max_duration.py run [--workers 4]
    python HAT_storm_max_duration.py diagnose
```

Notes that were in the code:

```text
THE CAUSE TEST. Barrier3D before 49fd069, the route_overwash axis-swap fix
of 2026-09-24, read Elevation[TS, i, d+1:d+10] (out of bounds whenever
i >= rows). Longer storms route for more steps. A detached worktree at the
commit before the fix (ce36866) runs the 72 h and 240 h series.
```

<details><summary>Function notes (the original docstrings)</summary>

**`builder_functions()`**

```text
load_data, calculate_r2_percent and create_storms, compiled from the
builder's source without running its module-level code.
```

**`trim_long_events()`**

```text
The nocap events, each longer than `limit` cut to the `limit` hours
above the berm centred on its peak TWL; Rhigh, Rlow, period and duration
recomputed on what is kept, exactly as the builder computes them.
```

**`_launch()`**

```text
Subprocess entry: the unchanged runner, with this period's storm file
pointed at a variant IN THIS PROCESS ONLY.
```

**`diagnose()`**

```text
For every run: did it stop, which domain drowned, in which model year,
by width or by height, and the storms of that year in both series.
```

**`compare()`**

```text
Each variant run against its matrix control (the committed 72 h series):
net shoreline change, overwash, skill against CoastSat, and the observed
overwash hit rate with each run dated by ITS OWN storm file.
```

</details>

#### HAT_storm_length_selection.py

Which storm duration rule should the hindcast use?

From the script's original header:

```text
HAT_storm_length_selection.py -- which storm duration rule should the hindcast use?
THE QUESTION (Hannah, 2026-09-28): pick the storm series, with overwash that
matches the observed record as the first priority.

Why it is open: the committed series (v3_72) DROPS every grouped event
longer than 72 h, which removes Isabel 2003, March 2018, Florence, Dennis
and Nor'Ida. The 72 h limit existed only because the pre-49fd069 Barrier3D
crashed on long storms (storms-overwash-and-dunes/2026-09-28-storm-max-duration).
Barrier3D also routes each storm at its PEAK Rhigh for its WHOLE duration, so
an event's length is a lever on how much overwash it makes, and the builder's
24 h grouping makes some events of 150-190 h out of strings of smaller surges.

THE CANDIDATES (all keep every event; "trimL" cuts an event longer than L to
the L hours around its peak TWL, recomputing Rhigh/Rlow/period on what is kept)
    drop72   the committed series: the control (the matrix runs)
    trim24, trim36, trim48, trim72, trim96, trim120, trim168
    full     no limit (= a 240 h limit: the longest event 1996-2024 is 193 h)

STAGE 1 -- overwash, against the imagery (8-overwash-analysis)
    Each candidate on both windows, full_management (the imagery is the managed
    island; also the scenario the ends are solved on) and natural. Edge rates
    as the matrix (solved on drop72). Overwash barely depends on the two end
    rates, so this is a fair screen.

    WHICH STORM MADE THE OVERWASH. Barrier3D records overwash per domain per
    YEAR. Every storm of a year is tested against one crest, the dune after
    the year's growth; that crest is recomputed exactly (Barrier3d.DuneGrowth
    on the saved dunes, which SeaLevel has already lowered in place), each
    storm's gaps come from the model's own DuneGaps, and the year's overwash
    is shared among the storms that reach a gap in proportion to the water
    they put through it (gap width x Qdune(Rexcess) x hours, Qdune as in
    Barrier3D). A storm's share is dated by its end less the observed record's
    7-day grace. `validate` checks the sharing against the storm replay.

    SCORES, per image x domain cell the imagery assessed, model overwash =
    shared overwash in the window since the previous image > THRESHOLD:
        POD    hit rate: observed overwash the model reproduces
        POFD   false-alarm rate: observed-absent cells the model overwashes
        PSS    Peirce skill score = POD - POFD (base-rate free; a series that
               overwashes everywhere scores ~0, not 1)
        timing r  per-image domain counts, observed vs model
        space  r  per-domain share of images with overwash, observed vs model
    "Model only" is not all error (washover fades; NC-12 is cleared), which is
    why PSS, not accuracy, is the headline.

STAGE 2 -- shoreline, fair (`ends`, then `score`)
    For the finalists and drop72: the two end rates re-solved on the
    candidate (full_management, Newton on LRR at GIS 1 and 90, the protocol
    of be_edge_domain_solve.py), then natural and managed runs at those ends:
    interior RMSE and bias against CoastSat, and overwash scored again.

NOTHING IN THE MAIN CODE CHANGES: series are built by the builder's own
functions into this folder; runs are the unchanged runner with the period's
storm file swapped in its own process (HAT_storm_max_duration._launch);
Barrier3D is the current one (49fd069).

WHERE: output/raw_runs/experiments/storms-overwash-and-dunes/2026-09-28-storm-length-selection/

USAGE
    python HAT_storm_length_selection.py build
    python HAT_storm_length_selection.py run [--workers 4]
    python HAT_storm_length_selection.py validate
    python HAT_storm_length_selection.py score
    python HAT_storm_length_selection.py ends --variants drop72 trimXX ...
    python HAT_storm_length_selection.py figures
```

<details><summary>Function notes (the original docstrings)</summary>

**`storm_shares()`**

```text
{(pad, storm row index): overwash m3/m} -- each domain-year's QowTS
shared among that year's storms that reached a dune gap.
```

**`overwash_cells()`**

```text
One row per assessed image x domain: observed, and model overwash
(shared m3/m in the window since the previous image).
```

**`validate()`**

```text
Share each storm's overwash as score() does, then replay the same
domain-years through Barrier3D (natural runs, where the replay is exact)
and compare per-storm volumes: does the sharing put the overwash on the
right storms, i.e. in the right image window?
```

**`solve_ends()`**

```text
Newton on LRR at GIS 1 and 90 (full_management), from the matrix ends,
each end on its own secant, as be_edge_domain_solve.py prints them.
```

</details>

#### HAT_storm_event_splitting.py

Should back-to-back storms be separate events?

From the script's original header:

```text
HAT_storm_event_splitting.py -- should back-to-back storms be separate events?
WHY (Hannah, 2026-09-29: "test splitting the events"). The storm builder joins
above-berm spells less than 24 h apart into one event (weather_grouping = 24).
At Hatteras the berm is overtopped at most high tides during an active spell,
so two storms a week apart chain into one long event: Edouard + Fran 1996
(one event, 119 h above the berm) and Jose + Maria 2017 (186 h). The adopted
series (v3_trim24) then keeps the 24 h around the event's single highest peak,
so Fran and Jose are not in the model at all, and 1996-2024 loses 56 spells of
>= 8 h above the berm inside merged events (41 peak above 2 m MHW).

THE CANDIDATES (every one trims each event to 24 h, as adopted)
    trim24   the adopted series (hindcast_storms/*_v3_trim24): the control.
             Rebuilt here from the builder's functions and checked identical.
    g12      the builder with weather_grouping = 12 h: the one-number change.
             It recovers Fran and Jose, but a merged storm's tidal fragments
             shorter than 8 h then fall under the minimum-duration rule and
             are dropped, so it has FEWER events and storm-hours than trim24.
    split12  the 24 h grouping kept to define a weather system, then the system
             split wherever the water stays below the berm for >= 12 h; a
             piece shorter than 8 h is folded into the piece before it (or
             after, for the first), so no hour the adopted series counts is
             lost. Each piece is then trimmed to 24 h and dated by its start.

RUNS: managed (full_management), both windows, edgeBE, the site config's end
rates (not re-solved, as in the trim-length check), the unchanged runner with
the period's storm file swapped in its own process (HAT_storm_max_duration).
SCORES: as the trim-length check -- overwash against the imagery (POD, POFD,
PSS, timing and space r), interior RMSE/bias against LOWESS-7, total overwash.

NOTHING IN THE MAIN CODE CHANGES.

    python HAT_storm_event_splitting.py build   # series + what each recovers
    python HAT_storm_event_splitting.py run     # 6 runs, 5 at a time
    python HAT_storm_event_splitting.py score

WHERE: output/raw_runs/experiments/storms-overwash-and-dunes/2026-09-29-event-splitting/
```

Notes that were in the code:

```text
The runner records the storm file relative to data/hatteras_init and
fails at the end of the run on a path outside it (2026-09-29: four
runs lost their metadata and shoreline matrix that way). A relative
path with ".." resolves to the same file and satisfies it.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_pieces()`**

```text
Split one system's above-berm hours (a sorted DatetimeIndex) at gaps
>= gap_h; fold a piece shorter than min_h into its neighbour.
```

**`split_series()`**

```text
systems: the untrimmed 24 h-grouped events (the builder's 'full' run).
gap_h None = no split (reproduces trim24).
```

</details>

#### HAT_storm_height_test.py

Are the storm water levels too low for Hatteras?

Note (2026-09-30): the `slope` action in the original header below no longer exists; the slope0p10 series is built by `build`, and `main()` accepts `gauges`, `build`, `run` and `score`.

From the script's original header:

```text
HAT_storm_height_test.py -- are the storm water levels too low for Hatteras?
WHY (storms-overwash-and-dunes/2026-09-28-dune-ceiling-per-domain): once the dunes
match the 2009 lidar, the low spots overwash where Irene really did. But
Irene also overwashed 47 of the 60 domains with dunes of 4.3 m and up, and the
model manages 10: a modelled Irene Rhigh of 3.80 m MHW cannot clear them.
The storm series takes water level from the Duck gauge (8651370), about 80 km
north of the reach, and adds Stockdon (2006) R2% run-up on WIS ST63228 waves
with one beach slope, 0.06.

PART 1 -- observations (no model runs)
    gauges   peak water level (m above MHW) at Duck against the gauges in or
             near the reach, for the named storms 1996-2024. The ocean-side
             record inside the reach is the Cape Hatteras Fishing Pier
             (8654400, historic); Oregon Inlet Marina (8652587) and USCG
             Station Hatteras (8654467) sit inside the inlets, on the sound
             side. Fetched from the NOAA CO-OPS API into this folder's data/.
    slope    the foreshore slope the run-up should use, measured from the
             domains' 10 m elevation profiles (first land cell to the dune
             toe), against the 0.06 in the storm builder
PART 2 -- sensitivity runs
    storm variants of the trim24 series: run-up slope 0.08 and 0.10
    (Stockdon recomputed, events re-found), and Rhigh/Rlow raised by 0.25 and
    0.5 m (a local surge Duck does not see; the events are unchanged)
    x dune ceilings uniform 5.5 m NAVD88 and per-cell (the two candidates)
    x managed, both windows. Controls: the trim24 runs of those two ceilings.

NOTHING IN THE MAIN CODE CHANGES (storm files and ceilings are swapped in
each run's own process, as in the earlier experiments).

WHERE: output/raw_runs/experiments/storms-overwash-and-dunes/2026-09-28-storm-height/

USAGE
    python HAT_storm_height_test.py gauges
    python HAT_storm_height_test.py slope
    python HAT_storm_height_test.py build
    python HAT_storm_height_test.py run [--workers 6]
    python HAT_storm_height_test.py score
```

<details><summary>Function notes (the original docstrings)</summary>

**`fetch()`**

```text
Hourly heights (historical stations: 'hourly_height'; otherwise
'water_level'), m above MHW. Cached under data/.
```

**`build()`**

```text
slope0p10: the builder's chain with run-up slope 0.10, events re-found
(every event kept), trimmed to 24 h. plusX: the trim24 series with Rhigh
and Rlow raised X m (the same events).
```

</details>

#### HAT_storm_rlow_duration.py

Do Rlow and duration, defined as Barrier3D expects them, improve the model?

```text
WHY (Hannah, 2026-10-01: "set up the Rlow + duration experiment", after a
review of the storm builder). Two storm columns do not mean what Barrier3D
assumes they mean.
  Rlow      the builder takes min(TWL) over the event. An event is TWL above
            the berm, so the median Rlow is 1.36 m MHW, the berm itself.
            Barrier3D uses Rlow only to choose the overwash regime
            (inundation when Rlow > a gap's mean crest, barrier3d.py:1520).
            Upstream MSSM (UNC-CECL/MSSM_RENCI, multivariateSeaStorm_NCB.m)
            sets hourly Rlow = TWL - S/2, Stockdon swash S, and the storm's
            Rlow is its maximum.
  duration  Barrier3D holds the peak Rexcess for every hour of `duration`
            (Discharge[:, 0, gap] = Qdune for the whole storm). Upstream
            yearly_storms halves the duration for that reason; the adopted
            series trims each event to the 24 h around its peak instead.

THE SERIES (all split12: 24 h grouping, split at >= 12 h below the berm)
    split12     the adopted series (v3_split12_trim24): the control. Rebuilt
                here and checked identical in both windows.
    rlow        Rlow = max(TWL - S/2) over the event; trim 24 h as adopted.
    dureq       Rlow as the control; duration = the hours at the peak level
                that pass the same Rexcess^1.5 flow over the berm as the
                whole event, sum((TWL - berm)^1.5) / (Rhigh - berm)^1.5,
                rounded, >= 1. No 24 h cap: the event is not trimmed.
    rlow_dureq  both.
The berm is the lowest crest a gap can have, so dureq is the longest the
equivalent duration can be; a higher crest makes the flow more peaked.

RUNS: full_management and natural, both windows, edgeBE, the site config's
end rates (not re-solved), the unchanged runner with the storm file swapped in
its own process (HAT_storm_max_duration._launch). 16 runs.
SCORES: as the event-splitting check (overwash vs imagery: POD, POFD, PSS,
timing and space r; interior RMSE/bias vs LOWESS-7; total overwash), plus the
share of storm-gaps Barrier3D routed as inundation (_InundationCount).

NOTHING IN THE MAIN CODE CHANGES.

    python HAT_storm_rlow_duration.py build   # 3 series + tables/series.csv
    python HAT_storm_rlow_duration.py run     # 16 runs, 5 at a time
    python HAT_storm_rlow_duration.py score   # tables/scores.csv

WHERE: output/raw_runs/experiments/storms-overwash-and-dunes/2026-10-01-rlow-and-duration/
```

#### HAT_storm_runup_by_slope.py

Does storm run-up from each domain's own beach slope improve the overwash skill?

```text
WHY (Hannah, 2026-10-01: "test the spatially varying runup with beach slope").
Every domain sees the same storm: one gauge, one WIS point, one slope (0.06).
After the storm input, duration and dune recovery were ruled out, what is left
is the height the water reaches at each place against the dunes there.

SLOPES: per domain on the 2009 lidar (0-elevation/2009-2014/1-gapfill-1m), along
each east-west row from the ocean-side MHW crossing (0.36 m NAVD88) to the first
cell at the berm (1.7 m), corrected for the shoreline's angle to the rows,
median over rows. Both windows read the 2009 lidar: the 1996 graft
(2009-2014-1996) left 19 domains unmeasurable and gave slopes to 0.72. Clipped
to Stockdon's range (0.01-0.16), then a 5-domain running median.

THE SERIES (same events and hours as the adopted v3_split12_trim24; only
Rhigh and Rlow change, recomputed hour by hour with the domain's slope)
    control   0.06 everywhere: checked identical to the adopted file, all 90
    measured  each domain's slope (median 0.075, 0.05-0.13)
    pattern   the measured slopes rescaled to median 0.06: the spatial pattern
              without the shift in the mean
Each domain's series is set on its Barrier3D object after build_cascade, in the
run's own process (CASCADE itself takes one storm file for all domains).

RUNS: full_management and natural, both windows, edgeBE, site-config ends. 12.
SCORES: hit rate, false-alarm rate, skill, false alarms and misses by zone,
alongshore r of overwash frequency, interior RMSE/bias vs LOWESS-7.

NOTHING IN THE MAIN CODE CHANGES.

    python HAT_storm_runup_by_slope.py slopes | build | run | score

WHERE: output/raw_runs/experiments/storms-overwash-and-dunes/2026-10-01-runup-by-slope/
```

#### HAT_trim_length_adopted.py

Does the storm trim length still matter once the dunes are realistic?

From the script's original header:

```text
HAT_trim_length_adopted.py -- does the storm trim length still matter once the dunes are realistic?
WHY (Hannah, 2026-09-28: "run the trim-length check first"). trim24 was chosen
in storms-overwash-and-dunes/2026-09-28-storm-length-selection on the OLD dunes (held
near 3 m MHW by the default Dmaxel), and 24 h was the shortest length tried.
Barrier3D applies a storm's peak Rhigh for its whole duration, so the length is
a lever on how much sand overwash moves. This re-asks the question on the
setup being adopted:
    Barrier3D  branch hatteras/adopted (worktree ../Barrier3D-adopted): the
               three overwash fixes + per-cell dune ceilings, switched on
               (DuneCeilingFromStart true, floor 0.5 m) in each run's process
    storms     every event kept, trimmed to 12, 24, 48, 72 h, or full length
               (12: the builder, --long-events trim --max-duration 12, written
               here; 24: the adopted hindcast_storms v3_trim24 files; 48, 72,
               full: the storm-length selection's verified files)
    runs       managed (full_management), both windows, the site config's
               current end rates (the LOWESS-7 solve)
SCORES: overwash against the imagery (as before), the 2010 dune crest against
the lidar, interior RMSE/bias against LOWESS-7 (run_registry.skill_vs_target).
Scoring runs under the same Barrier3D, so the storm sharing uses the fixed
DuneGaps and the per-cell DuneGrowth.

NOTHING IN THE MAIN CODE CHANGES.

WHERE: output/raw_runs/experiments/storms-overwash-and-dunes/2026-09-28-trim-length-adopted/
```

Notes that were in the code:

```text
The worktree these runs used was removed on 2026-09-28 when hatteras/adopted
was checked out in ../Barrier3D itself (the editable install); same branch.
```

---

## Model windows

Which start and end years make good calibration and test windows.

Each study drives the unchanged runner (`../../HAT_hindcast.py`); see `../README.md` for how
the experiments are run and filed, and how they import each other across folders.

| Script | What it asks |
|---|---|
| `HAT_candidate_windows.py` | What does the model give on the candidate windows 1996-2015 and 2010-2026? |
| `HAT_candidate_windows_plot.py` | Observed against modelled shoreline change rate on the candidate windows, both smoothed at the same widths. |
| `HAT_window_1996_2010.py` | The 1996-2010 window on today's setup: run through December 2010, both ends re-solved. |
| `HAT_full_window_1996_2025.py` | One run over the whole record, 1996 through December 2025 (30 model years), with its ends solved against CoastSat; `final --scenario` runs a rung of the management ladder. Moved here from `hatteras_ms/calibration/` 2026-10-23. |
