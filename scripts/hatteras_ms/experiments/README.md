# experiments

One-off studies, each asking one question of the hindcast. A study drives the
unchanged runner (`../HAT_hindcast_1984_2024.py`) as a subprocess with
`HAT_IGNORE_SETTINGS=1`, or in-process with one input swapped, and files its
runs under `output/raw_runs/experiments/<theme>/<date>-<name>/`. Each script's
header says what it asks and how to run it; `run` / `score` / `plot` actions
are separate so a finished sweep can be re-scored without re-running.

| Theme | Scripts |
|---|---|
| Wave climate | `HAT_wave_grid_smoothed_score.py` (+ `_plot`), `HAT_wave_grid_fixed_ends.py`, `HAT_wave_shortlist_ends_solved.py`, `HAT_wave_recommendation_figures.py`, `HAT_metres_2_wave_sensitivity.py` (+ `_plot`) |
| Island offset | `HAT_metres_1_offset_units.py` (+ `_plot`), `HAT_offset_source_comparison.py` (and its `_div10`, `_option_a` variants), `HAT_offset_source_shoreline_v2.py`, `HAT_offset_source_0922_figures.py` |
| End domains | `HAT_resolve_ends_metres.py`, `HAT_resolve_ends_on_position_change.py`, `HAT_position_change_ends_figure.py` |
| Storms and overwash | `HAT_storm_max_duration.py`, `HAT_storm_length_selection.py`, `HAT_storm_event_splitting.py`, `HAT_storm_height_test.py`, `HAT_trim_length_adopted.py`, `HAT_excess_overwash_diagnosis.py`, `HAT_dune_ceiling_rebuild.py`, `HAT_dune_ceiling_per_domain.py` |
| Code checks | `HAT_metres_3_overwash_fix.py` (+ `_plot_explained`), `HAT_barrier3d_gap_momentum_fix.py`, `HAT_adopt_dune_ceiling_check.py` |
| Topography and domains | `HAT_peaisland_extension.py`, `HAT_run_crest_experiment.py`, `HAT_plot_crest_experiment.py` |
| Road relocation | `HAT_relocation_comparison.py`, `HAT_relocation_period_compare.py`, `HAT_relocation_dune_position_check.py`, `HAT_score_relocation_timing.py`, `HAT_score_road_position.py`; conclusions in `RELOCATION_COMPARISON_RESULTS.md` |

Several scripts import another study's module and re-point its folders (the
`_div10` / `_option_a` variants, the storm studies' shared `MD` / `S`), so
rename or move one only together with the scripts that import it.

`superseded_20260907/` holds the retired 1984 seaward row-insert set, which
cannot be re-run; its `WHY.md` says why.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### HAT_adopt_dune_ceiling_check.py

Does the committed per-cell dune ceiling reproduce the experiment that proposed it?

From the script's original header:

```text
HAT_adopt_dune_ceiling_check.py -- does the committed per-cell ceiling reproduce the experiment?
Step 4 of adopting the per-cell dune ceiling (Hannah, 2026-09-28).

    reproduce  the four trim24 per-cell cases of
               storms-and-overwash/2026-09-28-dune-ceiling-per-domain
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

### HAT_barrier3d_gap_momentum_fix.py

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

### HAT_dune_ceiling_per_domain.py

Does a dune ceiling taken from each domain's own dunes keep the low spots storms break through?

From the script's original header:

```text
HAT_dune_ceiling_per_domain.py -- does a dune ceiling taken from each domain's own dunes keep the low spots storms break through?
WHY (storms-and-overwash/2026-09-28-dune-ceiling-and-rebuild): one island-
wide Dmaxel of 5.5 m NAVD88 matches the 2009 lidar on average and fixes
1996-2010 (PSS 0.60), but every dune grows toward it within a few years. That
erases the real low spots. In 2010-2024 the hit rate falls to 0.23, and Irene
overtops 8 of the 60 higher-dune domains against 47 observed.

THE CEILINGS (Hannah, 2026-09-28: "run the per-domain ceiling test"), each
built from the run's OWN starting dunes (DuneDomain[0]: the 1996-survey mosaic
for 1996-2010, the 2009 lidar for 2010-2024), each held at least FLOOR_M above
the berm, since a ceiling at the berm divides by zero in DuneGrowth:
    dom_median   each domain's Dmaxel = its median starting crest
    dom_p25      each domain's Dmaxel = its 25th-percentile starting crest
    cell         each dune CELL's ceiling = that cell's own starting crest
                 (Barrier3d.DuneGrowth wrapped in-process to take an array;
                 the scalar Dmax it returns, used by the flux limiter and
                 CASCADE's growth-rate reset, is the domain median)
    rebuild rule as now (it made no difference at realistic ceilings);
    storms trim24 and drop72; both windows; managed and natural.
Controls: the current model (3.4 m everywhere) and the uniform 5.5 m ceiling,
both from the earlier experiments.

NOTHING IN THE MAIN CODE CHANGES: in each run's process,
cascade_pipeline.hindcast.build_cascade is wrapped to set the ceilings on the
constructed model before the first step; the storm file is swapped as before.
Scoring applies the same DuneGrowth wrapper, so its crest reconstruction
matches the run.

WHERE: output/raw_runs/experiments/storms-and-overwash/2026-09-28-dune-ceiling-per-domain/

USAGE
    python HAT_dune_ceiling_per_domain.py run [--workers 6]
    python HAT_dune_ceiling_per_domain.py score
```

Notes that were in the code:

```text
The matrix controls these experiments compared against were archived on
2026-09-28 (archive/2026-09-28-loess10-ends/, when the runner's target moved
to LOWESS-7 and the ends were re-solved in a parallel session). Read from there.
```

<details><summary>Function notes (the original docstrings)</summary>

**`install_cell_growth()`**

```text
Barrier3d.DuneGrowth taking a per-cell ceiling where the object carries
one (_hat_cell_dmax, dam above the berm, shape (BarrierLength,)); every
other object runs the original. Same arithmetic as barrier3d.py.
```

**`shoreline_lowess7()`**

```text
Interior RMSE and bias against ONE target for every run (CoastSat LRR,
LOWESS-7, raw for GIS 1-10, as the runner builds it since 2026-09-28), so
runs made before and after the target change compare on equal terms.
```

</details>

### HAT_dune_ceiling_rebuild.py

Do a Hatteras dune ceiling and a rebuild that never lowers dunes fix the excess overwash?

Note (2026-09-30): the `figures` action in the original header below no longer exists; `main()` accepts `run` and `score`.

From the script's original header:

```text
HAT_dune_ceiling_rebuild.py -- does a Hatteras dune ceiling and a rebuild that never lowers dunes fix the excess overwash?
WHY (storms-and-overwash/2026-09-28-excess-overwash-diagnosis): the model's
dunes sit at ~3 m MHW while the 2009 lidar puts Hatteras foredunes at ~4.9 m.
Two settings hold them there:
    Dmaxel   never set for Hatteras, so Barrier3D's default 3.4 m NAVD88
             (3.04 m MHW, Virginia Coast Reserve) is the logistic growth
             ceiling, and every taller dune shrinks every year
    rebuild  when ANY front-row dune cell falls below 1.64 m MHW,
             roadway_manager.rebuild_dunes resets the WHOLE dune field to the
             design height, 3.0 m MHW, cutting down taller dunes

THE EXPERIMENT (Hannah, 2026-09-28: "run the Dmaxel and rebuild experiment")
    Dmaxel     3.4 (current), 5.5, 7.5, 9.0 m NAVD88
    rebuild    current   | as now
               nolower   | same trigger, but cells taller than the design
                           height keep their height (max(old, rebuilt))
               nolower43 | nolower, with the design height 4.3 m NAVD88
                           (3.94 m MHW; the NC-12 safe crest the roadway
                           manager's own docstring cites, Velasquez 2020)
    storms     drop72 (committed) and trim24, the files of
               2026-09-28-storm-length-selection
    windows    1996-2010 and 2010-2024; managed (full_management) every cell,
               natural for Dmaxel only (it has no rebuild)
    The two current-Dmaxel / current-rebuild managed cells and the current-
    Dmaxel natural cells are existing runs (the matrix and the length
    selection), reused.

NOTHING IN THE MAIN CODE CHANGES. In each run's own process, before the
unchanged runner executes: cascade.brie_coupler.set_yaml also writes Dmaxel
into the run's parameter copy, and rebuild_dunes is wrapped in both modules
that call it (roadway_manager, beach_dune_manager). The storm file is swapped
as in the length selection. Barrier3D is the current one (49fd069).

SCORES: overwash against the imagery (the length selection's POD/POFD/PSS/
timing/space, storm-dated), the model's 2010 dune crest against the 2009
lidar (1996-2010 runs), and interior RMSE/bias against CoastSat (at the matrix
edge rates, so not yet a fair shoreline comparison).

WHERE: output/raw_runs/experiments/storms-and-overwash/2026-09-28-dune-ceiling-and-rebuild/

USAGE
    python HAT_dune_ceiling_rebuild.py run [--workers 6]
    python HAT_dune_ceiling_rebuild.py score
    python HAT_dune_ceiling_rebuild.py figures
```

### HAT_excess_overwash_diagnosis.py

Why does the model overwash more than the imagery shows?

From the script's original header:

```text
HAT_excess_overwash_diagnosis.py -- why does the model overwash more than the imagery shows?
THE QUESTION (Hannah, 2026-09-28). Every storm series, the committed one
included, overwashes far more domains than the observed record in most image
windows (storms-and-overwash/2026-09-28-storm-length-selection, stage 1).

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

WHERE: output/raw_runs/experiments/storms-and-overwash/2026-09-28-excess-overwash-diagnosis/
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

### HAT_metres_1_offset_units.py

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
runner builds it (section 8 of HAT_hindcast_1984_2024.py), and every run's
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
The environment one run reads, built the way HAT_run_all builds it:
every HAT_* variable named here, none inherited from the shell.
```

**`coastsat_target()`**

```text
The CoastSat LRR target, GIS 1-90, as the runner builds it for a start
year (section 8 of the runner: LOWESS at 10 domains, the southern 10 raw).
Shared with scripts/hatteras_ms/experiments/HAT_metres_2_wave_sensitivity.py.
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

### HAT_metres_1_offset_units_plot.py

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

### HAT_metres_3_overwash_fix.py

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

### HAT_metres_3_overwash_fix_plot_explained.py

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

### HAT_offset_source_0922_figures.py

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

    python scripts/hatteras_ms/experiments/HAT_offset_source_0922_figures.py
```

Notes that were in the code:

```text
(b) 2010-2024 was never run in this study: the panel carries the
observation alone so the layout matches the others (full management x
both periods, Hannah 2026-09-28); no run is added to a record
```

### HAT_offset_source_comparison.py

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

    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison.py run
    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison.py score
    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison.py plot
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

### HAT_offset_source_comparison_div10.py

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

    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_div10.py run --jobs 4
    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_div10.py score
    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_div10.py plot
```

Notes that were in the code:

```text
2010: the dune line has a pre-metres build; the 2010 shoreline offset was
first built 2026-09-28, so asrun reads its v1. On the 1996 builds the two
differ only in the 15+15 buffer domains, not the real ones.
```

### HAT_offset_source_comparison_option_a.py

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

    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_option_a.py run --jobs 4
    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_option_a.py score
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

### HAT_offset_source_shoreline_v2.py

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

    python scripts/hatteras_ms/experiments/HAT_offset_source_shoreline_v2.py run --jobs 4
    python scripts/hatteras_ms/experiments/HAT_offset_source_shoreline_v2.py score
    python scripts/hatteras_ms/experiments/HAT_offset_source_shoreline_v2.py plot
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

### HAT_peaisland_extension.py

Does the buffer's orientation matter? The domain set extended north over Pea Island to GIS 115.

From the script's original header:

```text
HAT_peaisland_extension.py -- does the buffer's orientation matter?
THE QUESTION (Hannah, 2026-09-16). The hindcast models GIS 1-90 and pads
each end with 15 invented domains that extrapolate the local shoreline
slope, then bridge back to close BRIE's periodic ring; edgeBE pins GIS 1
and 90 to their observed rates with boundary source/sink terms that take up
whatever the buffer gets wrong (+32.2 and +10.0 m/yr in 1996-2010). What if the
buffer carried the REAL coast instead -- Pea Island north to GIS 115 --
with the ends re-solved there? (A one-domain southern extension was run and
removed the same evening: no domain polygon lies south of GIS 1.)

THE DESIGN
    geometries   n115 (GIS 1-115) against base
    offset modes asrun (the compressed planform every calibrated run uses)
                 and detrended (the planform at full strength; no calibrated
                 baseline exists, so base-detrended is solved here too)
    period       1996-2010, full_management, no groin, relocations off, Hs 2.5
    topography   1984-start CURRENT for GIS 1-90; the buffer profile beyond
    management   none on the extension (no road, fills, relocation, BE)
    stage 0      zeroBE on every member: the orientation effect with no
                 boundary term anywhere, against the zeroBE matrix run
    stage 1      the end domains solved by Newton steps (edgeBE, one probe
                 per step through HAT_BE_OVERRIDE), against the edgeBE
                 matrix run
    score        interior RMSE on GIS 2-89 against the SURVEYED target (the
                 runner's rmse_interior_m_yr, identical in meaning for every
                 geometry), the solved end values, and the rates on GIS 80-90

WHERE THINGS ARE
    inputs   2-brie-offset/1996/ext/<geometry>/          the offsets
             5-scr/3-rates/coastsat/lrr/1996_2010/ext/            the targets
    runs     output/raw_runs/experiments/topography-and-domains/2026-09-16-pea-island-domain-extension/<member>/
             one member per <geometry>-<mode>: its 1996_2010/zeroBE/ run is
             stage 0, its step<k>/ folders are the Newton probes, and SOLVED
             names the step that stands as the solved run
    logs     output/raw_runs/experiments/topography-and-domains/2026-09-16-pea-island-domain-extension/logs/<member>/
    answer   RESULTS.md and figures/ beside NOTE.md in that folder

USAGE
    python HAT_peaisland_extension.py stage0                  # the 5 zeroBE runs
    python HAT_peaisland_extension.py check                   # base geometry, edgeBE:
                                                              # must reproduce the matrix row
    python HAT_peaisland_extension.py probe --member n115-asrun --step 1 \
        --override "1=32.2,115=12.0"                          # one Newton probe
    python HAT_peaisland_extension.py next --member n115-asrun  # the next probe, from
                                                              # the runs so far
    python HAT_peaisland_extension.py score                   # RESULTS.md + figures

Each run is the ordinary hindcast runner driven through the environment,
exactly as HAT_run_all.py drives the matrix; nothing here reimplements a
run. Runs are never overwritten (pass --overwrite to redo one).
```

Notes that were in the code:

```text
the detrended planform has no calibrated 90-domain run to compare
against, so its baseline is solved here alongside the extensions
```

```text
the regression check: base geometry, compressed, edgeBE -- must land
on the matrix row's numbers exactly
```

```text
the zeroBE stage-0 run sits under zeroBE/, the probes under edgeBE/;
the solve script reads each run's preset off the index, so every run
is passed with its own tag and nothing else
```

```text
Per end: the solve script prints CONVERGED once the residual is inside
tolerance, and "no secant" once two probes imposed the same value --
which only happens after a converged end's step rounded to 0.0. Either
is done. Converged means every end is.
```

```text
An end the script gave no step for (converged, or no secant)
keeps its last imposed value: the runner needs a rate at every end.
```

```text
a probe row is "<run name> <imposed> <model> <residual>"; the
"next probe" line also starts with HAT_ and has no numbers
```

```text
A probe whose |value| exceeds this is not a boundary term any more; the
matrix values are tens of m/yr and Barrier3D's overwash router has been seen
to die silently under runaway progradation. Stop and say so instead.
```

<details><summary>Function notes (the original docstrings)</summary>

**`run_env()`**

```text
The environment one run reads, built the way HAT_run_all builds it:
every HAT_* variable named here, none inherited from the shell.
```

**`launch()`**

```text
One hindcast run, filed under experiments/<TAG>/<member>/ (stage 0)
or experiments/<TAG>/<member>/step<k>/ (a Newton probe).
```

**`solve_history()`**

```text
(run_name, preset, tag) of the member's stage-0 run and every probe
on disk, oldest first, read off the experiment folder.
```

**`_member_step()`**

```text
(member, step) from a run's tag: <TAG>/<member> is stage 0,
<TAG>/<member>/step<k> the k-th Newton probe. (None, None) otherwise.
```

**`mark_solved()`**

```text
A SOLVED file in the member folder naming the step that is the
solved run, for a reader who is not going to parse the index.
```

**`next_probe()`**

```text
The solve script over the member's runs so far, with the member's
geometry in the environment. Returns (override or None, converged,
text).
```

**`cmd_solve()`**

```text
Newton steps for one member until both ends converge or --max-steps
is reached: next probe, run it, repeat.
```

**`_target()`**

```text
The CoastSat target as the runner builds it: the surveyed GIS 1-90
table (what the interior score uses), or the extension's table over
GIS 0-115 (what an extended run's ends are solved against).
```

**`collect()`**

```text
Every run of the experiment plus the two matrix baselines, as rows:
member, step, preset, index row, and the per-domain LRR.
```

**`draw_compare()`**

```text
The figure to read first: the main approach (base geometry, edgeBE,
compressed planform, the matrix run) against the solved n115 extension
on GIS 1-90 only. One panel (Hannah, 2026-09-16: no difference or
residual panels, no southern extension).
```

**`draw()`**

```text
The whole extended reach: the alongshore LRR of the solved runs
against the target, a panel per offset mode, in house style.
```

</details>

### HAT_plot_crest_experiment.py

The three arms of the GIS 84-86 crest experiment, compared (frozen: two arms' runs are gone).

From the script's original header:

```text
Compares the three arms written by HAT_run_crest_experiment.py.

    pea1989base    v1 as shipped -- GIS 85 setback floored to 0
    pea1989keep    +N rows, the 1996 dune crest left standing in the interior
    pea1989lower   +N rows, that crest shaved to the backdune platform

    FROZEN 2026-09-07. The keep/lower run outputs were deleted with every run
    on modified topography (their topography had gone on 2026-09-03), so this
    script can no longer be re-run; output/raw_runs/experiments/topography-and-domains/2026-09-02-pea-island-row-insert-control/results/ is the
    record. Only pea1989base (v1) still exists under output/raw_runs/.

WHAT THE FIGURES ARE FOR
    Figure 1 is the test. It draws the road setback through time at GIS 84, 85
    and 86 with the prescribed 1989 relocation marked. The claim being checked
    is that the baseline's setback hits zero YEARS BEFORE 1989 and triggers an
    emergent relocation, so the prescribed displacement is then added to a
    synthetic base rather than to the evolved 1984 position. If the inserts
    work, their traces reach 1989 still positive.

    Figure 2 is the control. The insert touches three domains out of ninety, so
    the island-wide shoreline change rate should be unmoved everywhere else. A
    difference out at GIS 40 would mean the insert leaked through the alongshore
    coupling and the experiment is not isolating what it claims to.

READING THE SETBACK TRACE
    roadway_manager keeps `_road_setback_TS` in metres and rewrites it every
    year as `setback += dune_migrated`. A relocation shows as an upward jump: an
    EMERGENT one to the fixed `relocation_setback_m`, a PRESCRIBED one by that
    domain's own measured displacement. Which is which is the point, so both are
    marked rather than left to the eye.

USAGE
    python HAT_plot_crest_experiment.py
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
experiments/topography-and-domains/2026-09-02-pea-island-row-insert-control/<member>/ since 2026-09-16; the two
older layouts (arms/<arm>/ and the loose <arm>/) are tried after it.
```

### HAT_position_change_ends_figure.py

Figure 3 of the wave recommendation, redrawn on the ends solved on position change.

From the script's original header:

```text
Figure 3 redrawn on the position-change ends (2026-09-28).

Hannah: "redraw the figure with the position-change ends when done". The solve
(HAT_resolve_ends_on_position_change.py) runs full management only; this runs
the natural scenario at the same ends (option A waves), then draws
HAT_wave_recommendation_figures.fig3 on the two runs per window. The LRR-ends
figure is left as it is, for comparison.

    python scripts/hatteras_ms/experiments/HAT_position_change_ends_figure.py

Natural runs land under <study>/runs/natural_final/; the figure goes to the
study's figures/ folder.
```

### HAT_relocation_comparison.py

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
    python scripts/hatteras_ms/experiments/HAT_relocation_comparison.py
    python scripts/hatteras_ms/experiments/HAT_relocation_comparison.py --period 1996 --preset edgeBE
    python scripts/hatteras_ms/experiments/HAT_relocation_comparison.py --arm-a DIR --arm-b DIR
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
CSV and GIF in output/comparisons/relocation_1984_2004/<preset>/ (now relocation/1984_2004/<version>/<preset>/) and left
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

### HAT_relocation_dune_position_check.py

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
    python scripts/hatteras_ms/experiments/HAT_relocation_dune_position_check.py
    python scripts/hatteras_ms/experiments/HAT_relocation_dune_position_check.py --preset calibBE
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

### HAT_relocation_period_compare.py

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
    python scripts/hatteras_ms/experiments/HAT_relocation_period_compare.py
    python scripts/hatteras_ms/experiments/HAT_relocation_period_compare.py --presets zeroBE edgeBE --version v2
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

### HAT_resolve_ends_metres.py

Re-solve the two end domains under the metres offset.

From the script's original header:

```text
Re-solve the two end domains under the metres offset (2026-09-27).

Hannah, 2026-09-26/27: sweep the waves with the end source/sink terms FIXED,
"but we might need to resolve the ends now that we fixed the offset". The
stored values (HATTERAS_BE_EDGE_ONLY: 1996 +32.2 / +10.0, 2010 +72.6 / +31.3
m/yr) were solved at Hs 2.5 on the /10 offset. Chosen with Hannah:

    reference  the step-2 baseline wave climate, Hs 1.0 m, Tp 8 s,
               asymmetry 0.8, high-angle 0.45: a neutral setting, not the
               winner of either search, so the fixed ends do not pre-favour
               the sweep that uses them
    scenario   full management (road, beach and dune management, fills; no
               relocations, no groin), as the matrix end values always were;
               the pair is then used for natural and managed runs alike
    target     each window's CoastSat LRR: GIS 1 against the raw domain mean,
               GIS 90 against the LOWESS value (10 domains until 2026-09-28,
               7 since: common.SMOOTH_DOMAINS)
    solve      step 0 a fresh zeroBE run; then HAT_wave_shortlist_ends_solved's
               safeguarded step (secant capped at +-30 m/yr until a probe lies
               on each side of the target, then interpolation held inside the
               bracket), first step from the metres response measured on
               2026-09-26 (about 0.24 m/yr of residual per m/yr imposed at
               GIS 1, 0.20 at GIS 90); converged at |residual| <= 0.02 m/yr,
               at most MAX_STEPS probes
    output     tables/ends.json -- {period: {"1": rate, "90": rate}} -- read by
               HAT_wave_grid_fixed_ends.py. The config (HATTERAS_BE_EDGE_ONLY)
               is NOT changed.

WHERE: output/raw_runs/experiments/end-domain-boundaries/2026-09-27-ends-resolved-metres-offset/

    python scripts/hatteras_ms/experiments/HAT_resolve_ends_metres.py

    --tag   file the solve under another study (2026-09-28: the LOWESS-7 re-solve)
    --seed  "1996=4.8394,17.545;2010=18.8,24.535": step 1 probes these ends
            instead of the first-gain guess from zeroBE, so a re-solve near a
            known answer starts there
```

Notes that were in the code:

```text
Matched on the wave settings (fixed 2026-09-27): the folder holds the
step-0 run of every reference solved so far, and taking the first one
found gave the Hs-2 and Hs-2.5 solves the Hs-1 run's residuals as their
step 0. The values those solves adopted were each confirmed by a probe
run at the right waves, so only the solver's path was affected.
```

```text
Re-solve one window at other waves (2026-09-27, Hannah: "re-solve the
2010 ends at Hs 2"): --periods and the wave settings; the result is
merged into ends.json, the replaced value kept under "history". The
closest probe (smallest worst-end residual) is the answer, and --accept
is the residual the caller will act on (the 0.02 target is not reachable
at GIS 1 in 2010-2024, whose response is not monotonic below ~0.1).
```

### HAT_resolve_ends_on_position_change.py

Solve the two end domains on position change, not LRR.

From the script's original header:

```text
Solve the two end domains on position change, not LRR (2026-09-28).

Hannah, 2026-09-28: "calculate the edge source sink values we should be
matching based on the position change instead" -- option 1 of three: an
experiment, the solver and the config unchanged. The LRR-solved ends match
the rates to 0.15 m/yr but overshoot the observed end-minus-start change at
GIS 1 1996-2010 (+38 m) and GIS 90 2010-2024 (+21 m), because the observed
ends do not move linearly and the model does.

    target     the CoastSat total change, mean of the end year's images minus
               mean of the start year's (5-scr/3-rates/coastsat/total_change/
               <window>/smoothed/tables/domain_smoothed.csv): GIS 1 raw, GIS 90
               LOWESS over common.SMOOTH_DOMAINS (7) -- the same treatment the
               LRR target gives each end (Hannah: "use the smoothed value at
               GIS 90"; built by HAT_wave_recommendation_figures.
               observed_change_smoothed)
    model      end-minus-start position, the run's change_rate_m_yr x 14
    residual   (model - observed) / 14 years, in m/yr, so HAT_resolve_ends_
               metres' solver, gains and 0.02 m/yr tolerance (0.28 m of
               change) apply as they are
    waves      option A (Hs 2.0, Tp 7.5, asymmetry 0.6, high-angle 0.5),
               full management, as the LRR solve
    seed       the LOWESS-7 LRR-solved ends (2026-09-28-ends-resolved-lowess7)

WHERE: output/raw_runs/experiments/end-domain-boundaries/2026-09-28-ends-solved-on-position-change/

    python scripts/hatteras_ms/experiments/HAT_resolve_ends_on_position_change.py
```

<details><summary>Function notes (the original docstrings)</summary>

**`residuals()`**

```text
(model - observed) end-minus-start change, per year. `targets` (the
LRR target R.main builds) is not used.
```

</details>

### HAT_run_crest_experiment.py

Run the 1984-2004 hindcast once per crest-experiment arm.

From the script's original header:

```text
Runs the 1984-2004 hindcast three times -- baseline, insert-with-crest-kept,
insert-with-crest-shaved -- so the GIS 84/85/86 row insert can be judged on
model behaviour rather than on cross-sections.

WHY IT IS A SCRIPT AND NOT THREE COMMANDS
    The road-setback CSV is GLOBAL state: `hatteras_site_config.py:142`
    hardcodes its path, so selecting an arm means overwriting a file every other
    reader of this repo also sees. Doing that by hand leaves the tree pointing at
    an experiment arm the moment anything goes wrong. Here the restore is in a
    `finally`.

    THE TOPOGRAPHY IS NOT SELECTED THAT WAY, and the first version of this
    script got it wrong. It wrote `dune-topo/CURRENT` per arm -- and every arm
    still ran on v1, because `resolve_version` reads the EXTRACTOR's VERSION
    literal before it reads CURRENT. Two of three arms were silent duplicates of
    the control. The version now comes from HAT_TOPO_VERSION_1984_START, which
    outranks the extractor, is scoped to one product, and dies with the
    subprocess instead of persisting.

WHY ARMS RATHER THAN RUN NAMES
    All three produce the SAME run name -- the name is derived from the
    management switches, and those are identical by design. Each is filed as
    an EXPERIMENT (HAT_RUN_KIND=experiment, HAT_RUN_TAG=topography-and-domains/2026-09-02-pea-island-row-insert-control/<arm>)
    under raw_runs/experiments/, and the run index is keyed on
    (run_name, kind, tag), so nothing overwrites anything. The existing
    matrix run is never touched.

RELOCATIONS ON OR OFF -- TWO DIFFERENT QUESTIONS
    --relocations 1 asks: does an emergent relocation fire BEFORE the prescribed
        1989 Pea Island event and corrupt the base its displacement is added to?
        That is a question about whether the prescribed history is applied to the
        right island.

    --relocations 0 asks: left to itself, WHEN does the module relocate the road,
        and how close is that to 1989? That is a question about model skill, and
        it is the one you cannot ask with the events switched on -- prescribing
        the 1989 relocation and then checking whether the model produces 1989 is
        circular.

    The second needs a road that starts BEHIND the dune. A setback floored to 0
    relocates in year 1 by construction and carries no information about
    anything, which is why the baseline arm is a control here and not a
    prediction.

ARMS
    pea1989base    v1              + the setback CSV as shipped -- the ONE arm left

    RETIRED 2026-09-07 (Hannah: keep only unmodified topography). The insert
    arms islandv5 (as-built v5 = today's v4) and blocksv4 (as-built v4 =
    today's v3) went with the layers v3-v8 they ran on; their run outputs
    under output/raw_runs/ were deleted too, as were the outputs of the arms
    already retired on 2026-09-03 (blocksdate*, blocksdsas*, blocksduneline*,
    blocksminimum*, pea1989keep*, pea1989lower*). Sizes and reasons in
    data/hatteras_init/1-barrier3d-domains/archive_purge_20260907.csv.
    HAT_plot_crest_experiment.py's keep/lower comparison is therefore frozen at
    output/raw_runs/experiments/topography-and-domains/2026-09-02-pea-island-row-insert-control/results/.

USAGE
    python HAT_run_crest_experiment.py [--dry-run] [--arms a,b]
```

Notes that were in the code:

```text
Anchored by SEARCHING UPWARD for the project root rather than by
counting parent directories (2026-09-13). A counted depth is correct
only while the file stays where it was written, and these moved into
subfolders of hatteras_ms. Six files here already did it this way.
```

```text
The runner stays at the top of hatteras_ms; this driver moved into
experiments/ on 2026-09-13, so it names the folder rather than its own.
```

```text
Every insert and crest-edit arm is gone (2026-09-07, see the docstring):
islandv5 / blocksv4 with the layers v3-v8; blocksdate*, blocksdsas*,
blocksduneline*, blocksminimum*, pea1989keep*, pea1989lower* had lost
their topography on 2026-09-03 and lost their run outputs on 09-07.
`_arm` is kept so a future arm can be added in one line.

None = "leave the live forcing-tree CSV". That live file is the
v2-measured one (GIS 85/86 floored to 0), not the v1-era one this arm
was first defined against; v1/ carries its own v1-era CSV if that
pairing is wanted (HAT_run_row_insert_set.py's `original` arm uses it).
```

```text
Most of these versions were retired to dune-topo-experiments/ on
2026-09-02. Nothing outside dune-topo/ resolves through topo_dirs or
HAT_TOPO_VERSION_1984_START, so the run would otherwise fail deep in
the hindcast with a missing-array error that names no cause.
```

```text
NOTE: CURRENT is deliberately NOT written. See the docstring --
it loses to the extractor literal, so writing it here would look
like arm selection while doing nothing.
```

```text
Filed as raw_runs/experiments/topography-and-domains/2026-09-02-pea-island-row-insert-control/<member>/,
the member being the old arm name without its pea1989 prefix.
```

```text
ALWAYS put the tree back. An experiment arm left in CURRENT would make
every later road measurement and every later run silently read a
fabricated topography.
```

### HAT_score_relocation_timing.py

Scores each arm's predicted NC-12 relocation year against the 1989 and 1999 events.

From the script's original header:

```text
Scores each arm's PREDICTED NC-12 relocation year against the two documented
events: 1989 Pea Island (GIS 84-87) and 1999 inter-village (GIS 9-14).

ONLY MEANINGFUL WITH THE PRESCRIBED EVENTS OFF. With HATTERAS_ROAD_EVENTS on,
1989 and 1999 are inputs, and scoring the model against its own input is
circular. Pass arms built with --relocations 0.

WHAT IS BEING DISCRIMINATED
    The two independent measurements of the 1984 offset disagree by a factor of
    ~3 at GIS 85 -- the digitized dune line says 65.9 m, the DSAS shoreline
    record says 19.5 m. Neither can be preferred on its own terms. But they
    imply very different relocation dates, and the relocation dates are
    observed. So the timing test is the tie-breaker the two measurements cannot
    provide for each other.

THE HOLDOUT, AND WHY IT MATTERS
    There are TWO events, so an arm can be judged on one and tested on the
    other. `--holdout 1989` scores only the 1999 block, and vice versa. An N
    chosen because it reproduces 1989 has no claim on 1999, and that is the
    check worth having: fitting to both at once produces a better number and no
    way to know whether it means anything.

CENSORING
    A domain whose road never relocates inside 1984-2004 is RIGHT-CENSORED, not
    an error of +20 years. It is reported as ">2004" and excluded from the mean
    error, with the count stated -- averaging a censored value in would quietly
    reward an arm for never relocating anything.

USAGE
    python HAT_score_relocation_timing.py --arms pea1989base
    (the insert arms this compared -- blocksv4, blocksduneline, blocksdsas... --
     lost their run outputs on 2026-09-07; only unmodified topography is kept)
    python HAT_score_relocation_timing.py --holdout 1989
```

Notes that were in the code:

```text
Anchored by SEARCHING UPWARD for the project root rather than by
counting parent directories (2026-09-13). A counted depth is correct
only while the file stays where it was written, and these moved into
subfolders of hatteras_ms. Six files here already did it this way.
```

```text
blocksv4 was the default until 2026-09-07, when every insert arm's run
output was deleted; pea1989base(noreloc) is the one experiment arm left.
```

```text
NAMED FOR THE ARMS SCORED. A fixed filename meant each run
silently replaced the previous arms' scores -- the
blocksduneline/dsas/minimum comparison was lost that way and had
to be re-derived from the runs.
```

### HAT_score_road_position.py

Scores each arm's modelled 2004 road setback against the measured 2004 setback.

From the script's original header:

```text
Scores each arm's MODELLED 2004 road setback against the MEASURED 2004 setback.

WHY THIS EXISTS -- the timing test cannot do it
    HAT_score_relocation_timing.py scores the predicted relocation YEAR. It
    ranks the DSAS-derived N far above the dune-line-derived one (MAE 1.2 yr
    against 11.0). That ranking is not trustworthy on its own, because the
    relocation year is a function of TWO unknowns that trade off exactly:

        small initial setback + trigger at zero
        large initial setback + trigger at a maintenance buffer

    Both reproduce 1989. One observation cannot resolve two unknowns, so a good
    timing score is evidence about the PAIR, not about N.

    The road's POSITION breaks the degeneracy. Where the road physically ended
    up by 2004 is measured -- RoadOffset_2004_domains.csv, an independent
    same-year measurement on the 2004-start topography -- and the two N
    estimates predict different positions regardless of when the move happened.

REQUIRES THE PRESCRIBED RELOCATIONS ON
    The real NC-12 was moved by NCDOT in 1989 and 1999. A model run with the
    events off has not been given those moves, so its 2004 position answers a
    different question. Score arms built with --relocations 1.

WHAT A GOOD SCORE DOES AND DOES NOT MEAN
    Agreement here says the modelled road ends the period where the real one
    did. It does NOT validate the relocation year, the dune history, or the
    fabricated land -- those need their own observations. It is one number
    against one measurement, which is exactly why it is worth having alongside
    the timing test rather than instead of it.

USAGE
    python HAT_score_road_position.py --arms pea1989base
    (the insert arms this compared -- blocksv4, blocksduneline, blocksdsas... --
     lost their run outputs on 2026-09-07; only unmodified topography is kept)
```

Notes that were in the code:

```text
Anchored by SEARCHING UPWARD for the project root rather than by
counting parent directories (2026-09-13). A counted depth is correct
only while the file stays where it was written, and these moved into
subfolders of hatteras_ms. Six files here already did it this way.
```

```text
blocksv4 was the default until 2026-09-07, when every insert arm's run
output was deleted; pea1989base is the one experiment arm left.
```

### HAT_storm_event_splitting.py

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

WHERE: output/raw_runs/experiments/storms-and-overwash/2026-09-29-event-splitting/
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

### HAT_storm_height_test.py

Are the storm water levels too low for Hatteras?

Note (2026-09-30): the `slope` action in the original header below no longer exists; the slope0p10 series is built by `build`, and `main()` accepts `gauges`, `build`, `run` and `score`.

From the script's original header:

```text
HAT_storm_height_test.py -- are the storm water levels too low for Hatteras?
WHY (storms-and-overwash/2026-09-28-dune-ceiling-per-domain): once the dunes
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

WHERE: output/raw_runs/experiments/storms-and-overwash/2026-09-28-storm-height/

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

### HAT_storm_length_selection.py

Which storm duration rule should the hindcast use?

From the script's original header:

```text
HAT_storm_length_selection.py -- which storm duration rule should the hindcast use?
THE QUESTION (Hannah, 2026-09-28): pick the storm series, with overwash that
matches the observed record as the first priority.

Why it is open: the committed series (v3_72) DROPS every grouped event
longer than 72 h, which removes Isabel 2003, March 2018, Florence, Dennis
and Nor'Ida. The 72 h limit existed only because the pre-49fd069 Barrier3D
crashed on long storms (storms-and-overwash/2026-09-28-storm-max-duration).
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

WHERE: output/raw_runs/experiments/storms-and-overwash/2026-09-28-storm-length-selection/

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

### HAT_storm_max_duration.py

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

WHERE: output/raw_runs/experiments/storms-and-overwash/2026-09-28-storm-max-duration/
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

### HAT_trim_length_adopted.py

Does the storm trim length still matter once the dunes are realistic?

From the script's original header:

```text
HAT_trim_length_adopted.py -- does the storm trim length still matter once the dunes are realistic?
WHY (Hannah, 2026-09-28: "run the trim-length check first"). trim24 was chosen
in storms-and-overwash/2026-09-28-storm-length-selection on the OLD dunes (held
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

WHERE: output/raw_runs/experiments/storms-and-overwash/2026-09-28-trim-length-adopted/
```

Notes that were in the code:

```text
The worktree these runs used was removed on 2026-09-28 when hatteras/adopted
was checked out in ../Barrier3D itself (the editable install); same branch.
```

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

    python scripts/hatteras_ms/experiments/HAT_wave_grid_fixed_ends.py run all --jobs 8
    python scripts/hatteras_ms/experiments/HAT_wave_grid_fixed_ends.py score
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
    python scripts/hatteras_ms/experiments/HAT_wave_grid_smoothed_score.py run all
    python scripts/hatteras_ms/experiments/HAT_wave_grid_smoothed_score.py run coarse --periods 1996
    python scripts/hatteras_ms/experiments/HAT_wave_grid_smoothed_score.py score
    python scripts/hatteras_ms/experiments/HAT_wave_grid_smoothed_score.py rescore-step2
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

    python scripts/hatteras_ms/experiments/HAT_wave_grid_smoothed_score_plot.py [--source step2|grid|both]
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

    python scripts/hatteras_ms/experiments/HAT_wave_shortlist_ends_solved.py run [--jobs 8]
    python scripts/hatteras_ms/experiments/HAT_wave_shortlist_ends_solved.py score
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
