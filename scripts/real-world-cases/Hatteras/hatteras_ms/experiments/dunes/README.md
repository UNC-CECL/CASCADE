# Dunes

Kept on 2026-10-23: HAT_dune_ceiling_dem_growth.py, HAT_dune_ceiling_growth_grid.py, HAT_dune_ceiling_uniform.py, HAT_dune_crest_vs_surveys.py, HAT_dune_growth_low_rates.py, HAT_low_growth_test_period.py, HAT_stable_dune_figures.py, HAT_stable_dune_growth_rate.py, HAT_test_period_domain_figures.py, HAT_test_period_growth_figures.py, HAT_test_period_growth_rate.py. The other scripts this README describes were retired (see ../RETIRED_STUDIES.md; recover from commit a7c4358c).

The dune ceiling (Dmaxel), dune growth and recovery rates, and how the modelled crest compares with the DEMs and surveys, in the calibration and test windows.

Each study drives the unchanged runner (`../../HAT_hindcast.py`); see `../README.md` for how
the experiments are run and filed, and how they import each other across folders.

| Script | What it asks |
|---|---|
| `HAT_dune_ceiling_uniform.py` | Does one island-wide dune ceiling do as well as the per-cell ceilings, and at what height? |
| `HAT_dune_ceiling_domain_value.py` | Does one dune ceiling per domain, from that domain's own starting crests, beat the per-cell ceilings? |
| `HAT_dune_ceiling_callback.py` | Can a CASCADE-side callback reproduce the per-cell dune ceilings on base Barrier3D dune growth? |
| `HAT_dune_ceiling_range.py` | Across a range of island-wide dune ceilings, how much overwash does the model make, and is it where the record puts it? |
| `HAT_dune_ceiling_per_domain.py` | Does a dune ceiling taken from each domain's own dunes keep the low spots storms break through? |
| `HAT_dune_ceiling_rebuild.py` | Do a Hatteras dune ceiling and a rebuild that never lowers dunes fix the excess overwash? |
| `HAT_dune_ceiling_dem_growth.py` | With the dune ceiling at the measured dune height, can a slower growth rate keep the low dunes low? |
| `HAT_dune_ceiling_growth_grid.py` | Which dune ceiling and growth rate together reproduce the 1996-2009 dunes and overwash best? |
| `HAT_dune_design_height.py` | Does raising the dune design height from 3.0 to 3.94 m MHW (4.3 m NAVD88) change the calibration run? |
| `HAT_dune_growth_rate.py` | How much do the dune growth rate bounds (rmin, rmax) change the dunes and the overwash? |
| `HAT_dune_growth_low_rates.py` | Which dune growth rate, down to none, best follows the surveyed dunes through 1996-2009? |
| `HAT_dune_growth_alongshore_sensitivity.py` | Does an alongshore dune growth pattern, measured from the surveys, change the model's skill? |
| `HAT_dune_growth_alongshore_survey_fit.py` | Does the dune growth rate vary alongshore? Each domain's rate fitted straight from the surveyed crests, no model runs. |
| `HAT_dune_growth_beach_width_check.py` | Does the survey-fit dune growth rate follow beach width? Step 2 of the alongshore growth-rate question. |
| `HAT_dune_recovery_diagnosis.py` | Why do breached dunes never rebuild in the model? |
| `HAT_dune_recovery_rate.py` | Does a dune that recovers from flat, instead of regrowing in proportion to its height, improve the overwash skill? |
| `HAT_dune_recovery_slow_growth.py` | Does a recovery term for flattened dunes remove the chronic overwash that slow dune growth leaves at Cape Point? |
| `HAT_dune_crest_dem_vs_model.py` | Does the model build the dunes up in its first years? |
| `HAT_dune_crest_vs_surveys.py` | How well do the modelled dunes follow the surveyed ones through the calibration and test windows? |
| `HAT_forced_dune_crests.py` | Moved 2026-10-27 to `../../../dune-storm-response/`, with `HAT_road_dune_rebuild.py` and `HAT_dune_collision.py` |
| `HAT_dune_findings_figures.py` | Three summary figures of the dune ceiling (Dmaxel) and growth-rate studies. |
| `HAT_dune_findings_simple_figures.py` | Publication versions of the dune ceiling and growth-rate findings, one point per figure. |
| `HAT_stable_dune_growth_rate.py` | Which dune growth rate reproduces the stable dunes of 1996-2009 under Dmax 10.5 m above the berm? |
| `HAT_stable_dune_figures.py` | Why the stable dunes at GIS 59-62 and 20-22 can calibrate the dune growth rate r, one point per figure. |
| `HAT_test_period_growth_rate.py` | Did the dunes build faster after 2009? Which growth rate reproduces the 2009-2025 dunes under Dmax 10.5 m? |
| `HAT_test_period_growth_figures.py` | Why 0.02-0.04 is the dune growth rate for the 2009-2025 test period, one point per figure. |
| `HAT_test_period_domain_figures.py` | The seven calibration domains followed into the 2009-2025 test period, one folder per domain. |
| `HAT_low_growth_test_period.py` | Do no or very slow dune growth hold up in the 2009-2025 test window, scored on dunes, major events and shoreline? |
| `HAT_measured_ceiling_slow_growth.py` | Can the dune ceiling be the measured dune height, if the growth rate is slowed to the slow-growth candidate? |
| `HAT_zerobe_test_growth.py` | Does the test-period growth rate (0.02-0.04) improve 2009-2025 with no source/sink, full management and the groin? |
| `HAT_cape_growth_rate.py` | Does a faster dune growth rate at Cape Point alone stop the cape eroding where it accretes? |
| `HAT_cape_growth_rate_figure.py` | Did the faster Cape Point dune growth change anything? Cape rate 0.10-0.20 against the calibrated rate, both windows. |
| `HAT_shoreline_rescore_new_ends.py` | Is the shoreline worse because of the new dune settings, or because of the rest of the 2026-10-20 setup? |

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### HAT_dune_ceiling_uniform.py

Does one island-wide dune ceiling do as well as the per-cell ceilings, and at what height?

Asked because a per-cell ceiling taken from the starting DEM never lets a
surveyed low dune grow past its survey height. Each domain's per-cell ceilings
are replaced in-process with one value (5.0, 5.5, 6.0, 6.5 m NAVD88; the
starting crest medians are 5.73 in 1996 and 5.39 in 2009). Managed edgeBE,
both periods, design height 3.94 m MHW, no groin.

Results 2026-10-13 (`output/raw_runs/experiments/dunes/2026-10-13-dune-ceiling-uniform/tables/scores.csv`):

| | per-cell | 5.0 | 5.5 | 6.0 | 6.5 |
|---|---|---|---|---|---|
| 1996-2009 overwash skill | 0.575 | 0.625 | 0.594 | 0.045 | 0.015 |
| 1996-2009 hits / false alarms | 51 / 113 | 48 / 52 | 45 / 44 | 4 / 9 | 2 / 9 |
| 1996-2009 crest vs 2009 lidar, r | 0.55 | 0.12 | 0.22 | 0.21 | 0.19 |
| 1996-2009 dune rebuilds | 44 | 8 | 3 | 1 | 1 |
| 1996-2009 shoreline RMSE | 1.211 | 1.242 | 1.235 | 1.233 | 1.228 |
| 2009-2025 overwash skill | 0.195 | 0.136 | 0.120 | 0.125 | 0.111 |
| 2009-2025 hit rate | 0.41 | 0.22 | 0.19 | 0.19 | 0.17 |
| 2009-2025 shoreline RMSE | 1.745 | 1.662 | 1.657 | 1.664 | 1.664 |
| 2009-2025 cell-years flattened | 8.7% | 3.7% | 3.2% | 3.0% | 2.7% |

A uniform ceiling raises calibration skill only by cutting false alarms, has a
cliff between 5.5 and 6.0 m (overwash all but stops), loses the alongshore
dune pattern, and halves the test-period hit rate. It does end the rebuild
cycle and lowers the test shoreline RMSE by 0.08 m/yr. Not adopted at first;
adopted later the same day at 5.5 m to stay on the base model (see
HAT_dune_ceiling_range.py).

### HAT_dune_ceiling_domain_value.py

Does one dune ceiling per domain, from that domain's own starting crests, beat the per-cell ceilings?

Keeps alongshore variability at the 500 m domain scale while letting low cells
inside a domain grow to the domain's typical height. Each domain's ceiling is
the median (dom50) or 75th percentile (dom75) of its starting cell crests,
in-process only. Same setup as HAT_dune_ceiling_uniform.py, whose per-cell and
island-wide runs it scores alongside.

Results 2026-10-13 (`output/raw_runs/experiments/dunes/2026-10-13-dune-ceiling-domain-value/tables/scores.csv`):

| | per-cell | dom50 | dom75 | uniform 5.0 |
|---|---|---|---|---|
| 1996-2009 overwash skill | 0.575 | 0.460 | 0.376 | 0.625 |
| 1996-2009 hits / false alarms | 51 / 113 | 36 / 45 | 28 / 23 | 48 / 52 |
| 1996-2009 crest vs 2009 lidar, r | 0.55 | 0.56 | 0.54 | 0.12 |
| 1996-2009 dune rebuilds | 44 | 3 | 3 | 8 |
| 1996-2009 shoreline RMSE | 1.211 | 1.203 | 1.226 | 1.242 |
| 2009-2025 overwash skill | 0.195 | 0.111 | 0.124 | 0.136 |
| 2009-2025 hit rate | 0.41 | 0.24 | 0.23 | 0.22 |
| 2009-2025 shoreline RMSE | 1.745 | 1.729 | 1.706 | 1.662 |

The domain values keep the alongshore dune pattern but lose overwash hits in
both periods: raising a domain's low cells to its median closes the gaps the
observed overwash went through. Not adopted.

### HAT_dune_ceiling_callback.py

Can a CASCADE-side callback reproduce the per-cell dune ceilings on base Barrier3D dune growth?

The per-cell ceilings (DuneCeilingFromStart, DuneCeilingFloor) are a local
Barrier3D change (d343461 on hatteras/adopted); Hannah's rule is to build on
Barrier3D, not edit it. The callback turns them off, so Barrier3D grows dunes
on its own scalar ceiling, and each year sets each cell's growth rate to
r * (1 - D/K) / (1 - D/Dmax) during the Barrier3D step only, so the front dune
row grows logistically toward its own starting crest K. The plain rates are
restored before the road and beach managers run. growthparam is one value per
cell shared by both dune rows, so the back row is approximate.

- `callback` (scalar Dmax = the domain's median cell ceiling, as the fork
  reports): fails. Wherever a dune sits near the median, 1 - D/Dmax is near
  zero and the rescaled rate explodes; the back row reached 26.9 m (GIS 6,
  1996 run) and shorelines moved 340-410 m.
- `callback2x` (scalar Dmax = twice the domain's tallest cell ceiling; also
  corrects D for the RSLR SeaLevel() applies before growth): stable.

Results 2026-10-13 (`output/raw_runs/experiments/dunes/2026-10-13-dune-ceiling-callback/tables/scores.csv`):

| | fork per-cell | callback2x |
|---|---|---|
| 1996-2009 overwash skill | 0.575 | 0.553 |
| 1996-2009 shoreline RMSE | 1.2110 | 1.2083 |
| 1996-2009 dune rebuilds | 44 | 38 |
| 2009-2025 overwash skill | 0.195 | 0.205 |
| 2009-2025 shoreline RMSE | 1.7446 | 1.7297 |
| 2009-2025 dune rebuilds | 81 | 128 |
| 2009-2025 low-quartile crest at end | 4.18 m | 3.30 m |
| crest difference, every cell-year, P95 / max | | 0.43-0.56 m / 4.5-5.3 m |
| shoreline difference, max | | 8 m (1996), 23 m (2009) |

Island scores match closely; cell by cell they do not. Base Barrier3D uses its
one scalar Dmax for growth and also as the run-up overwash sediment-flux cap
(fluxLimit = Dmax), the dune restart cap and, through CASCADE, the managers'
growth shut-off. Twice the tallest ceiling raises that flux cap well above the
fork's median, which is the likely source of the extra test-period rebuilds and
the lower low dunes.

### HAT_dune_ceiling_range.py

Across a range of island-wide dune ceilings, how much overwash does the model make, and is it where the record puts it?

Run after the parameter file moved to one island-wide Dmaxel of 5.5 m NAVD88
(2026-10-13); the 5.5 runs are the unchanged runner on that file (metadata
records "uniform Dmaxel 5.50 m NAVD88"). Dmaxel set in-process to 4.5, 5.0,
5.25, 5.75 and 6.0. Managed edgeBE, both periods. space_r = correlation, over
GIS 1-90, of the share of assessed images with observed and with modelled
overwash. Per-domain shares: tables/overwash_by_domain.csv; alongshore figure:
figures/overwash_by_domain_dune_ceiling_range.png.

| | 4.5 | 5.0 | 5.25 | 5.5 | 5.75 | 6.0 | per-cell (fork) | domain median |
|---|---|---|---|---|---|---|---|---|
| 1996-2009 skill | 0.607 | 0.625 | 0.625 | 0.594 | 0.045 | 0.045 | 0.575 | 0.460 |
| 1996-2009 space_r | 0.26 | 0.25 | 0.25 | 0.29 | 0.14 | 0.14 | 0.64 | 0.48 |
| 1996-2009 modelled / observed share | 0.157 / 0.095 | 0.140 | 0.140 | 0.125 | 0.018 | 0.018 | 0.230 | 0.113 |
| 2009-2025 skill | 0.065 | 0.136 | 0.124 | 0.120 | 0.141 | 0.125 | 0.195 | 0.111 |
| 2009-2025 space_r | -0.01 | -0.01 | -0.05 | -0.06 | -0.03 | -0.01 | 0.10 | 0.01 |
| 2009-2025 road zone hit/miss/false alarm | 21/90/70 | 16/95/15 | 14/97/12 | 13/98/11 | 14/97/7 | 12/99/10 | 36/75/62 | |

With one ceiling the dunes converge to one height, so overwash in 1996-2009
is Cape Point plus one storm year that overtops almost every domain at once
(a flat line at ~0.12 in the figure). Above 5.5 m that storm no longer
overtops and overwash all but stops. No value in the range puts overwash where
the record has it; in 2009-2025 none does (space_r ~0), and the per-cell and
per-domain ceilings barely do either.

That storm is Isabel (2003-09-18), run-up 5.58 m NAVD88, the only storm over
5 m in 1997-2025 (2010-2025 peaks at 4.66 m, March 2018). A dune at the
ceiling is overtopped only by a storm whose run-up beats it, so the cliff sits
between 5.5 and 5.75.

Added the same day: the "DEM maximum" ceiling a colleague uses, the 2009 lidar
crest maximum 12.1 m NAVD88 and its P95 8.4 m. Both leave 2 of 68 observed
1996-2009 overwash events (modelled share 1.5% against 9.5% observed) and 25
of 157 in 2009-2025: no storm reaches dunes grown toward those heights. The
12 m maximum also sits at the 12 m clip of the 1996 DEM product.

Also added: each period's own starting-DEM median crest, 5.73 (1996) and 5.39
(2009), run for both periods. 1996-2009 at 5.73: skill 0.045, 4 of 68 hits
(as 5.75), because grown dunes stand above Isabel's 5.58 m run-up; at 5.39:
0.627, as 5.25. 2009-2025 at 5.39: 0.126, 31 hits, as 5.5. Period-specific
medians would break calibration and leave the test unchanged, so one Dmaxel
holds for both periods: set on the 1996-2009 plateau (4.5-5.5 m, skill
0.59-0.63; collapse above ~5.6 m) and held fixed for the test.

**Adopted 2026-10-13 (Hannah):** one island-wide Dmaxel of 5.5 m NAVD88 for
both periods, to stay on the base model. Source: 2026-09-28, where 5.5 m
brought the modelled crests within 0.2 m of the 2009 lidar median (dunes grow
to the ceiling in ~5 yr); confirmed here as the top of the 1996-2009 overwash
plateau (4.5-5.5 m) and held fixed for the test. "Median of both starting DEMs"
was dropped as the justification: it uses the 2009 survey, the end of the
calibration period;
DuneCeilingFromStart and DuneCeilingFloor removed from the parameter file. The
per-domain median (base-compatible, but CASCADE has no per-domain Dmaxel
input) was not taken: better 1996 placement (space_r 0.48) but lower skill.
Overwash location is reproduced only at Cape Point; say so when reporting.

### HAT_dune_ceiling_per_domain.py

Does a dune ceiling taken from each domain's own dunes keep the low spots storms break through?

From the script's original header:

```text
HAT_dune_ceiling_per_domain.py -- does a dune ceiling taken from each domain's own dunes keep the low spots storms break through?
WHY (storms-overwash-and-dunes/2026-09-28-dune-ceiling-and-rebuild): one island-
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

WHERE: output/raw_runs/experiments/storms-overwash-and-dunes/2026-09-28-dune-ceiling-per-domain/

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
WHY (storms-overwash-and-dunes/2026-09-28-excess-overwash-diagnosis): the model's
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

WHERE: output/raw_runs/experiments/storms-overwash-and-dunes/2026-09-28-dune-ceiling-and-rebuild/

USAGE
    python HAT_dune_ceiling_rebuild.py run [--workers 6]
    python HAT_dune_ceiling_rebuild.py score
    python HAT_dune_ceiling_rebuild.py figures
```

### HAT_dune_design_height.py

Does raising the dune design height from 3.0 to 3.94 m MHW (4.3 m NAVD88) change the calibration run?

The 3.0 m MHW design height has no recorded source. 4.3 m NAVD88 is the NC-12
dune crest that Velasquez-Montoya et al. (2021, Natural Hazards Review 22(2))
give for the road not to be overwashed; with MHW 0.36 m it is 3.94 m MHW,
2.6 m above the berm. The value is overridden on the `build_cascade` call
in-process; the runner is unchanged, so its metadata and printout still say
3.0 for the variant (the scorer reads the applied value off the saved
roadway managers). 1996-2009, full_management, edgeBE, no groin.

Results 2026-10-13 (`output/raw_runs/experiments/dunes/2026-10-13-dune-design-height/tables/scores.csv`):

| | control (3.0) | 3.94 |
|---|---|---|
| overwash skill | 0.548 | 0.575 |
| hits / misses / false alarms | 50 / 18 / 121 | 51 / 17 / 113 |
| crest end minus 2009 lidar (median) | -0.10 m | +0.03 m |
| dune rebuilds | 32 | 44 |
| shoreline RMSE interior | 1.2098 m/yr | 1.2110 m/yr |

The control reproduces the matrix run exactly. The extra rebuilds are in GIS
48, 56 and 57, where some dune cells have per-cell ceilings of 1.1-1.25 m
above the berm. A rebuild to 2.6 m puts them far above their ceiling, the
logistic growth term shrinks them below the trigger within a year, and they
are rebuilt again: an annual rebuild cycle that the 3.0 m design height
mostly avoided. Adopted 2026-10-13 (Hannah, option 1: 3.94 as is, the rebuild
cycle noted); runner and notebook set DUNE_DESIGN_ELEVATION_M = 3.94.

### HAT_dune_growth_rate.py

How much do the dune growth rate bounds (rmin, rmax) change the dunes and the overwash?

Barrier3D draws each dune cell's growth rate uniformly between rmin and rmax
with a fixed seed, so only the range changes between variants. rmin/rmax set
in-process to 0.25-0.65 (CASCADE default), 0.35-0.85 (Barrier3D default) and
0.75-1.15; the 0.55-0.95 control is HAT_dune_ceiling_uniform.py's per-cell run.
Managed edgeBE, both periods, design height 3.94 m MHW, per-cell ceilings.

Results 2026-10-13 (`output/raw_runs/experiments/dunes/2026-10-13-dune-growth-rate/tables/scores.csv`):

| | 0.25-0.65 | 0.35-0.85 | 0.55-0.95 (current) | 0.75-1.15 |
|---|---|---|---|---|
| 1996-2009 overwash skill | 0.542 | 0.530 | 0.575 | 0.570 |
| 1996-2009 hits / false alarms | 50 / 125 | 49 / 123 | 51 / 113 | 51 / 116 |
| 1996-2009 crest end minus lidar | +0.02 m | +0.04 m | +0.03 m | +0.06 m |
| 1996-2009 dune rebuilds | 29 | 29 | 44 | 44 |
| 1996-2009 shoreline RMSE | 1.2105 | 1.2101 | 1.2110 | 1.2110 |
| 2009-2025 overwash skill | 0.194 | 0.188 | 0.195 | 0.199 |
| 2009-2025 cell-years flattened | 9.5% | 8.9% | 8.7% | 8.1% |
| 2009-2025 dune rebuilds | 71 | 74 | 81 | 116 |
| 2009-2025 shoreline RMSE | 1.7512 | 1.7446 | 1.7446 | 1.7422 |

The dunes, overwash and shoreline barely move across a doubling of the mean
rate: the per-cell ceilings start every dune at its equilibrium, so the rate
only acts on the few cells a storm has lowered. The one outcome it moves is
the number of rebuilds, because a faster rate also shrinks a rebuilt dune
above its own ceiling faster (the rebuild cycle of HAT_dune_design_height.py).

### HAT_dune_recovery_diagnosis.py

Why do breached dunes never rebuild in the model?

```text
WHY (Hannah, 2026-10-01: "look into why the breached dunes never rebuild").
storms_vs_overwash_1996_2024.png showed false-alarm hotspots at GIS 1-6 and
78-84, where the model's lowest dune crest sits at the berm for most of a run.

READ-ONLY on the four matrix runs (1996/2010 x managed/natural); no runs.
Per domain: cells flattened (crest < 0.5 m above the berm), per-cell ceilings,
regrowth years from DuneRestart to a typical year's largest storm (logistic, as
Barrier3D), shoreline cell moves, road management and dune rebuilds, and the
managed run's false alarms (storms_vs_overwash.cells). Each domain is given the
dominant cause: regrowth slower than the storms, ceiling below a typical
storm, progradation resetting the dune rows, or retreat bringing low interior
into them.

    python HAT_dune_recovery_diagnosis.py

WHERE: output/raw_runs/experiments/storms-overwash-and-dunes/2026-10-01-dune-recovery-diagnosis/
(NOTE.md has the findings, tables/dune_recovery_by_domain.csv the numbers)
```

### HAT_dune_recovery_rate.py

Does a dune that recovers from flat, instead of regrowing in proportion to its height, improve the overwash skill?

```text
WHY (Hannah, 2026-10-01: "test the first two together", after
HAT_dune_recovery_diagnosis.py). A flattened dune cell restarts at 7.5 cm and
grows as r*D*(1 - D/ceiling), so it needs 8-9 years to clear a typical year's
storm and is reset every year. The second fix proposed with it, a ceiling
floored at the alongshore neighbours' median, was dropped before running (Hannah:
"recovery only"): the 2010 north-end lows are whole low stretches of the 2009
lidar, so a +-5..25-cell floor moved the share of cells capped below a typical
storm only from 15% to 10-12% island-wide.

THE CHANGE (in-process only, Barrier3d.DuneGrowth replaced in the run's own
process): the adopted growth plus A * max(0, 1 - D/ceiling) per year, A to the
front dune row and A/3 to the back row (the same Cf = 3 decay). It adds a
fixed amount each year that tapers to zero at the cell's ceiling; the volume
counts in Qdg like the rest of dune growth.
    control  A = 0 (the adopted model; must reproduce the matrix)
    rec015   A = 0.15 m/yr
    rec030   A = 0.30 m/yr
    rec050   A = 0.50 m/yr

RUNS: full_management and natural, both windows, edgeBE, site-config ends,
the adopted storm series. 16 runs.
SCORES: overwash hit rate, false-alarm rate and skill by
storms_vs_overwash's rule (the run's own QowTS swapped into overwash_vs_model),
false alarms by zone (Cape Point, community, road), the share of dune
cell-years flattened, interior RMSE/bias vs LOWESS-7.

NOTHING IN THE MAIN CODE CHANGES.

    python HAT_dune_recovery_rate.py run
    python HAT_dune_recovery_rate.py score

WHERE: output/raw_runs/experiments/storms-overwash-and-dunes/2026-10-01-dune-recovery-rate/
```

### HAT_forced_dune_crests.py

Moved 2026-10-27 to `scripts/real-world-cases/Hatteras/dune-storm-response/`, with
`HAT_road_dune_rebuild.py` and `HAT_dune_collision.py`; its notes and results are in that
folder's README. Its output moved on 2026-10-29 to `output/raw_runs/dune-storm-response/`.
