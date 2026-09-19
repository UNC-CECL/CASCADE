# model_vs_observed - the model against an observation, four windows

The four hindcast windows (1984-2004, 1996-2010, 2004-2024, 2010-2024), each
with the modelled shoreline change rate drawn over an observation: edgeBE,
full management, groin off. Two observations, the shoreline (CoastSat) and
the digitised dune line, and for each the run whose two end domains were
solved on it, so the comparison is fair at GIS 1 and 90 by construction.
Same layout, same colours, one y range across every panel.

Named `rate_windows/` until 2026-09-18, when every folder and file was renamed
to say what is compared with what (Hannah): `vs_shoreline/` not `coastsat/`,
`net_change` not `endpoint`, `model_vs_shoreline_and_duneline_` not `both_`.
The figures and tables themselves did not change; the tables were checked
byte-identical after the redraw. One tree since 2026-09-17; before that it
was `observed_vs_modeled_windows/` (CoastSat, 09-15) and
`duneline_vs_modeled_windows/` (the dune line, 09-16).

## Which figure answers which question

Start with `vs_shoreline_and_duneline/model_vs_shoreline_and_duneline_grid.png`.

| question | folder | observation | model line |
|---|---|---|---|
| Does the model follow the shoreline? | `vs_shoreline/domain_means/` | CoastSat per-domain mean LRR, ±1 std | OLS rate, ends solved on CoastSat |
| ...against the smoothed target the run index scores? | `vs_shoreline/smoothed/` | the scoring target (LOESS) as the fill, the means as dots | OLS rate, CoastSat-solved |
| Does the model reproduce the dune line's net change over the window? | `vs_duneline/net_change/` | two surveys differenced, per domain, over the survey interval | endpoint rate, ends solved on the dune line (mean3) |
| ...with the scatter smoothed out, as the CoastSat target is? | `vs_duneline/net_change_smoothed/` | the same, 10-domain LOESS north of D10, raw means D1-10 | endpoint rate, dune-solved |
| Which observation does the model follow, the shoreline or the dune line? | `vs_shoreline_and_duneline/` | both observations as NET CHANGE between the same dune-line dates, as lines with no fill: the CoastSat shoreline (`3-rates/coastsat/endpoint`) blue and the dune line red, each with the target's LOESS | two, both endpoint rate: solid = ends solved on CoastSat, dashed = ends solved on the dune line |

Since 2026-09-18, `vs_shoreline_and_duneline/` compares net change with net
change (Hannah). `skill.csv` scores every run against the CoastSat net change
as well (`cs_endpoint_raw`, `cs_endpoint_loess`). The CoastSat LRR stays the
scoring target under `vs_shoreline/`.

The dune line is scored on NET CHANGE only (2026-09-18, Hannah: "these should
not be lrr, they would just be endpoint, we are tracking net change"). The
old `duneline/lrr/` reading, an OLS through every dune line in the window, is
retired with its product (`5-scr/archive/duneline_lrr_retired_20260918/`).

## Naming

`model_vs_<feature>_<reading>_<start>_<end>.png` for a window, `_grid` for
the 2 x 2 by model period (1984-start left, 1996-start right, the earlier
window of each above the later). The feature is `shoreline` (CoastSat) or
`duneline`; the reading is `means`, `smoothed`, `netchange` or
`netchange_smoothed`. So `model_vs_shoreline_smoothed_1984_2004.png`,
`model_vs_duneline_netchange_grid.png`,
`model_vs_shoreline_and_duneline_2010_2024.png`. Every figure folder keeps
its PDFs and `CAPTIONS.md` under `supporting/`.

```
vs_shoreline/domain_means/         model_vs_shoreline_means_<w>.png
vs_shoreline/smoothed/             model_vs_shoreline_smoothed_<w>.png
vs_duneline/net_change/            model_vs_duneline_netchange_<w>.png
vs_duneline/net_change_smoothed/   model_vs_duneline_netchange_smoothed_<w>.png
vs_shoreline_and_duneline/         model_vs_shoreline_and_duneline_<w>.png
tables/                    domain_rates_<w>.csv   one row per domain: every
                                                  reading of both observations,
                                                  both estimators of every model
                                                  set, the residual against each
                           skill.csv              bias and RMSE (model minus
                                                  observation), GIS 2-89, per
                                                  window x model set x estimator
                                                  x target
runs_used.csv              one row per window per model set: run, arm, folder,
                           timestamp, commit, topography, offsets, and the dune
                           line's vintages, survey dates and interval
y_bounds.txt               the shared y range and the rule behind it
sensitivity/               see below
```

## The runs

Three model sets, named for where their two end domains were solved
(`model_ends` in `runs_used.csv` and `skill.csv`):

| set | what | where the runs are |
|---|---|---|
| `coastsat` | the matrix: ends solved against the CoastSat target | `matrix/` (1984-2004 from `versions/version-pair/v2`, topography v2 as asked; the calibration arm is on v1) |
| `dune-mean3` | the dune-line end solve, mean of the end domain and its two inward neighbours | `experiments/2026-09-16-dune-edgesolve/mean3/` |
| `dune-raw` | the same solve, the end domain's own value | `experiments/2026-09-16-dune-edgesolve/raw/` |

| window | run | matrix arm |
|---|---|---|
| 1984-2004 | `HAT_1984_2004_edgeBE_road_bdm_nogroin` | `version-pair/v2` |
| 1996-2010 | `HAT_1996_2010_edgeBE_road_bdm_nogroin` | calibration (offsets v1, the re-digitized 1997 line) |
| 2004-2024 | `HAT_2004_2024_edgeBE_road_bdm_nourish_nogroin` | calibration (nourishment on) |
| 2010-2024 | `HAT_2010_2024_edgeBE_road_bdm_nourish_nogroin` | calibration (run 09-16, nourishment on) |

The model line is `lrr_m_yr` from the run's `tables/shoreline_change_rate.csv`
(the OLS slope over the annual shorelines, the estimator CoastSat and the run
index use) against an OLS observation, and `change_rate_m_yr` (the endpoint
rate) against a two-survey observation. Each reading is paired with its own
estimator; the one pairing that mixes them is a sensitivity.

## The observations

**CoastSat** is the per-transect LRR through ~250 satellite dates, read from
the window's `transect_lrr_full.csv` through `hat_observed_rates.lrr_csv`;
the scoring target is `cascade_pipeline.hindcast.build_target_table` on it,
exactly as the runner makes it. The observed-only figure is
`data/hatteras_init/5-scr/4-comparisons/coastsat_windows/`, and its drawing is imported
here rather than copied.

**The dune line:**

| window | start line | end line | interval |
|---|---|---|---|
| 1984-2004 | 1984 (1984-09-19) | 2004 (2004-05-25) | 19.68 yr |
| 1996-2010 | 1997 (1997-10-12) | 2009 (2009-05-30) | 11.63 yr |
| 2004-2024 | 2004 (2004-05-25) | 2023 (2023-07-01, assumed) | 19.10 yr |
| 2010-2024 | 2009 (2009-05-30) | 2023 (2023-07-01, assumed) | 14.09 yr |

Vintages through `hat_topo_version.DUNE_LINE_FOR_YEAR`; stations from
`2-brie-offset/raw_offsets/<vintage>_duneline_offset_raw.csv` read as the
hindcast's end-year target loader reads them; seaward positive; dates from
`duneline_vs_coastsat.KNOWN_SURVEY_DATES`. The script READS all of this from
the stored product `5-scr/3-rates/duneline/endpoint/<window>/` (2026-09-18)
rather than computing it, so the figures and the stored numbers cannot
disagree. The 1997, 2009 and 2023 lines were re-digitized on 2026-09-18. Smoothing is the CoastSat target's exact
treatment through the same builder; the LOESS fraction is the 5 km window
over the reach and comes out the same (0.111 vs 0.110). The model line is
never smoothed.

## Sensitivities

```
sensitivity/ends-swapped/      each target against the OTHER solve:
    vs_shoreline/domain_means, smoothed
                                   the shoreline on the dune-solved (mean3) runs
    vs_duneline/net_change, net_change_smoothed
                                   the dune line on the CoastSat-solved matrix
                                   runs (the main level as it was before 09-17)
sensitivity/dune-raw-solve/    vs_duneline/* on the raw-reading dune solve
sensitivity/mixed-estimator/   model_ols_vs_duneline_netchange_<w>.png: the
                               two-survey observation against the model's OLS
                               rate, the one pairing that mixes estimators, kept
                               so the choice is visible
```

Every sensitivity draws the same y range and reads the same `runs_used.csv`
and `tables/` as the main level; `skill.csv` scores all three model sets
against every target, so a pairing without a figure still has its number.

## What it shows (each target against the runs solved on it)

Model minus observation, interior GIS 2-89, m/yr, bias / RMSE:

| window | CoastSat target, CoastSat-solved run, OLS | dune-line target (OLS), dune-solved run, OLS | dune endpoint (LOESS), dune-solved run, endpoint |
|---|---|---|---|
| 1984-2004 | +0.16 / 1.22 | +0.13 / 1.72 | +0.16 / 1.83 |
| 1996-2010 | +0.03 / 1.14 | +0.40 / 1.70 | +0.39 / 1.69 |
| 2004-2024 | -1.00 / 1.79 | -0.87 / 1.42 | -0.62 / 1.51 |
| 2010-2024 | -1.33 / 2.31 | -0.76 / 1.44 | -0.74 / 1.44 |

With each target paired with its own end solve, the dune line is the harder
target in the two early windows (RMSE 1.7 against 1.1-1.2) and the easier one
in the two nourished windows (1.4 against 1.8-2.3), where the waterline
advanced 1-3 m/yr from Cape Point to GIS 40 while the dune line held, and a
model with no beach cannot follow that. The dune line is a rougher
observation (five transects through two or three moments), so part of its
RMSE is scatter: the three-survey OLS is within 0.75 m/yr of the endpoint
everywhere, and smoothing removes 0.2-0.5 from the RMSE and nothing from the
bias. The y range (±12 m/yr) is set by GIS 2 in 1996-2010, where the
1997-2009 dune line moved landward at 10 m/yr; it now applies to the CoastSat
panels too, which sat at ±8 when they were their own tree. The full grid of
window x model set x estimator x target is in `tables/skill.csv`.

## Style choices on record

The model line is black, not the site config's model orange (Hannah,
2026-09-15): the observed line already carries two hues and a fill. In the
smoothed readings the TARGET takes the fill and a light outline, the
per-domain values are DOTS over it, and the black model line is the only line
on the panel: three kinds of mark for three things. Two cuts came before this
(2026-09-15): the target as a heavy sign-coloured line could not be told from
the heavy black one, so of two candidates rendered (target as fill, target
dashed) Hannah chose the fill; then the means as a thin line were a second
same-coloured line beside the fill's outline, doubled over D1-10 where the
target IS the raw mean, so a dots candidate was rendered and chosen.

Drawn by `scripts/analyze_output/compare_runs/HAT_rate_windows.py`, which
imports the panel drawing and the CoastSat reader from the 5-scr producer, the
survey dates and dune reader from
`scripts/input_prep/5-scr/duneline_vs_coastsat/duneline_vs_coastsat.py`, and
resolves run folders through `cascade_pipeline.run_registry` with the arm
named explicitly. `--no-sensitivity` draws the main level only.

## Tracked

This README, `runs_used.csv`, `y_bounds.txt`, `tables/*.csv` and every
figure folder's `supporting/CAPTIONS.md`. Images and PDFs are regenerable and
stay ignored.
