# Dune growth rate, 1996-2009 (calibration period)

Where the runner's 1996 dune growth rate comes from: `DUNE_GROWTH_RATES[1996] = (0.01, 0.03)` in
`HAT_hindcast_1984_2024.ipynb` / `.py`, with Dmaxel 12.2 m (Dmax 10.5 m above the berm) in
`Hatteras-CASCADE-parameters.yaml`. Adopted 2026-10-20. The 2009-2025 rate is in
[`../../2009-start/dune-growth-rate/`](../../2009-start/dune-growth-rate/README.md).

**Method (Hannah's advisor):** pick dunes that stayed stable 1996-2009 and choose the r that makes the
model reproduce their 2009 crest. **Stable domains:** GIS 59-62 (Mid-island / Wimble Shoals) and GIS 20-22
(Buxton-Avon / Avon), seven in all. They have a stable shoreline (CoastSat 1996-2009 LRR within +/-0.5 m/yr),
no overwash in any 1996-2009 image, and a dune present in all nine 1996-2009 surveys, changing gently without
dips. Their crests (6.0-7.3 m) sit above Isabel's run-up (5.58 m), so growth alone sets their change. The
selection is in `HAT_stable_dune_growth_rate.py` (`GROUPS`) and shown in `figures/stable_dune_fig1_where.png`.
Domains left out are shown in `stable_dune_fig3_excluded.png` (GIS 44 beach building out, GIS 65 shoreline
retreating).

| path | what it is |
|---|---|
| `tables/stable_dune_scores.csv` | model minus lidar median crest per survey, rate and domain group |
| `tables/stable_dune_summary.csv` | per rate: 2010 error against the 2009 lidar, and RMSE over the surveys |
| `tables/island_and_test_check.csv` | island-wide and 2009-2025 check of 0.01-0.02, 0.02 and 0.01-0.03 |
| `figures/1-selection/` | why these domains: `stable_dune_fig1_where` (shoreline, dune change and overwash per domain) and `stable_dune_fig3_excluded` (two domains left out) |
| `figures/2-observed-dunes/all-stable-domains/` | the same dune, grown in place, all seven domains together: `stable_dune_fig2_profiles` (1996 and 2009 DEM profiles at GIS 61 and 21), `stable_dune_fig4_surveys` (crests through all nine surveys: no storm loss), `stable_dune_fig8_alongshore` (crest along the shore in every survey; 1996 vs 2009 r 0.73-0.74), `stable_dune_fig9_crest_shift` (87% of 348 rows keep their crest cell within one 10 m cell) |
| `figures/2-observed-dunes/gis<N>/` | one folder per stable domain (GIS 20, 21, 22, 59, 60, 61, 62), every observed-dune figure for that domain alone: `stable_dune_gis<N>_dem` (the 1996 and 2009 1 m DEMs and their difference, north up, true scale, with the dune window and each row's picked crest; was fig7), `_profile`, `_surveys`, `_alongshore` and `_crest_shift` (fig2, fig4, fig8 and fig9 for that domain) |
| `figures/3-model-fit/stable-stretches/` | why 0.01-0.03: `stable_dune_fit_trajectory` (the adopted rate against every survey, median and middle half of the cells), `stable_dune_fit_error_by_rate` (every rate tested: 2009 error and RMSE over the surveys, from `tables/stable_dune_summary.csv`), `stable_dune_fit_spread` (crest change per cell: a range keeps the median and doubles the spread of a fixed 0.02) |
| `figures/3-model-fit/island-wide/` | how 0.01-0.03 does across the island, current runner's 1996 matrix runs: `stable_dune_fit_island` (survey median and 10th percentile through 1996-2009), `stable_dune_fit_island_alongshore` (domain-median crest along the island, 1 January 2010 model vs the 2009 lidar, full management; `_natural`: natural run, groin on; `_mean`, `_natural_mean`: the same with domain means, which fit about the same) |

Each figure folder keeps its PDFs and `CAPTIONS.md` in `supporting/`. Producers: `HAT_stable_dune_figures.py`
(fig1-4, fig8-9, the gis<N> sets, the fit figures; `dems` draws the gis<N> sets, fig8 and fig9 alone, `fit` the fit figures) and the tables, `HAT_stable_dune_growth_rate.py`.

Scripts: `scripts/hatteras_ms/experiments/HAT_stable_dune_growth_rate.py` (`run`, `score`, `check`) and
`HAT_stable_dune_figures.py`. The model runs (`.npz` states, one per rate) and launch logs stay in
`output/raw_runs/experiments/dunes/2026-10-20-stable-dune-growth-rate-topo-v3/`, because the run registry
files every run under `output/raw_runs/`. The first pass, on the 1996 start before the ALACE correction,
is `output/raw_runs/experiments/dunes/2026-10-20-stable-dune-growth-rate/` (superseded; its NOTE.md
sets out the selection in more detail). Moved here from the experiment folder on 2026-10-21.

## The fit (2026-10-20), on the corrected 1996 start (dune-topo v3)

The 1996 start rebuilt with the ALACE survey lowered 0.24 m
(HAT_dem_1984_mosaic.py OVERRIDE_VERTICAL_OFFSET_M) and the survey crests with 1996 / 1997 / 1999 lowered
0.24 / 0.34 / 0.32 m (their offsets from the 2009 lidar on NC-12; 2000 and 2001 left at +0.10 / -0.10).
Dmaxel 12.2 m, groin on, beach slope 0.06.

2010 model minus 2009 lidar, domain-group median crest (m); RMSE over all surveys in brackets:

| r | GIS 20-22 | GIS 59-62 | all seven |
|---|---|---|---|
| 0.01-0.02 (adopted now) | -0.30 (0.20) | -0.21 (0.11) | -0.25 |
| 0.02 | -0.12 (0.14) | +0.11 (0.10) | +0.01 |
| 0.02-0.03 | +0.04 (0.11) | +0.29 (0.17) | +0.18 |
| 0.03 | +0.23 (0.15) | +0.45 (0.24) | +0.36 |

The correction cut the roadway manager's year-1 bulldozing (median domain now none; 24 domains still > 0.1 m;
GIS 59-62 +0.17 m, was +0.41), the leftover from the 20 m road footprint reaching the road's higher edge.
With that leftover removed: 0.02 -> -0.10, 0.02-0.03 -> +0.06 (all seven).
Reading: r = 0.02 matches the stable dunes on the corrected start, as the model runs (+0.01 m) and with the
lowest survey RMSE; without the leftover bulldozing the answer is 0.02-0.025. 0.01-0.02 is now 0.25 m low.
Figures: figures/ (HAT_stable_dune_figures.py). A range about 0.02, 0.01-0.03, was adopted the same day (below).

## Fixed rate or a range? (same day; two wide ranges added to RATES)

Cell-to-cell spread of the 1996-2009 crest change on the 350 stable cells with a crest in 1996 and 2009:

| | median | IQR (width) | 2010 error 20-22 / 59-62 | survey RMSE 20-22 / 59-62 |
|---|---|---|---|---|
| observed | +0.69 | +0.25 to +1.11 (0.86) | | |
| r = 0.02 | +0.72 | +0.64 to +0.85 (0.22) | -0.12 / +0.11 | 0.14 / 0.10 |
| r 0.01-0.03 | +0.68 | +0.47 to +0.91 (0.44) | -0.17 / +0.10 | 0.15 / 0.10 |
| r 0.005-0.035 | +0.57 | +0.24 to +0.93 (0.69) | -0.21 / -0.11 | 0.17 / 0.10 |

Survey-to-survey noise (1997->1998 change) is SD 0.17 m, so most of the observed spread is real. A fixed rate
gets the median but a quarter of the spread; 0.01-0.03 keeps the median and doubles the spread; 0.005-0.035
widens it further but lowers the median. The model's spread is placed at random (Barrier3D draws r per cell),
not where the real fast and slow dunes are.

## Island-wide and 2009-2025 check (`check`; tables/island_and_test_check.csv); 0.01-0.03 adopted the same day

| | 0.01-0.02 | 0.02 | 0.01-0.03 (adopted) |
|---|---|---|---|
| 2010 median / P10 / % < 4.3 m (lidar 5.39 / 2.36 / 31) | 4.88 / 3.38 / 14 | 5.44 / 3.40 / 14 | 5.37 / 3.43 / 14 |
| 1996-2009 low-dune P10 RMSE (bias) | 0.90 (+0.82) | 0.94 (+0.85) | 0.94 (+0.85) |
| 1996-2009 Isabel image | 45/48, 26 FA | 45/48, 22 FA | 45/48, 22 FA |
| 1996-2009 overall skill (hits / FA) | 0.521 (47 / 110) | 0.532 (47 / 103) | 0.533 (47 / 102) |
| 1996-2009 shoreline RMSE | 1.355 | 1.346 | 1.340 |
| 2009-2025 low-dune P10 RMSE (bias) | 0.33 (+0.00) | 0.33 (+0.08) | 0.33 (+0.06) |
| 2009-2025 Irene / after Sandy / Mar 2024 | 30/72, 6/13, 1/42 | 30/72, 6/13, 1/42 | 30/72, 6/13, 1/42 |
| 2009-2025 overall skill (hits / FA) | 0.135 (47 / 88) | 0.123 (44 / 84) | 0.126 (45 / 86) |
| 2009-2025 shoreline RMSE | 2.111 | 2.092 | 2.097 |

Reading: 0.01-0.03 matches the island-wide 2009 median (5.37 vs 5.39) and is level with or best on every other
score; Isabel 45/48 (43 on the uncorrected start). Shoreline RMSE is still well above the old model's (1.23 /
1.66): the ends were solved on 5.5 m / 0.55-0.95, beach slope 0.04 and the uncorrected start; re-solve next.

## One rate island-wide, or one per domain? (2026-10-21; island-wide kept, Hannah)

Rates fitted per domain straight from the surveys, no model runs (`HAT_dune_growth_alongshore_survey_fit.py`
`--widened --intervals`: per cell, the log-odds gain over storm-free survey intervals divided by their years).
The method reproduces the adopted rates (stable domains 0.023, 1996-2009) and finds an alongshore pattern that
repeats between the windows once smoothed (7-domain running medians r 0.74, domain by domain 0.38). Nothing
explains it: beach width runs the wrong way (Spearman -0.38 / -0.32), CoastSat net change shows nothing, and low
dunes regrow faster than the logistic allows (-0.41 / -0.33), a property of the growth law, not of place
(`HAT_dune_growth_beach_width_check.py`). Rodanthe / the S-curves (GIS 76-88) read as rebuilding, not growth.

Sensitivity test (`HAT_dune_growth_alongshore_sensitivity.py`; full management, edgeBE, groin; the smoothed
pattern divided by its mean, times the adopted bounds, through the runner's `HAT_GROWTH_OVERRIDE`):

| | shoreline RMSE 1996 / 2009 | crest median RMSE 1996 / 2009 | overwash skill 1996 / 2009 |
|---|---|---|---|
| island-wide (adopted) | 1.267 / 2.056 | 0.20 / 0.16 | 0.535 / 0.123 |
| each window's own pattern | 1.250 / 2.057 | 0.21 / 0.20 | 0.519 / 0.105 |
| both windows' mean pattern | 1.258 / 2.051 | 0.20 / 0.20 | 0.536 / 0.108 |

Reading: the shoreline moves by under 0.02 m/yr, the dunes and overwash are level or worse, and the modelled
crests follow the surveyed ones alongshore no better. One rate per period, island-wide, is kept. Runs, tables and
figures: `output/raw_runs/experiments/dunes/2026-10-21-alongshore-dune-growth-survey-fit/` and
`.../2026-10-21-alongshore-growth-sensitivity/`.
