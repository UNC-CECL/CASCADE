# 2026-10-20 dune growth rate from stable dunes, Dmax 10.5 m above the berm

Producer: scripts/hatteras_ms/experiments/HAT_stable_dune_growth_rate.py (`run`, then `score`). 1996 start,
managed edgeBE, groin on, Dmaxel 12.2 m NAVD88 (berm 1.7 + 10.5), one run per growth rate.
Run setup on the day: beach slope 0.06 (changed that day in another session) and the working-tree storm code;
all eight rungs share it, so compare them with each other, not with earlier runs.

Method (Hannah's advisor): set r by reproducing a stable dune between 1996 and 2009 with the new ceiling.
Stable = no overwash in any 1996-2009 image, CoastSat 1996-2009 LRR within +/-0.5 m/yr, dune present in all nine
surveys and changing gently without dips. Chosen: GIS 59-62 (LRR -0.10 to +0.24) and GIS 20-22. Crests 6.0-7.3 m,
above Isabel's run-up. Score: model minus lidar domain-median crest, same cells, 1 Jan 2010 vs the 2009 lidar.

| mean r (range) | GIS 20-22 | GIS 59-62 | all seven |
|---|---|---|---|
| 0 | -0.49 | -0.09 | -0.26 |
| 0.005 (0-0.01) | -0.37 | +0.09 | -0.11 |
| 0.01 (0.01) | -0.15 | +0.25 | +0.08 |
| 0.015 (0.01-0.02) | -0.03 | +0.43 | +0.23 |
| 0.02 (0.02) | +0.19 | +0.59 | +0.42 |
| 0.025 (0.02-0.03) | +0.31 | +0.77 | +0.57 |
| 0.03 | +0.53 | +0.93 | +0.76 |
| 0.04 | +0.87 | +1.26 | +1.09 |

Reading: the seven domains together cross zero at a mean r of about 0.008, between 0-0.01 (-0.11 m) and a
fixed 0.01 (+0.08 m). The two stretches disagree: 59-62 wants r ~0.003, 20-22 r ~0.015. The adopted 0-0.01 is
within 0.11 m of the stable dunes; the storm-free estimate (0.003-0.022, median 0.015) was slightly high.
Tables: tables/stable_dune_scores.csv (per survey), stable_dune_summary.csv. Figure: figures/.

## 0.005-0.01 added (same day; run name c1220_rk005_010_s0_groin, thousandths naming added to the launcher)

| mean r (range) | GIS 20-22 | GIS 59-62 | all seven | RMSE over surveys 20-22 / 59-62 |
|---|---|---|---|---|
| 0.005 (0-0.01) | -0.37 | +0.09 | -0.11 | 0.22 / 0.29 |
| 0.0075 (0.005-0.01) | -0.25 | +0.18 | **-0.003** | 0.19 / 0.31 |
| 0.01 (0.01) | -0.15 | +0.25 | +0.08 | 0.17 / 0.34 |

0.005-0.01 reproduces the seven stable domains' 2009 median to within 3 mm and has no zero-growth cells. Not yet
run for 2009-2025 or scored island-wide.

## Island-wide and 2009-2025 check (`check`; tables/island_and_test_check.csv), same setup as the ladder

| | 0-0.01 (adopted) | 0.005-0.01 |
|---|---|---|
| 2010 median / P10 / % < 4.3 m (lidar 5.39 / 2.35 / 31) | 5.43 / 3.33 / 15 | 5.66 / 3.43 / 14 |
| 1996-2009 low-dune P10 RMSE (bias) | 1.00 (+0.92) | 1.03 (+0.96) |
| 1996-2009 Isabel image | 43/48, 21 FA | 43/48, 21 FA |
| 1996-2009 overall skill (hits / FA) | 0.515 (45 / 95) | 0.513 (45 / 96) |
| 1996-2009 shoreline RMSE | 1.284 | 1.276 |
| 2009-2025 low-dune P10 RMSE (bias) | 0.37 (-0.10) | 0.34 (-0.05) |
| 2009-2025 Irene / after Sandy / Mar 2024 | 30/72, 6/13, 2/42 | 30/72, 6/13, 2/42 |
| 2009-2025 overall skill (hits / FA) | 0.150 (50 / 90) | 0.137 (48 / 90) |
| 2009-2025 shoreline RMSE | 2.119 | 2.113 |

Reading: the two are indistinguishable island-wide and in the test window; 0.005-0.01 fits the stable dunes
(-0.003 vs -0.11 m) and has no zero-growth cells, 0-0.01 fits the island-wide 2009 median better (5.43 vs 5.66).

## First-year bulldozing artefact (found drawing stable_dune_fig5)

In the road-managed domains the roadway manager's bulldoze() moves 3,400-5,100 m3 per domain of sand lying above
its road elevation (RoadElevation.csv, 0.4-1.0 m MHW) onto the dunes in model year 1 only (GIS 59-62, 20; none
in later years). It raises the dunes 0.2-0.45 m once, at every growth rate, including r = 0: 74% of the 55
managed road domains rise > 0.1 m; community zones (bdm, no roadway manager) do not. Likely cause: the starting
DEM's road surface sits above the road elevation the manager uses, so the difference is cleared as "overwash".
Not diagnosed further; it is in every managed run, not only these.

Removing it (subtracting the r = 0 run's year-1 rise per group: 59-62 +0.41 m, 20-22 +0.08 m):

| mean r (range) | GIS 20-22 | GIS 59-62 | all seven |
|---|---|---|---|
| 0.005 (0-0.01) | -0.46 | -0.32 | -0.38 |
| 0.0075 (0.005-0.01) | -0.33 | -0.23 | -0.27 |
| 0.01 | -0.23 | -0.15 | -0.19 |
| 0.015 (0.01-0.02) | -0.12 | +0.02 | -0.04 |
| 0.02 | +0.11 | +0.19 | +0.15 |

With the artefact removed the two stretches agree, and the rate is 0.01-0.02. As the model runs now (artefact in),
0.005-0.01 is right; the bulldozing supplied ~0.3 m the growth rate would otherwise have to.

## Why the road elevation and the 1996 start DEM differ (diagnosed same day)

RoadElevation.csv = mean 2009 USACE lidar under the 2008 NC-12 line, m MHW (HAT_road_elevation.py); datum is
consistent end to end (file, runner loader and interior grid are all MHW-relative). bulldoze() flattens interior
rows setback .. setback + 20 m to road_ele. Checked at GIS 20, 45, 55, 59-62 on the 1 m clips:

| GIS | road file | 1 m DEM on the line, 2009 / 1996 start | 10 m road cells, 2009 / 1996 start | road cells from 1996 ALACE |
|---|---|---|---|---|
| 59 | 0.56 | 0.58 / 0.86 | 0.53 / 0.86 | 100% |
| 60 | 0.45 | 0.46 / 0.71 | 0.44 / 0.71 | 96% |
| 61 | 0.58 | 0.60 / 0.79 | 0.56 / 0.79 | 98% |
| 62 | 0.55 | 0.56 / 0.82 | 0.53 / 0.81 | 99% |
| 20 | 0.98 | 0.99 / 1.29 | 0.96 / 1.24 | 99% |
| 45 | 0.76 | 0.79 / 1.04 | 0.74 / 1.01 | 100% |
| 55 | 0.80 | 0.83 / 1.06 | 0.77 / 1.07 | 100% |
(m MHW.) The bulldoze footprint sits on the road (0-20 m off), and 10 m averaging changes < 0.05 m.

Cause: the 1996-start DEM carries the 1996 ALACE survey across the road, and ALACE reads the same pavement
0.20-0.30 m higher than the 2009 lidar. Pavement does not rise, and the offset is the same everywhere, matching
the DEM build audit's bias_2009_minus_1996 of about -0.2 m: a vertical offset between the surveys (instrument or
geoid model), not real change. The 2009-based road elevation is therefore ~0.25 m below the 1996 start's road,
and year 1 clears the difference onto the dunes.
Consequence beyond the road: if ALACE is ~0.25 m high, the 1996 start dunes are too, so the real 1996-2009 dune
growth is ~0.25 m more than measured, and the stable-dune r would rise again (~+0.007, toward 0.02).
Fix options (Hannah decides): (1) a 1996 road elevation sampled from the 1996 start DEM; (2) correct the ALACE
offset in the 1996-start DEM (fixes road and dunes, rebuilds the 1996 start); (3) skip the first-year bulldoze.
