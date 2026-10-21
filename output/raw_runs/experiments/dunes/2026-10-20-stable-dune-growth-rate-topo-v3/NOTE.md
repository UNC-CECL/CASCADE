# 2026-10-20 dune growth rate from stable dunes, corrected 1996 start (dune-topo v3)

Same method and driver as ../2026-10-20-stable-dune-growth-rate/ (HAT_stable_dune_growth_rate.py; the folder is
chosen by the 1996-start topography version), on the 1996 start rebuilt with the ALACE survey lowered 0.24 m
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
Figures: figures/ (HAT_stable_dune_figures.py). Not adopted; Hannah decides.

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
