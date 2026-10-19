# 2026-10-19 dune ceiling at the measured P99 crest, with slower growth

Question (Hannah, after the advisor's year-0/year-3 crest figure): set the dune ceiling to the measured
dune height and slow the growth rate instead, so the low dunes stay low? Producer:
scripts/hatteras_ms/experiments/HAT_dune_ceiling_dem_growth.py (run / score). Managed edgeBE,
full_management scenario (groin off in all 8 runs, control included). Ceiling = each period's own
starting-DEM crest P99: 9.03 m NAVD88 (1996), 9.50 m (2009). Tables in tables/.

| | control 5.5 m, r 0.55-0.95 | P99, r 0.05-0.15 | P99, r 0.10-0.20 | P99, r 0.15-0.30 |
|---|---|---|---|---|
| 1996-2009 overwash skill (hits / FA) | 0.594 (45 / 44) | 0.443 (37 / 65) | 0.409 (33 / 49) | 0.245 (18 / 13) |
| 1996 low dunes (P10) year 3, start 3.93 m | 5.27 | 4.57 | 4.78 | 5.11 |
| 2009 median crest, model minus lidar | +0.24 m | +1.61 m | +2.33 m | +3.01 m |
| 1996-2009 shoreline RMSE | 1.2383 | 1.2234 | 1.2310 | 1.2353 |
| 2009-2025 overwash skill (hits / FA) | 0.120 (30 / 38) | 0.166 (46 / 68) | 0.152 (43 / 65) | 0.146 (40 / 58) |
| 2009 low dunes (P10) year 3, start 3.05 m | 4.67 | 3.64 | 3.94 | 4.21 |
| 2009-2025 shoreline RMSE | 1.6581 | 1.7812 | 1.7916 | 1.7956 |

Reading: slower growth does keep the low dunes lower early on, and the test window gains hits (30 -> 46),
but with a ceiling near 9 m every dune keeps growing for the whole run. By 2009 the median crest is
1.6-3.0 m above the lidar, Isabel no longer overtops them, and calibration skill drops 0.59 -> 0.25-0.44.
Test false alarms nearly double and test shoreline RMSE worsens. The real dune field barely changed
(median crest 5.73 m in 1996, 5.39 m in 2009); a logistic toward a high ceiling cannot hold it steady.
Not adopted. Untested next step: slow growth with the 5.5 m ceiling kept.

## Follow-up: the 5.5 m ceiling kept, only the growth rate slowed (c55_ variants)

| | control r 0.55-0.95 | r 0.25-0.45 | r 0.15-0.30 | r 0.10-0.20 | r 0.05-0.15 |
|---|---|---|---|---|---|
| 1996-2009 overwash skill (hits / FA) | 0.594 (45 / 44) | **0.612** (47 / 51) | 0.584 (48 / 79) | 0.557 (48 / 96) | 0.540 (48 / 107) |
| 1996 low dunes (P10) year 3, start 3.93 m | 5.27 | 4.82 | 4.61 | 4.46 | 4.34 |
| 2009 median crest, model minus lidar (r) | +0.24 (0.22) | +0.39 (0.40) | +0.23 (0.44) | +0.09 (0.48) | -0.11 (0.52) |
| 1996-2009 shoreline RMSE | 1.2383 | 1.2390 | 1.2314 | 1.2249 | 1.2232 |
| 2009-2025 overwash skill (hits / FA) | 0.120 (30 / 38) | 0.150 (43 / 66) | 0.148 (43 / 67) | 0.163 (47 / 73) | 0.180 (55 / 91) |
| 2009 low dunes (P10) year 3, start 3.05 m | 4.67 | 4.04 | 3.61 | 3.45 | 3.39 |
| 2009-2025 shoreline RMSE | 1.6581 | 1.7088 | 1.7369 | 1.7421 | 1.7440 |

Reading: this is the advisor's mechanism. Slower growth under the 5.5 m ceiling keeps the low dunes near
their surveyed heights, keeps the tall ones (no overgrowth: 2009 crest within 0.4 m of the lidar), and
improves the alongshore crest pattern against the lidar (r 0.22 -> 0.40-0.52). Overwash hits rise in both
windows, but false alarms rise faster, so skill moves little: 0.25-0.45 is the only rate at or above the
control in calibration (0.612) and lifts the test from 0.120 to 0.150. Test shoreline RMSE worsens
0.05-0.09 m/yr. Most test misses are sound-side (10-16 diagnosis), out of reach of any dune setting.
Groin off in all runs. Not adopted; Hannah decides.
