# 3-rates/coastsat/total_change/2010_2024 - provenance

Written 2026-09-22 10:20 by scripts/input_prep/5-scr/3-rates/coastsat/total_change/coastsat_total_change.py (--product total_change).

## Which product this is

**Total shoreline change.** The rate is fitted on 2010-2024 and evaluated over 2010-2024 -- the SAME window -- so nothing is extrapolated and this is not a projection. `../../projected/2010_2024/` is the other reading of this window: the 1996-2024 rate carried onto it instead of its own.

**Total shoreline change** = the transect's 2010-2024 LRR (`../../lrr/2010_2024/transect_lrr_full.csv`) x 14 yr (2024 - 2010). Per domain the mean over its transects; it equals the LRR table's `mean_lrr` x 14 to 0.0e+00 m/yr.

**Observed** = mean CoastSat position over calendar 2024 minus the mean over calendar 2010, per transect, from the same time series the LRR is fitted to. No rate anywhere in it. Median positions per transect: 18 in 2010, 33 in 2024. 0 of 906 transects have no position in one of the two years and no observed change.

Seaward positive. `observed_minus_total_m` > 0 means the shoreline ended more seaward than the trend implies.

## Island summary

Domain mean total +16.0 m, observed +11.4 m. Landward in 30 of 90 domains, observed in 34 of 90. Observed minus total per domain: mean -4.6 m, range -30.0 to +28.2 m; r = 0.93.
