# 3-rates/coastsat/total_change/1996_2024 - provenance

Written 2026-09-22 10:20 by scripts/input_prep/5-scr/coastsat_total_change/coastsat_total_change.py (--product total_change).

## Which product this is

**Total shoreline change.** The rate is fitted on 1996-2024 and evaluated over 1996-2024 -- the SAME window -- so nothing is extrapolated and this is not a projection. There is no projected counterpart of this window: the rate IS the 1996-2024 one, so a projection of it onto itself would be these same numbers.

**Total shoreline change** = the transect's 1996-2024 LRR (`../../lrr/1996_2024/transect_lrr_full.csv`) x 28 yr (2024 - 1996). Per domain the mean over its transects; it equals the LRR table's `mean_lrr` x 28 to 0.0e+00 m/yr.

**Observed** = mean CoastSat position over calendar 2024 minus the mean over calendar 1996, per transect, from the same time series the LRR is fitted to. No rate anywhere in it. Median positions per transect: 13 in 1996, 33 in 2024. 0 of 906 transects have no position in one of the two years and no observed change.

Seaward positive. `observed_minus_total_m` > 0 means the shoreline ended more seaward than the trend implies.

## Island summary

Domain mean total +4.3 m, observed +1.3 m. Landward in 40 of 90 domains, observed in 47 of 90. Observed minus total per domain: mean -2.9 m, range -51.3 to +56.7 m; r = 0.92.
