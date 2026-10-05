# 3-rates/coastsat/total_change/1996_2015 - provenance

Written 2026-10-03 22:00 by scripts/input_prep/5-scr/3-rates/coastsat/total_change/coastsat_total_change.py (--product total_change).

## Which product this is

**Total shoreline change.** The rate is fitted on 1996-2015 and evaluated over 1996-2015 -- the SAME window -- so nothing is extrapolated and this is not a projection. There is no projected counterpart of this window: the rate IS the 1996-2024 one, so a projection of it onto itself would be these same numbers.

**Total shoreline change** = the transect's 1996-2015 LRR (`../../lrr/1996_2015/transect_lrr_full.csv`) x 19 yr (2015 - 1996). Per domain the mean over its transects; it equals the LRR table's `mean_lrr` x 19 to 0.0e+00 m/yr.

**Observed** = mean CoastSat position over calendar 2015 minus the mean over calendar 1996, per transect, from the same time series the LRR is fitted to. No rate anywhere in it. Median positions per transect: 13 in 1996, 17 in 2015. 0 of 906 transects have no position in one of the two years and no observed change.

Seaward positive. `observed_minus_total_m` > 0 means the shoreline ended more seaward than the trend implies.

## Island summary

Domain mean total -6.6 m, observed -7.8 m. Landward in 51 of 90 domains, observed in 52 of 90. Observed minus total per domain: mean -1.2 m, range -49.2 to +24.2 m; r = 0.88.
