# 3-rates/coastsat/total_change/1996_2010 - provenance

Written 2026-09-22 10:20 by scripts/input_prep/5-scr/3-rates/coastsat/total_change/coastsat_total_change.py (--product total_change).

## Which product this is

**Total shoreline change.** The rate is fitted on 1996-2010 and evaluated over 1996-2010 -- the SAME window -- so nothing is extrapolated and this is not a projection. `../../projected/1996_2010/` is the other reading of this window: the 1996-2024 rate carried onto it instead of its own.

**Total shoreline change** = the transect's 1996-2010 LRR (`../../lrr/1996_2010/transect_lrr_full.csv`) x 14 yr (2010 - 1996). Per domain the mean over its transects; it equals the LRR table's `mean_lrr` x 14 to 0.0e+00 m/yr.

**Observed** = mean CoastSat position over calendar 2010 minus the mean over calendar 1996, per transect, from the same time series the LRR is fitted to. No rate anywhere in it. Median positions per transect: 13 in 1996, 18 in 2010. 0 of 906 transects have no position in one of the two years and no observed change.

Seaward positive. `observed_minus_total_m` > 0 means the shoreline ended more seaward than the trend implies.

## Island summary

Domain mean total -4.5 m, observed -10.1 m. Landward in 52 of 90 domains, observed in 53 of 90. Observed minus total per domain: mean -5.6 m, range -43.1 to +16.8 m; r = 0.90.
