# 3-rates/coastsat/total_change/2010_2026 - provenance

Written 2026-10-03 22:00 by scripts/input_prep/5-scr/3-rates/coastsat/total_change/coastsat_total_change.py (--product total_change).

## Which product this is

**Total shoreline change.** The rate is fitted on 2010-2026 and evaluated over 2010-2026 -- the SAME window -- so nothing is extrapolated and this is not a projection. There is no projected counterpart of this window: the rate IS the 1996-2024 one, so a projection of it onto itself would be these same numbers.

**Total shoreline change** = the transect's 2010-2026 LRR (`../../lrr/2010_2026/transect_lrr_full.csv`) x 16 yr (2026 - 2010). Per domain the mean over its transects; it equals the LRR table's `mean_lrr` x 16 to 0.0e+00 m/yr.

**Observed** = mean CoastSat position over calendar 2026 minus the mean over calendar 2010, per transect, from the same time series the LRR is fitted to. No rate anywhere in it. Median positions per transect: 18 in 2010, 3 in 2026. 0 of 906 transects have no position in one of the two years and no observed change.

Seaward positive. `observed_minus_total_m` > 0 means the shoreline ended more seaward than the trend implies.

## Island summary

Domain mean total +19.5 m, observed +19.1 m. Landward in 30 of 90 domains, observed in 26 of 90. Observed minus total per domain: mean -0.3 m, range -36.6 to +33.3 m; r = 0.93.
