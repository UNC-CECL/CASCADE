# 3-rates/coastsat/projected/2010_2024 - provenance

Written 2026-09-22 10:20 by scripts/input_prep/5-scr/coastsat_total_change/coastsat_total_change.py (--product projected).

## Which product this is

**Projected shoreline change.** The rate is fitted on 1996-2024 but evaluated over 2010-2024, a window it was NOT fitted on. That is what makes it a projection: the long-term trend asked what this window should have looked like. The same window's OWN rate is in `../../total_change/2010_2024/`, and the difference between the two is how much the long-term rate misses this period by before the model is involved.

**Projected shoreline change** = the transect's 1996-2024 LRR (`../../lrr/1996_2024/transect_lrr_full.csv`) x 14 yr (2024 - 2010). Per domain the mean over its transects; it equals the LRR table's `mean_lrr` x 14 to 0.0e+00 m/yr.

**Observed** = mean CoastSat position over calendar 2024 minus the mean over calendar 2010, per transect, from the same time series the LRR is fitted to. No rate anywhere in it. Median positions per transect: 18 in 2010, 33 in 2024. 0 of 906 transects have no position in one of the two years and no observed change.

Seaward positive. `observed_minus_projected_m` > 0 means the shoreline ended more seaward than the trend implies.

## Island summary

Domain mean projected +2.1 m, observed +11.4 m. Landward in 40 of 90 domains, observed in 34 of 90. Observed minus projected per domain: mean +9.3 m, range -27.5 to +84.9 m; r = 0.79.
