# 3-rates/coastsat/lrr_projected/1996_2024 - provenance

Written 2026-09-19 10:56 by scripts/input_prep/5-scr/coastsat_lrr_projected/coastsat_lrr_projected.py.

**Projected** = the transect's 1996-2024 LRR (`../../lrr/1996_2024/transect_lrr_full.csv`) x 28 yr (2024 - 1996). Per domain the mean over its transects; it equals the LRR table's `mean_lrr` x 28 to 0.0e+00 m/yr.

**Observed** = mean CoastSat position over calendar 2024 minus the mean over calendar 1996, per transect, from the same time series the LRR is fitted to. Median positions per transect: 13 in 1996, 33 in 2024. 0 of 906 transects have no position in one of the two years and no observed change.

Seaward positive. `observed_minus_projected_m` > 0 means the shoreline ended more seaward than its long-term trend predicts.

## Island summary

Domain mean projected +4.3 m, observed +1.3 m. Projected landward in 40 of 90 domains, observed in 47 of 90. Observed minus projected per domain: mean -2.9 m, range -51.3 to +56.7 m; r(projected, observed) = 0.92.
