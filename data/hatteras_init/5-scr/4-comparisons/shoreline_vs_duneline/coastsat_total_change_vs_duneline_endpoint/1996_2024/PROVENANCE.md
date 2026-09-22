# 4-comparisons/shoreline_vs_duneline/coastsat_total_change_vs_duneline_endpoint/1996_2024 - provenance

Written 2026-09-22 09:40 by scripts/input_prep/5-scr/total_change_vs_duneline/total_change_vs_duneline.py.

- Shoreline: `3-rates/coastsat/lrr/1996_2024/transect_lrr_full.csv`, lrr_m_yr x 28 yr (the CALENDAR interval; the fit runs 1 January 1996 to 31 December 2024). Fitted and evaluated in the same window, so nothing is extrapolated - this is the fitted trend's net change, not a projection.
- Dune line: `3-rates/duneline/endpoint/1996_2024/` as stored (1997 and 2023 lines), 1997-10-12 to 2023-07-01 (ASSUMED), 25.72 yr.
- Beach-width change = shoreline change - dune-line change; positive = widened.
- Domain means only: the two use different transects (~10 CoastSat, ~5 dune per 500 m domain).

## The interval mismatch

The shoreline is carried over 28 calendar years, the dune line measured over 25.72, a difference of +2.28 yr (Hannah chose the calendar interval, 2026-09-21: the year is what the period means). So the beach-width gap holds that much shoreline drift on top of true width change. Nothing is corrected for it; the interval-matched value is carried beside the headline one in every table (`*_dune_interval_m` = the same rate x the dune line's own span).

Mean beach-width change: **+19.1 m** at 28 yr, +18.7 m at 25.72 yr - a 0.4 m term.

## Island summary

Domain mean shoreline +4.3 m, dune line -14.8 m, beach width +19.1 m (range -31.9 to +107.6); beach narrowed in 20 of 90 domains. r(shoreline, dune line) = 0.77, slope 0.64, RMSE 31.2 m; they agree in sign in 65 of 90. The shoreline is landward in 40, the dune line in 61. y axis ±100 m.
