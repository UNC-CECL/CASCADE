# 4-comparisons/shoreline_vs_duneline/coastsat_total_change_vs_duneline_endpoint/2010_2024 - provenance

Written 2026-09-22 09:40 by scripts/input_prep/5-scr/4-comparisons/shoreline_vs_duneline/total_change_vs_duneline.py.

- Shoreline: `3-rates/coastsat/lrr/2010_2024/transect_lrr_full.csv`, lrr_m_yr x 14 yr (the CALENDAR interval; the fit runs 1 January 2010 to 31 December 2024). Fitted and evaluated in the same window, so nothing is extrapolated - this is the fitted trend's net change, not a projection.
- Dune line: `3-rates/duneline/endpoint/2010_2024/` as stored (2009 and 2023 lines), 2009-05-30 to 2023-07-01 (ASSUMED), 14.09 yr.
- Beach-width change = shoreline change - dune-line change; positive = widened.
- Domain means only: the two use different transects (~10 CoastSat, ~5 dune per 500 m domain).

## The interval mismatch

The shoreline is carried over 14 calendar years, the dune line measured over 14.09, a difference of -0.09 yr (Hannah chose the calendar interval, 2026-09-21: the year is what the period means). So the beach-width gap holds that much shoreline drift on top of true width change. Nothing is corrected for it; the interval-matched value is carried beside the headline one in every table (`*_dune_interval_m` = the same rate x the dune line's own span).

Mean beach-width change: **+14.7 m** at 14 yr, +14.8 m at 14.09 yr - a 0.1 m term.

## Island summary

Domain mean shoreline +16.1 m, dune line +1.4 m, beach width +14.7 m (range -23.9 to +81.9); beach narrowed in 21 of 90 domains. r(shoreline, dune line) = 0.66, slope 0.48, RMSE 25.3 m; they agree in sign in 69 of 90. The shoreline is landward in 30, the dune line in 37. y axis ±100 m.
