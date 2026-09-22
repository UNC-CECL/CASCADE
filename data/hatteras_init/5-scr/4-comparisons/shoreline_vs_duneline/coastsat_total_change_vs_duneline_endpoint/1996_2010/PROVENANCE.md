# 4-comparisons/shoreline_vs_duneline/coastsat_total_change_vs_duneline_endpoint/1996_2010 - provenance

Written 2026-09-22 09:40 by scripts/input_prep/5-scr/4-comparisons/shoreline_vs_duneline/total_change_vs_duneline.py.

- Shoreline: `3-rates/coastsat/lrr/1996_2010/transect_lrr_full.csv`, lrr_m_yr x 14 yr (the CALENDAR interval; the fit runs 1 January 1996 to 31 December 2010). Fitted and evaluated in the same window, so nothing is extrapolated - this is the fitted trend's net change, not a projection.
- Dune line: `3-rates/duneline/endpoint/1996_2010/` as stored (1997 and 2009 lines), 1997-10-12 to 2009-05-30, 11.63 yr.
- Beach-width change = shoreline change - dune-line change; positive = widened.
- Domain means only: the two use different transects (~10 CoastSat, ~5 dune per 500 m domain).

## The interval mismatch

The shoreline is carried over 14 calendar years, the dune line measured over 11.63, a difference of +2.37 yr (Hannah chose the calendar interval, 2026-09-21: the year is what the period means). So the beach-width gap holds that much shoreline drift on top of true width change. Nothing is corrected for it; the interval-matched value is carried beside the headline one in every table (`*_dune_interval_m` = the same rate x the dune line's own span).

Mean beach-width change: **+11.7 m** at 14 yr, +12.4 m at 11.63 yr - a 0.8 m term.

## Island summary

Domain mean shoreline -4.5 m, dune line -16.2 m, beach width +11.7 m (range -47.5 to +94.3); beach narrowed in 30 of 90 domains. r(shoreline, dune line) = 0.59, slope 0.47, RMSE 25.3 m; they agree in sign in 59 of 90. The shoreline is landward in 52, the dune line in 73. y axis ±100 m.
