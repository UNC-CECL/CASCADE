# 2026-10-09-ends-solved-on-net-change-1996_2009

**Question.** What GIS 1 and GIS 90 rates does the calibration period need now that its target ends with calendar 2008 ±6 months?

**Setup.** As the [2026-10-05 solve](../2026-10-05-ends-solved-on-net-change-1996_2009/README.md), with one change: the target is the CoastSat mean over **2007-07-01 to 2009-06-30** (was 2008-08-17 to 2010-08-17) minus the mean over 1995-10-12 to 1997-10-12, raw at GIS 1 and LOWESS-7 at GIS 90. 13 model years (1996–2008), full management, no groin, relocations off, shoreline v2 start. Tolerance 0.26 m. The secant starts from the zeroBE matrix run.

**Driver.** `scripts/hatteras_ms/HAT_end_solve_net_change.py solve --period 1996` (`SOLVE_DATE` 2026-10-09). Step tables in `solve.txt`; logs `output/logs/driver/dem_to_dem/end_solve_2026-10-09_1996_2009_*.log`.

**Answer.** GIS 1 **+1.4041**, GIS 90 **+10.3439** m/yr after four steps; residuals −0.10 m and −0.14 m.

| step | GIS 1 imposed | GIS 1 model (target +9.42 m) | GIS 90 imposed | GIS 90 model (target +0.88 m) |
|---|---|---|---|---|
| zeroBE | 0.00 | −9.23 | 0.00 | −52.35 |
| 1 | +13.67 | +125.61 | +38.99 | +82.41 |
| 2 | +1.89 | +16.47 | +15.40 | +22.86 |
| 3 | +1.13 | +5.46 | +6.69 | −14.99 |
| 4 | +1.40 | +9.32 | +10.34 | +0.74 |

**Reading.** The targets at the ends moved little (GIS 1 +10.65 → +9.42 m, GIS 90 +2.17 → +0.88 m), so the rates moved little (−0.09 and −0.22 m/yr). Gains as before: about 14 m per m/yr at GIS 1, 4.3 at GIS 90. The interior scores are unchanged (edgeBE +2.4 / 17.0 / 0.53).

**Status.** Current. In `HATTERAS_BE_EDGE_ONLY` for 1996; `matrix/1996_2009/edgeBE/` re-run on it is bit-identical to step 4.
