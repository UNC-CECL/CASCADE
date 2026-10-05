# 2026-10-05-ends-solved-on-net-change-1996_2009

**Question.** What source/sink rates do the two locked end domains (GIS 1 and GIS 90) need so that the model's end-minus-start shoreline matches the observed net change over the calibration period?

**Setup.** The DEM-to-DEM plan has calibration from 1996 to 2009 (13 model years) and testing from 2009 to 2025. The target is `5-scr/3-rates/coastsat/net_change/1996_2009/`: the CoastSat mean over 2008-08-17 to 2010-08-17 minus the mean over 1995-10-12 to 1997-10-12, seaward positive. It is read raw at GIS 1 and at 7-domain LOWESS at GIS 90, as the runner reads its target; the model side is the end domain's own value. Runs use full management, no groin, relocations off, and the shoreline v2 start. The secant starts from the zeroBE matrix run (both ends 0). The convergence tolerance is 0.02 m/yr × 13 yr = 0.26 m.

**Driver.** `scripts/hatteras_ms/HAT_end_solve_net_change.py solve --period 1996`, which calls `be_edge_domain_solve.py --target net_change` after every step. The step-by-step tables are in `solve.txt`; the logs are in `output/logs/driver/dem_to_dem/end_solve_1996_2009_*.log`.

**Answer.** GIS 1 **+1.4981**, GIS 90 **+10.5659** m/yr. The final residuals are −0.001 m and +0.002 m.

| step | GIS 1 imposed | GIS 1 model (target +10.65 m) | GIS 90 imposed | GIS 90 model (target +2.17 m) |
|---|---|---|---|---|
| zeroBE | 0.00 | −9.23 | 0.00 | −52.35 |
| 1 | +14.56 | +128.60 | +39.94 | +84.40 |
| 2 | +2.10 | +20.32 | +15.92 | +24.86 |
| 3 | +0.99 | +3.42 | +6.77 | −14.59 |
| 4 | +1.46 | +10.16 | +10.66 | +2.58 |
| 5 | +1.50 | +10.65 | +10.57 | +2.17 |

**Reading.** GIS 1 responds strongly: about 14 m of net change per m/yr, ten times the nominal first-step gain, which is why step 1 overshot. GIS 90 gives about 4.4 m per m/yr. Both ends are much smaller than the LRR-solved 1996–2015 values (+3.39 / +37.60).

**Status.** Current. Written to `HATTERAS_BE_EDGE_ONLY` for 1996, and for 2009 unchanged (the test period runs on the calibration sources and sinks).
