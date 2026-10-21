# 2026-10-09-ends-solved-on-net-change-2009_2025

**Question.** What GIS 1 and GIS 90 rates does the edge-only test run need now that it runs through 2025?

**Setup.** As the [2026-10-05 solve](../2026-10-05-ends-solved-on-net-change-2009_2025/README.md), with one change: the run is **17 model years (2009–2025, ending 1 Jan 2026)**, was 16. The target is unchanged: the CoastSat mean over 2025-02-17 to 2026-02-17 minus the mean over 2008-08-17 to 2010-08-17, raw at GIS 1 and LOWESS-7 at GIS 90. Full management (fills on), no groin, relocations off, shoreline v2 start. Tolerance 0.34 m.

**Driver.** `scripts/hatteras_ms/HAT_end_solve_net_change.py solve --period 2009` (`SOLVE_DATE` 2026-10-09). Step tables in `solve.txt`; logs `output/logs/driver/dem_to_dem/end_solve_2026-10-09_2009_2025_*.log`.

**Answer.** GIS 1 **+32.2422**, GIS 90 **+21.5341** m/yr after six steps; residuals +0.29 m and +0.13 m.

| step | GIS 1 imposed | GIS 1 model (target +158.00 m) | GIS 90 imposed | GIS 90 model (target +36.04 m) |
|---|---|---|---|---|
| zeroBE | 0.00 | −19.46 | 0.00 | −70.29 |
| 1 | +99.42 | +256.21 | +59.57 | +128.43 |
| 2 | +64.00 | +207.08 | +31.87 | +60.23 |
| 3 | +28.63 | +146.37 | +22.05 | +38.11 |
| 4 | +35.40 | +163.94 | +21.13 | +34.70 |
| 5 | +33.11 | +159.64 | +21.49 | +35.89 |
| 6 | +32.24 | +158.29 | +21.53 | +36.17 |

**Reading.** One more year changes the rates by −0.46 (GIS 1) and +0.47 m/yr (GIS 90). GIS 1 behaved better than in the 10-05 solve (local gain 1.5 m per m/yr; residual +0.29 against +0.43). The interior scores are unchanged (edgeBE −11.6 / 25.1 / 0.42).

**Status.** Current. In `HATTERAS_BE_EDGE_ONLY` for 2009; `matrix/2009_2025/edgeBE/` re-run on it is bit-identical to step 6.
