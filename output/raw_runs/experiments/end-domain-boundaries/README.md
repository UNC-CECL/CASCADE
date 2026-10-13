# end-domain-boundaries

The source/sink rates locked at the two end domains (GIS 1 and 90): what they must carry against each target.

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| [`2026-10-05-ends-solved-on-1996_2025`](2026-10-05-ends-solved-on-1996_2025/README.md) | The ends for the 30-yr full-window run, 1996 through 2025, on the shoreline v2 start. | +144.2227 / +68.4160 (residuals +0.009 / +0.005). GIS 1 gain is only 0.038. | **current** for the full-window run; not in the site config |
| [`2026-10-05-ends-solved-on-net-change-1996_2009`](2026-10-05-ends-solved-on-net-change-1996_2009/README.md) | The ends for the DEM-to-DEM calibration period, 1996 to 2009, solved on the observed net change in metres (not the LRR). | +1.4981 / +10.5659 m/yr (residuals −0.001 / +0.002 m) after five secant steps from zeroBE. | record; superseded 2026-10-09 (new calibration end) |
| [`2026-10-05-ends-solved-on-net-change-2009_2025`](2026-10-05-ends-solved-on-net-change-2009_2025/README.md) | The ends for the DEM-to-DEM test period, 2009 to 2025, solved on its own net change (for the edge-only base run). | +32.7049 / +21.0679 m/yr after eight steps (residuals +0.43 / +0.01 m; GIS 1 responds noisily, ±1–2 m). | record; superseded 2026-10-09 (17-yr test) |
| [`2026-10-09-ends-solved-on-net-change-1996_2009`](2026-10-09-ends-solved-on-net-change-1996_2009/README.md) | The calibration ends on the calendar 2008 ±6-month target (2007-07-01 to 2009-06-30). | +1.4041 / +10.3439 m/yr after four steps (residuals −0.10 / −0.14 m). | record; superseded 2026-10-12 |
| [`2026-10-09-ends-solved-on-net-change-2009_2025`](2026-10-09-ends-solved-on-net-change-2009_2025/README.md) | The test ends for the 17-yr run (2009–2025, ending 1 Jan 2026). | +32.2422 / +21.5341 m/yr after six steps (residuals +0.29 / +0.13 m). | record; superseded 2026-10-12 |
| `2026-10-12-ends-solved-on-net-change-1996_2009` | The calibration ends on the ±6-month offset (v3) and the 2008 ±6-month end. | +1.6765 / +10.8869 m/yr. | record; superseded the same evening |
| `2026-10-12-ends-solved-on-net-change-2009_2025` | The test ends on the ±6-month offset (v3) and the 2025 ±6-month end. | +29.1320 / +21.7947 m/yr. | record; superseded the same evening |
| `2026-10-12-ends-solved-on-calendar-years-1996_2009` | The calibration ends on the calendar-year setup: start calendar 1996 (offset v4), storms 1997–2009, target calendar 2009 − calendar 1996. | +1.6727 / +10.7761 m/yr after five steps (residuals +0.005 / −0.035 m). | **current**; in the site config for 1996 |
| `2026-10-12-ends-solved-on-calendar-years-2009_2025` | The test ends on the calendar-year setup: start calendar 2009 (offset v4), storms 2010–2025, target calendar 2025 − calendar 2009. | +32.4325 / +21.6164 m/yr after eight steps (residuals −0.05 / +0.00 m; GIS 1 noisy). | **current**; edgeBE for 2009 in the site config |

**Status** — **current**: its answer is in use now. **record**: a finished check, kept so the number can be traced.

Studies from before the 2026-10-05 calibration/test plan were moved on 2026-10-07 to `D:\archive\CASCADE_offload\output\raw_runs\experiments\end-domain-boundaries/`, with this topic's full table in the README there (the LRR, dune-line, candidate-window and run-length solves).

Back to [the map](../README.md).
