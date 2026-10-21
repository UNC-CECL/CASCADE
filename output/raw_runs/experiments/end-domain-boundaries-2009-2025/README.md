# end-domain-boundaries-2009-2025

The source/sink rates locked at the two end domains (GIS 1 and 90) for the test window, 2009-2025: what they must carry against each target. `HAT_end_solve_net_change.py` files its solves here. Split out of `end-domain-boundaries/` by window on 2026-10-20 (`../RENAMES.md`); the other windows are in [`../end-domain-boundaries-1996-2009/`](../end-domain-boundaries-1996-2009/README.md), [`../end-domain-boundaries-1996-2025/`](../end-domain-boundaries-1996-2025/README.md).

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| [`2026-10-05-ends-solved-on-net-change-2009_2025`](2026-10-05-ends-solved-on-net-change-2009_2025/README.md) | The ends for the DEM-to-DEM test period, 2009 to 2025, solved on its own net change (for the edge-only base run). | +32.7049 / +21.0679 m/yr after eight steps (residuals +0.43 / +0.01 m; GIS 1 responds noisily, ±1–2 m). | record; superseded 2026-10-09 (17-yr test) |
| [`2026-10-09-ends-solved-on-net-change-2009_2025`](2026-10-09-ends-solved-on-net-change-2009_2025/README.md) | The test ends for the 17-yr run (2009–2025, ending 1 Jan 2026). | +32.2422 / +21.5341 m/yr after six steps (residuals +0.29 / +0.13 m). | record; superseded 2026-10-12 |
| `2026-10-12-ends-solved-on-net-change-2009_2025` | The test ends on the ±6-month offset (v3) and the 2025 ±6-month end. | +29.1320 / +21.7947 m/yr. | record; superseded the same evening |
| `2026-10-12-ends-solved-on-calendar-years-2009_2025` | The test ends on the calendar-year setup: start calendar 2009 (offset v4), storms 2010–2025, target calendar 2025 − calendar 2009. | +32.4325 / +21.6164 m/yr after eight steps (residuals −0.05 / +0.00 m; GIS 1 noisy). | superseded 2026-10-13 (dune settings) |
| `2026-10-13-ends-solved-on-dune-settings-2009_2025` | The test ends re-solved on the 10-13 dune settings (one Dmaxel 5.5 m, design 3.94 m MHW). | +29.66 / +20.76 m/yr (converged; `solve.txt`). | in the site config (+29.6576 / +20.7639), but solved on the pre-10-20 storms and dunes; **re-solve needed** |

**Status** — **current**: its answer is in use now. **record**: a finished check, kept so the number can be traced.

Every solve here ran on the old storms and dune settings, and its model states were deleted on 2026-10-20. Solves from before 2026-10-05 are on `D:\archive\CASCADE_offload\output\raw_runs\experiments\end-domain-boundaries\`, with that topic's full table in the README there (the LRR, dune-line, candidate-window and run-length solves).

Back to [the map](../README.md).
