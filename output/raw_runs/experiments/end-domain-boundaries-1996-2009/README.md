# end-domain-boundaries-1996-2009

The source/sink rates locked at the two end domains (GIS 1 and 90) for the calibration window, 1996-2009: what they must carry against each target. `HAT_end_solve_net_change.py` files its solves here. Split out of `end-domain-boundaries/` by window on 2026-10-20 (`../RENAMES.md`); the other windows are in [`../end-domain-boundaries-2009-2025/`](../end-domain-boundaries-2009-2025/README.md), [`../end-domain-boundaries-1996-2025/`](../end-domain-boundaries-1996-2025/README.md).

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| [`2026-10-05-ends-solved-on-net-change-1996_2009`](2026-10-05-ends-solved-on-net-change-1996_2009/README.md) | The ends for the DEM-to-DEM calibration period, 1996 to 2009, solved on the observed net change in metres (not the LRR). | +1.4981 / +10.5659 m/yr (residuals −0.001 / +0.002 m) after five secant steps from zeroBE. | record; superseded 2026-10-09 (new calibration end) |
| [`2026-10-09-ends-solved-on-net-change-1996_2009`](2026-10-09-ends-solved-on-net-change-1996_2009/README.md) | The calibration ends on the calendar 2008 ±6-month target (2007-07-01 to 2009-06-30). | +1.4041 / +10.3439 m/yr after four steps (residuals −0.10 / −0.14 m). | record; superseded 2026-10-12 |
| `2026-10-12-ends-solved-on-net-change-1996_2009` | The calibration ends on the ±6-month offset (v3) and the 2008 ±6-month end. | +1.6765 / +10.8869 m/yr. | record; superseded the same evening |
| `2026-10-12-ends-solved-on-calendar-years-1996_2009` | The calibration ends on the calendar-year setup: start calendar 1996 (offset v4), storms 1997–2009, target calendar 2009 − calendar 1996. | +1.6727 / +10.7761 m/yr after five steps (residuals +0.005 / −0.035 m). | superseded 2026-10-13 (dune settings) |
| `2026-10-13-ends-solved-on-dune-settings-1996_2009` | The calibration ends re-solved on the 10-13 dune settings (one Dmaxel 5.5 m, design 3.94 m MHW). | +0.98 / +10.80 m/yr (converged; `solve.txt`). | in the site config (+0.9769 / +10.8017), but solved on the pre-10-20 storms and dunes; **re-solve needed** |

**Status** — **current**: its answer is in use now. **record**: a finished check, kept so the number can be traced.

Every solve here ran on the old storms and dune settings, and its model states were deleted on 2026-10-20. Solves from before 2026-10-05 are on `D:\archive\CASCADE_offload\output\raw_runs\experiments\end-domain-boundaries\`, with that topic's full table in the README there (the LRR, dune-line, candidate-window and run-length solves).

Back to [the map](../README.md).
