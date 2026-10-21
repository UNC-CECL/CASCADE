# end-domain-boundaries-1996-2009

The source/sink rates locked at the two end domains (GIS 1 and 90) for the calibration window, 1996-2009: what they must carry against each target. `HAT_end_solve_net_change.py` files its solves here. Split out of `end-domain-boundaries/` by window on 2026-10-20 (`../RENAMES.md`); the other windows are in [`../end-domain-boundaries-2009-2025/`](../end-domain-boundaries-2009-2025/README.md), [`../end-domain-boundaries-1996-2025/`](../end-domain-boundaries-1996-2025/README.md).

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| `2026-10-13-ends-solved-on-dune-settings-1996_2009` | The calibration ends re-solved on the 10-13 dune settings (one Dmaxel 5.5 m, design 3.94 m MHW). | +0.98 / +10.80 m/yr (converged; `solve.txt`). | in the site config (+0.9769 / +10.8017), superseded 2026-10-20 |
| `2026-10-20-ends-solved-on-corrected-start-1996_2009` | The calibration ends re-solved on the v3_240 storms, beach slope 0.06, Dmaxel 12.2 m with growth 0.01-0.03, and the 1996 start with the ALACE offset corrected (dune-topo v3). | +3.6848 / +12.7633 m/yr after seven steps (residuals +0.03 / 0.00 m; `solve.txt`). | **current**; in the site config |

**Status** — **current**: its answer is in use now. **record**: a finished check, kept so the number can be traced.

Every solve here ran on the old storms and dune settings, and its model states were deleted on 2026-10-20. Solves from before 2026-10-05 are on `D:\archive\CASCADE_offload\output\raw_runs\experiments\end-domain-boundaries\`, with that topic's full table in the README there (the LRR, dune-line, candidate-window and run-length solves).

Studies here from before 2026-10-13 were moved on 2026-10-20 to `D:\archive\CASCADE_offload\output\raw_runs\experiments\end-domain-boundaries-1996-2009\`, with this topic's table as it stood in the README there.

Back to [the map](../README.md).
