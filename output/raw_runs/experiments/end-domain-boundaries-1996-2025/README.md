# end-domain-boundaries-1996-2025

The source/sink rates locked at the two end domains (GIS 1 and 90) for the 30-yr full-window run, 1996 through 2025. `HAT_full_window_1996_2025.py` files its solves here. Split out of `end-domain-boundaries/` by window on 2026-10-20 (`../RENAMES.md`); the other windows are in [`../end-domain-boundaries-1996-2009/`](../end-domain-boundaries-1996-2009/README.md), [`../end-domain-boundaries-2009-2025/`](../end-domain-boundaries-2009-2025/README.md).

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| [`2026-10-05-ends-solved-on-1996_2025`](2026-10-05-ends-solved-on-1996_2025/README.md) | The ends for the 30-yr full-window run, 1996 through 2025, on the shoreline v2 start. | +144.2227 / +68.4160 (residuals +0.009 / +0.005). GIS 1 gain is only 0.038. | superseded 2026-10-12; never in the site config |
| `2026-10-12-ends-solved-on-1996_2025` | The full-window ends on the calendar-year setup. | +137.3 / +68.2 m/yr, both ends converged after three runs (`solve.txt`). | superseded 2026-10-13 (dune settings) |
| `2026-10-13-ends-solved-on-dune-settings-1996_2025` | The full-window ends on the 10-13 dune settings. | Not converged: the last two runs imposed the same rate; next probe GIS 1 +137.3 (`solve.txt`). | unfinished |

**Status** — **current**: its answer is in use now. **record**: a finished check, kept so the number can be traced.

Every solve here ran on the old storms and dune settings, and its model states were deleted on 2026-10-20. Solves from before 2026-10-05 are on `D:\archive\CASCADE_offload\output\raw_runs\experiments\end-domain-boundaries\`, with that topic's full table in the README there (the LRR, dune-line, candidate-window and run-length solves).

Back to [the map](../README.md).
