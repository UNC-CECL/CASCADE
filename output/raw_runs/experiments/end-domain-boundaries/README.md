# end-domain-boundaries

The source/sink rates locked at the two end domains (GIS 1 and 90): what they must carry against each target.

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| [`2026-10-05-ends-solved-on-1996_2025`](2026-10-05-ends-solved-on-1996_2025/README.md) | The ends for the 30-yr full-window run, 1996 through 2025, on the shoreline v2 start. | +144.2227 / +68.4160 (residuals +0.009 / +0.005). GIS 1 gain is only 0.038. | **current** for the full-window run; not in the site config |
| [`2026-10-05-ends-solved-on-net-change-1996_2009`](2026-10-05-ends-solved-on-net-change-1996_2009/README.md) | The ends for the DEM-to-DEM calibration period, 1996 to 2009, solved on the observed net change in metres (not the LRR). | +1.4981 / +10.5659 m/yr (residuals −0.001 / +0.002 m) after five secant steps from zeroBE. | **current**; in the site config for 1996 |
| [`2026-10-05-ends-solved-on-net-change-2009_2025`](2026-10-05-ends-solved-on-net-change-2009_2025/README.md) | The ends for the DEM-to-DEM test period, 2009 to 2025, solved on its own net change (for the edge-only base run). | +32.7049 / +21.0679 m/yr after eight steps (residuals +0.43 / +0.01 m; GIS 1 responds noisily, ±1–2 m). | **current**; edgeBE for 2009 in the site config since 2026-10-06 |

**Status** — **current**: its answer is in use now. **record**: a finished check, kept so the number can be traced.

Studies from before the 2026-10-05 calibration/test plan were moved on 2026-10-07 to `D:\CASCADE_offload\output\raw_runs\experiments\end-domain-boundaries/`, with this topic's full table in the README there (the LRR, dune-line, candidate-window and run-length solves).

Back to [the map](../README.md).
