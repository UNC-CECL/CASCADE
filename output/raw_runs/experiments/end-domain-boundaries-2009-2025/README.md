# end-domain-boundaries-2009-2025

The source/sink rates locked at the two end domains (GIS 1 and 90) for the test window, 2009-2025: what they must carry against each target. `HAT_end_solve_net_change.py` files its solves here. Split out of `end-domain-boundaries/` by window on 2026-10-20 (`../RENAMES.md`); the other windows are in [`../end-domain-boundaries-1996-2009/`](../end-domain-boundaries-1996-2009/README.md), [`../end-domain-boundaries-1996-2025/`](../end-domain-boundaries-1996-2025/README.md).

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| `2026-10-13-ends-solved-on-dune-settings-2009_2025` | The test ends re-solved on the 10-13 dune settings (one Dmaxel 5.5 m, design 3.94 m MHW). | +29.66 / +20.76 m/yr (converged; `solve.txt`). | in the site config (+29.6576 / +20.7639), superseded 2026-10-20 |
| `2026-10-20-ends-solved-on-corrected-start-2009_2025` | The test ends re-solved on the v3_240 storms, beach slope 0.06, Dmaxel 12.2 m with growth 0.01-0.03 (the 2009 start is unchanged). | +58.3237 / +29.4069 m/yr, step 7 of eight (residuals -0.64 / 0.00 m). GIS 1 scatters about +/-0.7 m near the root, more than the 0.32 m tolerance; step 7 sits on the fitted root (+58.32). Step 8 is not used: the solver left the converged GIS 90 out of its override, so that run fell back to the preset's old +20.76 (fixed in `HAT_end_solve_net_change.py`). | **current**; in the site config |

**Status** — **current**: its answer is in use now. **record**: a finished check, kept so the number can be traced.

Every solve here ran on the old storms and dune settings, and its model states were deleted on 2026-10-20. Solves from before 2026-10-05 are on `D:\hatteras_ch1\model_runs\output\raw_runs\experiments\end-domain-boundaries\`, with that topic's full table in the README there (the LRR, dune-line, candidate-window and run-length solves).

Studies here from before 2026-10-13 were moved on 2026-10-20 to `D:\hatteras_ch1\model_runs\output\raw_runs\experiments\end-domain-boundaries-2009-2025\`, with this topic's table as it stood in the README there.

Back to [the map](../README.md).
