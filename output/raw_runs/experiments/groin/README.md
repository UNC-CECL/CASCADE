# groin

The Cape Point groin (GIS 5/6): how strong it must be, and in what form, to hold the observed fillet.

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| `2026-10-05-blocking-fit-dem-to-dem` (b0.3–0.9 × f0.1–0.8, plus dipole M12/f0.3 and M60/f0.6) | Which groin fits the 1996–2009 calibration period? | Blocking b0.6/f0.6: date RMSE 4.0 m vs 69.2 m with no groin; f is loose (0.5–0.8). Write-up in `hard-structures/groin/groin-module-test/1-dem-to-dem/2026-10-05-blocking-fit-calibration/`. | **current**; pinned in `output/calibration/groin/joint_fit.json` |
| `2026-10-08-schedule-refit` (instant2004 / instant1996 / ramp1996, each b0.3–0.9 × f0.1–0.8) | Does the groin fit better if it starts failing after the 1995 repair rather than at 2004? | On the annual CoastSat gap all three schedules score 9.8–10.7 m with a weak groin (b × f ≈ 0.36, the same as the pin, so 2009–2025 runs are identical); the photos still prefer the 2004 step. Write-up in `hard-structures/groin/groin-module-test/1-dem-to-dem/2026-10-08-schedule-refit/`. | **current**; decision pending |
| `2026-10-08-conserving-groin` (b0.60_f0.6 smoke run only) | Does a groin that conserves sand across the GIS 5\|6 face still fit? | Stopped: the pinned groin's imbalance offsets BRIE's own non-conservative solve, so the pin is the one that conserves sand island-wide. Write-up in `hard-structures/groin/groin-module-test/1-dem-to-dem/2026-10-08-conserving-groin/`. | record |

**Status** — **current**: its answer is in use now. **record**: a finished check, kept so the number can be traced.

Studies from before the 2026-10-05 calibration/test plan were moved on 2026-10-07 to `D:\CASCADE_offload\output\raw_runs\experiments\groin/`, with this topic's full table in the README there (the option-A and instant-2004 grids).

Back to [the map](../README.md).
