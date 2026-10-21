# groin-blocking

The blocking groin (`BlockingGroinCallback`) and the drift and conserving versions tried on the way to the wave-driven one. **Superseded 2026-10-14** by the wave-driven groin ([`../groin-wave-driven/`](../groin-wave-driven/README.md)). The blocking code was removed from `cascade/groin.py` that day; the scripts ran against commit `4d937b79`. Split out of `groin/` on 2026-10-20 (`../RENAMES.md`). Write-ups are in `hard-structures/groin/4-hindcast/2-blocking-1996-2025/`.

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| `2026-10-05-blocking-fit-dem-to-dem` (b0.3–0.9 × f0.1–0.8, plus dipole M12/f0.3 and M60/f0.6) | Which groin fits the 1996–2009 calibration period? | Blocking b0.6/f0.6: date RMSE 4.0 m vs 69.2 m with no groin; f is loose (0.5–0.8). Write-up `2026-10-05-blocking-fit-calibration/`. | superseded 2026-10-14; `output/calibration/groin/joint_fit.json` is no longer read |
| `2026-10-08-schedule-refit` (instant2004 / instant1996 / ramp1996, each b0.3–0.9 × f0.1–0.8) | Does the groin fit better if it starts failing after the 1995 repair rather than at 2004? | All three schedules score 9.8–10.7 m on the annual CoastSat gap with a weak groin (b × f ≈ 0.36); the photos still prefer the 2004 step. The schedules were carried into the wave-driven fit. | record |
| `2026-10-08-conserving-groin` (b0.60_f0.6 smoke run only) | Does a groin that conserves sand across the GIS 5\|6 face still fit? | Stopped: the pin's imbalance offsets BRIE's own non-conservative solve, so the pin is the one that conserves sand island-wide. | record |
| `2026-10-08-drift-groin` (b0.02–0.6 × f0.1–0.8) | Does a groin that blocks the waves' drift (BRIE untouched) fit 1996–2009? | Joint floor 15.8–16.7 m at b 0.08–0.10 (pin 20.9). With a correction for BRIE's leak at its bends it became the wave-driven groin. | record |

**Status**: **record** means a finished check, kept so the number can be traced. **superseded** means a later study replaced the answer.

Studies from before 2026-10-05 are on `D:\archive\CASCADE_offload\output\raw_runs\experiments\groin\`, with that topic's full table in the README there.

Back to [the map](../README.md).
