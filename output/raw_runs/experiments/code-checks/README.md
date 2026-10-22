# code-checks

Whether a code or model change moves the results: re-runs against stored runs.

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| `2026-10-15-groin-guard-check` | Does the groin's tested-range guard leave the Buxton hindcast alone, and stop the wave settings the stress test found broken? | Yes: silent and identical (0.0 m) on option A in both windows; stops asymmetry 0.8 (face 46.3°) and high-angle 0.3 (sand balance). | record |
| `2026-10-20-groin-cleanup-check` | Does the 10-20 cleanup of `cascade/groin.py` (renamed `GroinCallback`) change any result? | No: shoreline matrix bit-identical and diagnostics byte-identical in both windows (`old/` vs `new/`). | record |

**Status** — **current**: its answer is in use now. **record**: a finished check, kept so the number can be traced.

Studies from before the 2026-10-05 calibration/test plan were moved on 2026-10-07 to `D:\hatteras_ch1\model_runs\output\raw_runs\experiments\code-checks/`, with this topic's full table in the README there. Among them are the Barrier3D overwash fixes (09-24, 09-28), which every run still uses.

Studies here from before 2026-10-13 were moved on 2026-10-20 to `D:\hatteras_ch1\model_runs\output\raw_runs\experiments\code-checks\`, with this topic's table as it stood in the README there.

Back to [the map](../README.md).
