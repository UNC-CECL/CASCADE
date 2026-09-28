# code-checks

Whether a code or model change moves the results: re-runs against stored runs, and the Barrier3D route_overwash fix.

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| [`2026-09-14-calibrated-pair-rerun-current-code`](2026-09-14-calibrated-pair-rerun-current-code/NOTE.md) | Is the calibrated pair current under the 09-14 code? | Yes: 1984 identical, 2004 moves ≤ 0.0084 m/yr. | record |
| [`2026-09-14-relocation-arm-rerun-new-code`](2026-09-14-relocation-arm-rerun-new-code/NOTE.md) | What does rounding relocations to whole cells change in the 1984–2004 arm? | GIS 11 crosses the drowning line by one cell; otherwise code-only noise. | record |
| [`2026-09-14-relocation-rounding-probes`](2026-09-14-relocation-rounding-probes/NOTE.md) | Single probes around that rounding change. | Unrounded RMSE 0.5231 vs rounded 0.5229. | record |
| [`2026-09-14-site-config-split-check`](2026-09-14-site-config-split-check/NOTE.md) | Did splitting the site config change the model? | No: four runs bit-identical. | record |
| [`2026-09-24-metres-3-barrier3d-overwash-fix`](2026-09-24-metres-3-barrier3d-overwash-fix/NOTE.md) | How much does fixing Barrier3D's route_overwash axis swap move the results? | The bug caused the silent crashes; every run since 09-24 uses the fix (local Barrier3D branch). | **current** (the fix is in use) |

**Status** — **current**: its answer is in use now. **superseded**: a later study
re-asked it; follow the pointer. **record**: a finished check or a result from an
earlier set-up (÷10 offset, Hs 2.5 calibration), kept so the number can be traced.

Back to [the map](../README.md).
