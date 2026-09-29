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
| [`2026-09-28-barrier3d-overwash-gap-momentum-fix`](2026-09-28-barrier3d-overwash-gap-momentum-fix/NOTE.md) | How much do three further overwash fixes (DuneGaps cells, gap discharge slice, inundation momentum C) move the results? | They add 2–15% more overwash domain-years. Mean retreat grows by 0.1–2.6 m; the observed hit rate rises slightly and RMSE barely changes. Merged into Barrier3D `hatteras/adopted` (local, not pushed). | **current**: adopted 2026-09-28 |
| [`2026-09-28-per-cell-dune-ceiling-reproduces`](2026-09-28-per-cell-dune-ceiling-reproduces/) | Does the Barrier3D per-cell dune-ceiling feature reproduce the in-memory experiment, and change nothing when off? | Yes on both. | record |
| [`2026-09-28-adoption-end-to-end`](2026-09-28-adoption-end-to-end/) | Does the default runner on the adopted setup reproduce the storm and dune experiments? | Yes. | record |
| [`2026-09-29-default-storms-sandbag-fix`](2026-09-29-default-storms-sandbag-fix/NOTE.md) | Do the fixes to CASCADE's default storm file name and the `sandbag_management_on` broadcast change any Hatteras result? | No: the 2010 full-management run is bit-identical. The defaults run again, and the default-storms guard fires again. | record |
| [`2026-09-29-cruft-cleanup`](2026-09-29-cruft-cleanup/NOTE.md) | Does removing the cascade/ cruft (empty `res_manager.py`, local debug prints, formatting drift via black) change any Hatteras result? | No: the 2010 full-management run is bit-identical, and the defaults still run. | record |

**Status** — **current**: its answer is in use now. **superseded**: a later study
re-asked it; follow the pointer. **record**: a finished check or a result from an
earlier set-up (÷10 offset, Hs 2.5 calibration), kept so the number can be traced.

Back to [the map](../README.md).
