# Experiments, by topic

Each study answers one question and lives in `<topic>/<date>-<what it tested>/`. Its README or NOTE is the record. Each topic's README has a table of its studies: the question, the answer, and whether it is still current. A study's folder path under `experiments/` is also its runs' tag in `run_index.csv`. Folders moved on 2026-10-20 are listed in `RENAMES.md`.

Only studies under the 2026-10-05 calibration/test plan are kept here. Everything older is on `D:\archive\CASCADE_offload\output\raw_runs\experiments\` (see `OFFLOADED.md`).

**Model states deleted 2026-10-20.** Every `.npz` from runs on the old storms was deleted without a copy (`../README.md`). Only studies from 10-19 and 10-20 still have model states. The rest keep their tables, logs, figures and metadata, so their numbers can be traced but not reloaded.

## Topics

| topic | what it covers | current answer, in short |
|---|---|---|
| [`end-domain-boundaries-1996-2009/`](end-domain-boundaries-1996-2009/README.md) | The source/sink rates at the two end domains (GIS 1 and 90) for the calibration window, solved against each target. | +0.98 / +10.80 m/yr (solved 10-13, in the site config). Needs re-solving on the 10-20 settings. |
| [`end-domain-boundaries-2009-2025/`](end-domain-boundaries-2009-2025/README.md) | The same for the test window. | +29.66 / +20.76 m/yr (solved 10-13, in the site config). Needs re-solving on the 10-20 settings. |
| [`end-domain-boundaries-1996-2025/`](end-domain-boundaries-1996-2025/README.md) | The same for the 30-yr full-window run. | The 10-13 solve is unfinished; the last converged one (10-12) is +137.3 / +68.2 m/yr. |
| [`groin-wave-driven/`](groin-wave-driven/README.md) | The wave-driven groin (`GroinCallback`), the runner's only groin: strength per window, whole-life strength, face, fillet shape. | 0.16 for 1996–2009, 0.10 for 2009–2025, off from the 2017 fill; 0.25 for the whole life. Fitted on pre-10-20 settings. |
| [`groin-blocking/`](groin-blocking/README.md) | The blocking, drift and conserving groins, before the wave-driven one. | Superseded 2026-10-14. Record only. |
| [`dunes/`](dunes/README.md) | The dune ceiling, growth rate and design height, and the modelled dunes against the DEMs and surveys. | Dmaxel 12.2 m NAVD88 (Dmax 10.5 m above the berm), growth 0.005–0.01, design height 3.94 m MHW (10-20). |
| [`storms-overwash/`](storms-overwash/README.md) | The storm series and how modelled overwash is scored against the imagery. | Storms `v3_240` (notebook rule, 10-20). About a quarter to half of observed overwash is likely sound-side and out of the model's reach. |
| [`source-sink/`](source-sink/README.md) | The per-domain BE field (set 1) and whether it transfers to the test. | Set 1 fits calibration to 0.49 m; the test bias is the 2021 step; set 2 not derived. |
| [`management/`](management/README.md) | Whether the management input (fills, volumes, footprints) arrives as intended and helps the fit. | Reported fill footprints kept; the CoastSat-observed footprints made the test worse. |
| [`target-end-window/`](target-end-window/README.md) | Options for the start and end lines of the calibration and test targets. | Superseded by the 10-12 calendar-year setup. A ±6-month start changes scores by about 1 m; the test end window moves the target 1–3 m, far short of the model's gap. |
| [`revetment/`](revetment/README.md) | The simple revetment (repaired wall) on a straight coast. | Tested 10-12; module not adopted (write-up in `hard-structures/revetment/`). |
| [`code-checks/`](code-checks/README.md) | Whether a code or model change moves the results. | Records only. The groin guard and the 10-20 groin cleanup change nothing on the Buxton hindcast. |

## The model as configured now (2026-10-20)

- **Plan:** calendar-year setup since 2026-10-12. Calibrate 1996–2009 (storms 1997–2009, graded on calendar 2009), test 2009–2025 (storms 2010–2025, graded on calendar 2025). Resume document `scripts/hatteras_ms/DEM_TO_DEM_CALIBRATION.md`.
- **Waves:** option A (Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5), the `HAT_hindcast_config` / `hat_run.yaml` defaults. The study that chose them is on D: (`wave-climate/2026-09-27-wave-recommendation`).
- **Storms:** `v3_240` (`storms-overwash/2026-10-20-storm-notebook-rule`); Barrier3D `hatteras/adopted`.
- **Dunes and beach:** one Dmaxel of 12.2 m NAVD88, growth 0.005–0.01, design height 3.94 m MHW; beach slope 0.06 (`BEACH_SLOPE` in the runner).
- **Groin:** `GroinCallback`, `hat_run.yaml` `groin.strength_1996` 0.16 / `strength_2009` 0.10.
- **End rates:** `hatteras_site_config.HATTERAS_BE_EDGE_ONLY`, so a plain edgeBE run uses them with no override. To try others, pass `HAT_BE_OVERRIDE="1=<GIS 1>,90=<GIS 90>"` with the edgeBE preset.
- **Waiting on a re-run** after the 10-20 changes: the end solves, the groin fit, the matrix and the wave sweep.
- **Old folder names** in logs and run metadata: `RENAMES.md` here for 2026-10-20, `chains/RENAMES.md` on D: for older ones.
