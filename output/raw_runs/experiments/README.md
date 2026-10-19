# Experiments, by topic

Each study answers one question and lives in `<topic>/<date>-<what it tested>/`. Its README or NOTE is the record. Each topic's README has a table of its studies: the question, the answer, and whether it is still current. A study's folder path under `experiments/` is also its runs' tag in `run_index.csv`.

Only studies under the 2026-10-05 calibration/test plan are kept here. Everything older is on `D:\archive\CASCADE_offload\output\raw_runs\experiments\` (see `OFFLOADED.md`).

## Topics

| topic | what it covers | current answer, in short |
|---|---|---|
| [`end-domain-boundaries/`](end-domain-boundaries/README.md) | The source/sink rates at the two end domains (GIS 1 and 90), solved against each target and window. | 1996_2009 +1.4041 / +10.3439; 2009_2025 edgeBE +32.2422 / +21.5341 (re-solved 2026-10-09); 1996_2025 +144.2227 / +68.4160 m/yr, all on net change or the window's own target. |
| [`groin/`](groin/README.md) | The Cape Point groin: strength and form. | Blocking b0.6/f0.6 on 1996–2009, pinned 2026-10-05; the 10-08 schedule refit (failure from 1996, same b × f) is awaiting a decision. |
| [`source-sink/`](source-sink/README.md) | The per-domain BE field (set 1) and whether it transfers to the test. | Set 1 fits calibration to 0.49 m; the test bias is the 2021 step; set 2 not derived. |
| [`management/`](management/README.md) | Whether the management input (fills, volumes, footprints) arrives as intended and helps the fit. | Reported fill footprints kept; the CoastSat-observed footprints made the test worse. |
| [`target-end-window/`](target-end-window/README.md) | Options for the start and end lines of the calibration and test targets. | Calibration end is now calendar 2008 ±6 months (+3.6 m seaward of the DEM-centred line). Test end window moves the target 1–3 m, far short of the model's gap. A ±6-month start changes scores by about 1 m; ±1 yr kept. |
| [`code-checks/`](code-checks/README.md) | Whether a code or model change moves the results. | Records only. |

## The model as configured now

- **Plan:** calibrate 1996–2008 (ends 1 Jan 2009), test 2009–2025 (ends 1 Jan 2026, since 2026-10-09); resume document `scripts/hatteras_ms/DEM_TO_DEM_CALIBRATION.md`.
- **Waves:** option A (Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5), the `HAT_hindcast_config` / `hat_run.yaml` defaults. The study that chose them is on D: (`wave-climate/2026-09-27-wave-recommendation`).
- **End rates:** `hatteras_site_config.HATTERAS_BE_EDGE_ONLY`, so a plain edgeBE run uses them with no override. To try others, pass `HAT_BE_OVERRIDE="1=<GIS 1>,90=<GIS 90>"` with the edgeBE preset.
- **Storms:** `v3_split12_trim24`; per-cell dune ceilings and the Barrier3D overwash fixes (Barrier3D `hatteras/adopted`).
- **Old folder names** in logs and run metadata: `chains/RENAMES.md` on D:.
