# 2026-10-08-drift-groin

**Question.** Does the drift-blocking groin (option 1) fit the 1996-2009 gap across the groins in the full model, and with what blocking fraction b and failure fraction f?

**Why.** The straight-coast tests (`../../../2-module-tests/1-straight-coast/`) showed the pinned groin blocks BRIE's diffusive exchange, which misses the background drift: it traps nothing on a grid-aligned coast and on the wrong side on coasts turned with the drift. The drift groin reads the drift the waves actually carry. It leaves BRIE's sand leak at the bends it makes uncorrected, because the correction (module 6) also changes BRIE's behaviour at Buxton's natural bend (`../../../2-module-tests/3-real-planform/corrected_groin/`). BRIE and Barrier3D are not changed.

**Setup.** `drift_groin.py` drives the unchanged runner with `BlockingGroinCallback` replaced in a child process. Each year the groin reads the angle of the face just updrift of GIS 5|6, looks up the wave climate's net drift there from BRIE's own tables, and moves `b · Q_net · dt / (dy · depth)` from GIS 5 to GIS 6 through `x_s_dt`. The amounts are equal and opposite, so `applied_net_m` is 0. Everything else is the 10-05 setup: full management, the solved edgeBE ends, relocations off, instant failure at the 2004 step (b falls to b·f). Each run dir has `drift_groin.txt`, because the runner's report and run name (`groinblock`) do not show the patch. The diagnostics CSV gains `net_drift_m3_yr` and `approach_angle_deg`.

Grid: b 0.02-0.6 × f 0.1-0.8 (56 runs) under `output/raw_runs/experiments/groin/2026-10-08-drift-groin/`. b 0.6 moves about 56 m/yr at Buxton, about ten times the pin, so the grid reaches low.

**Scores.** `grid_scores.csv`: the 10-05 photo-date RMSE and the 10-08 annual CoastSat RMSE of the GIS 5|6 gap change, the groin's applied budget, seaward change at GIS 3, 4, 7, 8, and `island_vs_none_m`, the island's summed seaward change against the no-groin run (BRIE's leak at the groin's bends). The pinned arm (10-05 runs) and the no-groin baseline are scored beside it.

    python drift_groin.py smoke [--b 0.1 --f 0.6]
    python drift_groin.py grid [--streams 4]
    python drift_groin.py score
    python drift_groin.py test --b <b> --f <f>      # the chosen cell on 2009-2025

## Smoke cell, b 0.1 f 0.6 (calibration)

| run | annual RMSE | photo RMSE | island vs no groin | net applied |
|---|---|---|---|---|
| no groin | 41.6 | 69.2 | 0 | — |
| pinned b 0.6 f 0.6 (the 10-05 pin) | 29.2 | 4.1 | −22 m | +186 m |
| drift b 0.1 f 0.6 | 21.5 | 9.7 | +113 m | 0.0 m |

Drift at the groin 0.92 M m³/yr southward throughout. The approach face turns from −9.5° to −5.2° as the fillet grows. The groin moves 9.6 m/yr until the 2004 step, then 5.7 m/yr.

## Results (grid, 2026-10-08)

All 56 runs clean (`grid_run.log`, `score.log`). Joint score = root mean square of the two RMSEs.

| run | annual RMSE | photo RMSE | joint | gap change at the end, model / CoastSat | island vs no groin |
|---|---|---|---|---|---|
| no groin | 41.6 | 69.2 | — | −88 / −28 | 0 |
| pinned b 0.6 f 0.6 (the 10-05 pin) | 29.2 | 4.0 | 20.9 | −18 / −28 | −22 m |
| best pinned (b 0.4 f 0.8) | 10.7 | 26.9 | 20.5 | | −9 m |
| **drift b 0.1 f 0.4** | **19.0** | **11.7** | **15.8** | −23 / −28 | +107 m |
| drift b 0.05 f 0.8 (best annual) | 13.7 | 37.6 | 28.3 | −46 / −28 | +46 m |

- **The drift groin fits a little better than the pin, at about a sixth of its blocking fraction.** b 0.1 blocks 0.09 M m³/yr of the 0.92 M m³/yr arriving. The best cells all sit at b 0.1; f matters little there (joint 15.8-17.7 for f 0.1-0.6).
- **The fit is steep in b.** Between b 0.05 and 0.15 the annual RMSE goes 14 → 17 → 41 m. The grid is too coarse there; a refinement at b 0.06-0.14 is the next step.
- **The two scores still pull apart**, as for the pin: the annual CoastSat score prefers a weaker groin (b 0.05) and the three photo dates a stronger one (b 0.1-0.15).
- **Refinement, b 0.06-0.14 × f 0.2-0.8 (32 runs, `fine_run.log`).** The joint optimum does not move: a flat floor of 15.8-16.7 m from b 0.08 f 0.8 to b 0.10 f 0.4 (b 0.09 f 0.6: annual 16.9, photos 14.8, joint 15.9). Underneath it the two scores want different groins. The annual CoastSat score is best at b 0.06-0.07 with f 0.6-0.8 (9.9 m; photos 25-31 m). The photo score is best at b 0.12 f 0.4 (1.8 m; annual 30 m). No single (b, f) satisfies both; the joint cell is a compromise between them, as it was for the pin.
- **Test, 2009-2025, b 0.1 f 0.4 (2026-10-08).** Annual RMSE 26.2 m, photo RMSE 31.7 m (8 dates), against no groin 23.5 / 31.1 and the pin 31.0 / 36.3. Out of sample the drift groin does better than the pin but no better than no groin. It holds the gap too wide: the model's gap change at the end is −17.5 m against CoastSat's −50.4 m (no groin −57.8 m). The observed gap kept closing after 2009, faster than a groin at b·f = 0.04 lets it. This is the failure-schedule question already open in `../2026-10-08-schedule-refit/`, not specific to this module. The test's island change against no groin is +53 m. The pinned test run reads −416 m; that run predates this study and may not share the baseline's setup, so check it before quoting.
- **BRIE's leak grows with the groin, then saturates.** Island change against no groin is +19 m at b 0.02, +107 m at the best cell, and about +200 m from b 0.2 up: about 1.1 M m³ over the 13 years at the best cell. The pin's is −22 m.
