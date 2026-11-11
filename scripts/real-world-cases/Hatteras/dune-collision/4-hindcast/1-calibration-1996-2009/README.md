# 1-calibration-1996-2009 — the module fitted on the 1996-2009 lidar and CoastSat

`HAT_dune_collision.py` and `HAT_collision_calibration.py`. Managed edgeBE 1996-2009 with the groin, as the
headline matrix run, driving the unchanged runner (`../../../hatteras_ms/HAT_hindcast.py`) in-process.
Each script's header lists its actions and variants. Runs: `output/raw_runs/dune-storm-response/`
`2026-10-26-dune-collision/` (HAT_dune_collision.py) and `2026-10-29-collision-calibration/`
(HAT_collision_calibration.py scores the `cal_*` variants filed under the first).

`HAT_dune_collision.py` is also the hub the other folders use: its `toes`, `twl`, `calibrate` and
`empirical` actions are the measurements in `../../1-observations/`, and
`../../3-module-tests/1-rule-comparison/` imports it for the calibration cells. It stays in one piece
here until it is split (top README, Next).

The record below is the folder README as it stood on 2026-11-03, moved here unchanged. The year-end
rules (`width_*`, `larson_*`) ran outside the model and came before the module; the `storm_*` and
`cal_*` runs use `cascade/dune_collision.py`.

## The record (from 2026-10-27)

Interior LRR bias / RMSE against CoastSat (m/yr), managed edgeBE 1996-2009:

| Run | Bias | RMSE |
|---|---|---|
| control | -0.445 | 1.276 |
| forced to lidar crests, 1998, 2000, 2004 (`events`) | -0.865 | 1.578 |
| raise-only rebuild to 3.0 m NAVD88 on every road domain (`floor_3p0_anyroad`) | -0.366 | 1.241 |
| Larson collision, Cs 1.7e-4, sand to the beach (`larson_beach`) | -0.991 | 1.836 |

- Realistic (lower) dunes give more overwash, and the shoreline then retreats
  more than CoastSat shows, most at Rodanthe (GIS 80-87).
- Larson collision erosion puts the dune lowering in the right places (crest
  r 0.53 -> 0.71 in 2009) but lowers dunes on wide beaches where the lidar shows
  none, and flattened dunes never recover (p10 ~1.8 m against 2.4-3.1 m).
- Long, de Bakker & Plant (2014, GRL 41, 2899-2906) give no collision equation:
  the crest changes by less than ~10%, up to ~50% just below overtopping.
- In the Hatteras lidar, crest loss in collision is set by beach width, not
  R_high/D_high: dz/D_high = 0.50 max(0, 1 - W / 45 m), fit on Dennis and
  2004-05 (r 0.48 cells, 0.63 domain medians; R_high/D_high alone 0.29).
  `HAT_dune_collision.py empirical`, figure `collision_empirical.png`.
- Beach-width rule runs (2026-10-27/28), sand to the beach unless noted.
  Crest RMSE is the island median crest against the lidar over the surveys
  (m); stable is the 2009 median crest at GIS 20-22 and 59-62 less the lidar's (m):

  | Run | Bias | RMSE | GIS 80-89 | Overwash | Crest RMSE | Stable |
  |---|---|---|---|---|---|---|
  | control | -0.445 | 1.276 | 1.072 | 6521 | 0.197 | -0.27 |
  | every year the water reaches the toe (`width_beach`) | -1.056 | 1.680 | 1.147 | 16245 | 0.868 | -0.56 |
  | storm years only, R_high >= 3.1 m (`width_storms_beach`) | -0.568 | 1.425 | 0.886 | 9813 | 0.502 | -0.37 |
  | storm years, growth r 0.03-0.05 (`width_storms_r04`) | -0.453 | 1.352 | 0.878 | 8120 | 0.273 | +0.35 |
  | storm years, growth r 0.05-0.07 (`width_storms_r06`) | -0.392 | 1.296 | 0.890 | 7095 | 0.158 | +1.03 |
  | storm years, growth r 0.07-0.09 (`width_storms_r08`) | -0.346 | 1.269 | 0.972 | 6151 | 0.205 | +1.02 |
  | storm years, growth r 0.09-0.11 (`width_storms_r10`) | -0.294 | 1.254 | 1.084 | 5171 | 0.214 | +1.63 |

  Every year over-erodes (collision reaches the toe in 11 of 13 years). Storm
  years only gives the best crest pattern (r 0.82 against 0.53 in 2009) and the
  best Rodanthe, but at 0.01-0.03 the island median crest falls to 4.2 m by
  2009 (lidar 5.4). Refitting the growth rate restores the median (r 0.05-0.07
  best) but no single logistic rate fits both: growth scales with height, so
  the tall stable dunes overshoot by ~1 m while the collision-lowered dunes
  barely regrow (38% of cells below 4 m in 2009 against the lidar's 27%).
  Buxton (GIS 8-15) is worse than the control in every collision run (RMSE 1.29-1.47 against 1.04).

- The rules above lower the dunes after each model year, outside the model.
  Process-based since 2026-10-28: `cascade/dune_collision.py`
  (DuneCollisionManager) lowers, inside Barrier3D's storm loop and without
  editing Barrier3D, each cell a storm leaves in collision (lidar toe < peak
  Rhigh < crest) by Larson's rule per storm, 4 Cs (Rhigh - toe)^2 t / T / 20 m;
  later storms that year see the lowered dunes. Barrier3D's own overwash
  lowering is untouched (checked: identical on an overwash storm). Cs refit in
  the per-storm form (the peak Rhigh over the duration gives ~12x the hourly
  exposure): island 2.0e-5 (r 0.28, background growth +0.058 m/yr), Dennis
  4.2e-5. Sand lost.

  | Run | Bias | RMSE | GIS 8-15 | GIS 80-89 | Overwash | Crest RMSE | Stable |
  |---|---|---|---|---|---|---|---|
  | control | -0.445 | 1.276 | 1.037 | 1.072 | 6521 | 0.197 | -0.27 |
  | per storm, Cs 2.0e-5 (`storm_cs2e5`) | -0.786 | 1.603 | 0.971 | 0.833 | 10864 | 0.465 | -1.22 |
  | per storm, Cs 4.2e-5 (`storm_cs4e5`) | -1.028 | 1.759 | 0.872 | 1.137 | 13901 | 0.913 | -2.30 |
  | per storm, Cs 2.0e-5, growth 0.05-0.07 (`storm_cs2e5_r06`) | -0.676 | 1.479 | 0.990 | 0.811 | 8717 | 0.205 | +0.16 |
  | per storm, Cs 4.2e-5, growth 0.05-0.07 (`storm_cs4e5_r06`) | -0.923 | 1.678 | 0.866 | 1.022 | 11795 | 0.723 | -1.09 |

  `storm_cs2e5_r06` is the first run better than the control at both Buxton
  and Rodanthe, with the island median crest on the lidar (5.40 against 5.41 m
  in 2009) and the stable dunes within 0.16 m: with collision in the model,
  growth 0.05-0.07 no longer overbuilds them. Its island RMSE is worse, from
  the cape (GIS 2-7: 3.20 against 2.17) and the centre (GIS 26-75), which the
  control already erodes too much (bias -0.84). Cs 4.2e-5 brings GIS 10-12 to
  CoastSat but erodes everywhere else too much.
- Committee test (2026-10-28): does a sink where the model's dunes do not erode
  make Dennis lower them? A sink at GIS 9-14 (`BUXTON_SINK`, CoastSat LRR less
  the control's, -0.35 to -1.92 m/yr, through `HAT_BE_OVERRIDE`), without
  (`sink`) and with per-storm collision (`storm_cs2e5_r06_sink`):

  | Run | Bias | RMSE | GIS 8-15 | GIS 80-89 | Crest RMSE | Stable |
  |---|---|---|---|---|---|---|
  | control | -0.445 | 1.276 | 1.037 | 1.072 | 0.197 | -0.27 |
  | `sink` | -0.547 | 1.248 | 0.515 | 1.072 | 0.167 | -0.27 |
  | `storm_cs2e5_r06` | -0.676 | 1.479 | 0.990 | 0.811 | 0.205 | +0.16 |
  | `storm_cs2e5_r06_sink` | -0.783 | 1.461 | 0.611 | 0.811 | 0.231 | +0.16 |

  No for Dennis: the 1999 collision lowering at GIS 8-15 is identical with and
  without the sink, and the 1999-2000 crests differ by < 0.3 m. The sink does
  lower the dunes, but in 2003-04, once the shoreline has retreated a full
  10 m cell: Barrier3D drops the front dune row and fills the back row from the
  lower interior (GIS 12-13: 7.3 -> 4.6 m against the lidar's ~4 m from 1999,
  and GIS 14 to 4.5 m where the lidar keeps 5.6-6.1 m). The shoreline gain is
  the sink fitted to CoastSat, not a dune result (GIS 12-13 now overshoot,
  -3.0 against -2.0 m/yr). Crests per survey: `sink_crests_gis8_15.csv`.
- With the beach (2026-10-28): `DuneCollisionManager(beach_width=..., beach_width_limit=...)`
  scales Larson by max(0, 1 - W / 45 m), W each cell's median beach width over the
  surveys (1998 left out), held for the run. Offline, against the lidar crest
  change (domain medians, fit intervals): r 0.28 -> 0.48, Dennis 0.63 -> 0.86,
  and Cs 2.0e-5 -> 1.24e-4, on Larson et al.'s (2004) value for US East Coast
  storms (Birkemeier et al. 1988 data: 1.3e-4 fit, 1.7e-4 mean per case; their
  laboratory and short field sets give ~1.4e-3). A width from the
  survey before each storm fits worse (r 0.38).

  | Run | Bias | RMSE | GIS 8-15 | GIS 80-89 | Cape 2-7 | Centre 26-75 | Crest RMSE | Stable | Crest r 1999 / 2009 |
  |---|---|---|---|---|---|---|---|---|---|
  | control | -0.445 | 1.276 | 1.037 | 1.072 | 2.17 | 1.29 | 0.197 | -0.27 | 0.75 / 0.53 |
  | `storm_cs2e5_r06` | -0.676 | 1.479 | 0.990 | 0.811 | 3.20 | 1.45 | 0.205 | +0.16 | 0.82 / 0.72 |
  | `storm_w45` | -0.787 | 1.483 | 0.876 | 1.099 | 3.36 | 1.37 | 0.305 | -0.27 | 0.87 / 0.77 |
  | `storm_w45_r04` | -0.707 | 1.405 | 0.831 | 1.073 | 3.29 | 1.26 | 0.122 | +0.41 | 0.87 / 0.78 |
  | `storm_w45_r06` | -0.671 | 1.362 | 0.838 | 1.011 | 3.23 | 1.20 | 0.203 | +1.08 | 0.86 / 0.76 |

  Growth 0.03-0.05 (`storm_w45_r04`) is the dune fit: island median crest
  5.65 / 5.30 / 5.43 m in 1999 / 2004 / 2009 against the lidar's 5.53 / 5.38 /
  5.41, the best crest RMSE and crest pattern of any run, stable dunes +0.41 m.

  The beach puts the lowering in the right places (best crest pattern of any
  run) and the centre now beats the control. Not yet: GIS 10-12 drop in 1999
  but not far enough (5.1-6.4 m against the lidar's 1.7-4.0 m); the cape is no
  better, so collision is not what is wrong there; and neither growth rate
  fits, 0.01-0.03 leaves the island crest low, 0.05-0.07 overbuilds the stable
  dunes (behind wide beaches, now untouched by collision) by 1 m.
- Dennis at Buxton (2026-10-28, `storm_w45_r04`): the forcing is not short.
  The storm series holds Dennis's 3.65 m peak for 116 h; the hourly TWL is
  above 2.5 m for ~52 h, so the model gets ~6x the hourly record's exposure.
  The 1999 shortfall at GIS 10-12 (model lowers 0.6 / 1.3 / 0.9 m, lidar
  1.7 / 2.1 / 1.6 m) has two parts:
  1. The model enters Dennis too tall: the 1998 survey (flown 1998-09) is after
     Bonnie, which had already lowered GIS 11-12 by 1.3 / 0.8 m since 1997;
     the model's Bonnie (2.64 m, 64 h) barely reaches the toes (2.4-2.8 m), and
     GIS 12 stands at 7.3 m against the lidar's 5.9 m before Dennis.
  2. Neighbours with the same toe (2.4 m), beach (20-30 m) and crest respond
     differently: GIS 9 lost 0.1 m, GIS 10 1.7 m. No rule on these inputs can
     tell them apart, so the fit over-lowers GIS 9 (0.7 m) and under-lowers
     GIS 10-12.
- From 2026-10-28 `Hatteras-CASCADE-parameters.yaml` sets `StormStart: 1`
  (Hannah: storms always run from 1997). Every run above used 2. Effect, small
  (`control_ss1`, `storm_w45_r04_ss1`): bias -0.424 / -0.691, RMSE 1.279 /
  1.396, GIS 8-15 1.036 / 0.833, GIS 80-89 1.117 / 1.016 (against -0.445 /
  -0.707, 1.276 / 1.405, 1.037 / 0.831, 1.072 / 1.073 with 2). The 1997 storms
  (peak 2.45 m) lower the island's dunes 0.02 m on average and leave Buxton
  before Dennis unchanged: they do not explain the 1997-98 lidar loss.
- Hs 3.0 m (advisors, 2026-10-28; `control_hs3`, `storm_w45_r04_hs3`, StormStart 1,
  edgeBE ends still those solved at Hs 2.0). Hs is CASCADE's wave_height (BRIE
  alongshore diffusivity, shoreface depth 8.9 Hs = 26.7 m); the storm series,
  so the collision volumes, are unchanged. Shoreline RMSE (m/yr), Hs 2.0 -> 3.0:

  | Run | Island | Cape 2-7 | GIS 8-15 | Centre 26-75 | GIS 80-89 |
  |---|---|---|---|---|---|
  | control | 1.279 -> 1.174 | 2.17 -> 1.67 | 1.04 -> 1.13 | 1.28 -> 1.16 | 1.12 -> 1.32 |
  | beach collision, r 0.03-0.05 | 1.396 -> 1.224 | 3.27 -> 2.44 | 0.83 -> 0.94 | 1.25 -> 1.16 | 1.02 -> 0.97 |

  Dunes unchanged (crest RMSE 0.126, stable +0.42 m). Hs 3.0 helps where
  alongshore transport rules (cape, centre); with it the collision run is
  within 0.05 m/yr of the control island-wide and beats it at Buxton and
  Rodanthe. To confirm once the ends are re-solved at Hs 3.0.

  Ends re-solved at Hs 3.0 with StormStart 1 (`be_end_solve_net_change.py --hs 3`,
  `end-domain-boundaries-1996-2009/2026-10-28-ends-solved-on-obx-hs3-1996_2009`,
  step 8): GIS 1 +3.07, GIS 90 +2.06 m/yr (edgeBE at Hs 2.0: +3.94, +1.21).
  Runs `control_hs3_ends` / `storm_w45_r04_hs3_ends`: island RMSE 1.185 / 1.235,
  cape 1.74 / 2.51, GIS 8-15 1.13 / 0.94, centre 1.16 / 1.16, GIS 80-89
  1.35 / 0.98; within 0.07 m/yr of the runs on the Hs 2.0 ends. The groin face
  still over-erodes at Hs 3.0 (GIS 6 -3.59 against CoastSat -1.22; -2.6 at
  Hs 2.0): the groin's 0.10 strength, fitted at Hs 2.0, is the next thing to refit.

  Hs 2.5 (`control_hs2p5`, `storm_w45_r04_hs2p5`, on the Hs 2.0 ends). Shoreline
  RMSE (m/yr), collision run at Hs 2.0 / 2.5 / 3.0: island 1.40 / 1.29 / 1.23,
  cape 3.27 / 2.86 / 2.51, GIS 8-15 0.83 / 0.88 / 0.94, centre 1.25 / 1.19 /
  1.16, GIS 80-89 1.02 / 0.95 / 0.98, GIS 6 -2.59 / -3.04 / -3.59 (CoastSat
  -1.22). Control: island 1.28 / 1.19 / 1.19, GIS 80-89 1.12 / 1.07 / 1.35.
  Hs 3.0 wins island-wide only through the cape; 2.5 is as good or better
  everywhere else, and the groin face erodes less.

  Groin strength sweep, collision run, 0.04 / 0.06 / 0.08 / 0.10 (`storm_w45_r04_hs2p5_g*`,
  `storm_w45_r04_hs3_g*`; Hs 3.0 on its own ends). Island RMSE moves <= 0.02
  m/yr. A weaker groin makes GIS 6 worse, not better (Hs 3.0: -4.77 / -4.49 /
  -3.84 / -3.59; Hs 2.5: -4.03 / -3.92 / -3.21 / -3.04; CoastSat -1.22) and
  barely moves GIS 2-3. 0.10 is best of the four at both Hs: the groin-face
  over-erosion at higher Hs is not the groin holding back too much.

  Extended to 0.13 / 0.16 / 0.19. Island RMSE stays flat (Hs 3.0 1.23-1.25,
  Hs 2.5 1.29-1.30; the 0.02 wobbles are run noise). A stronger groin moves
  the error from GIS 6 to GIS 7: Hs 3.0 GIS 6 -3.59 -> -2.58, GIS 7 -1.52 ->
  -0.70 (CoastSat -1.22, -2.01); GIS 6-8 RMSE 1.57 -> 1.25 at Hs 3.0, 1.41 ->
  1.22 at Hs 2.5; Buxton 0.94 -> 0.91 / 0.88 -> 0.85; GIS 2-3 slightly worse.
  No optimum inside the range on the LRR; the groin's own fit metrics (gap,
  fillet, net GIS 2-12 / 4-9, `hard-structures/groin/.../2026-10-20-spread-relaxation`)
  were not rerun.
- Found on the way: the Hatteras runs keep Barrier3D's default `StormStart: 2`,
  so the 1997 storms (time 1; 5 storms, peak 2.45 m NAVD88) never run.

## 2026-11-04: the model-change steps, the road rule fix and the domain toes

- The module-off runs never got the 198 m road rule: `SHORELINE_GATE` was applied only in the module-on branch
  (`cal_control_g16_rd198` first ran without it). Fixed: `road_gate()` runs in both branches.
- `cal_step1_base` (module off, no road rule, Rodanthe 40%, growth 0.01-0.03, the rate before the module),
  `cal_step2_growth` (0.02-0.04, the 1996 default since 2026-10-29), `cal_step3_module`, `cal_step4_road`,
  `cal_step5_rodanthe90`, `cal_step6_starttoes` (per-cell toes; start-profile toes and widths in the last), then on
  one toe and one beach factor per domain from the 1996 and 1997 lidar (`TOE_DOMAIN_AVG`): `cal_step3_module_dt`,
  `cal_step4_road_dt`, `cal_step6_toeavg`; `cal_step6_toeavg_cw` with the per-cell beach widths. `RODANTHE_40`
  holds Rodanthe at 40% before step 5.
- Results and the settings adopted: `../../README.md`, Settings. The 1996 start-profile toes are 0.3 m below the
  1996 and 1997 survey grids' (0.5 m at Avon) and over-erode: crest RMSE 1.33 m against 0.71 m.

## 2026-11-04 (evening): Larson alone

- `cal_step6_toeavg_la`: `cal_step6_toeavg` with Larson alone (`LARSON_ALONE`: the beach factor 1 everywhere, Cs
  1.9e-5, its offline fit). Larson alone became the default (`hat_run.yaml` `dune_collision.beach: none`).
- `LA_CS_SWEEP`: the same at Cs 1e-5, 3e-5, 5e-5 and 8e-5 (`cal_la_cs10` ... `cal_la_cs80`), scored by
  `../3-on-vs-off/HAT_collision_cs_sweep.py`. 1.9e-5 kept.
- `cal_step3_module_la`, `cal_step4_road_la`: steps 3 and 4 with Larson alone (Rodanthe 40%), for the step figure.

## 2026-11-09 to 11: toes, the road fix, the berm toe and the row modes

What was found, and the figures: `../TESTING_2026-11-09_to_11.md`. The variants added here, all Larson alone, growth
0.02-0.04, groin 0.16 and, unless named otherwise, no road rule and Rodanthe 40% (the module-only setup):

- Automorph toes (`AUTOMORPH_TOES`, `input_prep/9-dune-toes/1-toes/automorph_domain_toes_1996.csv`):
  `cal_automorph_cs10`, `_cs19` (with the road rule and Rodanthe 90%, as `cal_step6_toeavg_la`);
  `cal_automorph_module_cs10` and `AUTOMORPH_CS_SWEEP` (`cal_automorph_module_cs30` ... `_cs200`), module only.
- Toes per cell (`automorph_cell_toes[_all]_1996.csv`): `cal_automorph_cells_9_15`, `cal_automorph_cellmed_9_15`
  (GIS 9-15 only); `AUTOMORPH_CELLS_ALL`, `cal_automorph_cellsall_cs10` ... `_cs200` (every domain).
- The road fix (`NO_FIRST_SCRAPE`, `_noscrape`): the road manager's first bulldozing kept off the dunes; since
  2026-11-09 the default in `cascade/roadway_manager.py`, so every later run has it whether named so or not.
  `cal_step2_growth_noscrape` (module off), `cal_automorph_module_cs*_noscrape`, `cal_automorph_cellsall_cs*_noscrape`.
- The toe at the berm (`TOE_AT_BERM`, 1.7 m NAVD88): `cal_toeberm_cs0p1`, `_cs0p3`, `_cs1` ... `_cs100`; lateral
  erosion `cal_toeberm_lat_cs1` ... `_cs250` (`COLLISION_ROWS`, the module's `rows` option).
- `COLLISION_ROWS`: which dune rows lose the collision volume per variant ("proportional", every variant not listed;
  "front", "equal", "lateral").
- Every run with the module now writes `logs/collision_cells_<variant>.csv` and its `tables/collision_cells.csv`
  (`collision_cell_table`: per domain and year the eroded volume, collision storms, cells in collision, front and
  back rows collapsed, lateral deficits); `python HAT_dune_collision.py collapse <variant>` draws them
  (`collision_collapse_<variant>.png`).
