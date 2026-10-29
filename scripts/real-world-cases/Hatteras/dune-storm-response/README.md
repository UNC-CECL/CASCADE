# dune-storm-response — how the dunes respond to storms, and how they are rebuilt

Active work (2026-10-26 on). The modelled dunes match the 1996 lidar (r 0.99) but
stop matching it after Dennis (1999): Barrier3D lowers a dune only where the water
overtops it (Rhigh > crest), so storms that only reach the dune face, the
collision regime of Sallenger (2000), leave the model's dunes untouched while
the lidar shows them dropping at GIS 2-3, 10-12 and 84-86.

The three scripts drive the unchanged runner (`../hatteras_ms/HAT_hindcast.py`),
patched in-process only: managed edgeBE 1996-2009 with the groin, as the
headline matrix run. Model output is under
`output/raw_runs/dune-storm-response/<tag>/` on the workspace (run kind
`dune-storm-response`; moved out of `experiments/dunes/` on 2026-10-29).

| Script | What it asks | Output tag |
|---|---|---|
| `HAT_forced_dune_crests.py` | If the dunes are set to the lidar crests after key storms, does the shoreline follow CoastSat? | `2026-10-26-forced-dune-crests` (holds the control run the others use) |
| `HAT_road_dune_rebuild.py` | Does a raise-only, lower post-storm dune rebuild along NC-12 change the shoreline? | `2026-10-26-road-dune-rebuild` |
| `HAT_dune_collision.py` | Collision-regime dune erosion: lidar toes, the hourly TWL, the Larson et al. (2004) and beach-width calibrations, and the runs | `2026-10-26-dune-collision` |

## Where things stand (2026-10-27)

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
  and Cs 2.0e-5 -> 1.24e-4, near Larson's laboratory 1.7e-4. A width from the
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
- Found on the way: the Hatteras runs keep Barrier3D's default `StormStart: 2`,
  so the 1997 storms (time 1; 5 storms, peak 2.45 m NAVD88) never run.

## Next

1. Put the eroded sand on the beach in the per-storm runs (the year-end runs
   did better with it), and wire DuneCollisionManager into `cascade.py` behind
   a switch, as the groin.
2. The cape and centre over-erosion, which collision adds to but does not cause.
3. `StormStart: 1` for the hindcasts.
4. Pair the chosen rule with the 3.0 m raise-only rebuild floor.
3. Once a rule is chosen: the toes, beach width and TWL go to `input_prep/`, the
   calibration's parameters to `data/hatteras_init/`, and the rule itself to a
   CASCADE module (`cascade/dune_collision.py`, as `cascade/groin.py`) switched
   on from `hat_run.yaml`.

## HAT_forced_dune_crests.py

Asked (Hannah, 2026-10-26) because the managed edgeBE 1996-2009 run misses the
CoastSat erosion at GIS ~8-15. Dennis (1999) runs up to 3.65 m NAVD88, the
modelled dunes there stand 5-7 m, and Barrier3D lowers a dune only where a
storm overtops it, so Dennis leaves them untouched while the lidar shows
GIS 10-12 dropping to ~1.7 m. Datums (MHW 0.36, berm 1.7 m NAVD88) and the
Duck gauge were checked and are not the cause: the Cape Hatteras Pier gauge
peaked lower than Duck during Dennis.

In-process, before the given storm year, each cell's two dune rows are scaled
so the crest equals the survey crest (`survey_dune_crests.csv`; beach landward
of the window counts as the berm). Managed edgeBE, groin, 1996-2009.
Variants: control, buxton (GIS 8-12, 1999), buxton_pre (8-12, 1998),
buxton_both (8-12, 1998 and 1999), buxton_8_15 (8-15, 1999), island (1-90, 1999),
events (1-90; 1998, 2000, 2004), buxton_events (8-15; 1998, 2000, 2004),
all_surveys (1-90; every survey).

    python HAT_forced_dune_crests.py run [variants...]
    python HAT_forced_dune_crests.py score
    python HAT_forced_dune_crests.py plot     # position change vs CoastSat
    python HAT_forced_dune_crests.py crests   # control vs lidar, GIS 8-15

Results 2026-10-26 (`output/raw_runs/dune-storm-response/2026-10-26-forced-dune-crests/`):
the control reproduces the matrix run (interior bias -0.445 m/yr). Forcing the
1999 crests at GIS 8-12 or 8-15 makes the erosion appear: GIS 11-12 match
CoastSat, GIS 8-10 overshoot (GIS 10 -4.0 against ~-0.9 m/yr). 1998 forcing
alone barely changes anything. GIS 13-15 do not move: their lidar crests
(5.7-6.9 m) stand above every storm's run-up, so the erosion there needs
another mechanism (background or alongshore loss). Interior bias: buxton
-0.538, buxton_8_15 -0.538, island -0.689 m/yr. Forcing every event makes it
worse (events -0.865, all_surveys -0.931; GIS 9-12 overshoot 2-4x).

The modelled and lidar crests agree through 1997 and part in 1998-1999,
when Bonnie and Dennis lowered the real dunes and not the modelled ones.
Island-wide, per run: `cascade_pipeline/plotting/dune_crest_comparison.py`
(drawn by `run_figures.py` into each run's `figures/dune_crests/`); in the
managed 1996-2009 run 17-21 domains end more than 1 m above the lidar after
1999. Not adopted: a diagnosis, not a model change.
