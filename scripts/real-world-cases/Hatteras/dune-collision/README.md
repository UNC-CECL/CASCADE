# dune-collision — dunes lowered by storms that reach them but do not overtop them

Everything about the dune collision module (`cascade/dune_collision.py`, `DuneCollisionManager`): the
problem in the lidar, how the module works, what it does on its own, how it was fitted to Hatteras, and
later, what it implies for the future. Laid out like `../hard-structures/groin/`: the folders are grouped
by **what kind of claim** they hold and numbered in the order the work builds on itself.

Named `dune-storm-response/` until 2026-11-03. The runs keep that name: they are under
`output/raw_runs/dune-storm-response/<tag>/` on the workspace (run kind `dune-storm-response`).
Notes, data and figures that are not model runs go in the matching folders of
`D:\research\UNC\hatteras_ch1\model_runs\dune-collision\`.

```
1-observations/            the problem, in the lidar and CoastSat; no model fit
    README.md                  lidar crest loss in collision, and what sets it (beach width)
    1-mismatch-along-island/   (figures only, in the workspace) the crest offset by domain and survey
    ...                        2-hotspots, 3-mechanism, 4-drivers: HAT_collision_motivation*.py here
    5-forced-dune-crests/      the lidar's crests forced into the model: the shoreline erosion appears
    6-toes-vs-dunelines/       the lidar toes against the dune traced from the imagery
2-module-explained/        how DuneCollisionManager works, step by step: regimes, beach factor, the model year
3-module-tests/            which rule, and what the module does on its own; no hindcast runs
    README.md                  the unit tests (tests/test_dune_collision.py) and what is still to test
    1-rule-comparison/         offline: which collision rule fits the lidar crest change best
4-hindcast/                the module in real model runs
    1-calibration-1996-2009/   Cs, beach limit, dune growth, Hs and groin fitted with the module on
    2-test-2009-2025/          the calibration carried unchanged to 2009-2025
5-forecast/                not started
```

| Script | Folder | Output tag |
|---|---|---|
| `HAT_forced_dune_crests.py` | `1-observations/5-forced-dune-crests/` | `2026-10-26-forced-dune-crests` (holds the control run the others use) |
| `HAT_collision_rule_comparison.py` | `3-module-tests/1-rule-comparison/` | none: offline, written to the workspace's `dune-collision/3-module-tests/1-rule-comparison/` (its `NOTE.md` lists every collision test) |
| `HAT_dune_collision.py` | `4-hindcast/1-calibration-1996-2009/` | `2026-10-26-dune-collision` |
| `HAT_collision_calibration.py` | `4-hindcast/1-calibration-1996-2009/` | `2026-10-29-collision-calibration` |
| `HAT_collision_test_2009.py` | `4-hindcast/2-test-2009-2025/` | `2026-10-29-collision-test-2009-2025` |

The road and village dune management studies (`HAT_road_dune_rebuild.py`, `HAT_measured_dune_heights.py`)
moved on 2026-11-03 to `../hatteras_ms/experiments/dunes/`.

## The problem

The modelled dunes match the 1996 lidar (r 0.99) but stop matching it after Dennis (1999): Barrier3D
lowers a dune only where the water overtops it (Rhigh > crest), so storms that only reach the dune face,
the collision regime of Sallenger (2000), leave the model's dunes untouched while the lidar shows them
dropping at GIS 2-3, 10-12 and 84-86.

## How the module works (`cascade/dune_collision.py`)

Inside Barrier3D's storm loop, without editing Barrier3D, each cell a storm leaves in collision
(lidar toe < peak Rhigh < crest) loses Larson et al.'s (2004) volume per storm,
4 Cs (Rhigh - toe)^2 t / T, taken as a crest lowering over the 20 m dune; later storms that year see the
lowered dunes. Barrier3D's own overwash lowering is untouched. Switched on in `../hatteras_ms/hat_run.yaml`
(`dune_collision`): Larson alone, Cs 1.9e-5, since 2026-11-04. Before that the volume was also scaled by
max(0, 1 - W / W0) for the beach in front of the dune (the beach factor): Cs 1.02e-4 with W0 25 m on the dry beach
(2026-11-03), and Cs 1.24e-4 with W0 45 m on the width to the crest. The beach factor was removed from the module on
2026-11-04; the last version with it is git tag `collision-with-beach-factor`. Toes from
`../input_prep/1-barrier3d-domains/3-dune-toes/`.

## Next

1. `2-module-explained/`: done 2026-11-03, eight figures and GIFs (regimes, beach factor, the model year, a
   year of storms on Buxton and on GIS 12, GIS 12 through the full model, one storm in cross-section, Dennis
   along the island), all on the final calibration.
2. `3-module-tests/`: the module on an idealized coast, as the groin's. Done 2026-11-03: stated bounds
   (`DuneCollisionOutsideTestedRange`), `2-long-term/`, `3-sensitivity/`, `4-bounds/`; left: the sand budget
   (needs the sand to go somewhere) and a real-domain long run. The calibration moved to the dry-beach
   widths the same day (Cs 1.02e-4, W0 25 m, checked in-model as `cal_g16_nocape_rd198_dry`); growth 0.02-0.04,
   Hs 3.0 and groin 0.16 were fitted on the crest-based width and not rechecked, and the 2009-2025 test
   window has not been rerun on the dry beach.
3. Split `HAT_dune_collision.py`: its `toes`, `twl`, `calibrate` and `empirical` actions into
   `1-observations/`, the runs and scores left in `4-hindcast/1-calibration-1996-2009/`.
4. The cape and centre over-erosion, which collision adds to but does not cause.
5. Pair the module with the 3.0 m raise-only road rebuild floor (`../hatteras_ms/experiments/dunes/`).

Older Next items, done: `StormStart: 1` for the hindcasts (2026-10-28); the module wired into
`cascade.py` behind a switch, its toes and beach widths in `input_prep/` (uncommitted as of 2026-11-03).
