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
    3-on-vs-off/               the figure scripts: module on and off, each test; figures on D: in 4-hindcast/<question>/
5-forecast/                not started
```

| Script | Folder | Output tag |
|---|---|---|
| `HAT_forced_dune_crests.py` | `1-observations/5-forced-dune-crests/` | `2026-10-26-forced-dune-crests` (holds the control run the others use) |
| `HAT_collision_rule_comparison.py` | `3-module-tests/1-rule-comparison/` | none: offline, written to the workspace's `dune-collision/3-module-tests/1-rule-comparison/` (its `NOTE.md` lists every collision test) |
| `HAT_dune_collision.py` | `4-hindcast/1-calibration-1996-2009/` | `2026-10-26-dune-collision` |
| `HAT_collision_calibration.py` | `4-hindcast/1-calibration-1996-2009/` | `2026-10-29-collision-calibration` |
| `HAT_collision_test_2009.py` | `4-hindcast/2-test-2009-2025/` | `2026-10-29-collision-test-2009-2025` |
| `HAT_collision_toe_comparison.py` | `3-module-tests/1-rule-comparison/` | none: offline, `tables/toe_comparison.csv` there |
| `HAT_collision_on_off_position.py`, `HAT_collision_model_steps.py` | `4-hindcast/3-on-vs-off/` | none: figures from the runs above, in the workspace's `dune-collision/4-hindcast/<question>/` (index: its README.md) |

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
4 Cs (Rhigh - toe)^2 t / T, taken as a crest lowering over the 20 m dune; later storms that year see the lowered
dunes. Barrier3D's own overwash lowering is untouched. Switched on in `../hatteras_ms/hat_run.yaml`
(`dune_collision`): Larson alone, Cs 1.9e-5, since 2026-11-04. Before that the volume was also scaled by
max(0, 1 - W / W0) for the beach in front of the dune (the beach factor): Cs 1.02e-4 with W0 25 m on the dry beach
(2026-11-03), and Cs 1.24e-4 with W0 45 m on the width to the crest. The beach factor was removed from the module on
2026-11-04; the last version with it is git tag `collision-with-beach-factor`. Toes from `../input_prep/1-barrier3d-domains/3-dune-toes/`, one per
domain since 2026-11-04 (below).

## Settings (2026-11-04, frozen)

> **Update 2026-11-11.** Superseded in part by the tests of 2026-11-09 to 11
> (`4-hindcast/TESTING_2026-11-09_to_11.md`): the toe is now the model's berm (1.7 m NAVD88, Laura's suggestion),
> which fits as well as the lidar toes with no toe input, and Cs about 1e-6 to 2e-6 fits it best island-wide; the road
> manager's first-year bulldozing, which put the starting topography on the dunes, is fixed (the default since
> 2026-11-09), so every earlier run has the old step; lateral erosion is under test. These have not yet been rerun
> with the road rule and Rodanthe 90%, nor on 2009-2025, and `hat_run.yaml` still sets the settings below.

The hindcast's collision setup, decided after the model-change steps below; not to be varied further without a
reason:

- the module: Larson et al. (2004) alone, per storm, no beach factor, Cs 1.9e-5, off at Cape Point (GIS 1-5);
- one toe per domain: the mean of its screened lidar toes in the surveys nearest the period's start (1996 and 1997;
  2004, 2005 and 2009, nothing after the start);
- the 198 m road-dune rule;
- Rodanthe village (GIS 77-83) returns 90% of its overwash to the shoreface, the other villages 40%
  (`site_layer/hatteras_site_config.py`, `zone_overwash_filter_pct`).

Why Larson alone (2026-11-04). Once the survey crests were re-read in retreat-only windows (2026-11-03), the beach
factor no longer earned its place. The old fixed windows read the beach as the crest where the dune had moved
landward out of the window, at exactly the narrow-beach domains (GIS 11: 3.83 to 1.68 m across Dennis, against
3.83 to 3.18 m now), which made beach width look like the control on collision loss. With the re-read crests,
offline, Larson alone and Larson x the beach factor tie (CV RMSE 0.44 m each, no erosion 0.45; Dennis r 0.49 each).
In the model the factor (domain, W0 25 m, Cs 1.02e-4) beat Larson alone by 0.02-0.06 m/yr on shorelines and
0.03-0.11 m on crests, at the cost of a beach-width input, a fitted limit and its own justification.

Cs 1.9e-5 is Larson alone's offline fit to the 1996-2009 lidar crest change (`3-module-tests/1-rule-comparison`,
`larson_storm`; 2.0e-5 on the old crests too). In the model (`4-hindcast/3-on-vs-off/HAT_collision_cs_sweep.py`,
Cs 1e-5 to 8e-5) more Cs is worse everywhere; 1e-5 and 1.9e-5 tie on the calibration crests (0.89, 0.90 m) and 1e-5
is 0.03 m/yr better on its shoreline. 1.9e-5 is kept: it has a source outside the model and gets the mean 1996-2009
crest right (bias +0.04 m; +0.15 at 1e-5). It is about ten times below Larson et al.'s published 1.3-1.7e-4
because the module holds each storm's peak runup over the storm's duration; the same fit with hourly water levels
(`larson_hourly`) gives 1.1e-4, inside the published range.

Each change added one at a time (`4-hindcast/3-on-vs-off/HAT_collision_model_steps.py`), RMSE, lower is better:

| Step | Calibration shoreline (m/yr) | Calibration crests (m) | Test shoreline (m/yr) | Test crests (m) |
|---|---|---|---|---|
| Base: module off, no road rule, villages 40% (calibration growth 0.01-0.03) | 1.14 | 0.91 | 2.24 | 1.41 |
| + dune growth 0.02-0.04 (calibration only) | 1.14 | 0.96 | | |
| + collision module (Larson alone, Cs 1.9e-5) | 1.19 | 0.91 | 2.48 | 1.60 |
| + 198 m road rule | 1.20 | 0.90 | 2.29 | 1.61 |
| + Rodanthe 90% | 1.21 | 0.90 | 2.16 | 1.61 |

What it says: in calibration the module takes back the crest error the faster dune growth adds (0.96 to 0.90 m) and
costs 0.05 m/yr of shoreline; in the test it leaves the dunes too low (bias -0.7 m), and the shoreline gain over the
base is the road rule and Rodanthe 90%. With the beach factor the same steps ended at 1.19 / 0.87 and 2.10 / 1.50
(`4-hindcast/_superseded/beach-factor-steps-until-2026-11-04/` in the workspace). Whether the module should stay
on is the open question.

The toes: per cell, the median over every survey of the period before 2026-11-04 (in 2009-2025, the 2014-2019
surveys). The start profiles' own toes fail in calibration (0.3 m below the 1996 and 1997 survey grids, 0.5 m at
Avon; the start is the 1 m survey, the others 3 m and 9 ft grids), and a low toe over-erodes (the erosion goes with
the square of the water above it). Offline (`3-module-tests/1-rule-comparison/HAT_collision_toe_comparison.py`)
domain toes fit the 1996-2009 crest change as well as per-cell toes (CV RMSE 0.443 against 0.439 m; no erosion at
all 0.452).

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
