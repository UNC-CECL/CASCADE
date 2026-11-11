# 3-module-tests — which rule, and what the module does on its own

No hindcast runs. Numbered from simplest to most complex, as `../../hard-structures/groin/3-module-tests/`.
Each script's header says what it runs and how; figures and their captions (`supporting/CAPTIONS.md`) are
in the workspace twin, `D:\research\UNC\hatteras_ch1\model_runs\dune-collision\3-module-tests\<folder>\figures\`.

```
1-rule-comparison/   offline against the lidar crest change: which collision rule fits best
2-long-term/         the module on Barrier3D's own test domain for 200 years: does it level off?
3-sensitivity/       how much the crest loss depends on Cs and the beach limit W0 (GIS 12, 1999)
4-bounds/            the module's stated limits against every storm the hindcast runs
```

- `1-rule-comparison/`: offline against the lidar crest change, which collision rule fits best.
- Unit tests, `tests/test_dune_collision.py` (repo root): the beach factor, the collision-cell test,
  Larson's volume by hand, the lowering never below the toe nor raising a row, zero Cs and a toe above
  every storm reproduce Barrier3D exactly, Barrier3D's own overwash lowering untouched, the beach
  following the shoreline, and the run-up-loss option. Since 2026-11-03 also the module's tested range:
  `DuneCollisionOutsideTestedRange` for Cs outside 0-2.5e-4, a toe below MHW, and, for a storm that
  lowers dunes, more than 200 h, a period outside 5-20 s or water more than 5 m above a toe (41 tests).
  Since 2026-11-10/11 (56 tests): the row modes (`rows`: "front" and "equal" by hand, both taking exactly the eroded
  volume, never below the toe nor raising a row, zero volume exact; an unknown mode rejected); lateral erosion
  (`erode_laterally` by hand: the front row keeps its height until its sand above the toe is gone, the rest carried
  to the back row and between storms, both rows gone stopping at the toe; `refill_deficit`: growth repays the deficit
  first; the deficit moved with the rows on dune migration, both ways; a full Barrier3D run with a valid deficit,
  pickled); and the per-cell series (`collision_volume_TS` adding up to the yearly `volume_TS`, no volume without a
  collision storm, a front-row collapse recorded, a flattened row topped up by a sliver not counted again). What the options do on Hatteras: `../4-hindcast/TESTING_2026-11-09_to_11.md`.

## 2-long-term (2026-11-03)

`collision_long_term_test.py` (about 2.5 min). Barrier3D alone on its default test domain
(`tests/test_coupling`, 50 cells, its 1000-year Virginia Coast Reserve storm series, peak water up to 4.0 m
above MHW), 200 years, module off / on with no beach / on behind a 25 m beach (W0 45 m), Cs 1.24e-4, each year
seeded alike. Dune migration with the shoreline is switched off, so only growth, overwash and collision act
(with it on, a one-cell shoreline retreat replaces the dunes with the low interior, and slow growth never
rebuilds them, module or not). Three domains: A, the test domain as it is (dune ceiling 1.5 m above the berm,
its own growth 0.55-0.95, toe 0.1 m above the berm); B, Hatteras-like (ceiling 6 m, dunes starting 4 m above
the berm, growth 0.02-0.04, toe 0.6 m above the berm, a Hatteras toe of 2.4 m NAVD88); C, B with A's low toe.
The equilibrium is the mean crest over years 100-190, before the series' largest storm (year 193, 4.04 m MHW
for 120 h).

| domain | module off | on, no beach | on, 25 m beach |
|---|---|---|---|
| A (low dunes, fast growth) | 1.45 m | 1.31 m (collision 0.064 m/yr) | 1.40 m (0.032 m/yr) |
| B (Hatteras-like) | 5.77 m | 3.20 m, still rising +0.68 m over the window | 5.06 m, rising +0.32 m |
| C (B, low toe) | 5.77 m | 0.07 m: at the toe by year ~33 | 0.07 m: at the toe by year ~90 |

- **The module does not run away on its own.** With fast growth (A) it takes 0.03-0.06 m/yr and the
  dunes sit 0.05-0.15 m below the module-off run, overwash setting the rhythm in all three.
- **With slow, Hatteras growth the outcome turns on the toe.** At a Hatteras toe (B) the dunes find a
  balance well above it (3.2 m with no beach, 5.1 m behind 25 m of beach, against 5.8 m without the
  module). Half a metre lower (C), collision outpaces growth and the dunes fall to the toe and stay there:
  Barrier3D's logistic growth is slowest on low dunes, so it cannot bring them back.
- **One large storm undoes it.** The year-193 storm (4.04 m for 120 h) takes B's unprotected dunes from
  ~3.8 m to 0.5 m (0.17 m by year 200), and the beach-backed ones to 3.7 m; slow growth has not rebuilt
  either by the end. Module off, B's dunes stand above that storm and do not move.
- **The limits work.** On the storm series as it is, the module stops the run in year 44 on a storm with a
  21.3 s period (`DuneCollisionOutsideTestedRange`); the runs above cap the periods at 20 s.

## 3-sensitivity (2026-11-03)

`collision_sensitivity_test.py` (seconds). GIS 12's 50 cells, 1999, offline through the module's own
functions (`../2-module-explained/collision_year_gif.py`'s `cells()` and `run_year()`, imported unchanged), at
Cs 0-2.5e-4 and W0 10-60 m on the dry beach (beach start to toe; GIS 12's dry beaches are 0-20 m), plus Larson's
laboratory 1.4e-3 (outside the bounds); the lidar's 1998-1999 mean crest change on the same cells (1.36 m) as
the reference. Rerun on the final calibration (Cs 1.02e-4, W0 25 m, `hat_run.yaml` since 2026-11-03).

- The year's mean lowering rises almost linearly with Cs until the lowest dunes reach their toes:
  0-1.83 m over Cs 0-2.5e-4 at W0 25 m. The hindcast's Cs 1.02e-4 and W0 25 m give 1.07 m, about 80% of
  the lidar's 1.36 m.
- W0 matters about as much as Cs on these narrow beaches: 0.64-1.25 m over W0 10-60 m at the hindcast's Cs.
- Larson's laboratory value would take 3.02 m, most dunes to their toes.
- History: on the crest-based widths (until the morning of 2026-11-03) the setting then in use, Cs 1.24e-4
  and W0 45 m, gave 0.84 m here; the widths moved to the dry beach and Cs and W0 were refit offline and
  checked in-model by another session (`cal_g16_nocape_rd198_dry`, `../4-hindcast/1-calibration-1996-2009/`).

## 4-bounds (2026-11-03)

`collision_bounds_figure.py` (seconds). All 317 storms of the 1996-2025 hindcast storm series against the
module's limits, the water taken above the lowest toe of the start each storm runs from (1.02 / 1.01 m
NAVD88), the most any cell sees.

- Every storm is inside the tested range, by little on two limits: the longest is 193 h (96% of 200 h);
  the highest, Isabel (5.58 m), puts the water 4.55 m above the lowest toe (91% of 5 m), though it overtops
  most dunes, where the module does not check it. Periods span 7.0-18.2 s, well inside 5-20 s.
- A longer storm series (a forecast) is the likeliest to reach a limit, first on duration.

## Still to do

- **The sand budget**, once the eroded sand goes somewhere (it is deleted now; see the top README's Next):
  sand removed from the dunes against sand placed on the beach, over a long run.
- The long-term test on a real Hatteras domain with dune migration on and the management modules, as the
  groin's `8-real-planform/`.
