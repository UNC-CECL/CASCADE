# 4-hindcast — the module in real model runs

- `1-calibration-1996-2009/`: the module's Cs and beach limit, then the dune growth rate, Hs and groin
  strength refit with it on. Best fit (2026-10-29): Hs 3.0, growth 0.02-0.04, groin ~0.16.
- `2-test-2009-2025/`: the 1996-2009 calibration run unchanged on 2009-2025.
- `3-on-vs-off/`: shoreline position change and fit as each model change is added, both windows
  (`HAT_collision_model_steps.py`), the Cs sweep for Larson alone (`HAT_collision_cs_sweep.py`), and the module on
  against off (`HAT_collision_on_off_position.py`).

The settings these settled on (2026-11-04): `../README.md`, Settings.

**2026-11-09 to 11:** the toe tests (Automorph, per cell, at the berm), the Cs sweeps on them, the road manager's
first-year bulldozing (found and fixed), and the new module options (row modes, lateral erosion, collapse tracking):
what, why, the numbers and the figures in `TESTING_2026-11-09_to_11.md`.
