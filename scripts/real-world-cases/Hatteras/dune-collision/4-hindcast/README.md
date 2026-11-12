# 4-hindcast — the module in real model runs

- `1-calibration-1996-2009/`: the module's Cs and beach limit, then the dune growth rate, Hs and groin
  strength refit with it on. Best fit (2026-10-29): Hs 3.0, growth 0.02-0.04, groin ~0.16.
- `2-test-2009-2025/`: the 1996-2009 calibration run unchanged on 2009-2025.
- `3-on-vs-off/`: the figure scripts for every test (module on against off, model steps, Cs, design, road, missed
  drops, storms, other forcings). Their figures are in the workspace's `dune-collision/4-hindcast/`, one folder per
  question since 2026-11-12, indexed in that folder's `README.md`.

The settings these settled on (2026-11-04): `../README.md`, Settings.

**2026-11-09 to 11:** the toe tests (Automorph, per cell, at the berm), the Cs sweeps on them, the road manager's
first-year bulldozing (found and fixed), and the new module options (row modes, lateral erosion, collapse tracking):
what, why, the numbers and the figures in `TESTING_2026-11-09_to_11.md`.
