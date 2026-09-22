# 2026-09-19-edgesolve-lrr1996_2024

**Question.** Hannah's CoastSat target is the FULL-PERIOD 1996-2024 LRR in both
model windows, not each window's own LRR (which the runner grades against and
the matrix end values were solved on). What do the two locked end domains
(GIS 1 and 90) carry when solved against it?

**Method.** `be_dune_edgesolve_loop.py --target coastsat --coastsat-window
1996_2024 --windows 1996 2010` (the `--coastsat-window` option was added to
the loop and to `be_edge_domain_solve.py` for this): the same Newton solve
as the matrix (model lrr_m_yr against the target table, GIS 1 raw, GIS 90
LOESS-10), with the 1996-2024 table in place of the window's own. Brackets are
the 09-18 matrix zeroBE and edgeBE full-management runs. Full management, no
groin, relocations off, Hs 2.5, offsets v1.

**Answer.** Converged in three steps each (`loop_log.csv`):

| window | GIS 1 / 90 (m/yr) | vs the matrix (sub-period) values |
|---|---|---|
| 1996-2010 | +28.5 / +24.5 | +32.2 / +10.0 |
| 2010-2024 | +37.1 / +25.4 | +72.6 / +31.3 |

The config (`HATTERAS_BE_EDGE_ONLY`) and the matrix are NOT changed; these
runs are what `output/comparisons/target_comparison/projected/`
(named `coastsat_full_period_lrr/` until 2026-09-21)
pairs with the full-period CoastSat target.

**Layout.** `coastsat/step<k>/<window>/edgeBE/<run>/`, `logs/`, `loop_log.csv`.
