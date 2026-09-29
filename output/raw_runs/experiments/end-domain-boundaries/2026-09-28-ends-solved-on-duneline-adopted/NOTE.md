# end-domain-boundaries/2026-09-28-ends-solved-on-duneline-adopted

**Question.** On the adopted model, what end rates at GIS 1 and 90 make the model match the dune line?

- Adopted model (2026-09-28): Barrier3D `hatteras/adopted` (the three overwash fixes and per-cell dune ceilings), storms `v3_trim24`. Otherwise option A: metres offset, Hs 2.0 / Tp 7.5 / asymmetry 0.6 / high-angle 0.5, full management, edgeBE, no groin.
- This is the 09-27 dune-line solve (`../2026-09-27-ends-solved-on-duneline-option-a/`), redone because the model changed (Hannah, 2026-09-28: "re-solve the dune line ends on the adopted setup"). Without it, `model_vs_observed/` and `target_comparison/` would put the adopted-model matrix beside pre-adoption dune-line runs.

**Method.** The 09-27 method, unchanged:

- `be_dune_edgesolve_loop.py --exp end-domain-boundaries/2026-09-28-ends-solved-on-duneline-adopted --windows 1996 2010 --smooth mean3 raw --tol 0.02`.
- Target: the dune-line endpoint change. The `raw` reading is the end domain's own value. The `mean3` reading, the main one, is the mean of the end domain and its two inward neighbours.
- Model side: the endpoint estimator.
- Start: the rebuilt matrix's zeroBE and edgeBE full-management runs, on the adopted model.
- `be_dune_edgesolve_results.py --solved 1996:mean3:3 1996:raw:2 2010:mean3:3 2010:raw:4` wrote `solved.csv`, `skill.csv` and `RESULTS.md`.

**Answer** (GIS 1 / GIS 90, m/yr):

| window | reading | step | ends | residuals |
|---|---|---|---|---|
| 1996–2010 | mean3 | 3 | −2.8 / +7.9 | −0.002 / −0.016, converged |
| 1996–2010 | raw | 2 | +0.6 / +10.6 | +0.010 / −0.015, converged |
| 2010–2024 | mean3 | 3 | +2.1 / +14.2 | **−0.032** / −0.014, accepted |
| 2010–2024 | raw | 4 | +1.9 / +14.0 | +0.010 / −0.000, converged |

For comparison, 09-27 found −3.0 / +7.6, +0.3 / +10.6, +3.4 / +15.0 and +3.1 / +14.8.

**2010 mean3 was accepted just outside the 0.02 tolerance, as on 09-27.**

- The solver writes probes to 0.1 m/yr, and GIS 1 sat 0.032 off at the closest reachable value.
- Steps 4–6 repeated that probe without a slope, and the loop stopped after 6.

**Reading.** As before, the dune-line ends barely change the interior. Interior RMSE over GIS 2–89 is within 0.06 m/yr of the CoastSat-solved run against both observations (`RESULTS.md`).

**Runs.** `mean3|raw/step<k>/<window>/edgeBE/<run>/` stay on disk only (`.gitignore`): the paths pass git's 260-character limit. Logs are in `logs/`, and every probe is in `loop_log.csv`.
