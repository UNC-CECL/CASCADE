# end-domain-boundaries/2026-09-27-ends-solved-on-duneline-option-a

**Question.** Under option A, what end rates at GIS 1 and 90 make the model match the dune line?
- Option A: metres offset; Hs 2.0 / Tp 7.5 / asymmetry 0.6 / high-angle 0.5; full management, edgeBE, no groin.
- This is the 09-18 dune-line solve (`../2026-09-18-end-domains-solved-on-redigitized-duneline/`), redone because the model changed.
- It was needed to redraw `output/comparisons/target_comparison/` and `model_vs_observed/` fairly (Hannah, 2026-09-27: "Now that we have the new wave climate/offset, dont we need to redo th analyses in here").

**Method.**
- Driver: `be_dune_edgesolve_loop.py --windows 1996 2010 --smooth mean3 raw --tol 0.02`.
- Target: the dune-line endpoint change, read two ways at each end: the end domain's own value (`raw`), and the mean of it and its two inward neighbours (`mean3`, the main reading).
- Model side: the endpoint estimator.
- Start: the option A matrix zeroBE and edgeBE full-management runs.
- `brackets()` in the driver gained the `offsetmetres` token for this solve.

**Answer** (`solved.csv`, `RESULTS.md`; every probe in `loop_log.csv`):

| window | reading | step | GIS 1 / 90 (m/yr) | residuals (m/yr) |
|---|---|---|---|---|
| 1996–2010 | mean3 | 4 | −3.0 / +7.6 | +0.018 / +0.010, converged |
| 1996–2010 | raw | 2 | +0.3 / +10.6 | +0.030 / −0.007 |
| 2010–2024 | mean3 | 5 | +3.4 / +15.0 | +0.075 / +0.003 |
| 2010–2024 | raw | 6 | +3.1 / +14.8 | −0.056 / −0.006 |

**Four chains were accepted just outside the 0.02 tolerance.**
- The solver writes probes to 0.1 m/yr.
- At the 1996 raw GIS 1 it printed no further step: there was no slope, because the last two probes were identical.
- The 2010 GIS 1 does not respond smoothly within about 0.1 m/yr (seen before in the option A end solve), so those chains bounced around their best value.
- The loop was stopped after step 6. The three step-7 probes it had started were killed with it, and their folders were deleted.
- These misfits are within the one accepted for the option A ends (0.045 m/yr at 2010 GIS 1).

**Reading.** The dune-line ends are far smaller than the CoastSat-solved ones (1996 +4.8 / +17.5; 2010 +18.8 / +24.5), but they barely change the interior:

| window | RMSE vs CoastSat, CoastSat-solved | RMSE vs CoastSat, dune-solved | RMSE vs dune line, CoastSat-solved | RMSE vs dune line, dune-solved |
|---|---|---|---|---|
| 1996–2010 | 1.05 | 1.08 | 2.04 | 1.98 |
| 2010–2024 | 2.25 | 2.28 | 1.65 | 1.67 |

All in m/yr, over GIS 2–89. This is the same reading as the /10 solve. `RESULTS.md` also lists 1984–2004 and 2004–2024 rows: those are the /10-era CoastSat brackets the results script always includes, not part of this solve.

**Runs.** `mean3|raw/step<k>/<window>/edgeBE/<run>/` are on disk only (`.gitignore`; the paths exceed git's 260-character limit here). Logs are in `logs/`.
