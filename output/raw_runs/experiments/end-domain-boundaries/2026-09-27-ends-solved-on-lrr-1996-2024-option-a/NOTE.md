# end-domain-boundaries/2026-09-27-ends-solved-on-lrr-1996-2024-option-a

**Question.** Under option A, what end rates at GIS 1 and 90 make the model match the full-period 1996–2024 CoastSat LRR in each window?
- Option A: metres offset; Hs 2.0 / Tp 7.5 / asymmetry 0.6 / high-angle 0.5; full management, edgeBE, no groin.
- This is the 09-19 solve (`../2026-09-19-end-domains-solved-on-lrr-1996-2024/`), redone because the model changed.
- It supplies the `ends_solved_on_coastsat` set of `output/comparisons/target_comparison/projected/` (Hannah, 2026-09-27).

**Method.**
- Driver: `be_dune_edgesolve_loop.py --windows 1996 2010 --target coastsat --coastsat-window 1996_2024 --tol 0.02`.
- The same solve as the matrix ends: the model's OLS rate against the target table, GIS 1 against the raw domain mean, GIS 90 against the LOWESS-10 value.
- The target is the 1996–2024 table rather than each window's own.
- Start: the option A matrix zeroBE and edgeBE runs.

**Answer** (`loop_log.csv`; `target_comparison.py` reads the last step):

| window | step | GIS 1 / 90 (m/yr) | residuals (m/yr) |
|---|---|---|---|
| 1996–2010 | 7 (same probe as step 4) | +4.5 / +27.5 | +0.022 / −0.000 |
| 2010–2024 | 7 | +4.8 / +20.5 | +0.046 / −0.005 |

The /10-era values were +28.5 / +24.5 and +37.1 / +25.4.

**How 2010 GIS 1 got there.**
- Its response was flat at first: cutting the imposed rate from 18.8 to 10 only moved the residual from +4.0 to +3.3.
- So the secant probed −32.9, which overshot to −19.5.
- That bracketed the root, and the chain closed in from there: 3.7 gave −1.38, 6.5 gave +1.97, 4.9 gave +0.17, 4.7 gave −0.097, and 4.8 gave +0.046.
- The loop was stopped after step 7. The step-8 probes (repeats) were killed and their folders deleted.
- 1996 GIS 1 stayed at +0.022 from step 4 on, because the solver printed no further GIS 1 step.

**Runs.** `coastsat/step<k>/<window>/edgeBE/<run>/` are on disk only (`.gitignore`; the paths exceed git's 260-character limit here). Logs are in `logs/`.
