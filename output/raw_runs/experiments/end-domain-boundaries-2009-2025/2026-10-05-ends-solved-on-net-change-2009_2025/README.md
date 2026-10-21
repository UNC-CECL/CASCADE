# 2026-10-05-ends-solved-on-net-change-2009_2025

**Question.** What GIS 1 and GIS 90 rates make the 2009–2025 edge-only run end on its own net-change target? The run had carried the calibration ends, which left GIS 1 at +9 m against +158 m observed and GIS 90 at −6 m against +36 m.

**Setup.** As in the 1996–2009 solve: `HAT_end_solve_net_change.py solve --period 2009`, full management, no groin, relocations off. The secant starts from the zeroBE matrix run. Target: `coastsat/net_change/2009_2025`, raw at GIS 1 and 7-domain LOWESS at GIS 90. Tolerance: 0.02 m/yr × 16 yr = 0.32 m. Run under the driver started 2026-10-05 (hence the folder date) and finished 2026-10-06.

**Answer.** GIS 1 **+32.7049**, GIS 90 **+21.0679** m/yr, from step 8 (the step limit). GIS 90 converged (+0.007 m). GIS 1 ended at +0.43 m against the target, but it does not respond smoothly: +33.19 gave +1.8 m and +27.61 gave −8.6 m, so ±1–2 m is its noise floor. Step tables: `solve.txt`.

**Result.** The matrix run `matrix/2009_2025/edgeBE/..._nogroin` was rerun on these values: GIS 1 +158.4 m (target +158.0), GIS 90 +36.0 (+36.0). Interior LOWESS-7 bias −10.4 m, RMSE 25.0 m, r 0.39 (with the calibration ends: −11.2 / 26.6 / 0.31).

**Status.** Current: edgeBE for 2009 in the site config. The test of the calibration field is still `domainBE`, which carries its own copy of set 1.
