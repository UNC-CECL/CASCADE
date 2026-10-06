# 2026-10-06-test-target-and-the-2021-step

**Question.** How much of the 2009–2025 test misfit is the island-wide 2020→2021 CoastSat step (about +17 m, near-uniform, present in never-nourished transects; `5-scr/1-observations/detrended_position/README.md`)?

**Method.** `test_target_2021_step.py` scores the four existing test runs (no new runs) over the interior GIS 2–89, with both sides smoothed over 7 domains, against three targets:
- **full**: the plan's target, 2025-08-17 ±6 months minus the 2009 start mean, against the 1 Jan 2025 model state.
- **pre_step**: 2019-08-17 ±1 yr (2018-08-17 to 2020-08-17, new `mean_shoreline/2018-08-17_2020-08-17/`) minus the same start, against the 1 Jan 2019 state (10 model years, the same 7.5-month offset).
- **no_step**: the full target minus each domain's measured 2020→2021 step (`step_2021_by_domain_prefill.csv`, smoothed). Those step values come from the detrended annual medians, so they carry about a year of each transect's trend (small).

| run | full: bias / RMSE / r | pre_step | no_step |
|---|---|---|---|
| zeroBE, no groin | −11.9 / 27.8 / 0.25 | +2.8 / 15.4 / 0.34 | +5.2 / 25.7 / 0.21 |
| edgeBE, no groin | −11.2 / 26.6 / 0.31 | +3.4 / 15.2 / 0.38 | +6.0 / 24.8 / 0.28 |
| edgeBE + groin | −11.5 / 27.2 / 0.25 | +3.1 / 15.7 / 0.29 | +5.7 / 25.4 / 0.21 |
| **domainBE + groin** | **−18.2 / 32.1 / 0.36** | **−1.4 / 16.1 / 0.37** | **−1.1 / 26.8 / 0.33** |

Observed interior mean: full +11.7 m, pre_step −0.9 m, no_step −5.5 m. The step itself averages +17.1 m.

**Reading.**
1. **The test bias is the step.** Without it, BE set 1 gets the island mean right (bias −1.1 m no_step, −1.4 m pre_step), and better than any run without it (+3 to +6 m). The calibration budget carries into the test period; the 2021 step is what it misses.
2. **The alongshore pattern does not carry.** RMSE barely changes with set 1 (pre_step 16.1 vs 15.2 m for edgeBE; no_step 26.8 vs 24.8), and r moves only from about 0.3 to 0.37. Where 1996–2009 eroded and 2009–2025 did not, set 1 is wrong locally.
3. **For BE set 2:** the full residual contains the step, about +17 m. Spread over 16 years as a rate, that is about +1 m/yr of accretion everywhere, and every forward scenario would repeat a one-time event every year. Deriving set 2 from the no_step residual avoids that.

**Status.** Diagnostic; nothing is adopted. Outputs are `scores.csv` and `domain_residuals.csv` (the domainBE residuals for each target, plus the step).
