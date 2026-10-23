# Split-sample validation of the groin's strength (2026-10-15)

**Question.** The window strengths, 0.16 for 1996-2009 and 0.10 for 2009-2025, were each fitted on their own
window (`../2026-10-14-window-refit/`), so neither window is an independent test of the groin. If the strength
fitted on one window is run unchanged on the other, does the groin still improve on no groin, and by how much?

**Setup.** The unchanged runner, with the module's tested-range guard (it stayed silent in every run), and the
setup of the window refit: edgeBE, full_management, no relocation, groin off from 2017. Each window is run at
no groin, its own strength, the other window's strength, and the whole-life 0.25 (fitted on 1967-2025, which
contains both windows, so not a test). Scored on the annual gap across the groins (GIS 5 minus GIS 6) against
CoastSat, with the window refit's scoring. Skill is against the no-groin run: 1 − MSE / MSE(no groin). The share
of the fit's improvement is how much of the gain from no groin to the window's own fit the transferred
strength keeps.

    python split_sample.py run   [--streams 4]
    python split_sample.py score                 # -> split_sample_scores.csv

Runs: `output/raw_runs/experiments/groin-wave-driven/2026-10-15-split-sample/<strength>/<window>/`. Logs: `logs/`. The
scores match the window refit's runs of the same cells to 0.0 m.

## Results

| window | strength | role | annual RMSE (m) | bias (m) | skill vs no groin | share of the fit's improvement |
|---|---|---|---|---|---|---|
| 1996-2009 | none | no groin | 47.9 | −44.5 | 0 | 0 |
| 1996-2009 | 0.16 | fitted on this window | 13.4 | +3.2 | 0.92 | 100% |
| 1996-2009 | **0.10** | **fitted on 2009-2025: the independent test** | **19.5** | −13.6 | **0.84** | **83%** |
| 1996-2009 | 0.25 | whole life (not independent) | 30.4 | +27.8 | 0.60 | 51% |
| 2009-2025 | none | no groin | 21.9 | −13.1 | 0 | 0 |
| 2009-2025 | 0.10 | fitted on this window | 9.7 | +1.3 | 0.81 | 100% |
| 2009-2025 | **0.16** | **fitted on 1996-2009: the independent test** | **14.1** | +9.2 | **0.59** | **64%** |
| 2009-2025 | 0.25 | whole life (not independent) | 28.1 | +20.5 | −0.65 | −51% |

## Findings

1. **The groin's strength transfers between windows with skill.** Calibrated on 1996-2009 and run unchanged on
   2009-2025, the groin cuts the error from 21.9 to 14.1 m (skill 0.59), keeping 64% of the improvement a fit
   on that window would give. The reverse transfer keeps 83% (skill 0.84). In both directions the groin with
   a borrowed strength is better than no groin.
2. **The transferred strength errs the way the groin aged.** The 1996 strength is too strong for 2009-2025
   (bias +9.2 m) and the 2009 strength too weak for 1996-2009 (−13.6 m), as a groin that weakened between the
   windows would make it (1995 last repair, 2003 damage, decay from 2015; `../../../1-observations/data/structure_history.md`).
   A strength falling from 0.16 to 0.10 is consistent with that record, though the window strengths also carry
   how fast the model relaxes the observed fillet (`../2026-10-14-window-refit/`, finding 4).
3. **The whole-life strength is not a substitute.** 0.25, fitted on one continuous run from 1967, overshoots
   from either observed start and is worse than no groin in 2009-2025 (skill −0.65), as the window refit found.
   Runs that start from an observed shoreline need a window strength.
4. **What this does and does not test.** Only the groin's strength is held out. The other tuned settings (the
   edgeBE ends, the dune settings, option A waves) were set with both windows in view, so this is a split-sample
   test of the groin, not of the whole model.

**For reporting.** Calibrated on 1996-2009 (strength 0.16) and validated on 2009-2025 without refitting, the
groin module reduces the RMSE of the gap across the Buxton groins from 21.9 m (no groin) to 14.1 m (skill 0.59
against no groin). Refitting on the validation window (0.10) lowers it to 9.7 m; the difference is consistent
with the groin's documented deterioration.
