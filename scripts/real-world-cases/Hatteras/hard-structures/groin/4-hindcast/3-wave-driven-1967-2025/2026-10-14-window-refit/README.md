# The groin's strength for each hindcast window (2026-10-14)

**Question.** The runner was set to the groin's whole-life strength, 0.25, fitted on one continuous run from 1967 (`../2026-10-13-pregroin-planform/`). The hindcast matrix does not run that way: each window starts from its observed shoreline (1996 or 2009), which already holds the fillet the groins built. What strength does a run from each observed start need, on the runner as it is now?

**Setup.** The unchanged runner, the corrected groin built in (commit `563a542e`), the dune settings of 2026-10-13 (design 3.94 m, island-wide Dmaxel). edgeBE, full management, no relocation, groin off from 2017: the setup of the 10-12 window fit. The strength is passed as `HAT_GROIN_STRENGTH`. Scored on the gap across the groins (GIS 5 minus GIS 6) against annual CoastSat and, for 1996-2009, the three photo dates, with the scoring functions of `../2026-10-09-corrected-groin/corrected_groin.py`.

    python window_refit.py run   [--streams 4]
    python window_refit.py score              # -> window_scores.csv

Runs: `output/raw_runs/experiments/groin-wave-driven/2026-10-14-window-refit/<strength>/<window>/`. Logs: `logs/`.

## Results

**1996-2009** (13 annual years; photos 1997, 2004, 2008)

| strength | annual RMSE (m) | bias (m) | photo RMSE (m) |
|---|---|---|---|
| no groin | 47.9 | -44.5 | 70.7 |
| 0.10 | 19.5 | -13.6 | 42.1 |
| 0.13 | 15.0 | -5.7 | 35.1 |
| **0.16** | **13.4** | +3.2 | 27.9 |
| 0.19 | 16.9 | +11.2 | 23.0 |
| 0.22 | 23.5 | +19.8 | 16.7 |
| 0.25 | 30.4 | +27.8 | 16.0 |

**2009-2025** (16 annual years)

| strength | annual RMSE (m) | bias (m) | 2010-2016 (m) | 2017-2025 (m) |
|---|---|---|---|---|
| no groin | 21.9 | -13.1 | 32.6 | 4.9 |
| 0.04 | 15.8 | -6.4 | 22.8 | 6.0 |
| 0.06 | 12.6 | -3.9 | 17.6 | 6.4 |
| 0.08 | 10.5 | -1.1 | 13.4 | 7.4 |
| **0.10** | **9.7** | +1.3 | 11.2 | 8.3 |
| 0.12 | 9.9 | +3.7 | 11.6 | 8.4 |
| 0.16 | 14.1 | +9.2 | 17.9 | 10.2 |
| 0.25 | 28.1 | +20.5 | 39.3 | 14.2 |

## Findings

1. **0.16 for 1996-2009 and 0.10 for 2009-2025**, the same as the 10-12 fit (`../2026-10-09-corrected-groin/`, findings 18-22), whose scores these match to 0.1 m. The 10-13 dune changes do not reach the groin.
2. **0.25 overshoots from either observed start**: by +28 m on average in 1996-2009, and in 2009-2025 it is worse than no groin (28.1 against 21.9 m). The start already holds the fillet, and the model's fillet is narrower than the real one (two domains against three or four), so a whole-life strength piles sand beside the face.
3. **The photo dates point higher** (0.22-0.25 on 1997, 2004, 2008), as on 10-12. Three dates against thirteen annual years: annual CoastSat decides, as in the 10-12 fit.
4. **The drop from 0.16 to 0.10 is consistent with the groin deteriorating**, but the window strengths measure how fast the model relaxes the observed fillet, not the groin itself: run as one history from 1967 they leave the gap 30-40 m below the record (`../2026-10-13-pregroin-planform/`).
5. **The 2009 value is fitted on the test window**, so for the groin 2009-2025 is no longer an independent test.

## What changed in the runner

`hat_run.yaml` `groin.strength_1996` 0.16 and `groin.strength_2009` 0.10, with an optional `groin.strength` that overrides both (and is required for any other start year with the groin on; 0.25 for a run that starts before the fillet exists, such as a forecast of the rebuilt groin). `HAT_run_matrix.py` keys each period's groin runs on that period's strength. At its defaults the runner reproduces the 0.16 and 0.10 runs here to 0.0 m (`output/raw_runs/experiments/groin-wave-driven/2026-10-14-window-refit-check/`).

## Figures

- `figures/matrix_groin_effect_full_management_edgeBE.png` (`matrix_groin_effect_figure.py`): net change along the island, groin on vs off, both windows.
- `figures/groin_simple_<calibration|test>.gif` (`groin_simple_gif.py`; git-ignored): shoreline change since the start around the groins, year by year, without and with the groin, over the fixed start and the observed change (CoastSat net-change target). Calibration 1996-2009 at 0.16 (GIS 1-15); test 2009-2025 at 0.10, off from 2017 (GIS 2-15). Matrix runs, full management, edgeBE.
