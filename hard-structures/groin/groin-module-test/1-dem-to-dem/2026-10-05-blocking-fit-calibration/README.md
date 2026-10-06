# 2026-10-05-blocking-fit-calibration

**Question.** Which blocking groin, with strength b and post-failure fraction f, reproduces the observed D5−D6 gap over the DEM-to-DEM calibration period (1996 → 2009)? And how does that groin do on the 2009 → 2025 test period, which it was not fitted to?

**Setup.** `blocking_fit.py` drives the unchanged runner:
- full management, the solved edgeBE ends (+1.4981 / +10.5659), relocations off;
- groin failure instant from the 2004 step;
- runs filed under `output/raw_runs/experiments/groin/2026-10-05-blocking-fit-dem-to-dem/`.

The score is the gap change at the wet/dry photo dates inside the window, compared with the observed change from the interpolated start value (the 2026-09-29 date RMSE). Calibration has only **three** dates: 1997, 2004 and 2008. The grid covers b 0.3–0.9 × f 0.1–0.8 (49 runs, all clean). Scores are in `grid_scores.csv`; logs are in `logs/`.

## Calibration, 1996–2009 (date RMSE, m)

No groin: 69.2 m (model −20 / −74 / −86 against observed +2 / +16 / −10).

| b \ f | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.8 |
|---|---|---|---|---|---|---|---|
| 0.3 | 46.1 | 44.9 | 44.0 | 43.0 | 42.0 | 41.0 | 39.1 |
| 0.4 | 37.5 | 35.9 | 34.1 | 32.6 | 31.2 | 29.6 | 26.9 |
| 0.5 | 28.6 | 26.2 | 23.9 | 21.6 | 19.4 | 17.4 | 13.6 |
| 0.6 | 21.5 | 18.1 | 14.8 | 11.1 | 7.5 | **4.0** | 5.1 |
| 0.7 | 19.3 | 15.9 | 13.2 | 11.3 | 11.1 | 12.6 | 19.9 |
| 0.8 | 26.7 | 25.2 | 24.6 | 25.1 | 26.8 | 29.8 | 38.5 |
| 0.9 | 39.1 | 38.7 | 39.2 | 40.8 | 43.5 | 47.4 | 58.0 |

**Best: b 0.6, f 0.6** (model +4 / +13 / −16). This is an interior minimum. b is well pinned down by the 2004 date. f rests on the 2008 date alone, so f 0.5–0.8 all sit within a few metres of it.

## Test, 2009–2025 (not fitted)

The observed gap changes at 2014–2023 are −5 / −19 / −11 / −42 / −61 / −40 / −29 / −49 m.

| run | all dates | 2014–2017 | 2018–2023 | 2017→2018 step, model vs observed |
|---|---|---|---|---|
| no groin | 31.1 | 45.6 | 17.6 | +31 vs −31 |
| b 0.6 f 0.5 | 31.0 | 16.0 | 37.1 | +32 vs −31 |
| **b 0.6 f 0.6** | 36.3 | **9.5** | 45.4 | +32 vs −31 |
| b 0.6 f 0.8 | 51.0 | 8.5 | 64.2 | +33 vs −31 |

**Reading.** Before 2018 the calibrated groin tracks the test period well: 9.5 m against 45.6 m with no groin. Every run, with or without a groin, then jumps about +32 m between 2017 and 2018, while the observed gap falls 31 m. The step is the same size with and without the groin, so it is not a groin effect. It lands in the year the Buxton 2017 fill (GIS 6–16) fires; the fill raises GIS 6 but not GIS 5. After that step the groin's higher gap only adds error. The leading candidate is how the model places the 2017 fill at the groin's downdrift cell. That is not tested yet.

**Status.** b 0.6 / f 0.6 is the calibration answer. Pinning it in the config is pending a decision.
