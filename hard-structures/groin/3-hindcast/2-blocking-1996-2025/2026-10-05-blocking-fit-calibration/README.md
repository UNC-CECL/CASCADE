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

## The 2017→2018 step: where the Buxton fill went (checked 2026-10-05)

The model places the 2017 fill on GIS 6–16 (southernmost groin ~200 m into GIS 6, as reported) and keeps it there. Observations put the south end's sand **south of the groin**:

| seaward change, m | GIS 4 | GIS 5 | GIS 6 | GIS 7 | GIS 8 | GIS 9 |
|---|---|---|---|---|---|---|
| wet/dry photos, 2017 → 2018 | +8.6 | **+19.8** | **−10.6** | +23.2 | — | — |
| CoastSat, 2017 H1 → 2018 H2 | +17.8 | **+31.7** | **−3.2** | +15.3 | +22.7 | +30.5 |
| model, 1 Jan 2017 → 1 Jan 2018 (no groin) | +1.4 | +0.6 | **+31.5** | +33.3 | +34.9 | +35.7 |
| model, same, b 0.6 f 0.6 | +0.6 | −0.1 | **+32.3** | +33.5 | +34.9 | +35.7 |

In CoastSat, GIS 8–9 gain first (late 2017). GIS 6 is up briefly in 2018 H1 (+15 m), then back to about −3 m by 2018 H2, while GIS 4–5 gain +18 to +32 m and keep it through 2019. So within about a year, the sand placed at the groin end moved south past the groin. In the model, almost nothing reaches GIS 4–5, even with no groin (+2 to +4 m by 2020). BRIE's annual diffusion is too slow at the cape to move it.

So the gap jump is a **fill placement/bypass mismatch at GIS 5|6**, not a groin effect. The 2022 Buxton fill uses the same footprint and will repeat it.

## End domains with the groin on (checked 2026-10-05)

The edgeBE ends were solved with no groin. With the pinned blocking groin, calibration GIS 1 ends at +10.61 m against the +10.65 m target (residual −0.04 m, inside the 0.26 m tolerance), and GIS 90 is unchanged (+2.17 m). The test period is unchanged at both ends (GIS 1 +8.96 m, GIS 90 −5.92 m). No re-solve is needed. The groin's reach stops at GIS 2–4: in calibration it halves their accretion (+4.6 / +12.1 m at GIS 3 / 4, against +9.8 / +26.0 with no groin). That is closer to CoastSat at GIS 3–4 (about +5 / +16).
