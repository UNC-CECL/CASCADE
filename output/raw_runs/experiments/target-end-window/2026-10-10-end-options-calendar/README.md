# 2026-10-10-end-options-calendar

**Question.** With no end tied to a DEM date, how does the selected end window (the last modelled year ±6 months) compare with the other calendar-based ends, in both periods?

**Why.** Hannah (2026-10-10): the end targets should not depend on the DEM date. The advisor prefers calibration's 1 Jul 2007 → 30 Jun 2009 window. This study applies the same rule to test (1 Jul 2024 → 30 Jun 2026) and keeps four alternatives for comparison. It replaces the 2026-10-07 option set, which mixed DEM-centred ends with a test model end that is now a year out of date.

**Method.** `end_options_calendar.py`, no new runs. Each end is defined from the last modelled year Y (calibration 2008, test 2025), using the same rule and colour in both periods:

| colour | key | window |
|---|---|---|
| red (drawn on top, same width) | `selected` | Y ±6 months: 1 Jul Y−1 → 30 Jun Y+1 |
| green | `year` | Y only |
| teal | `last_two` | Y−1 and Y |
| purple | `three_yr` | Y ±1 yr: Y−1 → Y+1 |
| orange | `end_state` | the model's final state (1 Jan Y+1) ±1 yr |

- The starts stay DEM-centred (the island offset's window): 1995-10-12 → 1997-10-12 and 2008-08-17 → 2010-08-17.
- CoastSat ends 2026-01-13, so the test windows that reach past 2025 hold only two weeks of 2026. The figures mark them with *.
- Six new mean shorelines were built for this study under `data/hatteras_init/5-scr/1-observations/mean_shoreline/`, all 906 transects at the default minimum of 10 positions. They are calendar 2008; 2007–2008; 2007–2009; 2024–2025; 2024–2026; and 2025-01-01 → 2027-01-01. The 2024-07-01 → 2026-06-30 window was built for `2026-10-10-test-end-calendar`.
- The runs are the current matrix: calibration 13 yr (final state 1 Jan 2009), test 17 yr (final state 1 Jan 2026).
- Scores use interior GIS 2–89, with both sides LOWESS-7.

## Results

**When the data falls** (`<period>/window_dates.csv`):

| end | calibration mean date (months before model end) | test mean date (months before model end) |
|---|---|---|
| selected | 2008-05-27 (−7.2) | 2025-04-23 (−8.3) |
| year | 2008-07-04 (−5.9) | 2025-07-12 (−5.7) |
| last_two | 2008-01-26 (−11.2) | 2025-01-27 (−11.1) |
| three_yr | 2008-05-15 (−7.6) | 2025-02-10 (−10.7) |
| end_state | 2008-10-28 (−2.1) | 2025-07-26 (−5.2) |

The selected window sits in a similar place in both periods: 64% of calibration's data falls in 2008 and 62% of test's in 2025. The 3-yr and end-state windows do not. Their halves after the last modelled year are full for calibration (2009) but almost empty for test (2026), so they mean different things in the two periods.

**How far each target sits from the selected one** (interior LOWESS-7, RMS): calibration 0.7–1.3 m, test 1.6–2.2 m. The largest local differences are at Buxton (GIS 6–8), GIS 28–31 and Wimble Shoals (GIS 64–68), up to about 7 m.

**Scores, edgeBE no groin** (bias / RMSE / r):

| end | calibration | test |
|---|---|---|
| selected | +2.4 / 17.0 / 0.53 | −11.4 / 24.4 / 0.44 |
| year | +3.2 / 17.3 / 0.53 | −11.9 / 25.2 / 0.43 |
| last_two | +1.7 / 16.5 / 0.54 | −9.9 / 23.1 / 0.45 |
| three_yr | +1.9 / 16.8 / 0.53 | −10.2 / 23.4 / 0.44 |
| end_state | +2.9 / 17.6 / 0.53 | −12.0 / 25.3 / 0.43 |

All three presets are in `<period>/scores.csv`. The edgeBE end rates were solved on the 2026-10-09 targets: the selected end for calibration, the DEM-anniversary end for test.

**Reading.** The options are within about 1.5 m in bias and 2 m in RMSE of each other in both periods, small against the misfit. The two windows weighted earliest (`last_two`, `three_yr` for test) score best because they sit before the 2025 accretion. That reflects when the data falls; it isn't evidence that those windows are better targets. The selected window is the one that falls the same way in both periods.

**Status.** Diagnostic. The calibration end in `NET_CHANGE_WINDOWS` is already the selected window. Adopting it for test means changing the test end there to 2024-07-01 → 2026-06-30 and re-solving the test end rates.

Figures: `<period>/figures/end_options_<period>_<preset>.png`, and both periods on one page in `figures/end_options_both_periods_<preset>.png` (captions in each `figures/supporting/CAPTIONS.md`). The two-period figure is on a ±55 m axis. It adds a second dashed line: the same run started from the ±6-month island offset (offset v3, from [`2026-10-09-6mo-starts`](../2026-10-09-6mo-starts/README.md)). Both runs are scored on the selected end with the ±1-yr start mean. edgeBE ±1 yr / ±6 months: calibration RMSE 17.0 / 16.8 m, bias +2.4 / +2.2 m; test RMSE 24.4 / 24.1 m, bias −11.4 / −11.1 m.

Back to [the topic](../README.md).
