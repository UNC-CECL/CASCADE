# 2026-10-10-test-end-calendar

**Question.** If the test end target stops depending on the DEM date, which calendar-year window should it be, and what does that do to the scores?

**Why.** Hannah (2026-10-10): the end targets should not be tied to the DEM date. Calibration's end already follows the calendar (calendar 2008 ±6 months, since 2026-10-09); the test end is still 17 Aug 2025 (the 2009 DEM date's anniversary) ±6 months. The starts stay DEM-centred (Hannah, same day), because they build the island offset the DEM topography starts from.

**Method.** `test_end_options.py`, no new runs. The 2009 start mean (17 Aug 2008 to 17 Aug 2010) is kept. Three ends are compared:

- `dem_2025`: the end in use, 2025-02-17 to 2026-02-17.
- `cal_2025`: calendar 2025 ±6 months, 2024-07-01 to 2026-06-30. This is calibration's rule. Its mean shoreline was built for this experiment (`data/hatteras_init/5-scr/1-observations/mean_shoreline/2024-07-01_2026-06-30/`). CoastSat ends 2026-01-13, so the window is off-centre.
- `year_2025`: calendar year 2025 only, 2025-01-01 to 2025-12-31.

The test runs are the current matrix (17 yr, 2009 through 2025, final state 1 Jan 2026). Scores use interior GIS 2–89, with both sides LOWESS-7.

## Results

**When the positions fall** (`window_dates.csv`). The mean date is measured from the model's final state (1 Jan 2026):

| end | positions | mean date | months before model end |
|---|---|---|---|
| in use (DEM date) | 32,133 | 2025-08-23 | −4.3 |
| option 1, calendar 2025 ±6 mo | 54,376 | 2025-04-23 | −8.3 |
| option 2, calendar year 2025 | 33,955 | 2025-07-12 | −5.7 |

For comparison, calibration's end (calendar 2008 ±6 months) has a mean date 7.2 months before its model end (1 Jan 2009).

**How the targets differ from the one in use** (interior, LOWESS-7):
- Option 1: mean −0.2 m, RMS 2.7 m, range −6.5 to +10.2 m. The largest differences are at GIS 6–8 (Buxton) and GIS 28–31.
- Option 2: mean +0.2 m, RMS 1.5 m.

**Scores** (`scores.csv`, bias / RMSE / r):

| run | in use | option 1 | option 2 |
|---|---|---|---|
| zeroBE, no groin | −13.4 / 28.3 / 0.27 | −13.2 / 27.6 / 0.29 | −13.6 / 28.4 / 0.27 |
| edgeBE, no groin | −11.6 / 25.1 / 0.42 | −11.4 / 24.4 / 0.44 | −11.9 / 25.2 / 0.43 |
| domainBE + groin | −20.2 / 34.0 / 0.36 | −20.0 / 32.9 / 0.41 | −20.4 / 33.9 / 0.38 |

The edgeBE end rates were solved on the end in use.

**Reading.** Either option moves the bias by less than 0.3 m and the RMSE by about 1 m. As with the start window, the fit can't choose. The choice is about the definition: option 1 uses the same rule as calibration; option 2 uses only the last modelled year.

**Status.** Diagnostic, pending Hannah's choice. Adopting either one means changing the test end in `NET_CHANGE_WINDOWS` and re-solving the test end rates.

Figure: `figures/test_end_options_edgeBE.png` (caption in `figures/supporting/CAPTIONS.md`).

Back to [the topic](../README.md).
