# 2026-10-07-end-window-options

**Question.** How much does each period's target change if its end line is the CoastSat mean over a different window? This covers the 1996–2009 calibration target and the 2009–2025 test target.

**Why.** Each run ends on 1 Jan of its end year: calibration on 1 Jan 2009 and test on 1 Jan 2025. The forcing of the two runs never overlaps. Each current end line, however, is centred on a DEM date (or its anniversary), not on the model's end. For calibration, that means the line includes Nor'Ida and the other storms of late 2009 and 2010, which the calibration model never sees.

**Method.** `target_end_window.py`. No new runs. Within a period, every target shares that period's start mean, so only the end changes. Each period's three existing matrix runs (zeroBE no groin, edgeBE no groin, domainBE + groin) are scored against each target over the interior GIS 2–89, with both sides smoothed over 7 domains.

```
2026-10-07-end-window-options/
├── README.md
├── target_end_window.py      writes both period folders
├── calibration/              1996–2009: 4 end lines
│   ├── domain_net_change.csv
│   ├── scores.csv
│   ├── window_dates.csv
│   └── figures/end_lines_vs_model_calibration_<zeroBE|edgeBE|domainBE>.png
└── test/                     2009–2025: 5 end lines
    └── (same files, figures end_lines_vs_model_test_<preset>.png)
```

In each `domain_net_change.csv`, the columns are: `<line>_m` and `<line>_lowess7_m` (observed net change), `diff_<line>_*` (minus the current line), and `model_<preset>_lowess7_m`. `scores.csv` and `window_dates.csv` keep a `period` column. Figure captions are in `<period>/figures/supporting/CAPTIONS.md`. Each figure overlaps all of the period's end lines (LOWESS-7) with one run's change and carries the house annotation layer.

> **2026-10-09:** the calendar 2008 ±6-month line (`mid_2009`) is now the calibration target, and the test runs through 2025 (17 years, ending 1 Jan 2026). The test numbers below are from the 16-year runs, now in `archive/2026-10-09-pre-test-2025/`; rerunning `target_end_window.py` would score the 17-year runs against a 1 Jan 2025 model end.

## Summary

| period | end lines vs the current one (interior, LOWESS-7) | effect on the scores | status |
|---|---|---|---|
| calibration 1996–2009 | mid_2009 +3.6 m seaward, jan_2009 +3.1 m, year_2009 +3.5 m; shape kept (r ≥ 0.988) | zeroBE/edgeBE bias 5.1/6.0 → 1.5/2.4 m (mid_2009), 1.6/2.5 m (year_2009) | diagnostic; not adopted |
| test 2009–2025 | −2.1 to +0.7 m on average; year_2025 is the same as the current line (+0.1 m, r 0.998) | bias moves 1–3 m, against a 10–18 m shortfall | diagnostic; not adopted |

## Calibration 1996–2009 (`calibration/`)

All targets start from the 1995-10-12 to 1997-10-12 mean.

| line | window | what it is centred on |
|---|---|---|
| dem_centred (current) | 2008-08-17 to 2010-08-17 | 2009 DEM date ±1 yr |
| mid_2009 | 2007-07-01 to 2009-06-30 | calendar 2008 ±6 months; new `mean_shoreline/2007-07-01_2009-06-30/` (906 transects, median 34 positions) |
| jan_2009 | 2008-01-01 to 2010-01-01 | model end (1 Jan 2009) ±1 yr; new `mean_shoreline/2008-01-01_2010-01-01/` (906 transects, median 32 positions) |
| year_2009 (added 10-09) | 2009-01-01 to 2009-12-31 | every position in the end year; new `mean_shoreline/2009-01-01_2009-12-31/` built with **min 5 positions, not 10** (878 of 906 transects, median 9, SE median 4.3 m) |

2009 is a thin imagery year. At the usual minimum of 10, only 438 transects qualify and nothing south of GIS 23 does; Buxton (GIS 1–22) has a median of 6 positions. The threshold was lowered for this line only, to cover all 90 domains, so its transect means are noisier than the other lines'.

**When the positions fall** (`window_dates.csv`). Imagery is uneven, so a window's mean date is not its centre.

| line | mean position date | vs model end (1 Jan 2009) | years from the start mean (1996-08-16) |
|---|---|---|---|
| dem_centred | 2009-07-29 | +6.9 months | 12.95 |
| mid_2009 | 2008-05-27 (64% of positions from 2008) | −7.2 months | 11.78 |
| jan_2009 | 2008-10-28 (69% of positions from 2008) | −2.1 months | 12.20 |
| year_2009 | 2009-07-18 (8,421 positions, a third of the others) | +6.5 months | 12.92 |

The model runs 13.0 yr. The current target is shifted about 7 months later than the model, but its span matches the model's length. The mid-2009 target is shifted about 7 months earlier and is 1.2 yr shorter. The jan_2009 window is centred on the model's end, but 2008 has more than twice as many positions as 2009, so its mean date still falls about 2 months early (0.8 yr shorter than the model).

**mid_2009 vs current.**
- Per transect: on average **+3.6 m seaward** (median +3.4, SD 6.0 m); 137 of 905 transects differ by more than 10 m.
- Per domain, interior: mean +3.7 m, RMS 6.5 m, range −10.3 to +19.4 m (largest at GIS 24–25 in Avon, and GIS 80–83 near the Rodanthe pier). Seaward in 70% of domains raw (88% at LOWESS-7).
- At LOWESS-7: mean +3.6 m, RMS 5.0 m, range −4.5 to +12.7 m, r = 0.988 with the current target.
- Ends: GIS 1 −1.2 m, GIS 90 +1.9 m.
- The interior observed mean is −9.8 m (current) and −6.1 m (mid-2009). At the island trend (about −0.75 m/yr), the 1.2 yr shorter span accounts for about 0.9 m of the 3.6 m. The rest is what late 2008 to mid-2010 did to the shoreline, beyond the trend.

**jan_2009 vs current.**
- Per transect: +2.9 m seaward on average (median +3.0, SD 4.8 m); 56 of 905 transects differ by more than 10 m.
- Per domain, interior: mean +3.0 m, RMS 5.2 m, range −8.0 to +15.8 m. Seaward in 72% of domains raw (83% at LOWESS-7).
- At LOWESS-7: mean +3.1 m, RMS 4.2 m, range −3.1 to +10.8 m, r = 0.991 with the current target.
- Ends: GIS 1 +0.1 m, GIS 90 +1.5 m.
- Interior observed mean −6.7 m. The jan_2009 line sits between the other two, nearer mid_2009 (about 0.5 m landward of it).

**year_2009 vs current.**
- Per transect: +3.5 m seaward on average (median +3.8, SD 4.3 m); 57 of 878 transects differ by more than 10 m.
- Per domain, interior: mean +3.4 m, RMS 4.9 m, range −4.7 to +12.6 m. Seaward in 81% of domains raw (89% at LOWESS-7).
- At LOWESS-7: mean +3.5 m, RMS 4.2 m, range −2.2 to +8.1 m, r = 0.994 with the current target.
- Ends: GIS 1 +1.2 m, GIS 90 +2.5 m.
- Interior observed mean −6.3 m. Its mean date (2009-07-18) is almost the same as the current line's (2009-07-29), so the 3.5 m is not timing or span. The current line's extra content is late 2009 to mid-2010: Nor'Ida and the 2010 storms. Their net effect is about 3.5 m of erosion, which the model never sees.

**Scores** (`scores.csv`, interior LOWESS-7, model minus observed, bias / RMSE / r)

| run | dem_centred (current) | mid_2009 | jan_2009 | year_2009 |
|---|---|---|---|---|
| zeroBE, no groin | +5.1 / 18.6 / 0.53 | +1.5 / 16.6 / 0.56 | +2.1 / 17.1 / 0.55 | +1.6 / 18.3 / 0.52 |
| edgeBE, no groin | +6.0 / 19.0 / 0.51 | +2.4 / 17.0 / 0.53 | +2.9 / 17.6 / 0.53 | +2.5 / 18.7 / 0.49 |
| domainBE + groin | +0.3 / 4.4 / 0.98 | −3.3 / 5.6 / 0.98 | −2.7 / 5.3 / 0.98 | −3.2 / 6.0 / 0.98 |

domainBE (BE set 1) was solved on the current target, so its lower score on that target is built in, not evidence. The runs that were not fitted to either target (zeroBE, edgeBE) sit about 3.6 m closer to the mid-2009 line, and their RMSE drops about 2 m.

**Reading.** Every alternative sits 3.1–3.6 m seaward of the current line. year_2009 shows this is not about timing: it has the same mean date as the current line, and the gap is what the storms of late 2009 and 2010 removed. The bias of the unfitted runs drops by the same ~3.5 m, but their RMSE barely moves on year_2009 (noisier line). The choice of end line shifts the target about 3.6 m seaward overall, by up to 10–19 m in a few domains, and leaves the alongshore pattern unchanged. If adopted, the change would show up mainly as a shift of about +0.28 m/yr in the per-domain BE field (3.6 m / 13 yr), concentrated in Avon and northern Rodanthe. The ends barely move. The test start (the 2009 offset) is unchanged; if the mid-2009 line were adopted, the calibration end and the test start would no longer be the same line.

## Test 2009–2025 (`test/`)

All targets start from the 2008-08-17 to 2010-08-17 mean, the line the test offset is built from. The model steps 1 Jan 2009 to 1 Jan 2025 (16 yr). New mean shorelines: `mean_shoreline/2023-07-01_2025-06-30/`, `2024-01-01_2026-01-01/`, `2024-07-01_2025-07-01/`, and (10-09) `2025-01-01_2025-12-31/` (906 transects each, min 10 positions; calendar 2025 has a median of 37).

| line | window | mean position date | vs model end (1 Jan 2025) | years from the start mean (2009-07-29) |
|---|---|---|---|---|
| dem_2025 (current) | 2025-02-17 to 2026-02-17, 17 Aug 2025 ±6 months | 2025-08-23 | +7.7 months | 16.07 |
| cal_2024 | 2023-07-01 to 2025-06-30, calendar 2024 ±6 months | 2024-06-25 | −6.2 months | 14.91 |
| jan_2025 | 2024-01-01 to 2026-01-01, model end ±1 yr | 2025-01-27 | +0.9 months | 15.49 |
| jan_2025_6mo | 2024-07-01 to 2025-07-01, model end ±6 months | 2024-12-24 | −0.3 months | 15.41 |
| year_2025 (added 10-09) | 2025-01-01 to 2025-12-31, every position in the end year | 2025-07-12 | +6.3 months | 15.95 |

17 Aug 2025 is the 2009 DEM date's anniversary, not a survey; there is no 2025 DEM.

**How far the alternatives sit from the current end** (interior GIS 2–89, raw per domain): cal_2024 −2.0 m (RMS 7.4), jan_2025 −2.1 m (RMS 4.4), jan_2025_6mo +0.7 m (RMS 5.8), year_2025 +0.1 m (RMS 1.2; per transect +0.02 m, none over 10 m). All keep the alongshore shape (r 0.96–0.99 at LOWESS-7). The south end moves most: GIS 1 differs by −17 to −41 m.

**Scores** (`scores.csv`, interior LOWESS-7, bias / RMSE / r)

| run | dem_2025 (current) | cal_2024 | jan_2025 | jan_2025_6mo | year_2025 |
|---|---|---|---|---|---|
| zeroBE, no groin | −11.9 / 27.8 / 0.25 | −10.6 / 24.6 / 0.32 | −10.2 / 25.8 / 0.28 | −13.2 / 27.0 / 0.31 | −12.1 / 27.9 / 0.26 |
| edgeBE, no groin | −10.4 / 25.1 / 0.39 | −9.0 / 21.9 / 0.45 | −8.7 / 23.1 / 0.42 | −11.6 / 24.3 / 0.44 | −10.6 / 25.2 / 0.40 |
| domainBE + groin | −18.2 / 32.1 / 0.36 | −16.9 / 28.4 / 0.47 | −16.5 / 29.6 / 0.42 | −19.5 / 31.1 / 0.45 | −18.5 / 32.0 / 0.38 |

**Reading.** In the test period the end window moves the target by 1–3 m on average, much less than the 10–18 m the model falls short by. Every window shows the island building out about 10–13 m on average while the runs sit near zero. Centring the end on the model's end date does not close that gap. The whole-year 2025 line is the current target in all but name: the two windows share 10.5 months.

## History

- 2026-10-07: calibration mid_2009 line, then jan_2009, then the test period, all in one folder `calibration-end-window/2026-10-07-mid-2009-end-line/`, with `calibration_end_window.py`.
- 2026-10-09: renamed to `target-end-window/2026-10-07-end-window-options/` because it covers both targets. The script became `target_end_window.py` and now writes one folder per period. The calibration figures gained a `calibration_` infix (`end_lines_vs_model_<period>_<preset>`), and the test model columns dropped theirs (`model_test_<preset>_*` → `model_<preset>_*`, now in their own file). The regenerated CSVs match the old ones exactly.
- 2026-10-09: added the whole-end-year lines year_2009 (min 5 positions) and year_2025.

Back to [the topic](../README.md).
