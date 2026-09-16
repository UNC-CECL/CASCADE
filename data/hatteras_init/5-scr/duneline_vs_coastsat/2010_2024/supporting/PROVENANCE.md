# Dune line vs CoastSat shoreline, 2010-2024

Written 2026-09-15 20:23 by `scripts/input_prep/5-scr/duneline_vs_coastsat/duneline_vs_coastsat.py`.

## Inputs

* dune lines: `2-brie-offset/raw_offsets/2009_duneline_offset_raw.csv`, `2023_duneline_offset_raw.csv` (first row per transect, domain mean, as `hindcast.load_absolute_dune_distance`). Both built by `duneline_to_raw_offsets.py`, so no GIS-vs-shapely metre between them.
* CoastSat LRR: `coastsat_lrr/2010_2024/transect_lrr_full.csv` (window 2010-01-01 to 2024-12-31, per-transect OLS).
* CoastSat endpoint: mean chainage within ±183 days of each survey date, per transect, from `coastsat_timeseries/`.
* transect → domain: `transect_domains/transect_domain_lookup.csv`.

## Survey dates

| line | date | source |
|---|---|---|
| 2010 | 2009-05-30 | Google Earth capture date; Hannah, 2026-09-15 — the 2009 line standing in for 2010 |
| 2024 | 2023-07-01 | **ASSUMED mid-year of 2023**; no flight date known for this line — the 2023 line standing in for 2024 |

Survey interval 14.09 yr. Sign: seaward positive in every column; a negative rate is retreat. Dune change is `-(ORIG_LEN_end - ORIG_LEN_start)`.

## Island-wide means (m/yr)

| dune line | CoastSat LRR | CoastSat endpoint |
|---|---|---|
| 0.09 | 1.15 | 0.50 |

## Agreement, per domain (y against dune rate x)

| y | n | r | slope | intercept | RMSE | bias (y − x) |
|---|---|---|---|---|---|---|
| CoastSat LRR | 90 | 0.51 | 0.76 | 1.08 | 2.00 | +1.06 |
| CoastSat endpoint | 90 | 0.51 | 0.70 | 0.44 | 1.64 | +0.41 |

CoastSat endpoint against CoastSat LRR (two estimators of the same series): r = 0.89, slope = 0.82, RMSE = 1.09, bias = -0.65 m/yr.

## Window occupancy

Median CoastSat observations per transect inside the start window: 9; inside the end window: 39. Transects with an empty window: start 0, end 0, of 906.

| window | asked for | median first obs | median last obs | truncated |
|---|---|---|---|---|
| start | 2008-11-28 to 2009-11-28 | 2008-12-04 | 2009-11-05 | no |
| end | 2022-12-30 to 2023-12-30 | 2023-01-21 | 2023-12-15 | no |

## Sensitivity of the endpoint rate to the assumed survey date

| centre shift | start | end | island mean endpoint (m/yr) | r vs dune | slope | RMSE |
|---|---|---|---|---|---|---|
| -6 mo | 2009-05-30 | 2022-12-30 | 0.80 | 0.38 | 0.51 | 1.92 |
| 0 | 2009-05-30 | 2023-07-01 | 0.50 | 0.51 | 0.70 | 1.64 |
| +6 mo | 2009-05-30 | 2023-12-30 | 0.33 | 0.58 | 0.78 | 1.50 |

## Read this before quoting it

* The dune line and the waterline are different features. A gap between them is beach-width change as much as it is disagreement.
* The LRR spans the calendar window; the endpoint spans the survey interval. They are not the same length of record.
* `n_obs_*` is the median per-transect count inside a one-year window. One storm inside a window moves that end.
