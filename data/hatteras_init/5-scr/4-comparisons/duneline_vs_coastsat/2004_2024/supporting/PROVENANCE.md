# Dune line vs CoastSat shoreline, 2004-2024

Written 2026-09-18 15:03 by `scripts/input_prep/5-scr/duneline_vs_coastsat/duneline_vs_coastsat.py`.

## Inputs

* dune lines: `2-brie-offset/raw_offsets/2004_duneline_offset_raw.csv`, `2023_duneline_offset_raw.csv` (first row per transect, domain mean, as `hindcast.load_absolute_dune_distance`). Both built by `duneline_to_raw_offsets.py`, so no GIS-vs-shapely metre between them.
* CoastSat LRR: `3-rates/coastsat/lrr/2004_2024/transect_lrr_full.csv` (window 2004-01-01 to 2024-12-31, per-transect OLS).
* CoastSat endpoint: mean chainage within ±183 days of each survey date, per transect, from `coastsat_timeseries/`.
* transect → domain: `transect_domains/transect_domain_lookup.csv`.

## Survey dates

| line | date | source |
|---|---|---|
| 2004 | 2004-05-25 | Google Earth capture date (the raw_GE frames carry none); Hannah, 2026-09-15 |
| 2024 | 2023-07-01 | **ASSUMED mid-year of 2023**; no flight date known for this line — the 2023 line standing in for 2024 |

Survey interval 19.10 yr. Sign: seaward positive in every column; a negative rate is retreat. Dune change is `-(ORIG_LEN_end - ORIG_LEN_start)`.

## Island-wide means (m/yr)

| dune line | CoastSat LRR | CoastSat endpoint |
|---|---|---|
| -0.09 | 0.49 | 0.41 |

## Agreement, per domain (y against dune rate x)

| y | n | r | slope | intercept | RMSE | bias (y − x) |
|---|---|---|---|---|---|---|
| CoastSat LRR | 90 | 0.25 | 0.23 | 0.51 | 2.06 | +0.58 |
| CoastSat endpoint | 90 | 0.37 | 0.34 | 0.44 | 1.89 | +0.50 |

CoastSat endpoint against CoastSat LRR (two estimators of the same series): r = 0.93, slope = 0.95, RMSE = 0.57, bias = -0.08 m/yr.

## Window occupancy

Median CoastSat observations per transect inside the start window: 20; inside the end window: 39. Transects with an empty window: start 0, end 0, of 906.

| window | asked for | median first obs | median last obs | truncated |
|---|---|---|---|---|
| start | 2003-11-24 to 2004-11-23 | 2003-11-29 | 2004-11-07 | no |
| end | 2022-12-30 to 2023-12-30 | 2023-01-21 | 2023-12-15 | no |

## Sensitivity of the endpoint rate to the assumed survey date

| centre shift | start | end | island mean endpoint (m/yr) | r vs dune | slope | RMSE |
|---|---|---|---|---|---|---|
| -6 mo | 2004-05-25 | 2022-12-30 | 0.63 | 0.34 | 0.30 | 1.98 |
| 0 | 2004-05-25 | 2023-07-01 | 0.41 | 0.37 | 0.34 | 1.89 |
| +6 mo | 2004-05-25 | 2023-12-30 | 0.29 | 0.38 | 0.35 | 1.87 |

## Read this before quoting it

* The dune line and the waterline are different features. A gap between them is beach-width change as much as it is disagreement.
* The LRR spans the calendar window; the endpoint spans the survey interval. They are not the same length of record.
* `n_obs_*` is the median per-transect count inside a one-year window. One storm inside a window moves that end.
