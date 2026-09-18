# Dune line vs CoastSat shoreline, 1984-2004

Written 2026-09-15 17:26 by `scripts/input_prep/5-scr/duneline_vs_coastsat/duneline_vs_coastsat.py`.

## Inputs

* dune lines: `2-brie-offset/raw_offsets/1984_duneline_offset_raw.csv`, `2004_duneline_offset_raw.csv` (first row per transect, domain mean, as `hindcast.load_absolute_dune_distance`). Both built by `duneline_to_raw_offsets.py`, so no GIS-vs-shapely metre between them.
* CoastSat LRR: `3-rates/coastsat/lrr/1984_2004/transect_lrr_full.csv` (window 1984-01-01 to 2004-12-31, per-transect OLS).
* CoastSat endpoint: mean chainage within ±183 days of each survey date, per transect, from `coastsat_timeseries/`.
* transect → domain: `transect_domains/transect_domain_lookup.csv`.

## Survey dates

| line | date | source |
|---|---|---|
| 1984 | 1984-09-19 | USGS 1984 aerial photo, Henderson release (`D:\Hatteras_GIS\Aerial\1984_henderson\1984_metadata`) |
| 2004 | 2004-05-25 | Google Earth capture date (the raw_GE frames carry none); Hannah, 2026-09-15 |

Survey interval 19.68 yr. Sign: seaward positive in every column; a negative rate is retreat. Dune change is `-(ORIG_LEN_end - ORIG_LEN_start)`.

## Island-wide means (m/yr)

| dune line | CoastSat LRR | CoastSat endpoint |
|---|---|---|
| -0.78 | -1.13 | -0.99 |

## Agreement, per domain (y against dune rate x)

| y | n | r | slope | intercept | RMSE | bias (y − x) |
|---|---|---|---|---|---|---|
| CoastSat LRR | 90 | 0.44 | 0.41 | -0.81 | 2.07 | -0.35 |
| CoastSat endpoint | 90 | 0.50 | 0.49 | -0.60 | 1.98 | -0.21 |

CoastSat endpoint against CoastSat LRR (two estimators of the same series): r = 0.94, slope = 0.98, RMSE = 0.68, bias = +0.14 m/yr.

## Window occupancy

Median CoastSat observations per transect inside the start window: 7; inside the end window: 20. Transects with an empty window: start 0, end 0, of 906.

| window | asked for | median first obs | median last obs | truncated |
|---|---|---|---|---|
| start | 1984-03-20 to 1985-03-20 | 1984-09-21 | 1985-02-28 | **yes** |
| end | 2003-11-24 to 2004-11-23 | 2003-11-29 | 2004-11-07 | no |

A truncated window does not average a full seasonal cycle. The CoastSat record begins 1984-09-21 on most transects (1984-05/06 on a few), so a window centred on the 1984-09-19 photo holds only the autumn and winter after it.

## Read this before quoting it

* The dune line and the waterline are different features. A gap between them is beach-width change as much as it is disagreement.
* The LRR spans the calendar window; the endpoint spans the survey interval. They are not the same length of record.
* `n_obs_*` is the median per-transect count inside a one-year window. One storm inside a window moves that end.
