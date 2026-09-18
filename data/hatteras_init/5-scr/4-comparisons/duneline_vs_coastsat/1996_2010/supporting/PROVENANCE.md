# Dune line vs CoastSat shoreline, 1996-2010

Written 2026-09-15 20:22 by `scripts/input_prep/5-scr/duneline_vs_coastsat/duneline_vs_coastsat.py`.

## Inputs

* dune lines: `2-brie-offset/raw_offsets/1997_duneline_offset_raw.csv`, `2009_duneline_offset_raw.csv` (first row per transect, domain mean, as `hindcast.load_absolute_dune_distance`). Both built by `duneline_to_raw_offsets.py`, so no GIS-vs-shapely metre between them.
* CoastSat LRR: `coastsat_lrr/1996_2010/transect_lrr_full.csv` (window 1996-01-01 to 2010-12-31, per-transect OLS).
* CoastSat endpoint: mean chainage within ±183 days of each survey date, per transect, from `coastsat_timeseries/`.
* transect → domain: `transect_domains/transect_domain_lookup.csv`.

## Survey dates

| line | date | source |
|---|---|---|
| 1996 | 1997-10-12 | USGS 1997 aerial photo, Henderson release (`D:\Hatteras_GIS\Aerial\1997_henderson`, Calendar_Date 19971012) — the 1997 line standing in for 1996 |
| 2010 | 2009-05-30 | Google Earth capture date; Hannah, 2026-09-15 — the 2009 line standing in for 2010 |

Survey interval 11.63 yr. Sign: seaward positive in every column; a negative rate is retreat. Dune change is `-(ORIG_LEN_end - ORIG_LEN_start)`.

## Island-wide means (m/yr)

| dune line | CoastSat LRR | CoastSat endpoint |
|---|---|---|
| -1.30 | -0.32 | -0.17 |

## Agreement, per domain (y against dune rate x)

| y | n | r | slope | intercept | RMSE | bias (y − x) |
|---|---|---|---|---|---|---|
| CoastSat LRR | 90 | 0.18 | 0.14 | -0.14 | 2.97 | +0.98 |
| CoastSat endpoint | 90 | 0.29 | 0.24 | 0.15 | 2.90 | +1.13 |

CoastSat endpoint against CoastSat LRR (two estimators of the same series): r = 0.90, slope = 0.96, RMSE = 0.91, bias = +0.15 m/yr.

## Window occupancy

Median CoastSat observations per transect inside the start window: 6; inside the end window: 9. Transects with an empty window: start 0, end 0, of 906.

| window | asked for | median first obs | median last obs | truncated |
|---|---|---|---|---|
| start | 1997-04-12 to 1998-04-12 | 1997-05-04 | 1998-03-04 | no |
| end | 2008-11-28 to 2009-11-28 | 2008-12-04 | 2009-11-05 | no |

## Read this before quoting it

* The dune line and the waterline are different features. A gap between them is beach-width change as much as it is disagreement.
* The LRR spans the calendar window; the endpoint spans the survey interval. They are not the same length of record.
* `n_obs_*` is the median per-transect count inside a one-year window. One storm inside a window moves that end.
