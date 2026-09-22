# Dune line vs CoastSat shoreline, 1984-2004 (net change)

Written 2026-09-22 09:13 by `scripts/input_prep/5-scr/4-comparisons/shoreline_vs_duneline/coastsat_vs_duneline.py`.

**Both sides are NET CHANGE between the same two dates** (2026-09-18, Hannah: the comparison is net position change on both sides). The CoastSat LRR, which this folder also drew until then, is not a like-for-like quantity for two surveys; it stays the model's scoring target in `3-rates/coastsat/lrr/`.

## Inputs (the stored products)

* dune line: `3-rates/duneline/endpoint/1984_2004/` (the 1984 and 2004 lines).
* CoastSat shoreline: `3-rates/coastsat/endpoint/1984_2004/` (mean position within ±6 months of each dune-line date, per transect, domain mean).

## Survey dates

| line | date | source |
|---|---|---|
| 1984 | 1984-09-19 | `coastsat_vs_duneline.KNOWN_SURVEY_DATES` |
| 2004 | 2004-05-25 | `coastsat_vs_duneline.KNOWN_SURVEY_DATES` |

Survey interval 19.68 yr. Seaward positive in every column.

## Island-wide means

| | net change (m) | as a rate (m/yr) |
|---|---|---|
| dune line | -15.4 | -0.78 |
| CoastSat shoreline | -19.4 | -0.99 |

## Agreement, per domain (shoreline y against dune x, m/yr)

| n | r | slope | intercept | RMSE | bias (y − x) |
|---|---|---|---|---|---|
| 90 | 0.50 | 0.49 | -0.60 | 1.98 | -0.21 |

## Window occupancy (CoastSat)

Median positions per transect inside the start window 7, the end window 20; empty windows: start 0, end 0, of 906 transects.

| window | asked for | median first obs | median last obs | truncated |
|---|---|---|---|---|
| start | 1984-03-20 to 1985-03-20 | 1984-09-21 | 1985-02-28 | **yes** |
| end | 2003-11-24 to 2004-11-23 | 2003-11-29 | 2004-11-07 | no |

A truncated window does not average a full seasonal cycle. The CoastSat record begins 1984-09-21 on most transects, so a window centred on the 1984-09-19 photo holds only the autumn and winter after it.

## Read this before quoting it

* The dune line and the waterline are different features; a gap between them is beach-width change as much as disagreement.
* One storm inside a ±6-month window moves that end.
