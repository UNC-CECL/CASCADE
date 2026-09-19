# Dune line vs CoastSat shoreline, 2010-2024 (net change)

Written 2026-09-19 11:19 by `scripts/input_prep/5-scr/duneline_vs_coastsat/duneline_vs_coastsat.py`.

**Both sides are NET CHANGE between the same two dates** (2026-09-18, Hannah: the comparison is net position change on both sides). The CoastSat LRR, which this folder also drew until then, is not a like-for-like quantity for two surveys; it stays the model's scoring target in `3-rates/coastsat/lrr/`.

## Inputs (the stored products)

* dune line: `3-rates/duneline/endpoint/2010_2024/` (the 2009 and 2023 lines).
* CoastSat shoreline: `3-rates/coastsat/endpoint/2010_2024/` (mean position within ±6 months of each dune-line date, per transect, domain mean).

## Survey dates

| line | date | source |
|---|---|---|
| 2010 | 2009-05-30 | `duneline_vs_coastsat.KNOWN_SURVEY_DATES` — the 2009 line standing in for 2010 |
| 2024 | 2023-07-01 | **ASSUMED 1 July**; no flight date known for this line — the 2023 line standing in for 2024 |

Survey interval 14.09 yr. Seaward positive in every column.

## Island-wide means

| | net change (m) | as a rate (m/yr) |
|---|---|---|
| dune line | +1.4 | +0.10 |
| CoastSat shoreline | +7.0 | +0.50 |

## Agreement, per domain (shoreline y against dune x, m/yr)

| n | r | slope | intercept | RMSE | bias (y − x) |
|---|---|---|---|---|---|
| 90 | 0.71 | 0.91 | 0.41 | 1.33 | +0.40 |

## Window occupancy (CoastSat)

Median positions per transect inside the start window 9, the end window 39; empty windows: start 0, end 0, of 906 transects.

| window | asked for | median first obs | median last obs | truncated |
|---|---|---|---|---|
| start | 2008-11-28 to 2009-11-28 | 2008-12-04 | 2009-11-05 | no |
| end | 2022-12-30 to 2023-12-30 | 2023-01-21 | 2023-12-15 | no |

## Sensitivity of the CoastSat net change to the assumed date

| centre shift | start | end | island mean (m/yr) | r vs dune | slope | RMSE |
|---|---|---|---|---|---|---|
| -6 mo | 2009-05-30 | 2022-12-30 | 0.80 | 0.60 | 0.77 | 1.63 |
| 0 | 2009-05-30 | 2023-07-01 | 0.50 | 0.71 | 0.91 | 1.33 |
| +6 mo | 2009-05-30 | 2023-12-30 | 0.33 | 0.73 | 0.92 | 1.25 |

## Read this before quoting it

* The dune line and the waterline are different features; a gap between them is beach-width change as much as disagreement.
* One storm inside a ±6-month window moves that end.
