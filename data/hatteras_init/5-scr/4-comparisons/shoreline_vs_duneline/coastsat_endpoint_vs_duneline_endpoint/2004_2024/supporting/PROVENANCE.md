# Dune line vs CoastSat shoreline, 2004-2024 (net change)

Written 2026-09-22 09:14 by `scripts/input_prep/5-scr/4-comparisons/shoreline_vs_duneline/duneline_vs_coastsat.py`.

**Both sides are NET CHANGE between the same two dates** (2026-09-18, Hannah: the comparison is net position change on both sides). The CoastSat LRR, which this folder also drew until then, is not a like-for-like quantity for two surveys; it stays the model's scoring target in `3-rates/coastsat/lrr/`.

## Inputs (the stored products)

* dune line: `3-rates/duneline/endpoint/2004_2024/` (the 2004 and 2023 lines).
* CoastSat shoreline: `3-rates/coastsat/endpoint/2004_2024/` (mean position within ±6 months of each dune-line date, per transect, domain mean).

## Survey dates

| line | date | source |
|---|---|---|
| 2004 | 2004-05-25 | `duneline_vs_coastsat.KNOWN_SURVEY_DATES` |
| 2024 | 2023-07-01 | **ASSUMED 1 July**; no flight date known for this line — the 2023 line standing in for 2024 |

Survey interval 19.10 yr. Seaward positive in every column.

## Island-wide means

| | net change (m) | as a rate (m/yr) |
|---|---|---|
| dune line | -1.7 | -0.09 |
| CoastSat shoreline | +7.8 | +0.41 |

## Agreement, per domain (shoreline y against dune x, m/yr)

| n | r | slope | intercept | RMSE | bias (y − x) |
|---|---|---|---|---|---|
| 90 | 0.37 | 0.34 | 0.44 | 1.89 | +0.50 |

## Window occupancy (CoastSat)

Median positions per transect inside the start window 20, the end window 39; empty windows: start 0, end 0, of 906 transects.

| window | asked for | median first obs | median last obs | truncated |
|---|---|---|---|---|
| start | 2003-11-24 to 2004-11-23 | 2003-11-29 | 2004-11-07 | no |
| end | 2022-12-30 to 2023-12-30 | 2023-01-21 | 2023-12-15 | no |

## Sensitivity of the CoastSat net change to the assumed date

| centre shift | start | end | island mean (m/yr) | r vs dune | slope | RMSE |
|---|---|---|---|---|---|---|
| -6 mo | 2004-05-25 | 2022-12-30 | 0.63 | 0.34 | 0.30 | 1.98 |
| 0 | 2004-05-25 | 2023-07-01 | 0.41 | 0.37 | 0.34 | 1.89 |
| +6 mo | 2004-05-25 | 2023-12-30 | 0.29 | 0.38 | 0.35 | 1.87 |

## Read this before quoting it

* The dune line and the waterline are different features; a gap between them is beach-width change as much as disagreement.
* One storm inside a ±6-month window moves that end.
