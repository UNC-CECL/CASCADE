# Dune line vs CoastSat shoreline, 1996-2024 (net change)

Written 2026-09-22 09:14 by `scripts/input_prep/5-scr/duneline_vs_coastsat/duneline_vs_coastsat.py`.

**Both sides are NET CHANGE between the same two dates** (2026-09-18, Hannah: the comparison is net position change on both sides). The CoastSat LRR, which this folder also drew until then, is not a like-for-like quantity for two surveys; it stays the model's scoring target in `3-rates/coastsat/lrr/`.

## Inputs (the stored products)

* dune line: `3-rates/duneline/endpoint/1996_2024/` (the 1997 and 2023 lines).
* CoastSat shoreline: `3-rates/coastsat/endpoint/1996_2024/` (mean position within ±6 months of each dune-line date, per transect, domain mean).

## Survey dates

| line | date | source |
|---|---|---|
| 1996 | 1997-10-12 | `duneline_vs_coastsat.KNOWN_SURVEY_DATES` — the 1997 line standing in for 1996 |
| 2024 | 2023-07-01 | **ASSUMED 1 July**; no flight date known for this line — the 2023 line standing in for 2024 |

Survey interval 25.72 yr. Seaward positive in every column.

## Island-wide means

| | net change (m) | as a rate (m/yr) |
|---|---|---|
| dune line | -14.8 | -0.58 |
| CoastSat shoreline | +5.0 | +0.20 |

## Agreement, per domain (shoreline y against dune x, m/yr)

| n | r | slope | intercept | RMSE | bias (y − x) |
|---|---|---|---|---|---|
| 90 | 0.81 | 0.94 | 0.74 | 1.16 | +0.77 |

## Window occupancy (CoastSat)

Median positions per transect inside the start window 6, the end window 39; empty windows: start 0, end 0, of 906 transects.

| window | asked for | median first obs | median last obs | truncated |
|---|---|---|---|---|
| start | 1997-04-12 to 1998-04-12 | 1997-05-04 | 1998-03-04 | no |
| end | 2022-12-30 to 2023-12-30 | 2023-01-21 | 2023-12-15 | no |

## Sensitivity of the CoastSat net change to the assumed date

| centre shift | start | end | island mean (m/yr) | r vs dune | slope | RMSE |
|---|---|---|---|---|---|---|
| -6 mo | 1997-10-12 | 2022-12-30 | 0.35 | 0.80 | 0.91 | 1.27 |
| 0 | 1997-10-12 | 2023-07-01 | 0.20 | 0.81 | 0.94 | 1.16 |
| +6 mo | 1997-10-12 | 2023-12-30 | 0.11 | 0.79 | 0.95 | 1.14 |

## Read this before quoting it

* The dune line and the waterline are different features; a gap between them is beach-width change as much as disagreement.
* One storm inside a ±6-month window moves that end.
