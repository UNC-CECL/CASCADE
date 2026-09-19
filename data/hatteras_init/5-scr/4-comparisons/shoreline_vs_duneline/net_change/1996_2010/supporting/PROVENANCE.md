# Dune line vs CoastSat shoreline, 1996-2010 (net change)

Written 2026-09-19 11:18 by `scripts/input_prep/5-scr/duneline_vs_coastsat/duneline_vs_coastsat.py`.

**Both sides are NET CHANGE between the same two dates** (2026-09-18, Hannah: the comparison is net position change on both sides). The CoastSat LRR, which this folder also drew until then, is not a like-for-like quantity for two surveys; it stays the model's scoring target in `3-rates/coastsat/lrr/`.

## Inputs (the stored products)

* dune line: `3-rates/duneline/endpoint/1996_2010/` (the 1997 and 2009 lines).
* CoastSat shoreline: `3-rates/coastsat/endpoint/1996_2010/` (mean position within ±6 months of each dune-line date, per transect, domain mean).

## Survey dates

| line | date | source |
|---|---|---|
| 1996 | 1997-10-12 | `duneline_vs_coastsat.KNOWN_SURVEY_DATES` — the 1997 line standing in for 1996 |
| 2010 | 2009-05-30 | `duneline_vs_coastsat.KNOWN_SURVEY_DATES` — the 2009 line standing in for 2010 |

Survey interval 11.63 yr. Seaward positive in every column.

## Island-wide means

| | net change (m) | as a rate (m/yr) |
|---|---|---|
| dune line | -16.2 | -1.39 |
| CoastSat shoreline | -2.0 | -0.17 |

## Agreement, per domain (shoreline y against dune x, m/yr)

| n | r | slope | intercept | RMSE | bias (y − x) |
|---|---|---|---|---|---|
| 90 | 0.66 | 0.73 | 0.85 | 2.02 | +1.22 |

## Window occupancy (CoastSat)

Median positions per transect inside the start window 6, the end window 9; empty windows: start 0, end 0, of 906 transects.

| window | asked for | median first obs | median last obs | truncated |
|---|---|---|---|---|
| start | 1997-04-12 to 1998-04-12 | 1997-05-04 | 1998-03-04 | no |
| end | 2008-11-28 to 2009-11-28 | 2008-12-04 | 2009-11-05 | no |

## Read this before quoting it

* The dune line and the waterline are different features; a gap between them is beach-width change as much as disagreement.
* One storm inside a ±6-month window moves that end.
