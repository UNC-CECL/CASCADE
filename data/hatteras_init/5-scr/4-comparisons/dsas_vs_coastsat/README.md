# 4-comparisons/dsas_vs_coastsat — do the two rate sources agree?

Built 2026-09-22. The question is whether DSAS (rates from the digitized
shorelines) and CoastSat (rates from the satellite shorelines) say the same
thing about the same 500 m of beach, **before** any smoothing enters.

```
calendar_windows/    CoastSat fitted over the DSAS calendar years.
                     BOTH windows: 1978-1997 and 1997-2019.
    dsas_vs_coastsat_raw.png
    slides/          the 3.4 in version
    supporting/      per-domain table, PDF, captions
survey_dates/        CoastSat anchored on the shoreline survey dates instead.
                     1997-2019 ONLY -- see "Why only one window" below.
    dsas_vs_coastsat_datematched.png    at +/-30 days and +/-6 months
    slides/
    supporting/
```

## Why only one window is date-matched

The date-matched method needs satellite positions around BOTH survey dates.
CoastSat's record here begins **1984-06-17**, so a window around a 1978
shoreline contains no imagery at all — there is nothing to average. The
1978-1997 pair cannot be date-matched by any window width, and its absence
from `survey_dates/` is a property of the satellite record, not an omission.

For the same reason the 1978-1997 panel in `calendar_windows/` is not a
like-for-like comparison: CoastSat's "1978-1997" rate is fitted from
1984-06-17 over 13.5 years against the DSAS 19.

Rebuild with
`python scripts/input_prep/5-scr/4-comparisons/dsas_vs_coastsat/dsas_vs_coastsat_raw.py`.

## What it says

| window | n | bias (CoastSat − DSAS) | RMSE | r | largest gap |
| --- | --- | --- | --- | --- | --- |
| 1978–1997 | 90 | −1.79 m/yr | 2.48 m/yr | 0.69 | 5.73 m/yr at GIS 65 |
| 1997–2019 | 83 | +0.70 m/yr | 0.98 m/yr | 0.90 | 2.67 m/yr at GIS 78 |

**Read the two panels differently.** 1997–2019 is a like-for-like comparison
and the sources agree well. 1978–1997 is not: CoastSat imagery begins in 1984,
so its rate there is fitted from 1984-06-17 over 13.5 years against the DSAS
19, missing the first six entirely. The two panels are not two measurements of
one disagreement — the first is largely a different period, not a different
method, and the 1.79 m/yr bias should not be quoted as a method offset.

Seven domains have no DSAS rate in 1997–2019 (n 83 of 90).

## Date-matched: anchoring CoastSat on the survey dates

`dsas_vs_coastsat_datematched.png` (+ `_slide`) asks the same question a
second way. Instead of fitting CoastSat over a calendar window, it takes the
mean satellite position within ±W days of each shoreline survey date and
differences them — the endpoint method `coastsat_endpoint.py` uses against the
dune line — so both sources describe motion between the same two moments.

Anchors **1997-09-27** and **2019-09-07**, 21.95 years apart.

| window | n domains | bias | RMSE | r | coverage |
| --- | --- | --- | --- | --- | --- |
| ±30 days | 48 | +1.09 m/yr | 1.28 | 0.92 | 457/906 transects; median 1 and 4 positions per end |
| ±6 months | 83 | +0.74 m/yr | 0.89 | **0.95** | 906/906; median 6 and 23 positions |
| (calendar window, for reference) | 83 | +0.70 m/yr | 0.98 | 0.90 | — |

**Date matching does improve the agreement**: r 0.90 → 0.95 and RMSE 0.98 →
0.89 at ±6 months, against the calendar-window comparison on the same
shorelines. Use the ±6-month row. The ±30-day row rests on a median of **one**
satellite position at the 1997 end and covers half the transects, so its
higher-looking r is not better evidence — it is a thinner one.

**On the 1997 date.** The survey date is **1997-09-27**, which also appears
commented into `coastsat_domain_lrr_specific_dates.py`. The 1997-01-01 on all
23 of the 1997 features in `nc_shorelines.geojson` is a placeholder Hannah
entered when the date was not to hand (confirmed 2026-09-22) — not a competing
record, and not to be cited as one. The geojson still carries it; correcting
the source data is a separate decision from using the right date here.

**What the archive did not do.**
`5-scr/archive/coastsat_lrr_superseded_20260810/dsas_coastsat_specific_dates/`
is named for this analysis but was produced with `SURVEY_DATES = []`, so it
fell through to the continuous-range mode: its median `n_obs` (417) matches
the plain calendar fit (414.5), and its statistics are the calendar
comparison's to two decimals. It is a calendar comparison under a
date-matching name, and should not be cited as evidence of date matching.

## Not the smoothed version

`6-scr-smooth/dsas_vs_coastsat/` compares the same two sources with a LOESS on
each, and exists to argue about the smoothing; its figures date from
2026-09-02 and predate the current figure style. This folder is the comparison
without it.

## Where the CoastSat rates come from

**Refitted here** (2026-09-22) from the current time series and the current
`transect_domain_lookup.csv`: all 906 transects, fitted over the DSAS calendar
window with the same loader, date filter, OLS and 3-position minimum as the
live windows. Only the window differs from `3-rates/coastsat/lrr/`.

They were read from `5-scr/archive/coastsat_lrr_superseded_20260810/` until
that date. The archived fits are still loaded, and the table carries them
beside the refit so the two can be differenced: **mean |difference| 0.050 m/yr
(1978–1997) and 0.039 (1997–2019)**, so the archive was close at domain level,
but the largest single domain moved 3.10 and 2.88 m/yr. The refit lifts r for
1978–1997 from 0.64 to 0.69 and leaves 1997–2019 unchanged.

**The refit is deliberately not written into `3-rates/coastsat/lrr/`.**
`hat_observed_rates.windows()` builds its window list by scanning that
directory, so a `1978_1997/` folder there would become a window for every
caller — `rates_figures.py` would draw figures for it, and the shared y bound
of every window figure is the largest domain mean over all windows, so the
existing figures would change. These years are a comparison, not a rate
product. No run is graded against them; the scoring target is
`3-rates/coastsat/lrr/`.
