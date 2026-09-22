# 4-comparisons/dsas_vs_coastsat — do the two rate sources agree?

Built 2026-09-22. The question is whether DSAS (rates from the digitized
shorelines) and CoastSat (rates from the satellite shorelines) say the same
thing about the same 500 m of beach, **before** any smoothing enters.

```
dsas_vs_coastsat_raw.png       two panels, one per DSAS window: the plain
                               per-domain mean rate from each source, no LOESS
supporting/
    dsas_vs_coastsat_raw.csv   per domain and window: DSAS, the CoastSat
                               refit, their difference, and the archived
                               CoastSat fit the refit replaced
    dsas_vs_coastsat_raw.pdf
    CAPTIONS.md
```

Rebuild with
`python scripts/input_prep/5-scr/dsas_vs_coastsat/dsas_vs_coastsat_raw.py`.

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
