# forward_from_1996/b-every_transect — the whole island

The same nested sweep on every CoastSat transect in the lookup
(906 of them, 24462 fits), aggregated to the 90 model
domains. It exists to answer one question about `../a-eight_sites/`: were the eight
representative, or did an even spread of eight happen to pick the unsettled
ones?

The figure drawn from these tables is `../../years_needed_alongshore.png`.

```
domain_convergence_summary.csv          a row per GIS domain: medians and
                                        quartiles across its ~10 transects
convergence_summary_all_transects.csv   a row per transect
window_convergence_transects_all.csv    the full sweep, 24462 rows
```

`years_needed_abs`, `_abs50` and `_abs100` are the ±0.25, ±0.5 and ±1.0 m/yr
columns the figure draws.

The domain number is the **median** across its transects, never the mean: one
transect that never settles would drag a mean to the end of the record and
report that the whole domain behaved that way.

## What this run found

- **Headline (CI overlap).** The rate settles after a median of 27 years (range 5–30 across the 906 transects), i.e. the window 1996–2022.
- Against the reference fit's 95% band, 67 of 906 transects (7%) have their 1996–2015 rate already inside the band; the median convergence window is 1996–2024.
- Against ±0.25 m/yr, 141 of 906 transects (16%) have their 1996–2015 rate already inside the band; the median convergence window is 1996–2022.
- Against ±20%, 120 of 906 transects (13%) have their 1996–2015 rate already inside the band; the median convergence window is 1996–2023.
- Against 3× the reference band, 221 of 906 transects (24%) have their 1996–2015 rate already inside the band; the median convergence window is 1996–2021.
- Against ±0.50 m/yr, 279 of 906 transects (31%) have their 1996–2015 rate already inside the band; the median convergence window is 1996–2021.
- Against ±1.00 m/yr, 580 of 906 transects (64%) have their 1996–2015 rate already inside the band; the median convergence window is 1996–2015.
- Against CI overlap, 238 of 906 transects (26%) have their 1996–2015 rate already inside the band; the median convergence window is 1996–2022.
- 686 of 906 transects (76%) enter the band and leave it again before settling, so the FIRST crossing is not the convergence window.
- 236 of 906 transects (26%) change SIGN between the 1996–2015 window and the reference — erosional over one and accretional over the other.
- The largest 1996–2015 error is GIS 1 (usa_NC_0032_0021) at -3.67 m/yr against a reference of +4.81 m/yr.
