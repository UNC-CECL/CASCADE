# backward_from_2020/all_transects — the whole island

The same nested sweep on every CoastSat transect in the lookup
(906 of them, 19026 fits), aggregated to the 90 model
domains. It exists to answer one question about `../sites/`: were the eight
representative, or did an even spread of eight happen to pick the unsettled
ones?

```
convergence_alongshore_backward_from_2020.png       three panels on the domain axis
tolerance_comparison_backward_from_2020.png         what "close enough" costs: all five
                                        tolerances on one error curve, and
                                        the island settled-by-record-length
domain_convergence_summary.csv          a row per GIS domain: medians and
                                        quartiles across its ~10 transects
convergence_summary_all_transects.csv   a row per transect
window_convergence_transects_all.csv    the full sweep, 19026 rows
supporting/                             the PDF and CAPTIONS.md
```

The domain number is the **median** across its transects, never the mean: one
transect that never settles would drag a mean to the end of the record and
report that the whole domain behaved that way.

## What this run found

- **Headline (CI overlap).** The rate settles after a median of 19 years (range 5–23 across the 906 transects), i.e. the window 2002–2020.
- Against the reference fit's 95% band, 71 of 906 transects (8%) have their 2010–2020 rate already inside the band; the median convergence window is 1999–2020.
- Against ±0.25 m/yr, 90 of 906 transects (10%) have their 2010–2020 rate already inside the band; the median convergence window is 2000–2020.
- Against ±20%, 73 of 906 transects (8%) have their 2010–2020 rate already inside the band; the median convergence window is 1999–2020.
- Against 3× the reference band, 215 of 906 transects (24%) have their 2010–2020 rate already inside the band; the median convergence window is 2003–2020.
- Against ±0.50 m/yr, 189 of 906 transects (21%) have their 2010–2020 rate already inside the band; the median convergence window is 2002–2020.
- Against ±1.00 m/yr, 351 of 906 transects (39%) have their 2010–2020 rate already inside the band; the median convergence window is 2006–2020.
- Against CI overlap, 298 of 906 transects (33%) have their 2010–2020 rate already inside the band; the median convergence window is 2002–2020.
- 739 of 906 transects (82%) enter the band and leave it again before settling, so the FIRST crossing is not the convergence window.
- 319 of 906 transects (35%) change SIGN between the 2010–2020 window and the reference — erosional over one and accretional over the other.
- The largest 2010–2020 error is GIS 81 (usa_NC_0036_0014) at +5.58 m/yr against a reference of -3.48 m/yr.
