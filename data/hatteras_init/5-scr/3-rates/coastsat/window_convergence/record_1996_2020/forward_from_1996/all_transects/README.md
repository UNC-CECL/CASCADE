# forward_from_1996/all_transects — the whole island

The same nested sweep on every CoastSat transect in the lookup
(906 of them, 19026 fits), aggregated to the 90 model
domains. It exists to answer one question about `../sites/`: were the eight
representative, or did an even spread of eight happen to pick the unsettled
ones?

```
convergence_alongshore_forward_from_1996.png       three panels on the domain axis
tolerance_comparison_forward_from_1996.png         what "close enough" costs: all five
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

- **Headline (CI overlap).** The rate settles after a median of 21 years (range 5–24 across the 906 transects), i.e. the window 1996–2016.
- Against the reference fit's 95% band, 86 of 906 transects (9%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2018.
- Against ±0.25 m/yr, 109 of 906 transects (12%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2017.
- Against ±20%, 94 of 906 transects (10%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2018.
- Against 3× the reference band, 252 of 906 transects (28%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2015.
- Against ±0.50 m/yr, 209 of 906 transects (23%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2015.
- Against ±1.00 m/yr, 445 of 906 transects (49%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2011.
- Against CI overlap, 274 of 906 transects (30%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2016.
- 638 of 906 transects (70%) enter the band and leave it again before settling, so the FIRST crossing is not the convergence window.
- 216 of 906 transects (24%) change SIGN between the 1996–2010 window and the reference — erosional over one and accretional over the other.
- The largest 1996–2010 error is GIS 81 (usa_NC_0036_0014) at -3.68 m/yr against a reference of -3.48 m/yr.
