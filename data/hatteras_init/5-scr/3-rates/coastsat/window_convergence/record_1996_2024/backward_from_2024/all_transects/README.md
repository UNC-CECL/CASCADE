# backward_from_2024/all_transects — the whole island

The same nested sweep on every CoastSat transect in the lookup
(906 of them, 22650 fits), aggregated to the 90 model
domains. It exists to answer one question about `../sites/`: were the eight
representative, or did an even spread of eight happen to pick the unsettled
ones?

```
convergence_alongshore_backward_from_2024.png       three panels on the domain axis
tolerance_comparison_backward_from_2024.png         what "close enough" costs: all five
                                        tolerances on one error curve, and
                                        the island settled-by-record-length
domain_convergence_summary.csv          a row per GIS domain: medians and
                                        quartiles across its ~10 transects
convergence_summary_all_transects.csv   a row per transect
window_convergence_transects_all.csv    the full sweep, 22650 rows
supporting/                             the PDF and CAPTIONS.md
```

The domain number is the **median** across its transects, never the mean: one
transect that never settles would drag a mean to the end of the record and
report that the whole domain behaved that way.

## What this run found

- **Headline (CI overlap).** The rate settles after a median of 22 years (range 5–27 across the 906 transects), i.e. the window 2002–2024.
- Against the reference fit's 95% band, 60 of 906 transects (7%) have their 2010–2024 rate already inside the band; the median convergence window is 2000–2024.
- Against ±0.25 m/yr, 114 of 906 transects (13%) have their 2010–2024 rate already inside the band; the median convergence window is 2001–2024.
- Against ±20%, 70 of 906 transects (8%) have their 2010–2024 rate already inside the band; the median convergence window is 2000–2024.
- Against 3× the reference band, 200 of 906 transects (22%) have their 2010–2024 rate already inside the band; the median convergence window is 2004–2024.
- Against ±0.50 m/yr, 239 of 906 transects (26%) have their 2010–2024 rate already inside the band; the median convergence window is 2004–2024.
- Against ±1.00 m/yr, 451 of 906 transects (50%) have their 2010–2024 rate already inside the band; the median convergence window is 2009–2024.
- Against CI overlap, 227 of 906 transects (25%) have their 2010–2024 rate already inside the band; the median convergence window is 2002–2024.
- 651 of 906 transects (72%) enter the band and leave it again before settling, so the FIRST crossing is not the convergence window.
- 150 of 906 transects (17%) change SIGN between the 2010–2024 window and the reference — erosional over one and accretional over the other.
- The largest 2010–2024 error is GIS 1 (usa_NC_0032_0021) at +5.91 m/yr against a reference of +3.16 m/yr.
