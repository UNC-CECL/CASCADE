# forward_from_1996/domain_means — the unit the model is graded on

The grading target is the domain MEAN of its transect rates
(`coastsat_domain_lrr.py` fits each transect and averages the slopes), so this
is the scale the answer actually has to be given at. `../sites/` and
`../all_transects/` work on single transects, which are noisier and therefore
an upper bound on the convergence window.

No refitting happens here: every window of every transect is already in
`../all_transects/window_convergence_transects_all.csv`, and this is that table
grouped to the 90 domains.

```
domain_mean_vs_transect_forward_from_1996.png   what averaging buys, both panels
window_convergence_domains_forward_from_1996.png
                                     the eight site domains as an ERROR, the
                                     like-for-like against ../sites/
tolerance_comparison_forward_from_1996.png      the five tolerances on one error curve
window_convergence_domains.csv       a row per domain per window
convergence_summary_domains.csv      a row per domain
supporting/                          the PDFs and CAPTIONS.md
```

## Two uncertainties, and the scoring uses the wider

`unc_m_yr` is the MEAN of the domain's transects' 95% half-widths.
`unc_if_independent_m_yr` is what the half-width of the mean would be if those
~10 transects were independent samples. They are not -- they are 10-metre-spaced
views of the same shoreline and move together -- so propagating them that way
divides the band by about √10 and manufactures a later convergence window
out of an assumption. Scoring uses the mean, which is **conservative**: it makes
convergence look earlier, not later. The independent figure is a column so the
size of the choice is visible.

## What this run found

- **Headline (CI overlap).** The rate settles after a median of 26 years (range 6–29 across the 90 domain means), i.e. the window 1996–2021.
- Against the reference fit's 95% band, 10 of 90 domain means (11%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2022.
- Against ±0.25 m/yr, 16 of 90 domain means (18%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2022.
- Against ±20%, 15 of 90 domain means (17%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2022.
- Against 3× the reference band, 28 of 90 domain means (31%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2021.
- Against ±0.50 m/yr, 30 of 90 domain means (33%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2020.
- Against ±1.00 m/yr, 48 of 90 domain means (53%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2014.
- Against CI overlap, 34 of 90 domain means (38%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2021.
- 69 of 90 domain means (77%) enter the band and leave it again before settling, so the FIRST crossing is not the convergence window.
- 26 of 90 domain means (29%) change SIGN between the 1996–2010 window and the reference — erosional over one and accretional over the other.
- The largest 1996–2010 error is GIS 36 at -3.57 m/yr against a reference of +0.01 m/yr.
