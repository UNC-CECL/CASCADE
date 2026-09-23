# backward_from_2020/domain_means — the unit the model is graded on

The grading target is the domain MEAN of its transect rates
(`coastsat_domain_lrr.py` fits each transect and averages the slopes), so this
is the scale the answer actually has to be given at. `../sites/` and
`../all_transects/` work on single transects, which are noisier and therefore
an upper bound on the convergence window.

No refitting happens here: every window of every transect is already in
`../all_transects/window_convergence_transects_all.csv`, and this is that table
grouped to the 90 domains.

```
domain_mean_vs_transect_backward_from_2020.png   what averaging buys, both panels
window_convergence_domains_backward_from_2020.png
                                     the eight site domains as an ERROR, the
                                     like-for-like against ../sites/
tolerance_comparison_backward_from_2020.png      the five tolerances on one error curve
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

- **Headline (CI overlap).** The rate settles after a median of 18 years (range 6–23 across the 90 domain means), i.e. the window 2002–2020.
- Against the reference fit's 95% band, 9 of 90 domain means (10%) have their 2010–2020 rate already inside the band; the median convergence window is 1999–2020.
- Against ±0.25 m/yr, 12 of 90 domain means (13%) have their 2010–2020 rate already inside the band; the median convergence window is 2000–2020.
- Against ±20%, 8 of 90 domain means (9%) have their 2010–2020 rate already inside the band; the median convergence window is 1999–2020.
- Against 3× the reference band, 21 of 90 domain means (23%) have their 2010–2020 rate already inside the band; the median convergence window is 2004–2020.
- Against ±0.50 m/yr, 18 of 90 domain means (20%) have their 2010–2020 rate already inside the band; the median convergence window is 2002–2020.
- Against ±1.00 m/yr, 40 of 90 domain means (44%) have their 2010–2020 rate already inside the band; the median convergence window is 2006–2020.
- Against CI overlap, 32 of 90 domain means (36%) have their 2010–2020 rate already inside the band; the median convergence window is 2002–2020.
- 65 of 90 domain means (72%) enter the band and leave it again before settling, so the FIRST crossing is not the convergence window.
- 32 of 90 domain means (36%) change SIGN between the 2010–2020 window and the reference — erosional over one and accretional over the other.
- The largest 2010–2020 error is GIS 35 at +4.50 m/yr against a reference of -0.30 m/yr.
