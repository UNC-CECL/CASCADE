# backward_from_2024/domain_means — the unit the model is graded on

The grading target is the domain MEAN of its transect rates
(`coastsat_domain_lrr.py` fits each transect and averages the slopes), so this
is the scale the answer actually has to be given at. `../sites/` and
`../all_transects/` work on single transects, which are noisier and therefore
an upper bound on the convergence window.

No refitting happens here: every window of every transect is already in
`../all_transects/window_convergence_transects_all.csv`, and this is that table
grouped to the 90 domains.

```
domain_mean_vs_transect_backward_from_2024.png   what averaging buys, both panels
window_convergence_domains_backward_from_2024.png
                                     the eight site domains as an ERROR, the
                                     like-for-like against ../sites/
tolerance_comparison_backward_from_2024.png      the five tolerances on one error curve
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

- **Headline (CI overlap).** The rate settles after a median of 22 years (range 9–26 across the 90 domain means), i.e. the window 2003–2024.
- Against the reference fit's 95% band, 5 of 90 domain means (6%) have their 2010–2024 rate already inside the band; the median convergence window is 2000–2024.
- Against ±0.25 m/yr, 9 of 90 domain means (10%) have their 2010–2024 rate already inside the band; the median convergence window is 2001–2024.
- Against ±20%, 9 of 90 domain means (10%) have their 2010–2024 rate already inside the band; the median convergence window is 2001–2024.
- Against 3× the reference band, 20 of 90 domain means (22%) have their 2010–2024 rate already inside the band; the median convergence window is 2004–2024.
- Against ±0.50 m/yr, 23 of 90 domain means (26%) have their 2010–2024 rate already inside the band; the median convergence window is 2004–2024.
- Against ±1.00 m/yr, 46 of 90 domain means (51%) have their 2010–2024 rate already inside the band; the median convergence window is 2009–2024.
- Against CI overlap, 21 of 90 domain means (23%) have their 2010–2024 rate already inside the band; the median convergence window is 2003–2024.
- 70 of 90 domain means (78%) enter the band and leave it again before settling, so the FIRST crossing is not the convergence window.
- 14 of 90 domain means (16%) change SIGN between the 2010–2024 window and the reference — erosional over one and accretional over the other.
- The largest 2010–2024 error is GIS 1 at +4.05 m/yr against a reference of +2.84 m/yr.
