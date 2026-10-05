# forward_from_1996/c-domain_means — the unit the model is graded on

The grading target is the domain MEAN of its transect rates
(`coastsat_domain_lrr.py` fits each transect and averages the slopes), so this
is the scale the answer actually has to be given at. `../a-eight_sites/` and
`../b-every_transect/` work on single transects. Averaging was expected to
settle sooner; it does not, so the disagreement between windows is real
shoreline behaviour rather than per-transect scatter.

No refitting happens here: every window of every transect is already in
`../b-every_transect/window_convergence_transects_all.csv`, and this is that table
grouped to the 90 domains.

```
window_convergence_domains.csv       a row per domain per window
convergence_summary_domains.csv      a row per domain
```

Tables only. The figure is `../../years_needed_alongshore.png`.

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

- **Headline (CI overlap).** The rate settles after a median of 26 years (range 6–30 across the 90 domain means), i.e. the window 1996–2021.
- Against the reference fit's 95% band, 7 of 90 domain means (8%) have their 1996–2015 rate already inside the band; the median convergence window is 1996–2024.
- Against ±0.25 m/yr, 14 of 90 domain means (16%) have their 1996–2015 rate already inside the band; the median convergence window is 1996–2022.
- Against ±20%, 11 of 90 domain means (12%) have their 1996–2015 rate already inside the band; the median convergence window is 1996–2023.
- Against 3× the reference band, 20 of 90 domain means (22%) have their 1996–2015 rate already inside the band; the median convergence window is 1996–2021.
- Against ±0.50 m/yr, 28 of 90 domain means (31%) have their 1996–2015 rate already inside the band; the median convergence window is 1996–2021.
- Against ±1.00 m/yr, 62 of 90 domain means (69%) have their 1996–2015 rate already inside the band; the median convergence window is 1996–2014.
- Against CI overlap, 22 of 90 domain means (24%) have their 1996–2015 rate already inside the band; the median convergence window is 1996–2021.
- 70 of 90 domain means (78%) enter the band and leave it again before settling, so the FIRST crossing is not the convergence window.
- 22 of 90 domain means (24%) change SIGN between the 1996–2015 window and the reference — erosional over one and accretional over the other.
- The largest 1996–2015 error is GIS 11 at -2.65 m/yr against a reference of +0.89 m/yr.
