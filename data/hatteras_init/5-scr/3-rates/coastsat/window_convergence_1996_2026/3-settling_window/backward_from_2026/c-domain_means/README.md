# backward_from_2026/c-domain_means — the unit the model is graded on

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

- **Headline (CI overlap).** The rate settles after a median of 25 years (range 10–28 across the 90 domain means), i.e. the window 2002–2026.
- Against the reference fit's 95% band, 2 of 90 domain means (2%) have their 2010–2026 rate already inside the band; the median convergence window is 2000–2026.
- Against ±0.25 m/yr, 12 of 90 domain means (13%) have their 2010–2026 rate already inside the band; the median convergence window is 2001–2026.
- Against ±20%, 8 of 90 domain means (9%) have their 2010–2026 rate already inside the band; the median convergence window is 2001–2026.
- Against 3× the reference band, 19 of 90 domain means (21%) have their 2010–2026 rate already inside the band; the median convergence window is 2004–2026.
- Against ±0.50 m/yr, 26 of 90 domain means (29%) have their 2010–2026 rate already inside the band; the median convergence window is 2004–2026.
- Against ±1.00 m/yr, 46 of 90 domain means (51%) have their 2010–2026 rate already inside the band; the median convergence window is 2010–2026.
- Against CI overlap, 22 of 90 domain means (24%) have their 2010–2026 rate already inside the band; the median convergence window is 2002–2026.
- 70 of 90 domain means (78%) enter the band and leave it again before settling, so the FIRST crossing is not the convergence window.
- 13 of 90 domain means (14%) change SIGN between the 2010–2026 window and the reference — erosional over one and accretional over the other.
- The largest 2010–2026 error is GIS 1 at +5.95 m/yr against a reference of +4.08 m/yr.
