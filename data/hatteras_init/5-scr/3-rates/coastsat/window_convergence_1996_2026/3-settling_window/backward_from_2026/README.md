# backward_from_2026 — how late can a window begin?

**When does each location settle on the long-term rate?** Fitted on the
CoastSat record **1996–2026**. The whole-profile companion questions are in `../../1-rate_profiles/` and `../../2-r_bias_rmse/`.

The END is pinned at 2026 and the START walks back, one year at a time: 2022–2026 through 1996–2026. It answers how recent a window can be and still recover the long-term rate — the other bracket on the same question.

Both directions converge on the same reference, the 1996–2026
rate, and that reference is the longest window of each sweep — fitted in
the same loop as every other window, so it cannot drift from a stored product.
Each direction marks a real model window: forward it is 1996–2015, the
first leg of the model chain, and backward it is 2010–2026, the second leg.

The figure for this sweep is `../years_needed_alongshore.png` (both
directions, one panel each).

```
a-eight_sites/     eight transects: shoreline position through time with the
                   fits drawn on it (the one figure here), plus tables
b-every_transect/  tables for all ~906 transects (the figure's source)
c-domain_means/    tables for the transects averaged to the 90 domains
```

The numbers below are the domain means under seven tolerances ("settles" = the
window's rate enters the tolerance and stays inside it for every longer window).
The figure uses only the three plain ones, ±0.25, ±0.5 and ±1.0 m/yr.

## The window this sweep gives

**DOMAIN MEANS settle at 2002–2026 (CI overlap, median over 90 domains); 24% have their 2010–2026 rate in that band. Single transects settle at 2002–2026**

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

Producer:
`scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_convergence.py`
`--direction backward`. Built 2026-10-02 by interview (Hannah).
