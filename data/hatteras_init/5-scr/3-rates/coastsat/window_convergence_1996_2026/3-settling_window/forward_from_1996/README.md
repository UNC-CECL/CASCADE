# forward_from_1996 — how much record do you need from 1996?

**When does each location settle on the long-term rate?** Fitted on the
CoastSat record **1996–2026**. The whole-profile companion questions are in `../../1-rate_profiles/` and `../../2-r_bias_rmse/`.

The START is pinned at 1996 and the END walks out, one year at a time: 1996–2000 through 1996–2026. It answers how much record the chain needs before the fitted rate stops depending on where it is cut off.

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

**DOMAIN MEANS settle at 1996–2021 (CI overlap, median over 90 domains); 24% have their 1996–2015 rate in that band. Single transects settle at 1996–2022**

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

Producer:
`scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_convergence.py`
`--direction forward`. Built 2026-10-02 by interview (Hannah).
