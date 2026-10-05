# backward_from_2024 — how late can a window begin?

**When does each location settle on the long-term rate?** Fitted on the
CoastSat record **1996–2024**. The whole-profile companion questions are in `../../1-rate_profiles/` and `../../2-r_bias_rmse/`.

The END is pinned at 2024 and the START walks back, one year at a time: 2020–2024 through 1996–2024. It answers how recent a window can be and still recover the long-term rate — the other bracket on the same question.

Both directions converge on the same reference, the 1996–2024
rate, and that reference is the longest window of each sweep — fitted in
the same loop as every other window, so it cannot drift from a stored product.
The marked year 2010 is a real model window in both: forward it is
1996–2010, the window the model is graded on, and backward it is
2010–2024, the second leg of the canonical chain.

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

**DOMAIN MEANS settle at 2003–2024 (CI overlap, median over 90 domains); 23% have their 2010–2024 rate in that band. Single transects settle at 2002–2024**

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

Producer:
`scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_convergence.py`
`--direction backward`. Built 2026-09-29 by interview (Hannah).
