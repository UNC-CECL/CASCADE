# window_convergence/record_1996_2020/forward_from_1996 — how much record do you need from 1996?

Fitted on the CoastSat record **1996–2020**. A sweep run on a
different record span lives under a different `record_` folder and is a
different experiment, not a version of this one: every window in it is scored
against a different reference.

The START is pinned at 1996 and the END walks out, one year at a time: 1996–2000 through 1996–2020. It answers how much record the chain needs before the fitted rate stops depending on where it is cut off.

Both directions converge on the same reference, the 1996–2020
rate, and that reference is the longest window of each sweep — fitted in
the same loop as every other window, so it cannot drift from a stored product.
The marked year 2010 is a real model window in both: forward it is
1996–2010, the window the model is graded on, and backward it is
2010–2020, the second leg of the canonical chain.

```
sites/          eight evenly spaced domains, one transect each, in full
all_transects/  every CoastSat transect, aggregated to the 90 domains
```

## The window this sweep gives

**DOMAIN MEANS settle at 1996–2015 (CI overlap, median over 90 domains); 33% have their 1996–2010 rate in that band. Single transects settle at 1996–2016**

- **Headline (CI overlap).** The rate settles after a median of 20 years (range 5–23 across the 90 domain means), i.e. the window 1996–2015.
- Against the reference fit's 95% band, 10 of 90 domain means (11%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2018.
- Against ±0.25 m/yr, 13 of 90 domain means (14%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2017.
- Against ±20%, 13 of 90 domain means (14%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2018.
- Against 3× the reference band, 27 of 90 domain means (30%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2015.
- Against ±0.50 m/yr, 23 of 90 domain means (26%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2015.
- Against ±1.00 m/yr, 48 of 90 domain means (53%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2011.
- Against CI overlap, 30 of 90 domain means (33%) have their 1996–2010 rate already inside the band; the median convergence window is 1996–2015.
- 63 of 90 domain means (70%) enter the band and leave it again before settling, so the FIRST crossing is not the convergence window.
- 21 of 90 domain means (23%) change SIGN between the 1996–2010 window and the reference — erosional over one and accretional over the other.
- The largest 1996–2010 error is GIS 1 at +3.01 m/yr against a reference of +0.22 m/yr.

Producer:
`scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_convergence.py`
`--direction forward`. Built 2026-09-23 by interview (Hannah).
