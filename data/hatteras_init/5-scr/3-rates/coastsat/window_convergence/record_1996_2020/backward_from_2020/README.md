# window_convergence/record_1996_2020/backward_from_2020 — how late can a window begin?

Fitted on the CoastSat record **1996–2020**. A sweep run on a
different record span lives under a different `record_` folder and is a
different experiment, not a version of this one: every window in it is scored
against a different reference.

The END is pinned at 2020 and the START walks back, one year at a time: 2016–2020 through 1996–2020. It answers how recent a window can be and still recover the long-term rate — the other bracket on the same question.

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

**DOMAIN MEANS settle at 2002–2020 (CI overlap, median over 90 domains); 36% have their 2010–2020 rate in that band. Single transects settle at 2002–2020**

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

Producer:
`scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_convergence.py`
`--direction backward`. Built 2026-09-23 by interview (Hannah).
