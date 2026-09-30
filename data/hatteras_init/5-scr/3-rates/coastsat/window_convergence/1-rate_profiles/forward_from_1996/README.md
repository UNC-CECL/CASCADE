# 1-rate_profiles/forward_from_1996 — when does the whole profile start to look like 1996–2024?

Every window of the nested family (1996–1997 … 1996–2024, the start pinned at
1996), drawn as a shoreline change rate profile along all 906 CoastSat
transects over the 1996–2024 reference, and scored by the alongshore Pearson r of
each window's profile against it. Built 2026-09-29 by interview (Hannah).

The companion to `../../2-settling_window/forward_from_1996/`, which scores each
location separately against tolerances, from five years. This one draws the
whole profile from **two** years, the minimum Hannah asked for.

```
window_profiles_forward_from_1996.png         (a) every window over the reference, (b) r per window
window_profiles_panels_forward_from_1996.png  one panel per window
window_profiles_transects.csv                     a row per transect per window (lrr, unc, n_obs, x_domain, diff)
window_profiles_correlation.csv                   a row per window: r_vs_reference, n_transects
```

## r against the reference

| window | years | r |
|---|---|---|
| 1996–1997 | 2 | -0.018 |
| 1996–1998 | 3 | 0.042 |
| 1996–1999 | 4 | -0.034 |
| 1996–2000 | 5 | 0.103 |
| 1996–2001 | 6 | 0.223 |
| 1996–2002 | 7 | 0.334 |
| 1996–2003 | 8 | 0.402 |
| 1996–2004 | 9 | 0.458 |
| 1996–2005 | 10 | 0.512 |
| 1996–2006 | 11 | 0.569 |
| 1996–2007 | 12 | 0.604 |
| 1996–2008 | 13 | 0.643 |
| 1996–2009 | 14 | 0.657 |
| 1996–2010 ← model window | 15 | 0.695 |
| 1996–2011 | 16 | 0.719 |
| 1996–2012 | 17 | 0.729 |
| 1996–2013 | 18 | 0.761 |
| 1996–2014 | 19 | 0.785 |
| 1996–2015 | 20 | 0.820 |
| 1996–2016 | 21 | 0.848 |
| 1996–2017 | 22 | 0.871 |
| 1996–2018 | 23 | 0.900 |
| 1996–2019 | 24 | 0.917 |
| 1996–2020 | 25 | 0.930 |
| 1996–2021 | 26 | 0.958 |
| 1996–2022 | 27 | 0.977 |
| 1996–2023 | 28 | 0.993 |
| 1996–2024 | 29 | 1.000 |

**Read r with care.** The windows are nested, so r reaches 1 at the
reference by construction. What it shows is where it gets there and how
steadily. r measures shape only: a window with every hotspot in the right
place at the wrong magnitude still scores high. Magnitude is scored by the
bias and tolerance tables in `../../2-settling_window/forward_from_1996/b-every_transect/`.

Producer:
`scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_profiles.py`
`--direction forward`.
