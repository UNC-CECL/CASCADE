# 1-rate_profiles/backward_from_2024 — when does the whole profile start to look like 1996–2024?

Every window of the nested family (2023–2024 … 1996–2024, the end pinned at
2024), drawn as a shoreline change rate profile along all 906 CoastSat
transects over the 1996–2024 reference, and scored by the alongshore Pearson r of
each window's profile against it. Built 2026-09-29 by interview (Hannah).

The companion to `../../2-settling_window/backward_from_2024/`, which scores each
location separately against tolerances, from five years. This one draws the
whole profile from **two** years, the minimum Hannah asked for.

```
window_profiles_backward_from_2024.png         (a) every window over the reference, (b) r per window
window_profiles_panels_backward_from_2024.png  one panel per window
window_profiles_transects.csv                     a row per transect per window (lrr, unc, n_obs, x_domain, diff)
window_profiles_correlation.csv                   a row per window: r_vs_reference, n_transects
```

## r against the reference

| window | years | r |
|---|---|---|
| 2023–2024 | 2 | 0.212 |
| 2022–2024 | 3 | 0.314 |
| 2021–2024 | 4 | 0.367 |
| 2020–2024 | 5 | 0.421 |
| 2019–2024 | 6 | 0.520 |
| 2018–2024 | 7 | 0.548 |
| 2017–2024 | 8 | 0.625 |
| 2016–2024 | 9 | 0.663 |
| 2015–2024 | 10 | 0.695 |
| 2014–2024 | 11 | 0.722 |
| 2013–2024 | 12 | 0.735 |
| 2012–2024 | 13 | 0.750 |
| 2011–2024 | 14 | 0.788 |
| 2010–2024 ← model window | 15 | 0.804 |
| 2009–2024 | 16 | 0.810 |
| 2008–2024 | 17 | 0.837 |
| 2007–2024 | 18 | 0.856 |
| 2006–2024 | 19 | 0.875 |
| 2005–2024 | 20 | 0.897 |
| 2004–2024 | 21 | 0.928 |
| 2003–2024 | 22 | 0.949 |
| 2002–2024 | 23 | 0.972 |
| 2001–2024 | 24 | 0.982 |
| 2000–2024 | 25 | 0.992 |
| 1999–2024 | 26 | 0.995 |
| 1998–2024 | 27 | 0.997 |
| 1997–2024 | 28 | 0.998 |
| 1996–2024 | 29 | 1.000 |

**Read r with care.** The windows are nested, so r reaches 1 at the
reference by construction. What it shows is where it gets there and how
steadily. r measures shape only: a window with every hotspot in the right
place at the wrong magnitude still scores high. Magnitude is scored by the
bias and tolerance tables in `../../2-settling_window/backward_from_2024/b-every_transect/`.

Producer:
`scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_profiles.py`
`--direction backward`.
