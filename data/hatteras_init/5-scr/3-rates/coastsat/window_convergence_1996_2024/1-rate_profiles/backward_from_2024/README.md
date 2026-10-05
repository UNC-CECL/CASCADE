# 1-rate_profiles/backward_from_2024 — what does each window's alongshore profile look like?

Every window of the nested family (2023–2024 … 1996–2024, the end pinned at
2024), drawn as a shoreline change rate profile along all 906 CoastSat
transects over the 1996–2024 reference. Windows from **two** years, the minimum
Hannah asked for. Built 2026-10-01.

```
window_profiles_overlay_backward_from_2024.png  every window over the reference (±10 m/yr)
window_profiles_panels_backward_from_2024.png   one panel per window, the reference last (±5 m/yr)
window_profiles_transects.csv                     a row per transect per window (lrr, unc, n_obs, x_domain, diff)
```

How close each profile is to the reference, as r, bias and RMSE with 95%
intervals: `../../2-r_bias_rmse/`. How many years each place needs:
`../../3-settling_window/backward_from_2024/`.

Producer:
`scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_profiles.py`
`--direction backward`.
