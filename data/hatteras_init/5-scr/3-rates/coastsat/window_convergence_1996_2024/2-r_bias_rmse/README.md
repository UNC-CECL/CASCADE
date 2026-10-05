# 2-r_bias_rmse — how close is each window's profile to 1996–2024?

Every nested window's alongshore shoreline change rate profile (all 906
CoastSat transects, from `../1-rate_profiles/<direction>/window_profiles_transects.csv`)
scored against the 1996–2024 profile three ways, by window length:

- **r**, alongshore Pearson correlation: are the hotspots in the same places?
  Shape only.
- **bias**, the mean of window rate minus 1996–2024 rate: is the overall level
  right? Negative = more erosional than the long-term record.
- **RMSE**, the root-mean-square of the same difference: the typical size of
  the miss at one transect, sign ignored (RMSE² = bias² + scatter²).

95% intervals from 1000 bootstrap resamples of the 90 domains.

```
window_profiles_r.png              r, both directions; r = 0.5 / 0.75 / 0.9 marked
window_profiles_bias_rmse.png      (a) bias, (b) RMSE, both directions
window_profiles_r_bias_rmse.csv    a row per direction per window length:
                                   r, bias_m_yr, rmse_m_yr, each with _lo95/_hi95
forward_from_1996/                 the two figures for one direction, with the
backward_from_2024/                calendar window on the top axis
```

At the 15-year model window: r 0.69 / 0.80, bias −0.49 / +1.00 m/yr, RMSE
1.51 / 1.55 m/yr (1996–2010 / 2010–2024). The windows are nested in the
reference, so every score reaches its perfect value at 29 years by
construction. Captions, with the intervals, in `supporting/CAPTIONS.md`.

Producer (reads the fits, no refitting; run `coastsat_window_profiles.py`
first after any change to them):
`scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_r_bias_rmse.py`
