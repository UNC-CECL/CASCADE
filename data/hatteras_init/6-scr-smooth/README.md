# 6-scr-smooth - smoothing the observed rates

What the LOESS smoothing does to the observed shoreline rates, and how the two
rate sources compare once smoothed.

```
method_comparison/    transect-based against domain-averaged smoothing, and
                      the inputs each produces
    01_transect_based/  02_domain_averaged/  03_cascade_inputs/  04_method_comparison/
dsas_vs_coastsat/     CoastSat against DSAS, both smoothed, on the retired
                      1978-1997 and 1997-2019 windows
```

Until 2026-09-18 these were `HAT_loess_method_comparison_output/` and
`HAT_loess_dsas_vs_coastsat_output/`, named after the scripts that made them.
They are named for what they compare now, and resolved through
`scripts/site_layer/hat_observed_rates.py` (`SMOOTH_METHOD_COMPARISON`,
`SMOOTH_DSAS_VS_COASTSAT`) like the 5-scr folders. Do not type the paths.
`method_comparison/03_cascade_inputs/` is read outside its producer, by
`8-overwash-analysis/HAT_overwash_vs_footprint.py`.

Both folders are gitignored: one run of their producer in
`scripts/input_prep/6-scr-smooth/` regenerates them.

The smoothing the model is actually graded against is not here: it is applied
in the runner, at transect resolution over along-coast distance, at a window of
10 domains with the southern 10 spliced in raw. This folder is the argument for
that choice, not the choice itself.
