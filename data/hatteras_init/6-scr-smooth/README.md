# 6-scr-smooth - smoothing the observed rates

What the LOESS smoothing does to the observed shoreline rates, and how the two
rate sources compare once smoothed.

```
HAT_loess_method_comparison_output/  transect-based against domain-averaged,
                                     and the inputs each produces
HAT_loess_dsas_vs_coastsat_output/   CoastSat against DSAS, both smoothed
```

**Both folders are named after the script that made them**, which is unlike the
rest of the tree and is on the list to fix.

The smoothing the model is actually graded against is not here: it is applied
in the runner, at transect resolution over along-coast distance, at a window of
10 domains with the southern 10 spliced in raw. This folder is the argument for
that choice, not the choice itself.
