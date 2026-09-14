# figure_making — figures drawn from finished runs and from the record

Grouped by **what the figure is of**, which is what you know when you go
looking. It was sixteen folders named after the script or the occasion, with no
map and 1212 output files mixed into the code.

```
island/        the island as the model starts it
management/    NC-12 and nourishment: the timeline, the rules, the investigation
shoreline/     shoreline change, observed and modelled
    dsas/        the DSAS rates
    chainage/    the raw CoastSat record, and the gifs of it
model_output/  generic plotting of a finished run's arrays
tools/         a colour picker and an .npy viewer; neither draws a figure
superseded_20260914/
```

## Where the figures go

**Not here.** Products land under `output/` — rule 1 of `ORGANIZATION.md`:

| Figures of | Land in |
|---|---|
| a manuscript or poster subject | `output/figures/<subject>/` |
| the observed record itself | `output/observations/` |
| a comparison across runs | `output/comparisons/` |

Every script here was writing beside itself, by an absolute path into one
machine's home directory. They are anchored on the repository now.

## superseded_20260914/

Three sets retired together, none deleted:

* `old_dsas_scripts/` — an earlier shoreline comparison, replaced by `shoreline/`
* `old_plot_tests/` — plots of the CASCADE package's own test scenarios, not of
  this site at all
* `old_rate_analysis/` — the rate analysis that `shoreline/dsas/` replaced

## What left this tree

The Pea Island scripts moved to `scripts/other_ms/pea_island_ms/`. They are a
**different site**, and `other_ms/` already held three manuscripts. Keeping them
here made the Hatteras figure tree look twice as large as it is.
