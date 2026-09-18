# figure_making — figures drawn from finished runs and from the record

Grouped by **what the figure is of**, which is what you know when you go
looking. It was sixteen folders named after the script or the occasion, with no
map and 1212 output files mixed into the code.

```
island/        the island as the model starts it; HAT_study_area_figures.py draws
               the generic site figures (study area, domain framework, one domain)
management/    NC-12 and nourishment: the timeline, the rules, the investigation
               (diagnose_road_drowning.py is STALE: it reads a 2009_v2
               topography layout and a HAT_1984_2004_base run that no longer
               exist, and its domains are the t=0 drowning the gap-filled DEM
               removed. Kept until it is either repointed or retired)
shoreline/     shoreline change, observed and modelled
    dsas/        the DSAS rates
    chainage/    the raw CoastSat record, and the gifs of it
model_output/  a finished run's arrays, and the cross-run figures: the
               hindcast result, the scenario grid, the GIS 11 drowning
               figure, the plan-view gif and the run re-renderer (from
               scripts/hatteras_ms/figures/, 2026-09-18)
tools/         HAT_figure_index.py (writes output/figures/README.md),
               HAT_clip_nc_coast.py (rebuilds the NC coast map layer), a
               colour picker and an .npy viewer
STYLE.md       the house style in words, written by write_style_sheet() in
               scripts/site_layer/hat_figure_style.py
superseded_20260914/
model_output/superseded_20260918/   the two 1978-1997 gif scripts; see WHY.md
```

## Where the figures go

**Not here.** Products land under `output/` — rule 1 of `ORGANIZATION.md`:

| Figures of | Land in |
|---|---|
| the site: reach, domains, one domain | `output/figures/site/` |
| the forcing record | `output/figures/forcing/` |
| NC-12 and nourishment | `output/figures/management/` (the investigation plots in `management/investigation/`) |
| shoreline change | `output/figures/shoreline/` |
| the island as the model starts it | `output/figures/initialization/` |
| the same figures for a projector | `output/figures/talk/<subject>/` |
| the house style, drawn | `output/figures/style/` |
| the observed record itself | `output/observations/` |
| a comparison across runs | `output/comparisons/` |

A cross-run figure that is finished for the manuscript is ALSO written to its
subject folder, in the house style, by the same run of the same script, so the
two copies cannot drift: `HAT_scenario_grid.py` writes
`output/figures/shoreline/scenario_grid.png` beside its `comparisons/` copy.
`HAT_hindcast_final_figure_loess.py` writes `shoreline/hindcast_<preset>.png`,
and `HAT_gis11_relocation_drown_figure.py` writes
`management/gis11_relocation_drown.png`. Both were redrawn for the house-style
column on 2026-09-18, and each has a `PUBLISH` switch to hold a figure back
while its layout is broken.

Resolve these through `scripts/site_layer/hat_figure_style.py`:
`figure_dir(subject, ...)`, `COMPARISONS_ROOT`, `OBSERVATIONS_OUT`.
`figure_dir()` raises on a subject not in the table above, so a new top-level
folder cannot appear by typo -- which is how
`output/figures/management_investigation/` came to sit beside `management/`
until 2026-09-18.

`output/figures/README.md` is a generated index of all of them: figure, what it
shows, and the script that draws it. Re-run
`scripts/figure_making/tools/HAT_figure_index.py` after adding a figure. A
figure with a dash in the "shows" column has no caption recorded, and one with
a dash under "drawn by" is an orphan nothing in this tree can reproduce.

Every script here draws in the house style: `apply_style()` from
`scripts/site_layer/hat_figure_style.py`, the vintage red/blue pair for two periods, and
nothing on the canvas that belongs in a caption. The five legacy scripts were
converted on 2026-09-17; before that each had its own typeface, font sizes and
palette.

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
