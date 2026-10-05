# scenario_grid

`scenario_grid_by_preset.png` (2026-10-03), written by
`scripts/figure_making/model_output/scenario_grid.py`; the same figure is published as
`output/figures/5-results/scenario_grid.png`.

Modelled shoreline change rate by domain for every management scenario, 1996–2015 and 2010–2026, under each source/sink preset. Periods are rows and presets (zeroBE, edgeBE) columns, in order of increasing correction. Each panel draws one line per scenario against the CoastSat LOWESS target in black; scenarios are a single-hue ramp ordered by management intensity, and each relocation arm is dashed in its non-relocation twin's colour because the two overlap at this scale. The y axis is shared across every panel. Reading across a row shows what the source/sink term does; the spread within a panel shows what management does. The model rate is lrr_m_yr, the OLS slope the runs are scored with, given the target's own 7-domain LOWESS (GIS 1–10 left raw) so both sides are smoothed alike; rates are m/yr, seaward positive. An empty panel is a run not yet made, not a result.
