# scenario_grid

`scenario_grid_by_preset.png` (2026-09-11), written by
`scripts/figure_making/model_output/scenario_grid.py`.

Modelled shoreline change rate (m/yr, seaward positive) by GIS domain for
every management scenario, on a 2 x 3 grid: periods as rows (1984-2004
above 2004-2024), source/sink presets as columns (zeroBE, edgeBE, calibBE, in
order of increasing correction). Each panel draws one line per scenario
against the CoastSat LOESS target in black; scenarios are a single-hue ramp
ordered by management intensity, and each relocation arm is dashed in its
non-relocation twin's colour because the two are indistinguishable at this
scale. The y axis is shared across all six panels. Reading across a row shows
what the source/sink term does; the spread within a panel shows what
management does. The model column is `lrr_m_yr`, the OLS slope the runs are
scored with.
