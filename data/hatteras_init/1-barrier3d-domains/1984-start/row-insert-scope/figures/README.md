# row-insert-scope / figures

Every figure for the 1984 seaward-row insert. All eight plotters now default
their `--out` here — resolved by `hat_topo_version.insert_figures_dir`, not
built by hand. They used to land loose in `dune-topo/`, which turned the version
folder into a figure dump.

## The scope: where the insert lands

| figure | drawn by |
|---|---|
| `HAT_row_insert_grid.png` | `HAT_report_row_insert_scope.py` — the Barrier3D grid, inserted cells in their own colour |
| `HAT_row_insert_plan.png` | same — the identical thing in plan view, on the DEM |
| `HAT_where_inserts_occur.png` | `HAT_plot_where_inserts_occur.py` — the v4→v5 scope change |

## The measurement: where N comes from

| figure | drawn by |
|---|---|
| `HAT_duneline_zoom_GIS83_87.png` | `HAT_plot_duneline_offset.py --zoom 83-87` |
| `HAT_how_N_determined_GIS84.png` / `_GIS85.png` | `HAT_plot_how_N_is_determined.py` — the feature/date split |
| `HAT_dunelines_on_DEM_GIS85.png` | `HAT_plot_dunelines_on_dem.py` |
| `HAT_dunelines_on_grid_GIS85.png` | `HAT_plot_dunelines_on_grid.py` |

## The fill: what the rows are made of

| figure | drawn by |
|---|---|
| `HAT_fill_options_grid_GIS85.png` | `HAT_plot_fill_options_grid.py` — candidates as Barrier3D domains, NC-12 at its real road elevation |
| `HAT_fill_options_GIS85.png` | `HAT_plot_fill_options.py` — the same as profiles |

## The result: base vs inserted

| figure | drawn by |
|---|---|
| `HAT_b3d_grid_v3_v5.png` / `_GIS85.png` | `HAT_plot_b3d_grid.py` — v3 against v5 |
|  `HAT_insert_three_scales_v3_v5.png` | `HAT_plot_insert_three_scales.py` |

## `frozen/` — cannot be regenerated

| figure | why |
|---|---|
| `HAT_b3d_grid_v1_v2.png` | `v2` was deleted 2026-09-03. The script survives; the topography does not |
| `HAT_v1_vs_v2.png` | **no drawing script exists** either. Origin unknown |

Kept as the only surviving picture of the pre-re-pick (`v1`/`v2`) insert. The
regenerable equivalent is `HAT_b3d_grid_v3_v5.png` on the current pick set.
**Do not try to rebuild these** — nothing on disk can.

## Deleted 2026-09-03

Five figures, all regenerable: two stale `fill_options` duplicates superseded by
renders in this folder, `HAT_fill_options_grid_GIS86.png` (old five-panel
format), and the `HAT_b3d_grid_v3_v4` pair superseded by `v3_v5`. Reasons and
sizes in `../../../archive_purge_20260903.csv`.
