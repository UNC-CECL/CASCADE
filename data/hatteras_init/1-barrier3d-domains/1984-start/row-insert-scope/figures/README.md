# row-insert-scope / figures

> Figure file names say `v3`/`v4`/`v5` in the AS-BUILT numbering: `v3` was
> today's `v2` (the base), `v4` today's `v3`, `v5` today's `v4`. Renumbered
> 2026-09-04; map in `../DUNE_TOPO_VERSION_GUIDE.md`. A re-plot writes the
> new names.

Every figure for the 1984 seaward-row insert. Steps are numbered in the order the argument runs (renumbered 2026-09-07:
measurement, footprint, fill, result). All plotters default
their `--out` here — resolved by `hat_topo_version.insert_figures_dir`, not
built by hand. They used to land loose in `dune-topo/`, which turned the version
folder into a figure dump.

## The measurement: where N comes from

| figure | drawn by |
|---|---|
| `HAT_duneline_zoom_GIS83_87.png` | `HAT_plot_duneline_offset.py --zoom 83-87` — labels carry the signed footprint rows (+ added / − removed) since 2026-09-07 |
| `HAT_how_N_determined_GIS84.png` / `_GIS85.png` | `HAT_plot_how_N_is_determined.py` — the feature/date split |
| `HAT_dunelines_on_DEM_GIS85.png` | `HAT_plot_dunelines_on_dem.py` |
| `HAT_dunelines_on_grid_GIS85.png` | `HAT_plot_dunelines_on_grid.py` |

## `2-footprint-1984/` — the symmetric 1984 footprint, both placements (2026-09-07)

The live analysis. One measurement (the paired 1984–1997 dune-line shift, `N = trunc(shift / 10 m)`,
rows added where the 1984 line lay seaward, removed where it lay landward) and TWO placements of the
same rows. Every figure reads `../footprint_1984_by_domain.csv`; captions in `CAPTIONS.md`.

**Root — placement-independent**

| figure | drawn by |
|---|---|
| `HAT_footprint_1984_rows.png` | `HAT_footprint_1984.py` — rows per domain, signed, communities banded |
| `HAT_footprint_1984_shift.png` | same — the paired shift with p10–p90 and the rows kept, island-wide (the canonical view) |
| `HAT_where_inserts_occur_blocks.png` | `HAT_plot_where_inserts_occur.py` — the two NC-12 relocation blocks zoomed |
| `HAT_row_insert_rows.png` | `HAT_report_row_insert_scope.py` — rows per domain, kept beside the report grid |

**`seaward/` — rows at the seaward edge, row 0 moves to the 1984 dune line, the setback becomes `setback_new_m`**

| figure | drawn by |
|---|---|
| `HAT_footprint_1984_grid.png` | `HAT_footprint_1984.py` — the grid in the CURRENT frame (row 0 fixed), elevation classes, added rows blank, removed rows hatched |
| `HAT_row_insert_grid.png` | `HAT_report_row_insert_scope.py` — the grid AS THE MODEL WOULD HOLD IT (dune fixed, interior pushed down) |
| `HAT_footprint_1984_plan.png` | `HAT_footprint_1984.py` — plan view in the dune-line 3-panel layout: hillshade, tickless, full domain boxes |
| `HAT_footprint_1984_setback.png` | same — road setback now and from the new row 0, all road domains |
| `HAT_where_inserts_occur_setback.png` | `HAT_plot_where_inserts_occur.py` — the same at the ten block domains |

**`behind-road/` — the same rows BEHIND NC-12 (advisor's placement): crest-to-road strip as measured, road and row 0 unmoved, today's setback kept**

| figure | drawn by |
|---|---|
| `HAT_row_insert_grid_behindroad.png` | `HAT_report_row_insert_scope.py --anchor road` — the model-frame grid with the blocks landward of the road (`insert_row_behind_road` in the table) |

`1-scope/` was retired the same day: every figure it held was redrawn on this footprint and lives here.
`HAT_row_insert_plan.png` (add-only) and `HAT_where_inserts_occur[_island].png` are gone for good.

## The fill: what the rows are made of

| figure | drawn by |
|---|---|
| `3-fill/HAT_insert_explainer_grid_GIS85.png` | `HAT_plot_insert_explainer_grid.py` — the same mechanics as Barrier3D plan-view grids: survey by year, v2 as the model sees it, the inserted domain, the five fills side by side (2026-09-04) |
| `3-fill/HAT_insert_explainer_GIS85.png` | `HAT_plot_insert_explainer.py` — the mechanics on one line: the survey by year, what the extraction keeps, where the added rows go, what each version writes into them (2026-09-04) |
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
