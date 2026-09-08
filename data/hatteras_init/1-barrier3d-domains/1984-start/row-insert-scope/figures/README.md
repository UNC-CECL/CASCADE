# row-insert-scope / figures

> Figure file names say `v3`/`v4`/`v5` in the AS-BUILT numbering: `v3` was
> today's `v2` (the base), `v4` today's `v3`, `v5` today's `v4`. Renumbered
> 2026-09-04; map in `../DUNE_TOPO_VERSION_GUIDE.md`. A re-plot writes the
> new names.

Every figure for the 1984 seaward-row insert. All eight plotters now default
their `--out` here — resolved by `hat_topo_version.insert_figures_dir`, not
built by hand. They used to land loose in `dune-topo/`, which turned the version
folder into a figure dump.

## The scope: where the insert lands

| figure | drawn by |
|---|---|
| `1-scope/HAT_row_insert_grid.png` | `HAT_report_row_insert_scope.py` — the grid AS THE MODEL WOULD HOLD IT: dune rows on top, added rows blank (red), removed rows hatched (blue), NC-12 marked. Symmetric since 2026-09-07; reads the footprint table |
| `1-scope/HAT_row_insert_rows.png` | same — rows per domain, signed, split out of the grid 2026-09-07 |
| ~~`HAT_row_insert_plan.png`~~ | retired 2026-09-07 (it was add-only). The plan view is `1-scope/HAT_footprint_1984_plan.png` below |
| `1-scope/HAT_where_inserts_occur_blocks.png` | `HAT_plot_where_inserts_occur.py` — the two relocation blocks zoomed: paired shift with spread and the rows kept. Rewritten 2026-09-07 on the footprint table (it used to compare two deleted layers) and split into three figures |
| `1-scope/HAT_where_inserts_occur_setback.png` | same — the NC-12 setback at the block domains, now and from the new row 0 |
| ~~`HAT_where_inserts_occur_island.png`~~ | retired 2026-09-07: it duplicated `HAT_footprint_1984_shift.png`, the canonical island-wide view |
| `1-scope/HAT_footprint_1984_grid.png` | `HAT_footprint_1984.py` (2026-09-07) — the grid with rows ADDED (blank, red) and REMOVED (hatched, blue), new row 0 and NC-12 marked; the symmetric, 10 m-rule footprint |
| `1-scope/HAT_footprint_1984_rows.png` | same — rows per domain, signed, on its own (split out of the grid figure 2026-09-07) |
| `1-scope/HAT_footprint_1984_plan.png` | same — plan view in the layout of the dune-line 3-panel figure: hillshade, tickless, full domain boxes, both dune lines and NC-12 1984 |
| `1-scope/HAT_footprint_1984_shift.png` | same — paired shift with p10–p90 and the rows kept, communities and relocation blocks marked (was panel (a) of `_N.png`, split 2026-09-07) |
| `1-scope/HAT_footprint_1984_setback.png` | same — road setback now and from the new row 0, per road domain (was panel (b) of `_N.png`) |

Captions for the three footprint figures are in `CAPTIONS.md` (written by the script); they carry no in-image titles.

## The measurement: where N comes from

| figure | drawn by |
|---|---|
| `HAT_duneline_zoom_GIS83_87.png` | `HAT_plot_duneline_offset.py --zoom 83-87` — labels carry the signed footprint rows (+ added / − removed) since 2026-09-07 |
| `HAT_how_N_determined_GIS84.png` / `_GIS85.png` | `HAT_plot_how_N_is_determined.py` — the feature/date split |
| `HAT_dunelines_on_DEM_GIS85.png` | `HAT_plot_dunelines_on_dem.py` |
| `HAT_dunelines_on_grid_GIS85.png` | `HAT_plot_dunelines_on_grid.py` |

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
