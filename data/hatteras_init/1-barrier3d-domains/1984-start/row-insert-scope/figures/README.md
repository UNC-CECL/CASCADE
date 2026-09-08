# row-insert-scope / figures

Every figure for the 1984 row footprint, in the order the argument runs:
**measure** the dune-line shift, turn it into a **footprint** of rows, argue the
**fill**. Each step is a numbered folder; every plotter defaults its `--out`
there through `hat_topo_version.insert_figures_dir`, never by hand. The words
that go under the live figures are in `CAPTIONS.md` (house style: no in-image
titles or footnotes). Reorganised 2026-09-08.

```
1-measurement/        where N comes from
2-footprint-1984/     the rows: placement-independent figures at the root,
    seaward/          ... rows at the seaward edge (row 0 moves to the 1984 line)
    behind-road/      ... the same rows behind NC-12 (the live placement)
3-fill/               what the rows contain: the copy fill
4-result/             the hindcast on the built version (dune-topo/v3) against v2
superseded-layers/    the record of the deleted layers v3-v8 (2026-09-04), regenerable
                      only if a layer is rebuilt
frozen/               pre-re-pick figures nothing on disk can rebuild
```

`4-result/` holds the first result: the hindcast on `dune-topo/v3`, the built
version of this footprint, against the same run on `v2`.

## `1-measurement/` — where N comes from

| figure | drawn by |
|---|---|
| `HAT_duneline_zoom_GIS83_87.png` | `HAT_plot_duneline_offset.py --zoom 83-87` — the two dune lines at true scale on GIS 83–87, each domain labelled with its offset and the signed footprint rows |
| `HAT_how_N_determined_GIS84.png` / `_GIS85.png` | `HAT_plot_how_N_is_determined.py` — the feature/date split of the 1984-line-vs-row-0 number |
| `HAT_dunelines_on_DEM_GIS85.png` | `HAT_plot_dunelines_on_dem.py` |
| `HAT_dunelines_on_grid_GIS85.png` | `HAT_plot_dunelines_on_grid.py` |

The measurement itself is `../duneline-shift/`; the footprint reads its two
per-profile files.

## `2-footprint-1984/` — the rows (2026-09-07)

One measurement (the paired 1984–1997 dune-line shift per profile, median per
domain, `N = trunc(shift / 10 m)`; rows added where the 1984 line lay seaward,
removed where it lay landward) and TWO placements of the same rows. Every
figure reads `../footprint_1984_by_domain.csv`.

**Root — placement-independent**

| figure | drawn by |
|---|---|
| `HAT_footprint_1984_rows.png` | `HAT_footprint_1984.py` — rows per domain, signed, communities banded |
| `HAT_footprint_1984_shift.png` | same — the paired shift with p10–p90 and the rows kept, island-wide (the canonical island-wide view) |
| `HAT_where_inserts_occur_blocks.png` | `HAT_plot_where_inserts_occur.py` — the two NC-12 relocation blocks zoomed |

**`seaward/` — rows at the seaward edge.** Row 0 moves to the 1984 dune line and
the setback becomes `setback_new_m`. The first placement; kept for comparison.

| figure | drawn by |
|---|---|
| `HAT_footprint_1984_grid.png` | `HAT_footprint_1984.py` — the grid in the CURRENT frame (row 0 fixed), elevation classes, added rows blank, removed rows hatched |
| `HAT_row_insert_grid.png` | `HAT_report_row_insert_scope.py` — the grid AS THE MODEL WOULD HOLD IT (dune fixed, interior pushed down) |
| `HAT_footprint_1984_plan.png` | `HAT_footprint_1984.py` — plan view in the dune-line 3-panel layout: hillshade, tickless, full domain boxes, the true-scale band off the 1997 line |
| `HAT_footprint_1984_setback.png` | same — road setback now and from the new row 0, all road domains |
| `HAT_where_inserts_occur_setback.png` | `HAT_plot_where_inserts_occur.py` — the same at the ten block domains |

**`behind-road/` — the same rows BEHIND NC-12** (advisor's placement; the live
one). The crest-to-road strip stays as measured; the block goes in directly
behind the model's two roadway rows AS PLACED under the 1984 setback,
`int(setback_new / 10) + 2` from row 0 (2026-09-08); row 0 and the dune do not move,
the road sits on measured cells at its 1984 distance from the crest, and the added
width is behind it. Domains without a model road (GIS 1–5, 8) anchor on the crest
row instead.

| figure | drawn by |
|---|---|
| `HAT_row_insert_grid_behindroad.png` | `HAT_report_row_insert_scope.py --anchor road` — the model-frame grid with the blocks blank behind the road |
| `HAT_footprint_1984_plan_behindroad.png` | `HAT_footprint_1984.py` — plan view with the true-scale band off NC-12's landward edge |

## `3-fill/` — what the rows contain: the copy fill (2026-09-07/08)

The block is a direct cell-by-cell copy of the N interior rows that follow the
insert point, so it fabricates no value. Audit in `../fill_copy_by_domain.csv`,
report in `../HAT_fill_copy_scope.txt`.

| figure | drawn by |
|---|---|
| `HAT_fill_copy_method_GIS<N>.png` (80, 85, 5, 49) | `HAT_fill_copy_scope.py` — THE METHOD in three stages, model frame: the domain as extracted (dune rows + interior), the N rows inserted blank, the copy fill with the source window and the copy drawn as an arrow |
| `HAT_fill_copy_grid_GIS<N>.png` (80, 85, 5, 49) | same — before and after only, two panels |
| `HAT_fill_copy_grid_island.png` | `HAT_report_row_insert_scope.py --anchor road --fill copy` — the whole island with every block filled by the copy rule and outlined in red; removals hatched |

Since v3 was built (2026-09-08) every after-panel is read from `dune-topo/v3` itself
(asserted equal to the rule) and shows NC-12 where v3's 1984 setback puts the model's
road, on measured cells with the block directly behind it; the pavement rows outlined.

GIS 85 shows the Pea Island case (the road at the dune, so the copied window is
the dune's landward face); GIS 80 the ordinary backbarrier case; GIS 5 the
no-road crest anchor; GIS 49 a three-row block mid-island.

## `4-result/` — the hindcast on v3 against v2 (2026-09-08)

| figure | drawn by |
|---|---|
| `HAT_footprint_v3_vs_v2_road.png` | `HAT_plot_footprint_result.py` — for every road domain, the setback each run started with and every year the model relocated NC-12, v2 (today's setbacks, floored) against v3 (the 1984 setbacks, no floor); the recorded 1989 and 1999 relocations marked. Reads the two runs' saved state: the calibration-tree run and arm `behindroad-copy` |
| `HAT_method_compare_v2_v3_GIS85.png` | `HAT_plot_method_compare.py` — the two methodologies side by side at one domain on the SAME pick set: v2 (measured road −15 m from row 0, floored to 0, road at the dune, no rows) against v3 (road at its 1984 setback of 44 m on measured cells, five rows copied in directly behind it). Only the block and the setback differ |
| `HAT_method_compare_v1_v3_GIS85.png` | same, `--base v1` — against the original extraction, so the 2026-09-02 re-pick (one cell at GIS 85) is part of the difference |

## `superseded-layers/` — the record of the deleted layers

Drawn 2026-09-03/04 for the add-only layers v3–v8 (rows at the seaward edge,
five fill rules), which were deleted 2026-09-07. The scripts survive with a
guard that names what is on disk; the figures can only be regenerated by
rebuilding a layer. File names use the AS-BUILT numbering (`v3` was today's
`v2`, `v5` today's `v4`; map in `../DUNE_TOPO_VERSION_GUIDE.md`).

| figure | drawn by, what |
|---|---|
| `HAT_insert_explainer_grid_GIS85.png` | `HAT_plot_insert_explainer_grid.py` — survey by year, v2 as the model sees it, the inserted domain, the five fills side by side |
| `HAT_insert_explainer_GIS85.png` | `HAT_plot_insert_explainer.py` — the same mechanics on one cross-shore line |
| `HAT_fill_options_grid_GIS85.png` | `HAT_plot_fill_options_grid.py` — the fill candidates as Barrier3D domain views, NC-12 at its road elevation |
| `HAT_fill_options_GIS85.png` | `HAT_plot_fill_options.py` — the same as profiles |
| `HAT_b3d_grid_v3_v5.png` / `_GIS85.png` | `HAT_plot_b3d_grid.py` — base against the measured+floor layer |
| `HAT_insert_three_scales_v3_v5.png` | `HAT_plot_insert_three_scales.py` — the layer at three scales |

## `frozen/` — cannot be regenerated

| figure | why |
|---|---|
| `HAT_b3d_grid_v1_v2.png` | the pre-re-pick `v2` was deleted 2026-09-03. The script survives; the topography does not |
| `HAT_v1_vs_v2.png` | **no drawing script exists** either. Origin unknown |

The only surviving picture of the pre-re-pick (`v1`/`v2`) insert. **Do not try
to rebuild these** — nothing on disk can.

## Retired

| figure | when, why |
|---|---|
| `1-scope/` (the folder) | 2026-09-07 — everything it held was redrawn on the symmetric footprint and lives in `2-footprint-1984/` |
| `HAT_row_insert_plan.png` | 2026-09-07 — add-only plan view; replaced by `2-footprint-1984/seaward/HAT_footprint_1984_plan.png` |
| `HAT_where_inserts_occur.png`, `_island.png` | 2026-09-07 — the four-panel original was split; the island-wide panel duplicated `HAT_footprint_1984_shift.png` |
| `HAT_row_insert_rows.png` | 2026-09-08 — identical to `HAT_footprint_1984_rows.png` |
| five stale layer figures | 2026-09-03 — regenerable duplicates; `../../../archive_purge_20260903.csv` |
