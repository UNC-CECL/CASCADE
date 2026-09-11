# `v3/figures` — the figures of a built version

Written 2026-09-10 19:17 by `HAT_plot_version_figures.py`
(`scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/5-build/`).

`v3` was **built** from `v2`, not extracted, so it does not get the
extractor's `qc/` and `gis_vs_processed/` sets (those compare the raw DEM
profile with the extraction and need the picks and the straightening frame).
It gets the counterparts a built version can answer, every number read from
its own arrays, its `RoadSetback_1984_dunestart.csv` and
`HAT_footprint_audit.csv`; nothing is re-measured.

On the figures `v2` is called "as extracted (1996 surface)" and `v3`
"1984 reconstruction"; the images carry no titles or statistics, so the captions
below do. Every PNG has a PDF beside it except the raster-only grid panels.

| figure | what |
|---|---|
| `grid/domain_NNN_grid_v3.png` (90) | one per domain, a single-column figure: (a) `v2` beside (b) `v3` as the model holds them — the two dune rows on top, drawn at berm + dune height, every interior row down the page, elevation classes (m MHW), NC-12's two rows at each version's setback, and the footprint: the inserted block outlined in red (dashed in the source: the rows it copies) or the removed rows hatched blue in the source and the seam marked in the version. The panel titles give the row count and the change. Unchanged domains have two identical panels. |
| `../HAT_dune_topo_summary_v3.png` | every domain on one page, v2 (grey) against v3 (purple), the villages banded: (a) interior rows added (red) or removed (blue) — 29 domains gain 73 rows, 23 lose 47, 38 are unchanged; (b) interior rows per domain; (c) the NC-12 setback, as measured on the 1996 surface and as the model receives it for 1984; (d) mean interior elevation over land cells, with the dune crest (berm + dune height, green), which the build does not change. Domain 1 is at Cape Point, 90 at Pea Island. The counterpart of the extractor's summary page. |
| `../HAT_dune_topo_island_planview_v3_<year>_{trimmed,padded}.png` | the island in plan view at the period's dune offsets (dune row plus interior; `padded` pads every domain to the extractor's cross-shore length, `trimmed` keeps each domain's own), NC-12 drawn where the **model** places it (the version's setback) rather than from the GIS mask, whose frame a built interior no longer shares. Villages as light bands over the water. |

Files:

- `HAT_dune_topo_summary_v3.png`
- `HAT_dune_topo_island_planview_v3_1984_trimmed.png`
- `HAT_dune_topo_island_planview_v3_1984_padded.png`
- `figures/grid/` — 90 domain panels
