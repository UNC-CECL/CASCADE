# `v3/figures` — the figures of a built version

Written 2026-09-09 15:07 by `HAT_plot_version_figures.py`
(`scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/5-build/`).

`v3` was **built** from `v2`, not extracted, so it does not get the
extractor's `qc/` and `gis_vs_processed/` sets (those compare the raw DEM
profile with the extraction and need the picks and the straightening frame).
It gets the counterparts a built version can answer, every number read from
its own arrays, its `RoadSetback_1984_dunestart.csv` and
`HAT_footprint_audit.csv`; nothing is re-measured.

| figure | what |
|---|---|
| `grid/domain_NNN_grid_v3.png` (90) | one per domain: `v2` beside `v3` as the model holds them — the two dune rows on top (berm + dune height), every interior row down the page, elevation classes (m MHW), NC-12's two rows at each version's setback, and the footprint: the inserted block outlined in red (add; dashed in the source: the rows it copies) or the removed rows hatched blue in the source and the seam marked in the version (remove). Unchanged domains have two identical panels. |
| `../HAT_dune_topo_summary_v3.png` | every domain on one page: rows added or removed, interior rows, the road setback the model receives, mean interior elevation and the dune crest, `v2` against `v3`, communities banded. The counterpart of the extractor's summary page. |
| `../HAT_dune_topo_island_planview_v3_<year>_{trimmed,padded}.png` | the island in plan view at the period's dune offsets, in the extractor's poster style, NC-12 drawn where the **model** places it (the version's setback) rather than from the GIS mask, whose frame a built interior no longer shares. |

Files:

- `HAT_dune_topo_summary_v3.png`
- `HAT_dune_topo_island_planview_v3_1984_trimmed.png`
- `HAT_dune_topo_island_planview_v3_1984_padded.png`
- `figures/grid/` — 90 domain panels
