# 2026-09-29 — output/figures/ before the numbered layout

The whole `output/figures/` tree as it stood on 2026-09-29, moved here intact
when Hannah asked for the figures to be easier to navigate and brought up to
date. `output/figures/` is gitignored, so this copy is the only record of it.

**What replaced it:** a numbered layout in the paper's order, with every figure
regenerated from its script on 2026-09-29:

    1-site/  2-observations/  3-model-inputs/  4-model-mechanics/  5-results/  talk/  style/

`output/figures/README.md` is the map, and `png_list.txt` beside this file
lists the 182 PNGs that were here. The old subject folders map across like this:

| old folder | new place |
|---|---|
| `site/` | `1-site/` |
| `shoreline/` observed (rates, dune lines, mean shoreline) | `2-observations/{shoreline,duneline,shoreline_vs_duneline,mean_shoreline}/` |
| `shoreline/hindcast_*`, `scenario_grid` | `5-results/` |
| `pipeline/<step>/`, `initialization/`, `forcing/`, `management/` | `3-model-inputs/<step>/` |
| `management/gis11_relocation_drown.png` | `5-results/` |
| `model/` | `4-model-mechanics/{barrier3d,brie,cascade,storm_routing}/` |

**Stale here, do not use:** the 09-17/09-18 site, forcing, management,
initialization and talk figures predate option-A waves (09-27), the Dmaxel
ceilings and v3_trim24 storms (09-28); the `mean_shoreline/*_1995_1997_*` and
`*_2009_2011_*` images use the calendar window that the DEM-centred window
replaced on 09-29.
