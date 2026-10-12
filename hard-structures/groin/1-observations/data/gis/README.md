# gis — the GIS layers the groin observations read

Moved here from `HAT-groin-gis-analysis/gis_data/` on 2026-10-08; it was `gis_data/` until 2026-10-12, when `1-observations/` split into `data/` and `analysis/`. Several layers are byte-identical
copies of a canonical file elsewhere in the repo; the groin scripts read the copy here, other code
reads the canonical one. Checked byte for byte on 2026-10-09.

| layer | read by | copy of |
|---|---|---|
| `groins_hatteras.geojson` | `../../analysis/2-gap_across_groins/coastsat_gap_across_groins.py`, `../../analysis/3-shoreline_rates_by_era/HAT_groin_shoreline_analysis.py`, `../../analysis/3-shoreline_rates_by_era/shoreline_plan_view_gif.py`, `scripts/input_prep/4-mgmt-forcings/beach_nourishment.py` | this is the only copy |
| `wet_dry_groin.geojson` | `HAT_groin_shoreline_analysis.py` | `data/hatteras_init/5-scr/1-observations/shoreline_inventory/wet_dry_groin.geojson` |
| `cascade_area.geojson` | `HAT_groin_shoreline_analysis.py` (study-area filter) | `data/hatteras_init/5-scr/1-observations/shoreline_inventory/cascade_area.geojson` |
| `HAT_domains.json` | `HAT_groin_shoreline_analysis.py` | `data/hatteras_init/5-scr/2-transect-frame/transect_domains/HAT_domains.json` (which `gap_across_groins_figures.py` and `shoreline_plan_view_gif.py` read) |
| `transects_100m.geojson` | `HAT_groin_shoreline_analysis.py` | `data/hatteras_init/2-brie-offset/transects/transects_100m.geojson` (which the 2-brie-offset stage reads) |
| `offshore_datum_line.geojson` | nothing | not byte-identical to the one `../../analysis/1-wetdry_photo_positions/HAT_geometric_distance_sanity_check.py` reads, `../../../3-hindcast/1-dipole-1967-2017/inputs/groin_init/island_offset/input/offshore_datum_line.geojson` |
| `CoastSat_shoreline_layer.geojson`, `CoastSat_polygon_layer.geojson` | nothing | CoastSat's global site layers (4,706 sites worldwide, EPSG:4326: median slope and mean trend per site), kept as reference |

`shapefiles/` holds the ArcGIS originals of four of these (`Offshore_Datum_Line_UTM`,
`Transects_100m`, `groins_hatteras`, `wet_dry_groin`). No code reads them, and they are git-ignored,
so they exist on this disk only. They sat loose beside the geojsons until 2026-10-09.
