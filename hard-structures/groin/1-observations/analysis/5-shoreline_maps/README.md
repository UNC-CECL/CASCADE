# shoreline_maps — where the shoreline was each year, as true-scale maps

**The map view.** The CoastSat shoreline around the Buxton groins and along the whole coast
(GIS 1-65), one frame per year from 1984 to 2025, drawn as real maps rather than rates or a single
gap. It measures nothing new: each line is the median CoastSat position per transect that year, from
the chainage table in `../4-shoreline_rates_by_era/tables/`. It is for seeing where, and when, the
shoreline moved; the numbers are in `../2-coastsat_gap/` (the gap each year) and
`../4-shoreline_rates_by_era/` (rates by era).

**Who reads it.** `../../../4-hindcast/2-blocking-1996-2025/2026-10-09-corrected-groin/plan_view_figure.py`
imports `shoreline_plan_view_gif.py` for its map layers, unchanged.

Until 2026-10-13 these were part of `../4-shoreline_rates_by_era/`: the script at its top level,
this folder's figure 1 as its figure 5 (`shoreline_rates_by_era_5_groin_window_four_years`), and the
GIFs in its `figures/` (loose until 2026-10-09, then `v1_plain/` and `v2_satellite/`).

```
shoreline_plan_view_gif.py   figure 1 and the map GIFs -> figures/; reads ../4-shoreline_rates_by_era/tables/
figures/                     figure 1; its PDF and CAPTIONS.md in supporting/
    v1_plain/                the map GIFs on a plain sand-and-sea basemap: whole reach, groin close-up,
                             shoreline_map_tiles/ (the coast in 1.5 km map GIFs)
    v2_satellite/            the same GIFs over satellite imagery
```

The GIFs are git-ignored (tens of MB); rerun to make them.

    python shoreline_plan_view_gif.py   # figure 1 and the map GIFs (--only static|whole|closeup|tiles,
                                        # --basemap plain|satellite); the tiles take a while

| # | figure | question |
|---|---|---|
| 1 | `shoreline_maps_1_groin_window_four_years` | how far, and when, did the shoreline around the groins move? The 1.5 km groin window as true-scale maps at 1985, 1998, 2010 and 2025, with its change since 1984-1986 and the fills marked (the static version of the close-up GIF, for papers and slides) |
| GIF | `v1_plain/shoreline_map_whole_reach` | where was the shoreline each year, 1984-2025? GIS 1-65 as one map |
| GIF | `v1_plain/shoreline_map_groin_closeup` | the same around the Buxton groins: a 1.5 km map (GIS 5-8), where it is on Hatteras, and how far the window's shoreline has moved since 1984-1986 |
| GIF | `v1_plain/shoreline_map_tiles/shoreline_map_tile_<nn>_GIS<a>-<b>` | the same 1.5 km view walked along the whole coast, GIS 1-65, numbered south to north and named by the domains each covers |

**The shoreline GIFs (2026-10-08).** True maps: UTM 18N, north up, drawn 1:1 with no stretching, a scale bar and north arrow on each. Every CoastSat observation sits at its real position, the transect's origin plus its chainage along that transect's own shore-normal direction from the CoastSat transect layer (all 1,255 lines start within 5 m of the chainage origins, checked in the script). One frame per year, 1984-2025: the median position per transect that year (at least 3 images) in black, the earlier years as thin lines coloured light to dark, the 1984-1986 shoreline dashed, the four groins as their real lines. The island is filled in sand up to that year's shoreline: the island outline (`site_layer.hat_map_layers.ISLAND_OUTLINE`, widened 250 m with square corners so it reaches every shoreline on these maps; they sit at most about 210 m seaward of it) minus everything seaward of the year's line, where a transect without images that year takes its 1984-1986 position. At true scale the year-to-year movement (tens of metres) is small next to kilometres of coast, so it is subtle by design. These replace three earlier versions tried the same day (the 100 m grid, a rotated frame, a change-from-1984 view), which either stretched the vertical scale or did not look like the coast, and the analysis script's own GIFs (`LEGACY_FIGURES = False` in `../4-shoreline_rates_by_era/HAT_groin_shoreline_analysis.py`). The 1.5 km views (the close-up and the tiles, centred every 1.4 km along the 1984-1986 shoreline, so neighbours overlap by 100 m) also draw the GIS domain bands (white lines, numbers on the land), a locator, and the window's median change since 1984-1986 per year (transects whose 1984-1986 position lies in the window), building up frame by frame. They replaced an eight-panel 4 km sections GIF, in which the change was too small to see at true scale.

**What else is on the maps.** NC-12 (the current line, `site_layer.hat_topo_version.ROAD_LINE_CURRENT`; it starts about 1 km north of the groins, so the groin window shows no road), village names at the middle of their GIS domains, Avon Pier, and the beach fills from the model's nourishment list (`data/hatteras_init/4-mgmt-forcing/nourishment/2-model-input/nourishment_projects.csv`: Buxton 2017 and 2022, GIS 6-16; Avon 2022, GIS 21-28) shaded on each window's change series where their footprint reaches it, with the fill's name and volume on the map in its year. The 1971 and 1973 Buxton fills predate CoastSat.

**v2, over satellite imagery (`--basemap satellite`, `figures/v2_satellite/`).** The same maps with Esri World Imagery behind them instead of the sand and ocean fill: the shorelines in white with a dark halo, earlier years from yellow to red, NC-12 in road yellow. The imagery is a single recent mosaic, so it shows today's beach under every year's line; it is context, not data. The 2025 line falls on the imaged beach, an independent check that the positions are placed correctly.
