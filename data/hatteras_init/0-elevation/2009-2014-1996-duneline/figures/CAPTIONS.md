# Figure captions

Written by `HAT_plot_duneline_offset.py` (`write_captions`) from the same table the figures draw from. Each heading names the figure's folder under `figures/` (island/, detail/, offset/). The figures carry no in-image titles or footnotes on purpose; use these under them.

## `detail/HAT_duneline_offset_simple.png`

The 1984 (red) and 1997 (blue) dune lines at true scale on four two-domain pairs (3–4, 19–20, 63–64, 79–80), one of them a control where the two agree. Equal aspect, nothing exaggerated in the map plane: each panel is 300 m cross-shore, centred on the 1984 line, by two 500 m domains alongshore. Grey relief is the 1 m gap-filled DEM shaded at 2.2× vertical exaggeration; the shading carries no readable elevation. Numbers beside each domain are its median offset, positive where 1984 lies seaward. Scale bar 50 m = 5 Barrier3D cells.

## `island/HAT_duneline_offset_simple_island.png`

Where the two dune lines are and where they disagree, over the whole island in three equal-aspect panels (south at left). Beside each map, the per-domain median offset as a bar, aligned row for row: red where the 1984 line lies seaward of 1997, blue where it lies landward; dashed guides at ±10 m, one Barrier3D cell. At this scale the two lines coincide within a line width nearly everywhere, which is what the bars are for. The ruler beside each bar is km north of the south end of domain 1, the same origin the ribbon figure uses. Shaded bands mark the pairs shown in the detail figure. Per-domain median offset (1984 minus 1997, seaward positive), from `duneline_offset_by_domain.csv`: island median +1.2 m, range -58.9 to +70.2 m; 52 of 90 domains differ by at least one 10 m Barrier3D cell, 47 of 90 are positive. Communities, village centres, piers and the groin are the project's own positions (`hatteras_site_config.HATTERAS_ANNOTATIONS`); structures are drawn seaward off the 1984 line at a fixed 380 m, a mark rather than a surveyed extent.

## `island/HAT_duneline_offset_simple_island_mean.png`

As the previous figure, with the per-domain MEAN offset on the bars in place of the median (island mean +3.6 m, range -60.3 to +67.8 m). The two differ where a domain's 1 m samples are skewed: a short stretch of large offset inside an otherwise quiet domain moves the mean and not the median.

## `island/HAT_duneline_offset_lines_island.png`

The two dune lines over the whole island as maps only, in 10-domain (~5 km) panels reading south (left) to north (right). Each panel is cropped at equal aspect to the envelope of the two lines plus 150 m either side, so the separation is visible on the map itself without exaggeration. Domain numbers every fifth domain on the landward edge; shaded bands mark the pairs in the detail figure. Communities, village centres, piers and the groin are the project's own positions (`hatteras_site_config.HATTERAS_ANNOTATIONS`); structures are drawn seaward off the 1984 line at a fixed 380 m, a mark rather than a surveyed extent.

## `island/HAT_duneline_offset_lines_island_3panel.png`

As the previous figure, in three 30-domain panels matching the map-and-bar figure's layout. At 15 km per panel a 50 m offset is about one line width, so the separation reads only where it is large. Communities, village centres, piers and the groin are the project's own positions (`hatteras_site_config.HATTERAS_ANNOTATIONS`); structures are drawn seaward off the 1984 line at a fixed 380 m, a mark rather than a surveyed extent.

## `offset/HAT_duneline_offset_ribbon.png`

(a) The 1984 and 1997 dune lines along the island at 1 m alongshore sampling, each drawn relative to their common 2 km boxcar-smoothed midline so the island's curvature drops out; the band between them is filled red where 1984 lies seaward and blue where it lies landward. (b) The difference, 1984 minus 1997, filled by sign; dashed guides at ±10 m, one Barrier3D cell. Grey bands mark the reaches shown in the DEM detail figure (17-21 (the quietest reach on the island); 62-68 (1984 line landward of 1997); 78-85 (1984 line seaward of 1997)). Domain numbers along the top; the distance axis starts at the south end of domain 1. Per-domain median offset (1984 minus 1997, seaward positive), from `duneline_offset_by_domain.csv`: island median +1.2 m, range -58.9 to +70.2 m; 52 of 90 domains differ by at least one 10 m Barrier3D cell, 47 of 90 are positive.

## `detail/HAT_duneline_offset_zooms.png`

The two dune lines at true scale on three reaches of five to eight domains: 17-21 (the quietest reach on the island); 62-68 (1984 line landward of 1997); 78-85 (1984 line seaward of 1997). Equal aspect; each panel is cropped to 300 m either side of the local 1984 line rather than the full 2000 m domain box. Grey relief is the 1 m gap-filled DEM shaded at 2.2× vertical exaggeration and carries no readable elevation; the domain boxes are outlined; the number beside each domain is its median offset, positive where 1984 lies seaward. Scale bar 50 m = 5 Barrier3D cells.

## `detail/HAT_duneline_offset_zoom_83_87.png`

As the previous figure, for domains 83–87: GIS 85 and its neighbours, the largest sustained positive run on the island. Where the 1984 footprint table is on disk the label also gives the number of Barrier3D rows the 1984 start would add (+) or remove (−) there, trunc(paired shift / 10 m): a row only once a full cell of change is measured, in either direction.

## `offset/HAT_duneline_offset_bydomain.png`

Median cross-shore offset between the 1984 and 1997 dune lines in each of the 90 domains, 1984 minus 1997 with seaward positive, with the interquartile range of the 1 m samples within the domain. Red markers: 1984 seaward; blue: 1984 landward. Dashed guides at ±10 m, one Barrier3D cell. Grey bands are the communities. Per-domain median offset (1984 minus 1997, seaward positive), from `duneline_offset_by_domain.csv`: island median +1.2 m, range -58.9 to +70.2 m; 52 of 90 domains differ by at least one 10 m Barrier3D cell, 47 of 90 are positive.

