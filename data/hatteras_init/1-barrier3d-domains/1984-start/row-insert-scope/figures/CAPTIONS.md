# Figure captions

Written by `HAT_footprint_1984.py` (`write_captions`) and by hand for the figures of the other two scripts. All of these live under `2-footprint-1984/` (root, `seaward/` or `behind-road/` — see README.md). The figures carry no in-image titles or footnotes on purpose; use these under them.

## `HAT_row_insert_grid.png`

(a–c) The 90 Barrier3D domains as the model would hold them, 50 alongshore cells each, every cross-shore row down the page from the two dune rows (right axis in metres; interiors run to 189 rows, white is off the array, cells are not square). Existing cells are shaded by elevation class (m above MHW), so the dune ridge, the backbarrier flat and the sound-side marsh are distinguishable. Where the 1984 dune line lay seaward of the 1997 line, N blank rows (red) are added between the dune and the existing interior, which is pushed down the page; where it lay landward, the existing rows 0 to |N|−1 are hatched (blue) and the interior the model would hold starts below them; the signed count is printed above each changed domain. The dark bar is NC-12 at its measured 1984 position, 20 m wide; it moves down with the interior where rows are added and stays put where rows are removed. Communities and villages along the bottom of each panel, from the site configuration. 29 domains gain 73 rows, 23 lose 47, 38 are unchanged.

## `HAT_row_insert_grid_behindroad.png`

(a–c) As `HAT_row_insert_grid.png`, with the same rows placed directly BEHIND THE ROADWAY ROWS instead of at the seaward edge. The roadway in the model is two straight rows at one setback per domain, so the block goes in at int(setback / 10) + 2 rows from interior row 0 and the strip from the dune crest through the road is kept exactly as measured; the added rows (red, blank) push only the backbarrier down the page, and removals (blue hatching) take backbarrier rows at the same index. The dark bar is NC-12 where the model holds it, at the setback the model receives today (floored at 0, so at GIS 85 and 86 the road is rows 0–1 and the block starts at row 2). Row 0 and the road do not move, so the model keeps today's setback. Domains without a model road (GIS 1–5, 8) fall back to the seaward placement. The insert row per domain is `insert_row_behind_road` in `footprint_1984_by_domain.csv`. N is identical to the seaward placement: 29 domains gain 73 rows, 23 lose 47, 38 are unchanged.

## `HAT_row_insert_rows.png`

Rows per domain under the 10 m rule, positive where rows are added and negative where existing rows are removed, with the communities banded along the axis. Identical in content to `HAT_footprint_1984_rows.png`; kept beside the report grid so the pair reads together.

## `HAT_where_inserts_occur_blocks.png`

(a) The inter-village NC-12 relocation block, GIS 9–14, and (b) the Pea Island block, GIS 84–87: per domain, the median of the 50 paired per-profile differences between the 1997 and 1984 dune-line crossings (points, p10–p90), positive where the 1984 line lies seaward, and the rows the 10 m rule keeps of it (bars, N × 10 m; red added, blue removed, grey unchanged), labelled with the signed row count. Dashed guides at ±1 cell. A green triangle marks N if the 13-year measurement were scaled to the 12-year 1984–1996 interval; recorded, not applied.

## `HAT_where_inserts_occur_setback.png`

The NC-12 setback at the ten relocation-block domains, in metres landward of interior row 0: as the model receives it today (grey; the v2 measurement floored at 0, so GIS 85 and 86 stand at 0) and from the new row 0 (coloured by what happens to the domain, with the p10–p90 over the 50 profiles), values printed above each bar. The new value is (road − row 0) + shift per profile, then the median, unrounded.

## `HAT_footprint_1984_grid.png`

(a–c) The 90 Barrier3D domains as the model indexes them, 50 alongshore cells each, every cross-shore row down the page from the CURRENT interior row 0 at 0 on every domain (topography v2; interiors run to 189 rows, white is off the array, cells are not square; right axis in metres). Existing cells are shaded by elevation class (m above MHW), so the dune ridge, the backbarrier flat and the sound-side marsh are distinguishable. Red rows above 0 are the rows that would be added to bring row 0 to the 1984 dune line, blank because no fill has been chosen; blue hatching marks existing rows 0 to |N|−1 that would be removed where the island has prograded since 1984; the signed count is printed above each changed domain. The black tick is where interior row 0 ends up; the dark bar is NC-12 at its measured 1984 position (seaward edge, 20 m wide), which does not move, so its distance to the tick is the new setback. Communities and villages along the bottom of each panel, from the site configuration. 29 domains gain 73 rows and 23 lose 47; 38 are unchanged. N = trunc(median paired shift / 10 m), so a row appears only once a full cell of change is measured; the largest are +7 at GIS 80 and -5 at GIS 63.

## `HAT_footprint_1984_rows.png`

Rows per domain under the 10 m rule, positive where rows are added (the 1984 dune line lay seaward of the 1997 line) and negative where existing rows are removed, with the communities banded along the axis. 29 domains gain 73 rows and 23 lose 47; 38 are unchanged. N = trunc(median paired shift / 10 m), so a row appears only once a full cell of change is measured; the largest are +7 at GIS 80 and -5 at GIS 63.

## `HAT_footprint_1984_plan.png`

(a–c) The same footprint in plan view, in three equal-aspect panels of thirty domains, south at left, each showing its domain boxes (2000 × 500 m) whole. Grey relief is the 2009-2014-1996 DEM hillshaded at 10 m and carries no readable elevation. The 1984 (red) and 1997 (blue) dune lines, NC-12 in 1984 (dashed), and each domain box shaded by N, red where rows are added and blue where they are removed, labelled on the landward edge. The solid band is the same N at true scale, the 1997 dune line offset seaward (add) or landward (remove) by N × 10 m, so one row is a 10 m sliver that needs a zoom; the band is anchored on the 1997 line as the map proxy for the existing array's seaward edge, which itself sits about 19 m seaward of interior row 0. The gap between the band's outer edge and the 1984 line is the truncation to whole cells. Communities as a bracket in the ocean margin, the pier and groin as seaward marks; scale bar 500 m = 50 cells. 29 domains gain 73 rows and 23 lose 47; 38 are unchanged. N = trunc(median paired shift / 10 m), so a row appears only once a full cell of change is measured; the largest are +7 at GIS 80 and -5 at GIS 63.

## `HAT_footprint_1984_shift.png`

Per domain, the median of the 50 paired per-profile differences between the 1997 and 1984 dune-line crossings (points, p10–p90 bars), positive where the 1984 line lies seaward, and the part of it the 10 m rule keeps as whole rows (filled bars, N × 10 m; red added, blue removed). Dashed guides at ±1 cell; points inside them become no rows, and the distance from a point to its bar is the truncation residual, always toward less change. Communities banded along the axis; the two NC-12 relocation blocks outlined. 29 domains gain 73 rows and 23 lose 47; 38 are unchanged. N = trunc(median paired shift / 10 m), so a row appears only once a full cell of change is measured; the largest are +7 at GIS 80 and -5 at GIS 63.

## `HAT_footprint_1984_setback.png`

For the 82 road domains, the NC-12 setback the model receives today (open circles, the v2 measurement floored at 0) and the setback from the new row 0 (filled, coloured by what happens to the domain, with the p10–p90 over the profiles), both in metres landward of interior row 0 on a symmetric-log axis. The new value is (road − row 0) + shift per profile, then the median, unrounded. GIS 85 and 86, floored at 0 today, become 44 and 12 m; GIS 16 falls to 1 m after losing four rows. No new setback is negative.
