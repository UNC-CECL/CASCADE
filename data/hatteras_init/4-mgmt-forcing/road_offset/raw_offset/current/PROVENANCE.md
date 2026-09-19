# nc12_current.geojson — NC-12 as it stands today

Extracted 2026-09-18 for the dune-line position figures
(`5-scr/4-comparisons/duneline_positions/`), which need a road line for the
2023 dune line. The 1978 and 2008 exports stand in for 1996 and 2010.

- **Source:** NCDOT route inventory, `D:\Hatteras_GIS\Roads\NCRoutes.gdb`,
  layer `NCRoutes`, route `30000012028` (NC-12, Dare County, milepost 0 to
  84.22). The layer's native CRS is EPSG:2264 (NC State Plane, US feet); the
  Z and M values were dropped.
- **Clip:** to the `nc12_2008` line's extent, ±1.5 km cross-shore and ±50 m
  alongshore. That leaves 60.9 km, against 58.4 km for the 2008 line.
- **Against 2008:** 86% of the line lies within 5 m of the 2008 alignment,
  and about 6 km sits more than 30 m from it (northing 3,902,700 to
  3,958,300, UTM 18N). Those are the post-2008 realignments, including the
  Pea Island bridge (2018) and the Rodanthe "Jug Handle" bridge (2022).
- **The date is not known.** The file carries no survey date; it is the
  alignment as of the database copy on D:, so it stands for about 2023.
- **Not a model input.** `hat_topo_version.ROAD_LINE_FOR_YEAR` does not list
  it, and no hindcast reads it.
