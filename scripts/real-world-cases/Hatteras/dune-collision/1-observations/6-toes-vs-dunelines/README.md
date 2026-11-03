# 6-toes-vs-dunelines — do the lidar toes sit on the dune Hannah traced?

`HAT_dune_toes_vs_dunelines.py` (2026-11-03). The module's toes are picked from the lidar by
`input_prep/1-barrier3d-domains/3-dune-toes/HAT_dune_toes.py`. This check puts every survey's toe, beach start and
crest back on the map and measures each against the dune line digitised from aerial imagery of the same year
(`data/hatteras_init/2-brie-offset/dunelines/`): 1997, 2004 and 2009 lines on the same year's lidar, the 2023 line
on the 2020 lidar (the last public survey). Distances run along the east-west profiles; positive is the line landward.

Outputs, in the workspace's `dune-collision/1-observations/6-toes-vs-dunelines/`: `dune_toes_vs_dunelines.png`
(offsets alongshore and their spread) and `dune_toes_on_imagery.png` (each survey's toes, crests and beach starts
with the traced line on the aerial image of the same months: the 1997 lidar on the 12 Oct 1997 USGS/USACE photo
with the 1997 line, the 2020 lidar on the Jan-Apr 2020 NOAA image with the 2023 line; GIS 12, 22, 34, 66, 85,
300 m squares); under `supporting/`, one row per survey cell (with map positions) and the summary by survey and
zone. `python HAT_dune_toes_vs_dunelines.py maps` redraws the imagery figure from the saved table.

Results (cells with a kept toe):

| Line vs lidar | Line landward of the toe, median | Middle half | Within 20 m |
|---|---|---|---|
| 1997 vs 1997 | +9 m | +3 to +14 m | 72% |
| 2004 vs 2004 | +17 m | -5 to +28 m | 43% |
| 2009 vs 2009 | +16 m | +8 to +25 m | 50% |
| 2023 vs 2020 | +7 m | +2 to +14 m | 83% |

- The traced line sits on the dune face, between toe and crest (4-12 m seaward of the crest, lidar 3.8-4.5 m
  NAVD88 at the line against toes of 2.3-2.9 m): the light/dark break the 1997 metadata names, likely the
  vegetation edge, not the toe. Same dune, a different point on it.
- 1997 and 2023/2020 agree to about one 10 m lidar cell; the Google Earth lines (2004, 2009) sit further back and
  scatter more.
- GIS ~32-36 and ~63-67 disagree in every year, the line 50-150 m seaward of the lidar toe: the lidar dune window
  likely holds a dune behind the traced one there. Checked with the 1997 lidar contoured on the 1997 photo
  (they register: house outlines and dune ridges fall on the photo's): the 1997 line follows the seaward edge of
  a broad, hummocky dune field, along the 2-3 m contours, and is right; the lidar toe lands 35-90 m inside the field,
  on a low swale between ridges, because the furthest-below-the-chord rule finds the deepest dip on a wide
  multi-ridge profile. There the dry beach the module reads is too wide. By 2020 the shore has moved back to a
  single foredune and the line and toes agree (+1 m at GIS 34, +5 m at GIS 66).
