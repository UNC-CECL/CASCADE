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
300 m squares; toes green where within 20 m of the line, two lidar cells, red where further, each panel titled
with the share that agrees); under `supporting/`, one row per survey cell (with map positions) and the summary by survey and
zone. `python HAT_dune_toes_vs_dunelines.py maps` redraws the imagery figure from the saved table.

Results (cells with a kept toe), with the toe picker as it stands (chord to the first dune, since 2026-11-03),
and in brackets the earlier picker (chord to the highest crest):

| Line vs lidar | Line landward of the toe, median | Middle half | Within 20 m | GIS 32-35, 63-67 |
|---|---|---|---|---|
| 1997 vs 1997 | +10 m (+9) | +5 to +15 m | 80% (72%) | +6 m (-68) |
| 2004 vs 2004 | +19 m (+17) | +4 to +30 m | 44% (43%) | -11 m |
| 2009 vs 2009 | +17 m (+16) | +10 to +26 m | 58% (50%) | +9 m |
| 2023 vs 2020 | +8 m (+7) | +2 to +14 m | 84% (83%) | +5 m (+4) |

- The traced line sits on the dune face, between toe and crest (4-12 m seaward of the crest, lidar 3.8-4.5 m
  NAVD88 at the line against toes of 2.3-2.9 m): the light/dark break the 1997 metadata names, likely the
  vegetation edge, not the toe. Same dune, a different point on it.
- 1997 and 2023/2020 agree to about one 10 m lidar cell; the Google Earth lines (2004, 2009) sit further back and
  scatter more.
- With the earlier picker GIS ~32-36 and ~63-67 disagreed in every year, the line 50-150 m seaward of the lidar
  toe. Checked with the 1997 lidar contoured on the 1997 photo
  (they register: house outlines and dune ridges fall on the photo's): the 1997 line follows the seaward edge of
  a broad, hummocky dune field, along the 2-3 m contours, and is right; the lidar toe lands 35-90 m inside the field,
  on a low swale between ridges, because the furthest-below-the-chord rule finds the deepest dip on a wide
  multi-ridge profile. There the dry beach the module reads is too wide. By 2020 the shore has moved back to a
  single foredune and the line and toes agree (+1 m at GIS 34, +5 m at GIS 66). Fixed in the picker
  (HAT_dune_toes.first_dune: the chord runs to the seaward-most peak with a 1.0 m dip behind it; a 0.5 m dip
  caught beach hummocks and overshot to +29 m): those stretches now agree like the rest. GIS 66 in 1997 still
  scatters, some cells still finding an inner ridge.

## The survey crest windows against the traced lines (2026-11-04)

`HAT_crest_windows_vs_dunelines.py`: for each survey with a traced line, is the line inside the crest window
(within 5 m), seaward or landward of it? Geometry only, for the start's fixed window, the shoreline-retreat window
(the default since 2026-11-03) and the follow-the-dune trial. Writes `supporting/crest_windows_vs_dunelines.csv`
and `_summary.csv`. Share of cells with the line inside the window, island:

| Survey | Fixed | Shoreline retreat | Follow |
|---|---|---|---|
| 1997 | 0.84 | 0.80 | 0.83 |
| 2004 | 0.74 | 0.71 | 0.69 |
| 2009 | 0.64 | 0.77 | 0.77 |
| 2020 (2023 line) | 0.74 | 0.79 | 0.81 |

The retreat window fixes 2009 at the hotspots (0.88-1.00 inside, against 0.08-0.83 fixed); 15-20% of cells keep
the line seaward of the window in every mode. Weak: GIS 75-77 in 2009 (0.53) and GIS 84-88 in 2020 (0.50).
