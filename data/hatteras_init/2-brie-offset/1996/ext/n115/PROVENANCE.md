# 1996 island offsets, extended reach `n115` (GIS 1 to 115)

Built 2026-09-18 15:02 by `scripts/input_prep/2-brie-offset/1-produce/build_island_offset.py --geometry n115` from `duneline_1997.geojson` (1997 imagery). An EXPERIMENT input (the Pea Island extension, 2026-09-16), not a version: `../../CURRENT` named `v3` (renumbered `v1` on 2026-09-19), which this build reproduces exactly on GIS 1-90 (max |diff| 0.000000 m, checked by step 2).

## What is different from the surveyed build

- The reach is GIS 1 to 115. Domains beyond GIS 1-90 are numbered on the surveyed reach's own line-intersects-polygon join onto the whole-island polygons (`1-barrier3d-domains/domain-geojson/domains_pea_hatteras_120.geojson`, `hat_extension_domains.join_lines`; a domain no polygon covers, GIS 0, by northing bin) and carry the dune line's measured offset; their topography is the shared buffer profile.
- 25 extension domains: 91..115.
- Padded to 145 (15 buffer domains each side, the same slope-and-bridge buffer as the surveyed build, now extrapolating from the extension's ends).

## Files

`Island_Dune_Offsets_1996_PADDED_145.csv` (read by the model when `HAT_GEOMETRY=n115`), `_CASCADE_Input.csv`, `_CASCADE_Input_unpadded.csv`, `_buffer_diagnostic.png`, and `1997_duneline_offset_raw_ext.csv` (the extension raw this was built from; the surveyed raw is the one `v3`, now `v1`, keeps).

## Rebuild

```
python scripts/input_prep/2-brie-offset/1-produce/build_island_offset.py --duneline duneline_1997.geojson --year 1996 --geometry n115
```

## Step output

### duneline_to_raw_offsets.py

```
17 transect(s) in no polygon skipped: LineID [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 127, 318, 594, 620, 621, 622]
  placed 155 transects by polygon join onto domains_pea_hatteras_120.geojson
Dune line : C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\dunelines\duneline_1997.geojson
Transects : C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\transects\transects_100m.geojson  (155 beyond GIS 1-90, numbered 91..121 by northing)
  transects with no crossing : 12  -> LineID [608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619]
  transects crossed >1 times : 0
  domains with <5 transects  : {119: 3, 120: 0, 121: 0}

Wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\raw_offsets\ext\1997_duneline_offset_raw_ext.csv  (155 rows)
```

### island_offset_hybrid.py

```
--- Processing 1996 ---
Input file: ['C:\\Users\\hanna\\PycharmProjects\\CASCADE\\data\\hatteras_init\\2-brie-offset\\raw_offsets\\1997_duneline_offset_raw.csv', 'C:\\Users\\hanna\\PycharmProjects\\CASCADE\\data\\hatteras_init\\2-brie-offset\\raw_offsets\\ext\\1997_duneline_offset_raw_ext.csv']
  115 domains processed.
  Baseline distance = 1953.198 m (min mean).

Geometry n115: GIS 1..115, 115 domains, 0 without a dune line []
  surveyed slice vs 1996/v3: max |diff| 0.000000 m

Unpadded offsets saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\ext\n115\Island_Dune_Offsets_1996_CASCADE_Input_unpadded.csv
Unpadded CASCADE-format file saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\ext\n115\Island_Dune_Offsets_1996_CASCADE_Input.csv

Padding summary (hybrid slope + bridge):
  Slope domains per side : 10
  Bridge domains per side: 5
  D1  boundary check     : left_buf[-1]=6360.76  D1=6245.11  diff=115.6455 m
  D90 boundary check     : right_buf[0]=4310.73  D90=4131.96  diff=178.7767 m
  Left  buffer range : 6360.76 – 7401.57 m
  Right buffer range : 4310.73 – 6512.46 m
  Padded length      : 145 (target: 145)

SUCCESS: Padded CASCADE input saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\ext\n115\Island_Dune_Offsets_1996_PADDED_145.csv

Diagnostic figure saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\ext\n115\Island_Dune_Offsets_1996_buffer_diagnostic.png
```
