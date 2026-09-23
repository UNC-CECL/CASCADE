# 1996 island offsets, v1

> **Renumbered 2026-09-19 (Hannah: the numbering restarts at the re-digitized lines).** This build was `v3` until 2026-09-19 and is now `v1`. The log below keeps the old names: its `v3` is this folder, its `v2` is `../superseded_20260919_pre-redigitized/v2/`, and the `_vs_` comparison files compare against that superseded build. The files are byte-identical to what they were under the old name.

> **Files renamed 2026-09-23 (Hannah: the files in a `v1` folder should say v1).** `offset_1996_v2_vs_v3.*` is now `offset_1996_superseded_v2_vs_v1.*` (CSV columns `model_v2_m/model_v3_m` etc. are now `model_superseded_v2_m/model_v1_m` etc.), and every `_buffer_diagnostic_v2` is now `_buffer_diagnostic_1to1`: that `_v2` was the second drawing of the figure, not a build. The step-output logs below keep the names the scripts printed. The comparison figure was redrawn the same day under its new name, so its legend and caption say `superseded_v2` and `v1`: `HAT_compare_offset_versions.py --year 1996 --a superseded_20260919_pre-redigitized/v2 --label-a superseded_v2 --b v1` with the same two raw files; the CSV came out identical.



Built 2026-09-18 15:01 by `scripts/input_prep/2-brie-offset/1-produce/build_island_offset.py` from `duneline_1997.geojson` (1997 imagery, standing in for the 1996 start through `DUNE_LINE_FOR_YEAR`).

## The line

| property | value |
|---|---|
| file | `2-brie-offset/dunelines/duneline_1997.geojson` |
| crs | EPSG:3725 |
| feature_type | Island orientation |
| year | 1997 |
| source_type | Aerial Imagery |
| method | Digitized from light/dark elevation break (no DEM available) |
| editor | H. Henry |
| edit_date | 1762214400000 |
| notes | n/a |

**No `imagery_date` property.** Add one to the geojson so the date does not have to be recovered later (the 1984 date came from a USGS metadata file, the 2004 date from memory).

## The intersection (step 1)

`duneline_to_raw_offsets.py` against `transects/transects_100m.geojson`: 450 transects in GIS 1-90, 0 with no crossing, 0 crossed more than once. Written to `raw_offsets/1997_duneline_offset_raw.csv` (the vintage's current raw) and copied here as `1997_duneline_offset_raw.csv`.

## The model input (step 2)

`island_offset_hybrid.py --year 1996 --version v3`: first row per transect, mean of the transects in each domain, zeroed on the most seaward domain (GIS 77, 1953.20 m from the datum), padded to 120 with the slope-and-bridge buffer. Files: `Island_Dune_Offsets_1996_PADDED_120.csv` (read by the model), `_CASCADE_Input.csv`, `_CASCADE_Input_unpadded.csv`, `_buffer_diagnostic.png`.

## Against v2 (step 3)

```
1996 v2 -> v3
  absolute (datum frame): 32 of 90 domains moved >= 0.5 m; mean -7.79 m, min -66.21 (GIS 35), max +0.19 (GIS 55)
  moved domains: {1: -16.9, 5: -1.2, 9: -1.5, 10: -3.2, 13: -1.5, 16: -16.0, 17: -1.9, 18: -1.9, 30: -6.0, 34: -16.9, 35: -66.2, 36: -25.8, 37: -7.0, 39: -1.0, 46: -41.7, 47: -4.6, 48: -7.7, 49: -17.0, 50: -37.2, 51: -12.9, 62: -25.6, 63: -54.1, 64: -49.5, 65: -42.2, 66: -61.2, 67: -41.0, 76: -23.5, 77: -19.9, 79: -34.8, 80: -37.7, 81: -8.0, 86: -15.0}
  baseline shift (v3 - v2 of the per-year minimum): +3.93 m
  model frame (each zeroed on its own minimum): mean -3.86 m, range -62.28 .. +4.12
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\v3\offset_1996_v2_vs_v3.csv
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\v3\offset_1996_v2_vs_v3.png
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\v3\supporting\offset_1996_v2_vs_v3.pdf
```

`offset_1996_superseded_v2_vs_v1.csv` and `.png` beside this file (PDF and caption under `supporting/`). In the fixed-datum frame a positive difference is the line moved LANDWARD.

## CURRENT

`../CURRENT` = `v1` (was `v3`) since this build.
 `hatteras_site_config._island_offset_file(1996)` resolves it; env `HAT_OFFSET_VERSION_1996` outranks the file for one run.

## Rebuild

```
python scripts/input_prep/2-brie-offset/1-produce/island_offset_hybrid.py --year 1996 --version v1 --raw-file data/hatteras_init/2-brie-offset/1996/duneline/v1/1997_duneline_offset_raw.csv
```

## Step output

### duneline_to_raw_offsets.py

```
Dune line : C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\dunelines\duneline_1997.geojson
Transects : C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\transects\transects_100m.geojson  (450 in GIS 1-90)
  transects with no crossing : 0
  transects crossed >1 times : 0

Wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\raw_offsets\1997_duneline_offset_raw.csv  (450 rows)
```

### island_offset_hybrid.py

```
--- Processing 1996 ---
Input file: C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\raw_offsets\1997_duneline_offset_raw.csv
  90 domains processed.
  Baseline distance = 1953.198 m (min mean).

Unpadded offsets saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\v3\Island_Dune_Offsets_1996_CASCADE_Input_unpadded.csv
Unpadded CASCADE-format file saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\v3\Island_Dune_Offsets_1996_CASCADE_Input.csv

Padding summary (hybrid slope + bridge):
  Slope domains per side : 10
  Bridge domains per side: 5
  D1  boundary check     : left_buf[-1]=6360.76  D1=6245.11  diff=115.6455 m
  D90 boundary check     : right_buf[0]=920.06  D90=844.29  diff=75.7712 m
  Left  buffer range : 5081.74 – 7401.57 m
  Right buffer range : 920.06 – 3921.83 m
  Padded length      : 120 (target: 120)

SUCCESS: Padded CASCADE input saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\v3\Island_Dune_Offsets_1996_PADDED_120.csv

Diagnostic figure saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\v3\Island_Dune_Offsets_1996_buffer_diagnostic.png
```

### HAT_compare_offset_versions.py

```
1996 v2 -> v3
  absolute (datum frame): 32 of 90 domains moved >= 0.5 m; mean -7.79 m, min -66.21 (GIS 35), max +0.19 (GIS 55)
  moved domains: {1: -16.9, 5: -1.2, 9: -1.5, 10: -3.2, 13: -1.5, 16: -16.0, 17: -1.9, 18: -1.9, 30: -6.0, 34: -16.9, 35: -66.2, 36: -25.8, 37: -7.0, 39: -1.0, 46: -41.7, 47: -4.6, 48: -7.7, 49: -17.0, 50: -37.2, 51: -12.9, 62: -25.6, 63: -54.1, 64: -49.5, 65: -42.2, 66: -61.2, 67: -41.0, 76: -23.5, 77: -19.9, 79: -34.8, 80: -37.7, 81: -8.0, 86: -15.0}
  baseline shift (v3 - v2 of the per-year minimum): +3.93 m
  model frame (each zeroed on its own minimum): mean -3.86 m, range -62.28 .. +4.12
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\v3\offset_1996_v2_vs_v3.csv
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\v3\offset_1996_v2_vs_v3.png
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\v3\supporting\offset_1996_v2_vs_v3.pdf
```
