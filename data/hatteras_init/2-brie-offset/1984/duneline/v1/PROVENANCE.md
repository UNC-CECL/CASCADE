# 1984 island offsets, v1

Built 2026-09-15 17:25 by `scripts/input_prep/2-brie-offset/build_island_offset.py` from `duneline_1984.geojson` (a 1984 survey, no stand-in).

## The line

| property | value |
|---|---|
| file | `2-brie-offset/dunelines/duneline_1984.geojson` |
| crs | EPSG:26918 |
| (none) | the file carries no properties |

**No `imagery_date` property.** Add one to the geojson so the date does not have to be recovered later (the 1984 date came from a USGS metadata file, the 2004 date from memory).

## The intersection (step 1)

`duneline_to_raw_offsets.py` against `transects/transects_100m.geojson`: 450 transects in GIS 1-90, 0 with no crossing, 0 crossed more than once. Written to `raw_offsets/1984_duneline_offset_raw.csv` (the vintage's current raw) and copied here as `1984_duneline_offset_raw.csv`.

## The model input (step 2)

`island_offset_hybrid.py --year 1984 --version v1`: first row per transect, mean of the transects in each domain, zeroed on the most seaward domain (GIS 78, 1934.26 m from the datum), padded to 120 with the slope-and-bridge buffer. Files: `Island_Dune_Offsets_1984_PADDED_120.csv` (read by the model), `_CASCADE_Input.csv`, `_CASCADE_Input_unpadded.csv`, `_buffer_diagnostic.png`.

## Against superseded_20260915_flat (step 3)

```
1984 superseded_20260915_flat -> v1
  absolute (datum frame): 90 of 90 domains moved >= 0.5 m; mean -1.28 m, min -1.70 (GIS 1), max -0.86 (GIS 82)
  moved domains: {1: -1.7, 2: -1.6, 3: -1.5, 4: -1.3, 5: -1.6, 6: -1.4, 7: -1.3, 8: -1.5, 9: -1.2, 10: -1.2, 11: -1.2, 12: -1.2, 13: -1.2, 14: -1.2, 15: -1.1, 16: -1.5, 17: -1.4, 18: -1.4, 19: -1.4, 20: -1.2, 21: -1.6, 22: -1.3, 23: -1.3, 24: -1.5, 25: -1.1, 26: -1.3, 27: -1.2, 28: -1.4, 29: -1.6, 30: -1.6, 31: -1.5, 32: -1.4, 33: -1.4, 34: -1.6, 35: -1.4, 36: -1.4, 37: -1.5, 38: -1.2, 39: -1.2, 40: -1.1, 41: -1.3, 42: -1.3, 43: -1.1, 44: -1.2, 45: -1.2, 46: -1.0, 47: -1.1, 48: -1.2, 49: -1.5, 50: -1.2, 51: -1.1, 52: -1.3, 53: -1.2, 54: -1.2, 55: -1.3, 56: -1.1, 57: -1.2, 58: -1.5, 59: -1.2, 60: -1.2, 61: -1.4, 62: -1.3, 63: -1.2, 64: -1.2, 65: -1.4, 66: -1.5, 67: -1.7, 68: -1.2, 69: -1.2, 70: -1.5, 71: -1.2, 72: -1.4, 73: -1.2, 74: -1.2, 75: -1.1, 76: -1.4, 77: -1.4, 78: -1.1, 79: -1.1, 80: -1.0, 81: -1.0, 82: -0.9, 83: -1.1, 84: -1.0, 85: -0.9, 86: -1.1, 87: -1.1, 88: -1.2, 89: -1.0, 90: -1.1}
  baseline shift (v1 - superseded_20260915_flat of the per-year minimum): +1.15 m
  model frame (each zeroed on its own minimum): mean -0.14 m, range -0.55 .. +0.28
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1984\v1\offset_1984_superseded_20260915_flat_vs_v1.csv
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1984\v1\offset_1984_superseded_20260915_flat_vs_v1.png
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1984\v1\supporting\offset_1984_superseded_20260915_flat_vs_v1.pdf
```

`offset_1984_superseded_20260915_flat_vs_v1.csv` and `.png` beside this file (PDF and caption under `supporting/`). In the fixed-datum frame a positive difference is the line moved LANDWARD.

## CURRENT

`../CURRENT` = `v1` since this build.
 `hatteras_site_config._island_offset_file(1984)` resolves it; env `HAT_OFFSET_VERSION_1984` outranks the file for one run.

## Rebuild

```
python scripts/input_prep/2-brie-offset/island_offset_hybrid.py --year 1984 --version v1 --raw-file data/hatteras_init/2-brie-offset/1984/duneline/v1/1984_duneline_offset_raw.csv
```

## Step output

### duneline_to_raw_offsets.py

```
reprojecting duneline_1984.geojson EPSG:26918 -> EPSG:3725
Dune line : C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\1-barrier3d-domains\2-brie-offset/dunelines\duneline_1984.geojson
Transects : C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\transects\transects_100m.geojson  (450 in GIS 1-90)
  transects with no crossing : 0
  transects crossed >1 times : 0

Wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\raw_offsets\1984_duneline_offset_raw.csv  (450 rows)
```

### island_offset_hybrid.py

```
--- Processing 1984 ---
Input file: C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\raw_offsets\1984_duneline_offset_raw.csv
  90 domains processed.
  Baseline distance = 1934.255 m (min mean).

Unpadded offsets saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1984\v1\Island_Dune_Offsets_1984_CASCADE_Input_unpadded.csv
Unpadded CASCADE-format file saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1984\v1\Island_Dune_Offsets_1984_CASCADE_Input.csv

Padding summary (hybrid slope + bridge):
  Slope domains per side : 10
  Bridge domains per side: 5
  D1  boundary check     : left_buf[-1]=6414.93  D1=6300.05  diff=114.8831 m
  D90 boundary check     : right_buf[0]=945.70  D90=864.08  diff=81.6232 m
  Left  buffer range : 5141.45 – 7448.88 m
  Right buffer range : 945.70 – 3987.74 m
  Padded length      : 120 (target: 120)

SUCCESS: Padded CASCADE input saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1984\v1\Island_Dune_Offsets_1984_PADDED_120.csv

Diagnostic figure saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1984\v1\Island_Dune_Offsets_1984_buffer_diagnostic.png
```

### HAT_compare_offset_versions.py

```
1984 superseded_20260915_flat -> v1
  absolute (datum frame): 90 of 90 domains moved >= 0.5 m; mean -1.28 m, min -1.70 (GIS 1), max -0.86 (GIS 82)
  moved domains: {1: -1.7, 2: -1.6, 3: -1.5, 4: -1.3, 5: -1.6, 6: -1.4, 7: -1.3, 8: -1.5, 9: -1.2, 10: -1.2, 11: -1.2, 12: -1.2, 13: -1.2, 14: -1.2, 15: -1.1, 16: -1.5, 17: -1.4, 18: -1.4, 19: -1.4, 20: -1.2, 21: -1.6, 22: -1.3, 23: -1.3, 24: -1.5, 25: -1.1, 26: -1.3, 27: -1.2, 28: -1.4, 29: -1.6, 30: -1.6, 31: -1.5, 32: -1.4, 33: -1.4, 34: -1.6, 35: -1.4, 36: -1.4, 37: -1.5, 38: -1.2, 39: -1.2, 40: -1.1, 41: -1.3, 42: -1.3, 43: -1.1, 44: -1.2, 45: -1.2, 46: -1.0, 47: -1.1, 48: -1.2, 49: -1.5, 50: -1.2, 51: -1.1, 52: -1.3, 53: -1.2, 54: -1.2, 55: -1.3, 56: -1.1, 57: -1.2, 58: -1.5, 59: -1.2, 60: -1.2, 61: -1.4, 62: -1.3, 63: -1.2, 64: -1.2, 65: -1.4, 66: -1.5, 67: -1.7, 68: -1.2, 69: -1.2, 70: -1.5, 71: -1.2, 72: -1.4, 73: -1.2, 74: -1.2, 75: -1.1, 76: -1.4, 77: -1.4, 78: -1.1, 79: -1.1, 80: -1.0, 81: -1.0, 82: -0.9, 83: -1.1, 84: -1.0, 85: -0.9, 86: -1.1, 87: -1.1, 88: -1.2, 89: -1.0, 90: -1.1}
  baseline shift (v1 - superseded_20260915_flat of the per-year minimum): +1.15 m
  model frame (each zeroed on its own minimum): mean -0.14 m, range -0.55 .. +0.28
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1984\v1\offset_1984_superseded_20260915_flat_vs_v1.csv
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1984\v1\offset_1984_superseded_20260915_flat_vs_v1.png
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1984\v1\supporting\offset_1984_superseded_20260915_flat_vs_v1.pdf
```
