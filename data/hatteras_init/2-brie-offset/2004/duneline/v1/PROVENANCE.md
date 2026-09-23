# 2004 island offsets, v1

Built 2026-09-15 17:25 by `scripts/input_prep/2-brie-offset/build_island_offset.py` from `duneline_2004.geojson` (a 2004 survey, no stand-in).

## The line

| property | value |
|---|---|
| file | `2-brie-offset/dunelines/duneline_2004.geojson` |
| crs | EPSG:3725 |
| (none) | the file carries no properties |

**No `imagery_date` property.** Add one to the geojson so the date does not have to be recovered later (the 1984 date came from a USGS metadata file, the 2004 date from memory).

## The intersection (step 1)

`duneline_to_raw_offsets.py` against `transects/transects_100m.geojson`: 450 transects in GIS 1-90, 0 with no crossing, 0 crossed more than once. Written to `raw_offsets/2004_duneline_offset_raw.csv` (the vintage's current raw) and copied here as `2004_duneline_offset_raw.csv`.

## The model input (step 2)

`island_offset_hybrid.py --year 2004 --version v1`: first row per transect, mean of the transects in each domain, zeroed on the most seaward domain (GIS 78, 1990.85 m from the datum), padded to 120 with the slope-and-bridge buffer. Files: `Island_Dune_Offsets_2004_PADDED_120.csv` (read by the model), `_CASCADE_Input.csv`, `_CASCADE_Input_unpadded.csv`, `_buffer_diagnostic.png`.

## Against superseded_20260915_flat (step 3)

```
2004 superseded_20260915_flat -> v1
  absolute (datum frame): 90 of 90 domains moved >= 0.5 m; mean -1.01 m, min -1.42 (GIS 30), max -0.75 (GIS 36)
  moved domains: {1: -1.0, 2: -1.2, 3: -1.1, 4: -1.2, 5: -1.2, 6: -1.0, 7: -0.9, 8: -1.1, 9: -0.9, 10: -1.2, 11: -0.9, 12: -0.9, 13: -1.0, 14: -1.0, 15: -1.1, 16: -1.1, 17: -1.1, 18: -1.0, 19: -0.9, 20: -0.9, 21: -1.0, 22: -1.4, 23: -0.9, 24: -1.1, 25: -0.8, 26: -0.9, 27: -0.8, 28: -1.0, 29: -1.2, 30: -1.4, 31: -1.0, 32: -1.1, 33: -1.1, 34: -0.9, 35: -1.2, 36: -0.8, 37: -1.0, 38: -1.1, 39: -0.9, 40: -1.0, 41: -1.0, 42: -1.1, 43: -1.0, 44: -0.9, 45: -1.1, 46: -1.1, 47: -1.0, 48: -1.1, 49: -1.0, 50: -0.9, 51: -1.0, 52: -1.0, 53: -0.8, 54: -0.9, 55: -0.9, 56: -1.0, 57: -0.9, 58: -1.1, 59: -1.0, 60: -1.0, 61: -0.9, 62: -0.8, 63: -0.8, 64: -0.9, 65: -1.4, 66: -0.8, 67: -1.1, 68: -1.1, 69: -1.1, 70: -1.0, 71: -1.1, 72: -0.8, 73: -1.1, 74: -1.0, 75: -0.9, 76: -1.1, 77: -1.0, 78: -0.9, 79: -0.9, 80: -1.1, 81: -1.0, 82: -1.1, 83: -1.2, 84: -1.1, 85: -1.1, 86: -1.1, 87: -0.8, 88: -0.9, 89: -1.1, 90: -0.9}
  baseline shift (v1 - superseded_20260915_flat of the per-year minimum): +0.95 m
  model frame (each zeroed on its own minimum): mean -0.06 m, range -0.47 .. +0.19
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2004\v1\offset_2004_superseded_20260915_flat_vs_v1.csv
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2004\v1\offset_2004_superseded_20260915_flat_vs_v1.png
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2004\v1\supporting\offset_2004_superseded_20260915_flat_vs_v1.pdf
```

`offset_2004_superseded_20260915_flat_vs_v1.csv` and `.png` beside this file (PDF and caption under `supporting/`). In the fixed-datum frame a positive difference is the line moved LANDWARD.

## CURRENT

`../CURRENT` = `v1` since this build.
 `hatteras_site_config._island_offset_file(2004)` resolves it; env `HAT_OFFSET_VERSION_2004` outranks the file for one run.

## Rebuild

```
python scripts/input_prep/2-brie-offset/island_offset_hybrid.py --year 2004 --version v1 --raw-file data/hatteras_init/2-brie-offset/2004/duneline/v1/2004_duneline_offset_raw.csv
```

## Step output

### duneline_to_raw_offsets.py

```
Dune line : C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\1-barrier3d-domains\2-brie-offset/dunelines\duneline_2004.geojson
Transects : C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\transects\transects_100m.geojson  (450 in GIS 1-90)
  transects with no crossing : 0
  transects crossed >1 times : 0

Wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\raw_offsets\2004_duneline_offset_raw.csv  (450 rows)
```

### island_offset_hybrid.py

```
--- Processing 2004 ---
Input file: C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\raw_offsets\2004_duneline_offset_raw.csv
  90 domains processed.
  Baseline distance = 1990.853 m (min mean).

Unpadded offsets saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2004\v1\Island_Dune_Offsets_2004_CASCADE_Input_unpadded.csv
Unpadded CASCADE-format file saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2004\v1\Island_Dune_Offsets_2004_CASCADE_Input.csv

Padding summary (hybrid slope + bridge):
  Slope domains per side : 10
  Bridge domains per side: 5
  D1  boundary check     : left_buf[-1]=6277.30  D1=6163.79  diff=113.5173 m
  D90 boundary check     : right_buf[0]=893.39  D90=817.26  diff=76.1234 m
  Left  buffer range : 5010.77 – 7298.96 m
  Right buffer range : 893.39 – 3866.68 m
  Padded length      : 120 (target: 120)

SUCCESS: Padded CASCADE input saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2004\v1\Island_Dune_Offsets_2004_PADDED_120.csv

Diagnostic figure saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2004\v1\Island_Dune_Offsets_2004_buffer_diagnostic.png
```

### HAT_compare_offset_versions.py

```
2004 superseded_20260915_flat -> v1
  absolute (datum frame): 90 of 90 domains moved >= 0.5 m; mean -1.01 m, min -1.42 (GIS 30), max -0.75 (GIS 36)
  moved domains: {1: -1.0, 2: -1.2, 3: -1.1, 4: -1.2, 5: -1.2, 6: -1.0, 7: -0.9, 8: -1.1, 9: -0.9, 10: -1.2, 11: -0.9, 12: -0.9, 13: -1.0, 14: -1.0, 15: -1.1, 16: -1.1, 17: -1.1, 18: -1.0, 19: -0.9, 20: -0.9, 21: -1.0, 22: -1.4, 23: -0.9, 24: -1.1, 25: -0.8, 26: -0.9, 27: -0.8, 28: -1.0, 29: -1.2, 30: -1.4, 31: -1.0, 32: -1.1, 33: -1.1, 34: -0.9, 35: -1.2, 36: -0.8, 37: -1.0, 38: -1.1, 39: -0.9, 40: -1.0, 41: -1.0, 42: -1.1, 43: -1.0, 44: -0.9, 45: -1.1, 46: -1.1, 47: -1.0, 48: -1.1, 49: -1.0, 50: -0.9, 51: -1.0, 52: -1.0, 53: -0.8, 54: -0.9, 55: -0.9, 56: -1.0, 57: -0.9, 58: -1.1, 59: -1.0, 60: -1.0, 61: -0.9, 62: -0.8, 63: -0.8, 64: -0.9, 65: -1.4, 66: -0.8, 67: -1.1, 68: -1.1, 69: -1.1, 70: -1.0, 71: -1.1, 72: -0.8, 73: -1.1, 74: -1.0, 75: -0.9, 76: -1.1, 77: -1.0, 78: -0.9, 79: -0.9, 80: -1.1, 81: -1.0, 82: -1.1, 83: -1.2, 84: -1.1, 85: -1.1, 86: -1.1, 87: -0.8, 88: -0.9, 89: -1.1, 90: -0.9}
  baseline shift (v1 - superseded_20260915_flat of the per-year minimum): +0.95 m
  model frame (each zeroed on its own minimum): mean -0.06 m, range -0.47 .. +0.19
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2004\v1\offset_2004_superseded_20260915_flat_vs_v1.csv
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2004\v1\offset_2004_superseded_20260915_flat_vs_v1.png
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2004\v1\supporting\offset_2004_superseded_20260915_flat_vs_v1.pdf
```
