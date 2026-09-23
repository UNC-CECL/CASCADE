# 2010 island offsets, v1

> **Renumbered 2026-09-19 (Hannah: the numbering restarts at the re-digitized lines).** This build was `v2` until 2026-09-19 and is now `v1`. The log below keeps the old names: its `v2` is this folder, its `v1` is `../superseded_20260919_pre-redigitized/v1/`, and the `_vs_` comparison files compare against that superseded build. The files are byte-identical to what they were under the old name.


Built 2026-09-18 15:01 by `scripts/input_prep/2-brie-offset/1-produce/build_island_offset.py` from `duneline_2009.geojson` (2009 imagery, standing in for the 2010 start through `DUNE_LINE_FOR_YEAR`).

## The line

| property | value |
|---|---|
| file | `2-brie-offset/dunelines/duneline_2009.geojson` |
| crs | EPSG:26918 |
| (none) | the file carries no properties |

**No `imagery_date` property.** Add one to the geojson so the date does not have to be recovered later (the 1984 date came from a USGS metadata file, the 2004 date from memory).

## The intersection (step 1)

`duneline_to_raw_offsets.py` against `transects/transects_100m.geojson`: 450 transects in GIS 1-90, 0 with no crossing, 0 crossed more than once. Written to `raw_offsets/2009_duneline_offset_raw.csv` (the vintage's current raw) and copied here as `2009_duneline_offset_raw.csv`.

## The model input (step 2)

`island_offset_hybrid.py --year 2010 --version v2`: first row per transect, mean of the transects in each domain, zeroed on the most seaward domain (GIS 76, 2004.07 m from the datum), padded to 120 with the slope-and-bridge buffer. Files: `Island_Dune_Offsets_2010_PADDED_120.csv` (read by the model), `_CASCADE_Input.csv`, `_CASCADE_Input_unpadded.csv`, `_buffer_diagnostic.png`.

## Against v1 (step 3)

```
2010 v1 -> v2
  absolute (datum frame): 50 of 90 domains moved >= 0.5 m; mean -6.73 m, min -65.61 (GIS 2), max +0.00 (GIS 7)
  moved domains: {1: -21.4, 2: -65.6, 3: -1.9, 4: -16.6, 5: -6.8, 6: -2.7, 10: -37.6, 11: -16.2, 13: -0.6, 15: -5.6, 16: -21.2, 20: -1.9, 27: -1.4, 28: -1.1, 29: -1.6, 30: -29.2, 31: -18.6, 32: -14.4, 33: -17.3, 34: -13.2, 35: -14.9, 38: -0.6, 42: -2.7, 43: -11.7, 45: -1.6, 46: -5.6, 48: -8.2, 49: -5.4, 50: -1.2, 51: -3.4, 52: -32.2, 53: -43.3, 54: -11.6, 55: -2.6, 60: -1.9, 61: -16.8, 62: -28.4, 63: -4.0, 65: -5.5, 68: -12.1, 69: -16.3, 70: -0.9, 71: -2.5, 73: -6.1, 74: -17.3, 75: -11.8, 84: -0.7, 85: -24.5, 86: -15.3, 87: -1.0}
  baseline shift (v2 - v1 of the per-year minimum): -0.00 m
  model frame (each zeroed on its own minimum): mean -6.73 m, range -65.61 .. +0.00
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\v2\offset_2010_v1_vs_v2.csv
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\v2\offset_2010_v1_vs_v2.png
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\v2\supporting\offset_2010_v1_vs_v2.pdf
```

`offset_2010_v1_vs_v2.csv` and `.png` beside this file (PDF and caption under `supporting/`). In the fixed-datum frame a positive difference is the line moved LANDWARD.

## CURRENT

`../CURRENT` = `v1` (was `v2`) since this build.
 `hatteras_site_config._island_offset_file(2010)` resolves it; env `HAT_OFFSET_VERSION_2010` outranks the file for one run.

## Rebuild

```
python scripts/input_prep/2-brie-offset/1-produce/island_offset_hybrid.py --year 2010 --version v1 --raw-file data/hatteras_init/2-brie-offset/2010/duneline/v1/2009_duneline_offset_raw.csv
```

## Step output

### duneline_to_raw_offsets.py

```
reprojecting duneline_2009.geojson EPSG:26918 -> EPSG:3725
Dune line : C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\dunelines\duneline_2009.geojson
Transects : C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\transects\transects_100m.geojson  (450 in GIS 1-90)
  transects with no crossing : 0
  transects crossed >1 times : 0

Wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\raw_offsets\2009_duneline_offset_raw.csv  (450 rows)
```

### island_offset_hybrid.py

```
--- Processing 2010 ---
Input file: C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\raw_offsets\2009_duneline_offset_raw.csv
  90 domains processed.
  Baseline distance = 2004.070 m (min mean).

Unpadded offsets saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\v2\Island_Dune_Offsets_2010_CASCADE_Input_unpadded.csv
Unpadded CASCADE-format file saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\v2\Island_Dune_Offsets_2010_CASCADE_Input.csv

Padding summary (hybrid slope + bridge):
  Slope domains per side : 10
  Bridge domains per side: 5
  D1  boundary check     : left_buf[-1]=6309.85  D1=6194.40  diff=115.4490 m
  D90 boundary check     : right_buf[0]=878.42  D90=806.49  diff=71.9269 m
  Left  buffer range : 5019.64 – 7348.89 m
  Right buffer range : 878.42 – 3855.01 m
  Padded length      : 120 (target: 120)

SUCCESS: Padded CASCADE input saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\v2\Island_Dune_Offsets_2010_PADDED_120.csv

Diagnostic figure saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\v2\Island_Dune_Offsets_2010_buffer_diagnostic.png
```

### HAT_compare_offset_versions.py

```
2010 v1 -> v2
  absolute (datum frame): 50 of 90 domains moved >= 0.5 m; mean -6.73 m, min -65.61 (GIS 2), max +0.00 (GIS 7)
  moved domains: {1: -21.4, 2: -65.6, 3: -1.9, 4: -16.6, 5: -6.8, 6: -2.7, 10: -37.6, 11: -16.2, 13: -0.6, 15: -5.6, 16: -21.2, 20: -1.9, 27: -1.4, 28: -1.1, 29: -1.6, 30: -29.2, 31: -18.6, 32: -14.4, 33: -17.3, 34: -13.2, 35: -14.9, 38: -0.6, 42: -2.7, 43: -11.7, 45: -1.6, 46: -5.6, 48: -8.2, 49: -5.4, 50: -1.2, 51: -3.4, 52: -32.2, 53: -43.3, 54: -11.6, 55: -2.6, 60: -1.9, 61: -16.8, 62: -28.4, 63: -4.0, 65: -5.5, 68: -12.1, 69: -16.3, 70: -0.9, 71: -2.5, 73: -6.1, 74: -17.3, 75: -11.8, 84: -0.7, 85: -24.5, 86: -15.3, 87: -1.0}
  baseline shift (v2 - v1 of the per-year minimum): -0.00 m
  model frame (each zeroed on its own minimum): mean -6.73 m, range -65.61 .. +0.00
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\v2\offset_2010_v1_vs_v2.csv
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\v2\offset_2010_v1_vs_v2.png
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\v2\supporting\offset_2010_v1_vs_v2.pdf
```
