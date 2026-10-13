# 1996 island offsets (shoreline), v4

Built 2026-10-12 18:58 by `scripts/input_prep/2-brie-offset/1-produce/build_island_offset.py --source shoreline` from `shoreline_mean_1996-01-01_1996-12-31.geojson` (the CoastSat mean shoreline over 1996-01-01 to 1996-12-31, the window `SHORELINE_WINDOW_FOR_YEAR` pairs with the 1996 start).

## The line

| property | value |
|---|---|
| file | `data/hatteras_init/5-scr/1-observations/mean_shoreline/1996-01-01_1996-12-31/shoreline_mean_1996-01-01_1996-12-31.geojson` |
| crs | EPSG:26918 |
| feature_type | Shoreline (CoastSat window mean) |
| year | 1996-01-01/1996-12-31 |
| imagery_date | 1996-01-01/1996-12-31 |
| source_type | Landsat 5/7/8 via CoastSat (coastsat.space) |
| method | Mean of the satellite positions in 1996-01-01 to 1996-12-31 per CoastSat transect (905 of 906 transects), geolocated as origin + chainage * unit vector and strung into one line in alongshore order. No outlier rejection, no smoothing. |
| editor | coastsat_mean_shoreline.py |
| edit_date | 2026-10-12 00:00:00 |
| notes | Window mean, not a survey. Median 13 positions per transect; see transect_means_1996-01-01_1996-12-31.csv and PROVENANCE.md. |

## The intersection (step 1)

`step1_line_to_raw_offsets.py` against `transects/transects_100m.geojson`: 450 transects in GIS 1-90, 0 with no crossing, 0 crossed more than once. Written to `raw_offsets/1996-01-01_1996-12-31_shoreline_offset_raw.csv` (the line's current raw) and copied here as `1996-01-01_1996-12-31_shoreline_offset_raw.csv`.

## The model input (step 2)

`step2_raw_offsets_to_island_offset.py --year 1996 --version v4 --source shoreline`: first row per transect, mean of the transects in each domain, zeroed on the most seaward domain (GIS 76, 1908.17 m from the datum), padded to 120 with the smooth wrap-around the model uses (`cascade_pipeline.hindcast.pad_offset_ring`), in metres: the padded file is exactly what offset_mode `metres` hands Cascade. Files: `Island_Shoreline_Offsets_1996_PADDED_120.csv` (read by the model), `_CASCADE_Input.csv`, `_CASCADE_Input_unpadded.csv`, `_buffer_diagnostic.png`.

## Against v3 (step 3)

```
1996 v3 -> v4
  model frame (each zeroed on its own minimum): mean +5.46 m, range +0.00 .. +10.72
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v4\offset_1996_v3_vs_v4.csv
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v4\offset_1996_v3_vs_v4.png
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v4\supporting\offset_1996_v3_vs_v4.pdf
```

`offset_1996_v3_vs_v4.csv` and `.png` beside this file (PDF and caption under `supporting/`). In the fixed-datum frame a positive difference is the line moved LANDWARD.

## CURRENT

`../CURRENT` = `v4` since this build.
 `hat_topo_version.offset_version(1996, 'shoreline')` resolves it; env `HAT_OFFSET_VERSION_1996_SHORELINE` outranks the file for one run.

## Rebuild

```
python scripts/input_prep/2-brie-offset/1-produce/step2_raw_offsets_to_island_offset.py --year 1996 --version v4 --source shoreline --raw-file data/hatteras_init/2-brie-offset/1996/shoreline/v4/1996-01-01_1996-12-31_shoreline_offset_raw.csv
```

## Step output

### step1_line_to_raw_offsets.py

```
reprojecting shoreline_mean_1996-01-01_1996-12-31.geojson EPSG:26918 -> EPSG:3725
Line      : C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\5-scr\1-observations\mean_shoreline\1996-01-01_1996-12-31\shoreline_mean_1996-01-01_1996-12-31.geojson
Transects : C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\transects\transects_100m.geojson  (450 in GIS 1-90)
  transects with no crossing : 0
  transects crossed >1 times : 0

Wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\raw_offsets\1996-01-01_1996-12-31_shoreline_offset_raw.csv  (450 rows)
```

### step2_raw_offsets_to_island_offset.py

```
--- Processing 1996 ---
Input file: C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\raw_offsets\1996-01-01_1996-12-31_shoreline_offset_raw.csv
  90 domains processed.
  Baseline distance = 1908.166 m (min mean).

Unpadded offsets saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v4\Island_Shoreline_Offsets_1996_CASCADE_Input_unpadded.csv
Unpadded CASCADE-format file saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v4\Island_Shoreline_Offsets_1996_CASCADE_Input.csv

Padding summary (smooth wrap-around, pad_offset_ring):
  Buffer domains per side : 15
  Real span               : 0.0 - 6220.9 m
  Buffer span             : 958.4 - 6291.2 m
  Largest angle, real     : 27.9 deg
  Largest angle, buffer   : 27.7 deg (BRIE goes anti-diffusive past ~42)

SUCCESS: Padded CASCADE input saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v4\Island_Shoreline_Offsets_1996_PADDED_120.csv

Diagnostic figure saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v4\Island_Shoreline_Offsets_1996_buffer_diagnostic.png
```

### HAT_compare_offset_versions.py

```
1996 v3 -> v4
  model frame (each zeroed on its own minimum): mean +5.46 m, range +0.00 .. +10.72
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v4\offset_1996_v3_vs_v4.csv
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v4\offset_1996_v3_vs_v4.png
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v4\supporting\offset_1996_v3_vs_v4.pdf
```
