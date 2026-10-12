# 1996 island offsets (shoreline), v3

> **Test build, not CURRENT (2026-10-09).** Built for `output/raw_runs/experiments/target-end-window/2026-10-09-6mo-starts/`, which runs the +/-6-month start against the +/-1-yr one (v2). `SHORELINE_WINDOW_FOR_YEAR` was switched to the +/-6-month window for the build only and then restored, so the table still pairs this start with the v2 window. The mean shoreline was built with **min 5 positions per transect, not 10**, because the 6-month window is thin (at 10, 713 transects and 81 domains; GIS 8-14, 76 and 81 missing). Runs read this build only with `HAT_OFFSET_VERSION_1996_SHORELINE=v3`.

Built 2026-10-09 17:55 by `scripts/input_prep/2-brie-offset/1-produce/build_island_offset.py --source shoreline` from `shoreline_mean_1996-04-12_1997-04-12.geojson` (the CoastSat mean shoreline over 1996-04-12 to 1997-04-12, the window `SHORELINE_WINDOW_FOR_YEAR` pairs with the 1996 start).

## The line

| property | value |
|---|---|
| file | `data/hatteras_init/5-scr/1-observations/mean_shoreline/1996-04-12_1997-04-12/shoreline_mean_1996-04-12_1997-04-12.geojson` |
| crs | EPSG:26918 |
| feature_type | Shoreline (CoastSat window mean) |
| year | 1996-04-12/1997-04-12 |
| imagery_date | 1996-04-12/1997-04-12 |
| source_type | Landsat 5/7/8 via CoastSat (coastsat.space) |
| method | Mean of the satellite positions in 1996-04-12 to 1997-04-12 per CoastSat transect (905 of 906 transects), geolocated as origin + chainage * unit vector and strung into one line in alongshore order. No outlier rejection, no smoothing. |
| editor | coastsat_mean_shoreline.py |
| edit_date | 2026-10-09 00:00:00 |
| notes | Window mean, not a survey. Median 10 positions per transect; see transect_means_1996-04-12_1997-04-12.csv and PROVENANCE.md. |

## The intersection (step 1)

`step1_line_to_raw_offsets.py` against `transects/transects_100m.geojson`: 450 transects in GIS 1-90, 0 with no crossing, 0 crossed more than once. Written to `raw_offsets/1996-04-12_1997-04-12_shoreline_offset_raw.csv` (the line's current raw) and copied here as `1996-04-12_1997-04-12_shoreline_offset_raw.csv`.

## The model input (step 2)

`step2_raw_offsets_to_island_offset.py --year 1996 --version v3 --source shoreline`: first row per transect, mean of the transects in each domain, zeroed on the most seaward domain (GIS 76, 1913.56 m from the datum), padded to 120 with the smooth wrap-around the model uses (`cascade_pipeline.hindcast.pad_offset_ring`), in metres: the padded file is exactly what offset_mode `metres` hands Cascade. Files: `Island_Shoreline_Offsets_1996_PADDED_120.csv` (read by the model), `_CASCADE_Input.csv`, `_CASCADE_Input_unpadded.csv`, `_buffer_diagnostic.png`.

## Against v2 (step 3)

```
1996 v2 -> v3
  model frame (each zeroed on its own minimum): mean +1.29 m, range -3.74 .. +7.56
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v3\offset_1996_v2_vs_v3.csv
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v3\offset_1996_v2_vs_v3.png
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v3\supporting\offset_1996_v2_vs_v3.pdf
```

`offset_1996_v2_vs_v3.csv` and `.png` beside this file (PDF and caption under `supporting/`). In the fixed-datum frame a positive difference is the line moved LANDWARD.

## CURRENT

Built with `--no-current`; `../CURRENT` still names `v2`.
 `hat_topo_version.offset_version(1996, 'shoreline')` resolves it; env `HAT_OFFSET_VERSION_1996_SHORELINE` outranks the file for one run.

## Rebuild

```
python scripts/input_prep/2-brie-offset/1-produce/step2_raw_offsets_to_island_offset.py --year 1996 --version v3 --source shoreline --raw-file data/hatteras_init/2-brie-offset/1996/shoreline/v3/1996-04-12_1997-04-12_shoreline_offset_raw.csv
```

## Step output

### step1_line_to_raw_offsets.py

```
reprojecting shoreline_mean_1996-04-12_1997-04-12.geojson EPSG:26918 -> EPSG:3725
Line      : C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\5-scr\1-observations\mean_shoreline\1996-04-12_1997-04-12\shoreline_mean_1996-04-12_1997-04-12.geojson
Transects : C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\transects\transects_100m.geojson  (450 in GIS 1-90)
  transects with no crossing : 0
  transects crossed >1 times : 0

Wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\raw_offsets\1996-04-12_1997-04-12_shoreline_offset_raw.csv  (450 rows)
```

### step2_raw_offsets_to_island_offset.py

```
--- Processing 1996 ---
Input file: C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\raw_offsets\1996-04-12_1997-04-12_shoreline_offset_raw.csv
  90 domains processed.
  Baseline distance = 1913.556 m (min mean).

Unpadded offsets saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v3\Island_Shoreline_Offsets_1996_CASCADE_Input_unpadded.csv
Unpadded CASCADE-format file saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v3\Island_Shoreline_Offsets_1996_CASCADE_Input.csv

Padding summary (smooth wrap-around, pad_offset_ring):
  Buffer domains per side : 15
  Real span               : 0.0 - 6217.7 m
  Buffer span             : 953.3 - 6287.6 m
  Largest angle, real     : 27.5 deg
  Largest angle, buffer   : 27.7 deg (BRIE goes anti-diffusive past ~42)

SUCCESS: Padded CASCADE input saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v3\Island_Shoreline_Offsets_1996_PADDED_120.csv

Diagnostic figure saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v3\Island_Shoreline_Offsets_1996_buffer_diagnostic.png
```

### HAT_compare_offset_versions.py

```
1996 v2 -> v3
  model frame (each zeroed on its own minimum): mean +1.29 m, range -3.74 .. +7.56
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v3\offset_1996_v2_vs_v3.csv
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v3\offset_1996_v2_vs_v3.png
  wrote C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\1996\shoreline\v3\supporting\offset_1996_v2_vs_v3.pdf
```
