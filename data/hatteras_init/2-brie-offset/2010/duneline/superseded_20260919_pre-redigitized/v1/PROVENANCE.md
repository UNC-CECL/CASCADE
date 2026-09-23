# 2010 island offsets, v1

Built 2026-09-15 20:08 by `scripts/input_prep/2-brie-offset/build_island_offset.py` from `duneline_2009.geojson` (2009 imagery, standing in for the 2010 start through `DUNE_LINE_FOR_YEAR`).

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

`island_offset_hybrid.py --year 2010 --version v1`: first row per transect, mean of the transects in each domain, zeroed on the most seaward domain (GIS 76, 2004.07 m from the datum), padded to 120 with the slope-and-bridge buffer. Files: `Island_Dune_Offsets_2010_PADDED_120.csv` (read by the model), `_CASCADE_Input.csv`, `_CASCADE_Input_unpadded.csv`, `_buffer_diagnostic.png`.

## CURRENT

`../CURRENT` = `v1` since this build.
 `hatteras_site_config._island_offset_file(2010)` resolves it; env `HAT_OFFSET_VERSION_2010` outranks the file for one run.

## Rebuild

```
python scripts/input_prep/2-brie-offset/island_offset_hybrid.py --year 2010 --version v1 --raw-file data/hatteras_init/2-brie-offset/2010/duneline/v1/2009_duneline_offset_raw.csv
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
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\v1\Island_Dune_Offsets_2010_CASCADE_Input_unpadded.csv
Unpadded CASCADE-format file saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\v1\Island_Dune_Offsets_2010_CASCADE_Input.csv

Padding summary (hybrid slope + bridge):
  Slope domains per side : 10
  Bridge domains per side: 5
  D1  boundary check     : left_buf[-1]=6333.51  D1=6215.78  diff=117.7300 m
  D90 boundary check     : right_buf[0]=878.37  D90=806.49  diff=71.8767 m
  Left  buffer range : 5045.95 – 7393.08 m
  Right buffer range : 878.37 – 3872.39 m
  Padded length      : 120 (target: 120)

SUCCESS: Padded CASCADE input saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\v1\Island_Dune_Offsets_2010_PADDED_120.csv

Diagnostic figure saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\v1\Island_Dune_Offsets_2010_buffer_diagnostic.png
```
