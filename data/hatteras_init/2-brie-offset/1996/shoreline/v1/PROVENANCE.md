# 1996 island offsets from the shoreline, v1

Built 2026-09-22 from the **CoastSat mean shoreline for calendar 1995–1997**,
not from a digitised dune line. The first build of this source, so `v1`
([[feedback-version-numbering-restarts]]: new source data restarts the count).

## The line

| property | value |
|---|---|
| file | `5-scr/1-observations/mean_shoreline/1995_1997/shoreline_mean_1995_1997.geojson` |
| crs | EPSG:26918 |
| feature_type | Shoreline (CoastSat window mean) |
| year | 1995-1997 — a **window**, not a survey date |
| source_type | Landsat 5/7/8 via CoastSat (coastsat.space) |
| method | mean of the positions in calendar 1995–1997 per CoastSat transect (905 of 906 transects), geolocated as `origin + chainage * unit_vector`, strung into one line in alongshore order; no outlier rejection, no smoothing |
| editor | `coastsat_mean_shoreline.py` |
| edit_date | 2026-09-22 |

Median 28 satellite positions per transect, within-window scatter (sd) median
12.9 m, **standard error of a transect mean median 2.4 m** — inside the 10 m
Barrier3D cell. One transect (`usa_NC_0032_0078`, domain 6, 6 positions) fell
below the ten-position minimum and was excluded; its domain kept its other
transects. Full accounting in that folder's `PROVENANCE.md`.

## The intersection (step 1)

`duneline_to_raw_offsets.py` against `transects/transects_100m.geojson`:
**450 transects in GIS 1-90, 0 with no crossing, 0 crossed more than once** —
the same clean pass a digitised line gives. Written to
`raw_offsets/1995_1997_shoreline_offset_raw.csv` and copied here.

The step is the dune-line script, used unchanged, and its `--duneline`
argument and `duneline_file` column therefore name a shoreline file. That is
the argument's name, not a claim about the feature; the `feature_type` and
`year` properties above travel into the raw file and say what it really is.

## The geolocation is validated by physics

Nothing in the chain checks that the geolocation is right — a wrong CRS, the
wrong end of a transect taken as the origin, or a sign error in the unit
vector would all produce a plausible-looking line. The beach does check it.
`ORIG_LEN` is measured from the **offshore** datum, so a smaller station is
more seaward, and the shoreline must sit seaward of the dune line:

```
dune station - shoreline station, over the 450 shared transects
  mean +37.4 m   median +31.8 m   sd 18.4 m   min +13.2 m   max +111.1 m
  shoreline seaward of the dune line on 450 of 450 transects (100%)
```

A median beach width of 32 m, positive everywhere, no domain reversed. A
geolocation error of any of the kinds above would scatter this or flip its
sign. (The dune line was flown 15.4 months after the centre of the
shoreline window, so part of this width is shoreline change rather than beach
— see the window note in the mean-shoreline `PROVENANCE.md`, and the
comparison in `../../comparisons/duneline_vs_shoreline/`.)

## The model input (step 2)

`island_offset_hybrid.py --year 1996 --source shoreline --version v1`: first
row per transect, mean of the transects in each domain, zeroed on the most
seaward domain (**GIS 76, 1907.607 m from the datum**), padded to 120 with the
same slope-and-bridge buffer the dune builds use. Nothing in the averaging,
zeroing or padding differs from a dune build — only the feature the stations
were measured to.

Files: `Island_Shoreline_Offsets_1996_PADDED_120.csv`, `_CASCADE_Input.csv`,
`_CASCADE_Input_unpadded.csv`, `_buffer_diagnostic.png`.

## Against the dune build

| | dune (`../../v1/`) | shoreline (here) |
|---|---|---|
| zeroed on | GIS 77, 1953.198 m | GIS 76, 1907.607 m |
| alongshore range | 6245.1 m | 6218.6 m |

Per domain, `shoreline − dune`: **mean +8.2 m, sd 17.2 m, range −51.0 to
+28.9 m**, correlation **r = 1.0000**.

**Read that carefully.** The near-perfect correlation is not the two features
agreeing about the beach — it is both of them being dominated by the same
thing, the 6.2 km of island curvature around Cape Hatteras. Against a 6.2 km
range, a 17 m sd is 0.3%. The physically interesting quantity is the
*residual*, ±20–50 m, which is roughly beach width and its alongshore
variation. Anyone comparing these two profiles should difference them, not
correlate them.

## What the model reads

Not this. `../PROVENANCE.md` says what would have to be checked first.

## Rebuild

```
python scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline.py

python scripts/input_prep/2-brie-offset/1-produce/duneline_to_raw_offsets.py \
    --duneline data/hatteras_init/5-scr/1-observations/mean_shoreline/1995_1997/shoreline_mean_1995_1997.geojson \
    --out 1995_1997_shoreline_offset_raw.csv

python scripts/input_prep/2-brie-offset/1-produce/island_offset_hybrid.py \
    --year 1996 --source shoreline --version v1 \
    --raw-file data/hatteras_init/2-brie-offset/1996/shoreline/v1/1995_1997_shoreline_offset_raw.csv
```

The last command reads the copy kept here, so this build is reproducible after
`raw_offsets/` moves on. Without `--raw-file` the step resolves the window's
current raw through `hat_topo_version.shoreline_raw_file_for_year(1996)`.

## Step output

### duneline_to_raw_offsets.py

```
  reprojecting shoreline_mean_1995_1997.geojson EPSG:26918 -> EPSG:3725
Transects : transects_100m.geojson  (450 in GIS 1-90)
  transects with no crossing : 0
  transects crossed >1 times : 0

Wrote raw_offsets\1995_1997_shoreline_offset_raw.csv  (450 rows)
```

### island_offset_hybrid.py

```
--- Processing 1996 ---
Input file: raw_offsets\1995_1997_shoreline_offset_raw.csv
  90 domains processed.
  Baseline distance = 1907.607 m (min mean).

Padding summary (hybrid slope + bridge):
  Slope domains per side : 10
  Bridge domains per side: 5
  D1  boundary check     : left_buf[-1]=6331.09  D1=6218.64  diff=112.4478 m
  D90 boundary check     : right_buf[0]=939.97  D90=863.80  diff=76.1767 m
  Left  buffer range : 5056.10 - 7343.12 m
  Right buffer range : 939.97 - 3912.59 m
  Padded length      : 120 (target: 120)
```
