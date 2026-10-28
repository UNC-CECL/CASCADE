# example — step 3 on 35 real transects

Run it from this folder:

```
python ../mean_shoreline_template.py --out expected_output
```

It uses the CONFIG default window, 1995 through 1997.

## What is in it

```
data/
    transects.geojson                    the 35 transects from ../../1-zone-join/example, EPSG:4326
    timeseries/usa_NC_0036_timeseries/   their 35 CoastSat CSVs, from ../../2-shoreline-change/example
expected_output/                         what the command above writes
    mean_shoreline_1995_1997.geojson       the line, EPSG:26918
    mean_shoreline_points_1995_1997.csv    one row per transect: n, mean, sd, se, x, y, used or why not
    mean_shoreline_1995_1997.png           the check map
```

Both inputs are copied unchanged from the earlier examples. This step doesn't
read the step-1 lookup, so all 35 transects are used, including `_0033` and
`_0067`, which the join leaves out.

## What you should see

- **35 used, 0 excluded.** Every transect has 26–28 images in 1995–1997, well
  over `MIN_OBS = 10`.
- **Standard error of a transect mean: median 2.1 m.** The within-window
  scatter (sd) is about 12 m, which is the beach moving with tides, waves and
  seasons; averaging 27 images brings the uncertainty of the mean down to ~2 m.
- **The map** shows one smooth line running from `_0033` (south) to `_0067`
  (north), crossing each grey transect once. A zig-zag would mean the ids are
  not in alongshore order.
- **The points match the full script.** For these 35 transects,
  `../../../1-observations/mean_shoreline/coastsat_mean_shoreline.py --window
  1995 1997` gives the same n and mean chainage, and x, y within 5 mm
  (checked 2026-10-28).

If your run differs from `expected_output/`, the difference is in your
environment (geopandas/pyproj version), not in the data.

Transects and time series from CoastSat (coastsat.space), tidally corrected as published.
