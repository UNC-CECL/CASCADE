# example — step 2 on 35 real transects

Run either script from this folder:

```
python ../shoreline_rates_template.py    --out expected_output
python ../shoreline_endpoint_template.py --out expected_output
```

Both use the CONFIG default window, 1996 through 2024.

## What is in it

```
data/
    timeseries/usa_NC_0036_timeseries/   35 CoastSat CSVs, one per transect, as downloaded
    transect_zones.csv                   the lookup from ../../1-zone-join/example
expected_output/                         what the two commands above write
```

The time series sit one folder down, the way a CoastSat download arrives; the
scripts search recursively. The lookup is copied unchanged from the step-1
example, so it holds 33 transects, not 35.

## What you should see

| zone | LRR (m/yr) | endpoint rate (m/yr) | net change (m) |
|---|---|---|---|
| 83 | -2.13 | -2.07 | -58.0 |
| 84 | -1.59 | -1.95 | -54.5 |
| 85 | -1.04 | -1.10 | -30.7 |

- **"2 not in the lookup".** `usa_NC_0036_0033` and `_0067` are fitted and
  kept in the per-transect tables with an empty `zone_id`, but averaged into no
  zone. Step 1 refused them, and step 2 respects that.
- **Zone 84's two rates differ** (-1.59 vs -1.95 m/yr). Zone 84 lost about
  30 m in 1996–2001, then retreated more slowly. The LRR is the trend through
  every pass, so the slower later years flatten it; the endpoint is the full
  distance moved, early loss included. That is the point of computing both:
  open `rates_per_transect_1996_2024.csv` and
  `endpoint_per_transect_1996_2024.csv` side by side to see which transects
  carry it.

If your run differs from `expected_output/`, the difference is in your
environment (pandas/scipy version), not in the data.

Time series from CoastSat (coastsat.space), tidally corrected as published.
