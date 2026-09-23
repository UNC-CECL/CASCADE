# mean_shoreline/1995_1997 -- the CoastSat window mean, as a line

Written by `scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline.py`
on 2026-09-23.

## What this is

The mean satellite shoreline position over **calendar 1995-1997**, per
CoastSat transect, geolocated and strung into one polyline. It is the
shoreline counterpart of a digitised dune line, and `2-brie-offset` consumes
it the same way.

It is **a window mean, not a survey on a date.** A single satellite pass
carries metres of tide, wave setup and cloud-edge noise; the mean over a
window of passes is the quantity a period can start from.

## The numbers

| | |
|---|---|
| CoastSat transects in the domain lookup | 906 |
| used (n_obs >= 10) | 905 |
| excluded | 1 |
| positions per transect | median 28, min 22, max 31 |
| within-window scatter (sd) | median 12.9 m |
| standard error of a transect mean | median 2.4 m |
| domains covered | 1-90 (90 of 90) |
| transects per domain | min 4, max 13 |
| vertex spacing along the line | median 50 m, p95 52 m, max 210 m |

The standard error is the number that matters for an island offset: at
~2 m it is inside the 10 m Barrier3D cell, so the alongshore shape of
this line is not sampling noise.

## Four things a reader should know

**1. The geolocation is the method, not a formality.** Aggregated to the 90
domains, raw mean chainage spans 124 m alongshore; the geolocated position
spans 6222 m, and the transect origins alone span 6169 m. The CoastSat
transect origins follow the shore around the cape, so a "shape" read straight
off chainage would be ~98% origin bookkeeping. Every other CoastSat product in
`5-scr` is a *difference* of chainage, where the origin cancels; this one is
not.

**2. The window is the calendar span, and it does not match the dune line.**
The 1996 period's dune line is digitised from imagery flown 1997-10-12,
**15.4 months** after the centre of this window (1996-07-01). Per
[[cascade-period-is-the-calendar-year]] the mismatch is reported here rather
than corrected by re-centring the window on the survey date. It matters when
this line is differenced against the dune line: part of any gap is those
fifteen months of shoreline change, not beach width.

**3. The years inside the window are not evenly sampled.** Landsat 7 does not
launch until 1999, so 1997 is thinner than the earlier years -- across a
sample of 80 transects, roughly 735 / 976 / 359 positions for 1995 / 1996 /
1997. The plain mean therefore leans slightly toward the early window. A
year-balanced mean was offered and not taken (Hannah, 2026-09-22); it would be
a small change to `window_means`.

**4. Tidal correction is probable but unrecorded.** The transect layer from
coastsat.space carries `beach_slope`, `cil` and `ciu` -- the per-transect
slope used for tidal correction and its confidence interval -- which is good
evidence these chainages are already tidally corrected. **Nothing in this
repository records a statement from the download, so it is not asserted
here.** The exposure is limited: `island_offset_hybrid.py` zeroes each build
on its own minimum, so a uniform tidal bias cancels completely and only the
alongshore variation in beach slope (0.04-0.06 here) survives, worth a few
metres.

## What was not done to the data

No outlier rejection -- the 8-19 m within-window scatter is the beach
moving, and `sd`, `se`, `n_obs` and the date span are in the CSV so a reader
can judge each transect. No smoothing of the line. No gap filling. A transect
below the minimum is excluded and named, never silently dropped.

### Excluded transects

| transect_id | domain_number | n_obs | excluded_because |
|---|---|---|---|
| usa_NC_0032_0078 | 6 | 6 | n_obs 6 < min_obs 10 |

## Files

| file | what it is |
|---|---|
| `shoreline_mean_1995_1997.geojson` | one LineString, EPSG:26918, with the metadata properties a dune line carries |
| `transect_means_1995_1997.csv` | per transect: n, mean, sd, se, date span, the geolocated point, domain, included/why not |
| `mean_shoreline_1995_1997.png` | the diagnostic figure |
| `mean_shoreline_1995_1997_island_outline.png` | panel (a) of the diagnostic alone, the line over the island outline |
| `on_imagery/` | the line and its ±1 sd band on the 1995, 1996 and 1997 USGS photographs at six sites, in two versions: `line_and_band/`, and `with_positions/` adding every satellite position behind the mean coloured by date; plus the whole island on the 1996 photographs in three segments (`line_and_band/..._island_1996.png`) and as one panel (`line_and_band/..._ribbon_1996.png`); written by `coastsat_mean_shoreline_on_imagery.py` (2026-09-23), sites in `on_imagery/supporting/sites.csv` |

Resolved through `hat_observed_rates.mean_shoreline_dir/_geojson/_csv`.
Never type these paths.
