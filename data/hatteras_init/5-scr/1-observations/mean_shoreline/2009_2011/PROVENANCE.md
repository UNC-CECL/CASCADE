# mean_shoreline/2009_2011 -- the CoastSat window mean, as a line

Written by `scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline.py`
on 2026-09-28.

## What this is

The mean satellite shoreline position over **calendar 2009-2011**, per
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
| used (n_obs >= 10) | 906 |
| excluded | 0 |
| positions per transect | median 48, min 22, max 65 |
| within-window scatter (sd) | median 15.3 m |
| standard error of a transect mean | median 2.2 m |
| domains covered | 1-90 (90 of 90) |
| transects per domain | min 4, max 13 |
| vertex spacing along the line | median 50 m, p95 52 m, max 152 m |

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
launch until 1999, so 2011 is thinner than the earlier years -- across a
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

No outlier rejection -- the 12-18 m within-window scatter is the beach
moving, and `sd`, `se`, `n_obs` and the date span are in the CSV so a reader
can judge each transect. No smoothing of the line. No gap filling. A transect
below the minimum is excluded and named, never silently dropped.

No transect was excluded.

## Files

| file | what it is |
|---|---|
| `shoreline_mean_2009_2011.geojson` | one LineString, EPSG:26918, with the metadata properties a dune line carries |
| `transect_means_2009_2011.csv` | per transect: n, mean, sd, se, date span, the geolocated point, domain, included/why not |
| `mean_shoreline_2009_2011.png` | the diagnostic figure |
| `mean_shoreline_2009_2011_island_outline.png` | panel (a) of the diagnostic alone, the line over the island outline |
| `on_imagery/` | the line and its ±1 sd band on the USGS photographs flown inside the window, at six sites and island-wide (three segments, and one ribbon panel); written by `coastsat_mean_shoreline_on_imagery.py`, which needs the D: drive (see its supporting/ folders) |

Resolved through `hat_observed_rates.mean_shoreline_dir/_geojson/_csv`.
Never type these paths.
