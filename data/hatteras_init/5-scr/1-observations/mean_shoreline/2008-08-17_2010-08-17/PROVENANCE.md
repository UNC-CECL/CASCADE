# mean_shoreline/2008-08-17_2010-08-17 -- the CoastSat window mean, as a line

Written by `scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline.py`
on 2026-09-29.

## What this is

The mean satellite shoreline position over **2008-08-17 to 2010-08-17, ±1 yr of the 2009 USACE NCMP topobathy lidar (CHARTS) (flown 2009-08-10 to 2009-08-24)**, per
CoastSat transect, geolocated and strung into one polyline. It is the
shoreline counterpart of a digitised dune line, and `2-brie-offset` consumes
it the same way.

It is **a window mean, not a survey on a date.** A single satellite pass
carries metres of tide, wave setup and cloud-edge noise; the mean over a
window of passes is the quantity a period can start from.

## The numbers

| | |
|---|---|
| window | 2008-08-17 to 2010-08-17, both days included |
| CoastSat transects in the domain lookup | 906 |
| used (n_obs >= 10) | 906 |
| excluded | 0 |
| positions per transect | median 28, min 10, max 37 |
| within-window scatter (sd) | median 14.6 m |
| standard error of a transect mean | median 2.8 m |
| domains covered | 1-90 (90 of 90) |
| transects per domain | min 4, max 13 |
| vertex spacing along the line | median 50 m, p95 52 m, max 153 m |

The standard error is the number that matters for an island offset: at
~3 m it is inside the 10 m Barrier3D cell, so the alongshore shape of
this line is not sampling noise.

## Four things a reader should know

**1. The geolocation is the method, not a formality.** Aggregated to the 90
domains, raw mean chainage spans 124 m alongshore; the geolocated position
spans 6222 m, and the transect origins alone span 6169 m. The CoastSat
transect origins follow the shore around the cape, so a "shape" read straight
off chainage would be ~98% origin bookkeeping. Every other CoastSat product in
`5-scr` is a *difference* of chainage, where the origin cancels; this one is
not.

**2. The window is centred on the start DEM's lidar survey, not on the calendar.**
The 2009 USACE NCMP topobathy lidar (CHARTS) was flown 2009-08-10 to 2009-08-24 (NOAA InPort item 54934, https://www.fisheries.noaa.gov/inport/item/54934); the window is ±1 yr of the middle of those flights, 2009-08-17. That survey is `0-elevation/2009-2014`, the 2010 start's topography (`2004-start`): 2009 USACE wherever it measured, 2014 Post-Sandy only in its nodata. The line becomes the 2010 shoreline island offset, a snapshot the model starts from beside that topography, so it is dated like the topography (Hannah, 2026-09-29) -- a deliberate exception to [[cascade-period-is-the-calendar-year]], which still governs the rates, total change and the scoring target. The 2010 period's dune line is digitised from imagery flown 2009-05-30 (the 2009 line), **2.6 months before** the centre of this window (2009-08-17).

**3. The years inside the window are not evenly sampled.** Positions behind the included means, by calendar year: 7091 / 8529 / 8441 for 2008 / 2009 / 2010. 2008 and 2010 are partial years: the window runs from 2008-08-17 to 2010-08-17. The plain mean leans toward the better-sampled stretch of the window. A year-balanced mean was offered and not taken (Hannah, 2026-09-22); it would be a small change to `window_means`.

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

No outlier rejection -- the 11-19 m within-window scatter is the beach
moving, and `sd`, `se`, `n_obs` and the date span are in the CSV so a reader
can judge each transect. No smoothing of the line. No gap filling. A transect
below the minimum is excluded and named, never silently dropped.

No transect was excluded.

## Files

| file | what it is |
|---|---|
| `shoreline_mean_2008-08-17_2010-08-17.geojson` | one LineString, EPSG:26918, with the metadata properties a dune line carries |
| `transect_means_2008-08-17_2010-08-17.csv` | per transect: n, mean, sd, se, date span, the geolocated point, domain, included/why not |
| `mean_shoreline_2008-08-17_2010-08-17.png` | the diagnostic figure |
| `mean_shoreline_2008-08-17_2010-08-17_island_outline.png` | panel (a) of the diagnostic alone, the line over the island outline |
| `on_imagery/` | the line and its ±1 sd band on the USGS photographs flown inside the window, at six sites and island-wide (three segments, and one ribbon panel); written by `coastsat_mean_shoreline_on_imagery.py`, which needs the D: drive (see its supporting/ folders) |

Resolved through `hat_observed_rates.mean_shoreline_dir/_geojson/_csv`.
Never type these paths.
