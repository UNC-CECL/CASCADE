# mean_shoreline/2025-02-17_2026-02-17 -- the CoastSat window mean, as a line

Written by `scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline.py`
on 2026-10-05.

## What this is

The mean satellite shoreline position over **2025-02-17 to 2026-02-17**, per
CoastSat transect, geolocated and strung into one polyline. It is the
shoreline counterpart of a digitised dune line, and `2-brie-offset` consumes
it the same way.

It is **a window mean, not a survey on a date.** A single satellite pass
carries metres of tide, wave setup and cloud-edge noise; the mean over a
window of passes is the quantity a period can start from.

## The numbers

| | |
|---|---|
| window | 2025-02-17 to 2026-02-17, both days included |
| CoastSat transects in the domain lookup | 906 |
| used (n_obs >= 10) | 906 |
| excluded | 0 |
| positions per transect | median 35, min 30, max 43 |
| within-window scatter (sd) | median 11.6 m |
| standard error of a transect mean | median 1.9 m |
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

**2. The window is given by dates, not calendar years.**
2025-02-17 to 2026-02-17, centred on 2025-08-18.

**3. The years inside the window are not evenly sampled.** Positions behind the included means, by calendar year: 29419 / 2714 for 2025 / 2026. 2025 and 2026 are partial years: the window runs from 2025-02-17 to 2026-02-17. The plain mean leans toward the better-sampled stretch of the window. A year-balanced mean was offered and not taken (Hannah, 2026-09-22); it would be a small change to `window_means`.

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

No outlier rejection -- the 8-18 m within-window scatter is the beach
moving, and `sd`, `se`, `n_obs` and the date span are in the CSV so a reader
can judge each transect. No smoothing of the line. No gap filling. A transect
below the minimum is excluded and named, never silently dropped.

No transect was excluded.

## Files

| file | what it is |
|---|---|
| `shoreline_mean_2025-02-17_2026-02-17.geojson` | one LineString, EPSG:26918, with the metadata properties a dune line carries |
| `transect_means_2025-02-17_2026-02-17.csv` | per transect: n, mean, sd, se, date span, the geolocated point, domain, included/why not |
| `mean_shoreline_2025-02-17_2026-02-17.png` | the diagnostic figure |
| `mean_shoreline_2025-02-17_2026-02-17_island_outline.png` | panel (a) of the diagnostic alone, the line over the island outline |
| `storm_check/` | were there big storms around this window? The storm record 3 yr either side, ranked in 1984-2024, and the mean without post-storm passes; written by `coastsat_mean_shoreline_storm_check.py` |

Resolved through `hat_observed_rates.mean_shoreline_dir/_geojson/_csv`.
Never type these paths.

## Why this window (2026-10-05)

This is the **test-period target** of the DEM-to-DEM plan (calibration 1996 → 2009, test 2009 → 2025). The test ends on 2025-08-17, the anniversary of the 2009 USACE DEM midpoint. The start and calibration-target means use ±1 yr around their DEM dates. This one uses **±6 months** (2025-02-17 → 2026-02-17), agreed with the advisor, because the CoastSat record ends on 2026-01-13 (2026-01-12 at sites 0044–0049). The portal had nothing newer on 2026-10-05. In practice the mean covers 2025-02-17 → 2026-01-13, so its centre sits about two weeks before 2025-08-17. That offset is reported here and not corrected.
