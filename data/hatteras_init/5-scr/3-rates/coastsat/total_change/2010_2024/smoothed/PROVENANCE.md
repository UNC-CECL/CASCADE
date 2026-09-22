# 3-rates/coastsat/total_change/2010_2024/smoothed - provenance

Written 2026-09-22 10:20 by scripts/input_prep/5-scr/3-rates/coastsat/total_change/coastsat_total_change.py (--product total_change), beside the raw comparison one level up.

**Total shoreline change** = the 2010-2024 LRR x 14 yr. The rate is fitted on the window it is evaluated over, so nothing is extrapolated.

## Why

The rate the model is graded against is not the raw rate: it is the raw domain mean over GIS 1-10 and a 10-domain alongshore LOESS of the transect values beyond (`cascade_pipeline.coastsat_loess`, imported here rather than re-implemented). The raw projected-vs-observed comparison therefore tests a quantity nobody feeds the model. This one tests the target as it is applied.

LOESS commutes with the x years multiply, so smoothing the RATE and smoothing the DISTANCE are the same operation; nothing here turns on the order. What matters is that both sides get the same pass at the same window, including the GIS 1-10 splice, so no residual is a smoothed quantity minus an unsmoothed one.

## The sweep

| LOESS window | n domains | bias (m) | RMS residual (m) | residual range (m) | sd total (m) | sd observed (m) | sign agreement | r |
|---|---|---|---|---|---|---|---|---|
| raw (none) | 90 | -4.6 | 11.1 | -30.0 to +28.2 | 27.3 | 27.6 | 91% | 0.93 |
| 3 domains / 1.5 km | 90 | -4.2 | 9.5 | -30.0 to +28.2 | 25.8 | 26.0 | 94% | 0.94 |
| 5 domains / 2.5 km | 90 | -4.1 | 9.0 | -30.0 to +28.2 | 24.5 | 24.9 | 93% | 0.95 |
| 10 domains / 5 km | 90 | -3.7 | 8.1 | -30.0 to +28.2 | 22.2 | 22.9 | 96% | 0.95 |

**Read the bias and the RMS residual, not r.** A symmetric smoother strips high-frequency variance that is uncorrelated between the two sides, so r rises with the window whether or not the smoothing is right; the column is named `r_inflated_by_smoothing` in `tables/residual_by_scale.csv` for that reason. The bias is close to smoothing-invariant and is the honest summary of whether the trend over- or under-predicts net change.

A departure that survives the 10-domain window is a place the 2010-2024 rate genuinely fails at the scale the model resolves. A departure that collapses between 3 and 10 domains was transect-scale scatter in the LRR estimate, the observed endpoint, or both.

## The figures

`total_change_smoothed_2010_2024_w<NN>.png` is one panel per window, the rate-derived distance against observed, both through that window's pass -- the per-place question, does the trend hold HERE.

`total_change_smoothed_2010_2024_overlay.png` (Hannah, 2026-09-21) puts every LOESS width's curve on one panel and drops the observed side, which answers the other question: what the window does to the target. It is the rate figure `3-rates/coastsat/lrr/2010_2024/smoothing_windows_2010_2024.png` in metres -- the LOESS commutes with the x 14 yr multiply, so the curves have the same shape and only the units differ. It is drawn because metres is the unit the model and the dune line are read in, not because it is a different field.

## Caveats

- The observed side is the thinner estimate: the LRR is fitted through a median 344 satellite positions per transect, while the endpoint uses 18 positions in 2010 and 33 in 2024. Most of the noise the LOESS is removing is probably observed-side, which is why both sides are smoothed.
- A 10-domain window is 5 km and the beach fills inside this record are 3-5 km wide (2014 at GIS 84-89, 2022 at GIS 6-15, 2022 at GIS 21-28). The fill signature in the residual is smeared at 10 domains and is clearest in the 3-domain panel; do not read its disappearance at 10 as evidence it was noise.
- GIS 1-10 are unsmoothed at every window, so the largest positive residual in the raw product (GIS 1) is carried through unchanged by construction.
