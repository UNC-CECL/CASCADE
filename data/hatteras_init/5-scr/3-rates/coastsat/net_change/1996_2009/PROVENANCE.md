# coastsat/net_change/1996_2009 -- the observed net change, calibration period

Written by `scripts/input_prep/5-scr/3-rates/coastsat/net_change/coastsat_net_change.py` on 2026-10-05.

## What this is

The **target** of the DEM-to-DEM plan. Each period gets sources and sinks fitted so that the
model's end-minus-start shoreline matches this, not the LRR. The LRR also carries storms, fills
and the groin.

| | |
|---|---|
| start window | `mean_shoreline/1995-10-12_1997-10-12/`, centred on 1996-10-12 |
| end window | `mean_shoreline/2008-08-17_2010-08-17/`, centred on 2009-08-17 |
| between centres | 12.85 yr |
| model run | 13 years, 1 Jan 1996 to 1 Jan 2009 |
| transects in both windows | 905 (dropped, in one window only: 1) |
| domain mean net change, GIS 1-90 | -9.2 m (range -73.7 to +44.2) |
| interior GIS 2-89, smoothed | mean -9.8 m, sd 21.1 m |
| median domain standard error | 1.2 m |

## Conventions

- Seaward positive. Net change = end mean chainage minus start mean chainage, per transect,
  then the mean of a domain's transects.
- `net_change_lowess7_m` is LOWESS over the 90 domain means at frac 7/90, with GIS 1-10
  left raw. This is what `matrix_vs_observed.smoothed` does to the model series, so both sides
  are smoothed alike.
- `se_net_change_m` combines the two window means' standard errors (in quadrature, then over
  the domain's transects). It is sampling noise only, not tide or datum error.
- The model's interval differs from the window-centre interval: the run starts on 1 Jan of the
  DEM year and steps whole years. This is reported here, not corrected.
