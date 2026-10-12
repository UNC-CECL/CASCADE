# coastsat/net_change/2009_2025 -- the observed net change, test period

Written by `scripts/input_prep/5-scr/3-rates/coastsat/net_change/coastsat_net_change.py` on 2026-10-12.

## What this is

The **target** of the DEM-to-DEM plan. Each period gets sources and sinks fitted so that the
model's end-minus-start shoreline matches this, not the LRR. The LRR also carries storms, fills
and the groin.

| | |
|---|---|
| start window | `mean_shoreline/2009-02-17_2010-02-17/`, centred on 2009-08-17 |
| end window | `mean_shoreline/2024-07-01_2026-06-30/`, centred on 2025-07-01 |
| between centres | 15.87 yr |
| model run | 17 years, 1 Jan 2009 to 1 Jan 2026 |
| transects in both windows | 905 (dropped, in one window only: 1) |
| domain mean net change, GIS 1-90 | +12.6 m (range -36.7 to +149.2) |
| interior GIS 2-89, smoothed | mean +10.7 m, sd 24.5 m |
| median domain standard error | 1.6 m |

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
- The end window runs past the CoastSat record, which stops on 2026-01-13. Its data
  covers 2024-07-01 to 2026-01-13.
