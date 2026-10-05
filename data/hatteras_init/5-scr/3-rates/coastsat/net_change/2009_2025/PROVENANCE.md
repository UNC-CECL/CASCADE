# coastsat/net_change/2009_2025 -- the observed net change, test period

Written by `scripts/input_prep/5-scr/3-rates/coastsat/net_change/coastsat_net_change.py` on 2026-10-05.

## What this is

The **target** of the DEM-to-DEM plan. Each period gets sources and sinks fitted so that the
model's end-minus-start shoreline matches this, not the LRR. The LRR also carries storms, fills
and the groin.

| | |
|---|---|
| start window | `mean_shoreline/2008-08-17_2010-08-17/`, centred on 2009-08-17 |
| end window | `mean_shoreline/2025-02-17_2026-02-17/`, centred on 2025-08-17 |
| between centres | 16.00 yr |
| model run | 16 years, 1 Jan 2009 to 1 Jan 2025 |
| transects in both windows | 906 (dropped, in one window only: 0) |
| domain mean net change, GIS 1-90 | +14.1 m (range -41.6 to +158.0) |
| interior GIS 2-89, smoothed | mean +11.7 m, sd 24.0 m |
| median domain standard error | 1.1 m |

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
- The end window is 2025-08-17 +/-6 months, not +/-1 yr, because CoastSat stops on
  2026-01-13. Its data covers 2025-02-17 to 2026-01-13.
