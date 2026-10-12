# coastsat/net_change/1996_2025 -- the observed net change, full period

Written by `scripts/input_prep/5-scr/3-rates/coastsat/net_change/coastsat_net_change.py` on 2026-10-12.

## What this is

Context, not a target: the full record from the calibration start window to the test end
window, for reading the full-window run (`scripts/hatteras_ms/HAT_full_window_1996_2025.py`)
against. Because the 2009 mean cancels between the two model periods, a transect's change here is its 1996-2009 change plus
its 2009-2025 change.

| | |
|---|---|
| start window | `mean_shoreline/1996-04-12_1997-04-12/`, centred on 1996-10-12 |
| end window | `mean_shoreline/2024-07-01_2026-06-30/`, centred on 2025-07-01 |
| between centres | 28.72 yr |
| model run | 30 years, 1 Jan 1996 to 1 Jan 2026 |
| transects in both windows | 905 (dropped, in one window only: 1) |
| domain mean net change, GIS 1-90 | +5.1 m (range -84.5 to +163.7) |
| interior GIS 2-89, smoothed | mean +2.6 m, sd 36.8 m |
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
- The end window runs past the CoastSat record, which stops on 2026-01-13. Its data
  covers 2024-07-01 to 2026-01-13.
