# coastsat/net_change/1996_2025/ext -- the full-period net change, GIS 1-120

Written by `scripts/input_prep/5-scr/3-rates/coastsat/net_change/coastsat_net_change_extension.py` on 2026-10-07.

## What this is

Context, not a target: the 1996-2025 net change of `../` carried north past the model domains to
Oregon Inlet. GIS 1-90 are the parent folder's transects unchanged. GIS 91-120 are
the extension domains of `2-transect-frame/transect_domains/transect_domain_lookup_ext.csv`.

| | |
|---|---|
| start window | 1995-10-12 to 1997-10-12, centred on 1996-10-12 |
| end window | 2025-02-17 to 2026-02-17, centred on 2025-08-17 |
| extension transects in both windows | 316 of 316 (dropped: 0) |
| domain mean net change, GIS 1-90 | +4.8 m |
| domain mean net change, GIS 91-120 | -16.3 m (range -75.8 to +69.5) |
| domain mean net change, GIS 1-120 | -0.5 m |

## How the extension was built

- No mean-shoreline window files exist for the extension transects. Their window means are
  computed in memory by `coastsat_mean_shoreline.window_means`, the function that writes
  `1-observations/mean_shoreline/`, with the same 10-position minimum. Nothing is
  written to `mean_shoreline/`.
- `net_change_lowess7_m` is LOWESS at 7 domains over all 120 domain means (frac 7/120),
  GIS 1-10 raw. Near GIS 90 it therefore differs slightly from the parent folder's
  curve, which stops at GIS 90.
- Seaward positive. `se_net_change_m` as in the parent folder: sampling noise only.
- The CoastSat transects south of GIS 1 (Cape Point, `usa_NC_0032_0001`-`0020`) have no domain,
  and the Hatteras Inlet to Cape Point shore (`usa_NC_0030`, `usa_NC_0031`) was never
  downloaded, so neither is here.
