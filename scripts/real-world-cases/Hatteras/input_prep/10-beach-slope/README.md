# 10-beach-slope — measuring Hatteras's beach slope

Data: `D:/research/UNC/ch1_hatteras/model_runs/data/hatteras_init/10-beach-slope/` (its README has the results).

The storm series' run-up (`3-env-forcings/3-storms/historical_storm_creation_v3_HAT.py`, Stockdon et al. 2006 R2%)
uses one slope, 0.06, from the NC state team. These scripts measure the island's own. The model does not read them
yet.

| Script | What | Writes |
|---|---|---|
| `HAT_beach_slope_mhw_berm.py` | foreshore slope, MHW 0.36 m to berm 1.7 m NAVD88, per 100 m transect and per domain; surveys 1996 and 1997 by default (`python ... 2000 2004` adds others, keeping those already measured) | `1-foreshore-mhw-to-berm/` |
| `HAT_beach_slope_figure.py` | the 1996 slope along the island and the island median by survey; `island_slope.csv` | `1-foreshore-mhw-to-berm/` |

Method (first pass, 2026-11-23): each transect sampled at 1 m from the ocean end on the survey's grid (2009 datum, as
`site_layer/hat_lidar_strips.py`), a 5 m median filter for spikes, the first point above MHW with real data just
seaward of it (so a data gap is not read as the water's edge), then the first point above the berm; slope = 1.34 m
over that distance, kept if 0.01-0.3. Domains in chunks of five so the 1 m grids stay small.

Decided (Hannah, 2026-11-23): this method; the per-domain values from the 1996 lidar; one island value for the storm
series, the median of the 1996 domain medians (0.099). The later surveys are the check through time.
