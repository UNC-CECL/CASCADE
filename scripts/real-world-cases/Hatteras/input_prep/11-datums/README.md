# 11-datums — checking the model's MHW and berm

Data and results: `D:/research/UNC/ch1_hatteras/model_runs/data/hatteras_init/11-datums/` (its README records the
sources: MHW 0.36 m NAVD88 confirmed against NOAA's Duck datums, epoch 1983-2001; the berm 1.7 m with no recorded
source; Barrier3D's `beta` not used in these runs).

| Script | What | Writes |
|---|---|---|
| `HAT_berm_elevation.py` | the berm crest per 100 m transect and per domain (first point 0.5 m above MHW where the beach flattens to 0.02 over 5 m; none if it climbs past 4 m first); surveys as arguments, 1996 by default, earlier surveys kept | `1-berm/` |

Profiles are read as `../10-beach-slope/HAT_beach_slope_mhw_berm.py` reads them (its `load`, 1 m samples, 5 m median
filter). First pass (2026-11-23): not yet checked against profiles; the model's berm is unchanged.
