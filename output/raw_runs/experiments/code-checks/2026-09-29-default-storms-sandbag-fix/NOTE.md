# 2026-09-29: default storm file and sandbag flag fixed. Hatteras results unchanged

Hannah, 2026-09-29: "fix 1a and 2 now" (`HATTERAS_CASCADE_CHANGES.md`, sections 5 and 6). Two latent defects in `cascade/cascade.py` and `cascade/cascade_groin.py`, fixed identically in both:

| defect | before | after |
|---|---|---|
| default storm file | `storm_file` defaulted to `cascade-default-original.npy` (`cascade.py`) / `cascade-default-old.npy` (`cascade_groin.py`), names a find-and-replace had made. Neither file exists; the file is `cascade-default-storms.npy`. The guard against using the default storms with another berm/MHW/slope compared against the broken name, so it never fired, and its message read "The default original only apply…". | the default, the guard and the message all say `cascade-default-storms.npy` / "default storms" again, as upstream |
| `sandbag_management_on` | stored as passed but indexed per domain in `update()`, so the default `False` raised `'bool' object is not subscriptable` in year 1 | a single value is broadcast to every domain (as upstream f7ad676b); a per-domain list is kept as given |

## Checks

1. **Defaults now run.** `Cascade("data/default_cascade_variables/", time_step_count=4, alongshore_section_count=2)` with no storm file and no sandbag argument: both classes run 3 years, and the sandbag flag is `[False, False]`. Before the fix, both would have failed at the storm file.
2. **The guard works again.** The same call with `berm_elevation=1.7` raises `CascadeError: The default storms only apply for a berm elevation=1.9 m NAVD88…`.
3. **Hatteras unchanged.** The 2010-2024 edgeBE full-management run was re-run here with the matrix's settings. The runner always passes its own storm file and a per-domain sandbag list. Against `raw_runs/matrix/2010_2024/edgeBE/HAT_2010_2024_edgeBE_offsetmetres_road_bdm_nourish_nogroin`:
   - shoreline matrix identical (max difference 0.0 m)
   - `shoreline_change_rate.csv`, `road_management.csv` and `nourishment_log.csv` identical

No re-run of the matrix is needed.

Run: `2010_2024/edgeBE/HAT_2010_2024_edgeBE_offsetmetres_road_bdm_nourish_nogroin/` (runner, `HAT_RUN_KIND=experiment`); log `output/logs/scratch/codecheck_default_storms_sandbag.log`.
