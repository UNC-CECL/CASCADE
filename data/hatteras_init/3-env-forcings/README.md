# 3-env-forcings — what the ocean does to the island

Sea level, storms, and the records they are derived from.

```
water_level/             the Duck gauge record, 1984-2024, m NAVD88
    noaa_cache_8651370/  the per-month download cache (git-ignored)
rslr/                    the mean-trend record, and one fit per PERIOD
storms/
    WIS_raw_data/        the wave record the run-up calculation needs
    hindcast_storms/     the model-facing series, one folder per WINDOW
        1984_2004/  1996_2010/  2004_2024/  2010_2024/
        1984_2024/       the spliced full-record series — NOT a period, it
                         exists for the full-span sweep
        old/             earlier attempts, not for use
    storm_check/         validation of a series against the record
    figures/             storm record and panel figures
storm_record/            the hurricane history documents this is checked against
wave_climate_duke/       the wave climate summary
```

## A window, not a year

`hindcast_storms/1996_2010/` is named for an **interval**, because a storm
series spans one. Surveys are named for their year instead. Same rule as the
rest of the init tree.

The end year is a **boundary**: the model spends `start..end-1`, so the
1996-2010 series drives a run of 1996 through 2009. The series carries one
extra year for that reason and the last step is never spent.

## Producers

All three live in `scripts/input_prep/3-env-forcings/`, and each writes here:

| Script | Writes |
|---|---|
| `NOAA_water_level/HAT_download_water_levels.py` | `water_level/` |
| `rslr/duck_rslr_analysis.py` | `rslr/` |
| `storm_creation_final/historical_storm_creation_v3_HAT.py` | `storms/hindcast_storms/<window>/` |

Both of the last two take the window as an argument and derive every path from
it, so the folder name and the file name cannot disagree about what was built:

```
python historical_storm_creation_v3_HAT.py --start-year 1996 --end-year 2010
```

## Moved here 2026-09-12

The gauge record, the sea level fits and the storm figures lived under
`scripts/input_prep/3-env-forcings/`, so more than five hundred data files sat
in the code tree. Three things came out of untangling it:

* **The download cache existed twice**, byte for byte, once in each tree. The
  scripts-side copy was removed; this one survives, and is regenerable by the
  downloader in any case.
* **The sea level script wrote to the current directory**, so where its figures
  landed depended on where it was launched from. It writes here now.
* **Several paths had never resolved on any machine** — drive-rooted literals
  and absolute paths into folder names that were renamed years ago. Those that
  belong to live producers are anchored on their own file now.

Three retired figure scripts under `storm_creation_final/from_Hannah/` still
point at `data/hatteras_init/storms/...`, a tree renamed before this work
began. They were already broken and are left alone rather than guessed at.
