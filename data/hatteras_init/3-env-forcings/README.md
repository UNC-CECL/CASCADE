# 3-env-forcings — what the ocean does to the island

Sea level, storms, and the records they are derived from.

```
1-records/                   what the forcings are built from, as downloaded
    water_level/             the Duck gauge record, 1984-2024, m NAVD88
        noaa_cache_8651370/  the per-month download cache (git-ignored)
    WIS_raw_data/            the wave record the run-up calculation needs
                             (git-ignored)
    storm_record/            the hurricane history documents the storm
                             series is checked against
    wave_climate_duke/       the wave climate (e_phi_0_OBX_yearly.nc, 381 MB,
                             git-ignored); no script reads it
2-rslr/                      sea level: the record, one fit per WINDOW, figures
    record/                  the NOAA monthly mean-trend download
    fits/                    duck_rslr_rates.csv (one row per window) and the
                             per-window monthly series with trend and residual
    figures/                 the three figures and their CAPTIONS.md
3-storms/
    hindcast_storms/         the MODEL-FACING series, one folder per WINDOW
        1984_2004/  1996_2010/  2004_2024/  2010_2024/
        1984_2024/           the spliced full-record series -- NOT a period, it
                             exists for the full-span sweep
    validation/<window>/     both validators' output, side by side
    figures/                 storm record and panel figures
    Notes                    the max-duration tests (72 h chosen)
archive/
    storms_superseded_20260914/   benton_storms/, testing_storms/ (the
                                  storm_check validators still read
                                  testing_storms/base_storms/)
    figures_roya/                 Roya's storm record figure
```

Grouped by job since 2026-09-18, like 5-scr and 8-overwash-analysis. Before
that the records sat beside the forcings built from them, the WIS export
inside `storms/`, and the storm validation in two places: `storms/storm_check/`,
whose two folders were named `storm_check_<window>.png` because a validator
wrote to a drive-rooted "folder" of that name, and a `validation/` folder
inside the `2004_2024` model-input window. The RSLR record stays in
`2-rslr/record/`, not `1-records/`: `rslr/` was laid out record -> fits ->
figures on 2026-09-15 and is kept whole.

Resolve every path through `scripts/site_layer/hat_env_forcings.py`
(`storm_series_file(start, end)`, `DUCK_GAUGE_FILE`, `RSLR_RATES_CSV`, ...);
`hatteras_site_config.py` builds each period's `storm_file` from it.

## A window, not a year

`3-storms/hindcast_storms/1996_2010/` is named for an **interval**, because a storm
series spans one. Surveys are named for their year instead. Same rule as the
rest of the init tree.

The end year is a **boundary**: the model spends `start..end-1`, so the
1996-2010 series drives a run of 1996 through 2009. The series carries one
extra year for that reason and the last step is never spent.

## Producers

All three live in `scripts/input_prep/3-env-forcings/`, and each writes here:

| Script | Writes |
|---|---|
| `1-records/HAT_download_water_levels.py` | `1-records/water_level/` |
| `2-rslr/duck_rslr_analysis.py` | `2-rslr/fits/`, `2-rslr/figures/` |
| `3-storms/historical_storm_creation_v3_HAT.py` | `3-storms/hindcast_storms/<window>/` |

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

Two retired figure scripts under `3-storms/from_Hannah/storm_creation/storm_figures/`
still point at `data/hatteras_init/storms/hindcast_storms/...`, a tree
renamed before this work began: `HAT_storm_record_figure.py` at
`fixed_storms/`, and `HAT_storm_record_figure_roya.py` at
`roya_storms_v2/1984_2004/`; they were already broken and are left alone rather
than guessed at. The third, `HAT_storm_panel.py`, reads the archived testing
storms through `hat_env_forcings.py` since 2026-09-18.
