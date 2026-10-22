# 2009-start

Barrier3D domains for the **2009 run** (test, 2009-2025). Named `2004-start/`
until 2026-10-19, when the 2004 start was retired; nothing in the arrays changed.

| | |
|---|---|
| DEM | `0-elevation/2009-2014/` - 2009 USACE, gaps filled from 2014 NOAA Post-Sandy |
| **what loads** | **`v1`** - `CURRENT` says so |
| picks | `1-extraction/picks/HAT_dune_search_windows_v1.json` (not regenerable) |

## Layout

```
1-extraction/
    npy-arrays/          domain_<N>.npy   m NAVD88, -10 nodata   <- extractor INPUT
    npy-arrays_survey/   domain_<N>.npy   provenance codes
    picks/
dune-topo/
    CURRENT              v1
    v1/                  topography/  dunes/  figures/  RUN_MANIFEST.txt
dune-growth-rate/        where the runner's growth rate for this start (0.02-0.04) comes from:
                         README, tables, figures (since 2026-10-21; runs in output/raw_runs/)
```

## History

`v1` is the former `2009-dune-topo/2009_v5`, moved and renumbered, not
rebuilt: the arrays are byte-identical to what every run before 2026-08-25
used. Files inside keep the names they were written with
(`HAT_dune_topo_settings_2009_v5.csv`, old paths in `RUN_MANIFEST.txt`); they
record a run that happened.

Runs published before 2026-08-25 used this topography for the 1984 period too,
because the runner had no per-period concept. Their provenance label reads
`2004-start`, which is where the arrays lived but not what those runs were.

`duneline-shift/` (dune line vs interior row 0, for the retired 2004 start) was
offloaded 2026-10-19 to
`D:\hatteras_ch1\model_runs\data\hatteras_init\1-barrier3d-domains\2004-start\`.
Lineage: `../LINEAGE.md`.
