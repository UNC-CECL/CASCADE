# 2010 island offsets (shoreline), v2

Built 2026-09-29 from the CoastSat mean shoreline over **2008-08-17 to 2010-08-17**, ±1 yr of the middle (2009-08-17) of the 2009 USACE NCMP topobathy lidar (CHARTS), flown 2009-08-10 to 2009-08-24 (NOAA InPort 54934) -- the lidar the 2010 start DEM is built on. v1 is the same CoastSat data over calendar 2009–2011. Decided by interview (Hannah, 2026-09-29): the offset is a snapshot the model starts from beside that DEM, so it is dated like the DEM. Same source data, a different averaging window, so the count continues (v2), not a new v1.

```
python scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline.py --centred-on usace_2009
python scripts/input_prep/2-brie-offset/1-produce/duneline_to_raw_offsets.py     --duneline data/hatteras_init/5-scr/1-observations/mean_shoreline/2008-08-17_2010-08-17/shoreline_mean_2008-08-17_2010-08-17.geojson     --out 2008-08-17_2010-08-17_shoreline_offset_raw.csv
python scripts/input_prep/2-brie-offset/1-produce/island_offset_hybrid.py --year 2010 --version v2 --source shoreline     --raw-file data/hatteras_init/2-brie-offset/raw_offsets/2008-08-17_2010-08-17_shoreline_offset_raw.csv
python scripts/input_prep/2-brie-offset/2-figures/HAT_compare_offset_versions.py --year 2010 --source shoreline --a v1 --b v2     --raw-a data/hatteras_init/2-brie-offset/2010/shoreline/v1/2009_2011_shoreline_offset_raw.csv     --raw-b data/hatteras_init/2-brie-offset/2010/shoreline/v2/2008-08-17_2010-08-17_shoreline_offset_raw.csv
```

The raw file is copied here, `2008-08-17_2010-08-17_shoreline_offset_raw.csv`. `hat_topo_version.SHORELINE_WINDOW_FOR_YEAR` names this window since 2026-09-29, so a build without `--raw-file` reads this raw file from `raw_offsets/`.

## Mean shoreline

906 of 906 CoastSat transects used; median 28 positions per transect (minimum 10, at usa_NC_0032_0078, GIS 6 -- exactly the minimum); median standard error of a transect mean 2.8 m. See `5-scr/1-observations/mean_shoreline/2008-08-17_2010-08-17/PROVENANCE.md`, and `5-scr/4-comparisons/mean_shoreline_windows/2010/` for the line against the calendar one.

## What validated it

The mean shoreline sits seaward of the dune line on 450 of 450 transects, median 43.4 m, range 10.0–111.0 m (the 2009 dune line, flown 2009-05-30, inside the window). Every transect crossed the line exactly once.

## Build

- Baseline distance: 1976.508 m (v1: 1976.7 m). Real span 0.0–6140.9 m.
- Largest angle between neighbouring domains: 26.1° on the real reach, 27.4° in the buffers; both under the ~42° anti-diffusive limit.
- The padded file is what offset_mode `metres` hands Cascade (`pad_offset_ring`).

## Against v1 (`offset_2010_v1_vs_v2.png`, `.csv`)

83 of 90 domains moved ≥ 0.5 m; datum-frame mean 1.75 m seaward; largest 11.4 m seaward at GIS 25 and 11.1 m landward at GIS 51. The zero domain is GIS 76 in both builds, which moved only 0.2 m, so the model frame tracks the datum frame: model-frame v2 − v1 mean −1.6 m, range −11.3 to +11.3 m (− seaward).

## Against the dune line (`../../comparisons/duneline_vs_shoreline/`)

Under `2009/comparisons/` since 2026-10-06, when the comparison against the CURRENT shoreline build moved up to the year; it sat in this build's `comparisons/` from 2026-09-29.

Drawn with `compare_offset_sources.py --year 2010 --shoreline-version v2` against the dune build `duneline/v1`. The shoreline is seaward of the dune line in 90 of 90 domains.

## CURRENT

`../CURRENT` = `v2` since 2026-09-29 (step 4). `hat_topo_version.SHORELINE_WINDOW_FOR_YEAR` names this build's window the same day, so `island_offset_hybrid.py --source shoreline` without `--raw-file` rebuilds v2. `HAT_OFFSET_VERSION_2010_SHORELINE=v1` selects the calendar build for one run. Until 2026-10-05 the shoreline source was an alternative arm and the default runner read the dune build; since then runs read the shoreline source by default (`hat_topo_version.RUN_OFFSET_SOURCE`), so **this build is what a run on this start reads**, and `HAT_ISLAND_OFFSET_SOURCE=duneline` puts a run back on the dune build. The beach-width check that was set as the condition for that switch is not recorded as done (see `../../../1996/shoreline/PROVENANCE.md`).
