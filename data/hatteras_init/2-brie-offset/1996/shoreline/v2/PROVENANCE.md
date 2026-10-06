# 1996 island offsets (shoreline), v2

Built 2026-09-29 from the CoastSat mean shoreline over **1995-10-12 to 1997-10-12**, ±1 yr of the middle (1996-10-12) of the 1996 fall East Coast NOAA/NASA ALACE lidar, flown 1996-10-09 to 1996-10-16 (NOAA InPort 48147) -- the lidar the 1996 start DEM is built on. v1 is the same CoastSat data over calendar 1995–1997. Decided by interview (Hannah, 2026-09-29): the offset is a snapshot the model starts from beside that DEM, so it is dated like the DEM. Same source data, a different averaging window, so the count continues (v2), not a new v1.

```
python scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline.py --centred-on alace_1996
python scripts/input_prep/2-brie-offset/1-produce/duneline_to_raw_offsets.py     --duneline data/hatteras_init/5-scr/1-observations/mean_shoreline/1995-10-12_1997-10-12/shoreline_mean_1995-10-12_1997-10-12.geojson     --out 1995-10-12_1997-10-12_shoreline_offset_raw.csv
python scripts/input_prep/2-brie-offset/1-produce/island_offset_hybrid.py --year 1996 --version v2 --source shoreline     --raw-file data/hatteras_init/2-brie-offset/raw_offsets/1995-10-12_1997-10-12_shoreline_offset_raw.csv
python scripts/input_prep/2-brie-offset/2-figures/HAT_compare_offset_versions.py --year 1996 --source shoreline --a v1 --b v2     --raw-a data/hatteras_init/2-brie-offset/1996/shoreline/v1/1995_1997_shoreline_offset_raw.csv     --raw-b data/hatteras_init/2-brie-offset/1996/shoreline/v2/1995-10-12_1997-10-12_shoreline_offset_raw.csv
```

The raw file is copied here, `1995-10-12_1997-10-12_shoreline_offset_raw.csv`. `hat_topo_version.SHORELINE_WINDOW_FOR_YEAR` names this window since 2026-09-29, so a build without `--raw-file` reads this raw file from `raw_offsets/`.

## Mean shoreline

905 of 906 CoastSat transects used (usa_NC_0032_0078, GIS 6, excluded: 5 positions); median 20 positions per transect (minimum 16); median standard error of a transect mean 2.4 m. See `5-scr/1-observations/mean_shoreline/1995-10-12_1997-10-12/PROVENANCE.md`, and `5-scr/4-comparisons/mean_shoreline_windows/1996/` for the line against the calendar one.

## What validated it

The mean shoreline sits seaward of the dune line on 450 of 450 transects, median 32.6 m, range 12.0–111.4 m (the 1997 dune line, flown 1997-10-12 -- the window's last day). Every transect crossed the line exactly once.

## Build

- Baseline distance: 1913.868 m (v1: 1907.607 m). Real span 0.0–6213.8 m.
- Largest angle between neighbouring domains: 28.0° on the real reach, 27.7° in the buffers; both under the ~42° anti-diffusive limit.
- The padded file is what offset_mode `metres` hands Cascade (`pad_offset_ring`).

## Against v1 (`offset_1996_v1_vs_v2.png`, `.csv`)

74 of 90 domains moved ≥ 0.5 m; datum-frame mean 0.6 m seaward; largest 9.0 m seaward at GIS 30 and 6.3 m landward at GIS 76. The zero domain is GIS 76 in both builds. It is also the domain that moved most landward (6.3 m), so in the model frame, where each build is zeroed on its own minimum, every other domain reads 6.3 m further seaward than its datum-frame change: model-frame v2 − v1 has mean −6.9 m and range −15.3 to 0.0 m (− seaward). BRIE sees only the shape, so the part that matters is the datum-frame change minus its mean.

## Against the dune line (`../../comparisons/duneline_vs_shoreline/`)

Under `1996/comparisons/` since 2026-10-06, when the comparison against the CURRENT shoreline build moved up to the year; it sat in this build's `comparisons/` from 2026-09-29.

Drawn with `compare_offset_sources.py --year 1996 --shoreline-version v2` against the dune build `duneline/v1`. The shoreline is seaward of the dune line in 90 of 90 domains.

## CURRENT

`../CURRENT` = `v2` since 2026-09-29 (step 4). `hat_topo_version.SHORELINE_WINDOW_FOR_YEAR` names this build's window the same day, so `island_offset_hybrid.py --source shoreline` without `--raw-file` rebuilds v2. `HAT_OFFSET_VERSION_1996_SHORELINE=v1` selects the calendar build for one run. Until 2026-10-05 the shoreline source was an alternative arm and the default runner read the dune build; since then runs read the shoreline source by default (`hat_topo_version.RUN_OFFSET_SOURCE`), so **this build is what a run on this start reads**, and `HAT_ISLAND_OFFSET_SOURCE=duneline` puts a run back on the dune build. The beach-width check that was set as the condition for that switch is not recorded as done (see `../../../1996/shoreline/PROVENANCE.md`).
