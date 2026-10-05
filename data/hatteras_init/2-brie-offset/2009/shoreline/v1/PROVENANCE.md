# 2010 island offsets (shoreline), v1

Built 2026-09-28 for the offset-source comparison over 2010–2024 (Hannah: the same ±1-year window as 1996). The same chain as `../../1996/shoreline/v1`:

```
python scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline.py --window 2009 2011
python scripts/input_prep/2-brie-offset/1-produce/duneline_to_raw_offsets.py \
    --duneline data/hatteras_init/5-scr/1-observations/mean_shoreline/2009_2011/shoreline_mean_2009_2011.geojson \
    --out 2009_2011_shoreline_offset_raw.csv
python scripts/input_prep/2-brie-offset/1-produce/island_offset_hybrid.py --year 2010 --version v1 --source shoreline
```

`hat_topo_version.SHORELINE_WINDOW_FOR_YEAR` gained `2010: (2009, 2011)` the same day, so the build and the runner resolve the raw file through the table. The raw file is copied here, `2009_2011_shoreline_offset_raw.csv`.

## Mean shoreline

Calendar 2009–2011; 906 of 906 CoastSat transects used, 0 excluded. There are a median of 48 positions per transect (minimum 22), and the median standard error of a transect mean is 2.2 m. See `5-scr/1-observations/mean_shoreline/2009_2011/PROVENANCE.md`.

## What validated it

The mean shoreline sits seaward of the 2009 dune line on **450 of 450** transects: median 42.2 m, range 0.4 to 106.5 m. In 1996 the same check gave 450 of 450, median 31.8 m. The 2009 dune line was flown 2009-05-30, inside the window.

## Build

- Baseline distance: 1976.7 m (the minimum domain mean). Real span 0.0–6144.7 m.
- Largest angle between neighbouring domains: 25.1° on the real reach and 27.4° in the buffers, both under the ~42° anti-diffusive limit.
- The padded file is what offset_mode `metres` hands Cascade: the smooth wrap-around from `pad_offset_ring`.

## CURRENT

`../CURRENT` named `v1` until 2026-09-29, when it moved to `v2` (the DEM-centred window); select this build for one run with `HAT_OFFSET_VERSION_2010_SHORELINE=v1`. The source is reached with `HAT_ISLAND_OFFSET_SOURCE=shoreline`. The 2010 dune build the runner reads by default is unchanged.

## Caveat carried from 1996

The Barrier3D interior topography is extracted in the dune-line frame. Before this offset becomes a default, check that the beach width is not counted twice.
