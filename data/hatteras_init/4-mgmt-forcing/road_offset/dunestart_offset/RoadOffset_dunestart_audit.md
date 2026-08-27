# NC-12 road offset measured from the extracted dune start

Generated 2026-08-27T08:13:03 by `HAT_road_offset_from_dune_start.py`.

Covers 1984 on 1984-start/v2, 2004 on 2004-start/v1. Each vintage is measured against interior row 0 of **its own period's extraction** -- the two are different islands, and 65 of 90 domains differ in interior shape between them.

| | |
|---|---|
| Extractor | `HAT_dune_topo_extractor.py` (ALONGSHORE_FLIP=True, STRAIGHTEN=True) |
| 1984 topography | `1984-start/v2` |
| 2004 topography | `2004-start/v1` |
| 1984 DEMs | `C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\1-barrier3d-domains\1984-start\npy-arrays` |
| 2004 DEMs | `C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\1-barrier3d-domains\2004-start\npy-arrays` |
| 1984 picked windows | `HAT_dune_search_windows_v2.json` |
| 2004 picked windows | `HAT_dune_search_windows_v1.json` |
| Reference | interior row 0 = picked dune crest + 1 cell |
| Road reference | seaward-most road cell per profile |
| Alongshore collapse | median over profiles with road |
| Elevation | median of non-water road cells, m MHW |
| Assumed road width | 20 m |

`road_setback` is metres landward of interior row 0 (`roadway_manager.py:99`). The road mask is put through the same orient / alongshore-flip / shear / trim chain as the topography, so the road and the dune line are compared on identical ground.

## Negative setbacks are floored, not measured away

A negative setback means the road falls SEAWARD of interior row 0. It cannot go into the model-facing CSV: `int(-50/10) = -5` and `xyz_interior_grid[-5:-3, :]` indexes from the landward end, so the road would be bulldozed into the bay with no error. `0.0` is not a "no road" sentinel either -- `build_roadway_management_on` uses the road span, not the value. The model-facing file therefore carries `max(setback, 0)`; the true signed value is in `RoadOffset_<year>_domains.csv` under `setback_dunestart_m`.

## Roadways that drown at initialisation are moved seaward

A roadway whose flanking rows are >20% water at t=0 is width-drowned by `roadway_manager.bulldoze` on the first call. `RoadwayManager` sets `_drown_break`, `cascade.py` never calls `update()` again, and that domain spends the whole hindcast as an **unmanaged barrier wearing a road label** -- no overwash removal, no dune rebuilding, no relocation. To keep the domain managed, the model-facing setback is moved to the nearest row seaward that passes bulldoze's own test.

The test transcribed here is bulldoze's: the rows checked are the NEIGHBOURS of the bulldozed band (`road_start - 1`, `road_end + 1`), never the band itself, and every cell counts. There is **no cap** on the distance moved -- the nearest viable row is taken however far it is. `setback_dunestart_m` is never altered; only `RoadSetback_<year>_dunestart.csv` carries the moved value, and every moved domain is flagged `MOVED_SEAWARD`.

### The assumption this rests on

**Nothing drowns at initialisation** on 1984-start/v2, 2004-start/v1 -- 0 domains in either year -- so no setback was moved and the rest of this section does not apply to this run.

That is the point of the gap-filled DEM. On `2009_v4` three domains per year drowned, and the audit for that version showed they drowned on **LiDAR coverage gaps, not measured water**: across the six flanking rows failing at GIS 78/79/80 there were 106 wet cells, of which 105 had never been surveyed and 1 was genuinely wet. The extractor writes no-data back as the water sentinel because Barrier3D has no representation for "unknown", so unsurveyed ground read as ocean and drowned the roadway. Filling those holes from the 2014 NOAA Post-Sandy DEM removes the cause, and the drown count goes to zero.

Those figures are quoted from the `2009_v4` audit and were measured there, not recomputed here -- see `dunestart_offset_ARCHIVE_2009_v4/`, deleted 2026-08-26 - see 1-barrier3d-domains/archive_purge_20260826.csv.

A move that lands inside the domain's own per-profile spread is a re-pick of the alongshore statistic. A move beyond it is a position no profile showed -- a different claim, flagged `BEYOND_MEASURED`.

## `delta_vs_legacy_m` is NOT the dune-line retreat

It is tempting to read the difference against the legacy `RoadSetback_<year>.csv` as the retreat between `<year>` and 2009. It is not, and the reported `corr(delta, retreat)` shows it: a like-for-like pair of measurements taken in two different years would correlate near +1 with the offset-derived retreat, and the measured correlation is strongly negative.

Three reasons the two files are not commensurable, from `HAT_setback_from_lines.py` (retired 2026-08-17, git blob `c37ec03e`):

1. **Different dune feature.** The legacy file computes `setback_i = road_cell_i - dune_cell_i` against a digitized same-year *dune-line geojson*. This script measures against the *DEM dune crest* found inside the picked search window. Those are different features, not the same feature at two dates.
2. **Different frame.** The legacy measurement is taken with both lines "raw, ocean-first" (`HAT_setback_from_lines.py:254`) -- i.e. unstraightened, so it still carries the obliquity smear this script removes.
3. **The legacy file is already floored.** It prints `FLOORED(x) raw setback was negative`, so at the domains where the two disagree most, the legacy value is itself a clamp, not a measurement.

So `delta_vs_legacy_m` is a *migration diagnostic* -- how much each domain's forcing changes if you adopt this method -- and nothing more. Closing the same-year comparison properly needs a same-year dune crest measured the same way, which needs a same-year DEM.

## The road reference is a buffered mask edge, not the road

`road_setback` is measured to the seaward-most cell of the **rasterized mask**, which `HAT_rasterize_road_to_domains.py` burns with `ROAD_BUFFER_M = 6.0` and `ALL_TOUCHED = True` -- roughly a 24 m mask for an ~8 m road. So the reference is not NC-12's centerline, and not its pavement edge either. `HAT_road_buffer_bias.py` measures the difference against the source geojson, in the original raster frame (orient / flip / shear / trim all cancel out of it, since interior row 0 is common to both):

| | 1984 | 2004 |
|---|---|---|
| median bias | -6.8 m | -6.8 m |
| p10 / p90 | -10.8 / -2.7 m | -10.8 / -2.8 m |

Negative means the mask edge sits **seaward** of the centerline, so every setback here runs about 7 m smaller than a centerline measurement would.

### That is close to right, and correcting it would be worse

`bulldoze` lays a 20 m block LANDWARD from `road_start`, so a block centred on the real road wants `road_start = centerline - 10 m`. The buffer supplies `centerline - 6.9 m`, leaving a residual misplacement of **3.1 m -- about a third of a cell**, and an order of magnitude below the 20-40 m p90 placement error the alongshore collapse to one scalar already causes (`HAT_method_comparison_figures.py`).

It is deliberately left alone. Per-profile setbacks are integer cell differences times the cell size, so domain medians land almost exactly on cell boundaries, and `bulldoze` truncates with `int()`. Applying the 3.1 m adjustment moves `road_start` a **full cell (10 m) seaward in 83% of 1984 domains and 90% of 2004 domains** -- trading a 3 m error for a 7 m error in the opposite direction. The buffer is accidentally doing very nearly the right job; the quantization is what would bite.

Two independent paths agree on the number: the raster-frame column comparison above, and the `road_x`/`road_y` columns in `RoadOffset_<year>_profiles.csv`, which are produced by inverting the full index chain back to map coordinates and whose shapely distance to the same geojson is 6.62-6.68 m median (max 12.6 m, no point beyond 60 m). That agreement also validates the index inversion -- an alongshore-flip error would mirror points within each 500 m domain and show up as distances in the hundreds of metres.

D8 is the one extreme outlier (-353 m in 1984). That is the Buxton bend, where NC-12 runs parallel to the raster rows so a "seaward-most cell" has no physical reading -- the same reason D8 is already `EXCLUDED_FROM_SPAN`. Its appearance here is a consistency check passing, not a new problem.

## 1984

- Road span: GIS 9-90 (82 with road, 7 without)
- Setback vs 2009 dune: median 210 m, range -10 to 590 m
- Delta vs the LEGACY `RoadSetback_1984.csv`: median 29 m, range -46 to 112 m
- Offset-derived dune-line retreat 1984->2004: median 43 m
- **corr(delta, retreat) = -0.397** -- see the caveat below.
- NEGATIVE, floored to 0: 1 domain(s)

| GIS | true setback (m) |
|---:|---:|
| 85 | -10 |
- DROWNS at initialisation: 0 domain(s)

### Flags

| GIS | flags |
|---:|---|
| 1 | NO_ROAD |
| 2 | NO_ROAD |
| 3 | NO_ROAD |
| 4 | NO_ROAD |
| 5 | NO_ROAD |
| 6 | NO_ROAD |
| 7 | NO_ROAD |
| 8 | SCATTER_SETBACK(160m),SCATTER_ELEV(0.70),WIDTH(750m),EXCLUDED_FROM_SPAN |
| 11 | SCATTER_ELEV(0.57) |
| 21 | SCATTER_SETBACK(51m) |
| 22 | SCATTER_SETBACK(50m) |
| 26 | SCATTER_SETBACK(60m) |
| 27 | SCATTER_SETBACK(80m) |
| 33 | SCATTER_SETBACK(60m) |
| 35 | SCATTER_SETBACK(41m) |
| 38 | SCATTER_SETBACK(50m) |
| 46 | SCATTER_SETBACK(50m) |
| 64 | SCATTER_SETBACK(51m) |
| 65 | SCATTER_SETBACK(50m) |
| 66 | SCATTER_SETBACK(70m) |
| 67 | SCATTER_SETBACK(70m) |
| 68 | SCATTER_SETBACK(102m) |
| 70 | SCATTER_SETBACK(100m) |
| 71 | SCATTER_SETBACK(51m) |
| 72 | SCATTER_SETBACK(51m) |
| 73 | SCATTER_SETBACK(51m) |
| 75 | SCATTER_SETBACK(60m) |
| 77 | SCATTER_SETBACK(50m) |
| 79 | SCATTER_SETBACK(71m) |
| 80 | SCATTER_SETBACK(41m) |
| 81 | SCATTER_SETBACK(163m) |
| 82 | SCATTER_SETBACK(111m) |
| 83 | SCATTER_SETBACK(50m) |
| 84 | SCATTER_SETBACK(60m),SCATTER_ELEV(1.00) |
| 85 | NEGATIVE(-10->floored 0),SCATTER_SETBACK(50m),SCATTER_ELEV(0.98) |
| 86 | SCATTER_ELEV(0.94) |
| 87 | SCATTER_SETBACK(80m) |
| 88 | SCATTER_SETBACK(51m) |

## 2004

- Road span: GIS 9-90 (82 with road, 7 without)
- Setback vs 2009 dune: median 190 m, range 20 to 600 m
- Delta vs the LEGACY `RoadSetback_2004.csv`: median 22 m, range -49 to 117 m
- NEGATIVE, floored to 0: 0 domain(s)
- DROWNS at initialisation: 0 domain(s)

### Flags

| GIS | flags |
|---:|---|
| 1 | NO_ROAD |
| 2 | NO_ROAD |
| 3 | NO_ROAD |
| 4 | NO_ROAD |
| 5 | NO_ROAD |
| 6 | NO_ROAD |
| 7 | NO_ROAD |
| 8 | SCATTER_SETBACK(166m),SCATTER_ELEV(0.71),EXCLUDED_FROM_SPAN |
| 21 | SCATTER_SETBACK(50m) |
| 26 | SCATTER_SETBACK(50m) |
| 27 | SCATTER_SETBACK(70m) |
| 33 | SCATTER_SETBACK(60m) |
| 48 | SCATTER_SETBACK(51m) |
| 63 | SCATTER_SETBACK(50m) |
| 67 | SCATTER_SETBACK(70m) |
| 68 | SCATTER_SETBACK(101m) |
| 70 | SCATTER_SETBACK(91m) |
| 71 | SCATTER_SETBACK(51m) |
| 72 | SCATTER_SETBACK(50m) |
| 73 | SCATTER_SETBACK(41m) |
| 75 | SCATTER_SETBACK(60m) |
| 80 | SCATTER_SETBACK(71m) |
| 81 | SCATTER_SETBACK(151m) |
| 82 | SCATTER_SETBACK(111m) |
| 83 | SCATTER_SETBACK(80m) |
| 84 | SCATTER_ELEV(0.54) |
| 85 | SCATTER_SETBACK(51m) |
| 87 | SCATTER_SETBACK(41m) |
| 88 | SCATTER_SETBACK(70m) |
