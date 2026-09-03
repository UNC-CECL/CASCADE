# Dune/topo extraction — how it changed, and what the version numbers mean

**Version numbers restart at `v1` per PRODUCT.** They are not globally unique.
`1984-start/v1` and `2004-start/v1` are different surfaces built from different
DEMs. A version number alone does not identify a topography — always record the
product with it.

    1-barrier3d-domains/
        1984-start/                 from DEM 2009-2014-1996   dune-topo/v1
        2004-start/                 from DEM 2009-2014        dune-topo/v1
        forecast/                   from a 2025 DEM, later
        buffer/                     shared; no product, no version, no tag
        domain-clips-1m/            LIVE INPUT - 1 m domain clips, read by
                                    HAT_rasterize_road_to_domains.py
        npy-arrays_2009_unfilled/   LIVE INPUT - un-gap-filled 2009 arrays
        control-picks/              LIVE INPUT - unstraightened control windows
        superseded/                 empty; see its README

Resolve paths through `scripts/hat_topo_version.py` — `topo_dirs(product)`,
`array_path(kind, gis, product)`, `domain_arrays(product)`. Never join these
strings by hand.

## Why the numbering was restarted (2026-08-26)

Before the tree went period-first there was ONE topography that both hindcast
periods read, named by the DEM year and a running version: `2009_v1` through
`2009_v5`. The restructure moved `2009_v5` to `2004-start/dune-topo/v5`
unchanged, which left a tree where one product started at v5 and the next
started at v1, and where run metadata said `2009_v5` for a directory called
`v5`. Three names for one thing.

The numbering was restarted so each product counts from 1. The cost is that
the version no longer encodes "fifth attempt" — that history is this file.

## The lineage

| former name | now | DEM | what changed |
|---|---|---|---|
| `2009_v1`, `2009_v2` | *deleted 2026-08-26* | 2009 unfilled | early passes, wrong grid |
| `2009_v3` | *deleted 2026-08-26* | 2009 unfilled | picked against the UNFILLED DEM |
| `2009_v4` | *deleted 2026-08-26* | 2009 unfilled | re-picked with the road drawn on the picker, so a dune-crest argmax could not lock onto the road embankment |
| `2009_v5` | **`2004-start/v1`** | 2009 + 2014 gap fill | a new DEM, not a re-pick. Road drowning at t=0 went 3/yr to 0 |
| — | *`1984-start/v1`, `v2` — deleted 2026-08-27* | 2009 + 2014 + 1996 ALACE | first pick set for the 1996-grafted DEM (v1, picked 2026-08-26) and its bridged post-process (v2). Cleared for a full restart, not superseded by a v3 — see below |

The deleted trees' OUTPUTS are gone; their PICK FILES survive in
`control-picks/`, so each is reproducible from its windows. See
`archive_purge_20260826.csv` for exactly what was removed.

`npy-arrays_2009_unfilled/` (the un-gap-filled 2009 arrays) and
`domain-clips-1m/` (the 1 m domain clips) were moved OUT of `superseded/` on
2026-08-26 — they are live inputs to `HAT_rasterize_road_to_domains.py`, which
builds the road masks the extractor requires, not archives.

## Where `domain-clips-1m/` came from (2026-08-26)

The clips were produced in ArcGIS on **2025-12-01 14:48:21** and imported from

    C:\Users\hanna\OneDrive - University of North Carolina at Chapel Hill\
        Ch1_CASCADE_hatteras\Hatteras_CASCADE_Input\domain_elevation\
        2009_domain_clipresample\domain_<N>\

That path lived only inside the per-raster `.tif.xml` ArcGIS metadata, which was
removed in the second purge of 2026-08-26 — it is recorded here instead so the
origin survives the sidecars. CRS is `NAD_1983_NSRS2007_UTM_Zone_18N` (EPSG 3725)
with `NAVD_1988` vertical (EPSG 5703), which the `.tif` files carry internally.

### What that purge removed, and what it did not

Removed: `.tif.ovr` pyramids and `.aux.xml`/`.tif.xml` sidecars — display and
statistics byproducts that GDAL rebuilds on demand — and clip domains **91–136**
plus unfilled arrays **91–131**, which sit outside the model reach. Both clip
consumers filter to domains 9–90 (`HAT_road_elevation.py:263`,
`HAT_rasterize_road_to_domains.py:184`), so nothing read them.

Kept: **`domain_111`**, because `buffer/README.md` cites it as the checkable
origin of the three buffer arrays, and every **`.tfw`** — world files are the
subject of the orientation proof at `HAT_dune_topo_extractor.py:374`.

No model input changed. See `archive_purge_20260826.csv`.

## What each DEM is

| product | DEM | built by |
|---|---|---|
| `2004-start` | `0-elevation/2009-2014/` — 2009 USACE, gaps filled from 2014 NOAA Post-Sandy | `HAT_dem_gap_fill.py` |
| `1984-start` | `0-elevation/2009-2014-1996/` — the above, plus 1996 NOAA/NASA ALACE overwriting measured ground wherever ALACE has data | `HAT_dem_1984_mosaic.py` |

The 1984-start DEM has **no road boundary** as of 2026-08-26: the 1996 override
used to be confined to the ocean side of the 1984 NC-12 line, and the landward
limit is now the ALACE swath's own edge. See that product's README for the
measurements behind the change.

## Array filenames carry no year

`domain_<N>_topography.npy`, `_dune.npy`, `_nodata.npy`. There is no tag.

They were `domain_<N>_topography_2009.npy` until 2026-08-26. The year was false
for both products — 2004-start is the 2009+2014 mosaic and 1984-start is
2009+2014+1996 — and a per-period tag was tried and reverted the same day: the
tag reached twelve scripts four different ways and no single search found them
all. The period lives in the DIRECTORY, which every reader must resolve anyway.
See the long note at the top of `scripts/hat_topo_version.py`.

## Reading a run's provenance

`run_metadata.json` and `run_index.csv` record **both** `topo_product` and
`topo_dune_version`. Rows written before 2026-08-26 have `topo_dune_version =
2009_v5` and were backfilled with `topo_product = 2004-start`, because
`2009_v5` IS the surface now called `2004-start/v1` — the same arrays under
their pre-restructure name.

That matters for period 1: a **1984** run tagged `2009_v5` read 2004-start's
surface, because before the restructure there was only one topography and both
periods shared it. Those runs are not comparable with runs on `1984-start/v1`.

## Picks are per version, and per product

`<product>/picks/HAT_dune_search_windows_<version>.json`. A version whose pick
file does not exist starts from `default_window()` — it does NOT inherit
another version's windows. Seeding is a deliberate file copy, and the extractor
records `prev_i0`/`prev_i1` when it happens.

## The 1984-start clear (2026-08-27)

`1984-start/dune-topo/v1` and `v2` were deleted outright, along with both of
their window sets, so the product restarts the pick/extract/audit/bridge process
from a blank slate at `v1` rather than carrying its history forward into a `v3`.
The extractor's `VERSION` literal was reset to `"v1"` in the same pass, so
`topo_dirs("1984-start")` raises until the new extraction exists — nothing can
run the 1984 period against a topography that is gone. `2004-start` is
untouched.

This is the same rule the purge above followed: **outputs are deleted, the
inputs that define them are kept.** Both window sets were copied to
`control-picks/HAT_dune_search_windows_1984-start_v{1,2}.json` first, so the
cleared extractions stay reproducible from their windows.

One addition to that rule. The 1996 aerial review — 58 holes of manual imagery
adjudication, the reference that actually decided the v1/v2 dropout verdicts —
was copied to `1984-start/aerial-review/`, **outside `dune-topo/`**. It is keyed
on `(domain, profile)`, which a cross-shore re-pick cannot move, so it survives
any number of re-extractions. Keeping it inside a version directory would have
made it collateral of the next clear, which is exactly the mistake
`superseded/README.md` was written about.

Sizes and file counts: `archive_purge_20260826.csv`, last three rows.

## 2026-09-03 — 1984-start tidy (19.3 MB)

Removed, logged in `archive_purge_20260903.csv`:

| path | files | why |
|---|---:|---|
| `1984-start/dune-topo/v2` | 184 | `v1 + rows, block scope`, built on the PRE-re-pick `v1`. Superseded by `v4` (same insert on the `v3` re-pick). |
| `1984-start/dune-topo-experiments/v1_D85_translate` | 183 | single-domain translate trial; predates the island-width fix; no script, no run |
| `1984-start/dune-topo-experiments/v1_blocks_none_dsas` | 182 | `--variant none`, DSAS shift source; no script, no run |
| `1984-start/dune-topo-experiments/v1_blocks_none_duneline` | 182 | `--variant none`, duneline shift source; no script, no run |

**Run outputs were not touched.** `output/raw_runs/blocksdate{,noreloc}` still
holds v2's results; only re-running from its inputs is now impossible.

Collateral edits, made BEFORE the deletion so nothing was ever left dangling:

- `HAT_run_crest_experiment.py` — `blocksdate` arm removed
- `HAT_plot_seaward_insert_compare.py` — default pair `v1`/`v2` → `v3`/`v5`
- `HAT_plot_fill_options.py`, `HAT_plot_fill_options_grid.py` — `INS_V` `v4` → `v5`
  (N verified identical at all ten block domains; figures byte-identical after)

Kept deliberately: `v1` (pre-re-pick extraction, still wired to the `pea1989*`
arms), `v4` (4.9 MB, the only way to redraw the v4→v5 scope figure), and all
five referenced variants in `dune-topo-experiments/`.

`dune-topo/CURRENT` was set to `v5` to record intent. **It is inert** — the
extractor's `VERSION` literal outranks it and `resolve_version()` returns `v3`.
See `1984-start/dune-topo/README.md`.

### Second pass, same day — consolidation (24.2 MB)

- `1984-start/duneline-shift/` **moved** to `1984-start/row-insert-scope/duneline-shift/`.
  Eight scripts repointed, including the two that write there. The path is now
  resolved once by `hat_topo_version.duneline_shift_dir(product)`; 2004-start
  keeps the plain layout, and that asymmetry lives in that function alone.
- `dune-topo-experiments/` **emptied** — all five remaining variants deleted
  with their arms (`blocksduneline`, `blocksdsas`, `blocksminimum`,
  `pea1989keep`, `pea1989lower`). Run outputs untouched; the scoring and
  plotting scripts read `output/raw_runs/`, so every comparison still works.
- `--shift-source dsas` **removed** from `HAT_insert_seaward_rows.py`. Its path
  was already broken (file had moved to `superseded/`), and the estimate
  measures the shoreline rather than the dune line. `HAT_measure_dsas_shift.py`
  was then DELETED (2026-09-03) - it is tracked in git, so recoverable. Its
  input `5-scr/scr-dsas-1978-2019/dsas_1978_1997_domain_means.csv` is
  untouched.
- `duneline-shift/superseded/` deleted; its record folded into the folder README.
