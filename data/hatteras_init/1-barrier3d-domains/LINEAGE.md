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

## 2026-09-04 — `1984-start/dune-topo/v6`: the fill decided

`v6` = `v3` + the same 98 rows at the same 38 domains as `v5`, with
`--fill median` instead of `--fill measured`. Every dry measured cell in the
inserted block is kept as measured; only the cells at or below MHW are given a
value, the median of the block's own dry cells. 4765 of 4900 inserted cells are
DEM (v5: 3732). Setbacks, land rows and N are identical to v5 in every domain;
mean land elevation moves 0.001 m median, −0.043 m at worst (GIS 85).

Decided in an interview under two constraints: the extracted interior is not
modified, and the simplest rule with the least fabrication wins. Matched
backdune, an alongshore analogue and a mass-conservative reconstruction were
considered and dropped — `row-insert-scope/HAT_fill_options.txt`, section
DECISION, has the argument.

**Accepted cost, recorded not corrected:** at GIS 85 and 86 the 1984 road sits
inside the 1996 dune crest that the block keeps, and `bulldoze()` hands that
crest to the dune in year 1 (+3.3 m and +2.5 m per dune cell). The fill reaches
the road at no other domain. Details in `1984-start/dune-topo/v6/README.md`.

`v5` is kept as the shipped measured-plus-floor reference. `dune-topo/CURRENT`
now says `v6`; it is inert for the reason recorded above. Nothing was wired into
the forcing tree: `4-mgmt-forcing/.../1984/RoadSetback_1984_dunestart.csv` still
carries the v3 setbacks (GIS 85 and 86 floored to 0).

## 2026-09-04, later — `v6` wired in; `CURRENT` now decides

- `hat_topo_version.resolve_version`: **`CURRENT` outranks the extractor's
  `VERSION` literal** (they were the other way round). The literal is what the
  extractor writes; `CURRENT` is what is read. `1984-start` resolves to `v6`;
  `2004-start` (`CURRENT` = `v1`) is unchanged. A fresh extraction is no longer
  adopted until `CURRENT` says so.
- `4-mgmt-forcing/.../1984/RoadSetback_1984_dunestart.csv` replaced by `v6`'s.
  The v3-measured file it replaced is saved as `v3/RoadSetback_1984_dunestart.csv`.
- `test_backdune` built (`--fill backdune`, same footprint) as the reference arm
  of the row-insert test. Not a candidate; deletable.
- The calibration-tree run `HAT_1984_2004_calibBE_road_bdm_groin` (v1, 2026-09-01)
  was re-run on `v6` in place. Its small outputs and index row are kept in
  `output/experiments/row_insert_test/prior_calibration_run_v1_20260901/`.
- The test itself: four arms (`islandv3`, `islandv5`, `islandv6`,
  `islandbackdune`, all `noreloc`) under `output/raw_runs/`, compared at GIS
  80-90 by `scripts/hatteras_ms/HAT_plot_row_insert_test.py` into
  `output/experiments/row_insert_test/`.

## 2026-09-04, evening — the six-fill set built; default reverted to `v3`

- **Default reverted.** `dune-topo/CURRENT` back to `v3` and the forcing-tree
  `RoadSetback_1984_dunestart.csv` back to the v3-measured file, until the set
  below has been run and compared. The calibration-tree calibBE road run made
  earlier today stays on `v6` (its metadata says so).
- **`v7`** flat backdune platform (`--fill backdune`) — replaces the deleted
  `test_backdune`, identical build. **`v8`** matched backdune, crest kept
  (`--fill matched-crest`, new rule: interior rows 0..N-1 copied seaward).
  **`v9`** matched backdune, crest skipped (`--fill matched-nocrest`, rows
  1..N). All three: same 98 rows / 38 domains as v5, bit-identical to v3 behind
  the block, setback CSV identical to v5's, 0 of 4900 inserted cells at their
  own coordinates.
- **Run registry:** `arm_component` accepts a two-level arm (`row-insert/median`)
  and `arms_holding` enumerates that level, so a set of arms files under one
  folder. Deeper nesting refused.
- **Tooling, not yet run:** `scripts/hatteras_ms/HAT_run_row_insert_set.py`
  (six arms, calibration settings, setback CSV swapped and restored) and
  `HAT_plot_row_insert_set.py` (island-wide skill, every road domain's
  relocations, GIS 80-90 detail) into `output/experiments/row_insert_set/`.
- The arm-tagged four-arm test from earlier today (`islandv*noreloc`,
  `output/experiments/row_insert_test/`) is superseded and will be deleted once
  the set has run.

## 2026-09-04, night — the 1984-start layers renamed to `<base>-<scope>-<fill>`

| old | new |
|---|---|
| `v4` | `v3-blocks-floor` |
| `v5` | `v3-island-floor` |
| `v6` | `v3-island-median` |
| `v7` | `v3-island-platform` |
| `v8` | `v3-island-matchedcrest` |
| `v9` | `v3-island-matchednocrest` |

`v1` and `v3` (extractions) unchanged. Folders renamed in place; each
`RUN_MANIFEST.txt` keeps its build-time header with a rename note appended.
Every script that names a layer by literal was repointed (the set driver and
plotter, the crest-experiment arms, the five input-prep plotters, the scope
report). Run metadata, `run_index.csv` rows and dated reports written before
the rename keep the old names — the map above and
`1984-start/row-insert-scope/DUNE_TOPO_VERSION_GUIDE.md` (new, the version
guide) translate them. Sections of this file above are history and keep the
names they were written with.

## 2026-09-04, last — renumbered to a plain sequence (supersedes the lineage names above)

Hannah's call: `v1` stays the original; everything else numbered from `v2` in
build order. The lineage names lasted about an hour.

| as built | interim | **final** | what |
|---|---|---|---|
| `v3` | `v3` | `v2` | re-pick base (extractor `VERSION = "v2"`; picks file renamed `_v2`) |
| `v4` | `v3-blocks-floor` | `v3` | blocks, measured + floor |
| `v5` | `v3-island-floor` | `v4` | island, measured + floor |
| `v6` | `v3-island-median` | `v5` | island, measured + median |
| `v7` | `v3-island-platform` | `v6` | island, flat platform |
| `v8` | `v3-island-matchedcrest` | `v7` | island, matched, crest kept |
| `v9` | `v3-island-matchednocrest` | `v8` | island, matched, crest skipped |

`CURRENT` = `v2`. Every functional literal repointed again; the v2 folder's
settings/figure files renamed `_v2`. **Sections above this one use the names
current when they were written**; run metadata and arm tags likewise. The map
lives in `1984-start/row-insert-scope/DUNE_TOPO_VERSION_GUIDE.md`.

## 2026-09-04 — the six-fill set RUN; the four-arm test deleted

All six arms exit 0 on the intended versions (`none`=v2, `measured-floor`=v4,
`median`=v5, `platform`=v6, `matched-crest`=v7, `matched-nocrest`=v8), under
`output/raw_runs/row-insert/<arm>/1984_2004/calibBE/`. Comparison in
`output/experiments/row_insert_set/` (README there). The earlier arm-tagged
test (`islandv3noreloc`, `islandv5noreloc`, `islandv6noreloc`,
`islandbackdunenoreloc`), its `run_index.csv` rows, its folder
`output/experiments/row_insert_test/` and `HAT_plot_row_insert_test.py` were
deleted as superseded; the v1 calibration-run snapshot moved to
`output/experiments/row_insert_set/prior_calibration_run_v1_20260901/`.
`CURRENT` is still `v2`; no fill has been adopted from the set.

## 2026-09-04 — relocation comparison by interior

Seventh arm `original` added to the set (`v1` + the v1-era setbacks, saved into
`v1/RoadSetback_1984_dunestart.csv` from the 1984start_v1 archive). Every arm
got a relocation-ON partner (same folder, `reloc` token) and
`HAT_relocation_comparison.py` ran on the seven pairs, groin on, calibBE.
Digest: `output/experiments/row_insert_set/relocation/`; numbers also in
`scripts/hatteras_ms/RELOCATION_COMPARISON_RESULTS.md`. Headline: the insert
moves every emergent relocation 3-10 years later; the fill moves it by at most
one year; `original` reproduces the published 0.30 / 0.40.

## 2026-09-07 — only unmodified topography: the layers `v3`–`v8` and every run on modified topography deleted

Hannah's decision, made in an interview: keep `v1` and `v2` (the two
extractions) and remove every version with inserted rows. Applied in full:

| removed | what | size |
|---|---|---|
| `1984-start/dune-topo/v3`–`v8` | the six layers (`v2` + rows; block scope, then five fills at 38 domains) | 29 MB |
| `output/raw_runs/row-insert/` | the six-fill set, controls `none` (v2) and `original` (v1) included, with relocation-ON partners — 14 runs | 4.2 GB |
| `output/raw_runs/blocksv4*`, `islandv5` | crest-experiment arms on as-built v4/v5 (today's v3/v4) | 0.9 GB |
| `output/raw_runs/blocksdate*`, `blocksdsas*`, `blocksduneline*`, `blocksminimum*`, `pea1989keep*`, `pea1989lower*` | arms whose topography went on 09-03; run outputs now gone too | 3.6 GB |
| `output/comparisons/relocation_1984_2004/row-insert/` | the seven relocation-comparison reports | 71 MB |
| `output/experiments/row_insert_set/` | the set comparison, gifs, logs, digest | 2.6 MB |
| 29 rows of `output/raw_runs/run_index.csv` | the runs above; pre-edit file kept as `run_index_archive_20260907_prepurge.csv` | — |

Kept: `v1`, `v2`, `CURRENT` (= `v2`), the forcing-tree CSV (v2-measured),
`pea1989base*` (v1), all of `row-insert-scope/` (the measurement of N, the
scope report, the fill argument, the figures, and the version guide — now a
record of deleted versions), and
`output/experiments/prior_calibration_run_v1_20260901/` (moved up one level out
of the deleted set folder, not deleted).

**The calibration-tree run `HAT_1984_2004_calibBE_road_bdm_groin` was re-run
on `v2`.** It had been re-run on the median fill (today's v5; metadata `v6`) on
09-04 and was the only modified-topography run in the tree. The rest of the
1984-2004 tree (all presets, all scenarios) still stands on `v1` from
2026-09-01, so that scenario is now the one tree run on `v2` — flagged, not
resolved.

**Scripts kept, arms retired.** `HAT_run_row_insert_set.py`,
`HAT_plot_row_insert_set.py` and `HAT_gif_domain_by_interior.py` keep only
`original` (v1) and `none` (v2); `HAT_run_crest_experiment.py` keeps
`pea1989base`; `HAT_plot_seaward_insert_compare.py` defaults to v1 vs v2; the
scoring scripts default to `pea1989base`. The insert plotters that name a layer
literal (`HAT_plot_b3d_grid`, `HAT_plot_fill_options[_grid]`,
`HAT_plot_insert_three_scales`, `HAT_plot_insert_explainer[_grid]`,
`HAT_plot_where_inserts_occur`) now call `hat_topo_version.require_version()`
first and exit naming what is on disk. `HAT_insert_seaward_rows.py` is
untouched: a layer can be rebuilt from `v2` and
`row-insert-scope/duneline-shift/duneline_retreat_1984_1997.csv` with the
recipe in the guide.

Sizes and reasons: `archive_purge_20260907.csv`.
