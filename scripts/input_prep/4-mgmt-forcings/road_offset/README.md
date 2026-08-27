# road_offset

Where NC-12 sits, per Barrier3D domain — the `road_setback` forcing CASCADE
spends. Two methods are kept, side by side, because the comparison between them
is a result in its own right.

Road *elevation* lives in a sibling folder,
`scripts/input_prep/4-mgmt-forcings/road_elevation/`. It is not an offset
product and carries no year — one elevation set serves both periods. That is
*not* because there is one DEM; there are two, and under the road they differ by
a median **+0.222 m**. That difference is the uncorrected 1996-vs-2009 survey
offset rather than a roadbed, so it is kept out of the forcing. See the note at
`HATTERAS_ROAD_ELEVATION_FILE` in `hatteras_site_config.py`.

```
1-produce/     each method: its producer, and the figures about THAT method
  old_method/    the legacy method, kept for comparison
2-audit/       read-only checks on one method's own output
3-figures/      per-domain views, either method
4-compare/     anything that spans BOTH methods
```

The split that matters is **4-compare/**: work about one method lives with that
method, work that decides *between* them lives on its own, because that is the
open question rather than a routine check. `2-audit/` asks "is this method's
output sound"; `4-compare/` asks "which method should we spend".

## The two methods

| | old | dune-start |
|---|---|---|
| Script | `1-produce/old_method/road_offset_pipeline.py` | `1-produce/HAT_road_offset_from_dune_start.py` |
| Zero point | same-year **digitised dune line** | **interior row 0**, the row `roadway_manager.py:99` indexes against |
| Samples/domain | 5 ArcGIS transects | up to 50 alongshore profiles |
| Statistic | `min(road) − min(dune)`, minima taken independently | per-profile difference, then median |
| Obliquity | uncorrected | mask sheared with the same per-profile shear as the topography |
| Output | `old_method_offset/<year>/` | `dunestart_offset/<year>/` |
| Drowns at t=0 | **0 (1984)**, 5 (2004) | **0, both years** |
| vs the rasterized road | +29 m / +22 m median | 0 m *by construction* |

Both are 2 rows × 82 cols, GIS IDs then metres, so switching is a drop-in.

**What the runner spends today** (`scripts/hatteras_site_config.py:96,110`) is
the **dune-start** method, switched 2026-08-18. `3-figures/HAT_road_domain_views.py`
defaults to `--method dunestart` to match; change both together, or the figures
stop describing the model you run.

**Each vintage is built on its OWN topography** — 1984 on `1984-start/v1`,
2004 on `2004-start/v1`. They are different islands: all 90 domains differ and
**65 differ in interior shape**. The pairing is defined once, as `YEAR_PRODUCT`
in `scripts/hat_topo_version.py`, and every script here resolves through it.
Nothing in this tree may call `topo_dirs()` without naming a product — that
default is what gave both vintages the 2004-start island until 2026-08-26.

## The workflow

```
# old method — one run per vintage
HAT_ROAD_YEAR=1984 python 1-produce/old_method/road_offset_pipeline.py
HAT_ROAD_YEAR=2004 python 1-produce/old_method/road_offset_pipeline.py

# dune-start method — masks first, then the measurement
# The masks are SHARED by both periods (they register to the resampled_*.tif
# grids, which no fill touches), so they are burnt once per road vintage.
HAT_ROAD_YEAR=1984 python 1-produce/HAT_rasterize_road_to_domains.py
HAT_ROAD_YEAR=2004 python 1-produce/HAT_rasterize_road_to_domains.py

# ONE run does BOTH vintages. It configures a separate extractor module per
# product, so you do not repoint HAT_dune_topo_extractor.py between them --
# and must not: this rewrites ONE audit markdown for the whole tree, so a
# half run publishes a document missing a period.
python 1-produce/HAT_road_offset_from_dune_start.py

# figures
python 1-produce/HAT_road_placement_on_domains.py     # per method, both years
python 1-produce/old_method/HAT_old_method_figures.py # how the old number is made

# audits — not optional polish
python 2-audit/HAT_check_geojson_vs_mask.py
python 2-audit/HAT_road_setback_audit.py

# comparison — which method to spend
python 4-compare/HAT_method_comparison_figures.py
python 4-compare/HAT_road_method_diagnostic.py
```

`HAT_check_geojson_vs_mask.py` is the only guard against silently re-burning the
wrong year, and `HAT_road_setback_audit.py` is the only thing that tells you a
domain is an unmanaged barrier wearing a road label for the whole hindcast.

## 1-produce/

| Script | Writes |
|---|---|
| `HAT_rasterize_road_to_domains.py` | `raster/<year>/masks/` — **the only script that masks the road**. Set `HAT_ROAD_YEAR`; do not save a per-year copy |
| `HAT_road_offset_from_dune_start.py` | `dunestart_offset/<year>/` — setback, elevation, per-domain and per-profile detail, audit |
| `HAT_road_placement_on_domains.py` | one PNG per method, into that method's folder |
| `old_method/road_offset_pipeline.py` | `old_method_offset/<year>/RoadSetback_<year>.csv` |
| `old_method/HAT_old_method_figures.py` | figures + audit for the old method's arithmetic |

`HAT_road_placement_on_domains.py` draws every method through one implementation
of `place_road` and the drown test, and `4-compare/` imports it rather than
reimplementing. Add a method by adding a dict entry, not by copying a file — a
drifted method comparison looks like a result.

## 2-audit/

Read-only. All three write tables or markdown, never a forcing file.

| Script | Isolates |
|---|---|
| `HAT_check_geojson_vs_mask.py` | the **rasterization** step alone, in the original raster frame — independent of every extractor transform |
| `HAT_road_setback_audit.py` | the setback against the interior CASCADE really runs, using `FindWidths` transcribed verbatim |

## 3-figures/

One script, three modes:

```
python 3-figures/HAT_road_domain_views.py --domains 52 --year 1984
python 3-figures/HAT_road_domain_views.py --domains drowning --mode map
python 3-figures/HAT_road_domain_views.py --browse --start 52

`--year` picks the topography as well as the setback: 1984 draws on
`1984-start`, 2004 on `2004-start`. It used to resolve one product at import,
before the CLI was parsed, so one of the two years was always drawn on the
other one's island.
```

`--mode map` is the land/water plan view, `--mode section` runs the **real**
`bulldoze()` on a copy of the real interior, `--browse` walks domains
interactively. `--method legacy|dunestart` picks which setback is drawn.

## 4-compare/

Cross-method only. Writes to `road_offset/` itself, not into either method's
folder, because neither output belongs to one method.

| Script | Answers |
|---|---|
| `HAT_method_comparison_figures.py` | where the two methods put the road, and how each compares to the road actually rasterized onto the grid |
| `HAT_road_method_diagnostic.py` | *why* each domain differs — legacy → dune-start split into FRAME (obliquity) and REFERENCE (which feature, which year) |

`HAT_method_comparison_figures.py` imports `place_road` and the drown test from
`1-produce/HAT_road_placement_on_domains.py` rather than reimplementing them, so
a per-method figure and a comparison figure cannot disagree about the same
domain. If that file ever moves, fix the `_PLACEMENT` path — do not copy the
functions across.

## Changelog

**2026-08-26 — each vintage on its own topography.** `1-barrier3d-domains` went
period-first on 2026-08-25 (`1984-start` / `2004-start`), and the setback
producer and the setback audit were updated for it. Five other scripts were
not. They looped over both years against a single module-level `topo_dirs()`
with no argument — `DEFAULT_PRODUCT`, i.e. `2004-start` — so every 1984 panel,
width and diagnostic was computed on the wrong island. Nothing errored.

The scale of it: **all 90 domains differ between the two products and 65 differ
in interior SHAPE** (GIS 11 is 165 rows on `1984-start`, 157 on `2004-start`).
That is four times the v3→v4 incident these files already carry warnings about.

* `YEAR_PRODUCT` is now defined **once**, in `scripts/hat_topo_version.py`, and
  imported by `hatteras_site_config.py` and every script here. It had been
  written out four times and omitted in three.
* `load_interiors()` **requires a year**. That is what surfaced two further
  scripts — `HAT_dunestart_modification_stages.py` and
  `HAT_oceanfloor_offset_check.py` — which had the same defect and no symptom.
* `HAT_road_domain_views.py` binds its topography from `--year` in `main()`
  instead of at import.
* The offset producer runs **both vintages in one invocation**, configuring a
  separate extractor module per product. The previous guard skipped a
  mismatched year, which made every run half a run — and since the audit
  markdown is rewritten whole, the 14:04 run on 2026-08-26 published a
  write-up with **no 2004 section at all** for a forcing that had not changed.
* Column renamed `setback_2009_m` → `setback_dunestart_m` (and `_floored_m`).
  Neither vintage is a 2009 measurement any more: `1984-start` is 2009+2014+1996
  and `2004-start` is 2009+2014.
* **2004 is byte-identical** after the rebuild — `2004-start/v1` is the renamed
  `2009_v5`. The 1984 files moved; the 2004 forcing did not.

**2026-08-26 (same day) — three stale `superseded/` clip paths.**
`domain-clips-1m` moved *into* `superseded/` on 08-25 and back *out* on 08-26.
Three scripts followed it in and not back out, each failing differently:

* `HAT_road_offset_from_dune_start.py` — wrote blank `road_x`/`road_y` for
  every profile, losing the shapely check that validates the index inversion.
* `HAT_check_geojson_vs_mask.py` — reported **0 profiles compared**, a median
  offset of `nan`, and then carried on and printed a report. This README calls
  it the only guard against re-burning the wrong year; it had been passing
  vacuously. It now exits non-zero rather than reporting on nothing. Repaired,
  it puts **100 % of centreline vertices inside the mask footprint**, both
  years, 4125 profiles each.
* `HAT_road_geojson_map.py` — printed "no clip tifs found" and drew nothing.

**2026-08-26 (same day) — road elevation is reproducible again.** `FILL_SOURCE`
was `2008_NOAA_IOCM`, a product deleted from disk and not rebuildable from this
repo, so the script raised on import. It now defaults to the live `2009-2014`
baseline and takes `HAT_ROAD_ELEV_FILL` from the environment. Rebuilding moved
**two domains** — GIS 78 by −0.006 m and GIS 79 by −0.015 m — which is the
bound the script's own header predicted. The `2009-2014-1996` product was
deliberately **not** used: it sits +0.222 m higher through the corridor, and
that is the uncorrected survey offset, not a road.


**2026-08-20 — rebuilt on `2009_v5` (gap-filled DEM).** The previous tree was
archived whole to `dunestart_offset_ARCHIVE_2009_v4/`, alongside the existing
v3 archive.

v5 is a **different DEM**, not just re-picked windows: `2009_pea_hatteras_filled`,
the 2009 base with its LiDAR holes filled from the 2014 NOAA Post-Sandy DEM.
**89 of 90 domains changed interior shape**, all deeper.

* **Drowns at initialisation: 3 per year → 0.** This is the v4 audit's own
  diagnosis acted on — those domains drowned on unsurveyed no-data read as
  water, not on measured water. No setback is moved seaward any more.
* 2004 negatives 0, 1984 negatives still 6 (GIS 10, 11, 12, 84, 85, 86).
* Setbacks moved in 21 (1984) / 26 (2004) of 82 domains, median 0 m. Largest:
  **GIS 79 400 → 500 m** and **GIS 80 440 → 490 m**, two of the three roadways
  the relocation logic acts on.
* **Masks were not re-burned and did not need to be.** They register to the
  `resampled_*.tif` grids, which the fill did not touch — filled and unfilled
  arrays have identical shapes in all 90 domains, the extractor ran with
  `REQUIRE_ROAD_MASKS = True` without erroring, and `HAT_check_geojson_vs_mask.py`
  re-confirms 100 % of centreline vertices inside the mask footprint.
* Fixed in `HAT_road_offset_from_dune_start.py`: the audit's "106 wet cells, 105
  never surveyed" paragraph was **hardcoded** and printed under whatever version
  ran, so on v5 it asserted a v4 measurement in a run where nothing drowns. It
  is now conditional on `n_drowned > 0`, and says where the number came from.
* Also corrected in this README: it claimed the runner spends the **old** method.
  It has spent dune-start since 2026-08-18 (`hatteras_site_config.py:96,110`).

**2026-08-17 (third pass)** — Added `4-compare/` and moved the two cross-method
scripts into it: `HAT_method_comparison_figures.py` (from `1-produce/`) and
`HAT_road_method_diagnostic.py` (from `2-audit/`). `2-audit/` now means "is this
method's own output sound"; `4-compare/` means "which method should we spend".
`HAT_road_method_diagnostic.py` also lost a hardcoded machine-specific
`PROJECT_ROOT`, which made it unrunnable anywhere else.

**2026-08-17 (second pass)** — Consolidated 17 scripts to 10, ~9,900 lines to
~5,900. Everything below is recoverable from git; the blob hash is given where
the file was never committed.

Deleted as superseded by `4-compare/HAT_method_comparison_figures.py` and
`1-produce/HAT_road_placement_on_domains.py`, after porting the two views they lacked (the
2004−1984 movement panel, and the per-profile error band):

* `3-figures/island_wide/HAT_plot_road_on_b3d_domains.py` — `8dd8cfac`
* `3-figures/island_wide/HAT_plot_road_method_on_island.py` — `1310eb42`
* `3-figures/island_wide/HAT_plot_road_placement_accuracy.py` — `e2d308e8`

Merged into `3-figures/HAT_road_domain_views.py` (drawing code ported verbatim;
config, CLI and setback selection are new):

* `3-figures/per_domain/HAT_plot_road_initialization.py` — `e9fd259c`
* `3-figures/per_domain/HAT_plot_barrier3d_grid.py` — `d7c2660b`
* `3-figures/per_domain/HAT_browse_road_domains.py` — `51cf5bb7`

Deleted `provenance/` entirely:

* `HAT_setback_from_lines.py` — `c37ec03e`
* `extractor_flipfalse/HAT_dune_topo_extractor.py` — `f5be4d3b`

> **The claim that folder rested on is false.** It was kept because the legacy
> setback "cannot be regenerated" — `duneline_1984.geojson` and
> `duneline_2004.geojson` really are absent from the repo. But
> `old_method/road_offset_pipeline.py` reproduces **both** shipped
> `RoadSetback_<year>.csv` files **byte-identically**, from ArcGIS point exports
> that are present. The pinned `FLIP = False` extractor copy existed only to
> support that script, and this README's own verification already showed the flip
> does not change setbacks (max difference 0.000000 m).

Also fixed: `road_offset_pipeline.py` had `OUTPUT_DIR = "/data/hatteras_init/…"`
and `HAT_rasterize_road_to_domains.py` had `PROJECT_ROOT = Path("/")` — both
stripped of their prefix, both resolving to the root of the current drive. The
rasterizer died immediately; the pipeline would have written the forcing
somewhere nothing reads. Both now derive their root from `__file__`.

**2026-08-17 (first pass)** — Reorganized into numbered folders. The setback
data folder `processed_offset/` was renamed `old_method_offset/`, which broke
`hatteras_site_config.py` and eight other files until they were repointed. Road
elevation moved out to `../road_elevation/`.
