# road_offset

Where NC-12 sits, per Barrier3D domain — the `road_setback` forcing CASCADE
spends. Two methods are kept, side by side, because the comparison between them
is a result in its own right.

Road *elevation* lives in a sibling folder,
`scripts/input_prep/4-mgmt-forcings/road_elevation/`. It is not an offset
product and carries no year: there is one DEM, so one elevation set serves both
periods.

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
| Drowns at t=0 | 7 (1984), 5 (2004) | **0, both years** |
| vs the rasterized road | +40 m / +23 m median | 0 m *by construction* |

Both are 2 rows × 82 cols, GIS IDs then metres, so switching is a drop-in.

**What the runner spends today** (`scripts/hatteras_site_config.py:96,110`) is
the **dune-start** method, switched 2026-08-18. `3-figures/HAT_road_domain_views.py`
defaults to `--method dunestart` to match; change both together, or the figures
stop describing the model you run.

Current outputs are built on **`2009_v5`** — the gap-filled DEM. See
`hat_topo_version.py`: the version is resolved once, from the extractor, so
bumping `VERSION` there moves this whole tree with it.

## The workflow

```
# old method — one run per vintage
HAT_ROAD_YEAR=1984 python 1-produce/old_method/road_offset_pipeline.py
HAT_ROAD_YEAR=2004 python 1-produce/old_method/road_offset_pipeline.py

# dune-start method — masks first, then the measurement
HAT_ROAD_YEAR=1984 python 1-produce/HAT_rasterize_road_to_domains.py
HAT_ROAD_YEAR=2004 python 1-produce/HAT_rasterize_road_to_domains.py
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
