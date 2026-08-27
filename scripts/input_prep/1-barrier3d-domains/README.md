# 1-barrier3d-domains - dune heights and interior topography

Takes the `domain_<N>.npy` elevation arrays that `0-elevation` produces and
turns them into the two things Barrier3D initialises from:

    topography/  domain_<N>_topography.npy   (200, 50)  dam
    dunes/       domain_<N>_dune.npy         (50,)      dam above berm
    (plus        domain_<N>_nodata.npy       coverage masks)

Input is metres NAVD88; output is Barrier3D's native decameters. That
conversion happens here and nowhere downstream.

## The stage is one script

```
HAT_dune_topo_extractor.py     the whole chain: pick -> extract -> figures
old/                           ancestors, kept for provenance only
```

Run it from anywhere. `MODE` selects `"pick"` (drag a cross-shore dune search
window per domain, saved to JSON after each one, safe to quit and resume),
`"run"` (extract using the saved windows), or `"pick_and_run"`.

## Which product it writes, and how everything else finds it

`TOPO_PRODUCT` and `VERSION` near the top of the extractor decide what gets
written. **Nothing else hardcodes that path.** `scripts/hat_topo_version.py`
parses those two names straight out of this file and every reader - the
runner, the groin sweep, the road tree, the figure scripts - resolves through
it:

    from hat_topo_version import topo_dirs, array_path, array_name

Version numbers restart at `v1` **per product**: `1984-start/v1` and
`2004-start/v1` are different surfaces from different DEMs. See
`data/hatteras_init/1-barrier3d-domains/LINEAGE.md`.

Because the resolver locates this file **by path**, moving or renaming it
breaks version resolution *silently* - `_extractor_state()` returns
`(None, None)` on a miss and resolution falls through to the `CURRENT` file.
If you move it, update `EXTRACTOR` in `hat_topo_version.py` and the
`parents[2]` in the `sys.path` line here, then check:

    cd scripts && python -c "import hat_topo_version as h; print(h.EXTRACTOR.is_file(), h._extractor_state())"

## Picks are frame-dependent - the one real footgun

A window picked on a *straightened* array is a perfectly valid index range on
an unstraightened one; it just points at different cells. `save_windows`
records `STRAIGHTEN` and the run pass refuses on a mismatch, but the pick
file name has to be set by hand to match. `STRAIGHTEN = True` writes
`{DEM_NAME}_straight`, `False` writes `DEM_NAME`.

`STRAIGHTEN = False` is also **how you get an oblique-uncorrected control** -
there is no separate script for it, and the unstraightened control windows
live in `data/.../1-barrier3d-domains/control-picks/`.

## The road overlay is display only

NC-12 is drawn on the picker so a dune-crest argmax cannot be dragged onto the
road embankment. It does not enter `find_dunes`, `build_interior` or
`straighten_profiles`: outputs are byte-identical with `SHOW_ROAD` either way,
only the figures and the road columns of the settings sheet change. The masks
it draws come from `4-mgmt-forcings/road_offset/`, which means that stage runs
*before* a pick pass, not after.

## old/

| file | what it was |
|---|---|
| `dune_topo_extractor_v3.py` | Lexi's original, the direct ancestor |
| `dune_topo_extractor_from_GIS.py` | the first Hatteras adaptation; source of the `-1.0` clamp and the fixed 8-px window the current file still explains |

Neither runs against the current tree - their paths (`hatteras_init/dunes/`,
`hatteras_init/topography/`, `hatteras_init/elevations/`) predate the numbered
reorganisation and no longer exist. They are here to explain the live file's
comments, not to be executed.

## What was removed 2026-08-26, and why

All recoverable at `7cd5af0` via `git show 7cd5af0:<path>`.

| removed | why |
|---|---|
| `topography_dunes/no_oblique_correction/HAT_dune_topo_extractor.py` | a fork of this script with the straightening deleted. `STRAIGHTEN = False` does the same job, and the fork predated the road overlay and per-product picks. Same filename and same docstring as the live file - the exact stale-read hazard `hat_topo_version.py` exists to prevent. |
| `gis-export-npy.py` | ArcGIS-era `.npy` exporter. Did not parse (`root_foimpolder`, a misindented `if`). Superseded by `0-elevation/2-produce/HAT_export_to_numpy.py`. Its `nodata_to_value=-10` convention is restated in `4-mgmt-forcings/road_offset/1-produce/HAT_rasterize_road_to_domains.py`, which is the only thing that still depended on it. |
| `export-rasters-to-npy.py` | another author's Ocracoke `Topo_2019` ArcGIS loop, hardcoded to a `C:\Users\frank\...` OneDrive path. Not Hatteras. |
| `buffers/buffer_creation.py` | one-shot copier for the three `buffer/` arrays. All three sources are gone and its destination was a root-absolute `/data/...`, so it could not run. Its provenance is now written up in `data/hatteras_init/1-barrier3d-domains/buffer/README.md`. |
| `topography_dunes/old/dune_topo_extractor_from_cascade.py` | seeded the next period from a finished run's `.npz` instead of a DEM. Dead as written - it targets `HAT_hindcast_1984_2024_updated.py` Section 8 and `TOPO_DUNE_INIT_YEAR`, neither of which exists. **Worth reviving for the `forecast/` product**, which otherwise needs a 2025 DEM; recover it before building that. |
