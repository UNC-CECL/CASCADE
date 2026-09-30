# 1-barrier3d-domains - dune heights and interior topography

Takes the `domain_<N>.npy` elevation arrays that `0-elevation` produces and
turns them into the two things Barrier3D initialises from:

    topography/  domain_<N>_topography.npy   (200, 50)  dam
    dunes/       domain_<N>_dune.npy         (50,)      dam above berm
    (plus        domain_<N>_nodata.npy       coverage masks)

Input is metres NAVD88; output is Barrier3D's native decameters. That
conversion happens here and nowhere downstream.

## The stage in two halves

```
1-extraction/                      DEM -> the domain arrays (v1, v2): the first half of the stage
    HAT_dune_topo_extractor.py     the whole chain: pick -> extract -> figures
    nodata_audit/                  the dropout bridge
    old_extractors/                ancestors, kept for provenance only
2-domain-reconstruction-1984/      v2 -> v3: the 1984 domains reconstructed from the 1996-based DEM,
                                   one subfolder per step, in the order the argument runs (2026-09-09)
    1-measurement/                 how far the dune line moved: the shift N comes from, and its per-domain plotters
    2-extent/                      how many rows, which domains: the footprint table, the scope report
    3-placement/                   WHERE the rows go: the road placement check; imagery-review/ (batch, window, summary)
    4-fill/                        what the rows contain: the copy fill, the fill plotters, the explainers
    5-build/                       HAT_build_footprint_version.py (-> dune-topo/v3)
    6-result/                      the hindcast on v3 against v2, and the version comparisons
```

**Two halves, one stage** (2026-09-09, Hannah). Both halves hand the runner the
same thing, a dune-topo version under `1984-start/`, so the reconstruction is
not a pipeline stage of its own: its data live here, its inputs include the
stage-4 road offset measured on v2, and once v3 exists that offset is measured
again on it. A stage between 1 and 2 would claim a linear order the method
does not have.

The footprint steps mirror `data/.../1984-start/2-domain-reconstruction-1984/<step>/` and
its `figures/<step>/`; paths are resolved through `hat_topo_version`
(`insert_scope_step`, `insert_figures_dir`, `duneline_shift_dir`), never built
by hand. Every script finds the repo root by walking up from its own file, so
the depth does not matter to it.

Run it from anywhere. `MODE` selects `"pick"` (drag a cross-shore dune search
window per domain, saved to JSON after each one, safe to quit and resume),
`"run"` (extract using the saved windows), or `"pick_and_run"`.

## Which product it writes, and how everything else finds it

`TOPO_PRODUCT` and `VERSION` near the top of the extractor decide what gets
written. **Nothing else hardcodes that path.** `scripts/site_layer/hat_topo_version.py`
parses those two names straight out of this file and every reader - the
runner, the groin sweep, the road tree, the figure scripts - resolves through
it:

    from site_layer.hat_topo_version import topo_dirs, array_path, array_name

Version numbers restart at `v1` **per product**: `1984-start/v1` and
`2004-start/v1` are different surfaces from different DEMs. See
`data/hatteras_init/1-barrier3d-domains/LINEAGE.md`.

Because the resolver locates this file **by path**, moving or renaming it
breaks version resolution *silently* - `_extractor_state()` returns
`(None, None)` on a miss and resolution falls through to the `CURRENT` file.
If you move it, update `EXTRACTOR` in `hat_topo_version.py` and the
`parents[2]` in the `sys.path` line here, then check:

    cd scripts && python -c "from site_layer import hat_topo_version as h; print(h.EXTRACTOR.is_file(), h._extractor_state())"

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

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### 1-extraction/HAT_dune_topo_extractor.py

Build each domain's Barrier3D dune and interior arrays from its DEM, with a hand-picked dune search window.

Notes that were in the code:

```text
HAT_dune_topo_extractor.py

Hatteras CASCADE dune & interior topography extractor with per-domain,
interactively selected dune search windows.

INPUT  : domain_#.npy, shape = (alongshore_rows, cross_shore_cols),
elevation in m NAVD88
OUTPUT : interior topography and dune height arrays, in decameters (dam)
+ a JSON of the per-domain dune search windows (re-runnable)

WORKFLOW
1. MODE = "pick"          -> step through domains, drag a dune search
window on the profile stack, saves to JSON
after every domain (safe to quit and resume)
2. MODE = "run"           -> extract using the saved JSON windows
3. MODE = "pick_and_run"  -> both in one pass

PICKER LAYOUT / KEYS
The picker draws the domain with the OCEAN AT THE BOTTOM: cross-shore runs
vertically (cell 0 = ocean, landward upward), alongshore runs left-right.
Left panel = elevation map, right panel = profile stack, shared y-axis.
This is display only; i0/i1 are still cross-shore indices from the ocean.

NC-12 is drawn on both panels (v4): filled road cells and a dashed centre
line on the map, a shaded cross-shore envelope on the profile stack, one
colour per road vintage. It is there to stop the window being dragged onto
the road embankment, which a dune-crest argmax will happily lock onto. The
road constrains nothing in code -- see ROAD OVERLAY in CONFIG.

drag vertically on either panel : set the cross-shore search window
enter / close window            : accept current window
r                               : reset to the default window
s                               : skip this domain (use DEFAULT_WINDOW_PX)
q / esc                         : quit picking, keep everything saved so far
```

```text
needed for a real, blocking, interactive picker window. Only force it if
tkinter actually exists, otherwise matplotlib fails later with a confusing
error at figure-creation time instead of here.
```

```text
Anchored 2026-09-14: absolute into a home directory, or into a tree
renamed since. Rule 5 of ORGANIZATION.md.
```

```text
Sized for a projected slide rather than a screen: at 19 in wide, 11 pt tick
labels are unreadable once the figure is scaled into a talk.
```

```text
--- MODE ---------------------------------------------------------------
v4 RE-PICKS every window with the road drawn on the picker (see ROAD OVERLAY).
The v3 windows were drawn without knowing where NC-12 sits, so a window could
sit landward of the road and take the road embankment for a dune crest. This
pass is an ADJUSTMENT, not a blind redraw: the picker opens on each domain's
saved v4 window (seeded from v3), so "r" resets to the v3 pick and accepting
without dragging keeps it.
```

```text
Domains to process. This filters BOTH passes -- picking AND running -- so it
is the modelled set, not just a picking subset.

list(range(1, 91))      the 90 domains CASCADE runs (D1 = Cape Point ->
D90 = Pea Island). The DEM folder holds 131; 91-131
are north of the study area and are not modelled,
so picking or extracting them is wasted work.
None                    every domain_*.npy found
[1, 5, 11, 33, 67, 74]  a trial set. Worth doing before committing to 90
hand-picks: 33 (25.6 deg) and 67 (22.8 deg) are the
worst obliquity on the island and two of the four
SCATTER flags; 11 (8.0 deg) is the domain the road
setbacks kept failing on; 74 (5.7 deg) is a
near-square control where panels 1 and 2 should look
almost identical; 1 and 5 are Cape Point, where a
LINEAR fit to the shoreline is most likely to fail.
```

```text
True for the v4 road re-pick. The v4 picks file is SEEDED from v3, so every
domain is already present and False would skip all 90 -- nothing would be
re-picked. Set back to False once the re-pick is finished, so a later run in
"pick_and_run" resumes instead of starting over.
```

```text
--- PATHS --------------------------------------------------------------
Everything for one settings variant lands in ONE run folder, so comparing
versions means comparing two directories:

data\hatteras_init\dune_topo\
picks\
HAT_dune_search_windows_2009_pea_hatteras.json  <- your picks (see below)
2009_v1\
RUN_MANIFEST.txt                     <- every setting that made this folder
HAT_dune_topo_settings_2009_v1.xlsx  <- per-domain sheet (+ .csv)
HAT_dune_topo_summary_2009_v1.png    <- all domains on one page
topography\  domain_7_topography_2009.npy   <- CASCADE reads these two
dunes\       domain_7_dune_2009.npy
figures\
gis_vs_processed\  domain_007_gis_vs_processed.png
qc\                domain_007_qc.png
2009_v2\   ... same shape, nothing shared, nothing overwritten

NOTE: topography\ and dunes\ moved. Point the hindcast runner at
RUN_DIR\topography and RUN_DIR\dunes for whichever version you're using.
v3 = v2 settings + ALONGSHORE_FLIP = True (see GEOMETRY below). Bumped rather
than reused so 2009_v2 survives as the unflipped reference to diff against.
v4 = v3 settings + the NC-12 overlay and a full re-pick of the dune windows
with the road visible. No processing setting changed, so a v4 run in "run"
mode on the v3 windows reproduces the v3 arrays byte for byte -- what makes
v4 different is the WINDOWS, which is exactly why it gets its own folder and
its own picks file rather than overwriting v3.
WHICH PRODUCT this run builds. One of the period folders under
data/hatteras_init/1-barrier3d-domains/ - see scripts/site_layer/hat_topo_version.py.
Added 2026-08-25 when the tree went period-first; before that there was only
one topography and both hindcast periods read it.

"1984-start"   from DEM 2009-2014-1996  (1996 ALACE, no road boundary)
"2004-start"   from DEM 2009-2014       (the baseline gap fill)
"forecast"     from a 2025 DEM, later
```

```text
scripts/ on the path, so the ARRAY NAMES come from the same resolver that
owns the directory layout - one definition, used by the writer here and by
every reader. hat_topo_version PARSES this file for TOPO_PRODUCT/VERSION and
never imports it, so importing it here creates no cycle.
```

```text
renumber). What this script WRITES.
Since 2026-09-04 it no longer decides what is
READ: dune-topo/CURRENT outranks it in
hat_topo_version.resolve_version. Do not point
this at v3-v8 (layers built ON v2) - a re-run
would overwrite their arrays. Guide: 1984-start/
2-domain-reconstruction-1984/DUNE_TOPO_VERSION_GUIDE.md.
To change the default, edit CURRENT; for one
run, HAT_TOPO_VERSION_1984_START.
```

```text
A LABEL ONLY. It is no longer written into any filename.

The arrays were domain_<N>_topography_2009.npy until 2026-08-26. The year was
false for both live products - 2004-start is the 2009+2014 mosaic, 1984-start
is 2009+2014+1996 - and the fix is not a better year but no year: the period
lives in the PRODUCT DIRECTORY, which every reader must resolve anyway. See
the long note in scripts/site_layer/hat_topo_version.py for why a per-period tag was
tried and reverted the same day.

What remains here is the text on two figures. It names the product, because
"2009 extracted crest" was captioning a 1996-grafted surface.
```

```text
v5 reads the GAP-FILLED arrays: 2009 base with its gaps filled from the 2014
NOAA Post-Sandy DEM (see data/hatteras_init/0-elevation/FIGURES.md).
25,591,292 cells filled at 1 m, 254,760 at this 10 m grid, mostly landward of
NC-12. The un-filled set is still at "2009_pea_hatteras" with its v4 picks
intact, so both can be run and compared.

PICKS ARE PER-VERSION and none exist for v5 yet - picks/ holds only
2009_v4, 2009_pea_hatteras and ..._straight. This run needs a pick pass
(MODE includes "pick"). Re-picking is warranted, not a formality: the filled
interior differs materially from what v4 was picked against.
```

```text
The product folder already says which run this is, so the run folder is the
VERSION alone: 2004-start/dune-topo/v5, not .../2009_v5. The migrated v5 keeps
its inner filenames (HAT_dune_topo_settings_2009_v5.csv) - those are outputs
of a run that happened, and renaming them would rewrite history.
```

```text
These two moved when data\hatteras_init was reorganized into the numbered
1-barrier3d-domains \ 2-brie-offset \ ... tree. They previously read
INIT_ROOT/"elevations"/DEM_NAME and INIT_ROOT/"dune_topo", neither of which
exists any more, so the script could not find its own inputs or the 2009_v2
run folder it wrote.
Period-first layout (2026-08-25). Was {DEM_YEAR}-raw/{DEM_YEAR}-npy-arrays/
{DEM_NAME} and {DEM_YEAR}-dune-topo/{RUN_NAME}, which keyed three path
segments on the DEM year and said nothing about which run the arrays were
for.
```

```text
Picks live OUTSIDE the run folder. They are the only artifact here you cannot
regenerate, and they describe where the dune sits in the DEM -- not which
settings variant you're testing -- so v2, v3... reuse them by default. Set
PICK_SET = RUN_NAME instead if you want a version to carry its own picks.
The picks are FRAME-DEPENDENT. A window picked on a straightened array is a
valid index range on an unstraightened one -- it just points at different
cells -- so the two frames cannot share a file. save_windows() stamps
"straightened" on every entry and the run pass refuses on a mismatch, but it
writes back to WINDOW_JSON on every domain: pointing this at the v1 set would
overwrite v1's unstraightened windows as you re-pick, and 2009_v1 would stop
being reproducible.
Set this by hand to match STRAIGHTEN below (it is defined further down, so it
cannot be referenced here):
STRAIGHTEN = True   ->  f"{DEM_NAME}_straight"
STRAIGHTEN = False  ->  DEM_NAME

v4 CARRIES ITS OWN PICKS -- PICK_SET = RUN_NAME, not the shared straight set.
This is the case the warning above describes. v4 re-picks all 90 windows, and
save_windows() writes back to WINDOW_JSON after EVERY domain, so pointing this
at f"{DEM_NAME}_straight" would destroy the v3 picks as you worked and 2009_v3
would stop being reproducible. The v4 file was seeded by copying the v3 one, so
every domain opens on its v3 window and an unchanged domain stays unchanged:

picks\HAT_dune_search_windows_2009_pea_hatteras_straight.json   v1-v3, FROZEN
picks\HAT_dune_search_windows_2009_v4.json                      v4, re-picked
```

```text
No TAG. Names come from hat_topo_version.array_name(), which is also what
every reader calls - see _gis_id() at the save site.
```

```text
--- ISLAND OFFSETS -----------------------------------------------------
Measured per-domain dune offsets used to place domains in a common cross-shore
frame. One value per GIS domain, header row = year.

CONVENTION (both inferred from the data -- change if your pipeline says otherwise):
OFFSET_ROW_ORDER = "D1_first"    row 0 of the CSV is domain 1 (Cape Point).
Check: this puts the largest raw_offset (6301 m) at Cape Point, which is what
a seaward-protruding headland should look like. Reversed puts the max at
Pea Island instead.
OFFSET_SEAWARD_POSITIVE = True   larger value = further seaward.
Check: 75/90 domains go negative 1984->2004, mean -2.07 m/yr. Seaward-
positive reads that as island-wide retreat at ~2 m/yr (right for Hatteras);
the other sign reads it as island-wide accretion (wrong).
```

```text
Each start's CURRENT build (2026-09-18). These named hindcast_<year>/
folders that no longer existed, and the fallback search found several
candidates per year, called it ambiguous and skipped the offsets.
```

```text
WHICH YEAR'S OFFSETS THIS PRODUCT IS PLOTTED AT (2026-08-27).

Both files above are still LOADED -- panel 3 of the island-offsets figure is
a 1984->2004 change rate and needs both. What PRODUCT_YEAR controls is the
PLAN VIEW, which used to be written once per offset year and so produced
four PNGs per run:

2004-start/dune-topo/v1/..._planview_v1_1984_padded.png
..._planview_v1_1984_trimmed.png
..._planview_v1_2004_padded.png
..._planview_v1_2004_trimmed.png

The topography is IDENTICAL in all four -- it is whatever TOPO_PRODUCT built.
Only off_cells differs, i.e. which canvas row each domain's row 0 lands on.
So the _1984_ pair above places the 2009+2014 surface at the 1984 measured
shoreline: an island that never existed, sitting in the 2004 product's folder
under a filename that reads like a 1984 initial condition.

That loop is older than the period-first tree. When ONE topography was read
by BOTH hindcast periods (see hat_topo_version.py), plotting it at both years
was the whole point. Now that 1984-start and 2004-start are separate DEM
products, each with its own offsets -- 1984-start/README.md line 149 pins
shoreline_offset to 2-brie-offset/1984/ -- the pairing is 1:1 and
the cross-product figures are noise.

The year is RESOLVED from the product through hat_topo_version.YEAR_PRODUCT,
never spelled here, so a third product cannot pick up a stale literal.
strict=False: "forecast" and "buffer" are not hindcast periods, and for those
PRODUCT_YEAR is None and the plan view falls back to every year it loaded --
the old behaviour, for the case where there is no relevant year to pick.
```

```text
Plan-view canvas, reproducing the ABSOLUTE placement of
initialization_figures.py (island_<year>_absolute.png) exactly:
offset_cells = round(offset_m / 10); each domain's topo row 0 (ocean side)
lands on canvas row = offset_cells; alongshore flipped with np.fliplr.
```

```text
ISLAND_FLIP_ALONGSHORE was REMOVED (2026-08-17). It applied np.fliplr to each
domain INSIDE the placement loop, so it did not mirror the island -- it reversed
the 50 cells of every 500 m block against the ascending domain order. That was a
workaround for the source arrays having the within-domain alongshore order
backwards, which ALONGSHORE_FLIP now fixes at load. Keeping both double-flipped:
measured seam/inner discontinuity ratio on the plotted canvas was

v2 arrays + flip  1.97 (right)      v2 arrays, no flip  21.15
v3 arrays + flip 21.15 (WRONG)      v3 arrays, no flip   1.97 (right)

so v2's plan view was only correct because two errors cancelled. The flip is
gone rather than defaulted False, per the decision to keep this code lean.
CONSEQUENCE: plotting a PRE-CORRECTION run (2009_v2 or earlier) through this
function now sawtooths, and _assert_alongshore_continuity below will say so.
To regenerate a v2 plan view, re-extract it with ALONGSHORE_FLIP = True instead.
```

```text
cells clip to the dark navy bottom of terrain,
so the model's cross-shore extent stays visible
and only outside-canvas NaN is light blue.
True = sentinel also renders light blue.
```

```text
Display only, does not touch the .npy files. ONE FIGURE PER YEAR PER MODE, so a
single run gives you both versions to compare:
"trimmed" = each domain only as tall as its own island, as stored
"padded"  = every domain given the same cross-shore extent (ISLAND_PAD_ROWS),
sentinel-filled landward, so the water behind each domain shows as
the navy wedge in the poster figure
If TRIM_INTERIOR_ROWS = False the arrays are already padded and the two match.
```

```text
200 = 2000 m, equal to TOPO_ROWS, so nothing is
cropped (the poster look)
100 = 1000 m, tighter; the script warns if that
crops real land off the sound side
```

```text
--- ISLAND SECTIONS ----------------------------------------------------
D1 = Cape Point (south) -> D90 = Pea Island (north). Labels the sheet, the
per-domain figures and the summary figure. Set to [] to disable.
```

```text
-3.0 keeps back-barrier marsh cells (Lexi's v3 edit)
-1.0 was the original dune_topo_extractor_from_GIS behavior
```

```text
--- NO DATA IS NOT WATER -----------------------------------------------
The source LiDAR carries -10.0 m NAVD88 where it has no return. Clamping used
to fold that into SENTINEL_WATER_M, so a cell the survey never saw became
indistinguishable from a cell measured below MHW. That conflation is not
cosmetic: roadway_manager.bulldoze drowns a roadway when >20% of the cells
BORDERING it sit at or below 0 m MHW, and a no-data cell satisfies that test.
In GIS 78/79/80 the row landward of NC-12 is 17-25 no-data cells and ZERO
genuinely wet ones, so all three roadways "width-drowned" at t=0 on the
strength of missing survey coverage. The giveaway is where the data stops:
those profiles end while the ground is still 0.5-0.7 m ABOVE MHW, whereas
their neighbours grade down through zero the way a real sound margin does.

So no-data is tracked separately and written to its own array. The topography
CASCADE reads is UNCHANGED -- no-data still lands on SENTINEL_WATER_M there,
because Barrier3D has no representation for "unknown" and inventing an
elevation would be worse. What changes is that the information now survives,
so any consumer that cares can ask.

<stem>_nodata.npy         bool, same shape as the topography array
True = this cell was never surveyed
```

```text
True reverses alongshore order after orienting. REQUIRED for Hatteras, because
the GIS row order and the domain numbering run in OPPOSITE directions:

1. Every resampled_domain_*.tfw has pixel Y size = -10 and zero rotation, so
the rasters are north-up and array ROW 0 IS THE NORTHERNMOST row. With
OCEAN_LOC="right" no rotation is applied, so axis 0 stays alongshore:
within a domain, index 0 = north, index 49 = south.
2. The .tfw upper-left northing increases monotonically with domain number
(D1 = 3,899,274 -> D90 = 3,944,002 m, 502.6 m per step), so DOMAIN NUMBER
INCREASES NORTHWARD. D1 = Cape Point / south, as OFFSET_ROW_ORDER and
SECTIONS already assume.

Unflipped, the assembled island is a 90-tooth sawtooth: each 500 m block is
internally mirrored against its neighbours. Measured on the 2009_v2 output,
mean jump at domain seams vs. mean jump within a domain:

as saved   flipped
island width        21.1x      1.97x     (134 m seam jumps -> 12.5 m)
dune height          3.45x     1.54x
mean interior elev   4.87x     1.41x

WHAT THIS DOES AND DOES NOT AFFECT. BRIE resolves ONE node per 500 m Barrier3D
domain (brie_coupler.py: dy=500, alongshore_section_count=ny) and exchanges a
scalar x_s, so it never sees within-domain cells -- the sawtooth does not
corrupt alongshore transport. Barrier3D's own 50-cell axis is near
mirror-symmetric (Q1/Q3 weighting and the i>0 / i<BarrierLength-1 boundaries
are symmetric; the router writes only into row d+1 so sweep order is
immaterial). The one asymmetry is a bug: DiffuseDunes loops
`range(2, BarrierLength)`, so alongshore cell 0 never exchanges sand with cell
1 while cell 49 does.

So this flag is about GEOGRAPHIC FIDELITY AND CROSS-INPUT CONSISTENCY, which is
where the real exposure is: NC-12 masks, community/nourishment zones and
setback CSVs are all built in GIS row order and must share ONE frame with the
topography, or the road sits at the mirrored alongshore position inside every
domain. Flipping HERE -- at load, before the shear -- is what keeps the picker,
the road masks, shear_like(), the QC figures and the saved .npy in that one
frame. Do NOT reverse only the arrays written to disk: that is what
PEA_dune_topo_extractor_..._alongshore_corrected.py does
(CORRECT_SAVED_ALONGSHORE), and it leaves every figure and mask mirrored
against the files CASCADE reads.

THE SAVED PICKS ARE UNAFFECTED. i0/i1 are scalar CROSS-SHORE indices, and the
flip commutes with straighten + water-trim: reversing alongshore negates the
polyfit slope but leaves ref[i] -> ref[n-1-i], hence shear[i] -> shear[n-1-i]
and an identical per-profile cross-shore frame. Verified on all 90 domains
(z_flipped == z_unflipped[::-1] exactly, c0 and n_cross_trimmed unchanged) by
HAT_alongshore_frame_check.py. No re-picking is needed.
```

```text
--- ROAD OVERLAY -------------------------------------------------------
NC-12 drawn on the picker and every per-domain figure. DISPLAY AND DIAGNOSTICS
ONLY: the road never enters find_dunes, build_interior or straighten_profiles,
so the arrays CASCADE reads are byte-identical with SHOW_ROAD either way. What
it changes is what YOU see while picking -- a search window that sits landward
of the road is picking the road embankment, not a dune crest.

THE MASKS ARE NOT MADE HERE. HAT_rasterize_road_to_domains.py burns the road
geojson onto each domain's resampled_*.tif affine, so they are cell-for-cell
aligned with the DEM .npy by construction. This script only reads them, and
refuses anything whose shape disagrees.

TWO VINTAGES ON A 2009 DEM. The road lines are 1978 and 2008 exports -- the
stand-ins for the 1984 and 2004 starts, paired in
hat_topo_version.ROAD_LINE_FOR_YEAR -- and the topography is 2009. That
mismatch IS the subject of RoadOffset_dunestart_audit.md, so both are drawn,
in different colours, with the LINE vintage in every label -- never read the
1978 line as 1984 topography, nor as 1978 topography.

RENAMED 2026-09-15. The masks were domain_<N>_road_1978.npy / _2004.npy and
now carry the line's true vintage. ROAD_YEARS are therefore LINE vintages,
not period starts, and they key every road_masks / road_stats dict below and
every "road ... <year>" column of the settings sheet.
```

```text
All 90 domains have masks for both years, so this
only ever fires if the raster tree moved or a
re-export changed a grid.
```

```text
D1-D7 (Cape Point) have ZERO road cells in both vintages -- NC-12 does not
reach the point. An empty mask is normal and silent; only a missing file or a
shape mismatch is an error.
```

```text
so a picked window can't wander onto the wet beach
where the shoreline curves within a domain
```

```text
dune in the domain -> alongshore alignment preserved.
False = each profile starts behind its own dune (v3's
original behavior; breaks alongshore alignment).
```

```text
True  = Lexi's v3: drop all-water rows, so each domain's
interior array is only as tall as its own island.
False = your dune_topo_extractor_from_GIS.py: every domain
is exactly (TOPO_ROWS, ALONG_COLS), padded landward
with sentinel. This is what produced the 2009_v1
arrays behind the poster figure.
```

```text
--- STRAIGHTEN ---------------------------------------------------------
Shear each alongshore profile so the shoreline runs HORIZONTALLY, before the
dune window is picked.

Why: the clip boxes are north-up (rot = 0.00 on all 131 rasters) while
Hatteras trends NNW, so cross-shore is due east-west and the shoreline crosses
each 500 m domain diagonally. Domain 11's dune runs cell 3 -> 13: 8 deg of
obliquity, 100 m of drift. Two consequences, and they are separate:

1. THE WINDOW has to be wide enough to span the diagonal. Domain 11's is
[3, 15] = 120 m -- also wide enough to catch a back-dune, a wooded ridge,
or a house. The two worst-obliquity domains (GIS 33 at 25.6 deg, GIS 67
at 22.8 deg) are two of the four SCATTER flags in the road setbacks.

2. THE WEDGE. USE_CONST_INTERIOR cuts horizontally at max(dune_loc) + 1,
throwing away everything seaward of that on every other profile.

Straightening fixes (1). USE_CONST_INTERIOR = False fixes (2). Set both --
even straightened, residual dune variability is real and the const cut still
costs 20-40 m of it.

NEITHER fixes the distance inflation: the profiles are still due east-west, so
cross-shore AND alongshore distances stay long by 1/cos(theta) -- 1% at 8 deg,
4% at 16 deg, 11% at 26 deg. A "500 m" domain spans 500/cos(theta) m of
shoreline. Only re-clipping with rotated boxes fixes that.

THE PICKS BECOME FRAME-DEPENDENT. A window picked straightened is a valid
index range on an unstraightened array; it just points at different cells.
save_windows records STRAIGHTEN and the run pass refuses on a mismatch. Use a
NEW WINDOW_JSON and a NEW RUN_NAME rather than overwriting a picked set.
```

```text
What to align on. "beach" = the first cell above BEACH_START_THR_M, which is
computed without a window -- no chicken-and-egg with the dune pick.
```

```text
"linear"  fit a straight line to start_beach, shear by that. Over 500 m the
island does not curve, so the obliquity IS linear: the fit removes
exactly the diagonal and leaves real alongshore variability in the
array. Immune to a few bad profiles.
"raw"     shear by start_beach itself. Flattens real structure too and folds
every noisy pixel into the geometry. Diagnostic only.
```

```text
Cells shifted in from beyond the original array were never surveyed
either, so they are no-data rather than water. They trim identically
(see water_col_bounds), so this changes the mask, not the topography.
```

```text
shear_like fills with np.zeros_like, so on a bool array the cells shifted
in from beyond the seaward end become False -- correct for a mask, which
is why the bool dtype has to survive this call.
```

```text
USE_CONST_INTERIOR: the cut is horizontal, so row 0 is the same cell on
every profile whether or not that profile found a dune.
```

```text
ROAD MASKS

Read-only consumers of HAT_rasterize_road_to_domains.py. Nothing here writes a
mask, and nothing here influences the dune search or the saved arrays.
```

```text
The mask is checked against the RAW DEM shape, before any orienting, so
a grid mismatch is reported against the thing the rasterizer actually
snapped to. No transpose or resize: that would hide a misregistration.
```

```text
Cap at ALONG_COLS, as the setback script does: profiles beyond the 50
CASCADE keeps are not part of the measurement.
```

```text
counted on the SEAWARD edge, because that is the value that gets
floored before it reaches the model
```

```text
sanity check that OCEAN_LOC is actually right (the original script's AUTO_ORIENT,
demoted to a warning so orientation stays an explicit, documented choice)
```

```text
No-data is identified on the RAW array, before the clamp folds it into the
water sentinel. It then rides the same shear/trim/slice path as z, marked
by a value far below the clamp, and is separated out again at save time.
```

```text
ORDER MATTERS: start_beach -> straighten -> water trim.
start_beach is found on the untrimmed array because the shear is what
defines the frame; c0 must then be measured on the array the window is
actually picked in, or the two disagree.
```

```text
LAST, because align_mask_to_topography needs the finished z, c0 and shear.
Loading the road cannot change any of them -- if it ever appears to, the
shape assert inside align_mask_to_topography is what will say so.
```

```text
--- SUGGESTED WINDOW ---------------------------------------------------
How far landward of the beach a foredune is allowed to be looked for, and how
much margin to leave either side of the crest the search finds.
```

```text
The road's seaward edge, so the suggestion can never propose a window
that reaches NC-12. road_masks is the aligned dict the caller passes
for display; if it is absent the hunt is simply unbounded landward.
```

```text
HALF-CELL OFFSETS, so the band covers the cells actually SEARCHED.

The window is half-open: find_dunes slices prof_arr[i, i0:i1], so
cell i1 is NOT searched. Drawn as axhspan(i0, i1) the band's edge
landed on the CENTRE of cell i1, so "drag to cover the crest"
excluded the crest -- and did so silently, since the argmax simply
pinned at i1-1. It cost the real crest at GIS 43, 64, 72, 85 and
86, up to 1.5 m at GIS 85 (3.44 m picked against 4.93 m actual).

Display only: no stored pick moves, no number changes. It just
makes the shaded band mean what a reader assumes it means.
```

```text
the suggestion, as an outline so it reads as a proposal rather
than a second selection
```

```text
argmax INSIDE the window. v3 used np.where(prof == dune_elev)[0][0] over
the whole profile, which can snap the dune onto an earlier cell of equal
elevation (common on a quantized DEM).
```

```text
Where SAVED interior row 0 sits on each profile, and the road measured
against it. Diagnostics only -- computed from dune_loc and the same
build_interior call the arrays came from, and written to nothing but the
settings sheet and the figures.
```

```text
Split no-data back out. The topography written is byte-identical to what
this script produced before the mask existed: every no-data cell goes back
to SENTINEL_WATER_M, because Barrier3D has no representation for
"unknown". The mask is what carries the distinction forward.
```

```text
the saved arrays are in the straightened frame; nothing about a .npy
says so, and a window or a mask from the other frame is silently wrong
```

```text
Shared elevation span so the three panels are comparable by eye -- just
labelled NAVD88 vs MHW. Scaling each to its own percentile lets the -10 m
shoreface swamp the ramp and the island reads as one flat colour.
```

```text
The road in the RAW frame: unsheared and untrimmed, so this is NC-12's real
diagonal across the north-up clip box. Compare it with the same road on
panel 2 -- that difference is what the shear removes, and it is the error
any raw cross-shore median of the road inherits.
```

```text
crop to the island so the diagonal is legible next to the straightened
panel; the full raw is mostly sound and shoreface
```

```text
The road re-indexed into the SAVED interior grid. This is the panel that
answers the question the other two cannot: is NC-12 inside the array
CASCADE reads, and at which interior row -- the same row
roadway_manager.bulldoze lands on with int(road_setback / dy). A road that
falls off the seaward edge here has a negative setback.
```

```text
the dune strip is a tenth the height of the map above it, so the map's
aspect would leave it a bar too short to read a tick off
```

```text
Road overlay. Recorded because the road columns in the sheet are
meaningless without knowing which vintage and which mask tree produced
them -- but note that NONE of these affect the saved arrays.
```

```text
centre-referenced, for continuity with the roya-style dune-to-road
number; NOT the value bulldoze indexes
```

```text
Profiles where the road is SEAWARD of interior row 0. Barrier3D cannot
represent that (int(negative/dy) indexes from the landward end), so it
is the flag the setback audit's NEGATIVE floor exists to handle.
```

```text
Union of keys in first-seen order, not rows[0].keys(). With
REQUIRE_ROAD_MASKS = False a domain whose mask is missing contributes no
road columns, and DictWriter raises on any row holding a key the header
does not -- so keying off the first row would turn one absent mask into a
crash after all 90 domains had been processed.
```

```text
NC-12's median cross-shore cell per domain, on the same axis as the search
window. Where the road line dips INTO the orange band, that domain's window
and the road overlap -- the crest argmax could be locking onto the road.
```

```text
save_windows() rewrites WINDOW_JSON after EVERY domain, so a picking run
against a SHARED pick set destroys the picks of every version pointing at
it, one domain at a time, with no prompt. That is how 2009_v3 would have
been lost to the v4 re-pick. Loud, but not fatal: sharing a pick set is
legitimate when you are only topping up domains that were never picked.
```

```text
The road in the SAVED interior frame, so it is padded, cropped and
placed by exactly the same rules as the topography it sits on. Built
before the pad/crop below so it goes through both with the grid.
```

```text
pad landward, matching where dune_topo_extractor_from_GIS.py left its
sentinel: interior row 0 is the ocean side, rows increase landward
```

```text
The road can be cropped off entirely here: NC-12 sits well
landward on the wide domains, so ISLAND_PAD_ROWS = 100 loses it
where it also loses real land. That is the same loss the
`cropped` warning already reports, not a separate bug.
```

```text
No per-domain flip here: the arrays already run south -> north within a
domain, matching the ascending domain order. See the removal note at
ISLAND_INCLUDE_DUNE.
```

```text
Same origin, same columns, same clip as the grid above -- the road is
placed by the topography's rule, not its own.
```

```text
NC-12 across the whole island, in the offset frame. At 45 km wide a
single 20 m road is around one pixel, so it is drawn with
pcolormesh on the same canvas rather than as a line: that way it
cannot drift relative to the topography it was placed against.
```

```text
back to the raw cross-shore axis: k + c0 + shear[i] (shear is zeros
when STRAIGHTEN is False, so this is the original `+ c0`)
```

```text
1) measured offsets, all domains

BOTH YEARS STAY. Panel 3 is a 1984->2004 change rate and is empty without
them, and the r() notes under panel 2 are only interpretable as a pair.
What changed 2026-08-27 is WEIGHT: the year this product is actually
built for (PRODUCT_YEAR) is drawn solid and heavy, the other dashed and
faded and labelled "reference". Before, the two read as equal candidates
for the initial condition, which is exactly the confusion the plan-view
split had -- see the PRODUCT_YEAR note.
```

```text
Open on the SAVED window when there is one, so a re-pick is an
adjustment: the v4 file was seeded from v3, so each domain shows
its v3 window, "r" resets to it, and accepting without dragging
keeps it. Falling back to default_window here -- as this did before
v4 -- would have made every one of the 90 re-picks a blind redraw.
```

```text
the frame this window was picked in. A window picked
straightened is a valid index range on an unstraightened
array; it just points at different cells. The run pass
refuses on a mismatch rather than quietly using it.
```

```text
What this window replaced, so the v3 -> v4 change is auditable
from the picks file alone. Equal values mean the road showed
nothing wrong with the old window and it was accepted as-is,
which is a result worth being able to see.
```

```text
absent == picked before straightening existed == False.
NOT "unknown, proceed": that would silently apply an
unstraightened window to a straightened array, which is a
valid index range pointing at the wrong cells.
```

```text
carried for island_plan_figure, which needs the road in the saved
interior frame and so needs the masks and the row0 line together
```

<details><summary>Function notes (the original docstrings)</summary>

**`water_col_bounds()`**

```text
First/last (inclusive) columns that are not entirely water.

`<=` rather than `==`: no-data now carries NODATA_SENTINEL_M, which is below
w_elev. With `==` a column of pure no-data would read as "not water" and
survive trimming, silently changing every array's shape. Water and no-data
are both "nothing to model here" for trimming purposes, so both trim.
```

**`remove_water_rows()`**

```text
Trim leading/trailing rows that are entirely water. `<=` for the reason
in water_col_bounds: no-data must trim like water or shapes change.
```

**`orient_ocean_right()`**

```text
Return an (alongshore, cross_shore) array with the ocean in the LAST column.

NOTE: v3's "top"/"left" branches used np.flip(arr), which flips BOTH axes and
silently reverses the alongshore order. These branches do not.
```

**`straighten_profiles()`**

```text
Shear each alongshore profile so the shoreline is horizontal.

Returns (z_sheared, start_beach_sheared, shear, obliquity_deg).

shear[i] = cells dropped from the SEAWARD end of profile i. Anything that
indexes the same grid afterwards -- an NC-12 mask, a dune-line mask -- has
to be sheared with this same array via shear_like(), or it points at
different ground than the topography does.

The fit's slope is cells cross-shore per cell alongshore. Both axes are
CELL_SIZE_M, so obliquity = atan(slope): the shoreline's angle to the grid's
east-west axis, measured directly rather than inferred from dune_loc's span.
```

**`shear_like()`**

```text
Apply an existing shear to another array on the same grid.

Must be the SAME shear straighten_profiles() returned for this domain. Used
by HAT_road_setback_extract.py to put the NC-12 and dune-line masks in the
frame the topography was saved in.
```

**`align_mask_to_topography()`**

```text
Put a raw GIS mask into the frame the topography was saved in.

MOVED HERE from HAT_road_offset_from_dune_start.py (2026-08-18). It lived
there as a private copy that re-derived this script's frame from outside it,
which meant two definitions of the same chain that had to be kept in step by
hand. It is now defined once, next to shear_like, and the setback script
imports it. Any change to load_profiles' ordering has to be reflected here or
the assert at the bottom fires.

The chain must match ``load_profiles`` operation for operation, or the mask
indexes different ground than ``dune_loc`` does:

    orient_ocean_right   same OCEAN_LOC and the same ALONGSHORE_FLIP
    [:, ::-1]            ocean-first, as load_profiles does to build ``raw``
    shear_like(shear)    the SAME per-profile shear, not a re-fit one
    [:, c0:c0+n_cross]   the SAME water-trim window

``load_profiles`` returns c0 but not c1; the trimmed width of ``z`` supplies
the rest, which is also a check that the two arrays end up the same shape.
```

**`interior_row0_line()`**

```text
Source cross-shore cell that becomes SAVED interior row 0, per profile.

MOVED HERE from HAT_road_offset_from_dune_start.py alongside
align_mask_to_topography, and for the same reason: the setback measures from
interior row 0, this script's figures and road columns measure from interior
row 0, and there must be exactly one definition of where that is.

``build_interior`` with USE_CONST_INTERIOR = False fills each column from
``prof_arr[i, dune_loc[i] + 1:]``, so interior row 0 is the cell one landward
of the crest. But TRIM_INTERIOR_ROWS = True then runs ``remove_water_rows``,
which drops leading AND trailing all-water rows -- so if interior row 0 were
all-water across every profile, the SAVED row 0 would be a different cell
and every setback would be off by that shift.

It is currently zero on all 90 domains, but it is computed rather than
assumed, because a change to the dune window or the water clamp could make
it nonzero without any other visible symptom.

TWO FIXES APPLIED IN THE MOVE, both of which are no-ops on the current
settings and both of which were latent bugs in the original:
  1. the all-water test is now `<= SENTINEL_WATER_M + 1e-9`, matching
     remove_water_rows. The original used `== SENTINEL_WATER_M`, which does
     NOT catch a leading row of pure no-data (NODATA_SENTINEL_M = -99 is
     below the sentinel, so `==` kept a row that the real trim dropped).
  2. lead_trim is only applied when TRIM_INTERIOR_ROWS is True. The original
     assumed it, so with TRIM_INTERIOR_ROWS = False it would have shifted
     row 0 by a trim that never happened.
```

**`load_road_masks()`**

```text
Load every ROAD_YEARS mask for one domain, in both frames.

Returns (road_raw, road_aligned), each {year: bool array}:

    road_raw      (n_along, n_raw)   ocean-first, UNSHEARED and UNTRIMMED --
                                    the frame panel 1 of the comparison
                                    figure draws, so the road's real diagonal
                                    across the domain stays visible
    road_aligned  (n_along, n_cross) the frame the topography and dune_loc
                                    live in, via align_mask_to_topography

Empty dicts if SHOW_ROAD is False. A domain with no road cells still gets an
all-False entry rather than being omitted -- D1-D7 are that case every run,
and callers should not have to distinguish "no road here" from "not loaded".
```

**`road_profile_positions()`**

```text
Seaward edge, landward edge and centre cell of the road, per profile.

NaN on profiles the road does not cross, which is a real state -- the road
leaves the domain, or the domain has no road at all.
```

**`processed_road_grid()`**

```text
Re-index road cells into the SAVED interior grid.

The interior is cut per profile at interior row 0, so a road cell's row in
the saved array is (source cell - row0[i]). This is what shows whether NC-12
is inside the array CASCADE actually reads, and at which interior row --
which is the same quantity roadway_manager.bulldoze indexes with
int(road_setback / dy).
```

**`road_offset_stats()`**

```text
Per-domain road geometry and setback from SAVED interior row 0.

`setback_median_m` IS `setback_dunestart_m` from RoadOffset_<year>_domains.csv --
same reference row, same frame, same shear, same edge, same statistic. It is
computed here independently so the two can be diffed; if they disagree, one
of the two frames has drifted.

MEASURED FROM THE SEAWARD EDGE OF THE ROAD BLOCK, NOT ITS CENTRE, and
reported as a MEDIAN over profiles. Both of those match
HAT_road_offset_from_dune_start.py, and neither is arbitrary:

  * the seaward edge is what `roadway_manager.bulldoze` indexes --
    `[road_start : road_start + road_width]` starts at the seaward edge, so
    that is the cell `int(road_setback / dy)` has to land on. Measuring the
    centre instead reads high by half the mask width, which on a ~24 m mask
    (ROAD_BUFFER_M = 6 plus all_touched on 10 m cells) is a systematic ~10 m.
  * the median, because NC-12 leaves some domains diagonally: the handful of
    profiles where the road clips a corner drag a mean by up to ~150 m while
    the median stays on the road proper.

The centre and the width are still reported, as geometry -- they are what
tells you the mask is a fat buffer around an 8 m road rather than the road.

Sign: POSITIVE = road LANDWARD of interior row 0 (the normal case, and the
only one Barrier3D can represent). Negative means the road sits seaward of
the dune line, which is what the audit's NEGATIVE floor is about.
```

**`add_road_plan_overlay()`**

```text
Draw exact road cells on an (alongshore x, cross-shore y) map panel.

The white contour is not decoration: the fill sits on a terrain colormap that
runs from near-black water to pale land, and no single fill colour is legible
against both. The outline is.
```

**`add_road_envelope()`**

```text
Road's cross-shore envelope on a profile-stack panel (elevation x, cell y).

A profile stack has no alongshore axis, so the road cannot be drawn cell by
cell -- what it can show is the band of cross-shore cells the road occupies
anywhere in the domain, which is what you compare the search window against.
```

**`load_profiles()`**

```text
Load one domain. Returns a dict:
    raw         : (n_along, n_cross_raw) RAW GIS elevation in m NAVD88,
                  oriented OCEAN-FIRST, untrimmed. For the comparison figure.
    z           : (n_along, n_cross) MHW-relative, clamped, water-trimmed,
                  OCEAN-FIRST (index 0 = ocean). This is what gets processed.
    start_beach : (n_along,) first index in z where z > BEACH_START_THR_M,
                  -1 if none
    c0          : cross-shore raw_offset of z within raw. WITHOUT
                  straightening, z[:, k] == raw column k + c0. WITH
                  straightening the mapping is per profile:
                      z[i, k] == raw[i, k + c0 + shear[i]]
    shear       : (n_along,) cells dropped from the seaward end of each
                  profile to make the shoreline horizontal. Zeros if
                  STRAIGHTEN is False. Any mask that has to index the same
                  grid must be put through shear_like() with THIS array.
    obliquity_deg : the shoreline's angle to the grid's east-west axis,
                  from the slope of the start_beach fit. 0.0 if not
                  straightened.
    road_raw    : {year: bool (n_along, n_cross_raw)} NC-12, ocean-first,
                  unsheared and untrimmed. Empty if SHOW_ROAD is False.
    road_masks  : {year: bool (n_along, n_cross)} NC-12 in the SAME frame as
                  z, via align_mask_to_topography. Empty if SHOW_ROAD is
                  False. Display and diagnostics only -- nothing downstream
                  of here lets the road affect the dune search or the arrays.
```

**`suggest_window()`**

```text
A crest-aware starting window: (i0, i1, median crest elevation).

WHY THE SEEDED WINDOW IS NOT ENOUGH. A re-pick seeded from the previous set
opens every domain on its previous window, which is precisely the thing a
re-pick exists to question -- and if that window clipped the crest, the
picker shows no sign of it. The argmax simply pins at i1-1 and looks
plausible. This proposes a window derived from the PROFILE instead.

METHOD. Per profile, take the unconstrained argmax over
[start_beach, start_beach + SUGGEST_REACH_PX) -- i.e. the highest ground in
the foredune zone, with no window imposed. Take the median of those
locations across the domain, then bracket it with margin. The result
contains the crest by construction, which the seeded window may not.

THE ROAD IS A WARNING, NOT A LIMIT -- and the first version of this got that
wrong. Bounding the hunt at NC-12's seaward edge looks prudent (the
embankment is a flat-topped ridge an argmax locks onto) but it fails exactly
where this project lives: at GIS 85 the road sits at cell 13 and the crest
at 14, so a road-bounded hunt excluded the crest BY CONSTRUCTION and
proposed a 0.54 m "dune" against a real 4.93 m one. Wherever the island has
migrated over the roadbed, the dune IS landward of the road.

So the hunt is unbounded landward within SUGGEST_REACH_PX, and
`road_seaward` is used only to tell the caller whether the proposed window
overlaps NC-12, which the picker surfaces as a warning for the eye.

Returns (i0, i1 EXCLUSIVE, median crest elevation, overlaps_road).
```

**`window_diagnostics()`**

```text
(median crest, % of profiles whose argmax pins at i1-1, is a higher cell
sitting just outside the landward edge?).

The pin fraction is the live tell for a clipped window: a crest that really
is the last cell in the window is fine, but a crest that pins there WHILE
cell i1 is higher means the window is cutting the dune off.
```

**`pick_window()`**

```text
Show the domain ocean-at-bottom and let the user drag a dune search window.

Display only: the cross-shore axis runs vertically with cell 0 (ocean) at the
bottom, landward upward. i0/i1 are still cross-shore indices from the ocean.

NC-12 is drawn on both panels (v4). It is there to tell you when a window has
wandered onto the road embankment: the road is a hard, flat-topped ridge a
dune-crest argmax will happily lock onto, and the v3 windows were picked
without being able to see it. The road never constrains the window in code --
only your eye.
```

**`_gis_id()`**

```text
'domain_11' -> '11'.

This script names arrays by the npy-arrays stem; hat_topo_version.
array_name() names them by GIS id. They MUST produce the same filename, so
the id is derived here rather than the name being spelled twice. An
unexpected stem raises instead of silently writing a file no reader will
look for - which is exactly how the old per-script tag went wrong.
```

**`comparison_figure()`**

```text
The whole chain, left to right: what came in, what you picked on, what
CASCADE gets.

RAW           the DEM as GIS exported it (m NAVD88, untrimmed). The clip
              boxes are north-up while the island trends NNW, so the
              shoreline crosses each domain diagonally -- and so does every
              overlay, because they are mapped back through the per-profile
              shear: raw[i, k + c0 + shear[i]]. The diagonal you see here is
              what the shear removes.

STRAIGHTENED  the array the picker showed you and the window was drawn on
              (m MHW, clamped, sheared, water-trimmed). Same overlays, now
              horizontal. If the beach-start line still slopes here, the
              linear fit did not capture that domain's shoreline.

PROCESSED     the .npy files CASCADE reads, converted back to m for display.

With STRAIGHTEN = False the first two panels are the same picture, minus the
trim and the datum shift. That is the point of showing both.
```

**`road_columns()`**

```text
Per-year NC-12 columns for the settings sheet.

`road setback <year> (m)` is the cross-check column: seaward edge, median
over profiles, metres from SAVED interior row 0, positive = landward. It
should equal `setback_dunestart_m` in RoadOffset_<year>_domains.csv to the metre.
See road_offset_stats for why the edge and the statistic are what they are.

Blank rather than 0 where the road is absent: D1-D7 have no NC-12 at all, and
a 0 there would read as "road exactly at interior row 0".
```

**`resolve_offset_file()`**

```text
Return the raw_offset CSV for a year. If the configured path is wrong, search
OFFSET_DIR rather than silently skipping -- the hindcast_* subfolder naming
is a convention this script only partly knows.
```

**`_island_norm()`**

```text
Terrain with 0 m pinned to colormap position 0.35.

A LOCAL copy, and no longer shared with anything. It was written to
match the initialization figure, which moved to the house elevation
classes on 2026-09-17 (hat_figure_style.elevation_cmap: a hard break at
0 m, one colour for water). This is a QC view inside the extractor, so
it was left on the ramp rather than changed in the same pass -- but it
is now the extractor's own choice, not a shared convention.
```

**`_assert_alongshore_continuity()`**

```text
Warn if the assembled alongshore axis is discontinuous at domain seams.

A per-domain alongshore reversal is almost invisible in a 45 km plan view --
it reads as roughness -- but it puts every 500 m block backwards. The signal
is unmistakable in numbers: compare the mean jump ACROSS domain seams with the
mean jump WITHIN a domain. A continuous island sits near 1; a per-domain
reversal drove this to 21 on the 2009_v3 arrays when the legacy plotting flip
was still applied, which is the bug this guard exists to catch.

Returns the ratio, or nan when it cannot be computed.
```

**`island_plan_figure()`**

```text
Plan view of the processed dune + interior for domains 1-90 at the measured
offsets, on the terrain ramp initialization_figures.py used until
2026-09-17 (see _island_norm); that figure is now in elevation classes.

ONE FIGURE PER CROSS-SHORE MODE, at PRODUCT_YEAR's offsets only. It used to
be one per offset YEAR per mode; see the PRODUCT_YEAR note for why the
off-year figure was a topography and a shoreline from different decades.
When PRODUCT_YEAR is None (a product with no hindcast year) it falls back
to every year loaded, which is the pre-2026-08-27 behaviour.

THE ROAD IS NOT RESTRICTED. Both NC-12 vintages stay on the canvas, in
their own colours -- that overlay is the subject of
RoadOffset_dunestart_audit.md and the reason SHOW_ROAD exists here at all.
The vintage mismatch it documents is a property of the ROAD LINES, not of
which shoreline the domains are placed against, so restricting the offsets
does not make it stale.
```

**`island_figure()`**

```text
All domains together: measured dune offsets vs the crest this run extracted.

The measured offsets and the extracted crest live in DIFFERENT frames -- the
offsets are in the model's common cross-shore frame, the extracted crest is in
the per-domain raw DEM array frame (m landward of that array's cell 0). They
are NOT differenced here. They're plotted on separate panels, and the script
reports the correlation between them so the frame relationship is testable
rather than assumed.
```

</details>

