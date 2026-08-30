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
python 3-figures/HAT_road_island_planview.py          # the whole island, both years

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

| Script | Draws | Writes |
|---|---|---|
| `HAT_road_island_planview.py` | the **whole island** at once — the road as the model builds it, coloured by road elevation | `dunestart_offset/HAT_road_island_planview_<year>.png` |
| `HAT_road_domain_views.py` | **one domain** at a time, either method | `road_offset/figures/` (not kept — see below) |
| `HAT_dunestart_modification_stages.py` | what flooring and seaward relocation moved | `dunestart_offset/modifications/` |
| `HAT_oceanfloor_offset_check.py` | the ocean-side floor, domain by domain | `dunestart_offset/modifications/` |
| `HAT_road_geojson_map.py` | the raw NC-12 lines on the DEM, no extractor transform | `raster/` |

Every one of these resolves its topography per year through `YEAR_PRODUCT`.

### HAT_road_island_planview.py

Both vintages in one invocation, no arguments. **One panel**, styled as the
extractor's plan view — same ocean background, poster `terrain` ramp with 0 m
pinned to position 0.35, axes rectangle, tick cadence, spine colours and title
format — so the two overlay and read as one pair.

**It draws the model's road, not the measured one.** From `roadway_manager.py:99`:

```
road_start = int(road_setback / dy)      # dy = 10 m, TRUNCATED not rounded
road_width = int(road_width / dx)        # 20 m / 10 m = 2 cells
new_road_domain = np.zeros((road_width, ncols)) + road_ele
```

So the road is a flat rectangle — two cells cross-shore, the domain's full 50
profiles alongshore, one constant row, one constant elevation. It steps between
domains rather than curving, and that step is the discretisation the model
imposes. An earlier version drew the per-profile survey line from
`RoadOffset_<year>_profiles.csv`; that is where the road *is*, not where CASCADE
*puts* it, which is the question this figure exists to answer.

Each segment is coloured by road elevation on a second scale, `RdPu` truncated
at 0.30 — magenta is the one strong hue `terrain` does not already spend, so the
road can never be misread as ground.

Both inputs are the **model-facing** files, for the same reason:

| | file | note |
|---|---|---|
| setback | `dunestart_offset/<year>/RoadSetback_<year>_dunestart.csv` | what `PERIOD["road_setback_file"]` resolves to; already floored and relocated |
| elevation | `road_elevation/RoadElevation.csv` | **one file for both vintages** |

That second row matters when reading the pair: the road colour is identical
between the 1984 and 2004 figures *by construction*. Any difference between them
is the setback moving or the island changing — never the roadbed, because the
model is not given a per-year roadbed.

One honest exaggeration: the true 2-cell band is about **0.8 pt** at this
figure's scale, too thin to hold a fill colour beside its own outline. It is
stroked at `ROAD_LW_PT = 1.6`, and every run prints the factor:

```
[road]  stroked 1.6 pt vs 0.8 pt true width -> x2.04
```

Position and alongshore extent are exact; only cross-shore thickness is
inflated, by about two. Raising `ROAD_LW_PT` past ~2 starts claiming extent the
model does not have. The extractor hits the same wall and solves it differently
— it draws its road with `pcolormesh` at true width, which it can afford
because its road carries no second variable.

### Figure conventions

Three files per vintage: `.png` (300 dpi), `.pdf` (vector), `_caption.txt`.

* **No title on the figure.** A journal sets the caption in the text, so a
  baked-in title duplicates it, cannot be copyedited, and has to be cropped by
  hand. The caption ships as the sidecar `.txt`, generated from the same
  constants the figure is drawn with, so it cannot drift from what is plotted.
* **A provenance footer instead.** Topo product/version, the three input
  filenames, the script, cell size, vertical exaggeration and the road stroke
  factor — enough that a PNG separated from this repo is still traceable.
* **Metric axes on top and right.** Domain number and canvas cell are what you
  need to trace a value back to a file, but neither is a length. Secondary axes
  give alongshore and cross-shore distance in km, at 1 cell = 10 m.
* **Vertical exaggeration is reported** (≈×1.5), computed from the axes as
  drawn rather than assumed. A plan view whose two axes sit at different scales
  has to say so, or the island reads as narrower than it is. It lives in the
  footer, not in an axes corner — every corner is a collision waiting on a
  different island shape, and the lower-left one is where the no-road hatch
  lands.
* **Domains with no road are hatched** along the lower frame and named in the
  legend (`no NC-12 in domain 1–8`). Without it, a domain with no line is
  indistinguishable from one whose line is hidden, and the reader cannot tell
  which.
* **Fonts embed as TrueType** (`pdf.fonttype = 42`) rather than matplotlib's
  default Type 3, which a number of journals reject at submission.

### The island ramp is switchable — `HAT_ISLAND_CMAP`

```
python 3-figures/HAT_road_island_planview.py                      # terrain, default
HAT_ISLAND_CMAP=oleron python 3-figures/HAT_road_island_planview.py
```

`terrain` is the default because it is what the extractor's plan view uses, and
these two are meant to read as one pair. `oleron` (Crameri, via `cmcrameri`) is
the perceptually uniform alternative: equal elevation steps look equal, and it
survives greyscale and colour-vision deficiency, neither of which `terrain`
does. Non-default ramps write to `_<ramp>`-suffixed filenames, so a comparison
keeps both rather than overwriting.

**The sea-level position differs between the two, and getting it wrong is
silent.** `terrain` has no built-in shoreline, so the extractor pins 0 m at ramp
position **0.35** by choice, to spend more of the ramp on land. `oleron` has its
land/sea break **built in at the exact middle**, so 0 m must be pinned at
**0.50**. Pin oleron at 0.35 and its own blue-to-green boundary lands at an
elevation that is not sea level — a shoreline drawn in the wrong place, which is
worse than the non-uniform map it replaced. `SEA_LEVEL_POS` holds the pair.

Switching the default here alone would break the pairing with
`HAT_dune_topo_island_planview_*`. Change both together or neither.
`cmcrameri` is an **optional** dependency, imported only when the ramp is asked
for, so the default path does not need it installed.

### HAT_road_domain_views.py

One script, three modes:

```
python 3-figures/HAT_road_domain_views.py --domains 52 --year 1984
python 3-figures/HAT_road_domain_views.py --domains drowning --mode map
python 3-figures/HAT_road_domain_views.py --browse --start 52
```

`--year` picks the topography as well as the setback: 1984 draws on
`1984-start`, 2004 on `2004-start`. It used to resolve one product at import,
before the CLI was parsed, so one of the two years was always drawn on the
other one's island.

`--mode map` is the land/water plan view, `--mode section` runs the **real**
`bulldoze()` on a copy of the real interior, `--browse` walks domains
interactively. `--method legacy|dunestart` picks which setback is drawn.

Its output folder, `road_offset/figures/`, is **not kept on disk** — 330 PNGs
that nothing reads and `.gitignore` never tracked. Regenerate what you need.

## 4-compare/

Cross-method only. Writes to `road_offset/method_comparison/`, not into either
method's folder, because neither output belongs to one method — and not to the
top level either, which is for the forcing, its source and its inputs. Both
scripts create the folder on demand, so it can be deleted whole.

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

**2026-08-28 (same day) — the data folder made self-describing.** Seven things
in `data/.../road_offset/` were not intentional. All are fixed; the folder now
carries its own `README.md` mapping every file to the script that writes it.

* **A name collision between a live audit and a dead one.**
  `RoadSetback_audit.{csv,md}` existed in *both* `dunestart_offset/` (current)
  and `old_method_offset/` (2026-08-17, orphaned). `HAT_road_setback_audit.py`
  hardcodes `dunestart_offset` at line 152, so nothing had written the second
  pair in eleven days — but the filename gave no way to tell them apart. The
  orphans are deleted; the legacy method keeps its correctly-named
  `RoadSetback_oldmethod_audit.md`.
* **`HAT_road_method_diagnostic.py` wrote into `dunestart_offset/`**, putting a
  cross-method result inside one of the two methods it compares and
  contradicting the rule this README states for `4-compare/`. Its sibling was
  always correct; only this one drifted.
* **All four cross-method outputs then moved to
  `road_offset/method_comparison/`.** They had been loose at the top level,
  sharing a directory listing with the forcing tree and its inputs with nothing
  to mark them as a different kind of thing — you read them when choosing a
  method, not when running the model. Nothing in the repo reads any of the four,
  so the move cost two `OUT_ROOT` lines and carried no risk. Verified by
  deleting the folder outright and letting both scripts rebuild it.

  `old_method_offset/` deliberately did **not** move under it. Six live code
  paths read that folder, one of them the producer of the model forcing via
  `read_two_row_csv()`, which returns `{}` rather than raising — and this exact
  directory was renamed once before, on 2026-08-17, taking nine files with it.
  It also remains a true sibling of `dunestart_offset/`: they are two methods,
  and nesting one under "comparison" would encode an importance ranking in the
  tree that both READMEs already state in words.
* Two stale references in the producer's own header: it cited the orphan audit,
  and named a column `setback_sameyear_m` that was renamed `setback_legacy_m`
  on 2026-08-26.
* `raster/<year>/RUN_MANIFEST.txt` recorded `CLIPRESAMPLE_ROOT` and
  `ELEV_NPY_DIR` under `1-barrier3d-domains/2009-raw/`, renamed by the 08-25
  restructure. The SETTINGS blocks are left unedited — they are what happened —
  and a dated footer gives the current names. The masks are unaffected.
* `raw_offset/<year>/nc12_<year>.csv.xml`, 722 KB of ArcGIS sidecar, removed
  after the one fact inside it (source `\\HANNAHS-LAPTOP\D$\Hatteras_GIS\Roads\`,
  created 2024-10-10) was written into the data README. Same treatment, same
  reasoning, as the `.tif.xml` purge of 2026-08-26.

Checked and found correct, so left alone: the 14-domain `raster/<year>/figures/`
subset is `QC_DOMAINS` by design, the `_rawframe` CSVs are the real
unstraightened control product, and no directory is empty.

**2026-08-28 (same day) — `HAT_road_island_planview.py` added.** There was no
view of the road against the *whole* island on the real elevation surface:
`HAT_dunestart_road_on_domains.py` draws a per-domain grid and
`HAT_road_domain_views.py` draws one domain at a time. This builds the
extractor's own canvas and its poster styling, so it overlays
`HAT_dune_topo_island_planview_<ver>_<year>_padded.png` and reads as its pair.
One figure per vintage, each on its own product.

It went through two wrong drafts worth recording. The first drew the
**per-profile survey line** and three stacked panels; that is where the road
is, not where CASCADE puts it, and the extra panels were not asked for. The road
is a flat 2-cell band per domain — `roadway_manager.py:99` truncates
`int(setback/10)` and fills 2 cells with one elevation — and the figure now draws
that, coloured by `RoadElevation.csv` on a second scale.

Writing it surfaced two defects in this README's `3-figures/` section, both now
fixed: it said "One script, three modes" while the folder held five scripts, and
the paragraph explaining `--year` was inside the code fence, so it rendered as
code rather than prose.

What the pair shows, worth knowing before it is read as a survey artefact: in
**2004 the road sits below the 1.34 m MHW berm nearly everywhere** (median
0.5-0.9 m MHW), while in **1984 domains 8-15 sit at or above it** (1.8-2.0 m).
Road elevation is one file for both vintages — `../road_elevation/RoadElevation.csv`,
sampled on the 2009-2014 baseline — so that contrast is the 1984 line sitting on
higher ground, not a change in roadbed height between the two dates.

**2026-08-28 — 1984 re-measured on the live `1984-start/v1`.** The 1984 numbers
in this tree had been measured against **`1984-start/v2`** using
**`HAT_dune_search_windows_v2.json`**. Both were deleted on 2026-08-27 when
`1984-start` was cleared back to a blank slate and re-extracted as a new `v1`
from the *same* DEM and the same `npy-arrays/` against a *new* pick set. Same
ground, different interior row 0 — which is the only thing this method measures
from.

Nothing errored, and nothing would have: `RoadOffset_dunestart_audit.md` named
`v2` in its own provenance table while `topo_dirs("1984-start")` resolved `v1`,
because after the clear there was exactly one version on disk for it to find.
This is the v3/v4 incident's failure mode with the versions renumbered.

* **15 of 83 road-bearing domains moved**, max 25 m, median 0. The road span is
  unchanged: no domain gains or loses a road, and 0 drown at initialisation
  either before or after.
* The movers are badly placed rather than numerous. **GIS 80: 565 → 540 m** —
  one of the three roadways the relocation logic acts on. **GIS 10-13:
  30→20, 30→10, 40→20, 40→30 m** — at 10 m cells that is `int(setback/10)`
  halving, so the road lands on a different row.
* **Cross-check that did not exist before.** The extractor now writes its own
  `road setback 1984 (m)` column, computed from the road *raster* against the
  interior it just saved. All **83 road domains agree to 0.00 m** with the
  re-measured `setback_dunestart_m`. The producer and the extractor are
  provably reading one interior — which is exactly what was untrue on 08-27.
* Against the rasterized road, the dune-start method is now median **+0 m**,
  mean |err| **0 m**, max |err| **10 m** for 1984 (one cell, worst case).
* **2004 is byte-identical again**, all five files, verified by diff against
  copies taken before the run. `2004-start/v1` never moved.
* `old_method_offset/` was **not regenerated** — it reads ArcGIS CSVs
  (`raw_offset/<year>/nc12_<year>.csv`, `2-brie-offset/raw_offsets/`), carries
  no topography dependency, and is not stale. Only its
  `HAT_old_method_road_on_domains.png` changed, because that figure draws the
  legacy setback on the new interior. Every legacy CSV and audit is untouched.
* `road_elevation/` and `road_relocation/` were not re-run: neither imports
  `hat_topo_version`, and both are frame-independent (an elevation and a
  line-to-line displacement).

Contrary to a note that circulated for a while, the **legacy method is
reproducible**. What is missing from the repo is the dune-line *geojsons*; the
pipeline never read them — it reads the exported offset CSVs, which are present
for both vintages.

**`figures/` is no longer kept on disk.** `HAT_road_domain_views.py` wrote 330
PNGs (70 MB) there; nothing read them, `.gitignore` line 158 (`*.png`) meant
none were tracked, and every 1984 panel in the set had been drawn on a deleted
topography until this rebuild. The folder was removed on 2026-08-28. The script
recreates it on demand — `out_dir.mkdir(parents=True, exist_ok=True)` runs
before every write — so regenerate whichever view you actually need:

```
python 3-figures/HAT_road_domain_views.py --domains 52 --year 1984   # one domain
python 3-figures/HAT_road_domain_views.py --domains all --year 1984  # the full set
```

`old_method_offset/` was **kept**. It is 1.9 MB and six tracked files, and
three live consumers read it: the producer at line 920 for `setback_legacy_m` /
`delta_vs_legacy_m`, `HAT_road_placement_on_domains.py` for its `old` entry,
and both `4-compare/` scripts. Deleting it would not raise —
`read_two_row_csv()` returns `{}` for a missing path — it would publish an
audit markdown whose migration diagnostic is silently all-NaN while the prose
explaining that diagnostic stays put. Removing the method is a code change, not
a folder delete.

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
