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
| Output | `old_method_offset/<year>/` | `dunestart_offset/measured/<year>/` |
| Drowns at t=0 | **0 (1984)**, 5 (2004) | **0, both years** |
| vs the rasterized road | +29 m / +22 m median | 0 m *by construction* |

Both are 2 rows × 82 cols, GIS IDs then metres, so switching is a drop-in.

**What the runner spends today** (`scripts/site_layer/hatteras_site_config.py:96,110`) is
the **dune-start** method, switched 2026-08-18. `3-figures/HAT_road_domain_views.py`
defaults to `--method dunestart` to match; change both together, or the figures
stop describing the model you run.

**Each vintage is built on its OWN topography** — 1984 on `1984-start/v1`,
2004 on `2004-start/v1`. They are different islands: all 90 domains differ and
**65 differ in interior shape**. The pairing is defined once, as `YEAR_PRODUCT`
in `scripts/site_layer/hat_topo_version.py`, and every script here resolves through it.
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
# HAT_ROAD_YEAR is the LINE vintage, 1978 or 2008 -- not the period. Which
# period reads which line is hat_topo_version.ROAD_LINE_FOR_YEAR (2026-09-15).
HAT_ROAD_YEAR=1978 python 1-produce/HAT_rasterize_road_to_domains.py
HAT_ROAD_YEAR=2008 python 1-produce/HAT_rasterize_road_to_domains.py

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
| `HAT_rasterize_road_to_domains.py` | `raster/<vintage>/masks/` — **the only script that masks the road**. Set `HAT_ROAD_YEAR` to the LINE vintage (1978 or 2008); do not save a per-vintage copy |
| `HAT_road_offset_from_dune_start.py` | `dunestart_offset/measured/<year>/` — setback, elevation, per-domain and per-profile detail, audit |
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
| setback | `dunestart_offset/measured/<year>/RoadSetback_<year>_dunestart.csv` | what `PERIOD["road_setback_file"]` resolves to; already floored and relocated |
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

**2026-09-15 — two axes, two integers; measured/ and derived/.** Hannah's
call, after an interview on how the forcings should be organised by start
year. Three things changed in `data/.../4-mgmt-forcing/`, and no number moved:

* **The lines and masks are filed by their TRUE vintage.** `raw_offset/1984/`
  → `raw_offset/1978/` (`nc12_1978.*`), `raw_offset/2004/` →
  `raw_offset/2008/`, `raster/1984/` → `raster/1978/` (masks
  `domain_<N>_road_1978.npy`), `raster/2004/` → `raster/2008/`, and the
  relocation measurement `road_relocation/1984_2004/` →
  `road_relocation/1978_2008/`. Before this the same integer meant a PERIOD
  under `dunestart_offset/` and a LINE under `raw_offset/` and `raster/`. The
  pairing lives once, in `hat_topo_version.ROAD_LINE_FOR_YEAR` (1984, 1996 →
  1978; 2004, 2010 → 2008), with `road_line_file()`, `road_mask_dir()` and
  `road_mask_file()` beside it; each refuses a period start where a line
  vintage is expected. `HAT_ROAD_YEAR`, `HAT_RELOC_FROM/TO` and the
  extractor's `ROAD_YEARS` are line vintages now; the rasterizer refuses 1984
  and 2004. The masks were renamed, not re-burnt (see the 2026-09-15 footer in
  each `RUN_MANIFEST.txt`, whose body is left as written).
* **`dunestart_offset/` is split into `measured/{1984,2004}` and
  `derived/{1996,2010}`.** Two of the four setback files were never
  measurements -- 1996 is 1984 plus the 1989 event, 2010 is a copy of 2004 --
  and nothing but a PROVENANCE.md said so. `hat_topo_version.ROAD_SETBACK_KIND`
  owns the split and `road_setback_dir()` / `road_setback_file()` /
  `road_setback_relpath()` resolve it; `HATTERAS_PERIODS[*]["road_setback_file"]`
  is built from the last of those rather than spelled out. Every script that
  formatted `dunestart_offset/{year}` now goes through the helper or spells
  `measured/` explicitly (the figure and compare scripts only ever loop over
  the two measured years).
* **The organising rule is written down**, in `data/.../4-mgmt-forcing/README.md`:
  state at year zero (the setback) is per start year; history (relocation
  events, nourishment projects) is one timeline the run window selects from;
  elevation has no year. 1996 stays derived -- no 1997 NC-12 line is planned.

Not done here, and done on 2026-09-18: four live scripts still spelled
`old_method_offset/`, a folder that became `superseded_20260911/` on 2026-09-11
(`HAT_road_offset_from_dune_start.py`, `HAT_road_placement_on_domains.py`,
`HAT_road_domain_views.py`, `HAT_road_method_diagnostic.py`). They failed soft,
as the data README warns, so from 09-11 to 09-18 the producer's
`setback_legacy_m` / `delta_vs_legacy_m` diagnostic was empty in any audit it
wrote; the setbacks themselves never read it. All four now resolve the folder,
at `road_offset/archive/superseded_20260911/`, through
`hat_topo_version.LEGACY_SETBACK_ROOT`, as every script now resolves every
4-mgmt-forcing path.

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

* `YEAR_PRODUCT` is now defined **once**, in `scripts/site_layer/hat_topo_version.py`, and
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

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### 1-produce/HAT_rasterize_road_to_domains.py

Burn a road geojson onto each domain's grid, so the setback is measured along the same profiles as the interior.

From the script's original header:

```text
Burn a road geojson onto the exact per-domain grids produced by the CASCADE
Clip & Resample Domains tool, so the setback can be measured along the same
profiles that define the interior.

THE ONLY SCRIPT THAT MASKS THE ROAD. Set ROAD_YEAR in CONFIG and run it:

    python HAT_rasterize_road_to_domains.py

The masks it writes are consumed by HAT_road_offset_from_dune_start.py (the
setback) and by the audit/ and figures/ scripts.

The road-side twin of gis-export-npy.py: same folder walk, same grids, same
array orientation -- rasterizing NC-12 instead of converting the DEM.

THE ONLY THING THAT MATTERS
Grid alignment. The road mask must be cell-for-cell identical to the elevation
array: same shape, same affine transform, same rotation (if any).

This script never guesses. It opens each domain's resampled_*.tif and takes the
affine straight from the raster, so alignment holds whether the Clip & Resample
tool produced north-up or shoreline-rotated grids.

OUTPUT LAYOUT
Everything for one road vintage lands in ONE folder, so comparing vintages
means comparing two directories:

  data\hatteras_init\4-mgmt-forcing\roads\raster\1978\
      RUN_MANIFEST.txt                        every setting that made this folder
      masks\
          domain_9_road_1978.npy              <- the setback scripts read these
          domain_10_road_1978.npy
          ...
      HAT_road_mask_diagnostics_1978.csv      per-domain numbers, incl. elevation
      HAT_road_mask_summary_1978.png          all domains on one page
      figures\
          domain_009_road_mask.png            per-domain map + profile

NOTE ON THE ELEVATION COLUMNS
road_elev_* is the elevation of the cells the mask landed on, MHW-relative
(the extractor subtracts MHW_M = 0.36 before anything else, so this matches
the frame CASCADE runs in).

Do NOT read a low value as "misregistered" without looking at the figures. Two
reasons the crown may not survive:

  1. ROAD_BUFFER_M = 6 plus all_touched=True gives a ~24 m wide mask on 10 m
     cells. NC-12 is ~8 m. Most of every "road" cell is shoulder and adjacent
     ground, and the median follows them.
  2. The DEM was resampled to 10 m. An 8 m road inside a 10 m cell is averaged
     with whatever else is in that cell.

So road_elev may be measuring the DEM's resolution rather than the road. The
LOW_ELEV flag is INFORMATIONAL. What actually proves registration is the
figures: does the mask trace the island, and does it sit where you digitized it.

The setback does not depend on any of this -- it comes from the geojson
geometry. The elevation only matters if you want a per-domain road_ele.

REQUIREMENTS
  rasterio, geopandas, shapely, numpy, matplotlib
```

Notes that were in the code:

```text
Derived from this file's own location, never hardcoded: this constant had been
reduced to Path("/") -- every path below resolved to the filesystem root and
the script died on "missing clipresample root". scripts/input_prep/
4-mgmt-forcings/road_offset/1-produce/<this file> -> parents[5].
```

```text
Same root gis-export-npy.py walks: one subfolder per domain, each holding
a resampled_*.tif. The REPO copy, not the OneDrive original, so a run does
not need the network drive.
MOVED 2026-08-25: 1-barrier3d-domains went period-first and the pre-90-domain
legacy went under superseded/. Path repointed so this script keeps reading
EXACTLY what it read before - no road number moves because of the reorg.
MOVED 2026-08-26 out of superseded/. These are LIVE INPUTS, not an
archive - four scripts read them, including the one that builds the
road masks the dune/topo extractor requires. Keeping them under a
directory called 'superseded' invited exactly the deletion this move
prevents.
```

```text
A GEOMETRY REFERENCE, NOT "THE ARRAYS THE EXTRACTOR READS" (corrected
2026-08-26). That is what this comment used to claim, and it stopped being
true when 1-barrier3d-domains went period-first: the extractor now reads
<product>/npy-arrays, one per period, and nothing reads this directory as an
input to a forcing.

What the cross-check needs is the GRID, not the elevations - the mask must
land on the same rows and columns the extractor will later shear and trim.
Every one of these arrays is (50, 200), identical in npy-arrays_2009_unfilled
and in both products, because all three are clipped from the same
resampled_*.tif boxes and the fills change values, never extents. So this
stays a valid alignment reference and is deliberately NOT repointed at a
product: doing that would make the road masks - which BOTH periods share -
depend on one period's DEM.

It held 131 arrays against the live 90; the extras were pre-90-domain legacy
and were purged on 2026-08-26, leaving domains 1-90. find_elev_npy() looks
up by domain id, so the count never mattered here either way.
```

```text
--- ROAD ---------------------------------------------------------------
THE ONE LINE TO CHANGE PER VINTAGE. Edit it and re-run; do NOT save a
per-year copy of this file. Everything below is derived from it, and the
output lands in its own road_offset\raster\<year>\ folder with a RUN_MANIFEST, so
the vintages stay separable without a second script.

This file used to be driven by HAT_rasterize_road_run.py, which patched these
globals at import. That driver is gone; its settings are folded in here, so
there is ONE rasterization implementation and one place for the year and the
paths. The rules below are the ones road_offset\raster\1978\RUN_MANIFEST.txt
records, so a 2008 run is made by the same rules as the 1978 masks on disk.

Overridable from the shell so building BOTH vintages needs no edit and no
second copy of this file -- which is the rule above, not an exception to it:
HAT_ROAD_YEAR=1978 python HAT_rasterize_road_to_domains.py
With the variable unset the constant below is what runs, so the file still
reads as "the one line to change".

ROAD_YEAR IS A LINE VINTAGE, NOT A PERIOD START (2026-09-15). The lines are
1978 and 2008 exports and are now filed under those years; which period reads
which is hat_topo_version.ROAD_LINE_FOR_YEAR. Passing 1984 or 2004 here is
refused rather than resolved to a folder that no longer exists.
```

```text
Per-domain map figures. None = every road domain (slow, ~110 figures);
a list = just those. Defaults to the domains this investigation cares about:
the 1999 relocation block, the crash suspect, the LOW_ELEV heartland, the
NODATA cases, and the 1989 Pea Island block.
```

```text
--- RASTERIZATION ------------------------------------------------------
Widen the road before burning. The geojson is a zero-width centerline, and a
zero-width line through a 10 m grid can skip cells diagonally, leaving gaps a
profile falls straight through. NC-12 is ~2 lanes plus shoulders.
Set 0.0 to burn the bare centerline: with ALL_TOUCHED the path stays
continuous, and road_elev then samples only the cells the road crosses. That
is the cheap test for whether the crown is resolved at 10 m.
```

```text
Display only. Must match OCEAN_LOC in HAT_dune_topo_extractor.py so the
figures show the same orientation the extractor works in.
```

```text
Check the inputs BEFORE mkdir. Carried over from the retired
HAT_rasterize_road_run.py driver, and it is not decoration: the mkdir
below is parents=True, so a wrong PROJECT_ROOT silently builds a whole
stray tree and the run looks normal. That is exactly how the leftover
C:\data\hatteras_init\ came to exist. Fail loudly instead.
```

<details><summary>Function notes (the original docstrings)</summary>

**`crs_to_metres()`**

```text
Metres per linear unit of src's CRS. EPSG:2264 (NC State Plane) is US
SURVEY FEET, so buffering by 6.0 there is 6 ft = 1.8 m, not 6 m.
```

**`ocean_first()`**

```text
(n_along, n_cross) with index 0 = ocean, matching what the extractor works
in after orient_ocean_right() -> [:, ::-1]. Display only.
```

**`domain_figure()`**

```text
Confirm, by eye, that the mask landed on the island where it should.

Left  : elevation map, ocean at the bottom, road mask overlaid.
Right : alongshore-mean cross-shore profile, with the road's cross-shore
        span shaded and the berm marked.
```

</details>

### 1-produce/HAT_road_offset_from_dune_start.py

NC-12 setback and road elevation per domain, measured from interior row 0, the reference CASCADE indexes against.

Notes that were in the code:

```text
HAT_road_offset_from_dune_start.py

NC-12 road setback and road elevation per Barrier3D domain, measured from the
SAME reference CASCADE indexes against: interior row 0 of the extracted
topography, which is one cell landward of the picked dune crest.

WHY THIS EXISTS
`roadway_manager.bulldoze` places the road with

road_start = int(road_setback / dy)                 # roadway_manager.py:99
old_road_domain = xyz_interior_grid[road_start:road_end, :]

so `road_setback` is metres landward of interior row 0, and
`cascade_pipeline/roadway.py:147` already documents it as "metres landward of
the dune line". This script is the first measurement that actually honours
that convention: it re-derives the dune line with the same picked windows and
the same straightening the topography was built with, then measures the road
against it, per alongshore profile.

It also removes a large error source. The clip boxes are north-up while the
island trends NNE, so NC-12 crosses each 500 m domain diagonally: on GIS 11
the road spans raw cross-shore cells 151-161, i.e. 110 m of cross-shore
extent for a 20 m road. Any method that reduces that to a raw cross-shore
median inherits the smear. Shearing the road mask with the SAME per-profile
shear as the topography (`shear_like`) collapses it.

TWO REFERENCE FRAMES, BOTH REPORTED
The topography is 2009; the road_offset are 1984 and 2004. Those disagree, and the
disagreement is the subject of old_method_offset/RoadSetback_oldmethod_audit.md.
This script does not pick a winner:
* setback_dunestart_m   measured directly here, road-vs-2009-dune. Internally
consistent with the grid CASCADE actually runs.
* setback_legacy_m  read from the EXISTING RoadSetback_<year>.csv, which
is a road-vs-same-year-dune measurement. No
extrapolation is performed to produce it.
* delta_m          the difference, which should track the dune-line
retreat between <year> and 2009. The 2-brie-offset
files are used ONLY to validate that, never as an input.

OUTPUTS (nothing existing is overwritten; new tree under dunestart_offset\)
dunestart_offset\measured\<year>\RoadSetback_<year>_dunestart.csv    2-row, model-facing
dunestart_offset\measured\<year>\RoadElevation_<year>_dunestart.csv  2-row, m MHW
dunestart_offset\measured\<year>\RoadOffset_<year>_domains.csv       per-domain detail
dunestart_offset\measured\<year>\RoadOffset_<year>_profiles.csv      per-profile detail
dunestart_offset\RoadOffset_dunestart_audit.md                       the write-up

measured\ because these ARE measurements: a digitised line against the
period's own extraction. The derived\ sibling (1996, 2010) is written by
HAT_road_setback_derived_vintages.py FROM these, never by this script.
Which line each year reads is hat_topo_version.ROAD_LINE_FOR_YEAR: 1984
reads the 1978 line's masks, 2004 the 2008 line's (2026-09-15).

WHERE THIS DEPARTS FROM "MEASURE, DON'T CORRECT" -- TWO PLACES, BOTH FLAGGED
Both act ONLY on the model-facing CSV. `setback_dunestart_m` in the _domains.csv
is always the measurement, and every adjusted domain carries a flag naming
what was done to it. (1) is the negative floor, below. (2) is the seaward
relocation of roadways that drown at initialisation, documented at
RELOCATE_DROWNING.

(1) NEGATIVE SETBACKS ARE FLOORED
A 1984 road measured against the 2009 dune line can land SEAWARD of interior
row 0, giving a negative setback. A negative cannot be written to the
model-facing file, because `int(-50/10) = -5` and
`xyz_interior_grid[-5:-3, :]` is valid Python that indexes from the LANDWARD
end -- the road would be bulldozed into the bay, silently. Nor is 0.0 a
"no road" sentinel: `build_roadway_management_on` decides which domains are
managed from the road SPAN, not the value, so a 0.0 inside the span still
means "managed, road at row 0".

So the model-facing CSV carries max(setback, 0.0) and the true signed value
is preserved in the _domains.csv and the audit, flagged NEGATIVE. That is a
floor, and it is called one everywhere it appears.

USAGE
python HAT_road_offset_from_dune_start.py
```

```text
Anchored 2026-09-14: absolute into a home directory, or into a tree
renamed since. Rule 5 of ORGANIZATION.md.
```

```text
There is now exactly ONE copy of HAT_dune_topo_extractor.py in the repo, and
this is it. Until 2026-08-26 there were four, only one of which had
ALONGSHORE_FLIP = True; importing a wrong copy measured the road in the
mirrored alongshore frame -- the exact mismatch this whole exercise is meant
to eliminate. The duplicates are deleted, so that particular trap is gone,
but keep resolving the extractor by this single path rather than by a
relative one.
```

```text
MASKS ARE KEYED BY LINE VINTAGE, NOT START YEAR (2026-09-15). raster/1978/
and raster/2008/ hold the 1978 and 2008 lines; the 1984 start reads the
first, the 2004 start the second, through hat_topo_version.ROAD_LINE_FOR_YEAR
and road_mask_dir() / road_mask_file(). Nothing here spells "raster/<year>".
```

```text
The existing same-year measurements, used as the second reference frame.
old_method_offset/ became a dated superseded folder on 2026-09-11 and this
kept naming it, so the second frame was silently absent; resolved 2026-09-18.
```

```text
Offset files, used ONLY to validate delta_m against measured retreat.
The flat <year>/ build this named was versioned on 2026-09-15, so the
check has been skipped since; each start's CURRENT build now (2026-09-18).
```

```text
EACH YEAR BELONGS TO A TOPOGRAPHY PRODUCT (2026-08-26).

The setback is measured from interior row 0, and row 0 is a property of the
EXTRACTION - which is now period-specific: 1984-start is built on the
1996-grafted DEM, 2004-start on the plain 2009+2014 one. Measuring the 2004
road against 1984-start row 0 and writing it to dunestart_offset/measured/2004/ with
nothing saying so is the class of error hat_topo_version.py exists to
prevent, and it is silent - the numbers look plausible.

FIRST FIX, SAME DAY, AND WHY IT WAS NOT ENOUGH. The first version of this
guard SKIPPED any year whose product did not match the extractor's own
TOPO_PRODUCT literal, telling you to repoint the extractor and re-run. That
made every run half a run - and write_audit() rewrites ONE markdown for the
whole tree, so the 14:04 1984-only run published an audit with no 2004
section at all, for a forcing that had not changed.

Now load_extractor(product) configures a SEPARATE extractor module per
product, so one invocation measures both vintages, each against its own row
0, and the audit is whole. The year -> product mapping is imported from
hat_topo_version.YEAR_PRODUCT rather than spelled here; a local literal is
how the figure scripts came to disagree with this one.
```

```text
--- ROAD SPAN ----------------------------------------------------------
Domains written to the model-facing CSVs. Matches the legacy files and the
rasterizer's own FIRST_ROAD_DOMAIN/LAST_ROAD_DOMAIN = 9, 90.

D8 is measured and reported but EXCLUDED, on evidence rather than convention:
at the Buxton bend NC-12 turns to run nearly east-west, PARALLEL to the raster
rows, so it crosses D8's corner rather than crossing the domain shore-normal.
HAT_check_geojson_vs_mask.py measures the consequence -- the road spans a
median of 75 cells per row there (102 max) against 2 cells in a normal domain,
and it touches only 25 of 50 profiles. A cross-shore "distance landward of the
dune" is not a meaningful quantity for a road running alongshore, so a scalar
setback for D8 would be a number without a physical reading.

Domains outside this span are still measured, still land in the _domains.csv,
and are flagged EXCLUDED_FROM_SPAN so the exclusion is visible rather than
silent.
```

```text
--- UNSTRAIGHTENED CONTROL PASS ----------------------------------------
A second measurement with STRAIGHTEN = False, so the old-vs-new difference can
be attributed instead of just reported. Holding the dune feature, the road
mask, the aggregation and this code fixed and changing ONLY the frame splits
the total change exactly:

new_straightened - legacy  =  (new_straightened - new_raw)   <- frame
+ (new_raw          - legacy)    <- dune feature

The control needs the pre-straightening pick set, because a window picked in
one frame is a valid index range in the other that points at different cells --
which is what the extractor's `straightened` guard exists to catch. That file
predates the flag (straighten: null), so the guard reads it as False and lets
the control through, which is correct.

Residual caveat, stated rather than hidden: the two passes necessarily use two
different pick files, so a small part of the "frame" component is really
window-choice. The windows were drawn to bracket the same dune in each frame,
so it is second-order, but it is not zero.

ALONGSHORE_FLIP is deliberately NOT varied. The setback is one scalar per
domain, collapsed by median over profiles; reversing the profile order leaves
that multiset unchanged, so the flip cannot move the setback at all. It decides
where the road sits WITHIN a domain, which is a different question.
```

```text
MOVED 2026-08-25: 1-barrier3d-domains went period-first and the pre-90-domain
legacy went under superseded/. Path repointed so this script keeps reading
EXACTLY what it read before - no road number moves because of the reorg.
```

```text
Assumed roadway width. RoadwayConfig.road_width_m is a single global 20.0
(cascade_pipeline/roadway.py:45), so a per-domain width is not consumable
today -- but cascade.py:43 does accept an array, so the MEASURED width is
reported as a diagnostic in case that changes.
```

```text
--- SEAWARD RELOCATION OF DROWNING ROADWAYS ----------------------------
The second departure from "measure, don't correct". Decided 2026-08-18.

WHAT IT DOES
A roadway whose flanking rows are >20% water at t=0 is width-drowned by
`roadway_manager.bulldoze` on the first call: RoadwayManager sets
_drown_break, cascade.py never calls update() again, and that domain spends
the entire hindcast as an UNMANAGED barrier wearing a road label -- no
overwash removal, no dune rebuilding, no relocation. Rather than lose the
domain, the road is moved to the nearest row SEAWARD that survives the test.

WHAT "NEAREST VIABLE" MEANS
The largest road_start < current for which bulldoze's own flanking test
passes -- `drown_test` below is that test, not an approximation of it. The
road's own band is NOT required to be dry, because bulldoze never looks at
it; requiring it would move roads further seaward than the model needs.
There is NO cap on the distance: the nearest viable row is taken however far
it is. See the caveat below for what that costs on GIS 79.

THE ASSUMPTION THIS RESTS ON, STATED PLAINLY
On the 2009_v3 topography the domains this fires on do NOT drown on measured
water. They drown on LiDAR coverage gaps. Across the six flanking rows that
fail at GIS 78/79/80 there are 106 wet cells, of which 105 were NEVER
SURVEYED and 1 is genuinely measured wet. The extractor writes no-data back
as SENTINEL_WATER_M because Barrier3D has no representation for "unknown",
and `wet_fraction` counts every cell (see the audit script -- that is
deliberate, so the figure reports what the run does).

So this relocation is NOT "the road was in water, we moved it out". It is
"the 2009 DEM has no data there, CASCADE reads no-data as water, so the road
is moved onto surveyed ground to keep the domain managed". Anyone reading a
managed-vs-unmanaged result at GIS 78-80 needs that sentence.

WHAT IT COSTS, PER DOMAIN
GIS 78  490 -> 470 m ( 20 m)  lands on the measured per-profile MINIMUM
GIS 80  490 -> 450 m ( 40 m)  inside the measured spread, ~p10 (already
flagged SCATTER_SETBACK(71m))
GIS 79  510 -> 400 m (110 m)  90 m SEAWARD OF THE SEAWARD-MOST PROFILE ever
measured in that domain. This one is not a
re-pick of the alongshore statistic; it is a
position no profile showed. Accepted knowingly
(no cap), and flagged BEYOND_MEASURED so it can
never be mistaken for a measurement.

Set RELOCATE_DROWNING = False to write the measured setbacks unchanged and let
the domains drown; everything is still measured and flagged either way.
```

```text
bulldoze's own constants. Not tuning knobs -- roadway_manager.py:50-51,125-165;
RoadwayManager re-asserts 0.2 at :531; drown_threshold = 0 at :766.
```

```text
The whole assignment line is replaced, trailing comment and all --
matching the quoted value would have to spell both quote styles, and
the comment on VERSION ("bump this per settings variant") is not
something this loader should try to preserve into a patched copy.
```

```text
FRAME ALIGNMENT

align_mask_to_topography() and interior_row0_line() USED TO LIVE HERE. They
were private copies that re-derived the extractor's frame from outside it --
two definitions of one chain, kept in step by hand. They now live in
HAT_dune_topo_extractor.py next to shear_like(), which is the only place that
knows the frame, and this script calls them through `ext`:

ext.align_mask_to_topography(raw_mask, dom)     (no `ext` first arg)
ext.interior_row0_line(prof_arr, dune_loc)

The extractor needed them anyway, to draw NC-12 on the picker and to report
road distances from interior row 0 in its settings sheet. Moving rather than
duplicating is what keeps this script's setback and the extractor's diagnostic
columns the same measurement.

TWO LATENT BUGS WERE FIXED IN THE MOVE, both no-ops on the current settings:
interior_row0_line's all-water test now matches remove_water_rows (`<=`, so a
leading row of pure no-data trims like the real thing instead of being kept),
and lead_trim is only applied when TRIM_INTERIOR_ROWS is True. lead_trim is 0
on all 90 domains either way, which is why the outputs did not move -- verified
by diffing the CSVs before and after the refactor.
```

```text
GEOLOCATION  (diagnostics only -- never feeds a setback)

Put the measured per-profile cells back on the map, so RoadOffset_*_profiles.csv
can be loaded into GIS and eyeballed against the road instead of being read as
bare cell indices. Borrowed in spirit from
roya_files/road_setback_from_road_json_to_interior_start_REVISED.py, which works
natively in map coordinates -- but NOT its method: that script reconstructs the
interior point by walking (c0 + interior_start + 0.5) * cell metres from the
domain polygon edge along a PCA axis, a formula with no shear term, which would
put back the obliquity this pipeline exists to remove. Here the affine comes
straight from the same resampled_domain_*.tif the mask was burned onto, and the
index chain is inverted exactly.
```

```text
MOVED TWICE, AND THE SECOND MOVE WAS MISSED UNTIL 2026-08-26.
2026-08-25: 1-barrier3d-domains went period-first and domain-clips-1m went
under superseded/ with the pre-90-domain legacy, so this path was repointed
there. 2026-08-26: it was moved back OUT - it is a LIVE input, read by four
scripts including the rasterizer that builds the road masks - and this path
was not updated with it.

Nothing errored, because the coordinate write-back is optional: every profile
simply got a blank road_x/road_y and one "[coords] no resampled tif" line.
Those two columns are what the audit's independent check of the index
inversion rests on - shapely distance from the reconstructed road points to
the digitised geojson, 6.62-6.68 m median - so losing them quietly loses the
check that would catch an alongshore-flip error.
```

```text
Only OCEAN_LOC = "right" is invertible with the simple chain below. The
"top"/"bottom" branches apply np.rot90, which swaps the axes -- inverting
that needs its own case, and getting it silently wrong would put points in
the wrong place on a map, which is worse than leaving them out.
```

```text
What the model-facing CSV actually carries, and why it differs.
Defaults live here rather than being added later so every record --
in-span, excluded, no-road, and the control pass -- has the same
fieldnames for write_csv's DictWriter.
```

```text
Same window the topography was extracted with, and the same frame guard the
extractor's run pass applies -- a window picked unstraightened is a valid
index range that points at different cells.
```

```text
road_setback locates the SEAWARD edge of the road block, because
bulldoze indexes [road_start : road_start + road_width]. Ocean-first
cross-shore means the seaward edge is the minimum index.
```

```text
Map coordinates of the two cells the setback is measured between, so
the pair can be plotted in GIS against the road. Diagnostics only --
setback_m above is computed from the indices, never from these.
```

```text
ext.TAG is gone: the arrays carry no year tag since 2026-08-26. The name
comes from the resolver, which is also what the extractor writes with.
```

```text
Past the end of the array: an IndexError mid-run, not a drowning.
Left alone here; the audit script is where that is adjudicated.
```

```text
A move inside the domain's own per-profile spread is a re-pick of the
alongshore statistic. A move beyond it is a position no profile
showed, which is a different claim and has to say so.
```

```text
First row per (domain, transect), then the mean of the transects in each
domain -- the same two steps island_offset_hybrid.py takes, so this and
the build differ only by the build's zeroing.
```

```text
ONE EXTRACTOR PER PRODUCT, both live for the whole run. Loading them up
front means a missing product fails before any file is written, rather
than after the first year has already been published.
```

```text
The picks are loaded up front for the SAME reason, and it is not
hypothetical: they used to be read inside the per-year loop, so a missing
2004 pick set was only discovered after 1984 had already been published
to the shared tree.
```

```text
Measured everywhere, written only within ROAD_SPAN. Flag the excluded
ones in the diagnostics so the gap between "has road" and "is forced"
is on the record.
```

```text
--- seaward relocation of roadways that drown at initialisation ---
Runs on the FLOORED value, because that is what CASCADE would index
with, and fills setback_model_m for every in-span domain.
```

```text
Stations grow LANDWARD from the offshore datum, so 2004 minus the
earlier year is the landward movement between them: >0 = retreat.
The operands were the other way round until 2026-09-22, which made
the printed median the negative of the retreat it named.
```

```text
Does delta_vs_legacy actually behave like retreat? If the legacy file
and this one were measuring the same thing in different years, this
correlation would be strongly POSITIVE. It is not -- see the audit.
```

```text
Does anything actually drown on THIS extraction? The no-data explanation
below is a measurement taken on 2009_v4, not a recomputation, so it must
not be printed under a version where nothing drowns -- it would assert a
cell count that was never measured on the arrays in front of it. The
gap-filled DEM (2009_v5) is exactly that case: 0 drowns, both years.
```

```text
Any vintage the run did not produce. The document is rewritten whole, so
a silent omission here reads as "this forcing does not exist".
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_extractor()`**

```text
Import the corrected extractor, configured for ONE topography product.

ONE RUN, BOTH VINTAGES (2026-08-26). This used to import the extractor once
and read TOPO_PRODUCT off it, so a year whose product did not match was
SKIPPED. That guard was right -- measuring 2004 against 1984-start row 0 is
silently wrong -- but it made a complete run impossible: whichever product
the extractor happened to sit on, the other year was skipped, and
`write_audit` rewrites ONE markdown file for the whole tree. A 1984-only run
on 2026-08-26 14:04 therefore published an audit with no 2004 section at
all, for a forcing that had not changed.

The extractor derives LOAD_PATH, PICKS_DIR, WINDOW_JSON and the save paths
from the TOPO_PRODUCT / VERSION literals at module level, so configuring it
means re-executing it with those two literals substituted. That is what this
does: the source is patched in memory, never on disk, and each product gets
its own module object. The alternative -- mutating TOPO_PRODUCT after import
-- would leave every derived path pointing at the old product, which is the
silent-wrong-row-0 failure again.

`product=None` keeps the file's own literals, which is what run_control's
unstraightened pass and any interactive import get.
```

**`domain_georeference()`**

```text
(transform, crs, n_rows, n_cols) for one domain, or None if unavailable.

Returns None rather than raising: the coordinates are a QC convenience, and a
missing tif or a missing rasterio must never stop a setback being measured.
```

**`cell_to_map()`**

```text
Aligned-frame (profile, cross-shore cell) -> (x, y) in the tif's CRS.

Inverts load_profiles' chain, in reverse order:

    aligned cell k  ->  ocean-first column   k + c0 + shear[profile]
    ocean-first     ->  oriented column      (n_cols - 1) - j
    ALONGSHORE_FLIP ->  original row         (n_rows - 1) - profile
    affine          ->  x, y at the CELL CENTRE (+0.5, +0.5)

The +0.5 is the same convention HAT_check_geojson_vs_mask.py documents: cell k
spans [k, k+1), so its centre is k + 0.5. Dropping it shifts every point half
a cell (5 m) northwest.
```

**`load_saved_interior()`**

```text
The interior array CASCADE actually initialises with.

Deliberately the SAVED .npy, not a re-derived interior: the drown test has
to run on the same bytes `HAT_road_placement_on_domains.py` and
`HAT_road_setback_audit.py` read, or the three scripts can disagree about
which domains drown. Paths come off the extractor rather than being
hardcoded, so bumping VERSION does not silently point this at stale arrays.
```

**`drown_test()`**

```text
bulldoze's width-drown test, transcribed. Returns (seaside, bayside, drowns).

The rows tested are the NEIGHBOURS of the bulldozed band -- road_start - 1
and road_end + 1 -- never the band itself. Interior values are decametres
MHW and bulldoze compares `grid * dz` against the threshold in metres, so
the * 10 is the model's, not a display convenience. Every cell counts,
no-data included: the extractor stores no-data as the water sentinel and
bulldoze reads the literal array.

None means the placement is not testable -- bulldoze indexes road_end + 1
with no bounds check, so that is an IndexError at t=0, not a drowning.
```

**`nearest_viable_seaward()`**

```text
Largest road_start strictly seaward of the current one that does not drown.

Scans landward-to-seaward and returns the first pass, so the road moves the
shortest distance that survives initialisation. None means no row seaward of
here is viable -- the domain keeps its measured setback and drowns.
```

**`relocate_drowning()`**

```text
Fill `setback_model_m` for every in-span domain, moving the drowned ones.

Mutates the records in place. `setback_dunestart_m` is never touched -- the
measurement stays recoverable from the _domains.csv, which is the whole
point of doing this here rather than in the measurement.
```

**`load_stations()`**

```text
Per-domain dune-line station from the SHARED offshore datum, metres,
LANDWARD-positive. 90 rows.

NOT the island-offset build (2026-09-22). This used to read
offset_file(year, "input"), and the retreat below differenced two of
those. Each build is zeroed on its OWN most seaward domain, and for these
two years those minima are 56.6 m apart (1984 zeroed at 1934.3 m, 2004 at
1990.9 m, both on GIS 78), so the difference carried a constant -56.6 m
and came out the wrong SIGN: a true median retreat of +13.8 m was reported
as -42.8 m before the sign convention below, i.e. 43 m of progradation.

A correlation is immune to a constant, so the corr(delta, retreat) result
this function exists to serve never moved -- only the median it printed
beside it, and that median reached RoadOffset_dunestart_audit.md.

The raw files share one offshore datum and have no such constant, so this
reads them, exactly as duneline_endpoint.py does for the same question.
```

**`load_windows()`**

```text
The dune-search windows for the product/version an extractor is set to.

The picks are the one input to this script that cannot be regenerated, and
the path to them is DERIVED, not configured:

    PICK_SET    = RUN_NAME = VERSION
    WINDOW_JSON = PICKS_DIR / f"HAT_dune_search_windows_{PICK_SET}.json"

So any process that bumps VERSION without writing a matching pick file
leaves this script pointing at a path that does not exist. That happened on
2026-08-26: nodata_audit/HAT_bridge_dropouts.py created 1984-start v2 and
moved the VERSION literal, and every re-run here died on a bare
FileNotFoundError naming a file nobody had ever created. The setbacks on
disk were fine -- they simply could not be re-measured.

The bridge script now carries picks forward as part of the bump. This
guard is the backstop for anything else that bumps a version, and it names
the fix rather than the missing file.
```

**`run_control()`**

```text
Re-measure with STRAIGHTEN = False so the frame effect is separable.

One control pass per vintage, each on that vintage's own extractor. The
flag is toggled and restored PER MODULE: there are two module objects now,
so a single saved/restored value would leave the other one unstraightened
for the rest of the process.
```

**`write_audit()`**

```text
The one write-up for the whole tree, covering every year in `audit`.

PROVENANCE IS PER VINTAGE NOW. This took a single `ext` and printed one
DEM path and one picks file for a document describing both years -- true
while a single extraction served every period, false since the tree went
period-first. Worse, the years are no longer guaranteed to be measured in
the same run, and this file is rewritten whole: a run that produced only
1984 published an audit whose 2004 section had simply vanished.

So: one provenance ROW per vintage, and a loud line naming any year that is
absent, rather than a document that quietly describes less than it claims.
```

</details>

### 1-produce/HAT_road_placement_on_domains.py

Where each setback method actually puts NC-12 on the Barrier3D interiors CASCADE starts from.

From the script's original header:

```text
Where each setback method actually puts NC-12 on the Barrier3D interiors CASCADE
initialises with -- the road placed exactly where `roadway_manager.bulldoze`
would put it from that method's model-facing RoadSetback CSV.

  a  1984 road on the 1984-start interiors (product+version resolved at run time)
  b  2004 road on the 2004-start interiors -- a DIFFERENT island, not the same
     one twice: 65 of 90 domains differ in interior shape between the products
  c  setback against island width, one line per period
  d  road movement between the two periods, under this method
  e  bulldoze's own drown test -- does CASCADE keep managing this roadway

ONE SCRIPT, BOTH METHODS, ON PURPOSE
Every method in METHODS is drawn by this same code and each figure lands beside
its own data. Two scripts would drift, and a drifted comparison is worse than no
comparison -- it looks like a result. Add a method by adding a dict entry, not
by copying this file.

  old        archive/superseded_20260911/<year>/RoadSetback_<year>.csv
  dunestart  dunestart_offset/measured/<year>/RoadSetback_<year>_dunestart.csv

Read-only with respect to the forcing: writes a PNG, a PDF beside it and a
CAPTIONS.md entry into each method's folder.

THE STYLE IS THE HOUSE STYLE, AND THIS MODULE RE-EXPORTS IT
Since 2026-09-10 every colour and every type size here comes from
scripts/site_layer/hat_figure_style.py, through apply_style(). The palette names this file
used to own -- LAND_CMAP, SURFACE, WATER, INK_MUTED, INK_SECOND, C_1984 /
C_2004 / C_YEAR, C_DROWN, SECTIONS -- are KEPT and now resolve to their house
equivalents, because HAT_dunestart_modification_stages.py,
HAT_method_comparison_figures.py, HAT_oceanfloor_offset_check.py and
HAT_road_geojson_map.py all read them off this module as `P.<name>`.

THE FRAME CAVEAT -- WHICH APPLIES TO ONE METHOD AND NOT THE OTHER
This figure always draws the setback landward of INTERIOR ROW 0 of that
period's own extraction, because that is the reference CASCADE applies it
against
(`roadway_manager.py:99`):

    road_start = int(road_setback / dy)      rows landward of interior row 0

For `dunestart` that is also the frame the number was MEASURED in, so drawing
and measuring agree.

For `old` they do not: that method measured against the same-year digitised dune
line, a different feature from a different year. That mismatch is the REFERENCE
component isolated in 2-audit/HAT_road_method_diagnostic.py. It is not a reason
to draw the figure differently -- a mis-referenced setback still lands wherever
CASCADE lands it -- but it is the reason the two figures differ, and it should be
read as a property of the method rather than of the island.

HOW THE ROAD IS PLACED  (transcribed from bulldoze, not approximated)
  road_start = int(setback / 10 m)                  truncation, not rounding
  band       = rows road_start .. road_start + 1    20 m, ROAD_WIDTH_CELLS = 2

The same int() is applied to every one of the 50 alongshore profiles in a
domain, so the road is a straight line across an island whose back edge is not.

WHAT THE INTERIOR ARRAYS ARE
domain_<d>_topography.npy, shape (rows, 50), values in DECAMETRES relative
to MHW. -0.30 dam (= -3 m) is the extractor's water sentinel, not an elevation.
Multiplied by 10 for display; masked, not ramped, where it is sentinel -- water
is a state, not a small magnitude.

REQUIREMENTS
  numpy, matplotlib
```

Notes that were in the code:

```text
The topography version is NOT hardcoded here any more. It was "2009_v3", and
when the dune windows were re-picked into 2009_v4 this kept drawing v3
interiors under v4 setbacks without erroring. See hat_topo_version.py.
parents[4] IS scripts/ -- hat_topo_version.py moved there 2026-08-20.
```

```text
The house style. Everything typographic and every colour comes from here now;
nothing in this file re-decides a font size or a grey. The names this module
used to define for its own palette survive below as aliases, because three
other scripts import this file for them.
```

```text
NOR IS THE PRODUCT (2026-08-26). Until today this line read

TOPO_DIR, DUNE_DIR, TOPO_RUN_NAME = topo_dirs()

at module level -- one topography, resolved once, drawn under BOTH panels of
a two-vintage figure. That was correct while a single extraction served every
period. It stopped being correct when the tree went period-first, and the
figure said so in its own caption: "the SAME interiors back both road panels,
since one topography serves every period."

They no longer do. All 90 domains differ between 1984-start and 2004-start
and 65 have a different interior SHAPE, so the 1984 panel was drawing a 1984
setback -- measured from row 0 of the 1984-start extraction -- against a
2004-start island. Same class of error as the v3/v4 one above, and just as
silent: the picture is plausible, it is simply of somewhere else.

Resolved per vintage now, through the one YEAR_PRODUCT mapping.
```

```text
array_name() is the single definition of these filenames - the same one
the extractor writes with. Nothing here spells a name.
```

```text
Each entry produces one figure, beside that method's own data.
setback  the MODEL-FACING file -- the one CASCADE reads, already floored
where the method floors it, because the question is where the model
puts the road, not what the method measured before clamping.
detail   optional per-domain CSV carrying a `flags` column, used only to
mark domains whose true value was negative and got floored to 0.
```

```text
measured/ only: YEARS below are the two measured starts. derived/
holds a copy (2010) and a one-event derivation (1996) of these,
which this figure would only draw twice.
```

```text
--- the drown test, transcribed from roadway_manager.bulldoze --------------
A roadway width-drowns when water cells BORDER it. Three details matter and
none of them is what you would guess from "the road is in water":

* the rows tested are the NEIGHBOURS of the bulldozed band, not the band --
road_start - 1 (seaside) and road_end + 1 (bayside);
* "water" is elevation <= 0 m MHW (drown_threshold=0), NOT the extractor's
-3 m sentinel. Real land that merely sits below MHW counts;
* it fires if EITHER side exceeds the fraction, strictly greater than.

roadway_manager.py:50-51 and 125-165. RoadwayManager re-asserts 0.2 at :531.
```

```text
--- the palette, re-pointed at the house one (2026-09-10) -------------------
These names are the module's public palette: HAT_dunestart_modification_
stages.py, HAT_method_comparison_figures.py, HAT_oceanfloor_offset_check.py
and HAT_road_geojson_map.py all read them off this module as `P.<name>`. They
are kept, and each one now RESOLVES to its house equivalent, so the figures in
this folder and the figures elsewhere in the project are one palette.

The earlier vintage is the house red and the later the house blue, which is
the vintage pair everywhere in this project; the pair used here before
(#2a78d6 blue 1984 / #eb6834 orange 2004) had 1984 blue, i.e. exactly
inverted against every other two-vintage figure in the repo.
```

```text
The drowned state. It was a dark crimson picked to separate from the old
blue/orange pair; against the house vintage RED it would now read as "1984".
ACCENT purple is the one house colour that is neither vintage, neither
reference nor fabricated ground, and it separates from both poles on hue and
from BASE on luminance. Drowned domains also carry a marker and a count, so
the state is never colour-alone.
```

```text
INK_SECOND is used for TEXT by the importers, INK_MUTED for rules, outlines
and grid; the house rule puts text in INK and rules in INK_MUTED.
```

```text
Elevation is drawn in CLASSES, not a ramp (house rule): a linear ramp over
0-4 m renders the whole back-barrier as one tone. LAND_CLASS_* is what the
panels here use. LAND_CMAP survives as a CONTINUOUS ramp in the same house
colours, because the two out-of-scope importers pair it with a plain
Normalize(LAND_VMIN, LAND_VMAX); handing them the discrete list under a
linear norm would paint 0-0.5 m land in the water colour.
```

```text
The alongshore reaches. Kept because HAT_oceanfloor_offset_check.py reads
them off this module to draw its own dividers. The figures BELOW no longer
use them: the house `town_bands()` shades the three village spans from
hatteras_site_config, which is the one authority on where the towns are.
```

```text
bulldoze indexes road_end + 1 with no bounds check, so a road this
far back is not "drowned", it is an IndexError at t=0.
```

```text
Only road_start is drawn. The 20 m band's landward edge sat 2 cells away,
which at this vertical scale read as line weight, not as width.

Segments are split by drown state so a failing domain is crimson in place,
rather than being annotated off to one side -- the question "where does
this road fail" is answered on the island, not in a caption.
```

```text
No marker on the plan views: the accent segment IS the signal, and a
triangle on top of a 1-domain-wide line obscured the thing it pointed at.
The state stays recoverable without colour -- panel (e) names every
failing domain on the same x-axis, and the caption carries the count.
```

```text
The three village spans, from hatteras_site_config -- a named strip
against the landward edge rather than a full-height wash, which on an
image panel would hide the island it is meant to locate.
```

```text
Each vintage is placed on ITS OWN interiors. Passing one shared dict
here is the bug this signature exists to make impossible.
```

```text
Anything cropped out of the plan view that is actually road would make the
picture a lie, so check rather than assume.
```

```text
--- (c) setback against the island it has to fit inside ----------------
ONE BAND PER VINTAGE. This was a single shaded band, drawn from the one
shared interiors dict, with both setback curves over it -- which read as
"here is the island, here are two roads on it". There are two islands.
The 1984 and 2004 widths differ by up to 80 m in places, and a setback
that fits one can run off the other.
```

```text
Lines, not a shaded band under each width. Two 10%-alpha bands, one
red and one blue, overprint to a purple wash across the whole panel
-- and purple is the drowned colour everywhere else in this figure.
```

```text
--- (d) where the road moved between the two periods -------------------
Ported from the retired 3-figures/island_wide/HAT_plot_road_on_b3d_domains
.py, which drew this for one method only.
```

```text
--- (e) what the bulldozed band actually lands on ----------------------
The series is a PERCENTAGE, so the threshold has to be scaled too -- at
DROWN_PCT it would sit on 0.2% and read as zero.
```

```text
ASCII hyphen, not U+2212: the Windows console is cp1252 and a
unicode minus raises UnicodeEncodeError here. Figure text is fine.
```

```text
The comparison the two figures exist to support, stated as numbers so it
does not have to be eyeballed off two PNGs.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_interiors()`**

```text
The Barrier3D interiors THIS vintage's setbacks were measured against.

`year` is required. It used to be absent, and a caller that cannot name a
year is a caller drawing two vintages on one island -- see topo_for_year.
```

**`load_years()`**

```text
{year: {interiors, canvas, shown}}, plus the crop shared by all panels.

The panels are stacked and share an x axis, so they must share a y extent
too -- otherwise a taller 1984 island would read as a wider barrier rather
than a taller canvas. The crop is therefore the max over vintages, and a
shorter vintage's canvas is NaN-padded up to it rather than drawn shorter.
```

**`_wet_fraction()`**

```text
Fraction of a border row at or below the threshold -- EVERY cell counted.

A no-data mask (domain_<d>_nodata.npy) exists beside the topography and
was briefly used here to drop never-surveyed cells from both the numerator
and the denominator. That is not the test bulldoze runs. bulldoze reads the
literal array, where the extractor has already written no-data back as the
water sentinel, so an unsurveyed cell IS a water cell as far as the model is
concerned -- and a road CASCADE would stop managing has to show as one here.
Excluding no-data flipped GIS 14/77/78/79/80 to pass while the model still
drowns them. Deliberately reverted: the figure reports what the run does.
```

**`place_road()`**

```text
bulldoze's own placement, and bulldoze's own drown test, per domain.

The drown test is transcribed rather than approximated -- it decides whether
CASCADE gives up managing this roadway, so an approximation of it would be a
figure about a model we are not running. Interior values are decametres MHW
and the model compares `grid * dz` against the threshold in metres, so the
* 10 here is bulldoze's, not a display convenience.
```

</details>

### 1-produce/HAT_road_setback_derived_vintages.py

Road setback files for the 1996 and 2010 starts, derived rather than measured.

Notes that were in the code:

```text
HAT_road_setback_derived_vintages.py

A road setback file for each of the two hindcast periods added 2026-09-11,
1996-2010 and 2010-2024.

NEITHER IS A MEASUREMENT. No NC-12 line was digitised for either vintage --
the repo holds two road lines, exported for 1978 and 2008 and labelled 1984
and 2004 (see "Road line vintages" in hatteras_site_config.py). So each new
file is DERIVED from a measured one, and the derivation lives here rather
than in a hand-edited CSV, because a hand-edit leaves no record of what was
added to what.

1996  =  RoadSetback_1984_dunestart.csv  +  the 1989 Pea Island relocation

The 1989 relocation (GIS 84-87) has ALREADY HAPPENED by 1996 and the
1999 one has not, so the 1996 road is the 1984 road with exactly one
event applied. The displacements are not retyped here: they are read
from HATTERAS_ROAD_EVENTS, so this file and the event the model fires
in a 1984-start run cannot disagree about how far the road moved.

WHAT THIS INHERITS. The 1984 setbacks are measured against interior
row 0 of the 1984-start extraction, so these are too, and the 1996
period must read that same product. Everything the 1978-vintage line
gets wrong about 1984 it also gets wrong about 1996, eighteen years
later rather than six.

2010  =  RoadSetback_2004_dunestart.csv, unchanged

The 2010 period reads the same topography product as the 2004 period
and the same road line, and no relocation in the record falls between
the two dates -- the next event is the 2022 bridge, which is inside the
run rather than at year zero. So the measurement is identical and this
is a copy under the period's own name, written rather than referenced
so every period names its own file (Hannah, 2026-09-11).

OUTPUTS
dunestart_offset/derived/1996/RoadSetback_1996_dunestart.csv
dunestart_offset/derived/2010/RoadSetback_2010_dunestart.csv
dunestart_offset/derived/<year>/PROVENANCE.md   what was derived from what

derived/ rather than beside the measured folders (2026-09-15): the address
says these are not measurements. The sources are read from
dunestart_offset/measured/, both through hat_topo_version.road_setback_file().

Nothing existing is read for writing, and both outputs are new files.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
```

<details><summary>Function notes (the original docstrings)</summary>

**`setback_path()`**

```text
measured/<year>/ for the sources, derived/<year>/ for the outputs --
hat_topo_version.ROAD_SETBACK_KIND decides which, not this script.
```

**`event_displacements()`**

```text
The displacement dict of the relocation event dated `year`.

Read from HATTERAS_ROAD_EVENTS rather than from the measurement CSV
directly, so this applies the SAME rounded, whole-cell displacement the
model applies when the event fires.
```

</details>

### 2-audit/HAT_check_geojson_vs_mask.py

Does the rasterized road mask sit where the source geojson says the road is?

Notes that were in the code:

```text
HAT_check_geojson_vs_mask.py

Does the rasterized road mask actually sit where the source geojson says the
road is? Read-only QC on the rasterization step alone.

WHY THIS IS SEPARATE FROM THE PLACEMENT FIGURE
HAT_method_comparison_figures.py scores modelled-vs-real road position,
which folds together three things: rasterization, the alongshore collapse to
one scalar, and cell quantization. This script isolates the FIRST one, and it
does so in the ORIGINAL raster frame -- no orientation, no alongshore flip,
no shear, no trim -- so it is independent of every transform the extractor
applies. If this passes, a registration error downstream is not the
rasterizer's fault.

WHAT IT CHECKS
1. CONTAINMENT  in each profile the geojson crosses, does it fall inside the
mask's cell footprint? The mask is a 6 m buffer of that same line with
all_touched, so containment should be essentially universal. Anything else
means a CRS, snap-raster or extent mismatch.
2. OFFSET       signed distance from the centerline to the mask's centre, in
metres. Should sit near zero with a spread of well under a cell.
3. COVERAGE     profiles where one source has road and the other does not.
4. YEAR IDENTITY  do the 1978 and 2008 masks actually differ, and only where
NC-12 was relocated? This is the check that the 2004 rasterization used
the 2004 line -- a patched driver could silently re-burn 1984.

USAGE
python HAT_check_geojson_vs_mask.py
```

```text
Anchored 2026-09-14: absolute into a home directory, or into a tree
renamed since. Rule 5 of ORGANIZATION.md.
```

```text
MOVED 2026-08-25: 1-barrier3d-domains went period-first and the pre-90-domain
legacy went under superseded/. Path repointed so this script keeps reading
EXACTLY what it read before - no road number moves because of the reorg.
```

```text
Documented relocation blocks: the 1999 Buxton-Avon work and the 1989 Pea
Island move. Used only to interpret check 4, never to gate it.
```

```text
densify so a long segment crossing the domain still lands in every
row it passes through
```

```text
A MISSING TIF USED TO BE INDISTINGUISHABLE FROM "no road here".
On 2026-08-26 domain-clips-1m moved out of superseded/ and TIF_FMT
was not updated with it. Every domain took this branch, so the check
reported "0 profiles with both line and mask" for BOTH years and a
median offset of nan -- and then carried on to Check 4 and printed a
result. This is the only guard against re-burning the wrong year;
a guard that passes vacuously is worse than no guard, so the tif
being absent is now counted and reported rather than skipped.
```

```text
`c` is a FRACTIONAL column from the inverse affine, where cell k
occupies [k, k+1). A cell's centre is therefore k + 0.5, so the
mask's centre is cells.mean() + 0.5. Comparing against
cells.mean() alone injects a spurious +0.5 cell (+5 m) bias into
every profile.
```

```text
Magnitude is what separates a real relocation from digitizing noise: a
moved road changes hundreds of cells, a redrawn line changes a handful.
```

```text
A VACUOUS PASS IS NOT A PASS. When domain-clips-1m moved out of
superseded/ on 2026-08-26 and TIF_FMT was not updated, every domain
was skipped for want of a tif: this printed "0 profiles", a median
offset of nan, then carried on to Check 4 and produced a report. The
README calls this the only guard against re-burning the wrong year,
so it has to fail loudly when it has compared nothing.
```

<details><summary>Function notes (the original docstrings)</summary>

**`geojson_cols_by_row()`**

```text
Centerline position per raster row, in ORIGINAL (row, col) cell units.

The geojson is EPSG:2264 in US survey feet; the domains are UTM 18N metres,
so this reprojects first. Vertices are binned to integer rows and averaged,
which is what a near-north-south line through a north-up grid needs.
```

</details>

### 2-audit/HAT_road_buffer_bias.py

How far the setback's buffered-mask edge sits from the real road, and whether it matters.

Notes that were in the code:

```text
HAT_road_buffer_bias.py

The setback is measured to the SEAWARD EDGE OF A BUFFERED MASK, not to the road.
How far apart are those, and does it matter?

WHY THIS EXISTS
HAT_rasterize_road_to_domains.py burns NC-12 with ROAD_BUFFER_M = 6.0 and
ALL_TOUCHED = True, giving a ~24 m mask for an ~8 m road.
HAT_road_offset_from_dune_start.py then takes that mask's seaward-most cell as
`road_start`, because `bulldoze` indexes [road_start : road_start + width].
So every setback inherits however far the buffer pushed that edge seaward, and
nothing else in this tree measures it -- the geojson is used to rasterize, to
check registration (HAT_check_geojson_vs_mask.py) and to draw maps, never to
produce a setback.

WHAT IT MEASURES
bias = (mask seaward cell) - (geojson centerline position), in metres,
NEGATIVE meaning the mask edge sits SEAWARD of the true centerline and the
reported setback is therefore SMALLER than a centerline measurement.

Measured in the ORIGINAL raster frame. `setback = (seaward - row0) * cell` and
row0 is identical either way, so orient / flip / shear / trim all cancel out of
the difference -- no transform is needed, and the result is independent of the
extractor.

WHAT IT FOUND (2026-08-19, both vintages)
median bias -6.8 m, p10 -10.8, p90 -2.7.

That is NOT a 6.8 m placement error, because of what it is compared against.
`road_width_m = 20.0` and bulldoze lays a 2-cell block LANDWARD from
road_start, so a block centred on the real road wants
road_start = centerline - 10 m. The buffer supplies centerline - 6.8 m. The
residual misplacement is 3.2 m -- about a third of a cell, and an order of
magnitude below the 20-40 m p90 placement error that the alongshore collapse
to one scalar already causes (HAT_method_comparison_figures.py).

DO NOT "CORRECT" IT. Per-profile setbacks are integer cell differences times
10 m, so domain medians land almost exactly on cell boundaries, and bulldoze
truncates with int(). Applying the 3.2 m adjustment moves road_start a FULL
cell (10 m) seaward in 83% of 1984 domains and 90% of 2004 domains, replacing
a 3 m error with a 7 m error in the other direction. The buffer is accidentally
doing very nearly the right job; the quantization is what would bite.

Cross-validated independently: the road_x/road_y columns
HAT_road_offset_from_dune_start.py writes into RoadOffset_<year>_profiles.csv
are produced by a completely different path (full index inversion to map
coordinates), and their shapely distance to the same geojson is 6.62-6.68 m
median -- the same number to ~0.2 m.

D8 is the one extreme outlier (-353 m in 1984). That is the Buxton bend, where
NC-12 runs parallel to the raster rows, so a "seaward-most cell" is not a
meaningful quantity -- the same reason D8 is already EXCLUDED_FROM_SPAN. Its
appearance here is a consistency check passing, not a new problem.

USAGE
python HAT_road_buffer_bias.py

REQUIREMENTS
numpy, geopandas, rasterio   (same set HAT_check_geojson_vs_mask.py needs)
```

```text
The centerline-to-raster-column projection is already solved next door, half-cell
convention and all. Importing it rather than re-deriving it keeps ONE definition
of "where is the centerline in this raster".
```

```text
bulldoze's modelled road block, from cascade_pipeline/roadway.py (road_width_m)
and roadway_manager.bulldoze: road_end = road_start + int(road_width / dx).
```

```text
Domains excluded from the model-facing span, reported separately rather than
dropped. Must match ROAD_SPAN in HAT_road_offset_from_dune_start.py.
```

```text
Ocean is at the RIGHT in the raw frame (OCEAN_LOC = "right"), so
ocean-first index 0 is the LAST raw column and the SEAWARD-most mask
cell is the MAXIMUM raw column.
```

```text
`c` is a fractional column from the inverse affine, where cell k
spans [k, k+1); a point therefore sits at cell-centre units c - 0.5.
Comparing against the raw integer injects half a cell (5 m).
```

### 2-audit/HAT_road_setback_audit.py

Audit the NC-12 setbacks, initial and after relocation, against the interiors CASCADE runs on.

From the script's original header:

```text
Audit the NC-12 setbacks -- initial AND prescribed-relocation -- against the
interior arrays CASCADE will actually run on, and write down where the road
lands off the island and what the model does about it.

This script is DIAGNOSTIC. It never modifies a setback, never writes a setback
file, and never touches the topography.

WHY AN AUDIT AND NOT A FIX
The setback MEASUREMENT was tested and holds up. Road and dune are digitised as
points along the SAME 1 m transects (FID_Transect_Points_1m / LineID /
ORIG_SEQ), from the same origin, in the same units, for the SAME YEAR -- so the
shoreline retreat cancels in the subtraction. road_offset_pipeline.py's one
shortcut is taking min(ORIG_LEN) over each layer INDEPENDENTLY, so the two
minima land on different transects. Measured on 2004 (the only vintage whose
raw road CSV survives), pairing road and dune on the SAME transect changes the
answer by a median of +1 m (worst -49/+54; |delta| > 25 m in 14 of 82 domains,
> 100 m in none). One grid cell. Not worth re-deriving.

WHAT ACTUALLY WENT WRONG IS THE FIT CHECK
HAT_setback_from_lines.py (retired 2026-08-17, git blob c37ec03e) bounds the
setback with

    interior_rows_of(d) -> np.load(...).shape[0]

which is Barrier3D's DomainWidth: the ROW COUNT of the array. Barrier3D never
uses that to decide where land is. barrier3d.py FindWidths() does:

    width = next((i for i, v in enumerate(InteriorDomain[:, bl])
                  if v <= SL), DomainWidth) - 1

-- per profile, the contiguous run of land from row 0 to the FIRST cell at or
below sea level. remove_water_rows() only trims LEADING and TRAILING all-water
rows, so one wide profile sets DomainWidth for all 50 and everything behind the
other 49 is sentinel padding. Median gap between the two numbers across GIS
9-90: 17 rows. Worst: 86 rows -- 860 m of setback the old check would approve
onto cells that are not there.

find_widths() below is FindWidths transcribed verbatim. The definition is
adopted, not invented; nothing in Barrier3D changes.

WHAT "DROWNS" ACTUALLY COSTS -- READ THIS BEFORE READING THE TABLES
roadway_drown is not a warning. Traced through the code, one drown does this:

  1. bulldoze() MUTATES barrier3d.InteriorDomain IN PLACE before it checks
     anything: xyz_interior_grid[road_start:road_end, :] = new_road_domain.
     The grid is passed by reference, so the write lands even though the
     function returns early afterwards. On GIS 52 at the 1984 setback, 53 of
     100 cells in the road rows were water (-0.3 dam) and all 100 came out at
     0.145 dam -- the model gains a 20 m ribbon of 1.45 m land across open
     water.

  2. RoadwayManager sets _drown_break = 1 and returns.

  3. cascade.py (~line 625) sees drown_break on EVERY later year and never
     calls update() again. _road_break[iB3D] = 1, dune growth rates reset to
     natural.

So from that year on there is no overwash removal, no dune rebuilding, no road
relocation, and _road_overwash_volume / _dunes_rebuilt_TS /
_rebuild_dune_volume_TS stay at zero for the rest of the run. A domain flagged
DROWNS_T0 is an UNMANAGED barrier wearing a road label for the whole hindcast.
If any result contrasts managed against unmanaged shoreline, those domains are
silently in the unmanaged group from year zero.

(group_roadway_abandonment is None in the runner, so this stays per-domain. Set
it, and one drowned domain abandons its entire group.)

THE ROAD'S OWN CELLS ARE NEVER CHECKED
bulldoze() tests the rows FLANKING the road -- road_end + 1 on the bay side,
road_start - 1 on the sea side -- and never looks at the road footprint itself.
It also skips row road_end: the road occupies [road_start:road_end], and the
bayside test is at road_end + 1, so the cell immediately behind the asphalt is
only used for np.size(). This audit reports pct_road_cells_water separately for
exactly that reason: a road can be substantially in the bay while the flanking
test is still under threshold.

PRESCRIBED RELOCATIONS BYPASS BOTH OF THE MODEL'S GUARDS
CASCADE has two guards on moving a road landward:

    get_road_relocation_elevation()  rebuilds the road at grade and refuses if
                                     mean(road_domain) <= 0 m MSL
    road_relocation_checks()         refuses if
                                     setback + 2*width > average_barrier_width
                                     (InteriorWidth_Avg, i.e. FindWidths)

Both live on the MODEL-DRIVEN path, taken when the dune migrates over the road.
The runner's HISTORICAL_ROAD_EVENTS path does not go through either -- it
assigns rm._road_setback = new_sb directly, leaves road_ele at its
SLR-decremented value rather than rebuilding at grade, and lets the next
bulldoze() flatten whatever is there. This audit replicates both guards on the
prescribed setbacks so the report can say which the model would have refused.

INPUTS
  SCENARIOS below: 2-row RoadSetback_<year>.csv files and/or literal
  {domain: setback_m} maps taken from HISTORICAL_ROAD_EVENTS
  TOPO_DIR/domain_<N>_topography.npy

OUTPUTS
  RoadSetback_audit.csv   every scenario x domain, machine readable
  RoadSetback_audit.md    the tracking document
  console verdict, exit code 1 if any scenario hits a hard wall

REQUIREMENTS
  numpy
```

Notes that were in the code:

```text
The method the runner spends (hatteras_site_config.py:78,91). Switched to
dune-start on 2026-08-18; the legacy tree is still on disk under
road_offset/archive/superseded_20260911/ for the method comparison.
```

```text
The runner uses ONE topography for BOTH hindcast periods, so the 1984 setbacks
and the 2004 setbacks are both spent against a 2009 grid. For 1984 that is a
25-year anachronism; it is the reason this audit exists.

This USED to say "must match TOPO_DUNE_VERSION in the runner", and it was
pinned to 2009_v2, then to 2009_v3. Pinning by hand is what went wrong: when
the dune windows were re-picked into 2009_v4 on 2026-08-19 this kept auditing
v3 interiors against v4 setbacks and said PASS, on a grid the model no longer
uses. 18 domains had different interiors and 10 a different SHAPE -- including
D79 and D80, two of the three roadways the relocation logic acts on.

So the version now comes from the extractor that PRODUCED the arrays, and a
version missing from disk is a loud error rather than a stale read. To audit
an older run deliberately, pass override= (see hat_topo_version.py).

The runner (HAT_hindcast_1984_2024.ipynb, HAT_hindcast_1984_2024.py,
HAT_groin_sweep_worker.py) still says 2009_v2, a path that does not exist, so
the mismatch with the runner remains REAL and is printed below rather than
hidden. Repointing the runner is a separate job.
parents[4] IS scripts/ -- hat_topo_version.py moved there 2026-08-20.
```

```text
--- the relocation bounds, DERIVED (2026-08-26) ----------------------------

These were hardcoded literals - {9: 73.0, 10: 97.0, ...} - described as
"1984 setback + 1978->1997 displacement". They were computed from a 1984
setback file that has since been re-measured against 1984-start/v1, so they
silently stopped meaning what their own note said. Worse, they drove a
DROWNS verdict.

They are now computed, from the two sources the MODEL uses:

* the current model-facing 1984 setback, RoadSetback_1984_dunestart.csv
* HATTERAS_ROAD_EVENTS[].displacement_m, which hatteras_site_config reads
from road_relocation_1978_2008.csv (mean_signed_landward_m) and which
the runner applies as new_setback = rm._road_setback + displacement

Evaluated at ZERO modelled retreat, so this stays the upper bound / worst
case it always claimed to be: any real retreat moves the road seaward of what
is audited here, and a domain clean at the bound is clean for every retreat.

THE EVENT DOMAIN LISTS COME FROM THE CONFIG, NOT FROM HERE. The 1999 event
covers GIS 9-14, not 9-15: hatteras_site_config excludes GIS 15 because the
two digitised lines CROSS inside it, so its mean and median displacement
disagree in sign. The old scenario audited 15 at a literal 106.0 m, a number
the measurement cannot support. Deriving the list removes that by
construction.
```

```text
EACH SCENARIO NAMES ITS OWN TOPOGRAPHY PRODUCT (2026-08-26).

This used to be a single module-level topo_dirs() with no argument, i.e.
DEFAULT_PRODUCT, and EVERY scenario was audited against it. Since the tree
went period-first that is wrong by construction: the 1984 setbacks are
measured from row 0 of the 1984-start extraction and the 2004 setbacks from
row 0 of 2004-start, so auditing both against one product checks at least one
of them against interiors the model will never see.

That is the SAME failure this file's header describes at v3/v4 - "18 domains
had different interiors and 10 a different SHAPE" - just with product in
place of version. Resolved per scenario, and both are printed and written to
the CSV so the pairing is on the record rather than assumed.
```

```text
Resolved from the product, not a literal. The tag became period-specific
on 2026-08-26 and was reverted the same day (see hat_topo_version.py); with
argument resolves the same DEFAULT_PRODUCT topo_dirs() just used, so the
directory and the filename can never disagree about which period this is.
array_name() is the single definition of these filenames - the same one
the extractor writes with. Nothing here spells a name.
```

```text
There is no longer a runner version to disagree with. Until 2026-08-20 the
runner carried its own hardcoded TOPO_DUNE_VERSION and this file kept a copy
of it so a mismatch could be warned about -- a copy that itself went stale
(it said 2009_v2 while the runner said 2009_v3 and the setbacks were built on
v5). Both now call topo_dirs(), so they resolve to the same extraction by
construction and the warning has nothing left to check.
The resolved product/version pairs are printed by main(), after SCENARIOS
exists. There is no module-level banner any more: there is no longer a single
topography to name.
```

```text
--- what to audit ------------------------------------------------------
"initial"    a 2-row RoadSetback_<year>.csv, the whole GIS 9-90 window
"relocation" a literal {gis: setback_m} lifted from HISTORICAL_ROAD_EVENTS in
HAT_hindcast_1984_2024.py. Relocation scenarios also
get the two model guards replicated.
```

```text
THE TWO DERIVED VINTAGES, added 2026-09-11 with the 1996-2010 and
2010-2024 periods. Neither was measured from a road line of its own, so
auditing the file the model actually loads matters more here, not less.

1996  the 1984 file with the 1989 Pea Island displacement applied. Its
other 78 domains are the 1984 values on the same grid, so the new
information is GIS 84-87 -- which the "1989 relocation" scenario
below already audits at the same numbers, by construction. Both
are kept: that scenario audits an EVENT, this audits a FILE, and
a derived file that stopped matching its own derivation is
exactly the failure worth catching.
2010  byte-identical to the 2004 file on the same product, so it is
NOT audited separately -- a second scenario would report the same
82 rows under a different label and invite them to be read as
independent agreement.
```

```text
The relocation events now carry a DISPLACEMENT, applied to whatever
setback the model is carrying at the event year:

new_sb = rm._road_setback + relocation_displacement_m[gis]

so the audited value depends on how much the model has retreated by then.
The setbacks below are that expression evaluated at ZERO modelled retreat,
i.e. 1984 setback + displacement -- which is exactly the old absolute
post_relocation_setback_m. So this scenario is the WORST CASE / upper
bound of the corrected form, and any real retreat moves the road seaward
of what is reported here. Audited at the bound deliberately: if a domain
is clean here it is clean for every retreat.
```

```text
--- Barrier3D / roadway_manager constants ------------------------------
Not tuning knobs. Every one read from model source:
ROAD_WIDTH_M   HAT_hindcast_1984_2024.py:2206
DX/DY/DZ       roadway_manager.bulldoze() defaults, dam
DROWN_THR      roadway_manager.py:766   drown_threshold = 0
PCT_WATER      roadway_manager.py:531   percent_water_cells_touching_road
SL             barrier3d FindWidths sea level, 0.0 dam at t = 0
```

```text
(GIS, old absolute post_relocation_setback_m, RoadSetback_2004.csv value).
The evidence that the old absolute values double-counted the retreat: both
events precede 2004, so the 2004 same-year measurement already IS the
post-relocation position. Reported in the markdown; nothing computes from it.
```

```text
--- the flanking rows, which is all bulldoze tests ----------------
Every cell counts, no-data included -- bulldoze reads the literal array
and the extractor stores no-data as the water sentinel. See wet_fraction.
```

```text
Barrier3D's two guards do not use the same notion of "island", and in a
holed domain they disagree:
FindWidths stops at the FIRST water cell, so land behind a gap is not
part of the island. Feeds InteriorWidth_Avg -> average_barrier_width,
which road_relocation_checks() uses to decide if there is room.
bulldoze reads the LITERAL cell at road_end + 1; land behind a gap
counts perfectly well.
So the road can be founded on ground the relocation logic believes is off
the back of the barrier. Counted, not corrected.
```

```text
No glob fallback. It existed to tolerate whatever year tag the arrays
happened to carry; there is no tag now, and a glob that matched more than
one file would have picked arbitrarily. A missing file is a missing file.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_find_project_root()`**

```text
Walk up until a directory holds data\hatteras_init.

NOT parents[N]. This file has moved, and the old parents[4] resolved to
scripts\ instead of the project root -- so every path below it was wrong
and nothing raised until a glob came back empty.
```

**`find_widths()`**

```text
barrier3d.py FindWidths(), verbatim -- including the `- 1` and the clamp at
zero, so it cannot drift from the model.

InteriorWidth[bl] is the contiguous run of land from row 0 of profile bl to
the first cell at or below sea level. It is the only width Barrier3D uses
to decide where the island is; land behind a water cell is invisible to it.
```

**`land_behind_first_water()`**

```text
Per profile: land cells BEHIND the first water cell, and the width of the
intervening gap.

Separates 'the island ends here' from 'FindWidths stopped at a 2-cell hole
and discarded 33 cells of barrier'. Decides nothing -- with no bathymetry
the two are indistinguishable by value -- but marks which domains carry the
doubt.
```

**`predict_bulldoze()`**

```text
Reproduce roadway_manager.bulldoze()'s indexing, its drown test, and what
it does to the road's own footprint -- without running it.

`wall` is None unless the setback would CORRUPT the run rather than merely
drown the road:

    NEGATIVE  road_start = int(-110/10) = -11, and numpy WRAPS a negative
              index to the back of the array. The road is bulldozed onto
              the sound side and the run finishes looking normal.
    OVERRUN   xyz_interior_grid[road_end + 1, :] is past the array ->
              IndexError mid-run.

Datum note: bulldoze compares (interior * dz) against drown_threshold = 0
described as "m MSL", but the extractor's arrays are MHW-RELATIVE (MHW_M =
0.36 subtracted first). The effective test is "at or below MHW". A
pre-existing inconsistency in the model -- reported, not fixed.
```

**`replicate_relocation_guards()`**

```text
What CASCADE's own relocation guards would say about this setback.

get_road_relocation_elevation(): rebuilds the road at grade and refuses if
    mean(road_domain) * dz <= 0 -- "Roadway cannot be relocated ... b/c the
    road would be at or below MSL".

road_relocation_checks(): refuses if
    road_relocation_setback + 2 * road_relocation_width
        > average_barrier_width
where average_barrier_width is barrier3d.InteriorWidth_AvgTS[-1] * 10,
i.e. the MEAN of FindWidths' InteriorWidth, in metres.

Both are on the model-driven relocation path. The runner's
HISTORICAL_ROAD_EVENTS path assigns rm._road_setback directly and takes
neither. Replicated here on the t=0 interior; the real 1989/1999 grid will
have evolved, so treat these as the initialisation-time verdict.
```

**`largest_setback_that_fits()`**

```text
Largest setback that neither drowns nor overruns. REFERENCE ONLY -- never
applied to a setback. Reported so the document can say how far the aerials'
position is from anything the modelled barrier could hold.
```

**`load_setback_csv()`**

```text
Read a 2-row CASCADE file and verify the ID row is what the runner assumes.

The runner does np.loadtxt(..., skiprows=1) and fills road_setbacks_full[]
BY POSITION, discarding the IDs. A missing, duplicated or out-of-order ID
silently shifts every domain north of it. Checked here because this is the
only place it can be.
```

**`wet_fraction()`**

```text
Fraction of a flanking row at or below the threshold -- EVERY cell counted.

A no-data mask (domain_<d>_nodata_<TAG>.npy) exists beside the topography
and was briefly used to drop never-surveyed cells from both the numerator
and the denominator. That is not the test bulldoze runs. bulldoze reads the
literal array, in which the extractor has already written no-data back as
the water sentinel, so an unsurveyed cell IS a water cell to the model.
Excluding it made GIS 78/79/80 read FITS while the run still drowns them --
exactly the silent pass this audit exists to catch. Deliberately reverted.
```

</details>

### 3-figures/HAT_dunestart_modification_stages.py

The dune-start setback before each of its two modifications, drawn stage by stage.

From the script's original header:

```text
The dune-start setback carries TWO modifications on top of the measurement, one
at each edge of the island. This script draws the island BEFORE each of them, so
the progression can be read stage by stage instead of taken on trust.

    setback_dunestart_m           raw measurement            1984: 1 negative, 1 drown
       |  OCEAN-SIDE MOVE -- negative setbacks floored to interior row 0
    setback_dunestart_floored_m   ocean-side applied         1984: 0 negative, 0 drown
       |  BAY-SIDE MOVE -- roadways drowning on a wet bayside row moved seaward
    setback_model_m          both applied, MODEL-FACING 1984: 0 negative, 0 drown

  Counts are from the run: 1984 on 1984-start/v1, 2004 on 2004-start/v1, and
  2004 is 0 negative / 0 drown at every stage. They have moved twice. The
  gap-filled DEM took the bayside drowns that the BAY-SIDE move existed for to
  zero, leaving every stage-0 drown a negative being tested on its wrapped row;
  then the 1984 vintage was re-measured on its OWN product (2026-08-26) and the
  1984 negatives went 6 -> 1, GIS 85 alone. The bay-side move is therefore a
  no-op under both topographies rather than a step with work to do -- a property
  of the DEMs, so the figure counts it from the data rather than asserting it.

  stage 0  HAT_dunestart_stage0_raw.png          before BOTH moves
  stage 1  HAT_dunestart_stage1_ocean_floor.png  ocean-side only, before the
                                                 bay-side move
  stage 2  ../HAT_dunestart_road_on_domains.png  the existing figure, both
                                                 applied -- not redrawn here

HOW FAR each ocean-side move actually is, cropped to the two stretches where it
fires, is HAT_oceanfloor_offset_check.py -> HAT_dunestart_oceanfloor_check.png.
At 90 domains these stage figures cannot show a 10 m move; that one can.

WHY THIS IMPORTS RATHER THAN COPIES
The drown test, the interiors, the canvas and the palette all come from
HAT_road_placement_on_domains.py by import. The whole value of a stage figure is
that stage 2 is the SAME test as stages 0 and 1; a transcribed copy that drifted
by one row would make the progression a fiction. Nothing about the test is
re-implemented here -- only the negative-setback case, which stage 2 cannot
contain by construction.

HOW A NEGATIVE SETBACK IS DRAWN  (stage 0 only)
A negative setback has two positions, and only ONE of them is drawn on the plan
view:

  TRUE      where the road was measured -- seaward of interior row 0, out in the
            dune/beach that the interior array does not cover. Drawn on a sand
            band below the island, on a y-axis extended past 0. This is the only
            position the plan view shows.
  WRAPPED   where CASCADE would actually put it. `int(-70/10) = -7` and
            `xyz_interior_grid[-7:-5, :]` is valid Python indexing from the
            LANDWARD end, so the road is bulldozed into the bay with no error
            raised. NOT drawn on the plan view -- a second mark for one road,
            in a place no measurement supports, reads as two roads rather than
            as one road and its consequence.

The wrap is still what decides those domains' drown state, so it is not lost:
panel C reports their percentages from the WRAPPED rows -- that is what the
model would test -- and marks them apart so the number is never read as a
measurement of the true position. The header text says so on the figure itself.

REQUIREMENTS
  numpy, matplotlib
```

Notes that were in the code:

```text
The house style, through the module that already resolves it. P.apply_style()
has run at import, so this file only needs the helpers.
```

```text
Sand, for the strip seaward of interior row 0 that the interior array does not
cover. Deliberately not the water colour -- a road measured out there is on
the beach, not in the sound, and the two must not read the same. The house
ADDED_FILL is that sand: it means ground that is not in the surveyed array.
```

```text
`png` is fixed: other files in this tree reference these names. The
LABELS are not -- "raw" was working vocabulary and says nothing about
what was or was not done to the number.
```

```text
Exactly what numpy does with these indices -- negative subscripts
count from the landward end. Not a simulation of it.
```

```text
The strip seaward of interior row 0. Only drawn when a road is out there,
so stage 1 keeps the same axes as stage 2 and the two stack cleanly.
Tested against -CELL/2, not 0: that is where the interior image starts, so
a stage with no negatives has floor_m == -CELL/2 and must draw no band.
```

```text
No in-band caption: any text here sits on top of the drowned roads
this band exists to show. The figure legend names the colour instead.
```

```text
The wrapped position is NOT drawn on the plan view. It was, and it put a
second mark for one road in a place no measurement supports, which read as
two roads rather than one road and its consequence. The wrap still governs
the drown state of these domains -- that is where it belongs, and the
lower panel marks them. The counts that used to sit in a box on this
panel are in the caption: they are statistics, not picture.
```

```text
Negatives are tested where numpy actually lands them, at the landward
end of the array, so they are marked apart -- the number is real, but
it does not describe the measured position.
```

```text
The landward-index marker carries no legend entry -- the lower panel
names it, and the entry was the longest item in the row. The band
keeps its swatch.
```

```text
PER VINTAGE (2026-08-26). This called P.load_interiors() with no
argument -- one island under both panels, which the caption used to state
outright. load_interiors() now requires a year, so this could not survive
the change silently.
```

<details><summary>Function notes (the original docstrings)</summary>

**`read_stage()`**

```text
One stage's setback per domain, from the per-domain diagnostics CSV.

Filtered on `setback_model_m` being finite, which is exactly the in-span
set the model-facing file carries -- so all three stages cover the same 82
domains and the progression compares like with like.
```

**`place_stage()`**

```text
Where each road sits at this stage, and what bulldoze's test says about it.

Non-negative setbacks are handed to the imported `place_road` unchanged, so
stages 0/1 and the stage-2 figure cannot drift apart. Only the negative case
is handled here, because stage 2 has none by construction.
```

</details>

### 3-figures/HAT_oceanfloor_offset_check.py

How big the ocean-side move max(setback, 0) is, domain by domain, drawn where it fires.

From the script's original header:

```text
The OCEAN-SIDE move -- `max(setback, 0)` -- is the one modification that changes
a MEASURED number before the model reads it. This script reports how big that
change is, domain by domain, and draws it zoomed in on the only two stretches of
island where it fires.

  setback_dunestart_m          what the profiles measured        (can be negative)
     |  OCEAN-SIDE MOVE
  setback_dunestart_floored_m  max(setback, 0)                   (never negative)

The whole-island stage figures (HAT_dunestart_stage0_raw.png, _stage1_...) show
WHERE the move fires. At 90 domains across a 900 m cross-shore window a 10 m
move is a line width, so they cannot show HOW FAR. This one crops to the
affected domains plus two either side, and takes the y-axis down to the measured
positions, so the displacement is drawn at a scale where it can be read.

WHAT COUNTS AS "DIFFERENT"
Three numbers, because they answer three different questions:

  delta_m       floored - raw, i.e. |raw|. How far the road was moved in metres.
  delta_rows    int(0/10) - int(raw/10). How far it moved in MODEL rows, which
                is what CASCADE actually sees. Truncation toward zero means
                these are not always delta_m / 10 -- a raw of -5 m is a 5 m move
                and a ZERO-row move, because int(-5/10) == 0 already.
  pct_seaward   share of the domain's own profiles that were themselves seaward
                of interior row 0. The domain value is a median; a domain at
                52% is a coin-flip that landed negative, and a domain at 100% is
                a road genuinely out in front of the dune start.

WHAT THE MOVE IS NOT
It is not a correction of the measurement. The road really was measured seaward
of interior row 0 -- that is what a road sitting in the dune/beach strip looks
like when the dune start is the reference. The move exists because a negative
value cannot be handed to `bulldoze`: `int(-60/10) = -6` and
`xyz_interior_grid[-6:-4, :]` is valid Python indexing from the LANDWARD end, so
the road silently lands in the bay. The counterfactual wrap row is reported here
per domain so the alternative to the move is on the record too.

OUTPUTS  (all under dunestart_offset\modifications\)
  HAT_dunestart_oceanfloor_check.png         the figure, captioned to stand alone
  HAT_dunestart_oceanfloor_check_panels.png  panels (a) and (b) and the key only,
                                             for a document that carries its own
                                             caption -- same axes, same draw code
  HAT_dunestart_oceanfloor_check.csv         the same numbers, machine-readable

REQUIREMENTS
  numpy, matplotlib
```

Notes that were in the code:

```text
The counterfactual: what bulldoze's drown test would say about the
wrapped road. Guarded because a raw setback deeper than the island
would index out of the array entirely.
```

```text
And the destination: row 0 is where the move puts it, so the move is
only defensible if the road survives there. Same test, same rows.
```

```text
Seaward of interior row 0 the interior array holds nothing -- this is the
dune/beach strip the road was measured on. Sand, not the water grey: a
road out here is on the beach, and the two must not read the same.
```

```text
The model-facing band, all 20 m of it, drawn as a band rather than a
line -- at this crop the road's own width is legible, and the move is
only meaningful against it. Unaffected neighbours get the same band,
faded, so the size of the move can be read against a setback the
method did not have to touch.
```

```text
Every profile, at its own alongshore position, ON TOP of the band --
the positive profiles of a floored domain sit inside the band's 20 m,
and under it they were invisible exactly where they matter.
```

```text
Wrapped here rather than with matplotlib's wrap=True, which measures
against the FIGURE width and would run the caption under the colorbar.
```

```text
PER VINTAGE (2026-08-26). measure_move() already took a year and was
handed ONE interiors dict for both -- so the 1984 rows were measured
against the 2004-start island. load_interiors() now requires the year.
```

```text
The background strip is context for a per-domain measurement, not the
measurement -- drawn from the FIRST vintage's interiors and labelled as
such rather than pretending to be both.
```

<details><summary>Function notes (the original docstrings)</summary>

**`read_domains()`**

```text
Every in-span domain for a year, with both sides of the ocean-side move.

In-span is `setback_model_m` finite -- the same filter the model-facing file
applies -- so a domain missing here is missing from the run, not from the
move.
```

**`read_profiles()`**

```text
Per-profile setbacks, keyed by domain, in profile order.

The domain number is a median of these. Without them the figure would report
a move on a single value and say nothing about whether the domain as a whole
was seaward of row 0 or straddling it.
```

**`context_rows()`**

```text
Every in-span domain's model-facing setback and profiles, for the year.

The zoom panels carry two unaffected neighbours either side, and a neighbour
with nothing drawn on it is not context. These are what the SAME method puts
on a domain it did not have to floor -- 30-120 m landward here -- which is
the only thing that makes "+60 m" a size rather than a number.
```

**`near_misses()`**

```text
Domains the move did NOT touch that still hold profiles seaward of row 0.

The domain setback is a median, so the move fires on a vote, not on a
physical boundary. A domain at 48% negative is one profile away from being
floored and a domain at 52% is one profile past it -- both exist here, three
domains apart. That is the sensitivity of this modification, and it belongs
beside the six moves rather than in a footnote nobody computes.
```

**`cluster()`**

```text
Affected domains grouped into the stretches of island they sit on.

Not hardcoded to the two stretches that fire today. Re-pick the dune windows
and a third could appear; it should get its own panel rather than be
swallowed by a fixed x range.
```

**`draw_caption()`**

```text
A per-domain table and a caption, in the register of a figure caption.

Four columns, not nine. The counterfactual wrap row and both drowning tests
were columns and are now one sentence -- they are a single result (the
constraint is what keeps these roadways on the island), and a column per
quantity made the reader assemble that result themselves. Every quantity
dropped from the table is still in HAT_dunestart_oceanfloor_check.csv.
```

**`build_figure()`**

```text
The check figure; `panels_only` drops everything but (a), (b) and the key.

Both versions come out of this one function rather than a second script.
The panels ARE the result and they get reused -- in a document that carries
its own caption, in a slide -- so the standalone title, the methods
paragraph and the table are the parts that go, not parts that get re-drawn
somewhere else and drift. The legend stays in both: it is not commentary,
it is what says which band is measured and which is applied.
```

</details>

### 3-figures/HAT_road_domain_views.py

One domain at a time: what CASCADE does to its grid with the road on it.

From the script's original header:

```text
One domain at a time: what CASCADE actually DOES to this grid.

Merged from three scripts that drew the same geometry three ways
(HAT_plot_road_initialization.py, HAT_plot_barrier3d_grid.py,
HAT_browse_road_domains.py). They shared a config block, a geometry() and a
drown test, and each copy was a place for those to drift. The drawing code below
is ported VERBATIM; only the config, the CLI and the setback selection are new.

MODES
  --mode map        land/water plan view. The road as ONE row index applied to
                    all 50 profiles, against an island whose landward edge is
                    not straight, plus the two rows bulldoze tests.
  --mode section    the assembled Barrier3D cross-section -- shoreface, beach
                    wedge, dune rows, interior, bay -- with the REAL bulldoze()
                    run on a copy of the real interior.
  --mode both       both figures per domain (default).
  --browse          interactive walk instead of writing files (n/p/w/q).
                    Needs a GUI backend; everything else runs headless.

WHICH SETBACK IS DRAWN -- the figures are only honest if this matches the run
  --method legacy      archive/superseded_20260911/<year>/RoadSetback_<year>.csv
                       what hatteras_site_config.py spends today
  --method dunestart   dunestart_offset/measured/<year>/RoadSetback_<year>_dunestart.csv
                       measured landward of interior row 0, the reference
                       roadway_manager.py:99 actually uses

Both files are 2 rows x 82 cols, GIS IDs then metres, so this is a drop-in.
Change the default the SAME DAY you change hatteras_site_config.py:78,91.

USAGE
    python HAT_road_domain_views.py --domains 52 --year 1984
    python HAT_road_domain_views.py --domains drowning --mode map
    python HAT_road_domain_views.py --browse --year 2004 --start 52
```

Notes that were in the code:

```text
Backend is chosen in main(), not at import: --browse needs an interactive one
and forcing Agg here would silently kill the window.
```

```text
Topography version resolved from the extractor, not hardcoded -- it was
"2009_v3" and kept drawing v3 interiors under v4 setbacks after the re-pick,
with no error. See hat_topo_version.py.
parents[4] IS scripts/ -- hat_topo_version.py moved there 2026-08-20.
```

```text
BOUND FROM --year, NOT AT IMPORT (2026-08-26).

This was `TOPO_DIR, DUNE_DIR, TOPO_RUN_NAME = topo_dirs()` at module level:
one topography, resolved before the CLI was parsed, so `--year 1984` drew a
1984 setback on 2004-start interiors and `--year 2004` drew a 2004 setback on
the same ones. One of those two was always wrong, and after the tree went
period-first it was the 1984 one -- silently, since a plausible island is
still an island.

Deferred to bind_topography(), called from main() once the year is known.
`SETBACK` is already rebound from --method the same way, so this follows the
pattern the file already uses rather than inventing one.
```

```text
One file, both periods. NOT because there is one surface -- there are two,
and they differ by a median +0.222 m in the road corridor -- but because that
difference is the uncorrected 1996-vs-2009 survey offset rather than a
roadbed. See hatteras_site_config.py, HATTERAS_ROAD_ELEVATION_FILE.
```

```text
Matches hatteras_site_config.py:78,91. Change both together, or these
figures stop describing the model you run.
```

```text
Exact name, not a glob. The old "…_topography_*.npy" needed a year tag
that no longer exists, and would have matched twice if one ever strayed
in from the other period - picking arbitrarily.
```

```text
Direct labels, placed INSIDE the axes. The sea-side colour WARNed on
contrast against a light surface in the validator, so it carries a label
rather than relying on the legend -- and the label sits on an opaque
white plate so the pale field underneath cannot erode it further.
(Placing these outside the axes clipped them against the right panel.)
```

```text
---- right: how many profiles END at each row ----
y is the SAME cross-shore row as the map (sharey), so a bar here lines up
with the row it describes. Counting profiles rather than tracing the edge
again keeps this panel from restating the map's black line.
```

```text
The road label sits BELOW the road box and to the left of the test-row
labels; placing it above collided with the rotated sea-side label.
```

```text
Count over the ACTUAL rewritten block (2 road rows x 50 profiles), not
the median profile -- summing metres down a median is not a quantity.
```

```text
Declared up front: --method's default READS SETBACK_METHOD, and Python
rejects a global statement that comes after a use in the same scope.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_find_project_root()`**

```text
Walk up until a directory holds data\hatteras_init.

NOT parents[N]. These files have moved twice, and the old parents[4]
silently resolved to scripts\ -- every data path below it was wrong, and so
was the sys.path entry the `import cascade.roadway_manager` depends on.
```

**`geometry()`**

```text
Reproduce bulldoze()'s indexing exactly.

Superset of the two originals: the section view used only road_start /
road_end / border / sea / bay / drowns, the map view also needs land, edge
and water_frac. ONE implementation, so the drown verdict cannot differ
between two figures of the same domain.
```

**`ensure_interactive_backend()`**

```text
Get a backend that opens a REAL window and blocks in plt.show().

PyCharm is the reason this exists. With "Show plots in tool window" on
(Settings > Tools > Python Scientific), PyCharm swaps the backend for
'module://backend_interagg', which draws into the SciView panel instead of a
window. plt.show() then returns IMMEDIATELY and no key_press_event is ever
delivered -- so the browser renders the first domain, never receives a
keystroke, and exits. That looks like "it only plots GIS 9".

Jupyter's 'module://matplotlib_inline.backend_inline' behaves the same way.

Testing for the literal string 'agg' does not catch either of them, so match
on the module:// prefix instead and switch to a windowing backend.
```

</details>

### 3-figures/HAT_road_geojson_map.py

The 1978 and 2008 NC-12 lines on the 2009 DEM in map coordinates: a reference sheet for the whole island.

From the script's original header:

```text
The 1978 and 2008 NC-12 geojsons (the 1984 and 2004 periods' roads) drawn on the 2009 DEM, in MAP coordinates, for
the whole modelled island. A reference sheet: find a GIS domain on the real
island, and see where each road line actually runs across it.

  data/hatteras_init/4-mgmt-forcing/road_offset/raster/
      HAT_road_geojson_on_2009_dem.png

WHY THIS IS A MAP AND NOT A DOMAIN-FRAME FIGURE
Every other road figure in this tree lives in the Barrier3D domain frame, which
is reached through orient -> alongshore flip -> shear -> water trim. That chain
is exactly what makes those figures comparable to the model, and exactly what
makes them useless for checking the model's inputs against the world.

This figure applies NO transform. The DEM is read in its native UTM 18N grid and
the geojson is reprojected onto it with `to_crs`, the same single step
HAT_check_geojson_vs_mask.py uses for the same reason. If the road line and the
road visible in the LiDAR agree here, the source data is right; whether the
DOMAIN arrays are right is a separate question that the domain-frame figures
answer.

WHY THE STRIPS RUN LEFT TO RIGHT
Hatteras runs very nearly north-south: GIS 1-90 spans 8.0 km east-west and 45.2
km north-south, 5.7:1. Cut into six 15-domain segments at TRUE NORTH and TRUE
SCALE, every segment is portrait (h/w 1.7-2.7) -- the island is simply taller
than it is wide at any segmentation. Stacking portrait panels vertically would
give a figure four feet tall, so the six run SOUTH (left) to NORTH (right)
instead. North is up in every panel and the scale is true in both axes; only the
reading order is unusual, and the panel titles carry it.

WHAT THE COLOURS MEAN, AND THE ONE DISTINCTION ONLY THIS FIGURE CAN DRAW
The source DEM stores NoData as exactly -10.0 m NAVD88, and everywhere else in
this project that value has already been folded into the water sentinel --
Barrier3D has no representation for "unknown", so by the time the topography is
saved, a LiDAR hole and a genuinely wet cell are the same number. This figure
reads the RAW tif, so it is the one place the two can be told apart, and it
draws them differently: never surveyed is a distinct grey from water.

That matters for reading GIS 78-80, where the roadway relocation fires. Those
domains drown on coverage gaps, not on measured water, and here you can see it.

REQUIREMENTS
  numpy, matplotlib, rasterio, geopandas
```

Notes that were in the code:

```text
MOVED TWICE. 2026-08-25 it went under superseded/ with the pre-90-domain
legacy and this path followed it there. 2026-08-26 it moved back OUT -- it is
a live input read by four scripts -- and this path did NOT follow, so the
script printed "no clip tifs found" and produced nothing. Repointed
2026-08-26; the same stale path was in the offset producer and in
2-audit/HAT_check_geojson_vs_mask.py, where it silently zeroed the check.
```

```text
Never-surveyed is the PALEST thing on the sheet, deliberately. Most of each
2000 m clip box is off-island and therefore NoData, so at any real saturation
it becomes the largest block of colour in the figure and the island reads as
mostly unsurveyed. It is absence of information and should recede; the holes
that matter sit INSIDE the island, where pale against green still reads.
```

```text
The two lines coincide except in the relocation blocks, so drawing them at the
same weight hides 1984 completely. 1984 goes underneath and wide, 2004 dashed
on top: coincident stretches read as a dashed orange line on a blue casing,
and the places they genuinely diverge are the places you see two lines.
```

```text
Extent from the TRIMMED size, not the file's, so the image lands on the
ground it was read from -- the trim discards up to s-1 rows/cols.
```

```text
The tifs carry a COMPOUND CRS (UTM 18N + NAVD88 height). If pyproj
declines to transform onto it, fall back to its horizontal part --
the vertical component is irrelevant to a plan-view reprojection.
```

```text
Three states, drawn as three layers rather than one ramp: water and
"never surveyed" are states, not small elevations, and the whole point
of reading the raw tif is that they can be separated here.
```

```text
No section names on this figure. The strips are only ~3 km wide, so a
boxed label sits on the island rather than beside it and hides the DEM
and the road lines underneath. The GIS numbers already locate a domain,
and the domain-frame figures carry the section bands.
```

```text
True scale in both axes means the panel widths are set by the data, not
chosen. Width ratios come from each strip's own bbox.
```

<details><summary>Function notes (the original docstrings)</summary>

**`constant_from_extractor()`**

```text
Read one numeric constant out of the extractor without importing it.

Importing a 2000-line module for two floats risks its import-time side
effects; copying the numbers risks them drifting apart. Parsing the
assignment is neither.
```

**`read_domain()`**

```text
One domain at full 1 m, block-reduced to the display resolution.

Reduced here rather than by rasterio's decimated read, because the two
quantities need OPPOSITE reductions and rasterio can only apply one per
read (and refuses `Resampling.min` on reads at all):

    elevation   mean of the SURVEYED cells in the block
    no-data     any() -- a block containing any hole is a hole

Averaging a hole together with real ground would invent an elevation and
quietly shrink the coverage gaps, which are the thing this figure exists to
show honestly. NoData is identified on the RAW NAVD88 values before the
datum shift, the same order the extractor uses.
```

</details>

### 3-figures/HAT_road_island_planview.py

Where NC-12 sits on the island CASCADE runs, one plan view per vintage: the model's road, not the measured one.

From the script's original header:

```text
Where NC-12 sits on the island CASCADE actually runs, one figure per vintage.

THIS DRAWS THE MODEL'S ROAD, NOT THE MEASURED ROAD
An earlier version of this script drew the per-profile measurement from
`RoadOffset_<year>_profiles.csv` -- a line that wanders with the real road's
obliquity. That is the survey, and it is not what Barrier3D runs. From
`roadway_manager.py:99`:

    road_start = int(road_setback / dy)          # dy = 10 m, TRUNCATED
    road_width = int(road_width / dx)            # 20 m / 10 m = 2 cells
    road_end   = road_start + road_width
    new_road_domain = np.zeros((road_width, ncols)) + road_ele

So the road the model spends is a **flat rectangle**: two cells cross-shore,
the domain's full 50 profiles alongshore, at one constant row, holding one
constant elevation. No obliquity, no scatter. This figure draws that, because
the question it exists to answer is "where does CASCADE put the road", not
"where did the survey find it".

The consequence is visible and intended: the road steps between domains rather
than curving, and the step is the discretisation the model imposes.

SOURCES ARE THE MODEL-FACING FILES, FOR THE SAME REASON
    setback    dunestart_offset/measured/<year>/RoadSetback_<year>_dunestart.csv
               the 2-row file hatteras_site_config.py resolves as
               PERIOD["road_setback_file"]. Already floored and already
               relocated seaward where a roadway would drown at t=0.
    elevation  road_elevation/RoadElevation.csv
               ONE file for both vintages -- HATTERAS_ROAD_ELEVATION_FILE.
               Not the per-year road_elev_mhw_median in the offset CSVs.

That second one matters for reading the pair: the road colour is identical
between the 1984 and 2004 figures by construction. Any difference you see
between them is the SETBACK moving, or the island underneath changing. It is
never the roadbed, because the model is not given a per-year roadbed.

THE CANVAS
Copied from `_build_island_canvas()` in `HAT_dune_topo_extractor.py`, by way of
`nodata_audit/HAT_plot_island_nodata.py`, so this overlays
`HAT_dune_topo_island_planview_<ver>_<year>_padded.png` cell for cell.
Importing the extractor instead would drag in its interactive picker.

    offsets  2-brie-offset/<year>/Island_Dune_Offsets_*.csv, metres,
             seaward positive, row 0 = domain 1. A 120-row file is stripped of
             its 15 buffer domains per end.
    origin   round(offset_m / 10) - the canvas row interior row 0 lands on
    padding  every domain padded landward to ISLAND_PAD_ROWS = 200 cells
    dune     canvas row origin - 1, one row
    columns  domains concatenated ascending, 50 profiles each, no per-domain
             flip - the arrays already run south to north
    y axis   origin="lower", so increasing row is LANDWARD. The dune sits at the
             bottom edge of the island band and the bay above it.

Each vintage is drawn on its own island - 1984 on `1984-start`, 2004 on
`2004-start`, through `hat_topo_version.YEAR_PRODUCT`. All 90 domains differ
between the products and 65 differ in interior shape.

INPUT   <product>/dune-topo/<version>/topography/domain_<N>_topography.npy  dam
        <product>/dune-topo/<version>/dunes/domain_<N>_dune.npy             dam
        2-brie-offset/<year>/Island_Dune_Offsets_<year>_*.csv      m
        dunestart_offset/measured/<year>/RoadSetback_<year>_dunestart.csv            m
        road_elevation/RoadElevation.csv                                    m MHW

OUTPUT  dunestart_offset/HAT_road_island_planview_<year>.png (and .pdf)
        dunestart_offset/CAPTIONS.md   the caption, keyed by file name
```

Notes that were in the code:

```text
Island elevation. The extractor's poster ramp, copied so the two plan views
are the same picture: `terrain` with 0 m pinned to colormap position 0.35, so
land gets 0.35-1.0 of the ramp and a 2 m dune is still readable against a
4 m maximum. Masked cells take the ocean colour rather than the axes white.
```

```text
Which island ramp. `terrain` is the extractor's, kept as the default so the
pair of plan views stays one picture. `oleron` is Crameri's perceptually
uniform topography map -- equal elevation steps look equal, and it survives
greyscale and colour-vision deficiency, which `terrain` does not.

HAT_ISLAND_CMAP=oleron python 3-figures/HAT_road_island_planview.py

THE SEA-LEVEL POSITION IS NOT THE SAME FOR THE TWO, and this is the whole
trap. `terrain` has no built-in shoreline, so the extractor pins 0 m at ramp
position 0.35 by choice, to spend more of the ramp on land. `oleron` has its
land/sea break BUILT IN at the exact middle of the ramp, so 0 m must be
pinned at 0.50 -- pin it anywhere else and the colormap's own blue-to-green
boundary lands at an elevation that is not sea level, which is worse than the
non-uniform map it replaced.
```

```text
Road elevation. `terrain` already spends blue, green, yellow, brown and white,
so the road ramp has to live in the one region it does not: magenta to deep
purple. The light end of RdPu is nearly white and would vanish against the
dune, so the ramp is truncated at 0.30 -- every road cell stays visibly
magenta, and the sequential light->dark = low->high reading survives.
```

```text
Stroke width for the road, in points. The TRUE band is 2 cells, which at this
figure's scale is ~1.0 pt -- so 1.6 is close to honest and still wide enough
to hold a fill colour beside its own outline. Raising this past ~2 starts
claiming cross-shore extent the model does not have; the printed
"road drawn at N x true width" line reports the factor every run.
```

```text
Embed fonts as TrueType in the PDF. Matplotlib's default is Type 3,
which a number of journals reject outright at submission and which no
vector editor can re-flow. Costs nothing to set. Not in the house
module; everything else here comes from it.
```

```text
The CURRENT build (2026-09-18). This took the first sorted match, which
for 1984 and 2004 is the superseded flat build, not the current one.
```

```text
A domain with no line is ambiguous on the figure -- "no road here" and
"the line is hidden under the dune" look identical. Carry the list so the
drawing can say which it is rather than leaving the reader to infer.
```

```text
Drawn at the width it will be printed: 190 mm, the house double column,
so the 9 pt type on it is 9 pt on the page. It was 20 in wide, where the
same type reduced to about 3 pt. The height is what the two colour scales
and the four labelled axes need; the vertical exaggeration that follows
from it is computed below and stated in the caption.
```

```text
Narrower than the extractor's 0.88 to leave a gutter for the cross-shore
distance axis; the colorbars move right by the same amount.
```

```text
The road is two cells on an 835-row canvas -- about one pixel, which is
why the extractor draws its own road with pcolormesh rather than a line.
Here the road also has to carry a COLOUR SCALE, and a one-pixel band
cannot: the hue would be invisible. So each domain is stroked as one flat
segment at the band's centre, thicker than 2 cells. Position and
alongshore extent are exact; only cross-shore thickness is exaggerated.
```

```text
Mark the domains that carry NO road, as a hatched strip along the bottom
of the frame. Without it, a domain with no line is indistinguishable from
one whose line is hidden -- and the reader has no way to tell which.
```

```text
---- real distance, alongside the model's own indices ----------------
The primary axes carry domain number and canvas cell, which is what you
need to trace a value back to a file. Neither is a length, so a reader
cannot judge scale from them. These add the metric axes: 1 cell = 10 m.
```

```text
Vertical exaggeration, computed from the axes actually drawn rather than
assumed. A plan view whose two axes are at different scales must say so,
or the island looks narrower than it is.
```

```text
Reported in the footer rather than inside the axes: the lower-left corner
is where the no-road hatch lands, and any in-axes corner is a collision
waiting on a different island shape.
```

```text
The endpoints used to be named here as "Cape Hatteras to Rodanthe", which
misplaces the north end by about 5 km: Rodanthe is domain 80, and domain
90 is Pea Island. The one house label carries no endpoints; they are in
the caption instead.
```

```text
NO in-figure title, and no footnote either. A journal sets the caption in
the text, so a title baked into the image duplicates it, cannot be
copyedited, and has to be cropped out by hand. The provenance line that
used to sit along the bottom edge is the same problem one size down, so
it has moved into the caption below with everything else.
```

```text
The caption lands in CAPTIONS.md beside the PNG, keyed by file name,
which is where every other figure in this project keeps its caption. It
was a sidecar <stem>_caption.txt until 2026-09-10.
```

```text
300 dpi is the usual raster floor for a figure submitted to a journal;
the PDF beside it is vector for everything except the rasterized mesh,
so text and the road stay sharp at any zoom.
```

<details><summary>Function notes (the original docstrings)</summary>

**`draw()`**

```text
The extractor's poster styling, with the road carrying a second scale.

Every frame value here is copied from the plan-view block at
HAT_dune_topo_extractor.py:2425 -- ocean background, axes rectangle, figure
aspect rule, colorbar geometry, tick cadence, spine colours, title format.
The one departure is the second colorbar, which the road elevation needs and
the extractor's version has no use for.
```

</details>

### 4-compare/HAT_method_comparison_figures.py

The two setback methods against each other, and both against where NC-12 actually is.

From the script's original header:

```text
The two setback methods against each other, and both against where NC-12
actually is.

  HAT_method_comparison_on_domains.png   both methods, each period on its own
                                         interiors
  HAT_method_vs_actual_road.png          both methods against the rasterized road

Written to road_offset/ itself rather than into either method's folder, because
neither figure belongs to one method.

ENCODING CHANGES HERE -- READ THIS FIRST
In the per-method figures, hue = YEAR. In these two, hue = METHOD and the year
is the panel:

    grey    the superseded method: the minimum road elevation minus the
            minimum dune elevation, taken independently per domain, against
            the same-year digitised dune line
    purple  the current method: per profile, measured landward from the dune
            start, then the domain median

That is the house BASE/ACCENT pair -- the input as it was against the change
under test -- and it leaves the vintage red/blue free. Because hue is NOT the
vintage in these two figures, the red pole is available, and it carries the
one thing that is a failure rather than a category: a roadway that drowns at
initialisation.

WHAT "ACTUAL ROAD" MEANS, AND WHAT IT DOES NOT
The reference band is the RASTERIZED NC-12 mask, per alongshore profile, in the
same frame both methods are drawn in: `road_seaward_cell - interior_row0_cell`,
read from dunestart_offset/measured/<year>/RoadOffset_<year>_profiles.csv. It is the road
as the model grid sees it, with all 50 profiles kept instead of collapsed.

It is NOT an independent check on the dune-start method. That method's setback
IS the median of this quantity, so the two agree by construction, differing only
by `int()` truncation and the negative floor. Read the band for two things it
does show:

  * how far the OLD method sits from the road actually burnt on the grid, which
    IS an independent comparison, because that method never saw this grid;
  * how much the road wanders WITHIN a domain -- the p10-p90 spread that any
    single scalar setback has to throw away, whichever method produced it.

REQUIREMENTS
  numpy, pandas, matplotlib
```

Notes that were in the code:

```text
The placement script owns load_years / load_interiors / place_road and the
transcribed drown test. Importing it keeps ONE implementation: a second copy
would drift, and a drifted method comparison looks like a result.
```

```text
Cross-method output goes in its own folder, NOT at road_offset/ level and NOT
inside either method's folder. The rule this satisfies is unchanged -- a
legacy-vs-dune-start result belongs to neither method -- but the top level is
for the forcing, its source and its inputs, and four loose comparison files
sitting beside them read as though they were part of the product.
```

```text
hue = METHOD here, not year. See the header. BASE is the input as it stood,
ACCENT the change under test. The rasterized road is the OBSERVATION both
methods are measured against, which is what C["REF"] means -- it was the
road ink for one draft, and a 42%-alpha near-black band is the same mid grey
as BASE, so the reference and the superseded method read as one thing. A
drowning roadway is the one failure state, and it takes the red pole, free
in this figure precisely because hue is not the vintage here.
```

```text
--- (C) error against the rasterized road ------------------------------
Both the per-DOMAIN error (line, against the domain's median road) and the
per-PROFILE spread of that error (band, against p10-p90 of the profiles).
The band is ported from the retired
3-figures/island_wide/HAT_plot_road_placement_accuracy.py: a method can sit
on the median road and still miss most individual profiles, and only the
band shows that.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_actual()`**

```text
The rasterized road per domain, kept as a spread rather than a scalar.

seaward_p10/p50/p90 are percentiles ACROSS the domain's profiles of the
road's seaward edge relative to interior row 0. width is the median
measured road width, so the band drawn is a real footprint.
```

**`load_placements()`**

```text
{method: {year: placed}} using the placement script's own logic.

`per` is the placement script's per-VINTAGE bundle, not one interiors dict.
It used to be the latter, which placed the 1984 road on 2004-start
interiors -- see the note at the top of HAT_road_placement_on_domains.py.
```

</details>

### 4-compare/HAT_road_method_diagnostic.py

What changes when the road forcing moves from the legacy method to the dune start, and why, per domain.

Notes that were in the code:

```text
HAT_road_method_diagnostic.py

What changes when the road forcing moves from the legacy method to measuring
against the extracted dune start -- and, crucially, WHY each domain changes.

THE DECOMPOSITION
The legacy-to-new difference confounds two things. A control run with
STRAIGHTEN = False -- same code, same DEM dune crest, same road mask, same
median aggregation, only the frame changed -- separates them exactly:

new_straightened - legacy  =  (new_straightened - new_raw)   FRAME
+ (new_raw          - legacy)    REFERENCE

FRAME      the obliquity correction. North-up clip boxes make NC-12 cross
each 500 m domain diagonally; straightening shears that out.
REFERENCE  everything else, and it is two effects that these files cannot
separate: the legacy method measures a digitized same-year
DUNE-LINE GEOJSON, this one measures the DEM DUNE CREST inside the
picked window, and those differ both in what feature they are and
in what year the feature dates from. Labelled honestly as one
component rather than split on an assumption.

Residual caveat: the two passes use two different pick files (a window picked
in one frame points at different cells in the other), so a little of FRAME is
really window choice. Second-order, not zero.

WHAT IT SHOWS
A  setback per domain, three methods, per year
B  the two components
C  |total change| as a FRACTION OF ISLAND WIDTH -- a 50 m shift is severe on
a 150 m island and minor on a 600 m one, so metres alone hide which
domains matter
D  road elevation, legacy vs new: a near-null result, kept because it
documents that the elevation was never the problem

USAGE
python HAT_road_method_diagnostic.py
```

```text
Derived, not hardcoded: a machine-specific absolute path here made the script
unrunnable for anyone else and silent about why. Two sibling scripts carried
the same defect in a form that resolved to the filesystem root.
```

```text
Topography version resolved from the extractor, not hardcoded -- it was
"2009_v3" and silently survived the re-pick into 2009_v4. See hat_topo_version.py.
parents[4] IS scripts/ -- hat_topo_version.py moved there 2026-08-20.
```

```text
PER VINTAGE, not once (2026-08-26). This was a module-level topo_dirs() with
no argument -- DEFAULT_PRODUCT -- and `island_width_m()` fed the SAME widths
into both years of the table, so `total_change_frac_width` normalised the
1984 change by the 2004-start island. All 90 domains differ between the two
products and 65 differ in interior SHAPE, so that denominator was wrong for
every 1984 row. Same failure as the v3/v4 one above, in product form.
```

```text
array_name() is the single definition of these filenames - the same one
the extractor writes with. Nothing here spells a name.
```

```text
Road elevation is NOT per-year, and it stays that way -- but the reason is no
longer "there is one 2009 DEM". There are two elevation products now, and in
the road corridor they disagree by a median +0.222 m (2009-2014-1996 minus
2009-2014, 54 of 82 domains beyond 0.05 m). That difference is the
uncorrected island-wide 1996-vs-2009 survey offset, not a roadbed, so it is
kept OUT of the forcing: one set, sampled on the 2009-2014 baseline, read for
both years in YEARS. The old per-year RoadElevation_<year>.csv pair does not
exist -- writing two files implied a measured change in roadbed height
between 1984 and 2004 that nothing supports.
See data/.../road_elevation/RoadElevation_audit.md and the note beside
HATTERAS_ROAD_ELEVATION_FILE in hatteras_site_config.py.
```

```text
4-compare output lands in method_comparison/, NOT inside either method's
folder -- this is a legacy-vs-dune-start comparison, so it belongs to neither.
It wrote into dunestart_offset/ until 2026-08-28, which put a cross-method
result under one of the methods it compares and contradicted the rule stated
in this folder's README. Its sibling, HAT_method_comparison_figures.py, was
always correct about that rule; only this one drifted.
```

```text
Hue is the VINTAGE here and line style is the method, so the house pair
applies: the earlier survey red, the later blue. Until 2026-09-10 this file
had them the other way round (1984 blue, 2004 orange), which put it at odds
with every other two-vintage figure in the project.
```

```text
Method is encoded by LINE STYLE, year by hue. That keeps the validated
two-colour palette (normal-vision dE 33.6, worst CVD dE 26.5) and means no
series is distinguished by colour alone.
The prose that used to be burned onto the canvas: a title sentence, a
three-line method paragraph and two statistics boxes. It goes to CAPTIONS.md
beside the PNG, with the numbers filled from this run.
```

```text
Honour the offset script's own span decision rather than
duplicating ROAD_SPAN here: domains it measured but excluded from
the model-facing files are not part of the forcing being compared.
```

```text
Right-aligned: the tall bars are at GIS 11-15, so a left-hand label sits
on top of them.
```

<details><summary>Function notes (the original docstrings)</summary>

**`island_width_m()`**

```text
Median land width per domain, from the interiors CASCADE will read.

`year` is required: the two periods start from different extractions, so
"the island width" is not one number per domain.
```

</details>
