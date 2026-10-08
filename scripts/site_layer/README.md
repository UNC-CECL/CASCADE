# `site_layer/` — what and where Hatteras is

The seam between the site-agnostic `cascade_pipeline/` package and the task
folders, which all consume Hatteras' own paths, place names and house style.
Six modules, each answering one question, each answering it in one place.

| Module | Answers | Imported by |
|---|---|---|
| `hatteras_site_config.py` | *What is this place* — domain geometry, town spans, periods, BE presets, road events, nourishment projects. | 87 |
| `hat_figure_style.py` | *What does a Hatteras figure look like* — one type scale, one palette, one elevation ramp, one way to letter a panel. | 102 |
| `hat_topo_version.py` | *Which* Barrier3D domains — which product (`1984-start` / `2004-start` / `forecast`) and which version inside it. | 77 |
| `hat_observed_rates.py` | *Where is the observed shoreline* — the CoastSat chainage, the per-window rate fits the model is graded against, the transect-to-domain lookup. | 22 |
| `hat_elevation_products.py` | *Which* elevation product and stage — `2009-2014` or `2009-2014-1996`, gapfilled 1 m or resampled 10 m. | 13 |
| `hat_extension_domains.py` | *Which alongshore reach* a run models, and the 500 m bins that number the coast beyond the 90 surveyed domains. | 11 |

## Why the indirection exists

**A hand-built path fails silently.** Four road scripts once hardcoded
`2009-dune-topo/2009_v3`; when the dune windows were re-picked into `v4` they
kept reading `v3` interiors while consuming `v4` setbacks, and 18 domains —
two of the three managed roadways among them — had their drown verdicts
computed on the wrong grid, with nothing raised. The same shape of bug hit
`HAT_road_elevation.py` when a fill source moved under `superseded/`.

Each location is resolved **once**, here, and a name that is not on disk is an
immediate, loud error listing what is. So never rebuild one of these paths by
hand: call `topo_dirs()`, `domain_arrays()`, `array_path()` or `product()` and
let it raise.

## How to import from it

```python
import sys
from pathlib import Path

REPO = next(p for p in Path(__file__).resolve().parents
            if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))

from site_layer.hat_topo_version import topo_dirs, array_name   # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS    # noqa: E402
```

Import the **module**, not the package. `__init__.py` deliberately re-exports
nothing: pulling the six in would make every consumer pay for all of them —
`hatteras_site_config` alone reads several CSVs at import — and would turn the
one real cycle in here (`hatteras_site_config` imports `hat_topo_version`;
`hat_figure_style` lazily imports `hatteras_site_config`) into an import-order
problem.

**Depth no longer matters.** Every module finds the repo root by searching
upward for `pyproject.toml`, not by counting `parents[N]`, which is Rule 5 of
`ORGANIZATION.md`. That is what made moving them off the `scripts/` root on
2026-09-18 a safe operation — the older `parents[1]` spelling would have
silently resolved to `scripts/` one level down.

## The two that also run

`hat_figure_style.py` writes the style sheet (`python scripts/site_layer/hat_figure_style.py`
→ `scripts/figure_making/STYLE.md` and `output/figures/style/`) and `hat_observed_rates.py` prints a
locations diagnostic. Running a file directly puts *its own folder* on
`sys.path` rather than `scripts/`, so both open with a
`if __package__ in (None, ""):` guard that adds `scripts/` back. Importing
them normally never takes that branch.

## Not for another site

`cascade_pipeline/` ships no site content by contract — a different study site
writes its own sibling of this package and never touches it. One leak already
exists (`cascade_pipeline/hindcast.py` imports `hat_topo_version`); folding
site content into the package would make it permanently Hatteras-only.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### __init__.py

site_layer: what and where Hatteras is, each answer resolved in one place.

From the script's original header:

```text
site_layer -- what and where Hatteras is.

The seam between the site-agnostic `cascade_pipeline/` package and the task
folders under `scripts/`, which all consume Hatteras' own paths, place names
and house style. Six modules, each answering one question, each answering it
in one place:

    hatteras_site_config     what is this place -- domain geometry, town
                             spans, periods, BE presets, road events,
                             nourishment projects
    hat_topo_version         which Barrier3D domains -- product and version
    hat_elevation_products   which elevation product and stage
    hat_extension_domains    which alongshore reach, and the 500 m bins
                             beyond the 90 surveyed domains
    hat_observed_rates       where the observed shoreline is, and the rate
                             fits the model is graded against
    hat_figure_style         what a Hatteras figure looks like

WHY THE INDIRECTION EXISTS
    A hand-built path fails silently. Four road scripts once hardcoded
    `2009-dune-topo/2009_v3`; when the dune windows were re-picked into `v4`
    they kept reading `v3` interiors while consuming `v4` setbacks, and 18
    domains -- two of the three managed roadways among them -- had their drown
    verdicts computed on the wrong grid, with nothing raised. Each location is
    resolved once, here, and a name that is not on disk is an immediate, loud
    error listing what is.

    So never rebuild one of these paths by hand. Call `topo_dirs()`,
    `domain_arrays()`, `array_path()` or `product()` and let it raise.

THIS PACKAGE IS DELIBERATELY NOT RE-EXPORTED
    Import the module you need, not the package:

        from site_layer.hat_topo_version import topo_dirs
        from site_layer.hatteras_site_config import HATTERAS_DOMAINS

    Pulling the six into this file would make every consumer pay for all of
    them -- `hatteras_site_config` alone reads several CSVs at import -- and
    would turn the one real cycle in here (site_config imports topo_version;
    figure_style lazily imports site_config) into an import-order problem.

NOT FOR ANOTHER SITE
    `cascade_pipeline/` ships no site content by contract: a different study
    site writes its own sibling of this package and never touches it. One leak
    already exists (`cascade_pipeline/hindcast.py` imports hat_topo_version);
    folding site content into the package would make it permanently Hatteras-
    only.
```

### hat_elevation_products.py

Which elevation product does a script read, and where does it live?

Notes that were in the code:

```text
hat_elevation_products.py

Which elevation product does a script read, and where does it live?

WHY THIS EXISTS
Six scripts used to build these paths by hand, each pasting its own copy of
.../0-elevation/{1-gapfill-1m,2-resampled-10m}/<TAG>/
and that is how HAT_road_elevation.py broke. It sets
FILL_SOURCE = "2008_NOAA_IOCM" and builds GAPFILL_1M_ROOT / FILL_SOURCE.
When 2008 was moved under superseded/ on 2026-08-25 the path stopped
resolving - and NOTHING RAISED. The globs simply returned nothing, the
script carried on, and every domain reported "no fill available".

Same failure as the one scripts/site_layer/hat_topo_version.py was written for: a
layout change that a hand-built path absorbs silently. Same fix. The
product is resolved ONCE, here, and a name that does not exist on disk is
an immediate, loud error listing what does.

THE LAYOUT IT RESOLVES  (product first, stage second - 2026-08-25)

data/hatteras_init/0-elevation/
2009-2014/                  the baseline DEM
1-gapfill-1m/           gapfill_audit.csv + clip_domain_*.tif
2-resampled-10m/        resample_audit.csv + resampled_domain_*.tif
figures/
2009-2014-1996/             the 1984-start DEM
1-gapfill-1m/           mosaic_1984_audit.csv + clip_domain_*.tif
2-resampled-10m/
figures/
superseded/<attempt>/       same shape, not for use
source-selection/           island-wide, belongs to no product
FIGURES.md                  figure design decisions, shared

BEFORE, it was stage first: 1-gapfill-1m/<TAG>/ and 2-resampled-10m/<TAG>/,
with every product's figures pooled in one figures/. The names were the
FILL SOURCE ("2014_NOAA_PostSandy"), which said what was added but not what
the product contained. Product folders are now named for their COMPOSITION,
deliberately not for a hindcast period: the 2009-2014 DEM currently serves
both the 1984 and 2004 periods, so naming it "2004-start" would assert
something that is not true.

USAGE
import sys; sys.path.insert(0, <repo>/scripts)
from site_layer.hat_elevation_products import product, PRODUCTS

p = product("2009-2014")
p.gapfill_1m / "clip_domain_7_filled.tif"
p.resampled_10m, p.figures, p.audit_1m

product() checks the directory exists and raises with the available names
if it does not. Pass check=False only when creating the product for the
first time, which is what the two producer scripts do.
```

```text
The survey codes a product's clip_domain_*_survey.tif may carry, MOST
SPECIFIC FIRST. This is the precedence HAT_dem_resample_clip.downsample_survey
resolves a mixed 2 x 2 block with, so the order is a decision. 2009 is the
base and 0 is "no survey saw it"; neither appears here.
```

```text
Nothing is superseded on disk right now. The 2008 NOAA IOCM attempt was
registered here until 2026-08-26 with the note "kept so its comparison
figures can be regenerated" - but its product folder had already been
deleted, so product("2008_NOAA_IOCM") raised rather than resolving. The
point-cloud path that built it is gone from HAT_dem_gap_fill.py too, so it
is not reproducible from this repo either. The machinery below stays: a
future superseded product registers here with superseded=True and lands
under 0-elevation/superseded/<name>/.
```

<details><summary>Function notes (the original docstrings)</summary>

**`product()`**

```text
Resolve a product by name.

Raises rather than returning a path that is not there. That is the whole
point of this module - see the note on HAT_road_elevation.py above.
Pass check=False when the product is about to be CREATED.
```

**`duneline_check_dir()`**

```text
The dune-line-vs-DEM check for one product: <product>-duneline/.

NOT a product, though it sits beside them: HAT_dem_duneline_coverage.py
and HAT_plot_duneline_offset.py measure where the digitised dune lines
fall on that product's surface, and write here (named 2026-09-18; the two
producers and the row-insert report each spelled the folder themselves).
```

</details>

### hat_env_forcings.py

Where are the sea-level and storm forcings, and the records they come from?

Notes that were in the code:

```text
hat_env_forcings.py

WHERE ARE THE SEA LEVEL AND STORM FORCINGS, AND THE RECORDS THEY COME FROM?

WHY THIS EXISTS
A dozen scripts typed data/hatteras_init/3-env-forcings/... themselves,
and hatteras_site_config spelled every period's storm file out in full.
Same answer as hat_observed_rates.py for 5-scr: resolved ONCE (2026-09-18).

THE LAYOUT, grouped by job (2026-09-18). Before that the records sat beside
the forcings built from them, the WIS export inside storms/, and the storm
validation in two places (storms/storm_check/ -- retired 2026-09-22 -- and a validation/ folder
inside a model-input window).

data/hatteras_init/3-env-forcings/
1-records/            what the forcings are built from, as downloaded
water_level/      Duck gauge 8651370, hourly, m NAVD88 (+ cache)
WIS_raw_data/     WIS station 63228 wave export (git-ignored)
storm_record/     the hurricane history this is checked against
wave_climate_duke/  e_phi_0_OBX_yearly.nc; no script reads it
2-rslr/               sea level: record/ fits/ figures/ (09-15 layout)
3-storms/
hindcast_storms/<window>/   the MODEL INPUTS, one per window,
plus 1984_2024/ (spliced, not a period)
validation/<window>/        both validators' output
figures/
archive/              retired storm series, and Roya's figure

The RSLR record sits in 2-rslr/record/, not 1-records/: rslr/ was laid out as
record -> fits -> figures on 2026-09-15 and is kept whole.

A WINDOW IS <start>_<end> (the end is a boundary; the model spends
start..end-1), the same naming rule as the rest of the init tree.
```

```text
benton_storms/ and testing_storms/, retired 2026-09-14; the storm_check
validators still read testing_storms/base_storms/.
```

```text
THE STORM SERIES THE HINDCAST RUNS ON (2026-09-29, Hannah adopted split12).
v3_split12_trim24: the builder's 24 h grouping, then each grouped event split
where the water stays below the berm >= 12 h (short pieces folded into a
neighbour), then every event cut to the 24 h above the berm around its peak.
The split puts back storms the grouping had chained to a larger one and the
trim then removed -- Fran 1996 (with Edouard), Jose 2017 (with Maria) -- and
changed no score (experiments/storms-overwash-and-dunes/2026-09-29-event-splitting).
Earlier defaults, both still on disk for reproducing older runs:
v3_trim24  2026-09-28 .. 09-29: every event kept, trimmed to 24 h, not split.
v3_72      until 2026-09-28: events over 72 h DROPPED (Isabel 2003, March
2018, Florence, Dennis), a limit that existed only because the
pre-49fd069 Barrier3D crashed on long storms.
Record: 3-storms/PROVENANCE.md; experiments/storms-overwash-and-dunes/.
```

<details><summary>Function notes (the original docstrings)</summary>

**`storm_series_file()`**

```text
The .npy CASCADE reads for a window: <window>_storms_<variant>.npy.
v3_split12_trim24 (the default since 2026-09-29): grouped events split at
>= 12 h below the berm, each trimmed to 24 h around its peak. v3_trim24
(2026-09-28): trimmed, not split. v3_72: events over 72 h dropped.
```

</details>

### hat_extension_domains.py

The alongshore reach a run models, by name, and the domain numbers beyond the 90 surveyed ones.

Notes that were in the code:

```text
hat_extension_domains.py

THE ALONGSHORE REACH A RUN MODELS, BY NAME -- and the 500 m bins that give
the coast beyond the 90 surveyed domains a domain number.

WHY THIS EXISTS (2026-09-16, the Pea Island extension experiment)
The hindcast models GIS 1-90, Cape Point to just north of Rodanthe, and
pads each end with 15 invented buffer domains that extrapolate the local
shoreline slope and then bridge back to close BRIE's periodic ring. The
edgeBE preset pins GIS 1 and 90 to their observed rates with source/sink
values that absorb whatever the buffer gets wrong.

Hannah's question: what if the buffer carried the REAL coast's orientation
instead -- Pea Island north to near Oregon Inlet, and the last kilometre
south to Cape Point -- with the end domains re-solved at the new ends?
The dune lines (every vintage) and the CoastSat record both reach Oregon
Inlet, so the coast is measured; only its domain numbers were missing.

THE NUMBERING
Hannah's whole-island domain polygons
(1-barrier3d-domains/domain-geojson/domains_pea_hatteras_120.geojson,
ID 1-121, EPSG:3725, 2000 x 500 m each) continue the surveyed numbering
north to Oregon Inlet; ID 1-90 are the surveyed polygons vertex for
vertex. An extension domain is numbered exactly as a surveyed one was: a
spatial join onto its polygon (join_lines, join_origins below). Nothing
about GIS 1-90 changes. An earlier numbering by northing bin (a 502.56 m
grid continued from GIS 90, the same evening) was removed once the
polygons existed; it had placed 4 dune-line and 40 CoastSat transects one
domain off, because the drawn polygons are not on a regular grid.

THE GEOMETRIES
"base"    GIS 1-90, the production reach. Every matrix run.
"n115"    GIS 1-115: 25 domains of Pea Island added, stopping ~3 km short
of Oregon Inlet where the rates turn inlet-dominated (Hannah,
2026-09-16: "lets do 115"). GIS 1 stays the southern end. No
southern extension: no polygon lies south of GIS 1, and the
1997 dune line ends 440 m south of it in any case.

HAT_GEOMETRY in the environment selects one; hatteras_site_config builds
HATTERAS_DOMAINS from it. Extension domains carry the shared buffer
topography (hat_topo_version.domain_arrays), the measured dune-line offset
(2-brie-offset/<year>/ext/<geometry>/), zero background erosion unless
HAT_BE_OVERRIDE names them, and no management of any kind.

USAGE
from site_layer.hat_extension_domains import GEOMETRIES, gis_bounds, join_lines
first, last = gis_bounds("n115")          # (1, 115)
join_lines(transects_gdf)                 # -> domain per transect line
```

```text
The surveyed reach. Topography, management inputs and the committed
CoastSat rate table all cover exactly this; anything outside it is an
extension domain.
```

```text
TWO JOIN RULES, because two different joins made the surveyed inputs and
each is reproduced exactly on GIS 1-90 (checked 2026-09-16):
join_lines    the 100 m dune-line transects: the LINE intersects the
polygon (450/450). A transect runs due west across the
island and meets one polygon; one in a sliver between
polygons meets none and is dropped, as ArcGIS dropped it.
join_origins  the CoastSat transects: the ORIGIN POINT (first vertex)
within the polygon (906/906), the rule of
coastsat_domain_mapping.py. (The centroid does not
reproduce it: 737/906.)
```

<details><summary>Function notes (the original docstrings)</summary>

**`domain_polygons()`**

```text
The whole-island polygons as a GeoDataFrame (gis, geometry) in
EPSG:3725, or None if the file is absent.
```

**`_join()`**

```text
gis per row of `gdf` (a GeoDataFrame in any CRS), NaN where no
polygon matches; the first match where several do.
```

</details>

### hat_figure_style.py

One typographic and colour standard for every Hatteras figure.

From the script's original header:

```text
One typographic and colour standard for every Hatteras figure.

WHY THIS EXISTS
    Each plotting script had been choosing its own font sizes, its own greys and
    reds, its own elevation ramp and its own way of labelling panels. Put two of
    them side by side in a document and they read as coming from different
    papers. Worse, the same quantity was drawn in different colours in different
    figures, so "grey" meant "v1" in one and "not applied" in another.

    Until 2026-09-10 there were TWO of these: this module (the row-insert era:
    DejaVu, a grey/dark-red base-accent pair, captions burned onto the canvas)
    and the STYLE block inside 0-elevation/3-figures/HAT_plot_duneline_offset.py
    (2026-09-04, Hannah: "more academic/professionally styled": Arial, the
    ColorBrewer RdBu poles, panel letters, nothing on the canvas that belongs in
    a caption). Hannah's call that day was that the 09-04 rules win wherever
    the two disagreed, and that the one module lives here, importable by every
    script the way hat_topo_version is, with the rules written out in
    scripts/figure_making/STYLE.md and the style drawn in output/figures/style/
    (both written by `write_style_sheet()`, or by running this file; they sat
    together in data/hatteras_init/9-figures/ until 2026-09-18).

THE RULES (the 2026-09-04 house style)
    typeface        Arial first (Helvetica, Liberation Sans, DejaVu Sans behind
                    it), 8-10 pt; text and axes in INK, secondary text in
                    INK_MUTED; thin 0.6 pt axes; hairline grid only when asked
    panels          a bold letter at the left of each panel title, `_title()`;
                    inside the corner when the title is wide, `_letter_inside()`
    maps            a north arrow and a scale bar on any map WITHOUT coordinate
                    ticks; a labelled UTM frame needs neither
    legends         frameless, outside the axes where the layout allows
                    (fig.legend(loc="outside ...") under constrained_layout)
    vintages        the EARLIER line or surface is C_1984 red, the LATER one
                    C_1997 blue, everywhere the two are drawn together; the
                    light fills are the bands between them
    canvas          NO title sentences, statistics lines or footnote
                    paragraphs on the image. That text goes in a CAPTIONS.md
                    beside the figure; `caption()` here does exactly that
    folders         a figure folder shows FIGURES: PNGs at the top, and the
                    PDFs, CAPTIONS.md, tables and provenance under
                    `supporting/` (`save()`, `record_caption()` and
                    `support_dir()` put them there)
    elevation       drawn in classes, not a ramp (`elevation_cmap()`); the
                    terrain colormap of HAT_plot_1984_mosaic is the one
                    deliberate exception, for the 1984-start DEM panels

COLOUR SEMANTICS -- do not reassign these locally
    C["BASE"]     the unmodified input (v1 / v2 as extracted)
    C["ACCENT"]   the modification under test (the insert, the built version)
    C["ROAD"]     NC-12
    C["ADDED"]    ground that was fabricated
    C["WATER"]    cells at or below sea level
    C["REF"]      a reference value: a median, a target, an observation
    C_1984 / C_1997 (and the _FILL pair)  the earlier / later vintage

ELEVATION IS DRAWN IN CLASSES, NOT A RAMP
    A continuous ramp is the wrong tool here. The back-barrier sits a few
    decimetres below MHW and the dune is five metres above it, so a linear ramp
    renders the entire island as one flat tone and hides the only distinction
    that matters -- which cells are land. `elevation_cmap()` returns a discrete
    scale with a hard break at 0 m.

USAGE
    from site_layer.hat_figure_style import apply_style, C, C_1984, C_1997, INK, _title
    apply_style()                       # first, before any figure is made
    ...
    caption(fig, "what the reader needs")   # lands in CAPTIONS.md on savefig

    python hat_figure_style.py          # (re)writes STYLE.md and the style sheet
```

Notes that were in the code:

```text
RUN AS A FILE, NOT IMPORTED. `python scripts/site_layer/hat_figure_style.py` puts this
file's OWN folder on sys.path, not scripts/, so a `site_layer.` import cannot
resolve and the __main__ block below would die on it. Importing the module
the normal way never takes this branch -- __package__ is "site_layer" then.
```

```text
WHERE FIGURES GO (2026-09-18). Every script in scripts/figure_making typed
output/figures/<subject> itself. The subjects are the ones
scripts/figure_making/README.md lists; a name outside them raises rather than
quietly starting a new top-level folder, which is how
output/figures/management_investigation/ came to sit beside management/.
```

```text

THE NUMBERED LAYOUT (2026-09-29, Hannah): top folders follow the paper's
order -- where, what was observed, what the model is fed, how it works, what
it gives. The old subjects (site, forcing, management, shoreline, model,
pipeline, initialization) were retired the same day and now raise, so a
script still using one fails loudly instead of rebuilding the old tree. The
pre-reorganisation tree is output/archive/2026-09-29_figures-pre-reorg/.
```

```text
The ColorBrewer RdBu poles: a warm/cool pair that stays distinct in
greyscale and under red-green colour deficiency. The SAME pair is used
wherever two vintages are drawn together, so red always means the earlier
line (1984), blue the later one (1997, 2004), and a fill the band between.
```

```text
Colour-blind safe: the base/accent pair is grey against a dark red, which
separates on luminance as well as hue, so it survives greyscale printing.
ACCENT was a dark red (#9e2a2b) until 2026-09-10, all but identical to the
vintage red C_1984, so "the modification under test" and "the 1984 line"
read as one colour across a document. It is now the PRGn purple pole:
distinct from the vintage pair, from REF green and from ADDED orange, and
still separated from BASE grey on luminance.
```

```text
An ORDERED variable gets a sequential ramp, not the vintage pair: the
alongshore smoothing width, drawn light (unsmoothed) to dark (widest). It is
anchored on the shoreline blue C_1997, already the CoastSat target's colour,
and stays clear of the amber shoals, the grey village bands and the black
model line, while surviving greyscale as a light-to-dark sequence. Added to
the style 2026-09-21, when the third script wanted the same four blues.
```

```text
COLUMN WIDTHS. A figure is drawn at the width it will be printed, so its
8-9 pt type is 8-9 pt on the page: 90 mm for a single column, 190 mm for a
double. Before 2026-09-10 figures were 11-19 in wide and their text shrank to
4 pt when reduced to a page. `figsize()` is the only way a figure should get
its size.
```

```text
ONE LABEL FOR THE ALONGSHORE AXIS. Three phrasings were in use ("GIS domain
(south -> north)", "CASCADE domain (1 at Cape Point, 90 at Pea Island)",
"domain (1 = south, Cape Hatteras)"); the endpoints belong in the caption.
```

```text
Only the spans this panel actually shows, and the label clamped to the
visible part of its span. The label sits at a DATA x, so on a panel
covering 30 of 90 domains an unrestricted label lands far off-axes and a
`bbox_inches="tight"` save then grows the figure to several times its
width (found 2026-09-10 on the footprint grid panels).
```

```text
Where a structure label may sit, tried in order: along the bottom of its
line, left of it then right; then along the top, under the village names.
```

```text
A line crosses the box BETWEEN its vertices (one vertex per domain,
a label a fraction of a domain wide), so test the segments, not
the vertices: 20 points along each.
```

```text
a 2 pt gap between the text and its line, so the halo that keeps
the text legible over data does not white out the line itself
```

```text
A FIGURE FOLDER SHOWS FIGURES. Everything a figure script writes beside the
PNG -- the vector copy, CAPTIONS.md, tables, provenance -- goes under this
subfolder, so opening the folder shows one image per figure and nothing
else (Hannah, 2026-09-15). `save()` and `record_caption()` do it for the PDF
and the captions; a script writing its own CSV or PROVENANCE.md uses
`support_dir(folder)` for the path.
```

```text
Added 2026-09-11 for the groin-sweep figures, which draw six error surfaces
over a parameter grid (the (M, f) heatmaps, the joint-fit surfaces, the
preset comparison). Each had chosen its own ramp -- magma_r, viridis,
RdYlBu_r -- so the same quantity was a different colour in adjacent figures,
and the saturated ramps collided with the annotations laid on top: magma's
purple end against the ACCENT marker, viridis's green against REF.

The rule is that a scalar error surface is drawn WITHOUT hue. It is the
background against which a best cell, a chosen pair, an iso-product curve or
a constraint is marked, and those marks are what the reader is meant to find;
reserving all colour for them is what makes them findable. Dark is worse,
which matches the convention of every error plot in the project.

Truncated at both ends: pure white reads as missing data (a sweep grid has
real holes, which `pcolormesh` leaves as the axes background) and pure black
hides a marker drawn on top of the worst cell.
```

```text
Figures in this project get read months later out of a folder, detached from
whatever conversation produced them, so each one needs to state its own
method somewhere. Until 2026-09-10 `caption()` wrote that text onto the
canvas; the house rule is that nothing on the image belongs in a caption, so
it now lands in a CAPTIONS.md next to the PNG, keyed by the file name, when
the figure is saved. Callers change nothing: `caption(fig, text)` then
`fig.savefig(path)` as before.
```

```text
Only the PNG gets an entry. `save()` writes a PDF beside it with
the same stem, and recording both put the same paragraph in
CAPTIONS.md twice (found 2026-09-10 during the restyle pass).
```

```text
A replaced entry's match stops at the newline before the next entry,
so the blank line between them was consumed on every re-run and the
file collapsed into one run-on block (found 2026-09-15). Put one
blank line back before every entry.
```

```text
A given out_dir takes both files (a trial); the default splits them: the
sheet is a figure, STYLE.md is documentation for whoever writes one.
```

```text
Three rows since 2026-09-11: the error ramp needs a strip of its own, and
nesting it under the elevation panel collapsed both to zero height under
constrained_layout.
```

<details><summary>Function notes (the original docstrings)</summary>

**`figsize()`**

```text
(w, h) in inches: `width` is "single", "double" or a number of inches;
the height is `height` or width * aspect, capped at a page.
```

**`town_bands()`**

```text
Village spans as light bands behind an alongshore axis, named once.
`spans` is {name: (first_gis, last_gis)}; default the site config's.
Call it AFTER the axis limits are set: spans outside the view are skipped
and a label is clamped to the visible part of its span.

`strip` draws the spans as a band of that fraction of the axes height
against the `where` edge, instead of a full-height wash. Use it when the
panel already shades something else: two full-height greys on one panel
cannot be told apart (the road alongshore figure, 2026-09-10).
```

**`_points_under()`**

```text
How many plotted points fall inside `txt`'s box (a little padded), over
every Line2D and LineCollection drawn in data coordinates.
```

**`_place_label()`**

```text
Put `name` along the line at `pos` in the first slot that covers no
plotted point, else the slot that covers fewest (Hannah, 2026-09-15:
labels were sitting on the data at Rodanthe Pier).
```

**`structures()`**

```text
The Buxton groin (solid hairline) and the two piers (dotted hairlines)
on an alongshore axis, each named once along its own line, reading upward
(Hannah, 2026-09-15). The lines stop short of the top so the village
labels there stay clear of them. A label goes at the bottom of its line
unless data is drawn there, in which case it moves to the other side of
the line or to the top -- so call this AFTER the data is plotted AND
after anything that resizes the axes at draw time (an outside legend, a
colourbar): the test is made in the layout as it stands.
`spans` is a site AnnotationConfig; default the Hatteras one. Lived in
coastsat_lrr_windows.py until 2026-09-15, when a third alongshore figure
wanted it.
```

**`support_dir()`**

```text
`<folder>/supporting/`, created. Where a figure script puts everything
that is not a PNG.
```

**`save()`**

```text
PNG at 300 dpi in the folder and, for `vector`, a PDF with the same
stem under `supporting/`, so a line or bar figure stays sharp in a
manuscript without the folder showing two files per figure. Returns the
paths.
```

**`apply_style()`**

```text
Idempotent. Called at the top of every figure so the style holds no
matter which entry point drew it, including an import from elsewhere.
```

**`error_cmap()`**

```text
The greyscale ramp for a scalar error or cost surface; dark is worse.

`reverse=True` for a surface where HIGH is better (a score, a share
explained), so that dark still means the outcome you do not want.
```

**`_np_linspace()`**

```text
numpy.linspace without importing numpy at module scope: this module is
imported by scripts that have not yet chosen a backend, and it has stayed
free of the numeric stack.
```

**`_letter_inside()`**

```text
The letter inside the top-left corner, for a panel whose title is wide
enough to run under a letter placed beside it.
```

**`panel_title()`**

```text
"(a) text" as one string, for a script that sets the title itself.
`_title()` is preferred: it puts the letter in bold at the left and the
title centred, which is what the house figures do.
```

**`_north_arrow()`**

```text
A north arrow in axes fraction. Only on maps WITHOUT coordinate ticks -
a labelled UTM frame is already north-up by construction.
```

**`_scalebar()`**

```text
A bar, because the coordinate ticks are gone. Drawn in data units, so
it scales with the panel and cannot disagree with it. White halo rather
than a box, so it sits on relief without blanking it. Under 1 km the label
also says how many Barrier3D cells that is (pass show_cells=False to
suppress it on a map that is not a model grid).
```

**`caption()`**

```text
Register `text` as the figure's caption. Written to CAPTIONS.md beside
the image by the figure's next savefig; nothing is drawn. `y` and `size`
are accepted for old callers and ignored.
```

**`_prune_captions()`**

```text
Drop entries whose figure is no longer beside the CAPTIONS.md.

`record_caption` replaces an entry by name but never removed one, so a
renamed or deleted figure left its caption behind forever (68 such
entries had accumulated by 2026-09-21, after the filename rename).
Pruning on every write keeps the file honest without anyone having to
remember.

Only an entry whose named file is MISSING is dropped, so a run that
redraws one figure never touches the captions of the others.
```

**`mark_offaxis()`**

```text
Mark values beyond +/-`half` at the axis edge, and say which they were.

A LINE that leaves the axis is simply cut by `ylim`, with nothing on the
canvas to say so, so a domain at +136 m reads as +100 (found 2026-09-22,
when the metre figures were fixed to +/-100 m). Scatter points already get
this treatment; this is the same for a line, and the caption clause below
names the values so nothing is lost silently.

Args:
    ax: the axes, already at its final ylim.
    x, y: the series as drawn. y beyond +/-half is what gets marked.
    half: the axis half-range.
    color: marker colour; the axes' ink by default.
    size: marker area in points squared.

Returns:
    [(x, y), ...] for the off-axis points, largest |y| first, for
    `offaxis_clause()`.
```

**`offaxis_clause()`**

```text
" Beyond +/-100 m, off the axis: the observed change +136 m at GIS 1."

Args:
    named: [(label, [(x, y), ...]), ...] as returned by `mark_offaxis`.
    half: the axis half-range, for the sentence.
    unit: the y unit.
```

**`compare_header()`**

```text
The 'what is being compared' line(s) above the panels.

Hannah, 2026-09-22: on a figure that puts two DIFFERENT measurements side
by side, which is which -- and over what dates -- has to be on the canvas,
not only in the caption. This is the one thing the style lets above the
panel titles, and it is a NAMING line, not a result: what each side is and
what interval it spans, never the numbers that came out. The summary stays
in `supporting/CAPTIONS.md`.

`target_comparison` carried this idea first (the source/sink line); it
is here so the dune-line comparisons place it identically.

Args:
    fig: the figure.
    lines: one string, or a sequence joined with newlines.
    size: point size; 8.5 sits just under the 10 pt panel titles.
```

**`record_caption()`**

```text
Write or replace the entry for `png_path.name` in the CAPTIONS.md under
`supporting/` beside it. Entries are '**`<file>`.** text' paragraphs;
other content is kept.
```

**`write_style_sheet()`**

```text
A swatch figure and a STYLE.md stating the rules, so the standard can be
seen and read without opening this file. Re-run after changing anything
above; both files say when they were written.
```

</details>

### hat_map_layers.py

Where are the shared map layers, and where does the house style live?

Notes that were in the code:

```text
hat_map_layers.py

WHERE ARE THE SHARED MAP LAYERS, AND WHERE DOES THE HOUSE STYLE LIVE?

WHY THIS EXISTS
Until 2026-09-18 these sat in data/hatteras_init/9-figures/, numbered like
a model-input stage although it was not one, and mixing three kinds of
thing: GIS layers (data), the written style (documentation for whoever
writes a figure script) and a rendered style sheet (a figure). They are
split by kind now:

data/hatteras_init/map_elements/     the map layers, un-numbered
hatteras_outline/                the island, NC 1:80k (subset)
nc_coast_80k/                    NC 1:80k clipped to the study area
natural_earth/                   the locator-map states
archive/                         the retired 1000 m domain polygons
scripts/figure_making/STYLE.md       the written house style
output/figures/style/                the rendered style sheet

Two of the layers were read straight off the D: GIS drive, so the figures
that drew them could not be made without it: the domain boxes
(D:/Hatteras_GIS/domains.geojson, identical in geometry and every
attribute to 5-scr's HAT_domains.json, which is what DOMAIN_BOXES points at)
and the NC 1:80k coastline (D:/Hatteras_GIS/Outlines/nc_80k/, 6 MB for the
whole state; NC_COAST is the 25 km window around the domains, rebuilt by
scripts/figure_making/tools/clip_nc_coast.py).

THE DOMAIN BOXES ARE NOT A MAP LAYER TO COPY. They are the model's frame and
are owned by 5-scr/2-transect-frame/; DOMAIN_BOXES re-exports that path so a
figure script needs one import, not a second copy that could drift.
```

```text
Hatteras Island only: the 208 NC 1:80k land polygons it overlaps, merged
into 52 (checked 2026-09-18: 88.85 km2 either way, zero difference).
```

```text
The coast around it -- sound shores, Ocracoke, Pea Island -- for maps whose
window reaches past the island (the overwash maps pad 7.2 km west).
```

### hat_observed_rates.py

Where is the observed shoreline data: the CoastSat chainage, the rate fits and the transect lookup?

Notes that were in the code:

```text
hat_observed_rates.py

WHERE IS THE OBSERVED SHORELINE DATA - the CoastSat chainage, the per-window
rate fits the model is graded against, and the transect-to-domain lookup that
ties them to Barrier3D domains?

WHY THIS EXISTS
Until 2026-09-12 these lived under `scripts/input_prep/5-scr/CoastSat/`.
They are DATA, and the rate CSV is a MODEL INPUT: section 8 of the hindcast
runner reads it on every run and the calibrated source/sink preset is fitted
against it. Everything else the model ingests comes from
`data/hatteras_init/`, so one input lived somewhere no one would look for
it, and anyone archiving or sharing the data tree shipped a model that
could not run.

Twenty-two files built that path by hand. Moving it once meant editing all
of them, which is the failure `hat_topo_version.py` was written to end for
topography. Same answer here: the location is resolved ONCE, and a window
that is not on disk is a loud error naming the ones that are.

THE LAYOUT, grouped by job (2026-09-18). Each group is numbered in the
order the work runs: observations, then the frame that ties them to
domains, then the rate fits built on both, then the comparisons.

data/hatteras_init/5-scr/
1-observations/             measured or digitized, not fitted by us
coastsat_timeseries/    raw per-transect chainage, by CoastSat site
dsas_1978_2019/         the DSAS rates, a different source
shoreline_inventory/    study-area and reference shorelines
2-transect-frame/
transect_domains/       the lookup, the transect layer, the domain
polygons, and the verification set
3-rates/                    the MODEL TARGETS: tables + one figure per
window (rates_figures.py).
Grouped by source since 2026-09-18.
coastsat/
lrr/<window>/       the OLS rate fits, one folder per WINDOW:
1984_2004 1996_2010 2004_2024 2010_2024,
and 1996_2024 (CONTEXT, not graded)
endpoint/<window>/  net change between +/-6-month means at
the dune-line dates (m and m/yr)
5yr_bins/<window>/  the OLS in successive 5-year bins,
1996_2010 2010_2024 1996_2024
window_convergence_1996_2024/ the same OLS on NESTED families of
1-rate_profiles/ windows, one pinned at each end: which
2-r_bias_rmse/   windows recover the long-term rate?
3-settling_window/
experiments/     See window_convergence_dir().
Seven tolerances scored, headline is
CI overlap
total_change/       the same rate as a DISTANCE: LRR(W) x the
<window>/       years of W, beside the observed change.
1996_2024 1996_2010 2010_2024
(was lrr_projected/ until 2026-09-21 --
nothing in it was ever a projection)
projected/          the 1996-2024 LRR carried onto a window it
<window>/       was NOT fitted on: x 14 yr over 1996_2010
and 2010_2024. The only projections here.
duneline/
endpoint/<window>/  net change between the two dune lines
(m and m/yr; replaced duneline_lrr/)
4-comparisons/
shoreline_vs_duneline/  (one question folder since 2026-09-19)
coastsat_endpoint_vs_duneline_endpoint/<window>/
(+ all_windows_stacked/)
coastsat_total_change_vs_duneline_endpoint/<window>/
(+ all_windows_stacked/, change_between_periods/)
duneline_positions/     where the 1997/2009/2023 dune lines sat:
maps, zooms, dune-road, beach width
(coastsat_windows/ and duneline_windows/ archived 2026-09-19 as
duplicates of 3-rates; trajectory_patterns/ and
two_period_comparison/ outputs deleted as stale)
archive/                    retired windows, old 5-year-bin runs,
the Rodanthe poster figures; not for use

Before 2026-09-18 all of these sat side by side at the top of 5-scr/, under
their old names (scr-dsas-1978-2019, coastsat_timeseries_lrr,
coastsat_lrr_windows, shoreline_change_patterns).

A WINDOW IS <start>_<end>, NOT A START YEAR. A rate fit spans an interval, so
it is named for one - unlike a dune line or a road alignment, which is a
survey at a moment and is named for its year. That split is the naming rule
the data tree follows (Hannah, 2026-09-12).

USAGE
from site_layer.hat_observed_rates import lrr_csv, transect_lookup
path = lrr_csv(1996, 2010)          # raises, listing windows, if absent
```

```text
RUN AS A FILE, NOT IMPORTED. `python scripts/site_layer/hat_observed_rates.py` puts this
file's OWN folder on sys.path, not scripts/, so a `site_layer.` import cannot
resolve and the __main__ block below would die on it. Importing the module
the normal way never takes this branch -- __package__ is "site_layer" then.
```

```text
3-rates is grouped BY SOURCE since 2026-09-18 (Hannah): coastsat/{lrr,
endpoint,5yr_bins} and duneline/endpoint. Until then the four sat flat as
coastsat_lrr/, coastsat_endpoint/, coastsat_5yr_bins/, duneline_endpoint/.
3-rates holds the tables and ONE house-style figure per window beside them
(scripts/input_prep/5-scr/3-rates/rates_figures.py); comparisons are in 4-comparisons.
```

```text
WHAT EACH WINDOW IS. Five window folders sit as peers under most products
and they are NOT equivalent: two chains, and one context window that is not
graded at all. Nothing in a folder name says so, which is how 1984_2004 gets
read as a current result (Hannah, 2026-09-21). This dict is the ONE
definition -- data/hatteras_init/5-scr/WINDOWS.md is written from it and
coastsat_total_change.py imports it rather than keeping its own copy.

All four periods are live in hatteras_site_config.HATTERAS_PERIODS; "legacy"
here means superseded as the MAIN chain, not deleted. See
[[cascade-canonical-periods]]: 1996 -> 2010 -> 2024 is the main chain.
```

```text
Shoreline vs dune line, ONE question folder since 2026-09-19 (Hannah):
endpoint_net_change/<window>/ (coastsat_vs_duneline.py, was 4-comparisons/
coastsat_vs_duneline/), endpoint_net_change/chains/ (net_change_vs_duneline.py, was
4-comparisons/net_change_1996_2024/), total_change/<window>/
(total_change_vs_duneline.py).
```

```text
Named for BOTH sides since 2026-09-22 (Hannah): the folder alone says what
was measured and how, with no lookup. The dune side is the same measured
endpoint in every comparison here -- the 2009 digitized line minus the
1997 one -- so what the old names failed to say is how the SHORELINE was
read, which is the only thing that differs.
was endpoint_net_change/  -> coastsat_endpoint_vs_duneline_endpoint/
was total_change/         -> coastsat_total_change_vs_duneline_endpoint/
```

```text
The stored dune-line observation (2026-09-18): the NET CHANGE between the
two lines that bound a window, per transect and per domain, in m and m/yr,
one folder per window. Written by
scripts/input_prep/5-scr/3-rates/duneline/duneline_endpoint.py. It replaced
duneline_lrr/ (an OLS through every line in the window, 2026-09-16; Hannah:
"we are tracking net change"); its archived copy,
archive/duneline_lrr_retired_20260918/, was deleted 2026-10-01 (recover
from git).
```

```text
The CoastSat counterpart (2026-09-18): net change in the CoastSat shoreline
between +/-6-month window means centred on the SAME dune-line survey dates,
so the two products difference like for like. Written by
scripts/input_prep/5-scr/3-rates/coastsat/endpoint/coastsat_endpoint.py.
```

```text
THE VOCABULARY (Hannah, by interview, 2026-09-21). An LRR turned into a
DISTANCE is named by the window it was FITTED on, not by the arithmetic:

TOTAL SHORELINE CHANGE   the rate is evaluated over the SAME window it was
fitted on. LRR(1996-2010) x 14 yr is total change,
so is LRR(1996-2024) x 28 yr. A shorter span
INSIDE the fit window (the 25.72 yr dune-line
interval) is still total change; the span is named
in the title, not in the folder.
PROJECTED SHORELINE      the rate is carried onto a window it was NOT
CHANGE                   fitted on. LRR(1996-2024) x 14 yr over 1996-2010
or 2010-2024 is the only case in the project.
OBSERVED CHANGE          no rate anywhere: the mean position over the end
calendar year minus the mean over the start one.

Before 2026-09-21 the folder `lrr_projected/` held all three windows built
from their OWN rate -- i.e. no projections at all -- which is what the rename
below fixes. Old paths: lrr_projected/ -> total_change/.
```

```text
TOTAL shoreline change (2026-09-19, Hannah's advisor; renamed 2026-09-21):
per transect lrr_m_yr x (end - start) years, the rate fitted on that same
window, beside the OBSERVED change between calendar-year mean positions at
the two ends. Windows 1996_2024, 1996_2010, 2010_2024. Written by
scripts/input_prep/5-scr/3-rates/coastsat/total_change/coastsat_total_change.py.
```

```text
PROJECTED shoreline change (2026-09-21, Hannah by interview): the 1996-2024
LRR carried over each 14-yr half, against the same observed change. Only
1996_2010 and 2010_2024 exist -- 1996_2024 would BE the total change above.
Written by the same script, --product projected.
```

```text
WINDOW CONVERGENCE (2026-09-23, Hannah by interview): which windows recover
the long-term rate, and which are too short? The SAME OLS as lrr/, run on two
families of NESTED windows that both converge on the 1996-2024 rate from
opposite sides -- forward pins the start at 1996 and walks the end out,
backward pins the end at 2024 and walks the start back. The pair brackets the
answer: forward gives a window 1996-YYYY, backward a window YYYY-2024. Run at
two scales, eight evenly spaced domains and every transect on the island.

WHY IT IS A PRODUCT AND NOT A COMPARISON. 4-comparisons reads stored tables
and never refits; this refits 200 times, so it belongs in 3-rates as another
reading of the CoastSat record (the rule 4-comparisons/README states).
Nothing is graded against it -- it says whether the 14-year grading window
is long enough, not what the target is.

A NESTED SWEEP CONVERGES BY CONSTRUCTION: the last window IS the reference,
so the result is the SHAPE of the approach and the year the curve settles,
never a yes/no match. Written by scripts/input_prep/5-scr/3-rates/coastsat/
window_convergence/coastsat_window_convergence.py.
```

```text
QUESTION FIRST (Hannah, 2026-09-29, "it is hard to understand"). The tree
splits by the question each product asks, numbered in reading order, and the
pinned direction sits under each question:

1-rate_profiles/<direction>_from_<year>/     what does each window's WHOLE
alongshore profile look like?
2-r_bias_rmse/[<direction>_from_<year>/]     how close is it to 1996-2024?
r, bias, RMSE (added 2026-10-01)
3-settling_window/<direction>_from_<year>/   when does each LOCATION
a-eight_sites/  b-every_transect/        settle on it?
c-domain_means/
experiments/record_cut_<end>/<direction>_from_<year>/
the settling sweep on a
truncated record

Until 2026-09-29 it was record_<start>_<end>/<direction>/{sites,all_transects,
domain_means,alongshore_profiles}/, which put two questions at one level and
the 2020-truncated experiment beside the main result as an equal.
```

```text
The three scales of the settling sweep, keyed as its --scale choices are.
Lettered so they list in reading order: the readable case, the island, and
the unit the model is graded on.
```

```text
The CoastSat and dune-line net change side by side, 1996-2024 and its
halves (net_change_vs_duneline.py, 2026-09-18): the net-change chain figure.
```

```text
TOTAL shoreline change against the dune line's measured net change, per
domain, one folder per window: 1996_2010, 2010_2024 and 1996_2024, each
window's rate fitted on that same window, x the CALENDAR span (14, 14, 28
yr). Plus chains/ and difference/. Written by
scripts/input_prep/5-scr/4-comparisons/shoreline_vs_duneline/total_change_vs_duneline.py.

ONE tree since 2026-09-21 (Hannah, by interview). It absorbed two folders
that were the same quantity under two names:
lrr_net_change/      the halves, each on its own rate -- renamed, it was
never a projection and its docstring said so.
projected/1996_2024/ the 1996-2024 rate x the 25.72 yr dune-line interval.
Not a projection either: same fit window, shorter
span. Its numbers were already carried here as the
*_dune_interval_m columns (verified identical to
0 m before the merge), so the headline figure is the
28 yr calendar span and the dune interval is a column
and a caption line. Retired to superseded_20260921/.
```

```text
The same comparison with the shoreline side PROJECTED instead: the
1996-2024 LRR carried onto each 14-yr half, against the dune line measured
over that half (2026-09-22, Hannah). 1996_2010 and 2010_2024 only -- over
the full period the rate window IS the change window, which is
TOTAL_CHANGE_VS_DUNELINE/1996_2024. Written by the same script,
total_change_vs_duneline.py --product projected.
```

```text
Where the dune line sat in 1997, 2009 and 2023: maps, imagery zooms, and
its distance to NC-12 and to the CoastSat shoreline (duneline_positions.py,
2026-09-18).
```

```text
The same period's mean shoreline over two windows -- the 3-yr calendar span
and the 2-yr span centred on the start DEM's lidar flights -- differenced per
transect and domain (coastsat_mean_shoreline_windows.py, 2026-09-29). One
folder per period start: mean_shoreline_windows/<period>/.
```

```text
THE DETRENDED POSITION (2026-09-23). Every CoastSat transect detrended against
its own 1996-2024 fit and averaged: the signal the whole island shares, and
the tests of what causes it. Built to explain a result in
3-rates/coastsat/window_convergence_1996_2024/ -- that the convergence window barely
varies between transects, which a per-transect cause cannot produce.

It is NOT a rate and NOT a model input, so it files under 1-observations with
the record it describes rather than under 3-rates with the targets. Resolved
here because two scripts share its matrix and the window_convergence README
points at it. Written by
scripts/input_prep/5-scr/1-observations/detrended_position/.
```

```text
Outputs deleted 2026-09-19 (stale, on the old 1984/2004 periods); the
scripts would regenerate here.
```

```text
THE MEAN SHORELINE (2026-09-22). One averaging window's mean satellite
shoreline, as a line on the ground: each CoastSat transect's chainage
averaged over the window and placed back in space, then the ~906 mean points
strung into a single polyline. Built by
scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline.py.

WHY IT IS A PRODUCT AND NOT AN INTERMEDIATE. Every other CoastSat product
here is a DIFFERENCE of chainage, in which each transect's arbitrary origin
cancels. This one is a POSITION, so the origin does not cancel and the
geolocation step is real work: aggregated to the 90 domains, raw chainage
spans 124 m alongshore while the geolocated position spans 6222 m, because
the transect origins follow the shore around the cape. The line is read
outside its producer -- 2-brie-offset turns it into the shoreline-derived
island offset -- so it is resolved here rather than typed there.
```

```text
The four windows drawn on one y axis (coastsat_lrr_windows.py). Since
2026-09-19 only its 2 x 2 is drawn, into 3-rates/coastsat/lrr/; the
per-window and halves figures duplicated 3-rates and were archived with
the old folder (archive/2026-09-19_4-comparisons_duplicates/, deleted
2026-10-01; recover from git).
```

```text
RETIRED 2026-09-19: duneline_windows.py duplicated 3-rates/duneline/
endpoint (window + chain figures); its output is archived.
```

```text
Two-window comparison figures (coastsat_two_period_comparison.py). Outputs
deleted 2026-09-19 (stale, old periods); the script would regenerate here.
```

```text
Retired windows (1978-1997, 1997-2019 and their specific-dates variants),
kept for the DSAS comparison in 6-scr-smooth, never a grading target.
Archived, but still read, so still resolved.
```

```text
6-scr-smooth: what the LOWESS smoothing does to the observed rates. Resolved
here too because its outputs are read outside their producers (2026-09-18,
when the two folders lost their HAT_*_output names).
method_comparison/   transect-based against domain-averaged smoothing;
03_cascade_inputs/ is read by the overwash work
dsas_vs_coastsat/    the two rate sources, both smoothed, on the
retired 1978-1997 / 1997-2019 windows
```

```text
Single files in transect_domains/ that scripts outside 5-scr read by name.
HAT_domains.json holds the 90 real domain boxes; the older 1000 m polygons,
now in map_elements/archive/, are NOT the model domains.
```

```text
Folders under coastsat_lrr/ that are not windows. Listed so windows() can
report what IS available without having to parse every directory name.
Everything but "custom" moved out on 2026-09-18; the names stay listed so a
stray copy restored from an old checkout is not mistaken for a window.
```

```text
THE EXTENSION (2026-09-16, the Pea Island extension experiment). The
transects beyond GIS 1-90, numbered by their whole-island domain polygon
(hat_extension_domains.join_origins), get their own lookup and their own rate
fit per window, beside the surveyed ones and never merged into them:

transect_domains/transect_domain_lookup_ext.csv
coastsat_lrr/<window>/ext/transect_lrr_full.csv       extension only
coastsat_lrr/<window>/ext/transect_lrr_with_base.csv  surveyed + extension

The with_base file is what an extended-geometry run loads as its active
dataset (one LOWESS over the whole reach); the surveyed file stays the
scoring table for GIS 2-89 so extended and base runs are graded alike.
Built by scripts/input_prep/5-scr/3-rates/coastsat/extension/coastsat_extension_lrr.py.
```

<details><summary>Function notes (the original docstrings)</summary>

**`window_convergence_dir()`**

```text
One settling sweep's folder: `forward_from_1996` or `backward_from_2024`.

NOT `<start>_<end>`, although rule 2 would ask for it: a span name claims
ONE interval and each of these folders holds twenty-five of them. What
they share is the end that is PINNED, so that is what the name gives.

The full 1996-2024 record files under `3-settling_window/`. Any other
record span is a DIFFERENT EXPERIMENT, not a version of the same product,
because every window in it is fitted against a different reference, so it
files under `experiments/record_cut_<end>/` (or `record_<start>_<end>/` if
the start moved too). Filing them apart also stops a truncated forward run
overwriting the full one, which shares its pinned year and so its folder
name (Hannah, 2026-09-23).

Two directions, because the pair brackets the answer (Hannah, 2026-09-23).
Forward pins 1996 and walks the end year outward: how much record do you
need from the start of the chain? Backward pins 2024 and walks the start
year back: how late can a window begin and still recover the long-term
rate? Both converge on the same 1996-2024 reference from opposite sides.
```

**`window_profiles_dir()`**

```text
One rate-profile family's folder, under `1-rate_profiles/`. Full
1996-2024 record only; the profiles have no truncated experiment.
```

**`window_scores_dir()`**

```text
The r / bias / RMSE scores of the rate profiles, under `2-r_bias_rmse/`:
the folder itself (both directions together) or, given a direction, that
direction's folder. Full 1996-2024 record only.
```

**`mean_shoreline_label()`**

```text
The token a window's folder and files are named with.

Calendar years give `1995_1997`. Since 2026-09-29 a window can also be
given by ISO dates -- the +/-1 yr windows centred on a start DEM's flights
-- and then the token is the dates, `1995-10-12_1997-10-12` (Hannah: name
them by exact dates, beside the calendar folders).
```

**`mean_shoreline_dir()`**

```text
One averaging window's folder: the line, the per-transect means, the
figure and PROVENANCE.md. Named for the span, per rule 2. `start`/`end`
are calendar years or ISO dates (see mean_shoreline_label).
```

**`mean_shoreline_geojson()`**

```text
The mean shoreline as ONE LineString, EPSG:26918, carrying the same
metadata properties a digitised dune line carries so the 2-brie-offset
intersection step reads it unchanged.
```

**`mean_shoreline_csv()`**

```text
Per CoastSat transect: the window mean, its scatter and count, the
dates it spans, and the geolocated mean point.
```

**`window_dir()`**

```text
The folder holding one window's rate fit.

Raises:
    FileNotFoundError: If that window has not been built, naming the ones
        that have. A window is built by
        scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_domain_lrr.py.
```

**`coastsat_endpoint_csv()`**

```text
Net CoastSat shoreline change for one window, between +/-6-month means
about the dune-line survey dates: `level` "transect" or "domain". Raises,
naming the producer, if absent.
```

**`dune_endpoint_csv()`**

```text
Net dune-line change for one window: `level` "transect" (per 100 m
transect) or "domain" (per GIS domain). Raises, naming the producer, if
absent.
```

**`lrr_csv_ext()`**

```text
The extension's rate fit for one window (with the surveyed reach by
default), raising if the extension has not been built for it.
```

</details>

### hat_overwash.py

Where does the observed overwash record live, and its figures and tables?

Notes that were in the code:

```text
hat_overwash.py

WHERE DOES THE OBSERVED OVERWASH RECORD LIVE, AND ITS FIGURES AND TABLES?

WHY THIS EXISTS
overwash_data.py already held the root for its three sibling scripts, but
the two model-comparison scripts in scripts/analyze_output/overwash/ typed
the workbook's path themselves, and both named a location it had left
(scripts/input_prep/8-overwash-analysis/), so neither could run. Same
answer as hat_observed_rates.py for 5-scr: resolved ONCE (2026-09-18).

THE LAYOUT, grouped by job (2026-09-18). Before that it was grouped by file
type -- figures/ and tables/ -- which split the footprint comparison across
figures/vs-footprint/ and a top-level vs-footprint/.

data/hatteras_init/8-overwash-analysis/
README.md
CAPTIONS.md               one entry per figure, headed by its folder
1-observations/           THE RECORD, and its long-form tables
Hatteras_Overwash_Data.xlsx
overwash_observations.csv  storms_by_image.csv
2-record/                 what the record shows
heatmaps/  map/
3-vs-footprint/           the record against the 1984 footprint:
figures here, their tables in tables/
4-vs-model/               the record against modelled overwash in
the hindcast runs (2026-09-27): figures
here, tables in tables/
```

```text
The label each figure's CAPTIONS.md entry carries, keyed by the short
folder name the scripts pass to overwash_data.upsert_caption. Its order is
the order the entries are sorted into.
```

### hat_source_sink.py

Where does the source/sink (background erosion) calibration live?

Notes that were in the code:

```text
hat_source_sink.py

WHERE DOES THE SOURCE/SINK (BACKGROUND EROSION) CALIBRATION LIVE?

WHY THIS EXISTS
Seven scripts under scripts/input_prep/7-source-sink/ built
data/hatteras_init/7-source-sink/... by hand, and two of them had to agree
on the same pass-0 backup file with a "keep the two in step" comment -- one
of which pointed at a file that did not exist. Same answer as
hat_observed_rates.py for 5-scr: the location is resolved ONCE (2026-09-18).

The field the MODEL runs on is not here: it is a dict literal in
hatteras_site_config.py. These are the fit's own products, its figures,
and the exported copy of the field.

THE LAYOUT, whose numbers match the script steps in
scripts/input_prep/7-source-sink/ (1-prepare writes into output/raw_runs/,
so it has no folder here)

data/hatteras_init/7-source-sink/
README.md                 generated by 4-export; describes the field
2-calibrate/
<pair>/               what one fit wrote: DOMAIN_BE_RATES*.txt,
be_zone_metrics.csv, cascade_base_lrr.csv,
convergence_history.json
prebe/                hatteras_site_config backups taken before
each apply pass
3-figures/
<pair>/1-field/  2-method/  3-limits/
4-export/                 be_rates_<period>.py and
be_calibration_domains.csv (default pair)
archive/                  superseded_20260914/ (_20260825/ and _20260902/
                          deleted 2026-10-01)

A PAIR IS TWO CALIBRATION PERIODS FITTED JOINTLY, named
<p1start>_<p1end>__<p2start>_<p2end>, e.g. 1984_2004__2004_2024. Until
2026-09-18 the default pair wrote to the ROOT of 2-calibrate/ and 3-figures/
and only other pairs got a folder, so what the unlabelled files belonged to
had to be known rather than read. Every pair has a folder now.
```

```text
Config backups written by be_apply_fit_to_config.py before each pass.
Shared by every pair: a backup is of the whole config, not of one fit.
```

```text
The field as it stood before PASS 1 of the frozen-zone masked iteration,
i.e. the one-shot (pass-0) solve. plot_be_zones.py and
export_be_calibration.py both split each final rate into pass-0 plus
what the iteration added, and must read the same file to agree.
```

### hat_topo_version.py

Which Barrier3D domains does a script read: which product, and which version within it?

Notes that were in the code:

```text
hat_topo_version.py

WHICH Barrier3D domains does a script read - which PRODUCT, and which VERSION
within it?

WHY THIS EXISTS
Four road scripts used to hardcode ".../2009-dune-topo/2009_v3/topography".
When the dune windows were re-picked into 2009_v4 (2026-08-19) they kept
reading v3 interiors while consuming v4 setbacks. Nothing errored. 18
domains had different interiors, 10 of them a different SHAPE - including
D79 and D80, two of the three roadways the relocation logic acts on. Every
drown verdict and placement number for those was computed on the wrong grid.

So the location is resolved ONCE, here, and a name that does not exist on
disk is an immediate, loud error rather than a silently stale read.

WHAT CHANGED 2026-08-25 - THERE ARE NOW TWO PRODUCTS
The tree was stage-keyed and held exactly one topography, which BOTH
hindcast periods read:

1-barrier3d-domains/2009-dune-topo/2009_v5/{topography,dunes}

It is now PERIOD-keyed, because the two periods start from different DEMs:

1-barrier3d-domains/
1984-start/    from DEM 2009-2014-1996   dune-topo/<version>/
2004-start/    from DEM 2009-2014        dune-topo/<version>/
forecast/      from a 2025 DEM, later
buffer/        shared
superseded/

2009_v5 became 2004-start/dune-topo/v1 (renamed from v5 on 2026-08-26 so
version numbers restart per product) - it was always built from
the 2009+2014 DEM, which is what the 2004 start uses. v3 and v4 are under
superseded/ (they were picked against the UNFILLED DEM).

WHY topo_dirs() STILL DEFAULTS TO A PRODUCT
Every existing caller - the road tree, the groin sweep, the poster script -
calls topo_dirs() with no arguments. Defaulting to DEFAULT_PRODUCT keeps
all of them resolving exactly what they resolved before the restructure, so
this change moves no road number and no published figure. The RUNNER is the
one caller that now passes a product, from HATTERAS_PERIODS[start]
["topo_product"].

SETTLED 2026-08-26. That note used to read "a live question": the road
scripts measured setbacks against 2004-start for BOTH periods. They no
longer do. Every vintage now resolves its own product through YEAR_PRODUCT
below, and RoadSetback_1984_dunestart.csv was re-measured on 1984-start.
What made it urgent rather than tidy: all 90 domains differ between the two
products and 65 have a different interior SHAPE.

HOW THE VERSION IS CHOSEN, in order
1. an explicit override= argument
2. HAT_TOPO_VERSION_<PRODUCT> in the environment, e.g.
HAT_TOPO_VERSION_1984_START=v5. Product-scoped, because
both products have a "v1" and a global override would silently resolve to
a real but wrong directory for the other period.
3. a CURRENT file in the product's dune-topo/ directory
4. the extractor's VERSION, but only if the extractor is currently pointed
at the SAME product
5. the only version present, if there is exactly one
otherwise: raise, listing what is on disk

CURRENT OUTRANKS THE EXTRACTOR LITERAL (swapped 2026-09-04). The extractor's
VERSION says what the extractor WRITES; CURRENT says what everyone READS.
They used to be one literal doing both jobs, which meant the only way to
make a layered version (the 1984-start layers v3-v8 of 2026-09-04 -- built
ON v2, never BY the extractor; DELETED 2026-09-07, only unmodified
extractions are kept) the
default was to edit the extractor to a name it would then overwrite on its
next run. So CURRENT existed, recorded intent, and was inert; the 1984-start
README carried a paragraph explaining that it did nothing. Now a product
with a CURRENT file reads that version, and a product without one still
follows the extractor ("bump VERSION and the tree follows" holds where no
CURRENT has been written - 2004-start today). A fresh extraction with
CURRENT still naming the old one is therefore NOT adopted until CURRENT is changed,
and that is the point: extracting and adopting are two decisions.

THE ENV RULE OUTRANKS THE EXTRACTOR (added 2026-09-02, after it cost a run).
A batch that selected its arm by writing CURRENT was ignored, because rule 3
fired first and the extractor was sitting on the same product. Two arms of a
three-arm experiment silently duplicated the control, exit code 0. CURRENT
is a persistent shared DEFAULT; per-run selection needs something that does
not mutate state every other reader sees.

USAGE
from site_layer.hat_topo_version import topo_dirs
TOPO_DIR, DUNE_DIR, RUN_NAME = topo_dirs()                  # 2004-start
TOPO_DIR, DUNE_DIR, RUN_NAME = topo_dirs("1984-start")
TOPO_DIR, DUNE_DIR, RUN_NAME = topo_dirs("2004-start", override="v1")
```

```text
Shared across products: the padding domains Barrier3D needs either side of
the 90 real ones. Not per-period - the buffer is not a survey of anything.
```

```text
The other shared inputs (2026-09-18). About twenty-five scripts typed these,
or re-joined the product paths the helpers below already return.
domain-clips-1m/<domain_N>/   clip_domain_N.tif (1 m, native) and
resampled_domain_N.tif (10 m, the Barrier3D
grid every measurement georeferences against)
control-picks/                window sets kept across a version clear
npy-arrays_2009_unfilled/     the pre-gap-fill arrays; still read by
HAT_rasterize_road_to_domains.py as a grid check
domain-geojson/               the 120-domain Pea-Hatteras polygons
```

```text
What topo_dirs() resolves when no product is named. See the note above: this
is the pre-restructure behaviour, kept so nothing that was not asked to move
moves.
```

```text
WHICH PRODUCT DOES A HINDCAST PERIOD READ. The single definition (2026-08-26).

Keyed by the period's START YEAR, which is also the year that labels every
road vintage in 4-mgmt-forcings: RoadSetback_1984_dunestart.csv is measured
from row 0 of the 1984-start extraction, RoadSetback_2004_dunestart.csv from
row 0 of 2004-start.

WHY IT LIVES HERE. It was written out three times and omitted five times.
HAT_road_offset_from_dune_start.py carried its own YEAR_PRODUCT literal,
HAT_road_setback_audit.py spelled the product per scenario, and
hatteras_site_config.py carries it per period as "topo_product" -- while
HAT_road_placement_on_domains.py, HAT_road_method_diagnostic.py and
HAT_road_domain_views.py looped over BOTH years against a single
module-level topo_dirs(), i.e. DEFAULT_PRODUCT for both.

That last one is not a cosmetic duplication. Between the two products ALL 90
domains differ and 65 have a different interior SHAPE -- GIS 11 is 165 rows
on 1984-start and 157 on 2004-start. It is the v3/v4 failure this file was
written for, four times larger, and it produced no error: a 1984 setback
scored against a 2004-start interior is simply a different island.

So the pairing is defined once, here, beside the resolver that consumes it.
hatteras_site_config.py imports it rather than repeating it; every road
script resolves through it rather than defaulting.
```

```text
The two periods added 2026-09-11 SHARE these products rather than owning
one each. A product is a DEM composition, not a period label: 1984-start
is the surface carrying the 1996 ALACE graft, which is the survey nearest
a 1996 start, and 2004-start is the 2009-plus-2014 surface, which is the
closest vintage match of any period to a 2010 start.

So a product name now names the period it was FIRST built for, not the
only one that reads it. Ask this mapping rather than inferring a product
from a year, and ask product_for_year() rather than indexing directly.
```

```text
WHICH NC-12 LINE A PERIOD'S ROAD IS MEASURED FROM. The single definition
(2026-09-15).

The two digitised centrelines under 4-mgmt-forcing/road_offset/raw_offset/
were exported off 1978 and 2008 imagery. Until 2026-09-15 their folders, the
rasterised masks under raster/, and the relocation measurement under
road_relocation/ were all named 1984 and 2004 -- the period starts they stand
in for -- which put two different axes under one integer: dunestart_offset/
1984 meant "the 1984 START", raw_offset/1984 meant "the 1978 LINE". They are
now named by the line's TRUE vintage, and this map is the only place a start
year is paired with a line. Ask road_line_for_year() rather than re-spelling
the pairing; a period start passed where a line vintage is expected is a
loud error (road_line_file), not a missing-file error later.

1996 and 2010 have no line of their own. 1996 is the 1978 line with the 1989
Pea Island relocation applied; 2010 is the 2008 line unchanged, because no
relocation falls between 2004 and 2010. See ROAD_SETBACK_KIND.
```

```text
WHICH SETBACK FOLDERS HOLD A MEASUREMENT AND WHICH A DERIVATION. The tree
under road_offset/dunestart_offset/ is split the same way (2026-09-15):

dunestart_offset/measured/<year>/   measured on the period's line against
row 0 of the period's own extraction
(HAT_road_offset_from_dune_start.py)
dunestart_offset/derived/<year>/    built FROM a measured file
(HAT_road_setback_derived_vintages.py)

so a derived file cannot be mistaken for a measurement by its address alone.
Before the split all four sat as flat siblings and only PROVENANCE.md said
which two were copies.
```

```text
THE REST OF 4-mgmt-forcing (2026-09-18). About thirty scripts typed these
themselves, and four still spelled old_method_offset/, which became a dated
superseded folder on 09-11 -- so they failed soft, comparing against nothing.
The layout is the one settled on 2026-09-15; this only names it once.
```

```text
Today's NC-12, from the NCDOT route inventory (2026-09-18), standing in for
2023/2024. NOT in ROAD_LINE_FOR_YEAR: no hindcast reads it; the dune-line
position figures do. See raw_offset/current/PROVENANCE.md.
```

```text
The old-method setbacks ("old_method_offset/" in older scripts), kept for the
method comparison: <year>/RoadSetback_<year>.csv.
```

```text
The 1984 dune-start setbacks as measured on 1984-start/v1, before the road
tree was re-measured on v2 ("dunestart_offset_ARCHIVE_1984start_v1").
```

```text
THE DUNE LINES, BY VINTAGE (2026-09-15). The same rule as the road lines:
a raw per-transect file under 2-brie-offset/raw_offsets/ is named for the
IMAGERY VINTAGE of the digitised line it came from, and this map is the only
place a period year is paired with a vintage. Until 2026-09-15 the 1996
start read a byte-identical COPY of the 1997 file filed under the name
1996_...; the copy is gone and the pairing lives here (Hannah: "a
DUNE_LINE_FOR_YEAR table, no copies"). When the 2010 and 2024 lines are
digitised, add them here under the year of their IMAGERY (2009 or 2023, if
that is what they are traced from), never under the period year.

A vintage may have more than one digitisation (duneline_1997.geojson and
duneline_1997_v2.geojson); the raw file under the vintage's name is the
CURRENT one, and each build under 2-brie-offset/<year>/v<n>/ keeps a copy of
the raw it was made from, so an older build is always reproducible.
```

```text
The rest of 2-brie-offset (2026-09-18): about fifteen scripts typed these,
four of them against layouts that no longer existed (the flat per-start
build, hindcast_<year>/ folders, a dunelines/ under 1-barrier3d-domains).
```

```text
WHICH FEATURE THE OFFSET WAS MEASURED FROM (2026-09-22). Until then every
build came from a digitised DUNE line, so the source was not worth naming.
The 1996 start now also has a build from the CoastSat satellite SHORELINE
(the 1995-1997 mean position; scripts/input_prep/5-scr/1-observations/
mean_shoreline/), which is a different feature, not a newer reading of the
same one -- so it is a separate source with its own v1, never a v2 of the
dune build (Hannah, 2026-09-22, and the rule in [[feedback-version-numbering-restarts]]).

The dune build KEEPS the flat layout it has always had, <year>/v<n>/, so
nothing the runner resolves moves. A non-default source nests one level
deeper, <year>/<source>/v<n>/, and because "shoreline" does not match the
v<n> pattern offset_version() scans for, adding it cannot disturb which
build the runner reads.
```

```text
A build's three files share one stem, named for the FEATURE the offset was
measured from, so a copy that leaves its folder still says what it is -- the
basename is the only thing a file carries with it.
```

```text
The default source keeps the plain env key it has always had, so an
override written before the sources were split still selects the build
it always selected. A non-default source gets its own key, so
overriding the shoreline arm cannot silently move what the runner reads.
```

```text
offset_YEAR_dir, not offset_start_dir: a comparison between two sources
belongs to neither, so it must not inherit one source's folder. It did
for a few minutes on 2026-09-22, when offset_start_dir started nesting
and this quietly followed it down into duneline/.
```

```text
THE SATELLITE SHORELINE, BY WINDOW (2026-09-22). The dune-line offset is
measured from a line digitised on ONE day, so DUNE_LINE_FOR_YEAR pairs a
period with a vintage YEAR. A CoastSat shoreline has no such day: a single
satellite pass carries metres of tide, wave setup and cloud-edge noise, so
the position a period starts from is a MEAN over a window of passes. The
pairing is therefore a period year -> a window, and this is the only place
it is spelled.

SINCE 2026-09-29 each window is +-1 yr of the middle of the lidar flights
the start DEM is built on (Hannah, by interview): the offset is a snapshot
the model starts from beside that DEM, so it is dated like the DEM.
1996: the 1996 ALACE lidar, flown 1996-10-09..16 -> 1995-10-12..1997-10-12
2010: the 2009 USACE NCMP lidar, flown 2009-08-10..24 -> 2008-08-17..2010-08-17
These are the windows of the shoreline v2 builds (CURRENT since the same
day). Only the OFFSET reads this table -- island_offset_hybrid's default raw
file -- so it departs from [[cascade-period-is-the-calendar-year]] for the
offset alone; observed change, rates and scoring keep their calendar windows
and never read it. The anchors live in coastsat_mean_shoreline.SURVEY_ANCHORS.

Until then: 1996 -> calendar 1995-1997 (2026-09-22), 2010 -> calendar
2009-2011 (2026-09-28); those built shoreline/v1, whose raw files are kept
in each v1 folder.
```

```text
a window is two calendar years (1995, 1997) or two ISO dates, and the
file is named as the mean_shoreline folder is (mean_shoreline_label)
```

```text
THE ARRAYS HAVE NO YEAR TAG (2026-08-26).

They are domain_<N>_topography.npy / _dune.npy / _nodata.npy. There is no
year in the name, and there should not be one.

A tag was tried, briefly, both ways. The name was the literal "2009" for a
long time, which was simply false - 2004-start is the 2009+2014 mosaic and
1984-start is 2009+2014+1996, so neither product is a 2009 DEM. Retagging by
period ("1984"/"2004") was then tried and reverted the same day, because the
tag turned out to be a DISTRIBUTED INVARIANT with no single place to fix and
no single grep to audit. Twelve live scripts build these paths, and the tag
reached them four different ways:

TOPO_DUNE_INIT_YEAR = "2009" then interpolated      5 scripts
the bare literal inline in the f-string             2 scripts
ext.TAG, imported from the extractor                1 script
globbed away as domain_*_topography_*.npy           2 scripts

The retag found the first form, missed the second, and broke both of those
scripts - they resolved the right DIRECTORY and then asked for a file that no
longer existed. That is the argument in one sentence.

WHAT THE TAG COULD NOT DO ANYWAY. It cannot catch a period mix-up. The tag
and the directory both derive from the same `product`, so a wrong product
gives a wrong directory AND a matching wrong tag - consistent, silent, no
error. What guards that is the runner's boot/run product assertion, not the
filename. And for the two scripts that glob the tag away, a stray file from
the other period would make the glob match twice and pick arbitrarily, which
is worse than no tag at all.

The period lives in the DIRECTORY, which every caller has to get right
regardless, and in each run's RUN_MANIFEST.txt. The buffer arrays have never
carried a tag and have never been confused.

CALLERS SHOULD NOT BUILD THESE NAMES. Use domain_arrays() or array_path()
below, so the directory and the filename come from one place.
```

```text
The one that has ALONGSHORE_FLIP = True. Three other copies of this file exist
in the repo and all are unflipped -- see the note in
HAT_road_offset_from_dune_start.py.
```

```text
THE TWO HALVES OF STAGE 1 (2026-09-09, Hannah). Under each product,
`1-extraction/` holds what the extractor reads and records (npy-arrays,
npy-arrays_survey, picks, the aerial review of the holes, the retired
experiments) and `2-domain-reconstruction-1984/` the 1984 reconstruction in
its six steps; `dune-topo/` stays at the product root because BOTH halves
write versions into it and it is what the runner loads. The scripts folder
scripts/input_prep/1-barrier3d-domains/ is split the same way.
```

```text
THE SEAWARD-ROW-INSERT FOLDER, and the paths that hang off it.

ONE definition, because eight plotting scripts and two measurement scripts
used to build these by hand - and two of them WRITE.

THE LAYOUT IS NOT SYMMETRIC BETWEEN PRODUCTS, deliberately. On 2026-09-03
everything belonging to the 1984-start seaward-row insert - the measurement of
N, the scope report, the fill comparison and every figure - was consolidated
under `2-domain-reconstruction-1984/`. 2004-start has no insert work and no such folder,
so its dune-line measurements stay at the product root.

The asymmetry is the price of that consolidation. It is contained here so a
caller cannot get it wrong, and so a future product does not inherit it by
accident: anything not listed gets the plain layout.
```

```text
The four sections of the insert figures folder, in the order the argument
runs: where the insert lands, where N came from, what the rows are made of,
and what the result looks like. Numbered so a directory listing reads in that
order, matching the numbered layout of data/hatteras_init itself.

WHY THIS IS HERE AND NOT IN THE PLOTTERS. Thirteen figures in one flat folder
is a dump, and a folder a plotter re-scatters on every run cannot be tidied by
moving files. Naming the section at the call site - and resolving it here - is
what makes the layout survive a re-plot. A section that is not one of these is
a typo, and raises rather than silently creating a new folder.
Renumbered 2026-09-07 to the order the argument runs: measure the dune-line
shift, turn it into a footprint of rows (two placements of the same rows:
seaward/, behind-road/; placement-independent figures at the step's root),
argue the fill, look at the result. 6-result is reserved: nothing has been
built and run on the footprint yet. The record figures of the deleted layers
sit in superseded-layers/ and the irreproducible pre-re-pick ones in frozen/;
neither is a section a plotter may write to (2026-09-08).
```

```text
Inside a section (2026-09-08, Hannah): island-wide figures go in `island/`
(or a named subfolder), and EVERY figure that shows one example domain goes
in `rows-added/` or `rows-removed/` by the sign of that domain's N in the
footprint table - never at the section root.
Six steps since 2026-09-09, in the order the ARGUMENT runs (Hannah): how far
the dune line moved, how many rows, WHERE they go (the two candidate
placements, the road check, and the imagery review that decides), what they
contain, the version built, what the model does. 5-build holds no figures.
```

```text
The DATA of the same steps (2026-09-09, Hannah: "organize this under
subfolders"): the tables and reports each step writes sit in a subfolder of
2-domain-reconstruction-1984/ named like its figure section, so a listing of the data
folder reads in the same order as figures/. Root keeps README.md,
DUNE_TOPO_VERSION_GUIDE.md and figures/. Resolved here for the same reason
the figure folders are: one definition, no hand-built paths in the scripts.
```

```text
"bridged" is written only by nodata_audit/HAT_bridge_dropouts.py: True where
an unsurveyed cell was filled by interpolation between measured neighbours.
It is a THIRD state, not a replacement for "nodata" - a bridged cell is still
a cell no survey saw, and the nodata mask keeps saying so.
```

```text
A domain outside the surveyed reach (an extended geometry, 2026-09-16)
has no array of its own and runs on the buffer profile, exactly as the
padding does; what makes it different from padding is its offset.
```

```text
THE ENVIRONMENT OUTRANKS THE EXTRACTOR LITERAL, and it has to.

Added 2026-09-02, after it cost a run. The crest experiment selected its
arm by writing dune-topo/CURRENT, ran, and reported dune-topo\v1 -- the
baseline -- because rule 2 below reads the extractor's VERSION literal
FIRST and the extractor happened to be sitting on the same product. The
CURRENT file was never consulted. Two arms of a three-arm experiment
were duplicates of the control, exit code 0, no warning.

That is this module's own failure mode, one level up: a caller that
cannot say "use THIS version" without editing a source file will end up
editing a source file, or will think it said it and be ignored. CURRENT
is not usable for that -- it is a persistent, shared default, and a batch
that sets it per arm is mutating global state for every other reader.
```

```text
CURRENT BEFORE THE EXTRACTOR LITERAL (2026-09-04) - see the header. The
literal is what the extractor writes; CURRENT is what is read. A product
without a CURRENT file behaves exactly as before.
```

<details><summary>Function notes (the original docstrings)</summary>

**`legacy_setback_file()`**

```text
The OLD-METHOD setback for a measured start (1984, 2004) -- superseded,
kept only to compare the two methods.
```

**`road_setback_relpath()`**

```text
road_setback_file() relative to INIT_ROOT, POSIX-style -- the form
HATTERAS_PERIODS[year]["road_setback_file"] carries.
```

**`offset_basename()`**

```text
The stem a build's three files share, e.g. Island_Dune_Offsets_1996.
island_offset_hybrid.py names its outputs from this rather than keeping a
second spelling of them (2026-09-22).
```

**`offset_start_dir()`**

```text
One source's folder for one period start: PROVENANCE.md, CURRENT,
v<n>/ builds, and any superseded_*/ and ext/.

EVERY source nests under its own name (2026-09-22, Hannah: "make it clear
with the naming of each folder what the original source was"). The dune
builds sat flat at <year>/ until then, from when they were the only kind,
which left a folder listing unable to say what <year>/v1/ was measured
from -- and left the two sources asymmetric once the shoreline arrived.

    <year>/duneline/v<n>/     from a digitised dune line
    <year>/shoreline/v<n>/    from a CoastSat window mean
    <year>/comparisons/       dune line vs the CURRENT shoreline build
    <year>/shoreline/<v>/comparisons/   dune line vs an older shoreline build
                              (both offset_source_comparison_dir; the
                              CURRENT pair moved up on 2026-10-06)

Nothing outside this module should join these parts by hand.
```

**`offset_year_dir()`**

```text
The period start's folder itself, which now holds only source folders
and comparisons/ -- no build of its own.
```

**`offset_version()`**

```text
Which build of a start's island offset every reader takes.

Order, mirroring topo_dirs(): HAT_OFFSET_VERSION_<year> in the
environment, then the CURRENT file, then the only v<n> directory present.
None for the flat layout (no v<n> directory at all). Several v<n> and no
CURRENT is an error, not a guess; so is naming a version that is not on
disk. Moved here from hatteras_site_config._island_offset_file on
2026-09-18 so the figure scripts that resolved it themselves share it.

A non-default source reads HAT_OFFSET_VERSION_<year>_<SOURCE> instead, so
overriding the shoreline arm cannot silently move the dune build the
runner reads (2026-09-22).
```

**`offset_comparison_dir()`**

```text
<year>/comparisons/<name>/ -- where a comparison BETWEEN builds lands.

A comparison is neither a version nor a source, so it does not belong in
either's folder. Until 2026-09-22 a version comparison was written into
the later version's own folder, which made a build's folder hold both the
build and a judgement about it; source comparisons are filed here from the
start. `name` says what was compared, e.g. "duneline_vs_shoreline".
```

**`offset_source_comparison_dir()`**

```text
Where a dune line vs shoreline comparison lands.

Against the CURRENT shoreline build: <year>/comparisons/<name>/, beside
both sources, since 2026-10-06 (Hannah: the shoreline and dune lines are
final, so the comparisons can sit at the year). Against any other
shoreline build: <year>/shoreline/<v>/comparisons/<name>/, filed with that
build as it has been since 2026-09-29, so a superseded comparison still
says which build it was drawn against. The builds used are written into
each comparison's README and caption.
```

**`offset_file()`**

```text
One file of a start's build. kind: padded (the model input), input
(90 values), or unpadded (90, with domain ids).
```

**`dune_line_for_year()`**

```text
The dune-line vintage a period start OR end year reads. With
`strict=False` an unknown year returns None instead of exiting, for the
end-year target, which is allowed to be missing.
```

**`dune_raw_file()`**

```text
raw_offsets/<vintage>_duneline_offset_raw.csv, the current build of
one LINE vintage's per-transect stations.
```

**`shoreline_raw_file()`**

```text
raw_offsets/<start>_<end>_shoreline_offset_raw.csv -- the per-transect
stations of one averaging window's mean shoreline, the shoreline
counterpart of dune_raw_file(). Written by duneline_to_raw_offsets.py from
the mean-shoreline geojson, which hat_observed_rates owns.
```

**`product_for_year()`**

```text
The topography product a period start year reads. Raises if unknown.

Use this in any loop over vintages. The failure mode it removes is a loop
body that resolves topo_dirs() once, outside the loop, and silently gives
every year the same interiors.
```

**`year_for_product()`**

```text
The period start year a topography product belongs to.

The inverse of product_for_year(), and it lives here for the same reason
the forward map does: the pairing is defined ONCE. A caller that needs
"which year is this product" -- the extractor's figure code does -- must
not re-spell {"1984-start": 1984} locally, because a third product added
to YEAR_PRODUCT would then update one direction and not the other.

strict=False returns None instead of raising, for products that are
legitimately not a hindcast period ("forecast", "buffer"). A caller using
it must have a defined behaviour for None -- see PRODUCT_YEAR in
HAT_dune_topo_extractor.py, which falls back to plotting every year.
```

**`_extractor_state()`**

```text
(TOPO_PRODUCT, VERSION) the extractor is currently configured for.

PARSED rather than imported: importing the extractor pulls in matplotlib
and a windowing backend, which is a heavy and fragile dependency for an
audit that needs neither. Only two literals are wanted.
```

**`insert_scope_step()`**

```text
The data folder of one step of the insert work (tables, reports), or a
named subfolder of it (3-placement has `imagery-review`). Created if
absent. `step` is one of INSERT_SCOPE_STEPS; anything else raises.
```

**`insert_figures_dir_for_domain()`**

```text
Where a figure of ONE example domain goes: the section's rows-added/ or
rows-removed/ (or unchanged/) by the sign of that domain's N, optionally
under a named subfolder (`under="road-check"`, `under="imagery-review"`). Created if absent.
```

**`insert_figures_dir()`**

```text
Where insert figures are written, created if it does not exist.

`section` is one of INSERT_FIGURE_SECTIONS. Omit it for the folder root,
which holds only the README, CAPTIONS.md and the `frozen/` figures no
script can rebuild - no plotter should write there. `sub` is a subfolder
of the section (INSERT_FIGURE_SUBFOLDERS: `island/`, a placement, or a
sign folder - use insert_figures_dir_for_domain for the latter); anything
else raises, for the same reason a wrong section does. Section roots hold
no figures either (2026-09-08).

IT MAKES THE DIRECTORY. Not a pure lookup, deliberately: none of the eight
plotters calls mkdir, so before the sections existed they all depended on
`figures/` already being on disk and every one of them would have failed on
a fresh checkout. Doing it in the one place that knows the layout is one
change instead of eight, and it cannot be forgotten by a ninth plotter.
```

**`array_path()`**

```text
Full path to one domain array, directory and name resolved together.

For the scripts that read a domain at a time - the road placement, the
setback audit, the units check. Use domain_arrays() when you want the
whole padded run.
```

**`domain_arrays()`**

```text
(elevation_paths, dune_paths) for a run, buffer-padded, as strings.

Replaces the build_domain_file_paths() that was copied verbatim into the
hindcast runner, the notebook and the groin sweep worker. Those three
copies each took an `init_year` and pasted it into the name; there is no
year any more and no name for a caller to build.

n_buffer is the padding on EACH side. The buffer profiles are shared by
every period and carry no tag, which is what these arrays now look like
too.
```

**`require_version()`**

```text
The version's directory, or a loud exit naming what IS on disk.

For scripts that carry a version LITERAL instead of resolving through
CURRENT -- the seaward-row-insert plotters name the layer they draw. On
2026-09-07 the 1984-start layers v3-v8 were deleted (Hannah's decision:
keep only unmodified topography), so a literal that was valid the day it
was written now names nothing. Failing here, before any array is opened,
says so in one line instead of a FileNotFoundError deep in a plotting loop.
```

**`env_override_name()`**

```text
The environment variable that pins ONE product's version.

Product-scoped on purpose. A bare HAT_TOPO_VERSION would apply to every
product, and BOTH products currently have a version called "v1" -- so a
global override set for one period would silently resolve to a real, wrong
directory for the other instead of failing.
```

**`current_topo_versions()`**

```text
{product: the dune-topo version it reads TODAY}, for every product.

What run_registry.rebuild_run_index needs to mark a run superseded: the
version each product resolves to now (CURRENT outranks the extractor
literal, see the rules above). A product with no dune-topo tree, or none
resolvable, is left out rather than guessed.
```

**`topo_dirs()`**

```text
(topography dir, dunes dir, run name), checked to exist.

Raises with what IS on disk rather than returning a path that will later
read as "no arrays for this domain".
```

</details>

### hatteras_site_config.py

Hatteras Island site config: the real domains, periods, forcing files, presets and management events.

From the script's original header:

```text
Hatteras Island site config: the real place names, spans, and labels.

This is application config, not library code -- it imports cascade_pipeline's
generic dataclasses and fills them in with Hatteras-specific content, the
same way any other CASCADE study site would. Nothing in cascade_pipeline itself
knows these values exist; a different site (Ocracoke, etc.) would write its
own sibling module in this same shape and never touch the package.

Import these presets from your run script / notebook:
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS, HATTERAS_ANNOTATIONS
```

Notes that were in the code:

```text
Sibling module in scripts/, which owns "where is data/hatteras_init". The
relocation cross-check below reads the period-2 setback file rather than
carrying a copy of its numbers, so this module needs the data root.
It also owns the year -> input pairings this file used to spell out:
YEAR_PRODUCT (topography), ROAD_LINE_FOR_YEAR (which digitised NC-12 line
a period's road is measured from) and road_setback_relpath() (where that
period's setback file is, measured/ or derived/). ROAD_LINE_FOR_YEAR is
imported here so a reader of this config sees the pairing by name.
```

```text
The other 4-mgmt-forcing paths this file names (road elevation, the
relocation measurement) come from the same module since 2026-09-18.
```

```text
Real domains GIS 1-90 (500 m each, south to north, Cape Point to Pea
Island), padded by 15 buffer domains on each side. These happen to match
DomainGeometry's own field defaults, but naming the instance explicitly
here (rather than relying on cascade_pipeline.domains.DEFAULT_DOMAINS) keeps
"this is Hatteras' geometry" visible at the call site.

THE REACH IS A NAMED GEOMETRY SINCE 2026-09-16 (the Pea Island extension
experiment). HAT_GEOMETRY in the environment picks one of
hat_extension_domains.GEOMETRIES; unset is "base", GIS 1-90, and every
matrix run. An extended geometry adds measured coast beyond GIS 90 (and
below GIS 1) on the shared buffer topography, with its own offset file
under 2-brie-offset/<year>/ext/<geometry>/ and its own CoastSat rows. Read
from the environment here, as HAT_OFFSET_VERSION_<year> is, because this
module is imported by scripts that never load HAT_hindcast_config; the
runner checks the two agree.
```

```text
WHICH FEATURE THE ISLAND OFFSET IS MEASURED FROM, for this run (2026-09-22).
"duneline" (the default, and what every run before this date used) or
"shoreline", the CoastSat window mean. Read from the environment for the
same reason HAT_GEOMETRY is: this module is imported by scripts that never
load HAT_hindcast_config.

It is a RUN-LEVEL choice, not a per-year one, because a run has one start
year. HAT_OFFSET_VERSION_<year> still picks the version WITHIN the source.

A run that moves this off the default is not a matrix run: it changes a
model input, so it needs HAT_RUN_KIND=experiment and a tag, or it derives
the same name as the matrix run it is being compared against and overwrites
it (the failure output/calibration/groin/README.md records for the rig
sweep). The runner refuses that combination rather than trusting it.
```

```text
The interior score is ALWAYS GIS 2-89 against the surveyed CoastSat table,
whatever the geometry, so an extended run and its 90-domain baseline are
graded on the same domains against the same target. The extended reach's
own end domains are excluded the way GIS 1 and 90 are in the base.
```

```text
The hindcast runs as one of two periods. Picking the start year resolves
every period-dependent forcing: run length, RSLR rate, storm series, the
BRIE island-offset file that sets the starting shoreline, the road setback
file, nourishment defaults, and which background-erosion preset applies.

Paths are relative to data/hatteras_init/ so this module stays independent
of where the repo is checked out; join them onto your own data base dir.

ROAD SETBACK: switched to the DUNE-START method on 2026-08-18.

The setback is metres landward of interior row 0 -- the row
roadway_manager.py:99 indexes against -- and the dune-start method is the
first one that measures it there. The legacy files measured against the
same-year digitised dune line instead, a different feature from a different
year, which put the road a median +40 m (1984) / +23 m (2004) landward of the
road actually rasterized onto the model grid, reaching 130 m.

TWO THINGS THAT MUST MOVE TOGETHER WITH THIS:
1. The topography must be the extraction these setbacks were measured
against -- they are metres landward of ITS interior row 0, so spending
them on another extraction measures from a row that does not exist.
This no longer has to be remembered: the runner, the sweep worker and
the road scripts all resolve the version through
scripts/site_layer/hat_topo_version.py, which reads VERSION out of the extractor.
Bump it there and everything moves together. (This comment used to
pin 2009_v3 by hand and went stale the day the setbacks moved to v4,
then again when it still said "Current: 2009_v5" after the tree went
period-first. It no longer names a version at all -- ask topo_dirs().)
There are now TWO extractions, one per period; see "topo_product".
2. The legacy files are still on disk under old_method_offset/ for the
method comparison. They are NOT interchangeable -- see
scripts/input_prep/4-mgmt-forcings/road_offset/README.md.

Both are 2 rows x 82 cols, GIS 9-90, so the swap itself is a drop-in.
```

```text
"topo_product" names the folder under
data/hatteras_init/1-barrier3d-domains/ that this period's Barrier3D domains
come from. It sits beside storm_file / island_offset_file / road_setback_file
because it is the same kind of thing: a per-period input the run ingests.

ADDED 2026-08-25. Before that the runner hardcoded ONE topography
("2009-dune-topo" / <version>) and BOTH periods read it, so a 1984 run and a
2004 run started from the same barrier. They no longer do:

1984  <- DEM 2009-2014-1996  (1996 ALACE overwriting measured ground
wherever ALACE has data; the landward limit
is the swath edge, near the dune toe. The
"ocean-side of the 1984 NC-12 line"
boundary this comment used to name was
dropped 2026-08-26 -- see that product's
README for the three-way measurement.)
2004  <- DEM 2009-2014       (the baseline gap-filled DEM)

The version WITHIN a product is still resolved, never pinned - see
scripts/site_layer/hat_topo_version.py.

NOT A LITERAL ANY MORE (2026-08-26). The value comes from YEAR_PRODUCT in
hat_topo_version.py, which is the same mapping every road script in
4-mgmt-forcings now resolves through. It was spelled out in four places and
omitted in three, and the three that omitted it - the placement figure, the
method diagnostic and the per-domain views - gave BOTH vintages 2004-start
interiors. Since 65 of 90 domains have a different interior shape between the
products, that is a different island, not a rounding difference. One mapping,
imported, so the runner and the forcing that feeds it cannot disagree.
WHICH BUILD OF THE ISLAND OFFSET A START READS (added 2026-09-15).

2-brie-offset/<year>/ used to hold one build. When the 1997 dune line was
re-digitised (duneline_1997_v2, local corrections) the 1996 start gained a
second build, and the two live side by side as 1996/v1/ and 1996/v2/ with a
CURRENT file naming the one every reader takes (2026-09-18: the 1997, 2009
and 2023 lines re-digitized; 1996/v3/ and 2010/v2/ were CURRENT, renumbered
1996/v1/ and 2010/v1/ on 2026-09-19 with the earlier builds under
<start>/superseded_20260919_pre-redigitized/) -- the same shape as
1-barrier3d-domains/<product>/dune-topo/. Resolved here, in one place, so
the runner cannot pin a path that a later re-digitisation silently leaves
stale. Order, mirroring hat_topo_version.topo_dirs():
1. HAT_OFFSET_VERSION_<year> in the environment (per-run selection that
does not mutate the shared default)
2. the CURRENT file in 2-brie-offset/<year>/
3. the only v* directory present, if exactly one
4. no v* directory at all: the flat layout, 2-brie-offset/<year>/<file>
(1984 and 2004 today; 2010 once it is built)
Several v* directories and no CURRENT is an error, not a guess.
```

```text
Every part of this path comes from hat_topo_version (2026-09-22). It was
built by hand here, which meant the source split -- <year>/v1/ becoming
<year>/duneline/v1/ -- would have left this resolving a path that no
longer exists, silently for every period at import.
```

```text
The version choice (env, CURRENT, the only v<n>; errors otherwise) lives
in hat_topo_version.offset_version since 2026-09-18, so the figure
scripts that used to repeat it read the same build the runner does.
```

```text
0.00391 m/yr fitted over 1984-2004 on the Duck gauge, stored to
0.001. The fits are in 3-env-forcings/2-rslr/fits/duck_rslr_rates.csv
(column config_m_yr is this rounding), written by
scripts/input_prep/3-env-forcings/2-rslr/duck_rslr_analysis.py.
```

```text
PAIRED WITH THE TOPOGRAPHY VERSION, AND NOTHING ENFORCES IT.
A setback is metres landward of interior row 0, so it belongs to the
extraction it was measured on. This file is the v2-era measurement;
the v1-era one it replaced is in
road_offset/archive/superseded_20260907/1984/.

Measured 2026-09-14: they differ at 27 of 82 domains, mean -12.4 m
and up to 205 m at GIS 35. Running v1 ARRAYS against these v2
setbacks moves every domain's rate by up to 0.005 m/yr -- small, but
it is a mismatch, not a result. Pinning HAT_TOPO_VERSION_1984_START
WITHOUT also pointing this at the matching file is the trap; it cost
a twelve-run comparison before it was noticed.
road_offset/dunestart_offset/measured/1984/RoadSetback_1984_dunestart.csv
-- MEASURED on the 1978 NC-12 line (ROAD_LINE_FOR_YEAR[1984]) against
row 0 of 1984-start.
```

```text
road_offset/dunestart_offset/measured/2004/RoadSetback_2004_dunestart.csv
-- MEASURED on the 2008 NC-12 line (ROAD_LINE_FOR_YEAR[2004]) against
row 0 of 2004-start.
```

```text
TWO MORE PERIODS, ADDED 2026-09-11. They OVERLAP the two above rather
than partitioning the record: four hindcast windows over one island, not
a timeline cut into quarters.

THE END YEAR IS A BOUNDARY, NOT A SIMULATED YEAR. run_cascade_simulation
spends start..start+run_years-1, so 1996 runs 1996-2009 and 2010 runs
2010-2023 -- the two new windows tile without overlapping each other, and
the 2010 survey is both the first one's target and the second one's start.

THEIR TOPOGRAPHY IS SHARED, NOT NEW (Hannah, 2026-09-11). 1996 reads
1984-start, whose DEM carries the 1996 ALACE graft -- the survey nearest
that start, and a better match for it than for 1984. 2010 reads
2004-start, built from the 2009 and 2014 lidar, which is the closest
vintage match of any period to its own start year.
```

```text
0.00402 m/yr fitted over 1996-2010 on the Duck gauge; the other three
periods are stored at this precision too. See rslr/fits/duck_rslr_rates.csv.
```

```text
DERIVED, NOT SURVEYED: built from the 1997 dune line, the nearest
island-wide survey (hat_topo_version.DUNE_LINE_FOR_YEAR[1996] ==
1997; the end-year target loader reads the same table). See
2-brie-offset/raw_offsets/PROVENANCE.md.
```

```text
DERIVED: the 1984 setbacks with the 1989 Pea Island relocation
applied, since that event precedes 1996 and the 1999 one does not.
No NC-12 line of 1996 vintage exists. See that folder's PROVENANCE.md.
road_offset/dunestart_offset/derived/1996/RoadSetback_1996_dunestart.csv
```

```text
0.00651 m/yr fitted over 2010-2024. It rounds up where 2004-2024
rounds down (0.00639), so the two differ by more in this table than
in the gauge record. See rslr/fits/duck_rslr_rates.csv.
```

```text
DERIVED, NOT SURVEYED: built 2026-09-15 from the 2009 dune line, no
2010 aerial imagery existing (DUNE_LINE_FOR_YEAR[2010] == 2009). So
the period starts from the island as surveyed a year EARLIER, the
mirror of the 1996 case. 2010/v1, CURRENT.
```

```text
A COPY of the 2004 file: same topography product, same road line, and
no relocation in the record between the two dates.
road_offset/dunestart_offset/derived/2010/RoadSetback_2010_dunestart.csv
```

```text
Sparse by design: a domain absent from a preset gets 0.0 m/yr. Sign follows
cascade/brie_coupler.py: (-) = erosion, (+) = accretion, in m/yr.

Three presets, one per hypothesis about where the alongshore sediment budget
is unresolved -- they are the source/sink axis of the run matrix:

zeroBE   no source/sink anywhere. Whatever the shoreline does is what
Barrier3D + BRIE produce unaided.
edgeBE   only the two end domains, which absorb the open-boundary artifact
at the ends of the modelled reach. Nothing is imposed on the
interior, so an interior misfit is the model's own.
calibBE  the full per-domain fit against the CoastSat LRR target rates.

Moved here from HAT_hindcast_1984_2024.py so the notebook and the
run script read the same numbers. data/hatteras_init/7-source-sink/4-export/
holds copies GENERATED from these by export_be_calibration.py; the older
partial copies that disagreed (the 2004 one truncated mid-dict) are under
7-source-sink/archive/. These are the ones that have been run.
```

```text
"zeroBE" preset -- empty rather than {1: 0, 90: 0}, because the sparse
contract already gives an absent domain 0.0 and an explicit zero reads like a
value that was solved for.
```

```text
GIS domains carrying the edge-only preset: the first and last REAL domains.
Not the padded buffers, which stay 0.0 in every preset. (1, 90) in the base
geometry; an extended geometry moves an end, and the value there comes from
HAT_BE_OVERRIDE while its solve is in progress (see below).
```

```text
"calibBE" preset -- the per-domain fit. HATTERAS_BE_RATES_EDGE is derived
from this below, so the two presets cannot disagree about the end domains.

FIT AGAINST THE LRR, NOT THE ENDPOINT DIFFERENCE. These rates are the
LOWESS-smoothed residual of an edgeBE base run against the CoastSat
target, and both sides of that residual are now the same estimator: the
target is a per-transect OLS slope (transect_lrr_full.csv, lrr_m_yr) and
the model side is read from the run's lrr_m_yr, its own OLS slope
through 21 annual states. Before 2026-08-22 the model side was
change_rate_m_yr -- (x[-1] - x[0]) / span -- so the residual carried the
gap between two estimators as if it were a sediment budget.

The refit moved 36 of 88 interior domains in 1984-2004 (mean 0.21, max
0.50 m/yr) and 65 of 88 in 2004-2024 (mean 0.39, max 1.60 m/yr). Period
2 moves further because that is the period with nourishment: a fill is
an instantaneous step in x_s, BRIE answers a step with a slowly-decaying
alongshore grid mode, and the 2022 Buxton and Avon fills are two years
from the end of the run -- so the endpoint difference was reading solver
ringing into the residual and calling it background erosion.

Regenerate with scripts/input_prep/7-source-sink/2-calibrate/
be_zone_residual_fit.py; GIS 1 and 90 are NOT taken from its
output, which writes them as 0.0 -- they stay the separately solved
buffer-cell values below.

THE END DOMAINS WERE RE-SOLVED AGAINST THE LRR ON 2026-08-23, for the
same reason the interior was: the values they replace were solved when
the model side was the endpoint difference, so they were fit to a
different estimator than the one now plotted. Previous values were
GIS 1 / 90 = -24.0 / +10.0 (1984) and +35.0 / +35.0 (2004).

HOW THEY WERE SOLVED. Two Newton steps against the base run below: a
first from the zeroBE/edgeBE secant, then a second on the local secant
through the two bracketing runs. The gain is SMALL -- d(LRR)/d(BE) is
0.092 to 0.123 across the four cases, i.e. only about a TENTH of an
imposed edge rate survives in that domain's own shoreline, the rest
being diffused alongshore by BRIE within a few domains. So each value
here is roughly
ten times the misfit it is there to close, and is a boundary-artifact
absorber, not a sediment budget: read as a real flux over the 22.25 m
shoreface, GIS 1 in 1984-2004 would be ~5e5 m^3/yr across one 500 m
domain. It is not capped for plausibility, because capping it would
just move the open-boundary artifact back into the figure.

THE GAIN IS NOT CONSTANT, which is why one step was not enough. It
STIFFENS with the imposed rate at GIS 1 in 1984-2004 (0.098 over
BE 0 to -24, but 0.123 over -24 to -46.4) and SOFTENS in 2004-2024
(0.097 to 0.092). At -42 m/yr for 20 years GIS 1 is being pushed
~840 m landward, far outside the range where a background rate acts
as a small perturbation, so a single global slope should not be
assumed if these are ever re-solved.

WHAT THEY WERE FIT TO. The value the rate-comparison figure DRAWS at
each end, so fit and figure cannot disagree:
GIS 1   raw per-domain transect mean -- LowessConfig.skip_southern_
domains is 10, so D1-D10 are drawn raw, not smoothed.
GIS 90  the LOWESS value of the primary window (7 domains since
2026-09-28, 10 before), which is what is drawn north of D10.
The two ends therefore use different estimators. That is deliberate:
it mirrors the splice the figure already makes.

BASE RUN. edgeBE / road_bdm / groin off, per period -- the same run
be_zone_residual_fit.py derives the interior residual from.

WHY GROIN OFF. Groin-off is a CHOICE, and it is the right one: the groin is
a structure whose
trapping is fitted against the same period-1 shoreline these BE rates are
fitted against, and letting both absorb the same misfit makes neither
identifiable. The BE fit goes first, groin off; the sweep then pins be1 at
the production edgeBE value and fits the groin on top.

THE GROIN FIT, RE-EXAMINED 2026-08-30. M = 60, f = 0.6 STANDS.

It was briefly marked void earlier that day. That was an over-correction and
is withdrawn. The reason given was true -- every period-1 sweep cell behind
the original fit was computed on the WRONG ISLAND, because
HAT_groin_sweep_worker.py resolved topography without naming a product and
DEFAULT_PRODUCT ("2004-start") answered, so a 1984 sweep built on the 2004
barrier. That bug is real and was fixed 2026-08-30; the drift guard now
reproduces its reference to exactly 0 m/yr. But re-scoring on the CORRECTED
island, at the re-calibrated be1 = -42.6, barely moved the answer:

D4-D8 demeaned profile RMSE, period 1, be1 = -42.6
no groin           15.20 m
M=50  f=0.8        11.44   <- best cell
M=40  f=1.0        11.44
M=50  f=1.0        11.52
M=60  f=0.8        11.58
M=60  f=0.6        11.58   <- the production value

11.58 against a best of 11.44 is 0.14 m across a 3.76 m improvement. The
production value is statistically indistinguishable from optimal, and the
topography bug moved it by less than the ridge it sits on.

BOTH M AND f ARE FITTED -- ON THE WINDOW THAT SPANS THE STRUCTURE'S LIFE.
Corrected 2026-08-30 after re-running the 1967-2018 rig. Earlier revisions of
this note said f was "not determined" or "prescribed rather than fitted".
That is true of the HINDCAST windows and false of the full-life window:

1967-2018 rig, 51 years, spans install (1969/70), the 1996 repair, the
2003 storm damage and 14 years of decline. RMSE at M = 60:

f     0.1     0.3     0.4     0.5     0.6     0.7     0.9
33.17   27.99   25.84   24.21   23.78   25.06   39.01
^ best, bracketed both sides

f = 0.6 is a clean INTERIOR minimum with steep curvature either side.

M IS NOT. Corrected 2026-08-30 -- an earlier revision of this note said
"So is M: 53.9 (M=20), 33.3, 27.9, 23.8 (M=60). Neither railed." The rig's
own CSV refuses that. M improves MONOTONICALLY to the last value that runs:

M      20      40      60      70          80          >=100
53.9    33.3    23.8    320-378     556-766     (blank)

The jump at M = 70 is a 13x discontinuity, not a fit degradation -- it is
the instability dipole_fit_notes.md already records ("M >= 70 went unstable and
M >= 100 drowned the barrier on the 41-domain rig"). M = 60 is the LARGEST
M THE RIG CAN HOLD, not the M where the fit stops improving. The two
documents had contradicted each other; dipole_fit_notes.md was right.

SO THE RIG CORROBORATES f, AND IS ONLY CONSISTENT WITH M. Do not write
"two independent routes agree on both parameters." Write: the rig brackets
f = 0.6 on both sides; it rails against a stability wall in M, so it cannot
distinguish M = 60 from any larger value. M = 60 stands on the PRODUCTION
fits (D4-D8 and the independent D3-D9), not on the rig.

"INDEPENDENT" IS ALSO GENEROUS. Since the 2026-08-30 repoint, both routes
use the same topography product (1984-start/v1), the same wave climate, the
same BRIE physics and the same wet/dry shoreline family. They are
independent WINDOWS, not independent EVIDENCE.

AND THE RIG RUNS 1967 OFF A 1984 ISLAND. RIG_TOPO_PRODUCT = "1984-start"
(in HAT_groin_hindcast_1967_2017.py) -- a deliberate 17-year anachronism
in the initial condition, accepted because the target is a shoreline
OFFSET rather than an elevation. It belongs in any methods description of
the rig. The rig also uses GROIN_INSTALL_YEAR = 1970 against the plan's
1969; one year, inside the build phase, but they are not the same number.

WHY THE HINDCAST WINDOWS CANNOT SEE f. They begin 15 years after installation
and end before or just after the collapse, so they contain almost none of the
1996-2003 deterioration ramp. Period-1 cumulative trapping is M(15.5 + 4.5f),
which f moves by only 29% across its whole range -- so period 1 reads high f
as simply "more trapping" and rails at 1.0. Period 2 is 20*M*f and prefers
f = 0. Neither is fitting the deterioration; they are fitting its absence.

f IS NOT A FREE KNOB EITHER WAY -- it encodes a maintenance record (installed
1969, last repaired 1995 (CSE 2013, p. 20), storm damage 2003, fillet peaks 2004). Making the
module a STATIC trapping rate was considered on 2026-08-30 and rejected on
measurement: at M = 60 setting f = 1 degrades period 2 from 14.90 to 17.87 m
(+20%), and at M = 95 from 14.90 to 20.38 (+37%). It also moves the modelled
D5-D6 differential the WRONG WAY -- observed is -2.47 m/yr, and the model goes
-0.55 at f = 0 to -0.18 at f = 1. Deterioration is doing real work.

M AND f ARE SET FROM DIFFERENT EVIDENCE, AND THAT IS DELIBERATE.
The authority for these values is hard-structures/groin/3-hindcast/1-dipole-1967-2017/dipole_fit_notes.md
(2026-08-24); this note summarises it and must not diverge from it. The
figures testing it, and what each one showed, are described in
scripts/hatteras_ms/groin-sweep/CALIBRATION_FIGURES.md -- the PNGs themselves
land in output/calibration/groin/figures/, which .gitignore does not track.

M from PERIOD 1.  f from the 1967 rig and from period 2.

WHY THEY CANNOT BOTH COME FROM PERIOD 1. Period-1 cumulative trapping is
M(15.5 + 4.5f), because period 1 mostly PRECEDES the 1996-2003 deterioration
ramp -- so f moves it by only 29% across its entire range. Period 2 is
20*M*f, where f = 0 gives zero: total leverage. Period 1 therefore fixes M
and barely sees f; period 2 fixes f and cannot see M at all.

THE RIDGE IS IN M(15.5 + 4.5f), NOT IN M*f. At f = 1.0 the best M is 50; at
f = 0.6 it is 70. So a poor score at M = 50, f = 0.6 is NOT evidence against
f = 0.6 -- it means M was set too low for that f.

A CORRECTION TO A CORRECTION, 2026-08-30. This note briefly claimed "only the
product M*f is identified" (taken from HAT_period1_top_n_figure.py's caption
without checking), and then, having tested that and found corr(RMSE, M*f) =
-0.07, claimed instead that M and f are separately and weakly constrained.
BOTH were wrong. The product test was right that M*f is not the invariant and
wrong about what is: the invariant is period-1 cumulative trapping. Quote M
and f as a pair, and cite dipole_fit_notes.md for why.

f = 0 IN THE PERIOD-2 SWEEPS IS THE RIGHT ANSWER, NOT A RAIL ARTEFACT. The
observations show the fillet declining after 2004, i.e. trapping ceased.
dipole_fit_notes.md records that considerable time was lost re-defining targets to
"fix" a result that was correct. Do not re-litigate it.

THE GROIN DOES REAL WORK BUT DOES NOT REPRODUCE THE SHAPE. 15.20 -> 11.44 m
is a 25% reduction, and it comes from matching the overall D4-D8 slope. The
observed profile has structure the model does not produce: a peak at D6
(+14 m observed, ~0 modelled), a dip at D7, a second peak at D8. Read the
residual as the split between what the groin explains and what the
source/sink calibration absorbs -- not as a successful shape fit.

THREE TARGETS THAT DO NOT WORK, so they are not retried:

1. THE FILLET (D5-D6 scalar). No admissible M can match it on this grid --
stated at HAT_groin_timeseries_check.py:29. Fitting it anyway on
2026-08-30 gave M = 95 at be1 = -42.6, with f railed at the grid bound
and the best M swinging 95 -> 160 between adjacent be1 values. Fitting
an unmatchable target is what produced the rail, not a real optimum.
2. THE FULL-PERIOD D1-D12 PROFILE. Ranks M = 0 best, monotonically. Not
because there is no groin: D2-D4 is Cape Point accretion the
parameterisation does not represent, and D6-D7 is an erosion trough
peaking one domain NORTH of the structure, which a groin actively
worsens by pushing D6 seaward. Together they swamp a ~17 m groin
signal. See HAT_fullperiod_windows.py.
3. NARROWING THAT WINDOW TO D5-D7 does not rescue it -- D6-D7 is inside
the narrow window too. D4-D8 DEMEANED is the window that works, because
demeaning removes the level offset the source/sink term owns and D1 is
excluded (the cape's 81-104 m change is ~5x the groin's signal).

PERIOD 2 IS NOT FITTED, BUT THE GROIN IS STILL ON FOR IT. GroinCallback
carries an ABSOLUTE CALENDAR timeline -- install 1969, deterioration onset
1996, linear ramp to 2003, then hold at M*f -- so no period-specific
configuration exists or is wanted. Running it in period 2 is right for
consistency of the structure's timeline, not because it explains that
period's shoreline.

What period 2 records is a RELEASE the module cannot produce: -76 m, of which
dipole_fit_notes.md attributes ~85% to the UPDRIFT side eroding once the structure
failed, not to impounded sand draining downdrift. Trapping is bounded at >= 0,
so the groin can stop adding sand but cannot drain the fillet. That -76 m is
carried by the source/sink calibration together with the Cape Point dynamics
the dipole does not represent.

WHAT THE GROIN EXPLAINS:  period 1  +17.2 m of an observed +52 m (33%).
period 2  ~0 of an observed -76 m.

So the 2026-08-30 joint fit railing at M = 160 was not a fitting failure to
be repaired -- fitting period 2 is the wrong thing to attempt.

BUT IT WAS STILL SITTING IN THE FILE THE PIPELINE READS. Found 2026-08-30:
output/calibration/groin/joint_fit.json held the RANKING'S answer (edgeBE M = 160
f = 0.8, zeroBE M = 140 f = 1.0), and HAT_run_all.py stage 6 passes whatever
that file holds to every groin run in the matrix. A stage-6 run would have
used M = 160. It prints "RAILED on M" while doing so, so the machinery knew.

The file is now PINNED by hand to M = 60, f = 0.6 for both presets, each
entry carrying its own superseded_ranking block, and the ranking is archived
at output/archive/2026-08-30_groin-railed-ranking/joint_fit_RAILED_ranking_20260830.json.
RE-RUNNING STAGE 5 OVERWRITES IT -- re-pin afterwards.

Two more places the pair is written, both corrected the same day:
scripts/hatteras_ms/hat_run.yaml carried M = 50 / f = 0.9 placeholders (a
2026-08-26 note left them stale on purpose pending the re-fit, which has now
happened and did not move the answer), and HAT_hindcast_config.py carried the
same values as fallback defaults. A MATRIX run reads joint_fit.json; a SINGLE
run reads the yaml. Check both before quoting a groin run's parameters.

AFFORDABILITY IS A SOFT BOUND, NOT A CEILING. M = 60 intercepts ~719,000
m3/yr against a 5-7e5 m3/yr littoral drift -- marginally above a LITERATURE
RANGE, which dipole_fit_notes.md is explicit is "not a hard limit". Earlier text
here treated it as one; it is a reason to prefer 60 over 70 (838k, ~1.3x the
drift), not a physical prohibition.

M IS NOT A SEDIMENT FLUX. It is an effective, grid-specific, FIELD-AGGREGATE
rate: the real fillet is ~190 m wide against a 500 m domain, and the four
Buxton groins span northings entirely inside D6, so one dipole expresses the
whole field. Do not divide M by four for a per-structure value, and do not
read it as a flux.

NOW MEASURED, 2026-08-30. groin_diagnostics.csv has logged the cumulative
displacement every year of every groin run since the module was written, and
nothing plotted it until HAT_groin_sediment_budget_figure.py. Over the rig's
50 years at M = 60, f = 0.6:

applied   2,400 m cumulative one-sided displacement  (+/- 28.7e6 m3)
realised     69 m fillet at the end   (peak 129 m)
retained    ~3%

BRIE's alongshore diffusion removes the rest. M is therefore the rate needed
to SUSTAIN a fillet against diffusion, not the rate at which sand is
impounded -- which is the quantitative version of "not a flux".

THIS CHANGES HOW THE AFFORDABILITY NUMBER READS. 719,000 m3/yr at M = 60 is a
GROSS RESTORING RATE set against a NET transport budget; they are not like
for like. The comparison is still a fair reason to prefer M = 60 over M = 95,
and it is a deliberate documented diagnostic in the run reports -- but
"marginally above the drift band" must NOT be read as "impounds more sand
than the coast carries." Figure: output/calibration/groin/figures/sediment_budget.png

HOW FAR THE VALUE IS ACTUALLY CONSTRAINED -- ROBUSTNESS, 2026-08-30.
Everything below was scored from the existing period-1 cells at
be1 = -42.6; no new runs. Read it before quoting M = 60 as "the fitted value".

WINDOW SENSITIVITY. D4-D8 is a CHOICE, and the answer moves with it:

window    best M    f     RMSE   no-groin   gain
D4-D8         60   1.0    10.09     12.19    2.10   <- production
D3-D9         60   1.0    13.93     15.98    2.05
D4-D7         70   1.0    10.54     13.63    3.09
D3-D8         95   0.8    12.91     16.80    3.89
D4-D9         40   1.0    12.23     12.65    0.41
D5-D7          0   0.0     6.18      6.18    0.00
D5-D8          0   0.0     6.18      6.18    0.00

THE SUPPORTING HALF: D3-D9, the one independent window that spans the
structure symmetrically, also returns M = 60 with the same gain (2.05 vs
2.10). M = 60 is therefore not an artefact of picking D4-D8 specifically.

THE QUALIFYING HALF: across defensible windows M ranges 40 to 95. The window
does real work in setting the answer. D4-D8 is justified by argument -- keep
D1's cape signal out, centre on the structure -- not forced by the data.

AND THE SHARPEST RESULT: D5-D7 and D5-D8, the windows TIGHTEST on the groin,
return M = 0 with a gain of exactly 0.00. The groin changes nothing where it
stands. The 2.10 m gain at D4-D8 comes from D4 and D8, the window's OUTER
edges. This is consistent with HAT_fullperiod_windows.py's finding that
D6-D7 carries an erosion trough peaking one domain north of the structure
that the module cannot produce -- but it should be said plainly: the domains
closest to the groin are the ones the groin explains least.

CONDITIONAL ON be1, NOT INDEPENDENT OF IT. The D4-D8 optimum trades off
monotonically with the edge rate: M = 50 at be1 = -46, 60 at -42.6, 80 at
-34, 125 at -10. Two things make this acceptable where the fillet's version
was not -- the surface is smooth and monotonic (the fillet swung 95 -> 160
non-monotonically between adjacent be1 values), and the GLOBAL minimum over
the whole (be1, M, f) space sits at be1 = -42.6, the independently calibrated
value. The groin fit prefers the same edge rate the source/sink fit arrived
at, which is a real if modest cross-check.

NOT PRESET-INDEPENDENT. Under zeroBE the D4-D8 fit improves monotonically all
the way to the grid edge (no-groin 21.9 -> M=140: 12.4, M=160: 12.4). With no
edge forcing the groin simply absorbs the missing source/sink term. A zeroBE
groin fit is meaningless, and M = 60 exists only in the presence of edgeBE.

PERIOD 2 CARRIES NO INFORMATION AT ALL. Scored on D4-D8, every M from 0 to
160 gives RMSE 14.90 (edgeBE) or 22.42 (zeroBE) -- identical to four
significant figures, not merely "weakly preferring zero". Confirms the plan:
by period 2 the deterioration schedule has trapping at ~0, so running the
groin there is a timeline-consistency choice and not a fitted one.

WHAT THIS ADDS UP TO. M = 60, f = 0.6 is the best-supported pair available,
and the support is asymmetric between the two parameters:

M comes from the PRODUCTION fits -- D4-D8 (11.58 vs 15.20 no-groin) and
the independent, structure-symmetric D3-D9, which returns M = 60 with
the same gain (2.05 vs 2.10). The 1967 rig does NOT confirm M: it rails
against a stability wall at M = 70 (see above).
f comes from the 1967 RIG, where 0.6 is bracketed on both sides, and from
period 2, where the observed fillet declines. The hindcast windows
cannot see f at all.

CONDITIONS THAT TRAVEL WITH THE PAIR. Conditional on be1 = -42.6 (the D4-D8
optimum moves 50 -> 125 as be1 goes -46 -> -10); void without edgeBE; carries
no information in period 2; and M varies 40-95 across defensible windows,
with the windows TIGHTEST on the structure (D5-D7, D5-D8) returning M = 0 at
gain 0.00. Quote it as a calibrated parameter pair with those conditions
attached -- not as a measured property of the structure.

DECIDED 2026-08-30: LOCKED. Single value, conditions stated, no ensemble over
M and no further runs. The 1967 window needs no new work -- the rig IS that
window (see below) and was re-run on the corrected topography that day.

BE CONVERGENCE, 2026-08-31: THREE PASSES, AND WHERE IT STOPS.
The field is converged on everything the protocol can converge, and the
remainder is the wavelength it deliberately does not chase.

pass    P1 domains>=0.5   max    mean      P2 domains>=0.5   max    mean
1          27          -1.70   0.35           19          -1.70   0.27
2           8          +1.00   0.10            7          -1.00   0.11
3           2          +0.90   0.05            3          -0.90   0.09

Passes 1->2 closed 41% of the max, matching the documented "one pass closes
42% (P1) and 57% (P2)". Pass 2->3 closed only 10%. That is not a stall to be
pushed through -- it is g, and the drop is diagnostic.

WHAT CLOSED, AND WHAT DID NOT. D83-D88 and D72-D74 -- contiguous, same-signed
-- dominated passes 1 and 2 and are now EXACTLY 0.0 in both periods. What
remains is D8 +0.9 against D10 -0.9 across a zero at D9, plus isolated single
domains at D22 and D32. be_zone_residual_fit.py:150 predicts precisely
this: "a contiguous same-signed block of corrections passes at g ~ 0.8-1.2,
while a pattern that alternates sign at the grid scale is damped to g ~ 0.1."
The measured 10% IS that g.

WHY MORE PASSES ARE THE WRONG ANSWER. At 10% per pass, 0.90 -> 0.50 needs six
more passes, each a full calibBE matrix rebuild -- ~12 hours to close a
feature the design already refused to close. The same comment records that
amplifying by 1/g was considered and REJECTED because "narrow features would
be amplified ~10x into rates that are indefensible read as sediment fluxes".
Iterating to the same place by brute force does not make them defensible. A
+-0.9 m/yr alternation between adjacent 500 m domains is more plausibly noise
in the model-observation comparison at that wavelength than a source/sink
signal worth encoding.

THE TOLERANCE TO QUOTE. The field is converged to 0.0 m/yr on contiguous
alongshore structure and to ~0.9 m/yr on isolated grid-scale features, with a
mean absolute residual of 0.05 (P1) and 0.09 (P2) m/yr. For scale, the
interior RMSE these runs are scored on is 1.1-2.8 m/yr, so the remainder sits
well inside the noise the model is fitted against. It does NOT meet the
script's own SIGNIFICANCE_THRESHOLD of 0.5 per zone, and should not be
described as if it did.

STATE ON DISK: config holds the PASS-2 field and all 22 calibBE runs are
built on it -- consistent. Pass 3's residual was measured and deliberately
NOT applied, because applying without a rebuild is what left the field two
days stale in the first place.

STALENESS AUDIT, 2026-08-31. NOTHING IS STALE, AND ONE WORRY WAS WRONG.
Run after the topography work, to check what the updated island invalidated.

55 production runs   topography CLEAN. Every 1984 run is on 1984-start/v1
and every 2004 run on 2004-start/v1 -- what hat_topo_version.py
resolves today. Zero mismatches.
55 production runs   BE field CLEAN. Every be_values_digest matches the
current config, including the 22 calibBE runs rebuilt that day.
period-1 sweep       CLEAN, by re-running cells and diffing.

THE WORRY THAT WAS WRONG. The period-1 edgeBE cells split across two dates --
427 at be1 = -10..-46 on 08-29 20:31, and 61 at be1 = -42.6 on 08-30 08:21 --
with the sweep worker's topography fix committed between them (562c75c,
08-30 13:01). Since result.json records NO topography, that looked like the
be1-sensitivity table above might be comparing two different islands, which
would have manufactured its own "global minimum sits at -42.6" result.

It does not. Re-running one cell from EACH batch reproduces the stored
result to 0.00e+00 -- bit for bit, across differential_m_yr, rmse_window,
bias_window and rates D4-D8. Both batches are on today's topography; the two
dates were two sweep sessions. The be1-sensitivity table stands as written.

PERIOD 2 IS NOT BIT-REPRODUCIBLE, AND PERIOD 1 IS. Re-running a period-2
cell twice gives two answers differing by 1.0e-04 -- MORE than either differs
from the stored value, which sits inside their spread. So the stored cell is
fine, but "re-run and compare" is only an exact staleness test for PERIOD 1.
For period 2 the tolerance is ~1e-4.

Not an unseeded RNG: roadway_manager seeds at 1973 and nothing else in
cascade/ draws randomly. Floating-point summation order is the likely
source, and period 2 accumulates more of it through the nourishment events.
1e-4 is five orders below anything that moves a ranking -- the RMSE
differences this calibration turns on are 0.1-4 m -- so it changes the TEST,
not any result.

THE fullperiod_1984_2024 SWEEP IS ALSO CLEAN -- checked 2026-08-31, and the
note claiming otherwise is withdrawn. CALIBRATION_FIGURES.md said panel (b)
of fig_three_targets ran on the PRE-FIX topography. It did when that was
written, and then the sweep was re-run 08-30 18:20 -- five hours AFTER the
worker fix at 13:01 -- and the figure rebuilt at 23:52 from those cells. The
note outlived the problem. Re-running cell M60_f0.50 reproduces to 8.3e-05,
the period-2 noise floor the 40-year window inherits. All 43 cells present.

STILL UNVERIFIED: only the 1967 rig, which files no run_index row at all.
Its two current runs were rebuilt 2026-08-31 at M = 60 / f = 0.6 and it now
lives in output/calibration/groin_rig/, away from production.

THE 1967 WINDOW HAS ALREADY BEEN RUN -- IT IS THE 41-DOMAIN RIG.
dipole_fit_notes.md recommends "fit on the 1967 window; apply in the hindcast",
because the hindcast windows begin 15 years after installation and record the
fillet's decay rather than its creation. That was checked on 2026-08-30 and
the answer is that the window exists already:

Change_from_wetdry_1967_D2_D12.csv covers D2-D12 -- ELEVEN real domains.
11 real + 15 buffer + 15 buffer = 41, the rig's exact domain count.

The rig was sized to the extent of the 1967 observations. Its sweep is at
hard-structures/groin/3-hindcast/1-dipole-1967-2017/results/sensitivity_sweep/.

RE-RUN 2026-08-30 ON 1984-start/v1; f MOVED ONTO THE PRODUCTION VALUE.
The rig had been resolving topo_dirs() with no product -- the same omission
that put the production sweep on the 2004 island -- plus two other stale
paths (domain_N_topography_2009.npy, and a "2009-buffer" directory that the
2026-08-25 restructure renamed). Repointed at 1984-start, the product nearer
the 1967 start in time and the one the production period-1 fit reads:

2026-08-24        2026-08-30 (1984-start/v1)
best cell          M=60, f=0.5       M=60, f=0.6
RMSE                    27.24              23.78

f moved ONTO the production value and the fit improved 13%. Read this as the
rig now supporting f, the parameter it can actually resolve -- NOT as
agreement on both. Its M = 60 sits against the stability wall documented
above, and would very likely rise if that wall were lifted.

ITS CACHE ALSO RESUMES ON (M, fraction, stage) WITH NO RECORD OF THE
TOPOGRAPHY, so the first re-run silently skipped all 67 stale cells and
reported the 2026-08-24 answer as fresh. Caught because the RMSE matched to
six decimals. Archive or delete HAT_groin_sweep_results.csv before any re-run
whose inputs have changed. Fourth instance of this bug class in this repo,
after the driver manifest key and the sweep worker's two call sites.

THE NUMBER STILL DOES NOT TRANSFER, the AGREEMENT does. The rig is a confined
41-domain array with its own structure parameters (install 1970, +25 yr onset
against the plan's 1969 / +27) and 1971/1973 nourishments the hindcast does
not carry. M is grid-specific. Its M >= 70 instability reappeared exactly as
documented -- RMSE 320-378 at M=70, and 11 of 39 attempts failed above it --
which is a rig artefact and not a production bound.

THE BLOCKER IS THE SHORELINE, NOT THE DEM. An earlier version of this note
said the 1967 window was impossible for want of a 1967 DEM. That was wrong:
the fillet is a SHORELINE quantity, so the 1984 DEM would serve perfectly
well for interior elevation. What does not exist is a 1967 SHORELINE for the
other 79 real domains -- the wet/dry table stops at D12.

WHY IT CANNOT BE LIFTED TO PRODUCTION GEOMETRY. Reconstructing 1967 offsets
for D2-D12 and holding D13-D90 at their 1984 values plants a discontinuity at
D12/D13. BRIE's diffusion length is ~3.2 km over 20 years, about six domains;
D12 to D8 is four. The artefact reaches the fit window before the build phase
finishes, contaminating the signal the exercise exists to measure. The only
alternative is inventing a 1967 shoreline for 79 domains.

M IS GRID-SPECIFIC, so the rig's value cannot be spent directly on the
120-domain grid -- a confined array preserves dipole amplitude that an open
one diffuses away. But AGREEMENT between the two is still evidence, and the
two routes agree: production-geometry period 1 on D4-D8, and a confined 1967
rig covering the build phase, both land on M = 60.

DECISION 2026-08-30: keep M = 60, f = 0.6. Revisit only if a pre-1984
shoreline record covering more than D2-D12 turns up outside this repo.

THE STABILITY BOUNDS DO NOT TRANSFER. "M >= 70 unstable, M >= 100 drowns" was
measured on the 41-domain rig. On the production 120-domain grid every cell
through M = 160 ran clean. Do not quote that ceiling for production.

The groin sits at GIS 5/6 and does not reach GIS 90 at all (on the pre-refit
runs D90 was identical to the 3 dp reported, groin-on vs groin-off, in
both periods); it moves GIS 1 by about 0.4 m/yr, so GIS 1 is worth
re-checking against the fitted groin.

STALE AS OF 2026-08-26 FOR PERIOD 1. Every rate below was derived from a
base run on the pre-restructure shared topography (now 2004-start/v1). The
1984-start product is a different surface, so the 1984 rates must be refit
before any calibBE or edgeBE run on it. Period 2 is unaffected -- its base
run read 2004-start, which has not changed.

D1 <-> D90 CROSS-TALK IS NEGLIGIBLE, so solving the two ends
independently is safe: with 15 buffer domains a side they are 30
domains -- 15 km -- apart around the ring, against a BRIE diffusion
length of sqrt(D*t) ~ 3.2 km over a 20 year period.
```

```text
Every value below was re-solved on 2026-08-28. The previous solution is in
git history and in output/archive/2026-08-28_full-tree/.

WHY. The 2026-08-23 values were derived on a base run against the
pre-restructure shared topography. `1984-start` is now a different surface
(re-extracted 2026-08-27 from the same DEM against a new pick set) and its
road setbacks were re-measured the same morning, moving 15 of 83 road-bearing
domains by up to 25 m. Period 2 was refit alongside it even though its inputs
had not changed, because the stage-5 groin joint fit intersects BOTH periods'
surfaces — mixed calibration vintages would make M and f mean two things.
(Period 2's zeroBE base run reproduced the archived one to 2e-4 m/yr, so its
inputs were verified unchanged rather than assumed.)

METHOD. Documented in be_zone_residual_fit.py: solve the two locked
ends by Newton steps on a secant, and the interior by iterated additive
passes (`be_apply_fit_to_config.py --add`) so each pass closes fraction g of whatever
misfit remains and no estimate of g is ever needed.

pass 0   interior from the edgeBE base runs           replace
pass 1-3 interior from the calibBE base runs          --add
GIS 90   re-solved after the interior settled          Newton, +3.0 probe
pass 4   final interior pass at the final edge values --add

INTERIOR RMSE OF THE BASE RUN (road_bdm, groin off), m/yr:

1984-2004   2004-2024
zeroBE               1.422       2.124
edgeBE               1.219       1.794
calibBE pass 0       0.721       0.763
calibBE pass 1       0.547       0.615
calibBE pass 2       0.527       0.583
+ GIS 90 re-solve    0.556       0.603     <- edges gained, interior gave back
calibBE final        0.517       0.563

Final mean interior bias: +0.058 (P1), +0.120 (P2).

STOPPING. `convergence_history.json` records the operative rule as "a pass
buys less than 5% of the standing RMSE". Interior passes hit that at pass 2->3
(3.8% P1, 5.3% P2). Note this is NOT the rule the LOWESS script's header
states ("no zone clears SIGNIFICANCE_THRESHOLD"); that one is unreachable
here, because a handful of domains carry residuals of 2-3 m/yr that are
narrower than the LOWESS window generating the correction and no smooth
alongshore field can close them. The 5% rule is the one that was met. The
residual it leaves — mean |smoothed residual| ~0.45 m/yr, with 19-28 domains
still above the 0.5 m/yr significance threshold — IS the tolerance this field
is converged to, and belongs in any methods paragraph built on it.

THE END-DOMAIN GAIN, MEASURED NOT ASSUMED:

GIS  1  edgeBE   0.109 (P1)  0.096 (P2)
GIS 90  edgeBE   0.104 (P1)  0.099 (P2)
GIS 90  calibBE  0.103 (P1)  0.094 (P2)

The last row settles a question this file used to leave open. The note below
says an isolated forced cell runs at g ~ 0.1, a forced ZONE at g ~ 1, and
that "GIS 90 sits on that join". Measured against a calibBE interior with a
coherent erosive zone at GIS 84-89 pressed against it, GIS 90's gain is
0.103 / 0.094 — indistinguishable from its own edgeBE gain, where every
neighbour is unforced. So GIS 90 does NOT sit on the join in practice: it
behaves as an isolated cell under both presets. The adjacent zone changes
WHERE its solution sits (calibBE +32.8 vs edgeBE +13.0 in P1) but not how
hard the cell is to move.

THE EDGES AND THE INTERIOR ARE WEAKLY COUPLED, and the loop contracts. GIS 90
drifted from converged (-0.06) to +0.62 over three interior passes as the
field at 84-89 grew, and re-solving it cost the interior 0.03 RMSE while
buying 0.61 at the edge — about 20:1, so it converges rather than oscillates.
The final interior pass was run AFTER the edge re-solve specifically so the
interior is fit against the final edge values and the table is self-consistent
as a set.

WHAT WAS NOT RE-SOLVED. edgeBE's GIS 90 (HATTERAS_BE_EDGE_D90) measured
converged throughout at -0.06 / -0.01 against the LOWESS target and was left
alone. GIS 1 in period 2 sits at -0.145 and was shrinking monotonically
(-0.231, -0.193, -0.145) but never cleared the threshold; it is a known small
residual, not an oversight.
```

```text
The 2004 calibrated preset falls back to the 1984 fit if it was never
solved separately. It has been, so this stays False -- but the flag is
what tells you whether a "calibrated 2004" run is really calibrated.
```

```text
GIS 90 IS THE ONE VALUE THE TWO PRESETS DO NOT SHARE.

GIS 1 is still sliced out of the calibrated fit below, so re-solving it
updates both presets at once -- that slicing is what stopped the old "base"
preset drifting into commented-out copies of these numbers, and it still
holds everywhere it can. GIS 1 can be shared because it lands within
0.22 m/yr of target in BOTH presets: nothing is forced near it (the
calibrated fit is 0.0 at GIS 2-4 in 1984-2004 and at GIS 2-7 in 2004-2024),
so its neighbourhood is the same either way.

GIS 90 CANNOT. Measured on the 2026-08-23 matrix, one shared value put
calibBE 1.136 m/yr (1984-2004) and 0.498 m/yr (2004-2024) off target there,
against 0.003 and 0.013 for edgeBE at the same number. The asymmetry is
structural, not a fitting error:

an ISOLATED forced cell is diffused away by its unforced neighbours, so
its gain d(LRR)/d(BE) is only about 0.1 and it needs roughly ten times
the misfit it closes. A forced ZONE moves together and is not diffused
away, so interior corrections run at a gain near 1.

GIS 90 sits on that join. In edgeBE it is isolated and the 0.1 gain holds.
In calibBE the fit puts a coherent erosive zone at GIS 84-89 (down to
-3.0 m/yr in 1984-2004) hard against it, and that zone drags GIS 90 down
by more than a metre a year -- an effect the per-domain LOWESS residual
never sees, because it fits each domain against a base run in which those
neighbours were unforced. So the two presets genuinely need different
numbers at GIS 90, and forcing one on them would mean one of them is
always wrong there.
```

```text
GIS 90 under edgeBE, solved on the edgeBE road_bdm base run. The calibBE
value for the same domain lives in HATTERAS_BE_RATES_CALIBRATED.
```

```text
PERIODS WITH NO CALIBRATED FIT TO SLICE GIS 1 OUT OF (2026-09-11).

The comprehension above takes edgeBE's GIS 1 from the calibrated preset, so
re-solving one moves both and they cannot drift apart. That works only where
a calibrated fit EXISTS. The periods added 2026-09-11 have none: they are
being solved edge-first, by the decision to do zero and edge before the
interior. Their two end values are therefore carried explicitly and merged
in below.

WHEN A CALIBRATED FIT ARRIVES for one of these, delete its entry here in the
same commit that adds the calibrated one. Leaving both would give GIS 1 two
homes, and the merge below would quietly win.
```

```text
(GIS 1, GIS 90), m/yr.

CURRENT VALUES: THE ADOPTED MODEL, 2026-09-28 evening (Hannah: "include
the overwash fixes, keep option A, go ahead"). Barrier3D hatteras/adopted
(the three overwash fixes + per-cell dune ceilings), storms v3_trim24,
option A waves, LOWESS-7 target, full management, edgeBE, no groin. The
LOWESS-7 values below were the seed:

1996-2010   GIS 1 +4.3509   GIS 90 +19.0935   residuals -0.006 / +0.007
2010-2024   GIS 1 +8.0      GIS 90 +21.2582   residuals +0.002 / -0.008

DUNE-CAP FIX, 2026-09-28 late (Hannah: "go with option 1, clip only the
bulldozed sand"): beach_dune_manager's 4 m cap now limits only the sand it
adds. 1996 was re-checked and still holds (-0.006 / +0.007). 2010 GIS 90
moved (+0.185 at +22.4937) and was re-solved by the secant, 4 steps from
+21.57 (experiments/end-domain-boundaries/2026-09-28-ends-resolved-dunecap/).
Before the fix it was +22.4937.

1996 by the safeguarded secant (4 steps). 2010 GIS 90 by the secant; 2010
GIS 1 by DIRECT PROBES: above ~+10 m/yr its response goes flat and even
reverses (+13.5 to +19 all leave ~+0.8 m/yr), which stalled the secant, so
it was mapped at -6, -2, +2, +6, +8, +8.8, +9.5, +10 and +8.0 closes it.
Record and every probe: output/raw_runs/experiments/end-domain-boundaries/
2026-09-28-ends-resolved-adopted/.

SPLIT12 STORMS, 2026-09-29 (Hannah: "use split12 as the storm series going
forward", then "re-run the matrix and re-solve the ends"). The storms are
now v3_split12_trim24 (grouped events split at >= 12 h below the berm, so
Fran 1996 and Jose 2017 are back). Seeded at the ends above, the secant
moved only GIS 1, by +0.04 in each window; GIS 90 held:

1996-2010   GIS 1 +4.3509 -> +4.3888   residuals +0.005 / -0.008   2 steps
2010-2024   GIS 1 +8.0    -> +8.0405   residuals +0.003 / -0.018   2 steps

(experiments/end-domain-boundaries/2026-09-29-ends-resolved-split12/)

--- LOWESS-7 option A values on the pre-adoption model, SUPERSEDED 2026-09-28 ---
OPTION A ON THE LOWESS-7 TARGET, 2026-09-28 (Hannah: "switch
the runner to 7 and re-solve the ends"). The runner's CoastSat target went
from LOWESS-10 to LOWESS-7 the same day; GIS 1 is graded against the raw
domain mean, so only GIS 90's target moved (+0.125 m/yr in 1996, -0.066 in
2010). Same waves, scenario and run as below; seeded at the LOWESS-10
values, converged to |residual| <= 0.02 m/yr:

1996-2010   GIS 1 +4.8394   GIS 90 +18.2545   residuals +0.006 / +0.016
2010-2024   GIS 1 +18.8657  GIS 90 +24.2358   residuals -0.005 / +0.001

Record and every probe: output/raw_runs/experiments/end-domain-boundaries/
2026-09-28-ends-resolved-lowess7/ (tables/ends.json).

--- LOWESS-10 option A values, SUPERSEDED 2026-09-28 ---
OPTION A, ADOPTED 2026-09-27 (Hannah). Re-solved on the
METRES island offset at the option A wave climate -- Hs 2.0 m, Tp 7.5 s,
asymmetry 0.6, high-angle 0.5, the HAT_hindcast_config defaults since the
same day -- on the edgeBE full-management nogroin run, against each
window's CoastSat LRR (GIS 1 raw, GIS 90 LOWESS-10), converged to
|residual| <= 0.05 m/yr:

1996-2010   GIS 1 +4.8394   GIS 90 +17.545   residuals +0.006 / -0.009
2010-2024   GIS 1 +18.8     GIS 90 +24.535   residuals -0.045 / -0.008

2010 GIS 1 was set by direct probes: its response is not monotonic
within ~0.1 m/yr. These values hold ONLY at those four wave settings;
change one and re-solve (HAT_resolve_ends_metres.py). Record and every
probe: output/raw_runs/experiments/end-domain-boundaries/
2026-09-27-ends-resolved-metres-offset/ (tables/ends.json). Option B is
HATTERAS_WAVE_OPTION_B below.

The two solves recorded below are the /10-OFFSET values (1996 +32.2 /
+10.0, 2010 +72.6 / +31.3, at Hs 2.5 / Tp 8 / asym 0.7 / ahf 0.1). They
are kept for the record; the matrix runs made on them are in
output/raw_runs/archive/2026-09-24-pre-metres/.

--- /10-offset solve, 1996, SUPERSEDED 2026-09-27 ---
SOLVED 2026-09-11, three Newton steps.

Fit exactly as the 1984 and 2004 end values were: model lrr_m_yr against
the target table's target_lrr_m_yr, on the edgeBE / road_bdm / groin-off
base run, with the same estimator on both sides of the residual. GIS 1 is
graded against the RAW per-domain mean and GIS 90 against the LOWESS-10
value, because that is the splice the rate figure draws.

step   GIS 1                          GIS 90
0    imposed  0.0  residual -2.930  imposed  0.0  residual -1.316
1    imposed 27.9  residual -0.480  imposed 12.5  residual +0.288
2    imposed 33.4  residual +0.125  imposed 10.3  residual +0.034
3    imposed 32.2  residual +0.003  imposed 10.0  residual +0.003

THE GAINS SIT INSIDE THE RANGE THE OTHER FOUR CASES SPAN. The local
secant through steps 0-1 gave d(LRR)/d(BE) = 0.088 at GIS 1 and 0.128 at
GIS 90, against 0.092 to 0.123 for the solved 1984 and 2004 cases. So the
tenth-of-what-you-impose behaviour holds here too, and these values are
again about ten times the misfit they close rather than a sediment
budget. Read them as boundary-artefact absorbers, nothing more.

A THIRD STEP WAS TAKEN, where 1984 and 2004 stopped at two. Not because
two was wrong: step 2 left +0.125 and +0.034, comparable to the -0.145
that stands at GIS 1 in period 2. It cost two minutes and removed the
question. Do not read the tighter convergence as a better-determined
number -- the target it converged onto carries the same uncertainty.

BOTH SIGNS ARE POSITIVE, unlike 1984, whose GIS 1 is -42.6. The southern
boundary needs sand added over 1996-2010 where it needed sand removed
over 1984-2004. That follows the observed target, which is +3.23 m/yr at
GIS 1 here; it is not evidence about Cape Point, which the model does not
represent (the calibrated fits are 0.0 through GIS 2-7).

SOLVED ON 1984-start/v2, WHILE THE 1984 PRESET IS ON v1 -- and that
turns out not to matter, which was worth one run rather than an
assumption either way (2026-09-12). Spending this pair on v1, the island
every published 1984 number was fitted on:

GIS 1 residual   GIS 90 residual
v2          +0.003           +0.003
v1          +0.048           +0.003

0.048 m/yr is a third of the -0.145 that stands unclosed at GIS 1 in
period 2, so these values transfer across the version within the
tolerance the published periods already accept. They do NOT need
re-solving when the 95 remaining v1 runs are moved to v2.

This says nothing about the INTERIOR. 65 of 90 domains differ in shape
between the two PRODUCTS, and this is a statement about two VERSIONS of
one product at two boundary cells. Do not generalise it.

Solve reproduced with:
scripts/input_prep/7-source-sink/2-calibrate/
be_edge_domain_solve.py --period 1996
```

```text
--- /10-offset solve, 2010, SUPERSEDED 2026-09-27 ---
SOLVED 2026-09-16, three Newton steps, the same protocol as 1996. Base
run: the 2010 matrix zeroBE / full_management / nourish / nogroin run,
2004-start v1, island offset 2010/v1, Hs 2.5.

step   GIS 1                          GIS 90
0    imposed  0.0  residual -8.430  imposed  0.0  residual -3.565
1    imposed 80.3  residual +0.919  imposed 34.0  residual +0.331
2    imposed 72.4  residual -0.022  imposed 31.1  residual -0.024
3    imposed 72.6  residual -0.003  imposed 31.3  residual +0.001

THE LARGEST END VALUE OF ANY PERIOD, at GIS 1. It is not a larger
artefact: the gains were 0.116 / 0.115 on the first secant, inside
the 0.09-0.13 every other case gave, so the value is again about ten
times its misfit. The misfit itself is what is large: the 2010-2024
CoastSat target at GIS 1 is +6.90 m/yr, against +3.23 in 1996-2010,
and the model produces -1.53 there unaided. Cape Point accreting at
that rate is not something the model represents (the calibrated fits
are 0.0 through GIS 2-7 in the two published periods), so this is the
boundary term supplying observed accretion the reach cannot make.
Read as a real flux it would be absurd; it is not one.

GIS 90 is +31.3 against +10.0 in 1996-2010, for the same reason:
the 2010-2024 target there is +2.22 m/yr, and the n115 extension
experiment already showed the north end grows when the observed
accretion at Pea Island is what it has to supply (+41.7 at GIS 115).

The probes are output/raw_runs/experiments/2026-09-16-edgesolve-2010/
(SOLVED names step3). Reproduced with:
be_edge_domain_solve.py --period 2010 --kind experiment
--run <base> --tag 2026-09-16-edgesolve-2010/base
--run <step> --tag 2026-09-16-edgesolve-2010/step<k> ...
```

```text
OPTION B, RECORDED, NOT WIRED (2026-09-27, Hannah). The one-parameter
alternative to option A: 2010-2024 at Hs 2.5 m (Tp, asymmetry and high-angle
as A), with its own ends, solved the same way; 1996-2010 is unchanged from A.
Raw share of the alongshore variation explained in 2010-2024: managed -122%
(A -135%), natural -538% (A -699%). No other single-parameter change helps
both scenarios. Neither window is fitted well; B only narrows the misfit.

NOTHING READS THIS. It holds the numbers so that running B never means
digging them out of a study folder. To run B, override the one wave field and
the two ends; the run name gains `waveHs2p5`, so it can never overwrite an
option A run:

HAT_START_YEAR=2010 HAT_HS=2.5 HAT_SOURCE_SINK_PRESET=edgeBE \
HAT_BE_OVERRIDE="1=8.0,90=40.399" python scripts/hatteras_ms/HAT_hindcast_1984_2024.py

Record: output/raw_runs/experiments/wave-climate/2026-09-27-wave-recommendation/
("Same vs period-specific"); ends in ends.json under `history` (Hs 2.5).
```

```text
AN EXTENDED GEOMETRY (2026-09-16) keeps only the standing values at domains
that are STILL ends -- GIS 90 is interior under n115 and must not carry
+10.0 -- and cannot pass the two checks below until its new end is solved:
the value arrives through HAT_BE_OVERRIDE, one Newton step at a time, and
section 4.3 of the runner refuses an edgeBE run whose end has none.
```

```text
An edge domain present but ZERO is the same failure the check above
exists to prevent, and the absence test does not catch it: edgeBE would
be byte-identical to zeroBE at that domain, which makes the two presets
indistinguishable there and collapses the edge secant to 0/0 with no
error. Found on 2026-08-28 while zeroing the interior of the calibrated
table -- GIS 1 is SLICED from that table, so emptying it would have
silently disarmed the edge preset.
```

```text
The split is only defensible while each side is actually solved against its
own preset. If a period ever gains an edgeBE D90 with no calibBE counterpart
-- or the two silently converge back to one number, which would mean the
split has stopped earning its complexity -- say so rather than let it pass.
```

```text
Canonical presets. Keys are the tokens that land in run-directory names, so
each one states which hypothesis was run -- "base" did not.
```

```text
Deprecated spellings, kept so older scripts keep running. resolve_be_preset()
maps these to the canonical key, and it is the canonical key that reaches the
run name -- an alias can never put a stale token in a directory name.
```

```text
GIS domains carrying NC-12. Domains 1-8 (Cape Point) have no road in the
modelled span.
```

```text
Permanent community zones. Roadway management is OFF here: inside a village
the road is a street network that is maintained rather than relocated, so
CASCADE's relocate-or-abandon logic does not describe it.

These are the PERMANENT settlement footprints only, deliberately narrower
than the beach-nourishment footprints, which extend past the villages.
```

```text
Per-domain road elevation, meters MHW-RELATIVE: the MEAN of the 2009 LiDAR
1 m cells in a 7 m corridor under the digitised 2004 NC-12 alignment.
Written by scripts/input_prep/4-mgmt-forcings/road_elevation/
HAT_road_elevation.py; see RoadElevation_audit.md beside the file.

NOT period-dependent, and that is deliberate -- but the REASON changed on
2026-08-26 and the old one is worth recording because it is now false.

It used to read: "there is one topography (2009) for every period, so one
elevation set serves both". There are two products now, and they do NOT agree
under the road. Sampling the 2004 alignment in a 3.5 m corridor on both:

2009-2014-1996 minus 2009-2014, corridor mean:  median +0.222 m
54 of 82 domains move more than 0.05 m; cell counts identical

Identical cell counts means ALACE REPLACED measured 2009 pavement rather
than filling holes in it. And +0.222 m is not a road: it is the island-wide
1996-vs-2009 survey offset, which mosaic_1984_audit.csv reports per domain at
median +0.255 m (p10 +0.14, p90 +0.33) and HAT_dem_1984_mosaic.py leaves
UNCORRECTED on purpose ("bias correction OFF, feathering OFF").

So a per-period road elevation built from each period's own DEM would push
1984 road_ele up ~0.22 m island-wide, and that increment would be the survey
offset, not a measured roadbed. A higher road is buried by overwash less
often, so it would reach the model. One file, built on the 2009-2014
baseline, keeps that offset out of a forcing.

Hence still no year in the filename -- a per-year name would imply a measured
change in roadbed height that no data supports.

The 2004 alignment is used rather than 1984 because it is the only digitised
line that lands on roadbed everywhere on the 2009 surface. Sampled along the
1984 line, the relocated domains (GIS 9-15, 84-87) return 2.37 m mean with a
within-domain sigma up to 1.70 m -- the abandoned corridor now lies UNDER the
foredune, so that sample is dune, not road. The 2004 value is both the better
measurement and the lower of the two.
```

```text
Historical NC-12 management events.

Relocations carry a DISPLACEMENT, not an absolute setback. CASCADE adds it to
whatever setback the model is carrying at the event year, which already
reflects the modelled shoreline retreat -- so the retreat is counted once. An
absolute setback referenced to the 1984 dune line would count it twice,
because the topography's own dune line has already moved landward by that
amount.

WHERE THE DISPLACEMENTS COME FROM -- measured, not typed.

Until 2026-08-20 these were eleven hand-entered literals attributed to a
1978->1997 cross-shore offset digitised in ArcGIS Pro. The 1997 line is not in
the repo -- only nc12_1978.geojson and nc12_2008.geojson are -- so those
numbers could not be re-derived, checked, or corrected. They are now read from
HAT_road_relocation_distance.py's per-domain measurement of the two lines that
ARE on disk.

WHICH COLUMN, AND WHY THE SIGNED ONE
mean_relocation_m        unsigned nearest-distance old line -> new line.
mean_signed_landward_m   the same displacement vectors projected onto the
landward normal. THIS ONE.
The unsigned column is dragged toward zero wherever the two lines share
vertices: the 2004 line was digitised by editing a copy of the 1984 one, so an
unedited stretch contributes an exact 0.000 m that is an editing artefact, not
a measurement. On GIS 9 and 87 that is 23% and 33% of samples, and the signed
mean is correspondingly LARGER (18.0 vs 13.8 m, 17.5 vs 11.8 m). Sign also
matters on its own terms: a relocation is landward by definition, so a column
that cannot express direction cannot contradict that claim.

ROUNDED TO WHOLE CELLS (Hannah, 2026-09-10: "a road can't relocate 17 m,
each grid cell is 10 x 10"). The measurement is a distance between two
digitised lines and is kept as measured in the CSV; what the model is FORCED
with is that distance rounded to the nearest 10 m, the Barrier3D cell, so a
prescribed move is a whole number of rows. roadway_manager floors the setback
to whole cells at placement anyway (int(road_setback / 10)); rounding here
makes the forcing say the same thing instead of carrying a fraction the grid
cannot hold. Nearest, not floor: 17.97 m is two cells' worth closer to 20
than to 10. No measured value sits on a 5 m boundary, so ties do not arise.
9  18.0 -> 20     10  46.9 -> 50     11  77.1 -> 80     12  68.4 -> 70
13 50.9 -> 50     14  24.7 -> 20     84  39.9 -> 40     85 108.7 -> 110
86 93.7 -> 90     87  17.5 -> 20

THE VINTAGE GAP. The lines were digitised off 1978 and 2008 imagery, so the
measured interval brackets BOTH events rather than either one. That is safe
only because the two events are disjoint in space -- 1989 moves GIS 84-87,
1999 moves GIS 9-15 -- so no domain's displacement is claimed twice. Do not
add a third relocation overlapping either span without revisiting this.
```

```text
Per-domain measurement written by
scripts/input_prep/4-mgmt-forcings/road_relocation/HAT_road_relocation_distance.py
Named by the two LINE vintages it was measured between (renamed from
1984_2004 on 2026-09-15, with the lines themselves), not by the periods those
lines stand in for.
```

```text
WHY THE 1999 EVENT STOPS AT GIS 14. The measurement classifies GIS 15
'relocated', but cannot say by how much or in which direction: the two
digitised lines CROSS inside that domain, so 56% of samples read landward and
44% seaward (sign_agreement 0.56), and the mean (+1.8 m) and median (-1.0 m)
disagree in SIGN -- that mean is the residue of two opposing populations, not
a displacement. Its mean magnitude, 4.0 m, is below the 5 m re-digitising
threshold; only its 12.2 m maximum is above, which is the sole reason it
classified 'relocated' rather than 'redigitized' at all. GIS 15 is therefore
left out of the 1999 event entirely rather than carried with a forced 0.0 m.
```

```text
WHAT THE ROUNDING COST, measured 2026-09-14. Rounding to whole cells moves
the 1999 event's landing at GIS 11 by +2.93 m (77.07 -> 80.0), and GIS 11 had
ONE CELL of margin against the drowning check. Re-running the 1984-2004
relocation arm under current code drowns NC-12 there in all eight calibBE and
edgeBE cells, where the runs stored on 2026-09-01 report it managed
throughout; the four zeroBE cells, which impose no background erosion and so
retreat less, do not drown. Forcing the unrounded displacements back reverses
it exactly: 0 drowned again.

The rounding is still right -- a prescribed move the model cannot represent at
sub-cell resolution should not pretend to -- but it is worth knowing that it
is what moved that domain over the line, and that GIS 11's margin is one cell
either way. See output/raw_runs/arms/recode-20260914/ for the comparison.
```

```text
Independent check on a relocated setback: the MEASURED 2004 setback for the
same domain. Both relocation events precede 2004, so a correctly displaced
setback should land near the 2004 same-year measurement.
Reporting only -- nothing reads this to decide anything.

READ FROM THE FILE PERIOD 2 ACTUALLY RUNS, not typed out here. This used to
hold eleven literals copied from the LEGACY old-method RoadSetback_2004.csv
(GIS 9 = 89 m, 10 = 83 m, 11 = 81 m ...). That copy survived the move to the
dune-start method on 2026-08-18 and both topography bumps after it, so the
printed check was comparing a dune-start setback against a number measured to
a different feature -- median +23 m landward, and on a topography the run no
longer uses. The same domains read 50 / 40 / 20 m in the file period 2
actually loads. A check that disagrees by construction is worse than none:
it invites a correct relocation to be read as a failed one.

Deriving it has the same motive as slicing HATTERAS_BE_RATES_EDGE out of the
calibrated preset -- the two can no longer drift apart, and bumping the
topography version moves this with the setbacks it is checking.
```

```text
Source: Hatteras_Management_Timelines.xlsx -> Nourishment_Timeline sheet.
Volumes are the reported project totals in cubic yards; cascade_pipeline
converts to m^3/m and spreads them evenly across each project's domains.

All three projects fall in 2004-2024, so a 1984-2004 run builds an empty
schedule from this same list -- no period-keying needed here.

The project extents are wider than HATTERAS_COMMUNITY_ZONES on purpose:
these are engineering footprints, and both Buxton and Rodanthe extend past
the settlement into the NC-12 corridor. That is what puts GIS 9-15 and
85-88 under both the roadway and beach-dune managers.
VOLUMES AND EXTENTS ARE FROM THE PROJECT RECORD, and both halves matter:
CASCADE applies volume/length, so a right volume on a short footprint is as
wrong as a wrong volume. Checked against the record 2026-08-22; Buxton was
already correct, Avon and Rodanthe were not.

project    was              now              record
Rodanthe   619 m^3/m        408 m^3/m        ~380 m^3/m
Buxton     183 m^3/m        183 m^3/m        ~197 m^3/m   (unchanged)
Avon       841 m^3/m        191 m^3/m        190-216 m^3/m

Avon was wrong on BOTH axes -- 2.2x the volume placed, over 55% of the
footprint -- which compounded to 3.9-4.4x the real fill density and a 68 m
instantaneous shoreline step. That step is what BRIE's Crank-Nicolson solve
rings on; see compute_lrr and HAT_hindcast_methods.md section 12.
```

```text
USACE/NCDOT emergency fill at the Mirlo Beach "S-curves", ~2 miles
(3.2 km) immediately NORTH of Rodanthe village, whose community zone
ends at GIS 83. Six domains is 3.0 km (1.86 mi) rather than the
6.4 domains 2 miles would take: GIS 90 is the locked north-end
domain carrying the source/sink boundary value (tens of m/yr), and
putting fill into a domain whose rate was pinned rather than
modelled would make the fill and the boundary condition
indistinguishable. The 7% shortfall is the cost of that exclusion
and is stated rather than absorbed into the volume.
```

```text
2.9 mi (4.7 km) north from the oceanfront groin at the lighthouse,
which sits at GIS 5.5 -- so 6-15, which is what was already here.
Left untouched: it is the one project whose configured density
(183 m^3/m) already matched the record (197 m^3/m).
```

```text
Placed from 3,000 ft north of Avon Pier (Due East Road) south to
Askins Creek North Drive at the village's southern boundary. The
pier is GIS 26 (HATTERAS_ANNOTATIONS.piers), 3,000 ft is 1.8
domains, and the southern village boundary is GIS 21 -- so 21-28,
4.0 km, matching the ~2.5 mi the record gives. The previous 23-26
was the MIDDLE of this footprint, missing ~2 km to the south and
~1 km to the north. Still entirely inside the Avon community zone
(21-31), so the beach_dune_manager footprint is unchanged by this.
```

```text
Overwash filtering percent for developed ground. CASCADE's overwash_filter
is a PERCENT -- filter_overwash divides by 100 -- and Rogers et al. (2015)
give 40-90%, residential to commercial. Hatteras' villages are residential,
so 40 is the low end of that range.

Earlier versions of this pipeline passed 0.4 here, on the assumption it was
a fraction. That filtered 0.4% of overwash, which is indistinguishable from
no filtering at all; BeachDuneConfig now rejects the fraction scale.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_island_offset_file()`**

```text
Path of the padded offset file, relative to INIT_ROOT.

120 domains in the base geometry. An extended geometry reads its own
build under <year>/ext/<geometry>/ (island_offset_hybrid.py --geometry),
which is not a version and does not move CURRENT; the path is returned
unchecked because every period resolves here at import and only the
period being run needs the file to exist -- the runner checks that.
```

**`island_offset_version()`**

```text
The version segment of the offset file this period resolves to.

"duneline/v2" since 2026-09-22, when every build moved under the SOURCE it
was measured from; "v2" before that, and "flat" for the unversioned layout
older still. Recorded in run metadata and run_index.csv (2026-09-15)
because a v1 run and a v2 run are otherwise identical on disk: the run name
carries no offset token and the file name is the same in every version
folder. Resolves through _island_offset_file so it can never disagree with
the file that was read.

READING AN OLD RUN: a token with no "/" is from before the sources were
split, and every build then was dune-derived -- "v1" means "duneline/v1".
```

**`resolve_be_preset()`**

```text
Resolves a source/sink preset name to its canonical key and rates.

Args:
    name: A canonical preset key or a deprecated alias.

Returns:
    A (canonical_name, rates_by_period) tuple. rates_by_period maps start
    year to a sparse {gis_id: rate_m_yr} dict.

Raises:
    ValueError: If the name is neither a canonical key nor an alias.
```

**`be_rates()`**

```text
The per-domain rates one preset supplies for one period.

WHY THIS IS NOT A PLAIN DICT LOOKUP. A period can be WIRED before it is
CALIBRATED, and since 2026-09-11 two of them are: all four periods resolve
their forcing files, but only 1984 and 2004 have solved edge and calibrated
presets. The other two carry zeroBE alone until their end domains are
re-solved against their OWN CoastSat target.

Indexing the preset dict directly gives a bare KeyError on the year, which
reads like a typo in the settings file rather than what it is -- a fit that
has not been done yet. This says which, and what would fix it.

Args:
    name: A canonical preset key or a deprecated alias.
    start_year: The period's start year.

Returns:
    A sparse {gis_id: rate_m_yr} dict. An absent domain is 0.0 m/yr.

Raises:
    ValueError: If the preset is unknown, or is not solved for that period.
```

**`_measured_displacements()`**

```text
Reads the measured landward displacement for one relocation event.

Args:
    gis_domains: The GIS domains the event moves. Hand-specified per event,
        NOT taken from the file: the measurement classifies ~40 domains
        'relocated' across the whole island, and which of them belong to
        which historical event is a fact about NC-12, not about the lines.

Returns:
    A {gis: displacement_m} dict, holding mean_signed_landward_m from the
    measurement rounded to the nearest CELL_M (see ROUNDED TO WHOLE CELLS
    above); `measured_displacement_m` returns the unrounded value.

Raises:
    ValueError: If a domain is absent from the measurement, or the
        measurement classifies it 'no_edit' / 'redigitized'. Either means
        the CSV cannot support a relocation there, and forcing one anyway
        would put a number the data does not contain into the model.
```

**`measured_displacement_m()`**

```text
The UNROUNDED measured displacement for one domain, for labels and
reports that want to show the measurement beside the forcing.
```

**`_relocation_check()`**

```text
Reads the measured setback for every domain a relocation event moves.

The domain list is taken from HATTERAS_ROAD_EVENTS rather than retyped, so
adding a relocation adds its cross-check automatically.

Args:
    period_start: Start year whose road_setback_file supplies the measured
        values. Both relocation events precede 2004, so period 2 is the
        first same-year measurement that postdates them.

Returns:
    A {gis: measured_setback_m} dict.

Raises:
    ValueError: If a relocated domain is absent from the setback file.
        load_road_setbacks fills an absent domain with 0.0, which would
        read here as "the road sits on the dune line" -- a plausible-
        looking number for a measurement that was never made.
```

</details>
