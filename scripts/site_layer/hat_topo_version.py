# ==============================================================================
# hat_topo_version.py
#
# WHICH Barrier3D domains does a script read - which PRODUCT, and which VERSION
# within it?
#
# WHY THIS EXISTS
#   Four road scripts used to hardcode ".../2009-dune-topo/2009_v3/topography".
#   When the dune windows were re-picked into 2009_v4 (2026-08-19) they kept
#   reading v3 interiors while consuming v4 setbacks. Nothing errored. 18
#   domains had different interiors, 10 of them a different SHAPE - including
#   D79 and D80, two of the three roadways the relocation logic acts on. Every
#   drown verdict and placement number for those was computed on the wrong grid.
#
#   So the location is resolved ONCE, here, and a name that does not exist on
#   disk is an immediate, loud error rather than a silently stale read.
#
# WHAT CHANGED 2026-08-25 - THERE ARE NOW TWO PRODUCTS
#   The tree was stage-keyed and held exactly one topography, which BOTH
#   hindcast periods read:
#
#       1-barrier3d-domains/2009-dune-topo/2009_v5/{topography,dunes}
#
#   It is now PERIOD-keyed, because the two periods start from different DEMs:
#
#       1-barrier3d-domains/
#           1984-start/    from DEM 2009-2014-1996   dune-topo/<version>/
#           2004-start/    from DEM 2009-2014        dune-topo/<version>/
#           forecast/      from a 2025 DEM, later
#           buffer/        shared
#           superseded/
#
#   2009_v5 became 2004-start/dune-topo/v1 (renamed from v5 on 2026-08-26 so
#   version numbers restart per product) - it was always built from
#   the 2009+2014 DEM, which is what the 2004 start uses. v3 and v4 are under
#   superseded/ (they were picked against the UNFILLED DEM).
#
# WHY topo_dirs() STILL DEFAULTS TO A PRODUCT
#   Every existing caller - the road tree, the groin sweep, the poster script -
#   calls topo_dirs() with no arguments. Defaulting to DEFAULT_PRODUCT keeps
#   all of them resolving exactly what they resolved before the restructure, so
#   this change moves no road number and no published figure. The RUNNER is the
#   one caller that now passes a product, from HATTERAS_PERIODS[start]
#   ["topo_product"].
#
#   SETTLED 2026-08-26. That note used to read "a live question": the road
#   scripts measured setbacks against 2004-start for BOTH periods. They no
#   longer do. Every vintage now resolves its own product through YEAR_PRODUCT
#   below, and RoadSetback_1984_dunestart.csv was re-measured on 1984-start.
#   What made it urgent rather than tidy: all 90 domains differ between the two
#   products and 65 have a different interior SHAPE.
#
# HOW THE VERSION IS CHOSEN, in order
#   1. an explicit override= argument
#   2. HAT_TOPO_VERSION_<PRODUCT> in the environment, e.g.
#      HAT_TOPO_VERSION_1984_START=v5. Product-scoped, because
#      both products have a "v1" and a global override would silently resolve to
#      a real but wrong directory for the other period.
#   3. a CURRENT file in the product's dune-topo/ directory
#   4. the extractor's VERSION, but only if the extractor is currently pointed
#      at the SAME product
#   5. the only version present, if there is exactly one
#   otherwise: raise, listing what is on disk
#
#   CURRENT OUTRANKS THE EXTRACTOR LITERAL (swapped 2026-09-04). The extractor's
#   VERSION says what the extractor WRITES; CURRENT says what everyone READS.
#   They used to be one literal doing both jobs, which meant the only way to
#   make a layered version (the 1984-start layers v3-v8 of 2026-09-04 -- built
#   ON v2, never BY the extractor; DELETED 2026-09-07, only unmodified
#   extractions are kept) the
#   default was to edit the extractor to a name it would then overwrite on its
#   next run. So CURRENT existed, recorded intent, and was inert; the 1984-start
#   README carried a paragraph explaining that it did nothing. Now a product
#   with a CURRENT file reads that version, and a product without one still
#   follows the extractor ("bump VERSION and the tree follows" holds where no
#   CURRENT has been written - 2004-start today). A fresh extraction with
#   CURRENT still naming the old one is therefore NOT adopted until CURRENT is changed,
#   and that is the point: extracting and adopting are two decisions.
#
#   THE ENV RULE OUTRANKS THE EXTRACTOR (added 2026-09-02, after it cost a run).
#   A batch that selected its arm by writing CURRENT was ignored, because rule 3
#   fired first and the extractor was sitting on the same product. Two arms of a
#   three-arm experiment silently duplicated the control, exit code 0. CURRENT
#   is a persistent shared DEFAULT; per-run selection needs something that does
#   not mutate state every other reader sees.
#
# USAGE
#     from site_layer.hat_topo_version import topo_dirs
#     TOPO_DIR, DUNE_DIR, RUN_NAME = topo_dirs()                  # 2004-start
#     TOPO_DIR, DUNE_DIR, RUN_NAME = topo_dirs("1984-start")
#     TOPO_DIR, DUNE_DIR, RUN_NAME = topo_dirs("2004-start", override="v1")
# ==============================================================================

from __future__ import annotations

import os
import re
from pathlib import Path

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents
                    if (_p / "pyproject.toml").exists())
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
DOMAIN_ROOT = INIT_ROOT / "1-barrier3d-domains"

# Shared across products: the padding domains Barrier3D needs either side of
# the 90 real ones. Not per-period - the buffer is not a survey of anything.
BUFFER_DIR = DOMAIN_ROOT / "buffer"

# What topo_dirs() resolves when no product is named. See the note above: this
# is the pre-restructure behaviour, kept so nothing that was not asked to move
# moves.
DEFAULT_PRODUCT = "2004-start"

PRODUCTS = ("1984-start", "2004-start", "forecast")

# WHICH PRODUCT DOES A HINDCAST PERIOD READ. The single definition (2026-08-26).
#
# Keyed by the period's START YEAR, which is also the year that labels every
# road vintage in 4-mgmt-forcings: RoadSetback_1984_dunestart.csv is measured
# from row 0 of the 1984-start extraction, RoadSetback_2004_dunestart.csv from
# row 0 of 2004-start.
#
# WHY IT LIVES HERE. It was written out three times and omitted five times.
# HAT_road_offset_from_dune_start.py carried its own YEAR_PRODUCT literal,
# HAT_road_setback_audit.py spelled the product per scenario, and
# hatteras_site_config.py carries it per period as "topo_product" -- while
# HAT_road_placement_on_domains.py, HAT_road_method_diagnostic.py and
# HAT_road_domain_views.py looped over BOTH years against a single
# module-level topo_dirs(), i.e. DEFAULT_PRODUCT for both.
#
# That last one is not a cosmetic duplication. Between the two products ALL 90
# domains differ and 65 have a different interior SHAPE -- GIS 11 is 165 rows
# on 1984-start and 157 on 2004-start. It is the v3/v4 failure this file was
# written for, four times larger, and it produced no error: a 1984 setback
# scored against a 2004-start interior is simply a different island.
#
# So the pairing is defined once, here, beside the resolver that consumes it.
# hatteras_site_config.py imports it rather than repeating it; every road
# script resolves through it rather than defaulting.
YEAR_PRODUCT = {
    1984: "1984-start",
    2004: "2004-start",
    # The two periods added 2026-09-11 SHARE these products rather than owning
    # one each. A product is a DEM composition, not a period label: 1984-start
    # is the surface carrying the 1996 ALACE graft, which is the survey nearest
    # a 1996 start, and 2004-start is the 2009-plus-2014 surface, which is the
    # closest vintage match of any period to a 2010 start.
    #
    # So a product name now names the period it was FIRST built for, not the
    # only one that reads it. Ask this mapping rather than inferring a product
    # from a year, and ask product_for_year() rather than indexing directly.
    1996: "1984-start",
    2010: "2004-start",
}

# WHICH NC-12 LINE A PERIOD'S ROAD IS MEASURED FROM. The single definition
# (2026-09-15).
#
# The two digitised centrelines under 4-mgmt-forcing/road_offset/raw_offset/
# were exported off 1978 and 2008 imagery. Until 2026-09-15 their folders, the
# rasterised masks under raster/, and the relocation measurement under
# road_relocation/ were all named 1984 and 2004 -- the period starts they stand
# in for -- which put two different axes under one integer: dunestart_offset/
# 1984 meant "the 1984 START", raw_offset/1984 meant "the 1978 LINE". They are
# now named by the line's TRUE vintage, and this map is the only place a start
# year is paired with a line. Ask road_line_for_year() rather than re-spelling
# the pairing; a period start passed where a line vintage is expected is a
# loud error (road_line_file), not a missing-file error later.
#
# 1996 and 2010 have no line of their own. 1996 is the 1978 line with the 1989
# Pea Island relocation applied; 2010 is the 2008 line unchanged, because no
# relocation falls between 2004 and 2010. See ROAD_SETBACK_KIND.
ROAD_LINE_VINTAGES = (1978, 2008)
ROAD_LINE_FOR_YEAR = {
    1984: 1978,
    1996: 1978,
    2004: 2008,
    2010: 2008,
}

# WHICH SETBACK FOLDERS HOLD A MEASUREMENT AND WHICH A DERIVATION. The tree
# under road_offset/dunestart_offset/ is split the same way (2026-09-15):
#
#     dunestart_offset/measured/<year>/   measured on the period's line against
#                                         row 0 of the period's own extraction
#                                         (HAT_road_offset_from_dune_start.py)
#     dunestart_offset/derived/<year>/    built FROM a measured file
#                                         (HAT_road_setback_derived_vintages.py)
#
# so a derived file cannot be mistaken for a measurement by its address alone.
# Before the split all four sat as flat siblings and only PROVENANCE.md said
# which two were copies.
ROAD_SETBACK_KIND = {
    1984: "measured",
    2004: "measured",
    1996: "derived",
    2010: "derived",
}

MGMT_ROOT = INIT_ROOT / "4-mgmt-forcing"
ROADS_ROOT = MGMT_ROOT / "road_offset"


def road_line_for_year(year: int) -> int:
    """The NC-12 line vintage (1978 or 2008) a period start year reads."""
    try:
        return ROAD_LINE_FOR_YEAR[int(year)]
    except (KeyError, TypeError, ValueError):
        raise SystemExit(
            f"\nno NC-12 line known for period start {year!r}. "
            f"Known: {', '.join(str(y) for y in sorted(ROAD_LINE_FOR_YEAR))}\n"
            f"Add it to ROAD_LINE_FOR_YEAR in {__file__}.\n")


def _check_line_vintage(vintage) -> int:
    v = int(vintage)
    if v not in ROAD_LINE_VINTAGES:
        raise SystemExit(
            f"\n{vintage!r} is not a road LINE vintage. The digitised lines are "
            f"{ROAD_LINE_VINTAGES}; a period start year goes through "
            f"road_line_for_year() first. (Before 2026-09-15 the lines were "
            f"filed under 1984 and 2004, the periods they stand in for.)\n")
    return v


def road_line_file(vintage: int) -> Path:
    """The digitised NC-12 centreline of one LINE vintage."""
    v = _check_line_vintage(vintage)
    return ROADS_ROOT / "raw_offset" / str(v) / f"nc12_{v}.geojson"


def road_mask_dir(vintage: int) -> Path:
    """Where HAT_rasterize_road_to_domains.py put one line vintage's masks."""
    v = _check_line_vintage(vintage)
    return ROADS_ROOT / "raster" / str(v) / "masks"


def road_mask_file(vintage: int, domain: int) -> Path:
    """One domain's rasterised road mask for one line vintage."""
    v = _check_line_vintage(vintage)
    return road_mask_dir(v) / f"domain_{int(domain)}_road_{v}.npy"


def road_setback_dir(year: int) -> Path:
    """The folder holding a period start year's setback files."""
    try:
        kind = ROAD_SETBACK_KIND[int(year)]
    except (KeyError, TypeError, ValueError):
        raise SystemExit(
            f"\nno road setback known for period start {year!r}. "
            f"Known: {', '.join(str(y) for y in sorted(ROAD_SETBACK_KIND))}\n"
            f"Add it to ROAD_SETBACK_KIND and ROAD_LINE_FOR_YEAR in "
            f"{__file__}.\n")
    return ROADS_ROOT / "dunestart_offset" / kind / str(int(year))


def road_setback_file(year: int) -> Path:
    """The model-facing 2-row RoadSetback_<year>_dunestart.csv for a period."""
    return road_setback_dir(year) / f"RoadSetback_{int(year)}_dunestart.csv"


def road_setback_relpath(year: int) -> str:
    """road_setback_file() relative to INIT_ROOT, POSIX-style -- the form
    HATTERAS_PERIODS[year]["road_setback_file"] carries."""
    return road_setback_file(year).relative_to(INIT_ROOT).as_posix()


# THE DUNE LINES, BY VINTAGE (2026-09-15). The same rule as the road lines:
# a raw per-transect file under 2-brie-offset/raw_offsets/ is named for the
# IMAGERY VINTAGE of the digitised line it came from, and this map is the only
# place a period year is paired with a vintage. Until 2026-09-15 the 1996
# start read a byte-identical COPY of the 1997 file filed under the name
# 1996_...; the copy is gone and the pairing lives here (Hannah: "a
# DUNE_LINE_FOR_YEAR table, no copies"). When the 2010 and 2024 lines are
# digitised, add them here under the year of their IMAGERY (2009 or 2023, if
# that is what they are traced from), never under the period year.
#
# A vintage may have more than one digitisation (duneline_1997.geojson and
# duneline_1997_v2.geojson); the raw file under the vintage's name is the
# CURRENT one, and each build under 2-brie-offset/<year>/v<n>/ keeps a copy of
# the raw it was made from, so an older build is always reproducible.
BRIE_ROOT = INIT_ROOT / "2-brie-offset"
RAW_OFFSET_DIR = BRIE_ROOT / "raw_offsets"
DUNELINE_DIR = BRIE_ROOT / "dunelines"
DUNE_LINE_FOR_YEAR = {
    1984: 1984,
    1996: 1997,   # no 1996 imagery; the nearest island-wide survey
    2004: 2004,
    2010: 2009,   # no 2010 aerial imagery (Hannah, 2026-09-15); the 2009 line
    2024: 2023,   # the 2023 NOAA imagery (D:\Hatteras_GIS\Aerial3); end year of 2004-2024 and 2010-2024
}


def dune_line_for_year(year, strict: bool = True):
    """The dune-line vintage a period start OR end year reads. With
    `strict=False` an unknown year returns None instead of exiting, for the
    end-year target, which is allowed to be missing."""
    try:
        return DUNE_LINE_FOR_YEAR[int(year)]
    except (KeyError, TypeError, ValueError):
        if not strict:
            return None
        known = ", ".join(f"{y} -> {v}" for y, v in sorted(DUNE_LINE_FOR_YEAR.items()))
        raise SystemExit(
            f"\nno dune line known for period year {year!r}. Known: {known}\n"
            f"Add it to DUNE_LINE_FOR_YEAR in {__file__} under the imagery vintage.\n")


def dune_raw_file(vintage) -> Path:
    """raw_offsets/<vintage>_duneline_offset_raw.csv, the current build of
    one LINE vintage's per-transect stations."""
    return RAW_OFFSET_DIR / f"{int(vintage)}_duneline_offset_raw.csv"


def dune_raw_file_for_year(year, strict: bool = True):
    """The raw file a period year reads, through DUNE_LINE_FOR_YEAR."""
    v = dune_line_for_year(year, strict=strict)
    return None if v is None else dune_raw_file(v)


def duneline_geojson(vintage, version: str | None = None) -> Path:
    """2-brie-offset/dunelines/duneline_<vintage>[_<version>].geojson."""
    suffix = f"_{version}" if version else ""
    return DUNELINE_DIR / f"duneline_{int(vintage)}{suffix}.geojson"


def product_for_year(year: int) -> str:
    """The topography product a period start year reads. Raises if unknown.

    Use this in any loop over vintages. The failure mode it removes is a loop
    body that resolves topo_dirs() once, outside the loop, and silently gives
    every year the same interiors.
    """
    try:
        return YEAR_PRODUCT[int(year)]
    except (KeyError, TypeError, ValueError):
        raise SystemExit(
            f"\nno topography product known for year {year!r}. "
            f"Known: {', '.join(str(y) for y in sorted(YEAR_PRODUCT))}\n"
            f"Add it to YEAR_PRODUCT in {__file__} -- and to HATTERAS_PERIODS "
            f"in hatteras_site_config.py, which imports this mapping.\n")


def year_for_product(product: str, strict: bool = True):
    """The period start year a topography product belongs to.

    The inverse of product_for_year(), and it lives here for the same reason
    the forward map does: the pairing is defined ONCE. A caller that needs
    "which year is this product" -- the extractor's figure code does -- must
    not re-spell {"1984-start": 1984} locally, because a third product added
    to YEAR_PRODUCT would then update one direction and not the other.

    strict=False returns None instead of raising, for products that are
    legitimately not a hindcast period ("forecast", "buffer"). A caller using
    it must have a defined behaviour for None -- see PRODUCT_YEAR in
    HAT_dune_topo_extractor.py, which falls back to plotting every year.
    """
    for year, prod in YEAR_PRODUCT.items():
        if prod == product:
            return int(year)
    if not strict:
        return None
    raise SystemExit(
        f"\nno period year known for topography product {product!r}. "
        f"Known: {', '.join(sorted(YEAR_PRODUCT.values()))}\n"
        f"Add it to YEAR_PRODUCT in {__file__}.\n")

# THE ARRAYS HAVE NO YEAR TAG (2026-08-26).
#
# They are domain_<N>_topography.npy / _dune.npy / _nodata.npy. There is no
# year in the name, and there should not be one.
#
# A tag was tried, briefly, both ways. The name was the literal "2009" for a
# long time, which was simply false - 2004-start is the 2009+2014 mosaic and
# 1984-start is 2009+2014+1996, so neither product is a 2009 DEM. Retagging by
# period ("1984"/"2004") was then tried and reverted the same day, because the
# tag turned out to be a DISTRIBUTED INVARIANT with no single place to fix and
# no single grep to audit. Twelve live scripts build these paths, and the tag
# reached them four different ways:
#
#     TOPO_DUNE_INIT_YEAR = "2009" then interpolated      5 scripts
#     the bare literal inline in the f-string             2 scripts
#     ext.TAG, imported from the extractor                1 script
#     globbed away as domain_*_topography_*.npy           2 scripts
#
# The retag found the first form, missed the second, and broke both of those
# scripts - they resolved the right DIRECTORY and then asked for a file that no
# longer existed. That is the argument in one sentence.
#
# WHAT THE TAG COULD NOT DO ANYWAY. It cannot catch a period mix-up. The tag
# and the directory both derive from the same `product`, so a wrong product
# gives a wrong directory AND a matching wrong tag - consistent, silent, no
# error. What guards that is the runner's boot/run product assertion, not the
# filename. And for the two scripts that glob the tag away, a stray file from
# the other period would make the glob match twice and pick arbitrarily, which
# is worse than no tag at all.
#
# The period lives in the DIRECTORY, which every caller has to get right
# regardless, and in each run's RUN_MANIFEST.txt. The buffer arrays have never
# carried a tag and have never been confused.
#
# CALLERS SHOULD NOT BUILD THESE NAMES. Use domain_arrays() or array_path()
# below, so the directory and the filename come from one place.

# The one that has ALONGSHORE_FLIP = True. Three other copies of this file exist
# in the repo and all are unflipped -- see the note in
# HAT_road_offset_from_dune_start.py.
EXTRACTOR = (PROJECT_ROOT / "scripts" / "input_prep" / "1-barrier3d-domains" / "1-extraction"
             / "HAT_dune_topo_extractor.py")


def _literal(name: str, text: str) -> str | None:
    """Read `NAME = "value"` out of the extractor source, or None."""
    m = re.search(rf'^{name}\s*=\s*["\']([^"\']+)["\']', text, re.MULTILINE)
    return m.group(1) if m else None


def _extractor_state() -> tuple[str | None, str | None]:
    """(TOPO_PRODUCT, VERSION) the extractor is currently configured for.

    PARSED rather than imported: importing the extractor pulls in matplotlib
    and a windowing backend, which is a heavy and fragile dependency for an
    audit that needs neither. Only two literals are wanted.
    """
    if not EXTRACTOR.is_file():
        return None, None
    text = EXTRACTOR.read_text(encoding="utf-8")
    return _literal("TOPO_PRODUCT", text), _literal("VERSION", text)


def product_dir(product: str) -> Path:
    if product not in PRODUCTS:
        raise SystemExit(
            f"\nunknown topography product {product!r}. Known: "
            f"{', '.join(PRODUCTS)}\n")
    return DOMAIN_ROOT / product


def dune_topo_root(product: str) -> Path:
    return product_dir(product) / "dune-topo"


# THE TWO HALVES OF STAGE 1 (2026-09-09, Hannah). Under each product,
# `1-extraction/` holds what the extractor reads and records (npy-arrays,
# npy-arrays_survey, picks, the aerial review of the holes, the retired
# experiments) and `2-domain-reconstruction-1984/` the 1984 reconstruction in
# its six steps; `dune-topo/` stays at the product root because BOTH halves
# write versions into it and it is what the runner loads. The scripts folder
# scripts/input_prep/1-barrier3d-domains/ is split the same way.
EXTRACTION_SUB = "1-extraction"


def extraction_dir(product: str) -> Path:
    """The extraction half of a product: arrays, picks and the audit records."""
    return product_dir(product) / EXTRACTION_SUB


def npy_dirs(product: str) -> tuple[Path, Path]:
    """(elevation arrays, survey/provenance arrays) the extractor reads."""
    d = extraction_dir(product)
    return d / "npy-arrays", d / "npy-arrays_survey"


def picks_dir(product: str) -> Path:
    """Where the dune-search windows of a product live (one JSON per version)."""
    return extraction_dir(product) / "picks"


# THE SEAWARD-ROW-INSERT FOLDER, and the paths that hang off it.
#
# ONE definition, because eight plotting scripts and two measurement scripts
# used to build these by hand - and two of them WRITE.
#
# THE LAYOUT IS NOT SYMMETRIC BETWEEN PRODUCTS, deliberately. On 2026-09-03
# everything belonging to the 1984-start seaward-row insert - the measurement of
# N, the scope report, the fill comparison and every figure - was consolidated
# under `2-domain-reconstruction-1984/`. 2004-start has no insert work and no such folder,
# so its dune-line measurements stay at the product root.
#
# The asymmetry is the price of that consolidation. It is contained here so a
# caller cannot get it wrong, and so a future product does not inherit it by
# accident: anything not listed gets the plain layout.
_INSERT_SCOPE = {"1984-start": "2-domain-reconstruction-1984"}


def insert_scope_dir(product: str) -> Path:
    """The seaward-row-insert folder. Raises for a product that has none."""
    sub = _INSERT_SCOPE.get(product)
    if sub is None:
        raise SystemExit(
            f"\n{product!r} has no 2-domain-reconstruction-1984 folder. Only "
            f"{', '.join(_INSERT_SCOPE)} carries the seaward-row insert.\n")
    return product_dir(product) / sub


# The four sections of the insert figures folder, in the order the argument
# runs: where the insert lands, where N came from, what the rows are made of,
# and what the result looks like. Numbered so a directory listing reads in that
# order, matching the numbered layout of data/hatteras_init itself.
#
# WHY THIS IS HERE AND NOT IN THE PLOTTERS. Thirteen figures in one flat folder
# is a dump, and a folder a plotter re-scatters on every run cannot be tidied by
# moving files. Naming the section at the call site - and resolving it here - is
# what makes the layout survive a re-plot. A section that is not one of these is
# a typo, and raises rather than silently creating a new folder.
# Renumbered 2026-09-07 to the order the argument runs: measure the dune-line
# shift, turn it into a footprint of rows (two placements of the same rows:
# seaward/, behind-road/; placement-independent figures at the step's root),
# argue the fill, look at the result. 6-result is reserved: nothing has been
# built and run on the footprint yet. The record figures of the deleted layers
# sit in superseded-layers/ and the irreproducible pre-re-pick ones in frozen/;
# neither is a section a plotter may write to (2026-09-08).
INSERT_FIGURE_SECTIONS = ("1-measurement", "2-extent", "3-placement", "4-fill", "5-build", "6-result")
# Inside a section (2026-09-08, Hannah): island-wide figures go in `island/`
# (or a named subfolder), and EVERY figure that shows one example domain goes
# in `rows-added/` or `rows-removed/` by the sign of that domain's N in the
# footprint table - never at the section root.
# Six steps since 2026-09-09, in the order the ARGUMENT runs (Hannah): how far
# the dune line moved, how many rows, WHERE they go (the two candidate
# placements, the road check, and the imagery review that decides), what they
# contain, the version built, what the model does. 5-build holds no figures.
SIGN_SUBFOLDERS = ("rows-added", "rows-removed", "unchanged")
INSERT_FIGURE_SUBFOLDERS = {
    "1-measurement": SIGN_SUBFOLDERS,
    "2-extent": ("island",),
    "3-placement": ("seaward", "behind-road",
                    "road-check", *(f"road-check/{s}" for s in SIGN_SUBFOLDERS),
                    "imagery-review", "imagery-review/island", *(f"imagery-review/{s}" for s in SIGN_SUBFOLDERS)),
    "4-fill": ("island", *SIGN_SUBFOLDERS),
    "5-build": (),
    "6-result": ("island", *SIGN_SUBFOLDERS),
}

# The DATA of the same steps (2026-09-09, Hannah: "organize this under
# subfolders"): the tables and reports each step writes sit in a subfolder of
# 2-domain-reconstruction-1984/ named like its figure section, so a listing of the data
# folder reads in the same order as figures/. Root keeps README.md,
# DUNE_TOPO_VERSION_GUIDE.md and figures/. Resolved here for the same reason
# the figure folders are: one definition, no hand-built paths in the scripts.
INSERT_SCOPE_STEPS = INSERT_FIGURE_SECTIONS


def insert_scope_step(product: str, step: str, sub: str | None = None) -> Path:
    """The data folder of one step of the insert work (tables, reports), or a
    named subfolder of it (3-placement has `imagery-review`). Created if
    absent. `step` is one of INSERT_SCOPE_STEPS; anything else raises."""
    if step not in INSERT_SCOPE_STEPS:
        raise SystemExit(
            f"\n{step!r} is not an insert step. Use one of {', '.join(INSERT_SCOPE_STEPS)}.\n")
    d = insert_scope_dir(product) / step
    if sub is not None:
        d = d / sub
    d.mkdir(parents=True, exist_ok=True)
    return d


def rows_sign_sub(n_cells: int) -> str:
    """The sign subfolder for a domain whose footprint is `n_cells` rows."""
    return "rows-added" if n_cells > 0 else ("rows-removed" if n_cells < 0 else "unchanged")


def footprint_n_cells(product: str, domain: int) -> int:
    """N for one domain from the footprint table (footprint_1984_by_domain.csv)."""
    import csv
    p = insert_scope_step(product, "2-extent") / "footprint_1984_by_domain.csv"
    if not p.is_file():
        raise SystemExit(f"\n{p} not found - run HAT_footprint_1984.py first\n")
    for r in csv.DictReader(open(p, newline="")):
        if int(r["domain"]) == int(domain):
            return int(r["n_cells"])
    raise SystemExit(f"\ndomain {domain} is not in {p}\n")


def insert_figures_dir_for_domain(product: str, section: str, domain: int,
                                  under: str | None = None) -> Path:
    """Where a figure of ONE example domain goes: the section's rows-added/ or
    rows-removed/ (or unchanged/) by the sign of that domain's N, optionally
    under a named subfolder (`under="road-check"`, `under="imagery-review"`). Created if absent."""
    s = rows_sign_sub(footprint_n_cells(product, domain))
    return insert_figures_dir(product, section, f"{under}/{s}" if under else s)


def insert_figures_dir(product: str, section: str | None = None,
                       sub: str | None = None) -> Path:
    """Where insert figures are written, created if it does not exist.

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
    """
    base = insert_scope_dir(product) / "figures"
    if section is not None:
        if section not in INSERT_FIGURE_SECTIONS:
            raise SystemExit(
                f"\n{section!r} is not an insert-figure section. Use one of "
                f"{', '.join(INSERT_FIGURE_SECTIONS)}.\n")
        base = base / section
    if sub is not None:
        allowed = INSERT_FIGURE_SUBFOLDERS.get(section, ())
        if sub not in allowed:
            raise SystemExit(
                f"\n{sub!r} is not a subfolder of {section!r}. "
                f"{'Use one of ' + ', '.join(allowed) if allowed else 'It has none'}.\n")
        base = base / sub
    base.mkdir(parents=True, exist_ok=True)
    return base


def duneline_shift_dir(product: str) -> Path:
    """The duneline-shift directory for a product. Read AND write path."""
    sub = _INSERT_SCOPE.get(product)
    base = product_dir(product)
    # the measurement of N is step 1 of the insert work (moved 2026-09-09)
    return (base / sub / "1-measurement" / "duneline-shift") if sub else (base / "duneline-shift")


# "bridged" is written only by nodata_audit/HAT_bridge_dropouts.py: True where
# an unsurveyed cell was filled by interpolation between measured neighbours.
# It is a THIRD state, not a replacement for "nodata" - a bridged cell is still
# a cell no survey saw, and the nodata mask keeps saying so.
ARRAY_KINDS = ("topography", "dune", "nodata", "bridged")


def array_name(kind: str, gis_id: int | str) -> str:
    """The filename for one domain array. The ONLY place this is spelled."""
    if kind not in ARRAY_KINDS:
        raise SystemExit(
            f"\nunknown array kind {kind!r}. Known: {', '.join(ARRAY_KINDS)}\n")
    return f"domain_{gis_id}_{kind}.npy"


def array_path(kind: str, gis_id: int | str,
               product: str | None = None,
               override: str | None = None) -> Path:
    """Full path to one domain array, directory and name resolved together.

    For the scripts that read a domain at a time - the road placement, the
    setback audit, the units check. Use domain_arrays() when you want the
    whole padded run.
    """
    topo, dune, _ = topo_dirs(product, override)
    parent = dune if kind == "dune" else topo
    return parent / array_name(kind, gis_id)


def domain_arrays(product: str | None = None,
                  override: str | None = None,
                  first_gis: int = 1,
                  last_gis: int = 90,
                  n_buffer: int = 0) -> tuple[list[str], list[str]]:
    """(elevation_paths, dune_paths) for a run, buffer-padded, as strings.

    Replaces the build_domain_file_paths() that was copied verbatim into the
    hindcast runner, the notebook and the groin sweep worker. Those three
    copies each took an `init_year` and pasted it into the name; there is no
    year any more and no name for a caller to build.

    n_buffer is the padding on EACH side. The buffer profiles are shared by
    every period and carry no tag, which is what these arrays now look like
    too.
    """
    topo, dune, _ = topo_dirs(product, override)
    buf_elev = str(BUFFER_DIR / "sample_1_topography.npy")
    buf_dune = str(BUFFER_DIR / "sample_1_dune.npy")

    # A domain outside the surveyed reach (an extended geometry, 2026-09-16)
    # has no array of its own and runs on the buffer profile, exactly as the
    # padding does; what makes it different from padding is its offset.
    from site_layer.hat_extension_domains import SURVEYED_GIS
    lo, hi = SURVEYED_GIS
    elev = [buf_elev] * n_buffer
    dunes = [buf_dune] * n_buffer
    for gis_id in range(first_gis, last_gis + 1):
        if lo <= gis_id <= hi:
            elev.append(str(topo / array_name("topography", gis_id)))
            dunes.append(str(dune / array_name("dune", gis_id)))
        else:
            elev.append(buf_elev)
            dunes.append(buf_dune)
    elev += [buf_elev] * n_buffer
    dunes += [buf_dune] * n_buffer
    return elev, dunes


def versions(product: str) -> list[str]:
    root = dune_topo_root(product)
    return sorted(p.name for p in root.iterdir()
                  if p.is_dir()) if root.is_dir() else []


def require_version(product: str, version: str, what: str = "") -> Path:
    """The version's directory, or a loud exit naming what IS on disk.

    For scripts that carry a version LITERAL instead of resolving through
    CURRENT -- the seaward-row-insert plotters name the layer they draw. On
    2026-09-07 the 1984-start layers v3-v8 were deleted (Hannah's decision:
    keep only unmodified topography), so a literal that was valid the day it
    was written now names nothing. Failing here, before any array is opened,
    says so in one line instead of a FileNotFoundError deep in a plotting loop.
    """
    d = dune_topo_root(product) / version
    if d.is_dir():
        return d
    raise SystemExit(
        f"\n{product}/dune-topo/{version} does not exist"
        + (f" ({what})" if what else "") + ".\n"
        f"versions present: {versions(product) or '(none)'}\n"
        f"The 1984-start layers v3-v8 were deleted 2026-09-07 -- see\n"
        f"  {DOMAIN_ROOT / 'archive_purge_20260907.csv'}\n"
        f"Rebuild one with HAT_insert_seaward_rows.py --dst-version <name>, or "
        f"pass a version that is on disk.\n")


def env_override_name(product: str) -> str:
    """The environment variable that pins ONE product's version.

    Product-scoped on purpose. A bare HAT_TOPO_VERSION would apply to every
    product, and BOTH products currently have a version called "v1" -- so a
    global override set for one period would silently resolve to a real, wrong
    directory for the other instead of failing.
    """
    return "HAT_TOPO_VERSION_" + product.upper().replace("-", "_")


def resolve_version(product: str, override: str | None = None) -> str:
    if override:
        return override

    # THE ENVIRONMENT OUTRANKS THE EXTRACTOR LITERAL, and it has to.
    #
    # Added 2026-09-02, after it cost a run. The crest experiment selected its
    # arm by writing dune-topo/CURRENT, ran, and reported dune-topo\v1 -- the
    # baseline -- because rule 2 below reads the extractor's VERSION literal
    # FIRST and the extractor happened to be sitting on the same product. The
    # CURRENT file was never consulted. Two arms of a three-arm experiment
    # were duplicates of the control, exit code 0, no warning.
    #
    # That is this module's own failure mode, one level up: a caller that
    # cannot say "use THIS version" without editing a source file will end up
    # editing a source file, or will think it said it and be ignored. CURRENT
    # is not usable for that -- it is a persistent, shared default, and a batch
    # that sets it per arm is mutating global state for every other reader.
    env_name = env_override_name(product)
    from_env = os.environ.get(env_name, "").strip()
    if from_env:
        return from_env

    # CURRENT BEFORE THE EXTRACTOR LITERAL (2026-09-04) - see the header. The
    # literal is what the extractor writes; CURRENT is what is read. A product
    # without a CURRENT file behaves exactly as before.
    current = dune_topo_root(product) / "CURRENT"
    if current.is_file():
        name = current.read_text(encoding="utf-8").strip()
        if name:
            return name
    ex_product, ex_version = _extractor_state()
    if ex_version and ex_product == product:
        return ex_version
    avail = versions(product)
    if len(avail) == 1:
        return avail[0]
    raise SystemExit(
        f"\nCannot decide which version of {product!r} to use.\n"
        f"versions present: {avail or '(none)'}\n"
        f"The extractor is currently set to product="
        f"{ex_product!r} version={ex_version!r}, which does not name this "
        f"product.\n"
        f"Fix by one of:\n"
        f"  - pass override=\"<version>\" to topo_dirs()\n"
        f"  - write the version into {current}\n"
        f"  - point the extractor at this product and re-run it\n")


def current_topo_versions() -> dict[str, str]:
    """{product: the dune-topo version it reads TODAY}, for every product.

    What run_registry.rebuild_run_index needs to mark a run superseded: the
    version each product resolves to now (CURRENT outranks the extractor
    literal, see the rules above). A product with no dune-topo tree, or none
    resolvable, is left out rather than guessed.
    """
    out = {}
    for product in PRODUCTS:
        try:
            out[product] = resolve_version(product, None)
        except (SystemExit, FileNotFoundError, OSError, ValueError):
            continue
    return out


def topo_dirs(product: str | None = None,
              override: str | None = None) -> tuple[Path, Path, str]:
    """(topography dir, dunes dir, run name), checked to exist.

    Raises with what IS on disk rather than returning a path that will later
    read as "no arrays for this domain".
    """
    product = product or DEFAULT_PRODUCT
    name = resolve_version(product, override)
    run = dune_topo_root(product) / name
    topo, dune = run / "topography", run / "dunes"

    if not topo.is_dir():
        raise SystemExit(
            f"\nTopography directory does not exist:\n    {topo}\n"
            f"product : {product}\n"
            f"versions present: {versions(product) or '(none)'}\n"
            f"Re-run HAT_dune_topo_extractor.py for that product and version, "
            f"or pass a different override to topo_dirs().\n")
    return topo, dune, name


def run_name(product: str | None = None) -> str:
    """Kept for callers that only want the label."""
    return topo_dirs(product)[2]
