"""
Which Barrier3D domains does a script read: which product, and which version within it?

    from site_layer.hat_topo_version import topo_dirs, domain_arrays
    TOPO_DIR, DUNE_DIR, VERSION = topo_dirs("2009-start")

Resolved once (environment, then CURRENT, then the extractor literal), with the
road, dune-line and offset paths beside it; a name not on disk raises. Details: scripts/site_layer/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-19
"""

from __future__ import annotations

import os
import re
from pathlib import Path

_HERE = Path(__file__).resolve()
# data/ and output/ live in the workspace on D:, not in the repo (hatteras_paths)
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
DOMAIN_ROOT = INIT_ROOT / "1-barrier3d-domains"

# Inputs shared by every product, under 1-barrier3d-domains/shared/ since 2026-10-19
SHARED_DIR = DOMAIN_ROOT / "shared"
# The buffer domains either side of the 90 real ones
BUFFER_DIR = SHARED_DIR / "buffer"
DOMAIN_CLIPS_DIR = SHARED_DIR / "domain-clips-1m"
CONTROL_PICKS_DIR = SHARED_DIR / "control-picks"
UNFILLED_2009_DIR = SHARED_DIR / "npy-arrays_2009_unfilled"
DOMAIN_GEOJSON_DIR = SHARED_DIR / "domain-geojson"
# Figures that compare products, so belong to neither
DOMAIN_FIGURES_DIR = DOMAIN_ROOT / "figures"


def domain_clip_file(domain: int, kind: str = "resampled") -> Path:
    """One domain's DEM clip: kind "resampled" (10 m) or "clip" (1 m)."""
    d = int(domain)
    if kind not in ("resampled", "clip"):
        raise ValueError(f"kind must be 'resampled' or 'clip', not {kind!r}")
    return DOMAIN_CLIPS_DIR / f"domain_{d}" / f"{kind}_domain_{d}.tif"

# What topo_dirs() resolves when no product is named
DEFAULT_PRODUCT = "2009-start"

# Named for the period that reads them since 2026-10-19; the DEMs and versions are unchanged
PRODUCTS = ("1996-start", "2009-start", "forecast")
# The names before 2026-10-19, still in run metadata and old scripts; resolved to the new ones
LEGACY_PRODUCTS = {"1984-start": "1996-start", "2004-start": "2009-start"}

# Which product a hindcast period reads, by start year: the single definition
YEAR_PRODUCT = {
    # First, so year_for_product() answers 1996 and 2009
    1996: "1996-start",
    2009: "2009-start",  # was key 2010 until 2026-10-05; the period starts in its DEM year
    # The retired 1984 and 2004 periods read the same DEMs; kept until their config entries go
    1984: "1996-start",
    2004: "2009-start",
}


def canonical_product(product: str) -> str:
    """The product name on disk today, translating a pre-2026-10-19 name."""
    return LEGACY_PRODUCTS.get(product, product)

# Which NC-12 line a period's road is measured from, by true vintage: the single definition
ROAD_LINE_VINTAGES = (1978, 2008)
ROAD_LINE_FOR_YEAR = {
    1984: 1978,
    1996: 1978,
    2004: 2008,
    2009: 2008,
}

# Which setback folders hold a measurement and which a derivation
ROAD_SETBACK_KIND = {
    1984: "measured",
    2004: "measured",
    1996: "derived",
    2009: "derived",
}

MGMT_ROOT = INIT_ROOT / "4-mgmt-forcing"
ROADS_ROOT = MGMT_ROOT / "road_offset"

# The rest of 4-mgmt-forcing
ROAD_LINE_ROOT = ROADS_ROOT / "raw_offset"            # <vintage>/nc12_<vintage>.geojson
# Today's NC-12 (NCDOT inventory), for the dune-line position figures; no hindcast reads it
ROAD_LINE_CURRENT = ROAD_LINE_ROOT / "current" / "nc12_current.geojson"
ROAD_RASTER_ROOT = ROADS_ROOT / "raster"              # <vintage>/masks/
ROAD_SETBACK_ROOT = ROADS_ROOT / "dunestart_offset"   # measured/ derived/ modifications/
ROAD_ARCHIVE = ROADS_ROOT / "archive"
# The old-method setbacks, kept for the method comparison
LEGACY_SETBACK_ROOT = ROAD_ARCHIVE / "superseded_20260911"
# The 1984 dune-start setbacks as measured on 1984-start/v1
SETBACK_1984_V1_DIR = ROAD_ARCHIVE / "superseded_20260907" / "1984"

ROAD_ELEVATION_DIR = MGMT_ROOT / "road_elevation"
# ONE file for every period; see HATTERAS_ROAD_ELEVATION_FILE for why.
ROAD_ELEVATION_FILE = ROAD_ELEVATION_DIR / "RoadElevation.csv"
ROAD_RELOCATION_ROOT = MGMT_ROOT / "road_relocation"
NOURISHMENT_DIR = MGMT_ROOT / "nourishment"            # figures only
# The source record behind HATTERAS_NOURISHMENT_PROJECTS.
MGMT_RECORD_XLSX = MGMT_ROOT / "Hatteras_Management_Timelines.xlsx"


def init_relpath(path: Path) -> str:
    """A path under INIT_ROOT as the POSIX string hatteras_site_config carries."""
    return Path(path).relative_to(INIT_ROOT).as_posix()


def legacy_setback_file(year: int) -> Path:
    """The OLD-METHOD setback for a measured start (1984, 2004) -- superseded,
    kept only to compare the two methods."""
    return LEGACY_SETBACK_ROOT / str(int(year)) / f"RoadSetback_{int(year)}.csv"


def road_relocation_dir(vintage_from: int, vintage_to: int) -> Path:
    """The measured displacement between two LINE vintages."""
    a, b = _check_line_vintage(vintage_from), _check_line_vintage(vintage_to)
    return ROAD_RELOCATION_ROOT / f"{a}_{b}"


def road_relocation_file(vintage_from: int, vintage_to: int) -> Path:
    a, b = _check_line_vintage(vintage_from), _check_line_vintage(vintage_to)
    return road_relocation_dir(a, b) / f"road_relocation_{a}_{b}.csv"


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
    return ROAD_LINE_ROOT / str(v) / f"nc12_{v}.geojson"


def road_mask_dir(vintage: int) -> Path:
    """Where HAT_rasterize_road_to_domains.py put one line vintage's masks."""
    v = _check_line_vintage(vintage)
    return ROAD_RASTER_ROOT / str(v) / "masks"


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
    return ROAD_SETBACK_ROOT / kind / str(int(year))


def road_setback_file(year: int) -> Path:
    """The model-facing 2-row RoadSetback_<year>_dunestart.csv for a period."""
    return road_setback_dir(year) / f"RoadSetback_{int(year)}_dunestart.csv"


def road_setback_relpath(year: int) -> str:
    """road_setback_file() relative to INIT_ROOT, POSIX-style -- the form
    HATTERAS_PERIODS[year]["road_setback_file"] carries."""
    return road_setback_file(year).relative_to(INIT_ROOT).as_posix()


# The dune lines, paired with a period year by imagery vintage: the only place, no copies
BRIE_ROOT = INIT_ROOT / "2-brie-offset"
RAW_OFFSET_DIR = BRIE_ROOT / "raw_offsets"
DUNELINE_DIR = BRIE_ROOT / "dunelines"
# The rest of 2-brie-offset
RAW_OFFSET_EXT_DIR = RAW_OFFSET_DIR / "ext"
RAW_OFFSET_SUPERSEDED = RAW_OFFSET_DIR / "superseded_20260915_gis_exports"
TRANSECT_DIR = BRIE_ROOT / "transects"
TRANSECT_FILE_100M = TRANSECT_DIR / "transects_100m.geojson"
TRANSECT_EXT_TABLE = TRANSECT_DIR / "transects_100m_ext.csv"

# Which feature the offset was measured from: the dune build keeps <year>/v<n>, others nest by source
OFFSET_SOURCES = ("duneline", "shoreline")
DEFAULT_OFFSET_SOURCE = "duneline"
# What a run reads unless HAT_ISLAND_OFFSET_SOURCE says otherwise: the DEM-centred shoreline since 2026-10-05
# (matrix runs before then are dune line). DEFAULT_OFFSET_SOURCE stays the build-naming default
RUN_OFFSET_SOURCE = "shoreline"

# A build's three files share one stem, named for the feature measured from
_OFFSET_STEMS = {
    "duneline": "Island_Dune_Offsets",
    "shoreline": "Island_Shoreline_Offsets",
}
_OFFSET_FILE_SUFFIXES = {
    "padded": "_PADDED_{n}.csv",            # what the model reads
    "input": "_CASCADE_Input.csv",          # 90 domains, one column
    "unpadded": "_CASCADE_Input_unpadded.csv",  # 90, with domain ids
}


def offset_basename(year: int, source: str = DEFAULT_OFFSET_SOURCE) -> str:
    """The stem a build's three files share, e.g. Island_Dune_Offsets_1996.
    step2_raw_offsets_to_island_offset.py names its outputs from this rather than keeping a
    second spelling of them (2026-09-22)."""
    return f"{_OFFSET_STEMS[_check_offset_source(source)]}_{int(year)}"


def _check_offset_source(source: str) -> str:
    s = str(source)
    if s not in OFFSET_SOURCES:
        raise ValueError(
            f"unknown island-offset source {source!r}; "
            f"known: {', '.join(OFFSET_SOURCES)}")
    return s


def offset_start_dir(year: int, source: str = DEFAULT_OFFSET_SOURCE) -> Path:
    """One source's folder for one period start: PROVENANCE.md, CURRENT,
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
    """
    return BRIE_ROOT / str(int(year)) / _check_offset_source(source)


def offset_year_dir(year: int) -> Path:
    """The period start's folder itself, which now holds only source folders
    and comparisons/ -- no build of its own."""
    return BRIE_ROOT / str(int(year))


def offset_version(year: int, source: str = DEFAULT_OFFSET_SOURCE):
    """Which build of a start's island offset every reader takes.

    Order, mirroring topo_dirs(): HAT_OFFSET_VERSION_<year> in the
    environment, then the CURRENT file, then the only v<n> directory present.
    None for the flat layout (no v<n> directory at all). Several v<n> and no
    CURRENT is an error, not a guess; so is naming a version that is not on
    disk. Moved here from hatteras_site_config._island_offset_file on
    2026-09-18 so the figure scripts that resolved it themselves share it.

    A non-default source reads HAT_OFFSET_VERSION_<year>_<SOURCE> instead, so
    overriding the shoreline arm cannot silently move the dune build the
    runner reads (2026-09-22).
    """
    y = int(year)
    s = _check_offset_source(source)
    base = f"2-brie-offset/{y}/{s}"
    # The default source keeps its plain env key; a non-default source gets its own
    env_key = (f"HAT_OFFSET_VERSION_{y}" if s == DEFAULT_OFFSET_SOURCE
               else f"HAT_OFFSET_VERSION_{y}_{s.upper()}")
    d = offset_start_dir(y, s)
    versions = sorted(p.name for p in d.iterdir()
                      if p.is_dir() and re.fullmatch(r"v\d+", p.name)) if d.is_dir() else []
    env = os.environ.get(env_key)
    current = d / "CURRENT"
    if env:
        version = env.strip()
    elif current.is_file():
        version = current.read_text(encoding="utf-8").strip()
    elif len(versions) == 1:
        version = versions[0]
    elif versions:
        raise RuntimeError(
            f"{base}/ holds {versions} and no CURRENT file; write one, or set "
            f"{env_key}.")
    else:
        return None
    if not (d / version).is_dir():
        raise FileNotFoundError(
            f"{base}/{version}/ does not exist (have {versions or 'no versions'}); "
            f"check CURRENT or {env_key}.")
    return version


def offset_build_dir(year: int, version: str | None = None,
                     source: str = DEFAULT_OFFSET_SOURCE) -> Path:
    """The folder of one build: the resolved version unless one is given."""
    s = _check_offset_source(source)
    v = version or offset_version(year, s)
    base = offset_start_dir(year, s)
    return base / v if v else base


def offset_comparison_dir(year: int, name: str) -> Path:
    """<year>/comparisons/<name>/ -- where a comparison BETWEEN builds lands.

    A comparison is neither a version nor a source, so it does not belong in
    either's folder. Until 2026-09-22 a version comparison was written into
    the later version's own folder, which made a build's folder hold both the
    build and a judgement about it; source comparisons are filed here from the
    start. `name` says what was compared, e.g. "duneline_vs_shoreline".
    """
    # offset_YEAR_dir, not offset_start_dir: a comparison of two sources belongs to neither
    return offset_year_dir(year) / "comparisons" / name


def offset_source_comparison_dir(year: int, name: str,
                                 shoreline_version: str | None = None) -> Path:
    """Where a dune line vs shoreline comparison lands.

    Against the CURRENT shoreline build: <year>/comparisons/<name>/, beside
    both sources, since 2026-10-06 (Hannah: the shoreline and dune lines are
    final, so the comparisons can sit at the year). Against any other
    shoreline build: <year>/shoreline/<v>/comparisons/<name>/, filed with that
    build as it has been since 2026-09-29, so a superseded comparison still
    says which build it was drawn against. The builds used are written into
    each comparison's README and caption.
    """
    v = shoreline_version or offset_version(year, "shoreline")
    if v == offset_version(year, "shoreline"):
        return offset_comparison_dir(year, name)
    return offset_build_dir(year, v, "shoreline") / "comparisons" / name


def offset_file(year: int, kind: str = "padded", total_domains: int = 120,
                version: str | None = None,
                source: str = DEFAULT_OFFSET_SOURCE) -> Path:
    """One file of a start's build. kind: padded (the model input), input
    (90 values), or unpadded (90, with domain ids)."""
    s = _check_offset_source(source)
    name = offset_basename(year, s) + _OFFSET_FILE_SUFFIXES[kind].format(
        n=int(total_domains))
    return offset_build_dir(year, version, s) / name
DUNE_LINE_FOR_YEAR = {
    1984: 1984,
    1996: 1997,   # no 1996 imagery; the nearest island-wide survey
    2004: 2004,
    2009: 2009,   # the 2009 line; key was 2010 until 2026-10-05 (no 2010 aerial imagery)
    2024: 2023,   # the 2023 NOAA imagery (D:\hatteras_ch1\gis_sources\Aerial\2023); end year of 2004-2024 and 2010-2024
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


# The CoastSat shoreline window per period: the DEM's calendar year (offset shoreline/v4, min 5
# positions), since 2026-10-12, with the storms starting the year after (advisor's proposal).
# Before: +/-6 months of the DEM's lidar flights (v3): 1996-04-12..1997-04-12, 2009-02-17..2010-02-17;
# +/-1 yr (v2): 1995-10-12..1997-10-12, 2008-08-17..2010-08-17
SHORELINE_WINDOW_FOR_YEAR = {
    1996: ("1996-01-01", "1996-12-31"),
    2009: ("2009-01-01", "2009-12-31"),
}


def shoreline_window_for_year(year, strict: bool = True):
    """The CoastSat averaging window a period start reads, as (start, end)."""
    try:
        return SHORELINE_WINDOW_FOR_YEAR[int(year)]
    except (KeyError, TypeError, ValueError):
        if not strict:
            return None
        known = ", ".join(f"{y} -> {a}-{b}" for y, (a, b)
                          in sorted(SHORELINE_WINDOW_FOR_YEAR.items()))
        raise SystemExit(
            f"\nno CoastSat shoreline window known for period year {year!r}. "
            f"Known: {known}\n"
            f"Add it to SHORELINE_WINDOW_FOR_YEAR in {__file__}.\n")


def shoreline_raw_file(window) -> Path:
    """raw_offsets/<start>_<end>_shoreline_offset_raw.csv -- the per-transect
    stations of one averaging window's mean shoreline, the shoreline
    counterpart of dune_raw_file(). Written by step1_line_to_raw_offsets.py from
    the mean-shoreline geojson, which hat_observed_rates owns."""
    # A window is two calendar years or two ISO dates, named as the mean_shoreline folder is
    a, b = (str(v) if "-" in str(v) else str(int(v)) for v in window)
    return RAW_OFFSET_DIR / f"{a}_{b}_shoreline_offset_raw.csv"


def shoreline_raw_file_for_year(year, strict: bool = True):
    """The raw file a period year reads, through SHORELINE_WINDOW_FOR_YEAR."""
    w = shoreline_window_for_year(year, strict=strict)
    return None if w is None else shoreline_raw_file(w)


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
    product = canonical_product(product)
    for year, prod in YEAR_PRODUCT.items():
        if prod == product:
            return int(year)
    if not strict:
        return None
    raise SystemExit(
        f"\nno period year known for topography product {product!r}. "
        f"Known: {', '.join(sorted(YEAR_PRODUCT.values()))}\n"
        f"Add it to YEAR_PRODUCT in {__file__}.\n")

# The arrays carry no year tag: the period lives in the directory; use domain_arrays() or array_path()

# The extractor with ALONGSHORE_FLIP = True (the other copies are unflipped)
from hatteras_paths import HATTERAS_DIR  # noqa: E402  (scripts live in the clone, not the workspace)
EXTRACTOR = (HATTERAS_DIR / "input_prep" / "1-barrier3d-domains" / "1-extraction"
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
    product = canonical_product(product)
    if product not in PRODUCTS:
        raise SystemExit(
            f"\nunknown topography product {product!r}. Known: "
            f"{', '.join(PRODUCTS)}\n")
    return DOMAIN_ROOT / product


def dune_topo_root(product: str) -> Path:
    return product_dir(product) / "dune-topo"


# The extraction half of a product; dune-topo/ sits beside it at the product root
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


# "bridged" (HAT_bridge_dropouts.py): a third state beside "nodata", not a replacement
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

    # A domain outside the surveyed reach runs on the buffer profile; only its offset differs
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


def env_override_name(product: str) -> str:
    """The environment variable that pins ONE product's version.

    Product-scoped on purpose. A bare HAT_TOPO_VERSION would apply to every
    product, and BOTH products currently have a version called "v1" -- so a
    global override set for one period would silently resolve to a real, wrong
    directory for the other instead of failing.
    """
    return "HAT_TOPO_VERSION_" + canonical_product(product).upper().replace("-", "_")


def _legacy_env_names(product: str) -> list[str]:
    """The pre-2026-10-19 variable names for a product (HAT_TOPO_VERSION_1984_START),
    still honoured so an old pin is not silently dropped."""
    new = canonical_product(product)
    return ["HAT_TOPO_VERSION_" + old.upper().replace("-", "_")
            for old, cur in LEGACY_PRODUCTS.items() if cur == new]


def resolve_version(product: str, override: str | None = None) -> str:
    if override:
        return override

    # The environment outranks CURRENT and the extractor literal, so a caller can pick a version
    for env_name in (env_override_name(product), *_legacy_env_names(product)):
        from_env = os.environ.get(env_name, "").strip()
        if from_env:
            return from_env

    # CURRENT before the extractor literal; without a CURRENT file, as before
    current = dune_topo_root(product) / "CURRENT"
    if current.is_file():
        name = current.read_text(encoding="utf-8").strip()
        if name:
            return name
    ex_product, ex_version = _extractor_state()
    if ex_version and canonical_product(ex_product or "") == canonical_product(product):
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

    What run_registry.scan_runs needs to mark a run superseded: the
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
    # Runs made before 2026-10-19 record the old names; answer for those too
    for old, new in LEGACY_PRODUCTS.items():
        if new in out:
            out[old] = out[new]
    return out


def topo_dirs(product: str | None = None,
              override: str | None = None) -> tuple[Path, Path, str]:
    """(topography dir, dunes dir, run name), checked to exist.

    Raises with what IS on disk rather than returning a path that will later
    read as "no arrays for this domain".
    """
    product = canonical_product(product or DEFAULT_PRODUCT)
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
