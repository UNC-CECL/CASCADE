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
#   2. the extractor's VERSION, but only if the extractor is currently pointed
#      at the SAME product - so "bump VERSION and the whole tree follows" still
#      holds for the product you are actively working on
#   3. a CURRENT file in the product's dune-topo/ directory
#   4. the only version present, if there is exactly one
#   otherwise: raise, listing what is on disk
#
# USAGE
#     from hat_topo_version import topo_dirs
#     TOPO_DIR, DUNE_DIR, RUN_NAME = topo_dirs()                  # 2004-start
#     TOPO_DIR, DUNE_DIR, RUN_NAME = topo_dirs("1984-start")
#     TOPO_DIR, DUNE_DIR, RUN_NAME = topo_dirs("2004-start", override="v1")
# ==============================================================================

from __future__ import annotations

import re
from pathlib import Path

_HERE = Path(__file__).resolve()
PROJECT_ROOT = _HERE.parents[1]
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
}


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
EXTRACTOR = (PROJECT_ROOT / "scripts" / "input_prep" / "1-barrier3d-domains"
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


def npy_dirs(product: str) -> tuple[Path, Path]:
    """(elevation arrays, survey/provenance arrays) the extractor reads."""
    d = product_dir(product)
    return d / "npy-arrays", d / "npy-arrays_survey"


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

    elev = [buf_elev] * n_buffer
    dunes = [buf_dune] * n_buffer
    for gis_id in range(first_gis, last_gis + 1):
        elev.append(str(topo / array_name("topography", gis_id)))
        dunes.append(str(dune / array_name("dune", gis_id)))
    elev += [buf_elev] * n_buffer
    dunes += [buf_dune] * n_buffer
    return elev, dunes


def versions(product: str) -> list[str]:
    root = dune_topo_root(product)
    return sorted(p.name for p in root.iterdir()
                  if p.is_dir()) if root.is_dir() else []


def resolve_version(product: str, override: str | None = None) -> str:
    if override:
        return override
    ex_product, ex_version = _extractor_state()
    if ex_version and ex_product == product:
        return ex_version
    current = dune_topo_root(product) / "CURRENT"
    if current.is_file():
        name = current.read_text(encoding="utf-8").strip()
        if name:
            return name
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
