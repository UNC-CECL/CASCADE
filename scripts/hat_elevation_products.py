# ==============================================================================
# hat_elevation_products.py
#
# Which elevation product does a script read, and where does it live?
#
# WHY THIS EXISTS
#   Six scripts used to build these paths by hand, each pasting its own copy of
#       .../0-elevation/{1-gapfill-1m,2-resampled-10m}/<TAG>/
#   and that is how HAT_road_elevation.py broke. It sets
#   FILL_SOURCE = "2008_NOAA_IOCM" and builds GAPFILL_1M_ROOT / FILL_SOURCE.
#   When 2008 was moved under superseded/ on 2026-08-25 the path stopped
#   resolving - and NOTHING RAISED. The globs simply returned nothing, the
#   script carried on, and every domain reported "no fill available".
#
#   Same failure as the one scripts/hat_topo_version.py was written for: a
#   layout change that a hand-built path absorbs silently. Same fix. The
#   product is resolved ONCE, here, and a name that does not exist on disk is
#   an immediate, loud error listing what does.
#
# THE LAYOUT IT RESOLVES  (product first, stage second - 2026-08-25)
#
#     data/hatteras_init/0-elevation/
#         2009-2014/                  the baseline DEM
#             1-gapfill-1m/           gapfill_audit.csv + clip_domain_*.tif
#             2-resampled-10m/        resample_audit.csv + resampled_domain_*.tif
#             figures/
#         2009-2014-1996/             the 1984-start DEM
#             1-gapfill-1m/           mosaic_1984_audit.csv + clip_domain_*.tif
#             2-resampled-10m/
#             figures/
#         superseded/<attempt>/       same shape, not for use
#         source-selection/           island-wide, belongs to no product
#         FIGURES.md                  figure design decisions, shared
#
#   BEFORE, it was stage first: 1-gapfill-1m/<TAG>/ and 2-resampled-10m/<TAG>/,
#   with every product's figures pooled in one figures/. The names were the
#   FILL SOURCE ("2014_NOAA_PostSandy"), which said what was added but not what
#   the product contained. Product folders are now named for their COMPOSITION,
#   deliberately not for a hindcast period: the 2009-2014 DEM currently serves
#   both the 1984 and 2004 periods, so naming it "2004-start" would assert
#   something that is not true.
#
# USAGE
#     import sys; sys.path.insert(0, <repo>/scripts)
#     from hat_elevation_products import product, PRODUCTS
#
#     p = product("2009-2014")
#     p.gapfill_1m / "clip_domain_7_filled.tif"
#     p.resampled_10m, p.figures, p.audit_1m
#
#   product() checks the directory exists and raises with the available names
#   if it does not. Pass check=False only when creating the product for the
#   first time, which is what the two producer scripts do.
# ==============================================================================

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

_HERE = Path(__file__).resolve()
PROJECT_ROOT = _HERE.parents[1]          # scripts/ -> repo root
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
ELEVATION_ROOT = INIT_ROOT / "0-elevation"

GAPFILL_1M = "1-gapfill-1m"
RESAMPLED_10M = "2-resampled-10m"
FIGURES = "figures"

SUPERSEDED_DIR = "superseded"
SOURCE_SELECTION_DIR = "source-selection"


@dataclass(frozen=True)
class Product:
    """One elevation product: a composition of surveys, and where it lives."""
    name: str
    surveys: tuple[str, ...]
    audit_name: str
    built_by: str
    summary: str
    superseded: bool = False

    @property
    def root(self) -> Path:
        base = ELEVATION_ROOT / SUPERSEDED_DIR if self.superseded else ELEVATION_ROOT
        return base / self.name

    @property
    def gapfill_1m(self) -> Path:
        return self.root / GAPFILL_1M

    @property
    def resampled_10m(self) -> Path:
        return self.root / RESAMPLED_10M

    @property
    def figures(self) -> Path:
        return self.root / FIGURES

    @property
    def audit_1m(self) -> Path:
        return self.gapfill_1m / self.audit_name

    @property
    def audit_10m(self) -> Path:
        return self.resampled_10m / "resample_audit.csv"


# The survey codes a product's clip_domain_*_survey.tif may carry, MOST
# SPECIFIC FIRST. This is the precedence HAT_dem_resample_clip.downsample_survey
# resolves a mixed 2 x 2 block with, so the order is a decision. 2009 is the
# base and 0 is "no survey saw it"; neither appears here.
FILL_CODES = {
    "2009-2014": (2014,),
    "2009-2014-1996": (1996, 2014),
}

PRODUCTS: dict[str, Product] = {
    "2009-2014": Product(
        name="2009-2014",
        surveys=("2009 USACE (base)", "2014 NOAA Post-Sandy (gap fill)"),
        audit_name="gapfill_audit.csv",
        built_by="scripts/input_prep/0-elevation/2-produce/HAT_dem_gap_fill.py",
        summary="The baseline DEM. 2009 USACE with its nodata filled from the "
                "2014 NOAA Post-Sandy survey. Fills holes only - no measured "
                "2009 cell is ever changed.",
    ),
    "2009-2014-1996": Product(
        name="2009-2014-1996",
        surveys=("2009 USACE (base)", "2014 NOAA Post-Sandy (gap fill)",
                 "1996 NOAA/NASA ALACE (override, no road boundary)"),
        audit_name="mosaic_1984_audit.csv",
        built_by="scripts/input_prep/0-elevation/2-produce/HAT_dem_1984_mosaic.py",
        summary="The 1984-start DEM. As 2009-2014, plus the 1996 ALACE survey "
                "OVERWRITING measured ground wherever ALACE has data. No road "
                "boundary since 2026-08-26 - the landward limit is the ALACE "
                "swath edge, so the graft seam lands at the dune toe.",
    ),
}

# Nothing is superseded on disk right now. The 2008 NOAA IOCM attempt was
# registered here until 2026-08-26 with the note "kept so its comparison
# figures can be regenerated" - but its product folder had already been
# deleted, so product("2008_NOAA_IOCM") raised rather than resolving. The
# point-cloud path that built it is gone from HAT_dem_gap_fill.py too, so it
# is not reproducible from this repo either. The machinery below stays: a
# future superseded product registers here with superseded=True and lands
# under 0-elevation/superseded/<name>/.
SUPERSEDED: dict[str, Product] = {}


def _known() -> str:
    live = "\n".join(f"    {n:<18} {p.summary.split('.')[0]}."
                     for n, p in PRODUCTS.items())
    old = "\n".join(f"    {n:<18} (superseded)" for n in SUPERSEDED)
    return f"{live}\n{old}" if old else live


def product(name: str, check: bool = True) -> Product:
    """
    Resolve a product by name.

    Raises rather than returning a path that is not there. That is the whole
    point of this module - see the note on HAT_road_elevation.py above.
    Pass check=False when the product is about to be CREATED.
    """
    p = PRODUCTS.get(name) or SUPERSEDED.get(name)
    if p is None:
        raise SystemExit(
            f"\nunknown elevation product {name!r}. Known:\n{_known()}\n")
    if check and not p.root.is_dir():
        raise SystemExit(
            f"\nelevation product {name!r} is not on disk:\n"
            f"    {p.root}\n"
            f"Build it with:\n    {p.built_by}\n"
            f"then HAT_dem_resample_clip.py --product {name}\n")
    return p


def fill_codes(name: str) -> tuple[int, ...]:
    """Survey codes this product's rasters may carry, most specific first."""
    return FILL_CODES.get(name, ())


def source_selection_dir() -> Path:
    """Island-wide DEM-candidate scoring - belongs to no single product."""
    return ELEVATION_ROOT / SOURCE_SELECTION_DIR
