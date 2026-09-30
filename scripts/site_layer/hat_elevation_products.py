"""
Which elevation product does a script read, and where does it live?

    from site_layer.hat_elevation_products import product
    p = product("2009-2014"); p.gapfill_1m, p.resampled_10m, p.figures

Resolves data/hatteras_init/0-elevation/<product>/<stage>/ once; a product not
on disk raises, listing the ones that are. Details: scripts/site_layer/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(
    _p for _p in _HERE.parents
    if (_p / "pyproject.toml").exists())          # scripts/ -> repo root
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


# Survey codes a clip_domain_*_survey.tif may carry, most specific first (the resample precedence)
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

# Nothing is superseded on disk now; a future superseded product registers here
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


def duneline_check_dir(name: str) -> Path:
    """The dune-line-vs-DEM check for one product: <product>-duneline/.

    NOT a product, though it sits beside them: HAT_dem_duneline_coverage.py
    and HAT_plot_duneline_offset.py measure where the digitised dune lines
    fall on that product's surface, and write here (named 2026-09-18; the two
    producers and the row-insert report each spelled the folder themselves).
    """
    if name not in PRODUCTS and name not in SUPERSEDED:
        raise SystemExit(
            f"\nunknown elevation product {name!r}. Known:\n{_known()}\n")
    return ELEVATION_ROOT / f"{name}-duneline"
