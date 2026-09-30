"""
Where are the shared map layers, and where does the house style live?

    from site_layer.hat_map_layers import ISLAND_OUTLINE, NC_COAST, DOMAIN_BOXES

The map layers in data/hatteras_init/map_elements/; the domain boxes are
re-exported from 5-scr/2-transect-frame/, never copied. Details: scripts/site_layer/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-22
"""

from __future__ import annotations

from pathlib import Path

if __package__ in (None, ""):
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from site_layer.hat_observed_rates import DOMAIN_BOXES  # noqa: E402,F401

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents
                    if (_p / "pyproject.toml").exists())
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"

MAP_ELEMENTS = INIT_ROOT / "map_elements"
# Hatteras Island only: the NC 1:80k land polygons it overlaps, merged
ISLAND_OUTLINE = MAP_ELEMENTS / "hatteras_outline" / "HAT_island_outline.shp"
# The coast around it (sound shores, Ocracoke, Pea Island) for maps reaching past the island
NC_COAST = MAP_ELEMENTS / "nc_coast_80k" / "nc_80k_hatteras_window.geojson"
NC_COAST_SOURCE = Path("D:/Hatteras_GIS/Outlines/nc_80k/nc_80k.shp")
NC_COAST_PAD_M = 25_000.0
# Natural Earth 1:10m states for the regional locator inset.
NE_STATES = MAP_ELEMENTS / "natural_earth" / "ne_10m_states_southeast_us.geojson"
ARCHIVE = MAP_ELEMENTS / "archive"

STYLE_DOC = PROJECT_ROOT / "scripts" / "figure_making" / "STYLE.md"
STYLE_SHEET_DIR = PROJECT_ROOT / "output" / "figures" / "style"


if __name__ == "__main__":
    for name in ("MAP_ELEMENTS", "ISLAND_OUTLINE", "NC_COAST", "NE_STATES",
                 "DOMAIN_BOXES", "ARCHIVE", "STYLE_DOC", "STYLE_SHEET_DIR"):
        path = globals()[name]
        print(f"{'ok' if path.exists() else 'MISSING':8} {name:16} "
              f"{path.relative_to(PROJECT_ROOT).as_posix()}")
