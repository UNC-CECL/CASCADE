# ==============================================================================
# hat_map_layers.py
#
# WHERE ARE THE SHARED MAP LAYERS, AND WHERE DOES THE HOUSE STYLE LIVE?
#
# WHY THIS EXISTS
#   Until 2026-09-18 these sat in data/hatteras_init/9-figures/, numbered like
#   a model-input stage although it was not one, and mixing three kinds of
#   thing: GIS layers (data), the written style (documentation for whoever
#   writes a figure script) and a rendered style sheet (a figure). They are
#   split by kind now:
#
#     data/hatteras_init/map_elements/     the map layers, un-numbered
#         hatteras_outline/                the island, NC 1:80k (subset)
#         nc_coast_80k/                    NC 1:80k clipped to the study area
#         natural_earth/                   the locator-map states
#         archive/                         the retired 1000 m domain polygons
#     scripts/figure_making/STYLE.md       the written house style
#     output/figures/style/                the rendered style sheet
#
#   Two of the layers were read straight off the D: GIS drive, so the figures
#   that drew them could not be made without it: the domain boxes
#   (D:/Hatteras_GIS/domains.geojson, identical in geometry and every
#   attribute to 5-scr's HAT_domains.json, which is what DOMAIN_BOXES points at)
#   and the NC 1:80k coastline (D:/Hatteras_GIS/Outlines/nc_80k/, 6 MB for the
#   whole state; NC_COAST is the 25 km window around the domains, rebuilt by
#   scripts/figure_making/tools/clip_nc_coast.py).
#
# THE DOMAIN BOXES ARE NOT A MAP LAYER TO COPY. They are the model's frame and
# are owned by 5-scr/2-transect-frame/; DOMAIN_BOXES re-exports that path so a
# figure script needs one import, not a second copy that could drift.
# ==============================================================================

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
# Hatteras Island only: the 208 NC 1:80k land polygons it overlaps, merged
# into 52 (checked 2026-09-18: 88.85 km2 either way, zero difference).
ISLAND_OUTLINE = MAP_ELEMENTS / "hatteras_outline" / "HAT_island_outline.shp"
# The coast around it -- sound shores, Ocracoke, Pea Island -- for maps whose
# window reaches past the island (the overwash maps pad 7.2 km west).
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
