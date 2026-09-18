# ==============================================================================
# hat_extension_domains.py
#
# THE ALONGSHORE REACH A RUN MODELS, BY NAME -- and the 500 m bins that give
# the coast beyond the 90 surveyed domains a domain number.
#
# WHY THIS EXISTS (2026-09-16, the Pea Island extension experiment)
#   The hindcast models GIS 1-90, Cape Point to just north of Rodanthe, and
#   pads each end with 15 invented buffer domains that extrapolate the local
#   shoreline slope and then bridge back to close BRIE's periodic ring. The
#   edgeBE preset pins GIS 1 and 90 to their observed rates with source/sink
#   values that absorb whatever the buffer gets wrong.
#
#   Hannah's question: what if the buffer carried the REAL coast's orientation
#   instead -- Pea Island north to near Oregon Inlet, and the last kilometre
#   south to Cape Point -- with the end domains re-solved at the new ends?
#   The dune lines (every vintage) and the CoastSat record both reach Oregon
#   Inlet, so the coast is measured; only its domain numbers were missing.
#
# THE NUMBERING
#   Hannah's whole-island domain polygons
#   (1-barrier3d-domains/domain-geojson/domains_pea_hatteras_120.geojson,
#   ID 1-121, EPSG:3725, 2000 x 500 m each) continue the surveyed numbering
#   north to Oregon Inlet; ID 1-90 are the surveyed polygons vertex for
#   vertex. An extension domain is numbered exactly as a surveyed one was: a
#   spatial join onto its polygon (join_lines, join_origins below). Nothing
#   about GIS 1-90 changes. An earlier numbering by northing bin (a 502.56 m
#   grid continued from GIS 90, the same evening) was removed once the
#   polygons existed; it had placed 4 dune-line and 40 CoastSat transects one
#   domain off, because the drawn polygons are not on a regular grid.
#
# THE GEOMETRIES
#   "base"    GIS 1-90, the production reach. Every matrix run.
#   "n115"    GIS 1-115: 25 domains of Pea Island added, stopping ~3 km short
#             of Oregon Inlet where the rates turn inlet-dominated (Hannah,
#             2026-09-16: "lets do 115"). GIS 1 stays the southern end. No
#             southern extension: no polygon lies south of GIS 1, and the
#             1997 dune line ends 440 m south of it in any case.
#
#   HAT_GEOMETRY in the environment selects one; hatteras_site_config builds
#   HATTERAS_DOMAINS from it. Extension domains carry the shared buffer
#   topography (hat_topo_version.domain_arrays), the measured dune-line offset
#   (2-brie-offset/<year>/ext/<geometry>/), zero background erosion unless
#   HAT_BE_OVERRIDE names them, and no management of any kind.
#
# USAGE
#     from site_layer.hat_extension_domains import GEOMETRIES, gis_bounds, join_lines
#     first, last = gis_bounds("n115")          # (1, 115)
#     join_lines(transects_gdf)                 # -> domain per transect line
# ==============================================================================
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents
                    if (_p / "pyproject.toml").exists())
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"

DOMAIN_CRS = "EPSG:3725"

# The surveyed reach. Topography, management inputs and the committed
# CoastSat rate table all cover exactly this; anything outside it is an
# extension domain.
SURVEYED_GIS = (1, 90)

# name -> (first_gis, last_gis), both inclusive.
GEOMETRIES = {
    "base": (1, 90),
    "n115": (1, 115),
}
BASE_GEOMETRY = "base"
BUFFER_DOMAINS_PER_SIDE = 15
DOMAIN_SPACING_M = 500.0


def gis_bounds(name: str) -> tuple[int, int]:
    """(first_gis, last_gis) of a named geometry; raises naming the known ones."""
    try:
        return GEOMETRIES[name]
    except KeyError:
        raise KeyError(f"no geometry {name!r}; have {sorted(GEOMETRIES)}") from None


def is_extended(name: str) -> bool:
    return gis_bounds(name) != GEOMETRIES[BASE_GEOMETRY]


def extension_gis(name: str) -> list[int]:
    """The domain numbers a geometry adds beyond the surveyed reach."""
    first, last = gis_bounds(name)
    lo, hi = SURVEYED_GIS
    return [g for g in range(first, last + 1) if g < lo or g > hi]


def geometry_label(name: str) -> str:
    """One string for a run's metadata: the name and the reach it spans."""
    first, last = gis_bounds(name)
    return f"{name} (GIS {first} to {last})"


# =============================================================================
# THE WHOLE-ISLAND POLYGONS (Hannah, 2026-09-16 evening)
# =============================================================================
# TWO JOIN RULES, because two different joins made the surveyed inputs and
# each is reproduced exactly on GIS 1-90 (checked 2026-09-16):
#   join_lines    the 100 m dune-line transects: the LINE intersects the
#                 polygon (450/450). A transect runs due west across the
#                 island and meets one polygon; one in a sliver between
#                 polygons meets none and is dropped, as ArcGIS dropped it.
#   join_origins  the CoastSat transects: the ORIGIN POINT (first vertex)
#                 within the polygon (906/906), the rule of
#                 coastsat_domain_mapping.py. (The centroid does not
#                 reproduce it: 737/906.)
EXTENDED_DOMAIN_POLYGONS = (INIT_ROOT / "1-barrier3d-domains" / "domain-geojson"
                            / "domains_pea_hatteras_120.geojson")
EXTENDED_POLYGON_ID = "ID"


@lru_cache(maxsize=1)
def domain_polygons():
    """The whole-island polygons as a GeoDataFrame (gis, geometry) in
    EPSG:3725, or None if the file is absent."""
    import geopandas as gpd
    if not EXTENDED_DOMAIN_POLYGONS.is_file():
        return None
    g = gpd.read_file(EXTENDED_DOMAIN_POLYGONS)
    g = g[[EXTENDED_POLYGON_ID, "geometry"]].rename(columns={EXTENDED_POLYGON_ID: "gis"})
    g["gis"] = g["gis"].astype(int)
    return g.to_crs(DOMAIN_CRS)


def _join(gdf, predicate):
    """gis per row of `gdf` (a GeoDataFrame in any CRS), NaN where no
    polygon matches; the first match where several do."""
    import geopandas as gpd
    polys = domain_polygons()
    if polys is None:
        raise FileNotFoundError(f"no whole-island polygons at {EXTENDED_DOMAIN_POLYGONS}")
    left = gpd.GeoDataFrame({"_row": range(len(gdf))}, geometry=gdf.geometry.values,
                            crs=gdf.crs).to_crs(DOMAIN_CRS)
    joined = gpd.sjoin(left, polys, how="left", predicate=predicate)
    return joined.groupby("_row")["gis"].first().reindex(range(len(gdf))).values


def join_lines(gdf):
    """Domain per transect LINE, by intersection (the 100 m transects)."""
    return _join(gdf, "intersects")


def join_origins(gdf):
    """Domain per transect ORIGIN (first vertex) within a polygon (CoastSat)."""
    import geopandas as gpd
    from shapely.geometry import Point
    pts = gpd.GeoDataFrame(
        geometry=[Point(g.coords[0]) if g.geom_type == "LineString"
                  else (Point(g.geoms[0].coords[0]) if g.geom_type == "MultiLineString" else g)
                  for g in gdf.geometry], crs=gdf.crs)
    return _join(pts, "within")
