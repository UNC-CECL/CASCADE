"""
The alongshore reach a run models, by name, and the domain numbers beyond the 90 surveyed ones.

    from site_layer.hat_extension_domains import GEOMETRIES, gis_bounds, join_lines

"base" is GIS 1-90; an extended geometry adds Pea Island domains numbered by
the whole-island polygons. HAT_GEOMETRY in the environment selects one. Details: scripts/site_layer/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents
                    if (_p / "pyproject.toml").exists())
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"

DOMAIN_CRS = "EPSG:3725"

# The surveyed reach: topography, management inputs and the CoastSat table cover exactly this
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


# The whole-island polygons (Hannah, 2026-09-16 evening)

# Two join rules, each reproducing its surveyed input: join_lines (dune-line), join_origins (CoastSat)
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
