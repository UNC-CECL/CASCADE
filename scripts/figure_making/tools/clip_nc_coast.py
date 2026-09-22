"""
clip_nc_coast.py
==============================================================================
Rebuild data/hatteras_init/map_elements/nc_coast_80k/ from the NC 1:80k
coastline on the D: GIS drive.

    python scripts/figure_making/tools/clip_nc_coast.py

WHY
    The overwash maps drew the coast straight off
    D:/Hatteras_GIS/Outlines/nc_80k/nc_80k.shp -- 6 MB for the whole state --
    so they could not be made without the drive. They only ever use the land
    in a window a few km around the domain boxes, so the repository keeps that
    window: every polygon within NC_COAST_PAD_M (25 km) of the boxes, clipped,
    as a geojson (0.75 MB) in EPSG:4326 like its source.

    Re-run it only if the source changes; the output is what the figures read.
==============================================================================
"""

from __future__ import annotations

import sys
from pathlib import Path

import geopandas as gpd
from shapely.geometry import box

REPO = next(p for p in Path(__file__).resolve().parents
            if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))

from site_layer import hat_map_layers as ml  # noqa: E402


def main():
    if not ml.NC_COAST_SOURCE.exists():
        raise SystemExit(f"\n{ml.NC_COAST_SOURCE} is not reachable: connect "
                         f"the drive.")
    dom = gpd.read_file(ml.DOMAIN_BOXES).to_crs(26918)
    b = dom.total_bounds
    pad = ml.NC_COAST_PAD_M
    window = box(b[0] - pad, b[1] - pad, b[2] + pad, b[3] + pad)
    coast = gpd.read_file(ml.NC_COAST_SOURCE)
    clipped = gpd.clip(coast.to_crs(26918), window).to_crs(coast.crs)
    ml.NC_COAST.parent.mkdir(parents=True, exist_ok=True)
    clipped.to_file(ml.NC_COAST, driver="GeoJSON")
    print(f"  {len(clipped)} of {len(coast)} polygons within {pad / 1000:.0f} km "
          f"of the domain boxes -> {ml.NC_COAST.relative_to(REPO).as_posix()} "
          f"({ml.NC_COAST.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
