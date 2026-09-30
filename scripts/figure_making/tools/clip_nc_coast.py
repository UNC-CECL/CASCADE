"""
Rebuild the clipped NC 1:80k coastline the figures read, from the source on the D: GIS drive.

    python scripts/figure_making/tools/clip_nc_coast.py

Keeps every polygon within NC_COAST_PAD_M of the domain boxes; writes
data/hatteras_init/map_elements/nc_coast_80k/. Needs the D: drive. Details: scripts/figure_making/tools/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
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


# Run: clip the state coastline to the padded domain window and save it
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
