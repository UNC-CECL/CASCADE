"""
How far does the 1996 ALACE survey sit above the 2009 USACE lidar, measured on NC-12 pavement?

    python scripts/input_prep/0-elevation/1-source-selection/HAT_alace_vertical_offset.py

Pavement does not grow or erode, so the median ALACE-minus-2009 difference along
the 2008 NC-12 line, per domain, is the offset between the two surveys. Read from
the built 1 m clips of the 1996-start (2009-2014-1996) and 2009-start
(2009-2014) products, on pixels each product took from its own survey; the
offset the 1996 product was built with (its audit's alace_vertical_offset_m,
0 if absent) is added back, so the answer is the raw survey offset whether or
not the product is corrected. GIS
9-13 are left out of the island value: NC-12 was relocated there in 1999, so the
two surveys see different pavement. A low-flat comparison (2009 ground between
0.6 and 2 m NAVD88) is reported beside it as a looser check that includes real
change. Writes alace_vertical_offset.csv beside the 1996 product. Local data
only. Details: scripts/input_prep/0-elevation/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-20
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import rasterize

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(_HERE.parents[1] / "2-produce"))
from site_layer import hat_topo_version as tv  # noqa: E402
from site_layer.hat_elevation_products import product  # noqa: E402
import HAT_dem_1984_mosaic as M  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
PRODUCT_1996, PRODUCT_2009 = "2009-2014-1996", "2009-2014"
ROAD_LINE = tv.MGMT_ROOT / "road_offset" / "raw_offset" / "2008" / "nc12_2008.geojson"
RELOCATED_1999 = range(9, 14)      # GIS 9-13: NC-12 moved landward in 1999
MIN_ROAD_PX = 50
FLATS_NAVD = (0.6, 2.0)            # m NAVD88, 2009 ground for the looser check
MIN_FLAT_PX = 500
OUT = product(PRODUCT_1996).root / "alace_vertical_offset.csv"
# -----------------------------------------------------------------------------


def read(root, g, kind):
    with rasterio.open(root / f"clip_domain_{g}_{kind}.tif") as s:
        return s.read(1), s.nodata, s.transform, s.crs


# The ALACE offset a built 1996 product carries, from its audit (built before 2026-10-20: none)
def built_offset():
    a = pd.read_csv(product(PRODUCT_1996).audit_1m)
    if "alace_vertical_offset_m" not in a:
        return 0.0
    v = a.alace_vertical_offset_m.unique()
    if len(v) != 1:
        raise SystemExit(f"audit mixes ALACE offsets {v}: rebuild the product")
    return float(v[0])


def main():
    line = gpd.read_file(ROAD_LINE)
    built = built_offset()
    r96, r09 = product(PRODUCT_1996).gapfill_1m, product(PRODUCT_2009).gapfill_1m
    rows = []
    for g in range(1, 91):
        z96, nd96, tr, crs = read(r96, g, "filled")
        src96 = read(r96, g, "survey")[0]
        z09, nd09, _, _ = read(r09, g, "filled")
        src09 = read(r09, g, "survey")[0]
        z96 = z96.astype(float) - built                             # back to the raw ALACE return
        z09 = z09.astype(float)
        both = (src96 == M.OVERRIDE_YEAR) & (src09 == 2009) & (z09 != nd09) & (z09 > -50)
        road = rasterize([(x, 1) for x in line.to_crs(crs).geometry], out_shape=z96.shape, transform=tr,
                         fill=0).astype(bool) & both
        flats = both & (z09 > FLATS_NAVD[0]) & (z09 < FLATS_NAVD[1])
        d, f = z96[road] - z09[road], z96[flats] - z09[flats]
        rows.append(dict(gis=g, road_px=int(road.sum()),
                         road_offset_m=round(float(np.median(d)), 3) if road.sum() >= MIN_ROAD_PX else np.nan,
                         road_iqr_m=round(float(np.subtract(*np.percentile(d, [75, 25]))), 3)
                         if road.sum() >= MIN_ROAD_PX else np.nan,
                         relocated_1999=g in RELOCATED_1999, flats_px=int(flats.sum()),
                         flats_offset_m=round(float(np.median(f)), 3) if flats.sum() >= MIN_FLAT_PX else np.nan))
    t = pd.DataFrame(rows)
    t.to_csv(OUT, index=False)
    v = t[~t.relocated_1999].road_offset_m.dropna()
    print(f"ALACE minus 2009 on NC-12: median {v.median():+.3f} m over {len(v)} domains "
          f"(IQR {v.quantile(.25):+.3f} to {v.quantile(.75):+.3f}); product built with {built:+.2f} m, "
          f"mosaic configured {M.OVERRIDE_VERTICAL_OFFSET_M:+.2f} m")
    print(f"the built product sits {v.median() + built:+.3f} m from 2009 on NC-12")
    print(f"low flats: median {t.flats_offset_m.median():+.3f} m (includes real change)")
    print(OUT)


if __name__ == "__main__":
    main()
