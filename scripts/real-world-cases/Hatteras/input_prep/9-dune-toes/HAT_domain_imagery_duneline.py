"""
One domain on the 1996 aerial photo (2026-11-05): the photo, the domain outline and the dune line Hannah traced.

    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_domain_imagery_duneline.py [--gis 85]

The photo is the USGS georeferenced scan of the 14 Oct 1996 USACE flight (aerial_imagery/1996/usgs_henderson/frames,
NC State Plane ft), resampled to 0.5 m on the domains' CRS. There is no 1996 dune line: the nearest is the one traced
on the 12 Oct 1997 photo (2-brie-offset/dunelines/duneline_1997.geojson), a year later. GIS 85 and 89 are where that
line sits closest to the 1997 lidar toes (dune-collision/1-observations/6-toes-vs-dunelines: most cells within 10 m).

Writes 6-imagery/domain_gis<NN>_1996_photo_duneline.png beside the toe tables (PDF under supporting/).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-05
"""
from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import HAT_toe_common as DT  # noqa: E402

PHOTO_GLOB = DT.S.GIS_ROOT / "aerial_imagery" / "1996" / "usgs_henderson" / "frames" / "*.tif"
DUNELINE = DT.S.REPO / "data" / "hatteras_init" / "2-brie-offset" / "dunelines" / "duneline_1997.geojson"
PAD_M = 100.0             # alongshore, beyond the domain
LAND_M, SEA_M = 220.0, 160.0   # cross-shore view: landward and seaward of the dune line
RES_M = 0.5
LINE_C, DOMAIN_C = "#00e5ff", "#ffd60a"
OUT_DIR = DT.OUT / "6-imagery"


def photo(crs, bounds):
    """The 1996 frames under `bounds` (domains' CRS), RGB uint8; the first frame with pixels wins."""
    import rasterio
    from rasterio.transform import from_origin
    from rasterio.warp import reproject, transform_bounds, Resampling
    x0, y0, x1, y1 = bounds
    w, h = int((x1 - x0) / RES_M), int((y1 - y0) / RES_M)
    t = from_origin(x0, y1, RES_M, RES_M)
    out = np.zeros((3, h, w), np.uint8)
    for f in sorted(glob.glob(str(PHOTO_GLOB))):
        with rasterio.open(f) as src:
            sb = transform_bounds(crs, src.crs, *bounds)
            if sb[2] < src.bounds.left or sb[0] > src.bounds.right or sb[3] < src.bounds.bottom or sb[1] > src.bounds.top:
                continue
            a = np.zeros_like(out)
            reproject(rasterio.band(src, [1, 2, 3]), a, src_nodata=0, dst_transform=t, dst_crs=crs, dst_nodata=0,
                      resampling=Resampling.bilinear)
            empty = (out == 0).all(axis=0)
            out[:, empty] = a[:, empty]
        if not (out == 0).all(axis=0).any():
            break
    return out


def figure(gis):
    import geopandas as gpd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    doms = gpd.read_file(DT.S.GIS_ROOT / "domains.geojson")
    poly = doms[doms.domain_id == gis].geometry.iloc[0]
    line = gpd.read_file(DUNELINE).to_crs(doms.crs).geometry.union_all()
    a, b, c, d = poly.bounds
    from shapely.geometry import box
    seg = line.intersection(box(a, b, c, d))                  # the dune line inside the domain; the ocean is east
    lx = float(np.mean([x for part in getattr(seg, "geoms", [seg]) for x in part.xy[0]]))
    bounds = (lx - LAND_M, b - PAD_M, lx + SEA_M, d + PAD_M)
    img = photo(doms.crs, bounds)
    if (img == 0).all():
        raise SystemExit(f"no 1996 photo covers GIS {gis}")

    aspect = (bounds[3] - bounds[1]) / (bounds[2] - bounds[0])
    fig, ax = plt.subplots(figsize=fs.figsize("single", aspect=aspect), layout="constrained")
    ax.imshow(np.moveaxis(img, 0, -1), extent=(bounds[0], bounds[2], bounds[1], bounds[3]), origin="upper")
    for g in [poly.boundary]:
        for part in getattr(g, "geoms", [g]):
            ax.plot(*part.xy, color=DOMAIN_C, lw=1.4, label="Domain outline")
    view = line.intersection(box(*bounds))
    for k, part in enumerate(getattr(view, "geoms", [view])):
        if not part.is_empty:
            ax.plot(*part.xy, color=LINE_C, lw=1.4, label="Dune line (traced, Oct 1997 photo)" if k == 0 else None)
    ax.set_xlim(bounds[0], bounds[2])
    ax.set_ylim(bounds[1], bounds[3])
    ax.set_xticks([])
    ax.set_yticks([])
    fs.spines_for_image(ax)
    fs.scalebar(ax, length_m=50, show_cells=False)
    fs.north_arrow(ax)
    h, l = ax.get_legend_handles_labels()
    keep = dict(zip(l, h))
    ax.legend(keep.values(), keep.keys(), loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1, fontsize=7,
              frameon=False)
    ax.text(0.98, 0.98, f"GIS {gis}\n14 Oct 1996 photo", transform=ax.transAxes, ha="right", va="top", fontsize=8,
            fontweight="bold", color="white", path_effects=fs.MAP_HALO)
    png = OUT_DIR / f"domain_gis{gis:02d}_1996_photo_duneline.png"
    fs.save(fig, png, close=True)
    fs.record_caption(png, (
        f"GIS domain {gis} (yellow outline) on the 14 Oct 1996 USACE aerial photo (USGS georeferenced scan), with the "
        f"dune line traced on the 12 Oct 1997 photo (cyan), the nearest traced line to 1996. North up."))
    print(f"wrote {png}")
    return png


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gis", type=int, default=85)
    figure(ap.parse_args().gis)


if __name__ == "__main__":
    main()
