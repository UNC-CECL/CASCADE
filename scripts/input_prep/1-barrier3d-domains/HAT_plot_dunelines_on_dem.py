#!/usr/bin/env python3
r"""
HAT_plot_dunelines_on_dem.py
==============================================================================
The dune lines on the RAW DEM, in map coordinates, before any Barrier3D
processing -- and what the DEM says at each line, island-wide.

WHY THIS FIGURE EXISTS
    The 1984 dune line is SUBMERGED in the surveyed surface. That is not a
    problem with the measurement, it is the measurement: by 1996 the island had
    retreated past where its 1984 dune stood, so the ground at that position is
    now below MHW and no survey covers it. How far offshore the line sits is
    what sets how many interior rows have to be added.

    Drawing it on a processed Barrier3D grid cannot show this, because that grid
    starts at the water trim and is expressed in cells. This draws the geometry
    on the product's own raster, in metres, in the raster's own CRS.

    DEM: 0-elevation/2009-2014-1996 (2-resampled-10m) -- the 1996 ALACE graft,
    which is the product the 1984-start arrays are exported from. Resolved
    through hat_elevation_products, not by joining strings.

WHAT N IS, IN THESE TERMS
    N is NOT how far offshore the 1984 line is. That distance -- row 0 to the
    1984 line -- also contains the offset between a digitized line and the
    model's interior row 0, which the 1997 line measures separately:

        offshore distance (total)  =  N  +  (line vs row 0)

    Inserting the full offshore distance would put row 0 on the digitized line,
    a light/dark break at the dune toe, when row 0 is one cell landward of the
    crest. At GIS 85 that would over-insert by ~1.7 cells.

USAGE
    python HAT_plot_dunelines_on_dem.py [--domain 85]
==============================================================================
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import rasterio
from pyproj import Transformer
from shapely.geometry import LineString, shape
from shapely.ops import transform as sh_transform, unary_union


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
from hat_topo_version import insert_figures_dir  # noqa: E402
from hat_elevation_products import product                        # noqa: E402
from hat_figure_style import (apply_style, C, caption,            # noqa: E402
                              panel_title, spines_for_image)

DL = REPO / "data/hatteras_init/1-barrier3d-domains/raw-duneline-geojson"
PRODUCT = "2009-2014-1996"
MHW_NAVD = 0.36
NODATA_BELOW = -900.0
L84, L97 = C["ACCENT"], "#1b6ca8"


def load_line(year, dst_crs):
    gj = json.load(open(DL / "duneline_{}.geojson".format(year)))
    src = gj.get("crs", {}).get("properties", {}).get("name", "EPSG:26918")
    geom = unary_union([shape(f["geometry"]) for f in gj["features"]])
    tr = Transformer.from_crs(src, dst_crs, always_xy=True)
    return sh_transform(lambda x, y, z=None: tr.transform(x, y), geom)


def tif_for(domain):
    d = product(PRODUCT).resampled_10m
    for nm in ("resampled_domain_{}_filled.tif".format(domain),
               "resampled_domain_{}.tif".format(domain)):
        if (d / nm).is_file():
            return d / nm
    return None


def sample_along(line, dem, T, shape_, keep_cols=True):
    """DEM values where `line` crosses each raster row of this domain."""
    vals, xs, ys = [], [], []
    for r in range(shape_[0]):
        yy = (T * (0, r + 0.5))[1]
        cut = line.intersection(LineString(
            [(T.c - 200, yy), (T.c + shape_[1] * T.a + 200, yy)]))
        pts = [cut] if cut.geom_type == "Point" else list(getattr(cut, "geoms", []))
        for p in pts:
            col = int((p.x - T.c) / T.a)
            if 0 <= col < shape_[1]:
                vals.append(dem[r, col])
                xs.append(p.x)
                ys.append(yy)
    return np.array(vals), np.array(xs), np.array(ys)


def classify(v):
    v = np.asarray(v, dtype=float)
    ok = v > NODATA_BELOW
    return (int((~ok).sum()),
            int((ok & (v <= MHW_NAVD)).sum()),
            int((ok & (v > MHW_NAVD)).sum()))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=int, default=85)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    D = args.domain
    apply_style()

    tif = tif_for(D)
    if tif is None:
        raise SystemExit("no raster for domain {}".format(D))
    with rasterio.open(tif) as s:
        dem = s.read(1).astype(float)
        T, crs, shp = s.transform, s.crs, s.shape
        bounds = s.bounds
    g84, g97 = load_line(1984, crs), load_line(1997, crs)

    v84, x84, y84 = sample_along(g84, dem, T, shp)
    v97, x97, y97 = sample_along(g97, dem, T, shp)

    fig = plt.figure(figsize=(11.2, 9.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0],
                          width_ratios=[1.35, 1.0])

    # ---- (a) map view ---------------------------------------------------
    ax = fig.add_subplot(gs[0, :])
    shown = np.where(dem > NODATA_BELOW, dem, np.nan)
    im = ax.imshow(shown, cmap="terrain", vmin=-1.5, vmax=6.0,
                   extent=[bounds.left, bounds.right, bounds.bottom, bounds.top],
                   origin="upper", aspect="equal", interpolation="nearest")
    # the MHW contour: the present shoreline
    ax.contour(np.flipud(shown), levels=[MHW_NAVD], colors="#0d3b66",
               linewidths=1.1,
               extent=[bounds.left, bounds.right, bounds.bottom, bounds.top])
    for g, col, lab in ((g84, L84, "1984 dune line"), (g97, L97, "1997 dune line")):
        parts = [g] if g.geom_type == "LineString" else list(g.geoms)
        first = True
        for pp in parts:
            xa, ya = pp.xy
            ax.plot(xa, ya, "-", color=col, lw=2.2,
                    label=lab if first else None, zorder=6)
            first = False
    ax.set_xlim(bounds.left, bounds.right)
    ax.set_ylim(bounds.bottom, bounds.top)
    ax.set_xlabel("easting (m, EPSG:3725)")
    ax.set_ylabel("northing (m)")
    ax.set_title(panel_title("a", "GIS {} — the dune lines on the raw "
                                  "{} DEM (dark blue = MHW)".format(D, PRODUCT)))
    ax.legend(fontsize=8, loc="lower left", framealpha=0.9)
    spines_for_image(ax)
    cb = fig.colorbar(im, ax=ax, fraction=0.020, pad=0.012)
    cb.set_label("elevation (m NAVD88)", fontsize=7.6)
    cb.ax.tick_params(labelsize=7)

    # ---- (b) what the DEM says at each line -----------------------------
    ax2 = fig.add_subplot(gs[1, 0])
    bins = np.linspace(-3, 7, 41)
    ax2.hist(v84[v84 > NODATA_BELOW], bins=bins, color=L84, alpha=0.75,
             label="1984 line  (n={})".format(int((v84 > NODATA_BELOW).sum())))
    ax2.hist(v97[v97 > NODATA_BELOW], bins=bins, color=L97, alpha=0.75,
             label="1997 line  (n={})".format(int((v97 > NODATA_BELOW).sum())))
    ax2.axvline(MHW_NAVD, color="#0d3b66", ls="--", lw=1.3, label="MHW")
    ax2.set_xlabel("DEM elevation at the line (m NAVD88)")
    ax2.set_ylabel("profiles")
    ax2.set_title(panel_title("b", "The 1984 line is below MHW; the 1997 "
                                   "line is on the dune"))
    ax2.legend(fontsize=7.6)
    ax2.grid(alpha=.3)

    # ---- (c) island-wide: is the 1984 line submerged? -------------------
    ax3 = fig.add_subplot(gs[1, 1])
    gis, m84, m97, sub = [], [], [], []
    for dd in range(1, 91):
        t = tif_for(dd)
        if t is None:
            continue
        with rasterio.open(t) as s:
            de = s.read(1).astype(float)
            TT, cc, sh = s.transform, s.crs, s.shape
        a, _, _ = sample_along(load_line(1984, cc), de, TT, sh)
        b, _, _ = sample_along(load_line(1997, cc), de, TT, sh)
        if a.size == 0 or b.size == 0:
            continue
        ok_a, ok_b = a > NODATA_BELOW, b > NODATA_BELOW
        if not (ok_a.any() and ok_b.any()):
            continue
        gis.append(dd)
        m84.append(float(np.median(a[ok_a])))
        m97.append(float(np.median(b[ok_b])))
        sub.append(100.0 * float(np.mean(a[ok_a] <= MHW_NAVD)))
    gis = np.array(gis)
    ax3.plot(gis, m84, "-", color=L84, lw=1.4, label="at the 1984 line")
    ax3.plot(gis, m97, "-", color=L97, lw=1.4, label="at the 1997 line")
    ax3.axhline(MHW_NAVD, color="#0d3b66", ls="--", lw=1.2, label="MHW")
    for lo, hi in ((9, 14), (84, 87)):
        ax3.axvspan(lo - .5, hi + .5, color="#ffe9b0", alpha=.5, zorder=0)
    ax3.set_xlim(0, 91)
    ax3.set_xlabel("GIS domain")
    ax3.set_ylabel("median DEM elevation\nat the line (m NAVD88)")
    n_sub = int(np.sum(np.array(m84) <= MHW_NAVD))
    ax3.set_title(panel_title("c", "1984 line below MHW in {} of {} domains"
                              .format(n_sub, len(gis))))
    ax3.legend(fontsize=7.6)
    ax3.grid(alpha=.3)

    fig_h = fig.get_figheight()
    n84, w84, l84c = classify(v84)
    fig.suptitle("The 1984 dune line is offshore — which is what sets N",
                 fontsize=11, fontweight="bold", x=0.055, ha="left",
                 y=1 - 0.24 / fig_h)
    caption(fig,
            "Geometry drawn on the product's own 10 m raster in EPSG:3725, "
            "before any orient / shear / water-trim. GIS {}: the 1984 line "
            "crosses {} profiles — {} below MHW, {} on land,\n"
            "{} no-data — while the 1997 line is on the dune throughout. N is "
            "the offshore distance MINUS the offset between a digitized line "
            "and interior row 0; inserting the full offshore\n"
            "distance would put row 0 on the toe rather than behind the crest."
            .format(D, len(v84), w84, l84c, n84),
            y=0.075 / fig_h, size=7.4)
    fig.subplots_adjust(top=1 - 0.58 / fig_h, bottom=0.80 / fig_h,
                        left=0.085, right=0.975, hspace=0.42, wspace=0.26)

    out = Path(args.out) if args.out else (
        insert_figures_dir("1984-start")
        / "HAT_dunelines_on_DEM_GIS{}.png".format(D))
    fig.savefig(out)
    print("wrote {}".format(out))
    print("  GIS {}: 1984 line median {:.2f} m NAVD88 ({} below MHW, {} land, "
          "{} nodata)".format(D, float(np.median(v84[v84 > NODATA_BELOW])),
                              w84, l84c, n84))
    print("  1984 line below MHW in {} of {} domains island-wide"
          .format(n_sub, len(gis)))


if __name__ == "__main__":
    main()
