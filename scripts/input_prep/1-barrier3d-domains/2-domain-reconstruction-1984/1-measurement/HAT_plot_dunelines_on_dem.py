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
from site_layer.hat_topo_version import insert_figures_dir_for_domain  # noqa: E402
from site_layer.hat_elevation_products import product                        # noqa: E402
from site_layer.hat_figure_style import (apply_style, C, C_1984, C_1997,     # noqa: E402
                              DOMAIN_AXIS_LABEL, caption, elevation_cmap,
                              figsize, open_frame, save, spines_for_image,
                              town_bands, _north_arrow, _scalebar, _title)

from site_layer.hat_topo_version import DUNELINE_DIR as DL  # noqa: E402
PRODUCT = "2009-2014-1996"
MHW_NAVD = 0.36
NODATA_BELOW = -900.0
L84, L97 = C_1984, C_1997
C_MHW = C["REF"]
# the two relocation blocks, GIS 9-14 and 84-87: the modification under test
BLOCKS = ((9, 14), (84, 87))


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

    fig = plt.figure(figsize=figsize("double", aspect=0.72),
                     constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.62, 1.0],
                          width_ratios=[1.25, 1.0])

    # ---- (a) map view ---------------------------------------------------
    # Elevation in the house classes, relative to MHW, so the water break in
    # the colour scale is the MHW contour the panel is about.
    ax = fig.add_subplot(gs[0, :])
    shown = np.where(dem > NODATA_BELOW, dem - MHW_NAVD, np.nan)
    cmap, norm, cbounds = elevation_cmap()
    extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]
    im = ax.imshow(shown, cmap=cmap, norm=norm, extent=extent,
                   origin="upper", aspect="equal", interpolation="nearest")
    # the MHW contour: the present shoreline
    ax.contour(np.flipud(shown), levels=[0.0], colors=C_MHW,
               linewidths=0.9, extent=extent)
    for g, col, lab in ((g84, L84, "1984 dune line"), (g97, L97, "1997 dune line")):
        parts = [g] if g.geom_type == "LineString" else list(g.geoms)
        first = True
        for pp in parts:
            xa, ya = pp.xy
            ax.plot(xa, ya, "-", color=col, lw=1.6,
                    label=lab if first else None, zorder=6)
            first = False
    ax.plot([], [], "-", color=C_MHW, lw=0.9, label="MHW contour")
    ax.set_xlim(bounds.left, bounds.right)
    ax.set_ylim(bounds.bottom, bounds.top)
    ax.set_xticks([])
    ax.set_yticks([])
    _title(ax, 0, "dune lines on the 1996 surface, GIS {}".format(D))
    ax.legend(loc="upper left", ncol=3)
    spines_for_image(ax)
    _scalebar(ax, 500.0, show_cells=False)
    _north_arrow(ax, x=0.955, y=0.60, length=0.16)
    # an inset beside the map, so the bar is the map's height and not the
    # gridspec row's (the map is aspect-equal and shorter than its row)
    cax = ax.inset_axes([1.012, 0.0, 0.018, 1.0])
    cb = fig.colorbar(im, cax=cax, boundaries=cbounds[1:], ticks=cbounds[1:-1])
    cb.set_label("elevation (m MHW)")
    cb.outline.set_linewidth(0.6)

    # ---- (b) what the DEM says at each line -----------------------------
    ax2 = fig.add_subplot(gs[1, 0])
    bins = np.linspace(-3, 7, 41)
    ax2.hist(v84[v84 > NODATA_BELOW], bins=bins, color=L84, alpha=0.8,
             label="1984 dune line")
    ax2.hist(v97[v97 > NODATA_BELOW], bins=bins, color=L97, alpha=0.8,
             label="1997 dune line")
    ax2.axvline(MHW_NAVD, color=C_MHW, ls="--", lw=1.0, label="MHW")
    ax2.set_xlabel("elevation at the line (m NAVD88)")
    ax2.set_ylabel("profiles")
    _title(ax2, 1, "elevation at each line, GIS {}".format(D))
    ax2.legend(loc="upper right")
    ax2.grid(axis="y")
    open_frame(ax2)

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
    town_bands(ax3)
    for lo, hi in BLOCKS:
        ax3.axvspan(lo - .5, hi + .5, color=C["ACCENT_FILL"], alpha=.45,
                    lw=0, zorder=0)
    ax3.plot(gis, m84, "-", color=L84, lw=1.1, label="at the 1984 line")
    ax3.plot(gis, m97, "-", color=L97, lw=1.1, label="at the 1997 line")
    ax3.axhline(MHW_NAVD, color=C_MHW, ls="--", lw=1.0, label="MHW")
    ax3.set_xlim(0, 91)
    ax3.set_xlabel(DOMAIN_AXIS_LABEL)
    ax3.set_ylabel("median elevation at the line\n(m NAVD88)")
    n_sub = int(np.sum(np.array(m84) <= MHW_NAVD))
    _title(ax3, 2, "median elevation, all domains")
    ax3.set_ylim(-1.9, None)   # room for the legend under the data
    ax3.legend(loc="lower center", ncol=3, columnspacing=0.8)
    ax3.grid(axis="y")
    open_frame(ax3)

    n84, w84, l84c = classify(v84)
    caption(fig,
            "The 1984 dune line is offshore in the surveyed surface, which is "
            "what sets N. (a) The digitized dune lines (1984 red, 1997 blue) on "
            "the {} product's own 10 m raster for GIS {}, in EPSG:3725, before "
            "any orient / shear / water-trim; elevation in classes relative to "
            "MHW ({} m NAVD88), with the MHW contour in green. The 1984 line "
            "crosses {} raster rows: {} below MHW, {} on land, {} no-data; the "
            "1997 line ({} rows) is on the dune throughout. (b) The raster "
            "elevation where each line crosses each row of GIS {}. (c) The "
            "median of the same quantity per domain along the island (1 at "
            "Cape Point, 90 at north Pea Island); the 1984 line lies below MHW "
            "in {} of {} domains. Grey bands are the villages, purple bands the "
            "two relocation blocks (GIS 9-14 and 84-87). N is the offshore "
            "distance MINUS the offset between a digitized line and interior "
            "row 0; inserting the full offshore distance would put row 0 on "
            "the toe rather than behind the crest."
            .format(PRODUCT, D, MHW_NAVD, len(v84), w84, l84c, n84,
                    int((v97 > NODATA_BELOW).sum()), D, n_sub, len(gis)))

    out = Path(args.out) if args.out else (
        insert_figures_dir_for_domain("1984-start", "1-measurement", D)
        / "HAT_dunelines_on_DEM_GIS{}.png".format(D))
    save(fig, out)
    print("wrote {}".format(out))
    print("  GIS {}: 1984 line median {:.2f} m NAVD88 ({} below MHW, {} land, "
          "{} nodata)".format(D, float(np.median(v84[v84 > NODATA_BELOW])),
                              w84, l84c, n84))
    print("  1984 line below MHW in {} of {} domains island-wide"
          .format(n_sub, len(gis)))


if __name__ == "__main__":
    main()
