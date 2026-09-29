"""
dem_composition_figures.py
==============================================================================
Which survey supplies each part of the model's topography, and how the 1 m
surveys become the 10 m grid Barrier3D reads.

    python scripts/figure_making/pipeline/0-elevation/dem_composition_figures.py

Writes to output/figures/3-model-inputs/0-elevation/:

    dem_sources_alongshore.png   per GIS domain, the share of measured cells
                                 each survey supplies, in both products:
                                 2009-2014 (the 2004/2010 starts) and
                                 2009-2014-1996 (the 1984/1996 starts)
    dem_resample_one_domain.png  one domain (GIS 45): the survey each 1 m cell
                                 comes from, the 1 m surface, and the 10 m
                                 resample the domain arrays are cut from

Read only: the products' clip_domain_<N>_{filled,survey}.tif (1 m) and
resampled_domain_<N>_{filled,survey}.tif (10 m), resolved through
scripts/site_layer/hat_elevation_products.py. Nothing is regenerated.

The tiles are north-up UTM boxes, 500 m alongshore by 2000 m cross-shore;
through the reach the ocean is to the east, so it sits at the right.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import rasterio  # noqa: E402
from matplotlib.colors import ListedColormap, BoundaryNorm  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))

from site_layer.hat_elevation_products import product  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, C_1984, C_1997, DOMAIN_AXIS_LABEL, figsize, figure_dir, save, record_caption,
    _title, open_frame, elevation_cmap, town_bands,
)

OUT = figure_dir("inputs", "0-elevation")
MHW_NAVD = 0.36
EXAMPLE_GIS = 45
# survey code -> label, colour. The vintage pair: the earlier survey red, the
# later blue; the 2009 base neutral.
SURVEYS = [(1996, "1996 NOAA/NASA ALACE (graft)", C_1984),
           (2009, "2009 USACE lidar (base)", "#b8b8b8"),
           (2014, "2014 NOAA Post-Sandy (gap fill)", C_1997)]


def read(path):
    with rasterio.open(path) as r:
        return r.read(1), r.nodata, r.bounds


def shares(prod_name):
    p = product(prod_name)
    out = np.zeros((90, len(SURVEYS)))
    for g in range(1, 91):
        s, _, _ = read(p.gapfill_1m / f"clip_domain_{g}_survey.tif")
        z, nd, _ = read(p.gapfill_1m / f"clip_domain_{g}_filled.tif")
        dry = (z != nd) & (z > MHW_NAVD)                 # land above MHW only
        n = dry.sum()
        for j, (code, _, _) in enumerate(SURVEYS):
            out[g - 1, j] = ((s == code) & dry).sum() / max(n, 1)
    return out


def fig_sources():
    prods = [("2009-2014", "2009-2014: the 2004 and 2010 starts"),
             ("2009-2014-1996", "2009-2014-1996: the 1984 and 1996 starts")]
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=4.4), sharex=True,
                             constrained_layout=True)
    g = np.arange(1, 91)
    stats = {}
    for i, (name, title) in enumerate(prods):
        sh = shares(name)
        stats[name] = sh.mean(axis=0)
        ax = axes[i]
        bottom = np.zeros(90)
        for j, (_, lab, col) in enumerate(SURVEYS):
            ax.bar(g, sh[:, j], bottom=bottom, width=0.9, color=col, lw=0)
            bottom += sh[:, j]
        ax.set_ylim(0, 1.12)
        ax.set_xlim(0.5, 90.5)
        ax.set_yticks([0, 0.5, 1])
        ax.set_ylabel("share of land\nabove MHW")
        town_bands(ax, label=(i == 0))
        open_frame(ax)
        _title(ax, i, title)
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    fig.legend(handles=[Patch(color=c, label=l) for _, l, c in SURVEYS], loc="outside lower center",
               ncol=3, frameon=False)
    out = save(fig, OUT / "dem_sources_alongshore.png")
    plt.close(fig)
    a, b = stats["2009-2014"], stats["2009-2014-1996"]
    record_caption(out[0],
        "Which survey supplies the model's topography, domain by domain: the share of each domain's "
        "1 m land cells (above MHW) taken from each survey. (a) The 2009-2014 product, which the 2004 and 2010 "
        "starts are extracted from: the 2009 USACE lidar is the base and the 2014 NOAA Post-Sandy survey "
        f"only fills its holes, no measured 2009 cell is changed (reach mean {a[1]:.0%} 2009, {a[2]:.0%} 2014). "
        "(b) The 2009-2014-1996 product for the 1984 and 1996 starts: as (a), with the 1996 NOAA/NASA ALACE "
        "survey overwriting measured ground wherever its swath has data, so the ocean side of the island is "
        f"the 1996 surface and its landward edge is a seam at the swath limit (reach mean {b[0]:.0%} 1996, "
        f"{b[1]:.0%} 2009, {b[2]:.0%} 2014). Water cells are not counted (the 2014 survey also supplies most of the sound floor). "
        "Survey acquisition dates are not recorded in the repository and are given by year only. "
        "GIS 1 is Cape Point, GIS 90 Pea Island.")
    return out


def fig_resample(gis=EXAMPLE_GIS):
    p96 = product("2009-2014-1996")
    s1, _, b = read(p96.gapfill_1m / f"clip_domain_{gis}_survey.tif")
    z1, nd1, _ = read(p96.gapfill_1m / f"clip_domain_{gis}_filled.tif")
    z10, nd10, b10 = read(p96.resampled_10m / f"resampled_domain_{gis}_filled.tif")
    s10, _, _ = read(p96.resampled_10m / f"resampled_domain_{gis}_survey.tif")
    z1 = np.where(z1 == nd1, np.nan, z1.astype(float)) - MHW_NAVD
    z10 = np.where(z10 == nd10, np.nan, z10.astype(float)) - MHW_NAVD
    ext = (0, b.right - b.left, 0, b.top - b.bottom)
    ext10 = (b10.left - b.left, b10.right - b.left, b10.bottom - b.bottom, b10.top - b.bottom)

    codes = [0] + [c for c, _, _ in SURVEYS]
    scmap = ListedColormap(["white"] + [c for _, _, c in SURVEYS])
    snorm = BoundaryNorm(np.arange(len(codes) + 1) - 0.5, scmap.N)
    to_idx = {c: i for i, c in enumerate(codes)}
    si = np.vectorize(lambda v: to_idx.get(int(v), 0))(s1)
    si10 = np.vectorize(lambda v: to_idx.get(int(v), 0))(s10)
    cmap, norm, bounds = elevation_cmap()

    fig, axes = plt.subplots(4, 1, figsize=figsize("double", height=7.0), sharex=True, sharey=True,
                             constrained_layout=True)
    panels = [(si, ext, scmap, snorm, "Which survey each 1 m cell comes from"),
              (z1, ext, cmap, norm, "The 1 m surface (m MHW)"),
              (z10, ext10, cmap, norm, "Resampled to 10 m, the grid the domain arrays are cut from"),
              (si10, ext10, scmap, snorm, "The survey each 10 m cell is attributed to")]
    for i, (arr, e, cm, nm, t) in enumerate(panels):
        ax = axes[i]
        ax.set_facecolor("white")
        im = ax.imshow(arr, cmap=cm, norm=nm, extent=e, origin="upper", interpolation="nearest",
                       aspect="auto")
        ax.set_ylabel("alongshore (m)")
        ax.set_yticks([0, 250, 500])
        open_frame(ax)
        _title(ax, i, t)
        if cm is cmap:
            cb = fig.colorbar(im, ax=ax, ticks=bounds[1:-1], fraction=0.025, pad=0.01)
            cb.outline.set_linewidth(0.5)
    axes[-1].set_xlabel("west → east (m), bay at the left, ocean at the right")
    fig.legend(handles=[Patch(color=c, label=l) for _, l, c in SURVEYS] + [Patch(fc="white", ec="0.6", label="no survey")],
               loc="outside lower center", ncol=4, frameon=False, fontsize=7.5)
    out = save(fig, OUT / "dem_resample_one_domain.png", vector=False)
    plt.close(fig)
    record_caption(out[0],
        f"From surveys to the 10 m grid, one domain (GIS {gis}, the 2009-2014-1996 product; the 2009-2014 "
        "product is the same without the 1996 layer). The domain box is 500 m alongshore by 2000 m "
        "cross-shore, north up, bay at the left and ocean at the right. (a) The survey behind each 1 m cell: "
        "the 1996 ALACE survey across the ocean side, the 2009 USACE lidar behind it, and the 2014 NOAA "
        "Post-Sandy survey where 2009 had no data. (b) The resulting 1 m surface, in the model's elevation "
        f"classes (m above MHW = NAVD88 - {MHW_NAVD} m; blue below MHW, white where no survey measured). "
        "(c) The same surface resampled to the 10 m cells Barrier3D uses (HAT_dem_resample_clip.py). "
        "(d) The survey each 10 m cell is attributed to; a mixed block goes to the most specific survey in "
        "it (hat_elevation_products.FILL_CODES).")
    return out


def main():
    apply_style()
    for f in (fig_sources, fig_resample):
        print(f.__name__, "->", f()[0].relative_to(REPO))


if __name__ == "__main__":
    main()
