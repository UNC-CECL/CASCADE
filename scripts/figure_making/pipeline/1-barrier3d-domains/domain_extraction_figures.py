"""
domain_extraction_figures.py
==============================================================================
How one domain's Barrier3D arrays are cut from the 10 m DEM, step by step,
using the extractor's own functions on its own saved picks.

    python scripts/figure_making/pipeline/1-barrier3d-domains/domain_extraction_figures.py [--gis 45] [--product 2004-start]

Writes output/figures/3-model-inputs/1-domains/domain_extraction_gis<N>.png.

WHAT IT CALLS
    HAT_dune_topo_extractor.load_profiles (orient, MHW, clamp, beach start,
    straighten, trim), find_dunes (the crest inside the picked window) and
    build_interior (everything landward of the crest), with the window from
    the product's picks file. extract_domain() is NOT called: it writes the
    arrays. Instead the result is checked against the saved arrays of the
    product's CURRENT version, and the script stops if they differ, so the
    figure shows what the model reads. The road overlay is switched off.

Every cross-shore panel has the ocean at the RIGHT.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "1-barrier3d-domains" / "1-extraction"))

import HAT_dune_topo_extractor as ex  # noqa: E402
from site_layer import hat_topo_version as tv  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, C, C_1984, C_1997, INK, INK_MUTED, CELL_M, figsize, figure_dir, save,
    record_caption, _title, open_frame, elevation_cmap,
)

OUT = figure_dir("inputs", "1-domains")
PROFILES = (8, 25, 42)          # alongshore profiles drawn in panel (c)


def extract(gis, product):
    ex.SHOW_ROAD = False
    dom = ex.load_profiles(tv.npy_dirs(product)[0] / f"domain_{gis}.npy")
    root = tv.dune_topo_root(product)
    ver = (root / "CURRENT").read_text(encoding="utf-8").strip()
    picks = json.load(open(root.parent / "1-extraction" / "picks" / f"HAT_dune_search_windows_{ver}.json"))
    w = picks[f"domain_{gis}"]
    i0, i1 = int(w["i0"]), int(w["i1"])
    de, dl = ex.find_dunes(dom["z"], dom["start_beach"], i0, i1)
    dh = de - (ex.BERM_ELEV_NAVD_M - ex.MHW_M)
    dh[np.isfinite(dh) & (dh < 0)] = ex.MIN_DUNE_H_M
    topo, _ = ex.build_interior(dom["z"], dl)
    topo = ex.remove_water_rows(topo, ex.SENTINEL_WATER_M)
    topo = np.where(topo <= ex.NODATA_SENTINEL_M + 1e-9, ex.SENTINEL_WATER_M, topo)
    st = np.load(root / ver / "topography" / f"domain_{gis}_topography.npy")
    sd = np.load(root / ver / "dunes" / f"domain_{gis}_dune.npy").ravel()[:ex.ALONG_COLS]
    ok = (st.shape == topo.shape and np.allclose(st, topo * 0.1)
          and np.allclose(sd, np.nan_to_num(dh[:ex.ALONG_COLS], nan=ex.MIN_DUNE_H_M) * 0.1))
    if not ok:
        raise SystemExit(f"GIS {gis}: the extractor functions do not reproduce {product}/{ver}; not drawing")
    return dom, i0, i1, de, dl, dh, topo, st, sd, ver


def fig_extraction(gis, product):
    dom, i0, i1, de, dl, dh, topo, st, sd, ver = extract(gis, product)
    cmap, norm, bounds = elevation_cmap()
    raw = dom["raw"] - ex.MHW_M                                  # ocean first, m MHW
    raw = np.where(dom["raw"] <= ex.RAW_NODATA_MAX_NAVD, np.nan, raw)
    z = np.where(dom["z"] <= ex.NODATA_SENTINEL_M + 1e-9, np.nan, dom["z"])
    n_along = raw.shape[0]
    y_al = (np.arange(n_along) + 0.5) * CELL_M
    above = raw > ex.BEACH_START_THR_M
    sb_raw = np.where(above.any(axis=1), above.argmax(axis=1), -1)
    ok = sb_raw >= 0
    fit = np.polyval(np.polyfit(np.arange(n_along)[ok], sb_raw[ok], 1), np.arange(n_along))

    fig = plt.figure(figsize=figsize("double", height=8.6), constrained_layout=True)
    gs = fig.add_gridspec(4, 2, width_ratios=[1, 0.02], height_ratios=[1, 1, 1.15, 1.2])

    def plan(ax, arr, x_extent, i):
        im = ax.imshow(arr, cmap=cmap, norm=norm, origin="lower", aspect="auto", interpolation="nearest",
                       extent=(0, x_extent, 0, n_along * CELL_M))
        ax.set_facecolor("white")
        ax.invert_xaxis()
        ax.set_ylabel("alongshore (m)")
        ax.set_yticks([0, 250, 500])
        open_frame(ax)
        return im

    ax = fig.add_subplot(gs[0, 0])
    im = plan(ax, raw, raw.shape[1] * CELL_M, 0)
    ax.plot((sb_raw[ok] + 0.5) * CELL_M, y_al[ok], ".", color=C["ACCENT"], ms=3)
    ax.plot((fit + 0.5) * CELL_M, y_al, color=C["ACCENT"], lw=1.0)
    ax.set_xlabel("from the ocean edge of the tile (m)")
    _title(ax, 0, f"The 10 m tile, ocean at the right; beach start and its fit ({dom['obliquity_deg']:.1f}°)")
    cax = fig.add_subplot(gs[0:2, 1])
    cb = fig.colorbar(im, cax=cax, ticks=bounds[1:-1])
    cb.set_label("m MHW")
    cb.outline.set_linewidth(0.5)

    ax = fig.add_subplot(gs[1, 0])
    plan(ax, z, z.shape[1] * CELL_M, 1)
    ax.axvspan(i0 * CELL_M, i1 * CELL_M, color=C["ADDED"], alpha=0.25, lw=0)
    good = dl >= 0
    ax.plot((dl[good] + 0.5) * CELL_M, y_al[good], "o", color=C_1984, ms=2.5)
    ax.set_xlabel("from the first non-water cell, straightened (m)")
    _title(ax, 1, "Sheared straight and trimmed; the picked window and each profile's crest")

    ax = fig.add_subplot(gs[2, 0])
    berm = ex.BERM_ELEV_NAVD_M - ex.MHW_M
    x = (np.arange(z.shape[1]) + 0.5) * CELL_M
    cols = [C_1997, INK, C["BASE"]]
    for k, p in enumerate(PROFILES):
        ax.plot(x, z[p], color=cols[k], lw=1.0, label=f"profile at {p * CELL_M:.0f} m")
        if dl[p] >= 0:
            ax.plot((dl[p] + 0.5) * CELL_M, z[p, dl[p]], "o", color=C_1984, ms=4, zorder=4)
    ax.axvspan(i0 * CELL_M, i1 * CELL_M, color=C["ADDED"], alpha=0.25, lw=0)
    ax.axhline(berm, color=INK_MUTED, lw=0.8, ls=(0, (3, 2)))
    ax.axhline(0, color=C_1997, lw=0.5)
    ax.set_xlim(min(x[-1], (i1 + 45) * CELL_M), 0)
    ax.set_ylabel("elevation (m MHW)")
    ax.set_xlabel("from the first non-water cell (m)")
    open_frame(ax)
    ax.legend(handles=ax.get_legend_handles_labels()[0] + [
        Line2D([], [], color=C_1984, marker="o", ls="", ms=4, label="crest: the dune"),
        Line2D([], [], color=INK_MUTED, lw=0.8, ls=(0, (3, 2)), label=f"berm {berm:.2f} m")],
        frameon=False, fontsize=7, ncol=2, loc="upper left")
    _title(ax, 2, "Along a profile: the highest cell in the window is the dune; landward of it is interior")

    sub = gs[3, 0].subgridspec(1, 2, width_ratios=[1, 0.16], wspace=0.04)
    ax = fig.add_subplot(sub[0, 0])
    ext = (0, st.shape[0] * CELL_M, 0, st.shape[1] * CELL_M)
    ax.imshow(st.T * 10.0, cmap=cmap, norm=norm, origin="lower", aspect="auto", interpolation="nearest",
                    extent=ext)
    ax.invert_xaxis()
    ax.set_ylabel("alongshore (m)")
    ax.set_yticks([0, 250, 500])
    ax.set_xlabel("landward of the dune (m)")
    open_frame(ax)
    _title(ax, 3, f"domain_{gis}_topography.npy: {st.shape[0]} × {st.shape[1]}, dam MHW")
    axd = fig.add_subplot(sub[0, 1], sharey=ax)
    axd.barh(y_al, sd, height=CELL_M * 0.85, color=C["ADDED"], lw=0)
    axd.tick_params(labelleft=False)
    axd.set_xlabel("dam")
    open_frame(axd)
    axd.set_title("dune.npy", fontsize=8, loc="center")

    out = save(fig, OUT / f"domain_extraction_gis{gis}.png", vector=False)
    plt.close(fig)
    record_caption(out[0],
        f"How one domain's Barrier3D arrays are made: GIS {gis}, the {product} product, dune-topo {ver}, drawn "
        "by the extractor's own functions (HAT_dune_topo_extractor.py) on its saved window pick, and checked "
        "to reproduce the saved arrays exactly. Ocean at the right throughout. (a) The 10 m DEM tile, 50 "
        "profiles by 500 m, in m above MHW (NAVD88 - 0.36 m); purple points are each profile's beach start, "
        f"the first cell above {ex.BEACH_START_THR_M} m MHW, and the line their linear fit, whose "
        f"{dom['obliquity_deg']:.1f}° angle is the shoreline's obliquity to the grid. (b) Each profile shifted "
        "by the fit so the shoreline runs straight, with water beyond the island trimmed (water below "
        f"{ex.WATER_CLAMP_M:.0f} m is clamped); the amber band is the dune search window picked by hand for "
        f"this domain (cells {i0}-{i1}), and red points are the highest cell each profile has inside it. (c) "
        "Three profiles: the crest is the dune, whose height above the berm "
        f"({ex.BERM_ELEV_NAVD_M} m NAVD88) becomes the dune array; everything landward of the crest becomes "
        "the interior. (d) The two files the model reads: the interior, cross-shore rows from the dune "
        "landward in decametres above MHW (drawn here in m, the same classes as above), trailing all-water "
        "rows trimmed; and the dune height above the berm per profile, in decametres (a profile with no dune "
        f"or one below the berm gets {ex.MIN_DUNE_H_M} m).")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gis", type=int, default=45)
    ap.add_argument("--product", default="2004-start")
    a = ap.parse_args()
    apply_style()
    print("->", fig_extraction(a.gis, a.product)[0].relative_to(REPO))


if __name__ == "__main__":
    main()
