#!/usr/bin/env python3
r"""
HAT_plot_insert_explainer_grid.py
==============================================================================
The 1984 seaward-row insert explained as Barrier3D PLAN-VIEW GRIDS (rows
cross-shore, columns alongshore, colour = elevation), on one domain.

    top row     (a) the survey around the pick, in the extractor's frame,
                    referenced to v2's interior row 0; 2009-survey cells
                    hatched; the crest pick and the measured road marked
                (b) the v2 domain as the model sees it: two dune rows then
                    the interior; nothing seaward of the crest exists
                (c) the inserted domain (v5 as the example): dune rows, the N
                    added rows outlined, the old crest a ridge inside, and
                    the road on rows 4-5 of the new interior
    bottom row  the first rows of v4-v8 side by side: same dune rows, same
                interior from old row 0 on, only the N added rows differ

Everything is drawn from the arrays on disk (v2, v4-v8 topography and dune
files, the survey-year clip) - no run output.

USAGE
    python HAT_plot_insert_explainer_grid.py [--domain 85]
==============================================================================
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(HERE))
import HAT_insert_seaward_rows as ins  # noqa: E402
from hat_topo_version import dune_topo_root, insert_figures_dir  # noqa: E402
from hat_topo_version import require_version  # noqa: E402
from hat_figure_style import elevation_cmap  # noqa: E402

PRODUCT = "1984-start"
BASE = "v2"
# DELETED 2026-09-07 with every layer (only unmodified topography is kept);
# the literal is kept as the name of what this drew. require_version() in
# main() says so before any array is opened.
LAYERS = [("v4", "v4  measured + floor"), ("v5", "v5  measured + median"),
          ("v6", "v6  flat platform"), ("v7", "v7  matched, crest kept"),
          ("v8", "v8  matched, crest skipped")]
CELL = 10.0
BERM_M = 1.34            # Barrier3D BermEl for these runs (m MHW); dune file is height above it
ROAD_ROWS = 2


def survey_in_z_frame(ext, dom, D):
    s = np.load(dune_topo_root(PRODUCT).parent / "npy-arrays_survey" / f"domain_{D}.npy").astype(float)
    s = ext.orient_ocean_right(s, ext.OCEAN_LOC)[:, ::-1]
    c0 = int(dom["c0"]); sh = np.asarray(dom["shear"]).astype(int)
    out = np.zeros(dom["z"].shape, int)
    for i in range(out.shape[0]):
        for k in range(out.shape[1]):
            src = k + c0 + sh[i]
            if 0 <= src < s.shape[1]:
                out[i, k] = int(s[i, src])
    return out


def model_grid(version, D, n_interior):
    """(rows x along) m MHW: two dune rows then the first n_interior rows."""
    root = dune_topo_root(PRODUCT) / version
    topo = np.load(root / "topography" / f"domain_{D}_topography.npy") * CELL
    dune = np.load(root / "dunes" / f"domain_{D}_dune.npy") * CELL + BERM_M
    dune = np.atleast_2d(dune)
    if dune.shape[0] == 1:
        dune = np.vstack([dune, dune])
    return np.vstack([dune, topo[:n_interior]])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=int, default=85)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    for _v, *_ in LAYERS:
        require_version("1984-start", _v, "LAYERS")
    D = args.domain

    off = ins.load_offset_module()
    ext = off.load_extractor(PRODUCT)
    windows = json.load(open(ext.WINDOW_JSON))
    dom = ext.load_profiles(ext.LOAD_PATH / f"domain_{D}.npy")
    prof = ext.masked_profiles(dom["z"])
    w = windows.get(f"domain_{D}")
    i0, i1 = ((int(w["i0"]), int(w["i1"])) if w else ext.default_window(prof, dom["start_beach"]))
    _e, dl = ext.find_dunes(prof, dom["start_beach"], i0, i1)
    r0 = np.asarray(ext.interior_row0_line(prof, dl)[0]).astype(int)
    survey = survey_in_z_frame(ext, dom, D)

    # survey grid referenced to each profile's own row 0
    span = np.arange(-12, 26)
    zg = np.full((span.size, len(r0)), np.nan)
    yg = np.zeros((span.size, len(r0)), int)
    for i in range(len(r0)):
        if r0[i] < 0:
            continue
        for j, k in enumerate(span):
            src = r0[i] + k
            if 0 <= src < prof.shape[1]:
                zg[j, i] = prof[i, src]
                yg[j, i] = survey[i, src]

    aud = list(csv.DictReader(open(dune_topo_root(PRODUCT) / "v4" / "HAT_seaward_row_insert_audit.csv")))
    rec = next(r for r in aud if int(r["domain"]) == D)
    n = int(rec["n_rows_inserted"])
    raw_sb = float(rec["setback_raw_before_m"])
    k_b = max(int(raw_sb / CELL), 0)                 # v2 road row
    k_c = int((raw_sb + 10 * n) / CELL)              # inserted road row

    cmap, norm, bounds = elevation_cmap()
    n_show = 26
    g2 = model_grid(BASE, D, n_show)
    g5 = model_grid("v5", D, n_show + n)
    nal = g2.shape[1]

    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(2, 5, height_ratios=(1.25, 1.0), hspace=0.32, wspace=0.22,
                          left=0.05, right=0.985, top=0.9, bottom=0.17)
    a = fig.add_subplot(gs[0, 0:2])
    b = fig.add_subplot(gs[0, 2])
    c = fig.add_subplot(gs[0, 3:5])

    def road_band(ax, row_top, label, colour="k"):
        ax.add_patch(Rectangle((-0.5, row_top - 0.5), nal, ROAD_ROWS, fill=False, ec=colour, lw=2))
        ax.text(nal - 1, row_top - 0.5 - 0.35, label, ha="right", va="bottom", fontsize=8, color=colour,
                bbox=dict(fc="white", ec="none", alpha=.8, pad=1))

    # (a) survey
    a.imshow(zg, cmap=cmap, norm=norm, aspect="auto", origin="upper",
             extent=(-0.5, nal - 0.5, span[-1] + 0.5, span[0] - 0.5))
    hatch = np.where(yg == 2009, 1.0, np.nan)
    a.imshow(hatch, cmap="gray", alpha=0.0, aspect="auto", origin="upper",
             extent=(-0.5, nal - 0.5, span[-1] + 0.5, span[0] - 0.5))
    a.contourf(np.arange(nal), span, (yg == 2009).astype(float), levels=[0.5, 1.5],
               colors="none", hatches=["////"], extend="neither")
    a.axhline(-1, color="k", lw=1.2)
    a.text(0.5, -3.3, "crest pick, row −1 (→ Barrier3D dune rows)", fontsize=8.5, va="bottom",
           bbox=dict(fc="white", ec="none", alpha=.85, pad=1))
    a.axhline(-0.5, color="k", lw=.6, ls=":")
    a.text(0.5, 1.6, "v2 interior row 0 is the row below the line", fontsize=8.5, va="top",
           bbox=dict(fc="white", ec="none", alpha=.85, pad=1))
    road_band(a, raw_sb / CELL, f"NC-12 measured at {raw_sb:+.0f} m from row 0")
    a.set_ylabel("cross-shore cell, relative to v2 row 0  (negative = seaward)")
    a.set_xlabel("alongshore cell")
    a.set_title(f"(a) the survey around the pick, GIS {D}   (hatched = 2009 survey, rest 1996)", loc="left", fontsize=10)
    a.set_ylim(span[-1] + 0.5, span[0] - 0.5)

    # (b) v2 as the model sees it
    b.imshow(g2, cmap=cmap, norm=norm, aspect="auto", origin="upper",
             extent=(-0.5, nal - 0.5, n_show - 0.5, -2.5))
    b.axhline(-0.5, color="k", lw=.8)
    b.text(0.5, -2.2, "dune rows", fontsize=8.5, va="top", color="white")
    b.text(0.5, 0.05, "row 0", fontsize=8.5, va="top")
    road_band(b, k_b, f"setback floored to 0 → rows {k_b}–{k_b+1}")
    b.set_title("(b) v2 as the model sees it", loc="left", fontsize=10)
    b.set_ylabel("model row")
    b.set_xlabel("alongshore cell")

    # (c) the insert, v5 as example
    c.imshow(g5, cmap=cmap, norm=norm, aspect="auto", origin="upper",
             extent=(-0.5, nal - 0.5, n_show + n - 0.5, -2.5))
    c.add_patch(Rectangle((-0.5, -0.5), nal, n, fill=False, ec="#d17f00", lw=2.2, ls="--"))
    c.text(0.5, 0.0, f"{n} added rows = new rows 0–{n-1}", fontsize=8.5, va="top", color="#d17f00",
           bbox=dict(fc="white", ec="none", alpha=.85, pad=1))
    c.axhline(n - 0.5, color="k", lw=.8)
    c.text(nal - 1, n + 0.05, "old row 0 — the 1996 crest, now a ridge inside", fontsize=8.5, va="top", ha="right",
           bbox=dict(fc="white", ec="none", alpha=.8, pad=1))
    c.text(0.5, -2.2, "dune rows, moved in front", fontsize=8.5, va="top", color="white")
    road_band(c, k_c, f"setback {raw_sb:+.0f}+{10*n} = {raw_sb+10*n:+.0f} m → rows {k_c}–{k_c+1}")
    c.set_title("(c) the inserted domain (v5 shown): road does not move, row 0 does", loc="left", fontsize=10)
    c.set_ylabel("model row")
    c.set_xlabel("alongshore cell")

    # bottom: the five fills, zoomed to the block
    zoom = n + 8
    for k, (v, label) in enumerate(LAYERS):
        ax = fig.add_subplot(gs[1, k])
        g = model_grid(v, D, zoom)
        ax.imshow(g, cmap=cmap, norm=norm, aspect="auto", origin="upper",
                  extent=(-0.5, nal - 0.5, zoom - 0.5, -2.5))
        ax.add_patch(Rectangle((-0.5, -0.5), nal, n, fill=False, ec="#d17f00", lw=1.8, ls="--"))
        ax.axhline(n - 0.5, color="k", lw=.8)
        ax.add_patch(Rectangle((-0.5, k_c - 0.5), nal, ROAD_ROWS, fill=False, ec="k", lw=1.6))
        ax.set_title(label, loc="left", fontsize=9.5)
        ax.set_xlabel("alongshore cell")
        if k == 0:
            ax.set_ylabel("model row")
            ax.text(nal - 1, k_c - 0.7, "road", ha="right", va="bottom", fontsize=8,
                    bbox=dict(fc="white", ec="none", alpha=.8, pad=1))
        else:
            ax.set_yticklabels([])

    cax = fig.add_axes([0.35, 0.095, 0.3, 0.015])
    fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax, orientation="horizontal",
                 ticks=[bb for bb in bounds if abs(bb) < 50])
    cax.set_xlabel("elevation (m MHW); leftmost class is water (≤ 0)", fontsize=8.5)
    cax.tick_params(labelsize=7.5)

    fig.suptitle(f"GIS {D}  ·  the 1984 seaward-row insert as Barrier3D grids", fontsize=13,
                 fontweight="bold", x=0.02, ha="left")
    fig.text(0.02, 0.005,
             "(a) is the DEM in the extractor's frame, each profile referenced to its own pick; the cells above the crest line are "
             "the 1996 beach and seaward dune face the extraction drops. (b) is what Barrier3D receives: the picked crest as two dune "
             "rows, then the interior from row 0. (c) prepends N rows at the dropped coordinates and puts the dune rows in front, so "
             "row 0 moves seaward by N while NC-12 stays on the ground it was measured on; the old crest remains inside. The bottom "
             "row shows the five fills of those N rows; everything below the black line is identical in all five and identical to v2.",
             fontsize=8, wrap=True, va="bottom")
    out = Path(args.out) if args.out else insert_figures_dir(PRODUCT, "3-fill") / f"HAT_insert_explainer_grid_GIS{D}.png"
    fig.savefig(out, dpi=150)
    print("wrote", out)


if __name__ == "__main__":
    main()
