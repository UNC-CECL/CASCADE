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
REPO = next(
    _p for _p in HERE.parents
    if (_p / "pyproject.toml").exists())                  # domain-reconstruction-1984/4-fill/ (2026-09-09)
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(HERE))
import HAT_insert_seaward_rows as ins  # noqa: E402
from hat_topo_version import dune_topo_root, insert_figures_dir  # noqa: E402
from hat_topo_version import require_version  # noqa: E402
from hat_figure_style import (apply_style, C, INK, INK_MUTED,  # noqa: E402
                              caption, elevation_cmap, figsize, save,
                              spines_for_image, _letter_inside, _title)

PRODUCT = "1984-start"
BASE = "v2"
# DELETED 2026-09-07 with every layer (only unmodified topography is kept);
# the literal is kept as the name of what this drew. require_version() in
# main() says so before any array is opened.
# (folder on disk, the fill rule it holds, the rule as a panel title). The
# version token names the LAYER, which is working vocabulary: it is carried
# here and in the caption's mapping, never on the canvas (2026-09-10).
LAYERS = [("v4", "measured + floor", "measured\n+ floor"),
          ("v5", "measured + median", "measured\n+ median"),
          ("v6", "flat platform", "flat\nplatform"),
          ("v7", "matched, crest kept", "matched,\ncrest kept"),
          ("v8", "matched, crest skipped", "matched,\ncrest skipped")]
CELL = 10.0
BERM_M = 1.34            # Barrier3D BermEl for these runs (m MHW); dune file is height above it
ROAD_ROWS = 2


def survey_in_z_frame(ext, dom, D):
    s = np.load(dune_topo_root(PRODUCT).parent / "1-extraction" / "npy-arrays_survey" / f"domain_{D}.npy").astype(float)
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
    apply_style()

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

    fig = plt.figure(figsize=figsize("double", aspect=0.80),
                     constrained_layout=True)
    gs = fig.add_gridspec(2, 5, height_ratios=(1.3, 1.0))
    gsb = gs[1, :].subgridspec(1, 5, wspace=0.12)
    a = fig.add_subplot(gs[0, 0:2])
    b = fig.add_subplot(gs[0, 2])
    c = fig.add_subplot(gs[0, 3:5])

    def road_band(ax, row_top, label=None, below=False):
        ax.add_patch(Rectangle((-0.5, row_top - 0.5), nal, ROAD_ROWS,
                               fill=False, ec=C["ROAD"], lw=1.4, zorder=8))
        if label:
            y = (row_top + ROAD_ROWS - 0.5 + 0.3) if below else (row_top - 0.8)
            ax.text(nal - 1, y, label, ha="right",
                    va="top" if below else "bottom", fontsize=7, color=INK,
                    zorder=9,
                    bbox=dict(fc="white", ec="none", alpha=.85, pad=1))

    # (a) the survey around the pick
    a.imshow(zg, cmap=cmap, norm=norm, aspect="auto", origin="upper",
             extent=(-0.5, nal - 0.5, span[-1] + 0.5, span[0] - 0.5))
    a.contourf(np.arange(nal), span, (yg == 2009).astype(float),
               levels=[0.5, 1.5], colors="none", hatches=["////"],
               extend="neither")
    a.axhline(-1, color=INK, lw=1.1, zorder=6)
    a.text(0.5, -3.2, "crest pick", fontsize=7, va="bottom", color=INK,
           bbox=dict(fc="white", ec="none", alpha=.85, pad=1))
    a.axhline(-0.5, color=INK, lw=.6, ls=":", zorder=6)
    a.text(0.5, 1.5, "row 0", fontsize=7, va="top", color=INK,
           bbox=dict(fc="white", ec="none", alpha=.85, pad=1))
    road_band(a, raw_sb / CELL, "NC-12 as measured")
    a.set_ylabel("cross-shore cell relative to row 0\n(negative = seaward)")
    a.set_xlabel("alongshore cell")
    a.set_ylim(span[-1] + 0.5, span[0] - 0.5)
    a.set_xticks([0, 25, 49])
    spines_for_image(a)
    a.set_title("the survey around the pick", loc="center")
    _letter_inside(a, 0)

    # (b) the domain as extracted, as the model sees it
    b.imshow(g2, cmap=cmap, norm=norm, aspect="auto", origin="upper",
             extent=(-0.5, nal - 0.5, n_show - 0.5, -2.5))
    b.axhline(-0.5, color=INK, lw=.8, zorder=6)
    b.text(nal - 1, -2.2, "dune rows", fontsize=7, va="top", ha="right",
           color="white")
    b.text(0.5, 0.1, "row 0", fontsize=7, va="top", color=INK,
           bbox=dict(fc="white", ec="none", alpha=.85, pad=1))
    road_band(b, k_b, "NC-12", below=True)
    b.set_ylabel("model row")
    b.set_xlabel("alongshore cell")
    b.set_xticks([0, 25, 49])
    spines_for_image(b)
    b.set_title("as extracted", loc="center")
    _letter_inside(b, 1)

    # (c) the same domain with the rows inserted
    c.imshow(g5, cmap=cmap, norm=norm, aspect="auto", origin="upper",
             extent=(-0.5, nal - 0.5, n_show + n - 0.5, -2.5))
    c.add_patch(Rectangle((-0.5, -0.5), nal, n, fill=False, ec=C["ADDED"],
                          lw=1.8, ls=(0, (4, 2)), zorder=7))
    c.text(0.5, 0.1, "rows added", fontsize=7, va="top", color=C["ADDED"],
           bbox=dict(fc="white", ec="none", alpha=.85, pad=1))
    c.axhline(n - 0.5, color=INK, lw=.8, zorder=6)
    c.text(nal - 1, n + 0.1, "old row 0", fontsize=7, va="top", ha="right",
           color=INK, bbox=dict(fc="white", ec="none", alpha=.85, pad=1))
    c.text(nal - 1, -2.2, "dune rows", fontsize=7, va="top", ha="right",
           color="white")
    road_band(c, k_c, "NC-12")
    c.set_ylabel("model row")
    c.set_xlabel("alongshore cell")
    c.set_xticks([0, 25, 49])
    spines_for_image(c)
    c.set_title("the rows inserted", loc="center")
    _letter_inside(c, 2)

    # bottom: the five fill rules, zoomed to the block
    zoom = n + 8
    bottom = []
    for k, (v, _rule, title) in enumerate(LAYERS):
        ax = fig.add_subplot(gsb[0, k])
        bottom.append(ax)
        g = model_grid(v, D, zoom)
        ax.imshow(g, cmap=cmap, norm=norm, aspect="auto", origin="upper",
                  extent=(-0.5, nal - 0.5, zoom - 0.5, -2.5))
        ax.add_patch(Rectangle((-0.5, -0.5), nal, n, fill=False,
                               ec=C["ADDED"], lw=1.4, ls=(0, (4, 2)), zorder=7))
        ax.axhline(n - 0.5, color=INK, lw=.8, zorder=6)
        ax.add_patch(Rectangle((-0.5, k_c - 0.5), nal, ROAD_ROWS, fill=False,
                               ec=C["ROAD"], lw=1.4, zorder=8))
        ax.set_xticks([0, 25, 49])
        if k == 2:
            ax.set_xlabel("alongshore cell")
        ax.set_title(title, loc="center", fontsize=8.5)
        _letter_inside(ax, 3 + k)
        spines_for_image(ax)
        if k == 0:
            ax.set_ylabel("model row")
        else:
            ax.set_yticklabels([])

    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=bottom,
                      orientation="horizontal", location="bottom",
                      boundaries=bounds[1:], ticks=bounds[1:-1],
                      shrink=0.42, aspect=34, pad=0.02)
    cb.outline.set_linewidth(0.6)
    cb.ax.tick_params(labelsize=7, length=2)
    cb.set_label("elevation (m MHW); leftmost class is water (≤ 0)",
                 fontsize=7.6, labelpad=2)

    caption(fig,
            "GIS {D}: the 1984 seaward-row insert drawn as Barrier3D grids, "
            "rows cross-shore and columns alongshore, colour the elevation "
            "class. (a) is the DEM in the extractor’s frame, each profile "
            "referenced to its own crest pick; the hatched cells are 2009 "
            "survey and the rest 1996, and the cells above the crest line are "
            "the beach and seaward dune face the extraction drops. NC-12 is "
            "at its measured offset, {sb:+.0f} m from row 0. (b) is what "
            "Barrier3D receives: the picked crest as two dune rows, then the "
            "interior from row 0, with the negative setback floored to zero so "
            "the road lands on rows {kb}–{kb1}. (c) prepends {n} rows at the "
            "dropped coordinates and puts the dune rows in front, so row 0 "
            "moves seaward by {n} cells while NC-12 stays on the ground it was "
            "measured on ({sb:+.0f} + {adv} = {aft:+.0f} m, rows {kc}–{kc1}); "
            "the old crest remains inside the domain as a ridge. It is drawn "
            "with the “measured + median” fill. (d)–(h) are the five fill "
            "rules for those {n} rows, built as dune-topo layers v4–v8 in that "
            "order; everything below the black line is identical in all five "
            "and identical to the domain as extracted. Elevations are metres "
            "above mean high water; one cell is 10 m."
            .format(D=D, n=n, sb=raw_sb, adv=int(10 * n),
                    aft=raw_sb + 10 * n, kb=k_b, kb1=k_b + 1,
                    kc=k_c, kc1=k_c + 1))

    # NOTE 2026-09-08: this figure now lives in figures/superseded-layers/; the script is
    # guarded (no layer on disk), so nothing is written here until a layer is rebuilt.
    out = Path(args.out) if args.out else insert_figures_dir(PRODUCT, "4-fill") / f"HAT_insert_explainer_grid_GIS{D}.png"
    written = save(fig, out, vector=False)
    print("wrote", written[0])


if __name__ == "__main__":
    main()
