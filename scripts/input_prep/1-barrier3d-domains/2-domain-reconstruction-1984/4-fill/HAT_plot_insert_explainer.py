#!/usr/bin/env python3
r"""
HAT_plot_insert_explainer.py
==============================================================================
How the 1984 seaward-row insert is built, on one cross-shore line at one domain
(GIS 85 by default), drawn from the real arrays.

    (a) the survey along that line, cell by cell, shaded by survey year, with
        the crest pick and the measured 1984 road position
    (b) what the extraction keeps: the crest cell becomes Barrier3D's dune
        rows, everything landward of it is the interior (row 0 first), and
        everything seaward is dropped
    (c) the insert: N rows are added at the dropped coordinates, the dune rows
        move in front of them, the old crest stays inside, and the road does
        not move -- its setback grows by 10 N
    (d) what each version writes into the N cells, over the survey values that
        are there

Alongshore MEDIANS of the domain's 50 profiles; the cross-shore axis is metres
from v2's interior row 0 (negative = seaward), the frame every version shares.

USAGE
    python HAT_plot_insert_explainer.py [--domain 85]
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
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]                  # domain-reconstruction-1984/4-fill/ (2026-09-09)
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(HERE))
import HAT_insert_seaward_rows as ins  # noqa: E402
from hat_topo_version import dune_topo_root, insert_figures_dir  # noqa: E402
from hat_topo_version import require_version  # noqa: E402

PRODUCT = "1984-start"
BASE = "v2"
# DELETED 2026-09-07 with every layer (only unmodified topography is kept);
# the literal is kept as the name of what this drew. require_version() in
# main() says so before any array is opened.
LAYERS = [("v4", "v4  measured + floor", "#2166ac", "-"),
          ("v5", "v5  measured + median", "#b2182b", "-"),
          ("v6", "v6  flat platform", "#1b7837", "--"),
          ("v7", "v7  matched, crest kept", "#e08214", "-"),
          ("v8", "v8  matched, crest skipped", "#7b3294", "--")]
YEAR_COLOUR = {1996: "#e8c98a", 2009: "#c9c9c9", 2014: "#a9b8c9", 0: "white"}
CELL = 10.0


def survey_in_z_frame(ext, dom, D):
    """Survey year per cell in the extractor's z frame (n_along, n_cross)."""
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
    r0, _ = ext.interior_row0_line(prof, dl)
    r0 = np.asarray(r0).astype(int)
    survey = survey_in_z_frame(ext, dom, D)

    # everything referenced to each profile's own row 0, then alongshore median
    span = np.arange(-12, 26)                     # cells relative to row 0
    z_rel = np.full((len(r0), span.size), np.nan)
    y_rel = np.zeros((len(r0), span.size), int)
    for i in range(len(r0)):
        if r0[i] < 0:
            continue
        for j, k in enumerate(span):
            src = r0[i] + k
            if 0 <= src < prof.shape[1]:
                z_rel[i, j] = prof[i, src]
                y_rel[i, j] = survey[i, src]
    zmed = np.nanmedian(z_rel, axis=0)
    ymode = np.array([np.bincount(col[col >= 0]).argmax() if (col >= 0).any() else 0
                      for col in y_rel.T])
    x = span * CELL                                # metres, row 0 at 0

    # the layers' added cells and v2's interior, alongshore medians
    audit = ins.csv.DictReader(open(dune_topo_root(PRODUCT) / "v4" / "HAT_seaward_row_insert_audit.csv"))
    n = next(int(r["n_rows_inserted"]) for r in audit if int(r["domain"]) == D)
    base = np.load(dune_topo_root(PRODUCT) / BASE / "topography" / f"domain_{D}_topography.npy") * CELL
    base_med = np.median(base[:26], axis=1)
    blocks = {}
    for v, *_ in LAYERS:
        a = np.load(dune_topo_root(PRODUCT) / v / "topography" / f"domain_{D}_topography.npy") * CELL
        blocks[v] = np.median(a[:n], axis=1)

    # road: measured 1984 setback (raw) from v2's own audit column
    import csv
    raw_sb = None
    for r in csv.DictReader(open(dune_topo_root(PRODUCT) / "v4" / "HAT_seaward_row_insert_audit.csv")):
        if int(r["domain"]) == D and r["setback_raw_before_m"] not in ("", "nan"):
            raw_sb = float(r["setback_raw_before_m"])
    road_x = raw_sb if raw_sb is not None else None       # metres from v2 row 0

    # ------------------------------------------------------------------ figure
    fig, axes = plt.subplots(2, 2, figsize=(15, 9.5))
    (a, b), (c, d) = axes
    crest_x = -CELL                                        # the crest cell is r0 - 1
    xlim = (-130, 250)

    def cells(ax, xs, zs, colours, hatch=None, alpha=1.0, edge="0.3"):
        for xi, zi, col in zip(xs, zs, colours):
            if np.isnan(zi):
                continue
            ax.add_patch(Rectangle((xi - CELL / 2, min(0, zi)), CELL, abs(zi), facecolor=col,
                                   edgecolor=edge, lw=.4, hatch=hatch, alpha=alpha))

    def road(ax, x_left, label, y=-0.75):
        """A 2-cell road bar with its seaward edge at x_left (m)."""
        if road_x is None:
            return
        ax.add_patch(Rectangle((x_left, y - 0.2), 2 * CELL, 0.4, facecolor="k"))
        ax.text(x_left + CELL, y - 0.3, label, ha="center", va="top", fontsize=8)

    for ax in (a, b, c, d):
        ax.set_xlim(*xlim)
        ax.set_ylim(-1.75, 6.6)
        ax.axhline(0, color="#4a90c2", lw=.8)
        ax.set_xlabel("cross-shore, metres from v2 interior row 0  (negative = seaward)")
        ax.set_ylabel("elevation (m MHW)")
        ax.grid(alpha=.25)

    # (a) the survey
    cells(a, x, zmed, [YEAR_COLOUR.get(int(y), "white") for y in ymode])
    a.plot(crest_x, zmed[span == -1][0] + 0.25, "v", color="k", ms=9)
    a.text(crest_x, zmed[span == -1][0] + 0.55, "crest pick", ha="center", fontsize=9)
    a.text(-95, 0.35, "beach", fontsize=9, color="0.3")
    a.text(-45, 2.9, "seaward\ndune face", fontsize=9, color="0.3", ha="center")
    a.text(90, 2.3, "backdune", fontsize=9, color="0.3")
    for yr, xx in ((1996, 30), (2009, 200)):
        a.add_patch(Rectangle((xx, 5.5), 18, 0.5, facecolor=YEAR_COLOUR[yr], edgecolor="0.3", lw=.4))
        a.text(xx + 22, 5.75, f"{yr} survey", va="center", fontsize=8.5)
    road(a, road_x, f"NC-12, measured {raw_sb:+.0f} m from row 0")
    a.set_title(f"(a) the survey along GIS {D}, cell by cell (alongshore median)", loc="left", fontsize=10.5)

    # (b) what v2 keeps
    keep = span >= 0
    cells(b, x[keep], zmed[keep], ["#d9b77a"] * keep.sum())
    cells(b, x[~keep & (span != -1)], zmed[~keep & (span != -1)], ["white"] * ((~keep & (span != -1)).sum()),
          hatch="///", edge="0.6")
    b.add_patch(Rectangle((crest_x - CELL / 2, 0), CELL, zmed[span == -1][0], facecolor="#7a5c2e", edgecolor="k", lw=.8))
    b.text(crest_x, zmed[span == -1][0] + 0.3, "crest cell →\nBarrier3D dune rows", ha="center", fontsize=8.5)
    b.annotate("interior row 0", xy=(0, base_med[0]), xytext=(45, 5.6), fontsize=9,
               arrowprops=dict(arrowstyle="->", lw=.8))
    b.text(-95, 4.6, "dropped:\nBarrier3D grows\nits own beach", fontsize=8.5, color="0.4", ha="center")
    k_b = max(int(raw_sb / CELL), 0)
    road(b, (k_b - 0.5) * CELL, f"setback {raw_sb:+.0f} m floored to 0 → road on rows {k_b}–{k_b+1}")
    b.set_title("(b) what the v2 extraction keeps", loc="left", fontsize=10.5)

    # (c) the insert geometry
    cells(c, x[keep], zmed[keep], ["#d9b77a"] * keep.sum())
    blk = (span >= -n) & (span < 0)
    cells(c, x[blk], zmed[blk], ["white"] * blk.sum(), edge="0.6", hatch="///")
    c.add_patch(Rectangle((-n * CELL - CELL / 2, 0), n * CELL, 6.0, fill=False, edgecolor="#d17f00", lw=2, ls="--"))
    c.text(-n * CELL / 2 - CELL / 2, 6.15, f"{n} added rows = new interior rows 0–{n-1}", ha="center", fontsize=9, color="#d17f00")
    c.add_patch(Rectangle((-(n + 1) * CELL - CELL / 2, 0), CELL, zmed[span == -1][0], facecolor="#7a5c2e", edgecolor="k", lw=.8))
    c.text(-(n + 1) * CELL - 12, 2.6, "dune rows\nmoved\nin front", ha="right", fontsize=8.5)
    c.annotate("old crest stays,\nnow a ridge inside", xy=(crest_x, zmed[span == -1][0]), xytext=(60, 5.4),
               fontsize=8.5, arrowprops=dict(arrowstyle="->", lw=.8))
    k_c = int((raw_sb + 10 * n) / CELL)
    road(c, (k_c - n - 0.5) * CELL,
         f"setback {raw_sb:+.0f} + {10*n} = {raw_sb + 10*n:+.0f} m → rows {k_c}–{k_c+1} of the new interior")
    c.set_title("(c) the insert: rows added at the dropped coordinates", loc="left", fontsize=10.5)

    # (d) the fills
    xb = (np.arange(n) - n) * CELL                          # x of the added cells
    d.step(np.append(xb, 0) - CELL / 2, np.append(zmed[blk], zmed[blk][-1]), where="post",
           color="0.45", lw=1.2, ls=":", label="1996 survey at those cells")
    xi = np.arange(0, 26) * CELL
    d.step(np.append(xi, xi[-1] + CELL) - CELL / 2, np.append(base_med, base_med[-1]), where="post",
           color="0.2", lw=1.4, label="v2 interior (shared by every version)")
    for v, label, colour, ls in LAYERS:
        d.step(np.append(xb, 0) - CELL / 2, np.append(blocks[v], blocks[v][-1]), where="post",
               color=colour, lw=1.8, ls=ls, label=label)
    if road_x is not None:
        rs = int((raw_sb + 10 * n) / CELL)
        d.axvspan(-n * CELL + rs * CELL - CELL / 2, -n * CELL + (rs + 2) * CELL - CELL / 2,
                  color="k", alpha=0.08)
        d.text(-n * CELL + (rs + 1) * CELL - CELL / 2, 6.1, f"road, rows {rs}–{rs+1}", ha="center", fontsize=8.5)
    d.axvline(-CELL / 2, color="k", lw=.6)
    d.set_xlim(-n * CELL - 25, 120)
    d.legend(fontsize=8, loc="upper right")
    d.set_title(f"(d) what each version writes into the {n} cells", loc="left", fontsize=10.5)

    fig.suptitle(f"GIS {D}  ·  how the 1984 seaward-row insert is built, and what fills it",
                 fontsize=13, fontweight="bold", x=0.02, ha="left")
    fig.text(0.02, 0.005,
             "Cells shaded by the survey year the DEM took them from (a). The extraction (b) turns the picked crest cell into "
             "Barrier3D's dune rows and keeps everything landward as the interior; the seaward cells are dropped. The insert (c) "
             "adds N = round(1984–1997 dune-line shift / 10 m) rows at those dropped coordinates, so interior row 0 moves seaward and "
             "the road's setback grows by 10·N while the road itself stays where it was measured. The versions (d) differ only in "
             "what those N cells hold: the survey's own values (v4 floors the low ones to the backdune platform, v5 keeps them and "
             "fills water with the block's dry median), a flat platform (v6), or the interior's own near-dune rows copied seaward "
             "(v7 from the crest, v8 from one row behind it).",
             fontsize=8, wrap=True, va="bottom")
    fig.tight_layout(rect=(0, 0.07, 1, 0.96))
    # NOTE 2026-09-08: this figure now lives in figures/superseded-layers/; the script is
    # guarded (no layer on disk), so nothing is written here until a layer is rebuilt.
    out = Path(args.out) if args.out else insert_figures_dir(PRODUCT, "4-fill") / f"HAT_insert_explainer_GIS{D}.png"
    fig.savefig(out, dpi=150)
    print("wrote", out)


if __name__ == "__main__":
    main()
