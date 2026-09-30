#!/usr/bin/env python3
"""
The candidate interior fills as Barrier3D domain views: the grid the model would be handed under each.

    python scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/4-fill/HAT_plot_fill_options_grid.py [--domain 85] [--rows 26]

Guarded: the layers it draws were deleted 2026-09-07, so it stops
before drawing until one is rebuilt. Details: scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import argparse
import csv
import importlib.util as _iu
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle


# Walk up until a directory holds data/hatteras_init
def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer.hat_topo_version import array_name, dune_topo_root            # noqa: E402
from site_layer.hat_topo_version import insert_figures_dir  # noqa: E402
from site_layer.hat_topo_version import require_version  # noqa: E402
from site_layer.hat_figure_style import (apply_style, C, INK, caption,       # noqa: E402
                              elevation_cmap, figsize, save,
                              spines_for_image, _title)

# --- CONFIG ------------------------------------------------------------------
# INS_V supplies N and the post-insert setback ONLY
BASE_V, INS_V = "v2", "v4"   # base was "v3" and the layer "v5" until the 2026-09-04 renumber
OFFSET_SCRIPT = (REPO / "scripts/input_prep/4-mgmt-forcings/road_offset"
                 / "1-produce/HAT_road_offset_from_dune_start.py")
BERM_EL_M = 1.7
DUNE_ROWS = 2
ROAD_CELLS = 2
BACKDUNE_ROWS = 3
# -----------------------------------------------------------------------------

# THE ROAD ELEVATION CASCADE ACTUALLY USES, m MHW
from site_layer.hat_topo_version import ROAD_ELEVATION_FILE as ROAD_ELEV_FILE  # noqa: E402


# Per-domain road elevation, m MHW
def road_elevation(domain):
    rows = list(csv.reader(open(ROAD_ELEV_FILE)))
    ids = [int(float(x)) for x in rows[0]]
    vals = [float(x) for x in rows[1]]
    return dict(zip(ids, vals)).get(domain)


# The extractor module, for the 1984-start product
def load_ext():
    spec = _iu.spec_from_file_location("hat_off", OFFSET_SCRIPT)
    m = _iu.module_from_spec(spec)
    sys.modules["hat_off"] = m
    spec.loader.exec_module(m)
    return m.load_extractor("1984-start")


# The DEM cell at each added position, per profile
def real_cells(dom, row0, n, n_along):
    z = dom["z"]
    out = np.full((n, n_along), np.nan)
    for i in range(n_along):
        for k in range(n):
            src = int(row0[i]) - n + k
            if 0 <= src < z.shape[1]:
                out[k, i] = z[i, src]
    return out


# `road_ele` paints the road block at a SINGLE elevation, which is what the model holds
def draw(ax, topo, dune, setback, n_added, nrows, idx, title,
         road_ele=None):
    cmap, norm, _ = elevation_cmap()
    n_along = topo.shape[1]
    strip = np.tile(BERM_EL_M + dune[None, :n_along], (DUNE_ROWS, 1))
    shown = np.vstack([strip, topo[:nrows, :]])
    if road_ele is not None and setback >= 0:
        r0 = int(setback / 10.0) + DUNE_ROWS
        shown[r0:r0 + ROAD_CELLS, :] = road_ele
    ax.imshow(shown, cmap=cmap, norm=norm,
              aspect="auto", interpolation="nearest", zorder=1)

    ax.axhline(DUNE_ROWS - 0.5, color=INK, lw=0.9, zorder=4)
    if n_added:
        ax.add_patch(Rectangle((-0.5, DUNE_ROWS - 0.5), n_along, n_added,
                               fill=False, ec=C["ADDED"], lw=1.4,
                               ls=(0, (4, 2)), zorder=5))

    rs = int(setback / 10.0) + DUNE_ROWS
    ax.add_patch(Rectangle((-0.5, rs - 0.5), n_along, ROAD_CELLS, fill=False,
                           ec=C["ROAD"], lw=1.5, zorder=8))
    if setback < 0:
        ax.add_patch(Rectangle((-0.5, rs - 0.5), n_along, ROAD_CELLS,
                               facecolor="none", hatch="////", ec=C["ROAD"],
                               lw=0.0, zorder=7))

    ax.set_xlim(-0.5, n_along - 0.5)
    ax.set_ylim(nrows + DUNE_ROWS - 0.5, -0.5)
    ax.set_xticks([0, 25, 49])
    ax.set_yticks([0, 5, 10, 15, 20, 25])
    spines_for_image(ax)
    _title(ax, idx, title)


# Run: one grid panel per candidate fill
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=int, default=85)
    ap.add_argument("--rows", type=int, default=26)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    require_version("1984-start", INS_V, "INS_V, the layer that supplies N and the post-insert setback")
    D = args.domain
    apply_style()

    R = dune_topo_root("1984-start")
    v3 = np.load(R / BASE_V / "topography" / array_name("topography", D)) * 10.0
    dune = np.load(R / BASE_V / "dunes" / array_name("dune", D)) * 10.0
    aud = {int(r["domain"]): r for r in csv.DictReader(
        open(R / INS_V / "HAT_seaward_row_insert_audit.csv"))}
    n = int(aud[D]["n_rows_inserted"])
    setb = float(aud[D]["setback_model_after_m"])
    setb_meas = float(aud[D]["setback_raw_before_m"])

    ext = load_ext()
    dom = ext.load_profiles(ext.LOAD_PATH / "domain_{}.npy".format(D))
    prof = ext.masked_profiles(dom["z"])
    w = json.load(open(ext.PICKS_DIR / "HAT_dune_search_windows_{}.json"
                       .format(BASE_V)))["domain_{}".format(D)]
    _e, dl = ext.find_dunes(prof, dom["start_beach"], int(w["i0"]), int(w["i1"]))
    row0, _ = ext.interior_row0_line(prof, dl)
    n_along = v3.shape[1]

    plat = np.median(v3[1:1 + BACKDUNE_ROWS, :], axis=0)
    real = real_cells(dom, row0, n, n_along)
    dry = np.isfinite(real) & (real > 0.0)
    road_ele = road_elevation(D)
    if road_ele is None:
        raise SystemExit(
            "domain {} has no entry in {} - it carries no NC-12, so there is "
            "no road elevation to draw.".format(D, ROAD_ELEV_FILE.name))

    blocks = [
        # Matched backdune: the first n rows copied in front of themselves (0% measured)
        ("matched backdune", "copies the present near-dune profile to the "
         "1984 position: real cells, but measurements of a different place, "
         "so none of the block is taken from the DEM at these coordinates",
         v3[:n, :].copy(), 0.0),
        # KEEP DRY, FILL WET WITH THE BLOCK'S OWN MEDIAN
        ("measured + median", "keeps every dry cell as measured and fills only "
         "the cells at or below mean high water, with the median of the "
         "block’s own dry cells",
         np.where(dry, real, np.median(real[dry]) if dry.any() else plat.mean()),
         100.0 * dry.mean()),
        ("raw DEM (control)", "is the 1996 cells as they are, no floor and no "
         "dry-land test — a control, not a candidate, drawn to show what the "
         "two guards reject",
         np.where(np.isfinite(real), real, plat[None, :]),
         100.0 * np.isfinite(real).mean()),
    ]

    fig, axes = plt.subplots(2, 2, figsize=figsize("double", aspect=0.95),
                             sharex=True, sharey=True,
                             constrained_layout=True)
    axf = axes.ravel()
    rs_ref = int(setb_meas / 10.0)
    draw(axf[0], v3, dune, setb_meas, 0, args.rows, 0, "domain as extracted")

    rs = int(setb / 10.0)
    notes = []
    for k, (lab, why, blk, pct) in enumerate(blocks, start=1):
        full = np.vstack([blk, v3])
        under = np.median(full[rs:rs + ROAD_CELLS, :], axis=1)
        notes.append(
            "({}) “{}” {}; {:.0f}% of the block comes from the DEM, its "
            "mean elevation is {:.2f} m, and the ground under the road was "
            "{:.2f} / {:.2f} m before bulldoze() flattened it."
            .format("abcd"[k], lab, why, pct, float(np.mean(blk)),
                    under[0], under[1]))
        draw(axf[k], full, dune, setb, n, args.rows, k, lab,
             road_ele=road_ele)

    for ax in axes[-1]:
        ax.set_xlabel("alongshore cell")
    for ax in axes[:, 0]:
        ax.set_ylabel("cross-shore row\n(0 = first interior cell)")

    cmap, norm, bounds = elevation_cmap()
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm),
                      ax=axes.ravel().tolist(), orientation="horizontal",
                      location="bottom", boundaries=bounds[1:],
                      ticks=bounds[1:-1], shrink=0.45, aspect=32, pad=0.015)
    cb.outline.set_linewidth(0.6)
    cb.ax.tick_params(labelsize=7, length=2)
    cb.set_label("elevation (m MHW); leftmost class is water (≤ 0)",
                 fontsize=7.6, labelpad=2)

    fig.legend(handles=[
        Line2D([0], [0], color=C["ROAD"], lw=1.5,
               label="NC-12, at its measured offset and road elevation"),
        Patch(facecolor="white", edgecolor=C["ROAD"], hatch="////", lw=0.6,
              label="hatched: seaward of interior row 0"),
        Line2D([0], [0], color=C["ADDED"], lw=1.4, ls=(0, (4, 2)),
               label="rows added behind the dune"),
    ], loc="outside lower center", ncol=3, frameon=False, fontsize=7.6)

    caption(fig,
            "GIS {D}, the candidate interior fills as the Barrier3D domain the "
            "model would be handed: two dune rows over the interior domain, "
            "row 0 first. Panels (b)–(d) differ ONLY in the {n} added rows; "
            "everything landward of the dashed band is identical in all three "
            "and identical to (a) from its row 0 on. NC-12 is drawn at its "
            "measured offset AND at the elevation bulldoze() gives it — one "
            "value, {re:.2f} m, for every road cell, read from "
            "RoadElevation.csv and independent of the fill. The road does not "
            "move; the added rows move interior row 0 out from under it. "
            "(a) is the domain as extracted, where the measured offset of "
            "{sm:+.0f} m puts the road on row {rr}, SEAWARD of interior row 0 "
            "(hatched). That is the failure, not a plotting artefact: int() "
            "truncates toward zero, so roadway_manager would index the "
            "interior from its landward end and bulldoze the sound-side "
            "marsh, and the domain would not initialise. {notes}"
            .format(D=D, n=n, re=road_ele, sm=setb_meas, rr=rs_ref,
                    notes=" ".join(notes)))

    out = Path(args.out) if args.out else (
        # Superseded-layers figure: guarded, writes nothing until a layer is rebuilt
        insert_figures_dir("1984-start", "4-fill")
        / "HAT_fill_options_grid_GIS{}.png".format(D))
    written = save(fig, out, vector=False)
    print("wrote {}".format(written[0]))


if __name__ == "__main__":
    main()
