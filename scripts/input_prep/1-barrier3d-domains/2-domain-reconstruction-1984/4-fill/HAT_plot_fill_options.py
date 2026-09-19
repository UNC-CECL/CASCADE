#!/usr/bin/env python3
r"""
HAT_plot_fill_options.py
==============================================================================
What can the added interior rows be MADE of? Four candidates and a control,
one domain.

THE PROBLEM
    The rows added behind the dune stand where land existed in 1984 and had
    eroded away by 1996. No survey covers that ground. The DEM does have cells
    at those coordinates, but they are the 1996 surface -- a later, lower
    landform at the same place. So every option below is a different answer to
    "what was here in 1984", and none of them is a measurement of it.

FOUR CANDIDATES AND A CONTROL
    Each is named by its RULE, in the legend, on the bar axes and in the
    caption -- the b/c/d/e tags this figure carried until 2026-09-10 read as
    panel letters beside the house style's own (a)/(b)/(c), which is a trap
    worth closing. The grid figure names the same rules the same way.

    flat backdune     flat, at the median of interior rows 1-3
    matched backdune  today's near-dune PROFILE copied to the 1984 position
                     -- "the 1984 backdune looked like the present one, just
                     further seaward". Real cells, so it carries alongshore
                     texture the flat fill cannot.
    measured + floor  the real DEM cell where it is dry land, floored at the
                     backdune platform. THIS IS THE SHIPPED RULE -- v4 at the
                     ten block domains, v5 island-wide. Same rule, new scope.
    measured + median  every dry cell kept AS MEASURED; only the cells at or
                     below MHW filled, with the median of the block's own dry
                     cells. `--fill median`. One guard, one constant, and no
                     measurement is ever raised -- 91% measured at GIS 85
                     against the floored rule's 47%.
    raw DEM          the 1996 cells as they are, no floor, no dry-land test.
                     A CONTROL, NOT A CANDIDATE -- drawn to show what the two
                     guards actually reject, in numbers.

    Dropped 2026-09-03: `taper`, a linear platform-to-row-0 ramp. Fully
    invented, and it anchored on row 0 -- which at GIS 85 IS the mis-picked
    1996 crest, so it inherited a known-bad endpoint. `--fill taper` still
    exists in HAT_insert_seaward_rows.py: removing a build capability is a
    different decision from removing a figure panel.

    `matched backdune` is NOT a --fill choice in HAT_insert_seaward_rows.py.
    It is drawn as a candidate; building it needs a new fill rule.

    Still undrawn: an alongshore analogue from a neighbouring domain, and a
    mass-conservative reconstruction. Neither is a one-line variant of the
    others, and the second is the only one that would be DERIVED rather than
    asserted.

WHAT TO LOOK FOR
    Where NC-12 lands. The road is a fixed 2-cell block at a fixed setback, so
    the only thing that changes between options is the ground under it.
    `measured + floor` puts it on 3.17 / 4.96 m because those cells are the
    1996 DUNE FACE; the flat and matched backdunes put it on backdune, where a
    road behind a dune belongs.

USAGE
    python HAT_plot_fill_options.py [--domain 85]
==============================================================================
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


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer.hat_topo_version import array_name, dune_topo_root          # noqa: E402
from site_layer.hat_topo_version import insert_figures_dir  # noqa: E402
from site_layer.hat_topo_version import require_version  # noqa: E402
from site_layer.hat_figure_style import (apply_style, C, INK, INK_MUTED,  # noqa: E402
                              caption, figsize, open_frame, save, _title)

# INS_V supplies N and the post-insert setback ONLY; the figure draws
# blocks it builds itself. Repointed v4 -> v5 on 2026-09-03: N is
# identical at all ten block domains (verified), and v5 is the version
# taken forward, so v4 no longer has to exist for this figure to build.
# DELETED 2026-09-07 with every layer (only unmodified topography is kept);
# the literal is kept as the name of what this drew. require_version() in
# main() says so before any array is opened.
BASE_V, INS_V = "v2", "v4"   # base was "v3" and the layer "v5" until the 2026-09-04 renumber
OFFSET_SCRIPT = (REPO / "scripts/input_prep/4-mgmt-forcings/road_offset"
                 / "1-produce/HAT_road_offset_from_dune_start.py")
BACKDUNE_ROWS = 3
ROAD_CELLS = 2


def load_ext():
    spec = _iu.spec_from_file_location("hat_off", OFFSET_SCRIPT)
    m = _iu.module_from_spec(spec)
    sys.modules["hat_off"] = m
    spec.loader.exec_module(m)
    return m.load_extractor("1984-start")


def real_cells(ext, dom, row0, n, n_along):
    """The DEM cell at each added position, per profile. NaN where off-array."""
    z = dom["z"]
    out = np.full((n, n_along), np.nan)
    for i in range(n_along):
        for k in range(n):
            src = int(row0[i]) - n + k
            if 0 <= src < z.shape[1]:
                out[k, i] = z[i, src]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=int, default=85)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    require_version("1984-start", INS_V, "INS_V, the layer that supplies N and the post-insert setback")
    D = args.domain
    apply_style()

    R = dune_topo_root("1984-start")
    v3 = np.load(R / BASE_V / "topography" / array_name("topography", D)) * 10.0
    aud = {int(r["domain"]): r for r in csv.DictReader(
        open(R / INS_V / "HAT_seaward_row_insert_audit.csv"))}
    n = int(aud[D]["n_rows_inserted"])
    setb = float(aud[D]["setback_model_after_m"])

    ext = load_ext()
    dom = ext.load_profiles(ext.LOAD_PATH / "domain_{}.npy".format(D))
    prof = ext.masked_profiles(dom["z"])
    wj = json.load(open(ext.PICKS_DIR / "HAT_dune_search_windows_{}.json"
                        .format(BASE_V)))
    w = wj["domain_{}".format(D)]
    _e, dl = ext.find_dunes(prof, dom["start_beach"], int(w["i0"]), int(w["i1"]))
    row0, _ = ext.interior_row0_line(prof, dl)
    n_along = v3.shape[1]

    plat = np.median(v3[1:1 + BACKDUNE_ROWS, :], axis=0)          # (along,)
    real = real_cells(ext, dom, row0, n, n_along)                  # (n, along)

    # --- the candidate blocks, all (n, n_along) in m MHW ------------------
    flat = np.repeat(plat[None, :], n, axis=0)
    matched = v3[:n, :].copy()          # today's near-dune profile, moved
    measured = np.where(np.isfinite(real) & (real > 0.0),
                        np.maximum(real, plat[None, :]), plat[None, :])
    raw = np.where(np.isfinite(real), real, plat[None, :])

    dry = np.isfinite(real) & (real > 0.0)

    # TAGGED b/c/d/e TO MATCH THE GRID FIGURE'S PANEL LETTERS, where (a) is the
    # reference. The two figures are read side by side and single letters that
    # meant different things in each was a trap worth closing.
    #
    # THE "% FROM THE DEM" IS CARRIED HERE, beside the block it describes.
    # It used to be recomputed further down by an `if tag == "A"` chain, which
    # fell through to the wrong branch the moment these tags were renamed and
    # silently reported 47% for every option. A number that describes a block
    # belongs with the block.
    #
    # `matched backdune` counts as 0%: its cells ARE real measurements, but of
    # a different place, copied. Measured-ness here means "measured AT THESE
    # COORDINATES", which is the only sense that bears on whether the fill is
    # invented.
    opts = [
        # (key, legend label, axis label, block, colour, % taken from the DEM)
        #
        # Five house colours, one per rule; the same colour names the rule in
        # all three panels. "measured + floor" is the rule the layers were
        # built with, so it carries C["ACCENT"], the change under test.
        ("flat", "flat backdune", "flat\nbackdune", flat, C["LATE"], 0.0),
        ("matched", "matched backdune", "matched\nbackdune", matched,
         C["REF"], 0.0),
        ("floor", "measured + floor", "measured\n+ floor", measured,
         C["ACCENT"], 100.0 * (dry & (real >= plat[None, :])).mean()),
        ("median", "measured + median", "measured\n+ median",
         np.where(dry, real, np.median(real[dry]) if dry.any() else plat.mean()),
         C["EARLY"], 100.0 * dry.mean()),
        ("raw", "raw DEM (control)", "raw DEM\n(control)", raw, C["BASE"],
         100.0 * np.isfinite(real).mean()),
    ]
    cols = [o[4] for o in opts]
    ticks = [o[2] for o in opts]

    fig = plt.figure(figsize=figsize("double", aspect=0.88),
                     constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0])

    # ---- (a) the five profiles in a common frame ------------------------
    ax = fig.add_subplot(gs[0, :])
    x_body = np.arange(v3.shape[0]) * 10.0
    keep = x_body <= 260
    ax.axvspan(-n * 10 - 5, -5, color=C["ADDED_FILL"], zorder=0)
    ax.text(-n * 10 - 12, 5.6, "rows added", fontsize=7.5, color=INK_MUTED,
            va="top", ha="left")
    ax.plot(x_body[keep], np.median(v3, axis=1)[keep], "-", color=INK,
            lw=2.4, label="interior, unchanged by any option", zorder=2)
    xa = (np.arange(n) - n) * 10.0
    for _k, lab, _t, blk, col, _pct in opts:
        med = np.median(blk, axis=1)
        ax.plot(np.r_[xa, 0.0], np.r_[med, np.median(v3[0])], "-o", color=col,
                lw=1.3, ms=2.4, label=lab, zorder=4)
    ax.add_patch(Rectangle((setb - n * 10, -0.78), ROAD_CELLS * 10, 0.40,
                           color=C["ROAD"], zorder=7))
    ax.annotate("NC-12", xy=(setb - n * 10 + 10, -0.38), xytext=(45, -0.74),
                fontsize=7.5, color=INK, va="center",
                arrowprops=dict(arrowstyle="->", lw=0.8, color=INK))
    ax.axvline(0, color=INK, lw=1.0)
    ax.text(5, 5.6, "row 0 as extracted", fontsize=7.5, rotation=90,
            va="top", color=INK_MUTED)
    ax.axhline(0, color=C["WATER"], lw=0.9, ls=":")
    ax.set_xlim(-n * 10 - 18, 260)
    ax.set_ylim(-0.9, 6.7)
    ax.set_xlabel("metres landward of interior row 0 as extracted")
    ax.set_ylabel("median elevation (m MHW)")
    ax.grid(axis="y")
    open_frame(ax)
    ax.legend(loc="upper right", ncol=3, fontsize=7.5, handlelength=1.4,
              columnspacing=1.0, borderpad=0.35)
    _title(ax, 0, "candidate fills in profile")

    # ---- (b) what NC-12 ends up sitting on ------------------------------
    axb = fig.add_subplot(gs[1, 0])
    rs = int(setb / 10.0)
    unders = []
    for _k, _lab, _t, blk, _col, _pct in opts:
        full = np.vstack([blk, v3])
        unders.append(np.median(full[rs:rs + ROAD_CELLS, :], axis=1))
    xi = np.arange(len(opts))
    u = np.array(unders)
    axb.bar(xi - 0.19, u[:, 0], 0.38, color=cols)
    axb.bar(xi + 0.19, u[:, 1], 0.38, color=cols, alpha=0.5)
    axb.axhline(float(np.median(plat)), color=C["REF"], ls="--", lw=1.0,
                zorder=3)
    axb.set_xticks(xi)
    axb.set_xticklabels(ticks, fontsize=7)
    axb.set_ylim(0, 6.4)
    axb.set_ylabel("elevation under NC-12 (m MHW)")
    axb.grid(axis="y")
    open_frame(axb)
    axb.legend(handles=[
        Patch(facecolor=INK_MUTED, label="seaward road cell"),
        Patch(facecolor=INK_MUTED, alpha=0.5, label="landward road cell"),
        Line2D([0], [0], color=C["REF"], ls="--", lw=1.0,
               label="backdune platform"),
    ], loc="upper left", fontsize=7, handlelength=1.4, borderpad=0.35)
    _title(axb, 1, "ground under NC-12")

    # ---- (c) how much of each block is invented -------------------------
    axc = fig.add_subplot(gs[1, 1])
    frac = [o[5] for o in opts]
    axc.bar(xi, frac, 0.62, color=cols)
    axc.set_xticks(xi)
    axc.set_xticklabels(ticks, fontsize=7)
    axc.set_ylim(0, 108)
    axc.set_yticks([0, 25, 50, 75, 100])
    axc.set_ylabel("% taken from the DEM")
    axc.grid(axis="y")
    open_frame(axc)
    _title(axc, 2, "share taken from the DEM")

    caption(fig,
            "GIS {D}: what the {n} interior rows added behind the dune should "
            "be made of. No survey covers ground that had eroded away by 1996, "
            "so every option is an assertion about 1984, not a measurement of "
            "it. The cells the DEM does hold at those coordinates are the 1996 "
            "surface — a later landform in the same place — which is why "
            "“measured + floor”, the rule the layers were built with, puts "
            "NC-12 on a {ra:.1f}/{rb:.1f} m dune FACE rather than on backdune; "
            "flat and matched backdune put the road where a road behind a dune "
            "belongs. “matched backdune” copies the present near-dune profile "
            "to the 1984 position: real cells, but measurements of a different "
            "place, so they count as 0% taken from the DEM. "
            "“measured + median” keeps every dry cell as measured and fills "
            "only the cells at or below mean high water with the median of the "
            "block’s own dry cells. “raw DEM” is a control, not a candidate: it "
            "imports sub-MHW beach, and is drawn to show what the two guards of "
            "“measured + floor” reject. The dashed reference in (b) is the "
            "backdune platform, the median of interior rows 1–3, at "
            "{plat:.2f} m; the road block is fixed in place, so the only thing "
            "that changes between options is the ground under it."
            .format(D=D, n=n, ra=u[2, 0], rb=u[2, 1],
                    plat=float(np.median(plat))))

    out = Path(args.out) if args.out else (
        # NOTE 2026-09-08: this figure now lives in figures/superseded-layers/; the script is
        # guarded (no layer on disk), so nothing is written here until a layer is rebuilt.
        insert_figures_dir("1984-start", "4-fill")
        / "HAT_fill_options_GIS{}.png".format(D))
    written = save(fig, out)
    print("wrote {}".format(written[0]))
    print("\n  option   road sits on        mean added elev   % from DEM")
    for (key, _lab, _t, blk, _col, _p), un, fr in zip(opts, unders, frac):
        print("  {:<9}  {:>5.2f} / {:>5.2f} m      {:>6.2f} m        {:>5.0f}%"
              .format(key, un[0], un[1], float(np.mean(blk)), fr))


if __name__ == "__main__":
    main()
