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
    Tagged b/c/d/e to match the grid figure's panel letters, where (a) is
    the reference. The two are read side by side.

    b  flat backdune     flat, at the median of interior rows 1-3
    c  matched backdune  today's near-dune PROFILE copied to the 1984
                     position -- "the 1984 backdune looked like the present
                     one, just further seaward". Real cells, so it carries
                     alongshore texture the flat fill cannot.
    d  measured+floor  the real DEM cell where it is dry land, floored at the
                     backdune platform. THIS IS THE SHIPPED RULE -- v4 at the
                     ten block domains, v5 island-wide. Same rule, new scope.
    e  measured+median  every dry cell kept AS MEASURED; only the cells at or
                     below MHW filled, with the median of the block's own dry
                     cells. `--fill median`. One guard, one constant, and no
                     measurement is ever raised -- 91% measured at GIS 85
                     against (d)'s 47%.
    f  raw DEM       the 1996 cells as they are, no floor, no dry-land test.
                     A CONTROL, NOT A CANDIDATE -- drawn to show what (d)'s
                     two guards actually reject, in numbers.

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
    the only thing that changes between options is the ground under it. (d)
    puts it on 3.17 / 4.96 m because those cells are the 1996 DUNE FACE;
    (b) and (c) put it on backdune, where a road behind a dune belongs.

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
from matplotlib.patches import Rectangle


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
from hat_topo_version import array_name, dune_topo_root          # noqa: E402
from hat_topo_version import insert_figures_dir  # noqa: E402
from hat_figure_style import apply_style, C, caption, panel_title  # noqa: E402

# INS_V supplies N and the post-insert setback ONLY; the figure draws
# blocks it builds itself. Repointed v4 -> v5 on 2026-09-03: N is
# identical at all ten block domains (verified), and v5 is the version
# taken forward, so v4 no longer has to exist for this figure to build.
BASE_V, INS_V = "v3", "v5"
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
        ("b", "flat backdune", flat, "#4c72b0", 0.0),
        ("c", "matched backdune  [today's profile, moved seaward]",
         matched, "#55a868", 0.0),
        ("d", "measured + floor  [shipped: v4, v5]", measured, C["ACCENT"],
         100.0 * (dry & (real >= plat[None, :])).mean()),
        ("e", "measured + median fill  [--fill median]",
         np.where(dry, real, np.median(real[dry]) if dry.any() else plat.mean()),
         "#c44e52", 100.0 * dry.mean()),
        ("f", "raw DEM, no floor  [control]", raw, "#8172b3",
         100.0 * np.isfinite(real).mean()),
    ]

    fig = plt.figure(figsize=(11.5, 8.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0],
                          width_ratios=[1.55, 1.0], hspace=0.42, wspace=0.26)

    # ---- (a) the five profiles in a common frame ------------------------
    ax = fig.add_subplot(gs[0, :])
    x_body = np.arange(v3.shape[0]) * 10.0
    keep = x_body <= 260
    ax.axvspan(-n * 10 - 5, -5, color=C["ADDED_FILL"], zorder=0)
    ax.text(-n * 10, 5.55, "added rows", fontsize=8, color="#7a5200", va="top")
    ax.plot(x_body[keep], np.median(v3, axis=1)[keep], "-", color="#444444",
            lw=4.5, label="v3 interior (unchanged by any option)", zorder=2)
    xa = (np.arange(n) - n) * 10.0
    for tag, lab, blk, col, _pct in opts:
        med = np.median(blk, axis=1)
        ax.plot(np.r_[xa, 0.0], np.r_[med, np.median(v3[0])], "-o", color=col,
                lw=1.7, ms=3.5, label="{}  {}".format(tag, lab), zorder=4)
    ax.add_patch(Rectangle((setb - n * 10, -0.75), ROAD_CELLS * 10, 0.45,
                           color="k", zorder=7))
    ax.annotate("NC-12", xy=(setb - n * 10 + 10, -0.3), xytext=(70, -0.62),
                fontsize=8.5, arrowprops=dict(arrowstyle="->", lw=1.0))
    ax.axvline(0, color="k", lw=1.1)
    ax.text(3, 5.55, "v3 row 0", fontsize=8, rotation=90, va="top")
    ax.axhline(0, color="#3b6ea5", lw=0.8, ls=":")
    ax.set_xlim(-n * 10 - 18, 260)
    ax.set_ylim(-0.9, 5.8)
    ax.set_xlabel("metres landward of v3 interior row 0")
    ax.set_ylabel("median elevation (m MHW)")
    ax.set_title(panel_title(
        "a", "GIS {} — four candidates and a control for the {} added rows"
             .format(D, n)))
    ax.legend(fontsize=8, loc="upper right", ncol=2)
    ax.grid(alpha=.3)

    # ---- (b) what NC-12 ends up sitting on ------------------------------
    axb = fig.add_subplot(gs[1, 0])
    rs = int(setb / 10.0)
    names, unders = [], []
    for tag, lab, blk, col, _pct in opts:
        full = np.vstack([blk, v3])
        under = np.median(full[rs:rs + ROAD_CELLS, :], axis=1)
        names.append(tag)
        unders.append(under)
    xi = np.arange(len(opts))
    u = np.array(unders)
    axb.bar(xi - 0.19, u[:, 0], 0.38, color=[o[3] for o in opts],
            label="seaward road cell")
    axb.bar(xi + 0.19, u[:, 1], 0.38, color=[o[3] for o in opts], alpha=0.5,
            label="landward road cell")
    axb.axhline(float(np.median(plat)), color=C["REF"], ls="--", lw=1.3,
                label="backdune platform {:.2f} m".format(float(np.median(plat))))
    axb.set_xticks(xi)
    axb.set_xticklabels(names)
    axb.set_ylabel("elevation under NC-12 (m MHW)")
    axb.set_xlabel("fill option")
    axb.set_title(panel_title("b", "What the roadway sits on"))
    axb.legend(fontsize=7.5)
    axb.grid(alpha=.3, axis="y")

    # ---- (c) how much of each block is invented -------------------------
    axc = fig.add_subplot(gs[1, 1])
    frac = [o[4] for o in opts]
    axc.barh(np.arange(len(opts)), frac, color=[o[3] for o in opts])
    axc.set_yticks(np.arange(len(opts)))
    axc.set_yticklabels(names)
    axc.invert_yaxis()
    axc.set_xlim(0, 100)
    axc.set_xlabel("% of the block taken from the DEM")
    axc.set_title(panel_title("c", "How much is measured"))
    axc.grid(alpha=.3, axis="x")

    fig_h = fig.get_figheight()
    fig.suptitle("What should the added interior rows be made of?",
                 fontsize=11.5, fontweight="bold", x=0.055, ha="left",
                 y=1 - 0.26 / fig_h)
    caption(fig,
            "No survey covers ground that had eroded away by 1996, so every "
            "option is an assertion about 1984, not a measurement of it. The "
            "cells the DEM does\n"
            "hold at those coordinates are the 1996 surface — a later landform "
            "in the same place — which is why option (d), the shipped rule, puts "
            "NC-12 on a\n"
            "3.2/5.0 m dune FACE rather than on backdune. (b) and (c) put the "
            "road where a road behind a dune belongs.\n"
            "(e) keeps every dry cell as measured and fills only the sub-MHW ones. "
            "(f) is a CONTROL, not a candidate: it imports sub-MHW beach, and "
            "is drawn to show what (d) rejects.",
            y=0.09 / fig_h)
    # The caption runs to four lines now that F carries its own control note,
    # so the axes need more room under them than the five-option version did.
    fig.subplots_adjust(top=1 - 0.62 / fig_h, bottom=1.30 / fig_h,
                        left=0.075, right=0.98)

    out = Path(args.out) if args.out else (
        insert_figures_dir("1984-start") / "HAT_fill_options_GIS{}.png".format(D))
    fig.savefig(out)
    print("wrote {}".format(out))
    print("\n  option   road sits on        mean added elev   % from DEM")
    for (tag, lab, blk, col, _p), un, fr in zip(opts, unders, frac):
        print("  {:<7}  {:>5.2f} / {:>5.2f} m      {:>6.2f} m        {:>5.0f}%"
              .format(tag, un[0], un[1], float(np.mean(blk)), fr))


if __name__ == "__main__":
    main()
