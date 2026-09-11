#!/usr/bin/env python3
r"""
HAT_plot_where_inserts_occur.py
==============================================================================
Where the 1984 footprint changes the domain and where it does not -- zoomed on
the two relocation blocks. Both signs. TWO figures since 2026-09-07 (they
were panels of one):

    HAT_where_inserts_occur_blocks.png   (a) (b) the two relocation blocks,
        GIS 9-14 and 84-87: the measured paired shift per domain with its
        p10-p90, and the rows it becomes under the 10 m rule, +N added / -N
        removed / no change.
    HAT_where_inserts_occur_setback.png  the NC-12 setback at those domains,
        as the model receives it now and from the new row 0 (with p10-p90).
    (an island-wide panel was drawn here too, and RETIRED the same day: it
     showed the same quantity as HAT_footprint_1984_shift.png, which is the
     canonical island-wide view - Hannah's call, 2026-09-07)

REWRITTEN 2026-09-07. This used to contrast two layers (block scope v3 against
island scope v4, both add-only) through their audit CSVs. Those layers were
deleted that day and the footprint became symmetric, so the figures read ONE
table, `2-domain-reconstruction-1984/2-extent/footprint_1984_by_domain.csv`, written by
HAT_footprint_1984.py.

THREE REASONS A DOMAIN IS UNCHANGED - now only one
    Under the old block scope "unchanged" meant either "measured, no cell
    needed" or "never asked". With island-wide scope and a symmetric rule there
    is one reason left: |shift| is under a full 10 m cell.

THE INTERVAL CAVEAT, DRAWN
    The dune lines are 1984 and 1997 -- 13 years. The DEM surface at row 0 is
    1996 ALACE, so the interval wanted is 12 years. The green triangles in the
    blocks figure show N if the measurement were scaled 12/13; recorded, not
    corrected -- scaling would assume steady change across 13 storm years.

USAGE
    python HAT_plot_where_inserts_occur.py
==============================================================================
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "0-elevation" / "3-figures"))
from hat_topo_version import insert_figures_dir, insert_scope_step# noqa: E402
import HAT_plot_duneline_offset as off  # noqa: E402  the house style
from hat_figure_style import (C, DOMAIN_AXIS_LABEL, INK_MUTED, caption,   # noqa: E402
                              figsize, open_frame, save)

FOOTPRINT_CSV = insert_scope_step("1984-start", "2-extent") / "footprint_1984_by_domain.csv"
FIG_DIR = insert_figures_dir("1984-start", "2-extent", "island")     # the blocks: placement-independent
FIG_SEAWARD = insert_figures_dir("1984-start", "3-placement", "seaward")  # the new setback assumes it
BLOCKS = ((list(range(9, 15)), "inter-village", "GIS 9–14, inter-village", 1999),
          (list(range(84, 88)), "Pea Island", "GIS 84–87, Pea Island", 1989))
YEARS_MEASURED, YEARS_WANTED = 13.0, 12.0
CELL_M = 10.0
C_ADD, C_ADD_FILL = off.C_1984, off.C_1984_FILL
C_REM, C_REM_FILL = off.C_1997, off.C_1997_FILL
C_NONE, C_REF, INK = C["BASE_FILL"], C["REF"], off.INK
C_NOW, C_SAME = C["BASE_FILL"], C["BASE"]   # the setback as extracted; unchanged domains


def n_of(m: float) -> int:
    return int(np.trunc(m / CELL_M))


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def load():
    if not FOOTPRINT_CSV.is_file():
        raise SystemExit(f"\n{FOOTPRINT_CSV} not found - run HAT_footprint_1984.py first\n")
    T = {int(r["domain"]): r for r in csv.DictReader(open(FOOTPRINT_CSV, encoding="utf-8"))}
    D = {
        "shift": {d: _f(r["shift_m_median"]) for d, r in T.items()},
        "p10": {d: _f(r["shift_m_p10"]) for d, r in T.items()},
        "p90": {d: _f(r["shift_m_p90"]) for d, r in T.items()},
        "n": {d: int(r["n_cells"]) for d, r in T.items()},
        "sb_now": {d: _f(r["setback_model_now_m"]) for d, r in T.items()},
        "sb_new": {d: _f(r["setback_new_m"]) for d, r in T.items()},
        "sb_p10": {d: _f(r["setback_new_p10_m"]) for d, r in T.items()},
        "sb_p90": {d: _f(r["setback_new_p90_m"]) for d, r in T.items()},
    }
    D["colour"] = {d: (C_ADD if v > 0 else C_REM if v < 0 else C_NONE) for d, v in D["n"].items()}
    return D


# =============================================================================

def fig_blocks(D) -> Path:
    off.apply_style()
    shift, p10, p90, n, colour = D["shift"], D["p10"], D["p90"], D["n"], D["colour"]
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", aspect=0.48),
                             constrained_layout=True,
                             gridspec_kw=dict(width_ratios=[6, 4]))
    for j, (ax, (blk, short, name, yr)) in enumerate(zip(axes, BLOCKS)):
        d = [g for g in blk if g in shift]
        xi = np.arange(len(d))
        sm = np.array([shift[g] for g in d])
        ax.bar(xi, [n[g] * CELL_M for g in d], 0.62, color=[colour[g] for g in d], zorder=3)
        ax.errorbar(xi, sm, yerr=[sm - np.array([p10[g] for g in d]),
                                  np.array([p90[g] for g in d]) - sm],
                    fmt="o", ms=3.0, color=INK, ecolor=INK_MUTED, elinewidth=0.8, zorder=5)
        for k, g in enumerate(d):
            lab = f"{n[g]:+d} row{'s' if abs(n[g]) != 1 else ''}" if n[g] else "no change"
            top = max(sm[k], p90[g], n[g] * CELL_M)
            ax.text(k, top + 2.5, lab, ha="center", va="bottom", fontsize=7.5,
                    color=colour[g] if n[g] else INK_MUTED,
                    fontweight="bold" if n[g] else "normal")
            n12 = n_of(shift[g] * YEARS_WANTED / YEARS_MEASURED)
            if n12 != n[g]:
                y12 = shift[g] * YEARS_WANTED / YEARS_MEASURED
                ax.plot([k], [y12], "v", color=C_REF, ms=6, zorder=6)
                ax.text(k + 0.36, y12, f"→ {n12:+d}", ha="left", va="center",
                        fontsize=7.5, color=C_REF)
        for t in (-CELL_M, CELL_M):
            ax.axhline(t, color=INK_MUTED, ls="--", lw=0.7, zorder=2)
        ax.axhline(0, color=INK, lw=0.7, zorder=2)
        ax.set_xlim(-0.7, len(d) - 0.3)
        ax.set_xticks(xi)
        ax.set_xticklabels(d)
        ax.set_xlabel(DOMAIN_AXIS_LABEL)
        if j == 0:
            ax.set_ylabel("1984 line − 1997 line (m)\n+ seaward")
        lo = min(-15.0, float(np.nanmin([p10[g] for g in d])) - 8)
        hi = max(30.0, float(np.nanmax([p90[g] for g in d])) + 18)
        ax.set_ylim(lo, hi)
        ax.grid(axis="y")
        ax.set_axisbelow(True)
        open_frame(ax)
        off._title(ax, j, f"{short}, relocated {yr}")
    fig.legend(handles=[Line2D([0], [0], marker="o", ms=3.0, color=INK, linestyle="none",
                               label="median of 50 paired profiles, p10–p90"),
                        Patch(facecolor=C_ADD, label="rows added (N × 10 m)"),
                        Patch(facecolor=C_REM, label="rows removed (N × 10 m)"),
                        Patch(facecolor=C_NONE, label="under one cell, unchanged"),
                        Line2D([0], [0], color=INK_MUTED, linestyle="--", label="±1 cell"),
                        Line2D([0], [0], marker="v", color=C_REF, linestyle="none", ms=6,
                               label="N if scaled to the 12-year interval")],
               loc="outside lower center", ncol=3, frameon=False)
    caption(fig,
            "Where the 1984 footprint changes the domain, in the two relocation "
            "blocks: (a) GIS 9–14, the inter-village reach relocated in 1999, "
            "and (b) GIS 84–87 on Pea Island, relocated in 1989. Points are the "
            "median paired shift of the 1984 dune line seaward of the 1997 line "
            "over 50 profiles, with the p10–p90 range; bars are the rows that "
            "shift becomes under the 10 m rule (truncated to whole cells): red "
            "rows added landward of NC-12, blue rows removed, grey under one "
            "cell and unchanged. Dashed lines mark one cell either way. The "
            "dune lines are 1984 and 1997 (13 years) while the surface at row 0 "
            "is 1996 ALACE (12 years wanted); a green triangle shows N if the "
            "shift were scaled by 12/13, recorded, not applied, because scaling "
            "would assume steady change across 13 storm years.")
    p = FIG_DIR / "HAT_where_inserts_occur_blocks.png"
    save(fig, p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_setback(D) -> Path:
    off.apply_style()
    n, colour = D["n"], D["colour"]
    sb_now, sb_new, sp10, sp90 = D["sb_now"], D["sb_new"], D["sb_p10"], D["sb_p90"]
    d = [g for blk, *_ in BLOCKS for g in blk if g in sb_new and np.isfinite(sb_new[g])]
    xi = np.arange(len(d))
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.46), constrained_layout=True)
    # a gap between the two blocks
    gap = next(i for i, g in enumerate(d) if g > 20)
    xi = np.where(np.arange(len(d)) >= gap, xi + 0.8, xi).astype(float)
    ax.bar(xi - 0.2, [sb_now[g] for g in d], 0.4, color=C_NOW, zorder=3)
    ax.bar(xi + 0.2, [sb_new[g] for g in d], 0.4, zorder=3,
           color=[colour[g] if n[g] else C_SAME for g in d])
    ax.errorbar(xi + 0.2, [sb_new[g] for g in d],
                yerr=[[max(sb_new[g] - sp10[g], 0) for g in d], [max(sp90[g] - sb_new[g], 0) for g in d]],
                fmt="none", ecolor=INK, elinewidth=0.8, capsize=2, zorder=5)
    for k, g in enumerate(d):
        ax.text(xi[k] - 0.2, sb_now[g] + 1.0, f"{sb_now[g]:.0f}", ha="center", va="bottom",
                fontsize=6.5, color=INK_MUTED)
        ax.text(xi[k] + 0.2, max(sp90[g], sb_new[g]) + 1.0, f"{sb_new[g]:.0f}", ha="center",
                va="bottom", fontsize=6.5, color=colour[g] if n[g] else INK, fontweight="bold")
    ax.axhline(0, color=INK, lw=0.7)
    ax.set_ylim(min(0.0, float(np.nanmin([sp10[g] for g in d])) - 5),
                float(np.nanmax([sp90[g] for g in d])) + 14)
    # the two blocks, named once each, above their bars
    for lo, hi, (_, _s, name, yr) in ((0, gap - 1, BLOCKS[0]), (gap, len(d) - 1, BLOCKS[1])):
        ax.text((xi[lo] + xi[hi]) / 2, 0.985, f"{name}, relocated {yr}",
                transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=7.5, color=INK_MUTED)
    ax.set_xticks(xi)
    ax.set_xticklabels(d)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("NC-12 setback (m landward of row 0)")
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    open_frame(ax)
    fig.legend(handles=[Patch(facecolor=C_NOW, label="setback as extracted (1996 surface), floored at 0"),
                        Patch(facecolor=C_ADD, label="1984 reconstruction, rows added"),
                        Patch(facecolor=C_REM, label="1984 reconstruction, rows removed"),
                        Patch(facecolor=C_SAME, label="1984 reconstruction, unchanged"),
                        Line2D([0], [0], color=INK, lw=0.8, label="p10–p90 over the profiles")],
               loc="outside lower center", ncol=3, frameon=False)
    caption(fig,
            "The NC-12 setback at the domains of the two relocation blocks "
            "(GIS 9–14, the inter-village reach relocated in 1999; GIS 84–87 "
            "on Pea Island, relocated in 1989): as extracted from the 1996 "
            "surface and floored at 0 (light grey), and from the new interior "
            "row 0 of the 1984 reconstruction, coloured by what the footprint "
            "did to that domain (red rows added, blue rows removed, dark grey "
            "unchanged), with the p10–p90 range over the 50 profiles. Values "
            "are labelled above each bar in metres.")
    p = FIG_SEAWARD / "HAT_where_inserts_occur_setback.png"
    save(fig, p, bbox_inches="tight")
    plt.close(fig)
    return p


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.parse_args()
    D = load()
    for f in (fig_blocks(D), fig_setback(D)):
        print(f"wrote {f}")
    nn = np.array(list(D["n"].values()))
    print(f"  added {int((nn > 0).sum())} domains / {int(nn[nn > 0].sum())} rows; "
          f"removed {int((nn < 0).sum())} / {int(-nn[nn < 0].sum())}; unchanged {int((nn == 0).sum())}")
    for stale in (FIG_DIR / "HAT_where_inserts_occur.png",
                  FIG_DIR / "HAT_where_inserts_occur_island.png"):
        if stale.is_file():
            stale.unlink()
            print(f"  removed {stale.name} (retired 2026-09-07)")


if __name__ == "__main__":
    main()
