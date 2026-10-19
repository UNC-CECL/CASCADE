#!/usr/bin/env python3
"""
One year of BRIE at a 100 m headland: what the two sections at each edge book for the same handover.

    python one_year_figures.py   ->  figures/brie_sand_headland_one_year.png
                                     figures/brie_sand_one_boundary_books.png

Each section's rate comes from the boundary on its north side. At the headland's south edge section 14 books
3.7 m arriving while section 15 (rate 0) books nothing leaving; at the north edge 22.8 m is booked leaving and
nothing arriving. The plain version, then the same with the rule, BRIE's code and one rate per boundary beside it.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-15
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import *  # noqa: E402,F403  (BRIE runs, drawing helpers, constants; puts scripts/ on the path)
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
FIGS = HERE / "figures"

# -----------------------------------------------------------------------------


# 2. A headland, one year, BRIE's real numbers
def fig_headland():
    b, r, theta, s, edges = headland_numbers()
    real1, *_ = run(HEADLAND, years=1)
    fig, ax = plt.subplots(figsize=(13.5, 7.2))
    blank(ax)
    coast(ax, s, base=-160)
    ax.set_ylim(-160, 330)
    ax.text(LO - 0.4, 300, "SEA", fontsize=13, color="#5b86a8", fontweight="bold", va="top")
    ax.text(HI - 0.6, -150, "north →", fontsize=12, color=INK_MUTED, ha="right", va="bottom")
    ax.text(16, 125, "a 100 m headland\n(three 500 m sections)", ha="center", fontsize=11, color=INK_MUTED)
    # each section's rate: shaded strip and label under the shoreline
    for i in range(LO, HI):
        fast = r[i] > 1e-9
        ax.add_patch(Rectangle((i - 0.42, -120), 0.84, 32, fc=GAIN if fast else "white", ec=INK_MUTED, lw=0.8,
                               alpha=0.9 if fast else 1))
        ax.text(i, -104, f"{1000 * r[i]:.0f}" if fast else "0", ha="center", va="center", fontsize=9.5,
                color="white" if fast else INK_MUTED, fontweight="bold")
        if fast:
            arrow(ax, i - 0.05, -78, i + 0.45, s[i + 1] * 0.5 if s[i + 1] > s[i] else s[i] * 0.5, col=GAIN, lw=1.6,
                  ms=12)
    ax.text(LO - 0.45, -135, "each section's rate (per 1,000), set by the boundary on its NORTH side:\n"
            "gold = that boundary is steep, so sand moves; 0 = it faces straight into these waves, so BRIE's rate is 0",
            fontsize=9.5, color=INK, va="top")
    # the two edges
    (s14, _), (s17, _) = edges[14], edges[17]
    bubble(ax, 13.2, 225, f"section 14 (rate {1000 * r[14]:.0f}):\n\"{s14:.1f} m of sand\narrives from the headland\"",
           fs=10)
    bubble(ax, 15.4, 225, "section 15 (rate 0):\n\"nothing leaves me\"", fs=10)
    bubble(ax, 17.0, 225, f"section 17 (rate {1000 * r[17]:.0f}):\n\"{-s17:.1f} m of sand\nleaves me, north\"",
           fs=10)
    bubble(ax, 19.1, 225, "section 18 (rate 0):\n\"nothing arrives\"", fs=10)
    for x, y in ((13.6, 180), (15.2, 190)):
        ax.annotate("", xy=(14.5, 110), xytext=(x, y), arrowprops=dict(arrowstyle="-", color=INK_MUTED, lw=0.8))
    for x, y in ((17.2, 180), (18.8, 190)):
        ax.annotate("", xy=(17.5, 110), xytext=(x, y), arrowprops=dict(arrowstyle="-", color=INK_MUTED, lw=0.8))
    ax.text(14.5, 65, f"+{s14:.1f} m\nappears", ha="center", fontsize=13, fontweight="bold", color=GAIN,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=GAIN, lw=1.4))
    ax.text(17.5, 65, f"{s17:.1f} m\nvanishes".replace("-", "−"), ha="center", fontsize=13, fontweight="bold",
            color=BAD, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=BAD, lw=1.4))
    vol = abs(real1[-1]) * DY * depth(b)
    fig.text(0.5, 0.015, f"After the whole year's step (BRIE's own update): the island has {abs(real1[-1]):.0f} m of "
             f"shoreline less, about {vol / 1000:,.0f},000 m³ (≈ {round(vol / TRUCK_M3, -3):,.0f} truckloads), "
             "although nothing removed any sand.", ha="center", fontsize=12, fontweight="bold", color=BAD)
    fig.suptitle("Why: at a headland, the two sections at each edge disagree about the same handover",
                 fontsize=15, fontweight="bold", y=0.985)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.08)
    png = FIGS / "brie_sand_headland_one_year.png"
    save(fig, png, dpi=200, close=True)
    record_caption(png, (
        "**Tests:** nothing new; one year of an unmodified `brie.Brie` (alongshore step only, option A waves, "
        "500 m sections) at a 100 m headland three sections wide, shown in plan view with north to the right "
        "(a sketch, not to scale). **How to read:** under the shoreline, each section's rate, which BRIE sets "
        "from the angle of the boundary on that section's north side (brie.py 820-826, 1294-1309): only the two "
        "sections whose north boundary is a headland edge get one; flat boundaries face straight into these "
        "waves, where BRIE's rate is zero. The bubbles give what each section books for the handover across the "
        "headland's two edges in one year (2 r Δ, BRIE's own numbers). **Shows:** at the south edge section 14 "
        "books sand arriving but the headland's first section, rate 0, books nothing leaving, so it appears; "
        "at the north edge section 17 books sand leaving and section 18 books nothing arriving, so it vanishes. "
        f"Over BRIE's full year the island loses {abs(real1[-1]):.1f} m of shoreline (× 500 m × "
        f"{depth(b):.1f} m active depth = about {vol:,.0f} m³)."))


# Plan view of the step's cells, with each edge's two books
def books_panel(ax, s, r, shared, title, i):
    lo, hi = 12, 21
    cells = np.arange(lo, hi)
    xs = np.repeat(np.r_[cells - 0.5, cells[-1] + 0.5], 2)[1:-1]
    ys = np.repeat(s[cells], 2)
    ax.fill_between(xs, -60, ys, color=C["ADDED_FILL"], lw=0)
    ax.fill_between(xs, ys, 230, color=C["WATER"], alpha=0.45, lw=0)
    ax.plot(xs, ys, color=INK, lw=1.6)
    made_total = 0.0
    for e in (14, 17):                                 # the step's south and north edges: boundary e|e+1
        d = s[e + 1] - s[e]
        rs_, rn_ = (r[e], r[e]) if shared else (r[e], r[e + 1])
        south, north = 2 * rs_ * d + 0.0, -2 * rn_ * d + 0.0   # what each side books across this boundary in a year
        made = south + north
        made_total += made
        top = max(s[e], s[e + 1])
        ax.plot([e + 0.5, e + 0.5], [-55, 225], color=INK_MUTED, lw=0.8, ls=(0, (2, 2)))
        ax.text(e + 0.45, top + 18, f"south books\n{south:+.2f} m", ha="right", va="bottom", fontsize=6.4, color=INK)
        ax.text(e + 0.55, top + 18, f"north books\n{north:+.2f} m", ha="left", va="bottom", fontsize=6.4, color=INK)
        ok = abs(made) < 0.005
        ax.text(e + 0.5, top + 92, "balanced" if ok else f"{made:+.2f} m {'made' if made > 0 else 'lost'}",
                ha="center", va="bottom", fontsize=7, fontweight="bold", color=SHARED_C if ok else C["GROIN"],
                bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.9))
    ax.text(lo - 0.3, -45, f"whole island, year 1 (first order): {made_total:+.2f} m".replace("-", "−"), fontsize=7,
            fontweight="bold",
            color=SHARED_C if abs(made_total) < 0.005 else C["GROIN"], va="bottom")
    ax.set_xlim(lo - 0.5, hi - 0.5)
    ax.set_ylim(-60, 230)
    ax.set_xticks(cells, [str(c) for c in cells], fontsize=6.5)
    ax.set_xlabel("domain (500 m each), south → north")
    open_frame(ax)
    _letter_inside(ax, i)
    ax.set_title(title, loc="center", fontsize=8.5)


def fig_one_boundary():
    b = brie_on(SHAPES["sharp step"], years=1)
    x = b._x_s.copy()
    r, theta = rates(b, x)
    s = seaward(b, x)
    fig = plt.figure(figsize=figsize("double", height=5.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=(0.42, 1))
    ax = fig.add_subplot(gs[0, :])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3)
    ax.axis("off")
    box(ax, 0.05, 0.2, 3.0, 2.6, "1  What should happen\nsand leaving one domain across a\nboundary arrives in the "
        "next: both\nbooks agree, the total stays the same", "#eef3f8", fs=6.8)
    box(ax, 3.45, 0.2, 3.1, 2.6, "2  What BRIE does (brie.py 820-826, 1294-1324)\neach domain takes ONE rate, r, from "
        "the angle\nof its NORTH boundary, and uses it on BOTH:\nx[i] += r[i] (x[i+1] − 2 x[i] + x[i−1])",
        "#fbeed3", ec=C["ADDED"], fs=6.8)
    box(ax, 6.95, 0.2, 3.0, 2.6, "3  So at a bend\nthe two domains at a boundary use\ndifferent rates for it: "
        "their books\ndisagree, and the gap is sand made or lost", "#fde0dd", ec=C["GROIN"], fs=6.8)
    _letter_inside(ax, 0)
    books_panel(fig.add_subplot(gs[1, 0]), s, r, False, "BRIE as written: a 100 m step, year 1", 1)
    ax = fig.add_subplot(gs[1, 1])
    books_panel(ax, s, r, True, "The same step, one rate per boundary", 2)
    fig.axes[1].set_ylabel("Shoreline (m), + seaward")
    fig.suptitle("Why BRIE's alongshore step can make or lose sand at a bend", fontsize=10)
    png = FIGS / "brie_sand_one_boundary_books.png"
    save(fig, png, dpi=250, close=True)
    k = np.argmax(r)
    record_caption(png, (
        "**Tests:** nothing new; this figure explains the bookkeeping in BRIE's alongshore step at one boundary. "
        "**How to read:** (a) the rule. Sand crossing a boundary should leave one domain and arrive in the next in "
        "the same amount. BRIE computes each domain's angle from the boundary to its north neighbour (brie.py "
        "820-826), turns it into one rate r from its wave-climate diffusivity table (1294-1309), and applies that "
        "rate to both of the domain's boundaries (1311-1324). (b) a real Brie (option A waves) started from a "
        "100 m seaward step three domains wide. At each edge of the step, the domain to the south books "
        "2 r_south Δ across the boundary and the domain to the north books −2 r_north Δ, where Δ is the step "
        "across the boundary and 2 r Δ is the first-order amount in one year of BRIE's Crank-Nicolson step. The "
        "domain whose north boundary is the step gets a rate from the step's angle "
        f"(up to {1000 * r[k]:.1f} per 1,000), the domain beyond it gets the rate of a flat boundary, which under "
        "option A is zero (BRIE clips negative diffusivity to zero), so the books disagree. (c) the same step "
        "with one rate per boundary: both sides use it and the books agree. **Shows:** BRIE's sand budget "
        "balances wherever neighbouring boundaries have the same angle; at a bend it does not, because the "
        "two domains disagree about how much crossed the boundary between them."))

def main():
    apply_style()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig_headland()
    fig_one_boundary()


if __name__ == "__main__":
    main()
