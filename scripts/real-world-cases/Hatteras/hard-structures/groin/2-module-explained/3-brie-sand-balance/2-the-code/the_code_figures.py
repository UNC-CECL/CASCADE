#!/usr/bin/env python3
"""
BRIE's alongshore step in its own code: "it looks at both sides" (true of the positions, not of the rate), and the
three lines worked through by hand on six sections.

    python the_code_figures.py   ->  figures/brie_sand_both_sides_one_rate.png
                                     figures/brie_sand_equation_worked_example.png
                                     figures/brie_sand_equation_step_by_step.gif   (git-ignored)

The update x[i] += r[i] (x[i+1] - 2 x[i] + x[i-1]) takes the shoreline positions from both sides, but r[i] comes
from one angle only: brie.py 824 builds it from x_s[np.r_[1:ny, 0]] - x_s, which is x[i+1] - x[i], the north
neighbour. The worked example is a real Brie on six sections, x_s = 1000, 1000, 900, 900, 900, 1000 m (landward +);
every number is BRIE's own (angle 820-828, rate 1294-1309, explicit half of the update 1315-1324, the year
from Brie.update()). The code shown is condensed (`self._` dropped, the index clamping omitted).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-15
"""
from __future__ import annotations

# Non-code files beside this script live in the workspace on D: (hatteras_paths.workspace_twin)
import sys as _hat_sys  # noqa: E402
from pathlib import Path as _HatPath  # noqa: E402
_HAT_DIR = next(_p for _p in _HatPath(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import workspace_twin  # noqa: E402
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import *  # noqa: E402,F403  (BRIE runs, drawing helpers, constants; puts scripts/ on the path)
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
FIGS = workspace_twin(HERE) / "figures"
POS_C, RATE_C, ANG_C = C["LATE"], C["ADDED"], C["ACCENT"]   # positions, rate, angle
MONO = "DejaVu Sans Mono"
X0 = np.array([1000.0, 1000.0, 900.0, 900.0, 900.0, 1000.0])    # BRIE's x_s (m, landward +)
N = len(X0)
NEXT, PREV = np.r_[1:N, 0], np.r_[N - 1, 0:N - 1]
CODE = [
    ("820-828", "theta = 180 * np.arctan2(x_s[np.r_[1:ny, 0]] - x_s, dy) / np.pi", "angle"),
    ("1294-1309", "r = np.maximum(0, coast_diff[round(90 - theta)] * dt / 2 / dy**2)", "rate"),
    ("1315-1324", "RHS = x_s + r * (x_s[np.r_[1:ny, 0]] - 2*x_s + x_s[np.r_[ny-1, 0:ny-1]])", "update"),
]
# -----------------------------------------------------------------------------


def chip(ax, x, y, text, col, fs=17, ha="left"):
    return ax.text(x, y, text, fontsize=fs, ha=ha, va="center", family="DejaVu Sans Mono", color=INK,
                   bbox=dict(boxstyle="round,pad=0.25", fc=col if col else "none", ec="none", alpha=0.25 if col else 0))


# 5. The formula, and one boundary with one rate everywhere and with BRIE's rates
def fig_both_sides():
    b, r, theta, s, edges = headland_numbers()
    fig = plt.figure(figsize=(13.5, 9.0))
    gs = fig.add_gridspec(3, 2, height_ratios=(0.62, 1, 0.36), hspace=0.32, wspace=0.08, left=0.03, right=0.97,
                          top=0.9, bottom=0.03)
    ax = fig.add_subplot(gs[0, :])
    blank(ax)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 10)
    ax.text(1, 9.3, "BRIE's update for one section, every year (brie.py 1311-1324):", fontsize=13, va="top")
    parts = [("x[i]  +=  ", None), ("r[i]", RATE_C), ("  × (  ", None), ("x[i+1]", POS_C), ("  − 2·x[i]  +  ", None),
             ("x[i−1]", POS_C), ("  )", None)]
    widths = [12.8, 5.0, 8.2, 7.6, 19.0, 7.6, 3]
    x = 6
    centres = {}
    for (t, col), w in zip(parts, widths):
        chip(ax, x, 5.8, t, col)
        centres[t] = x + w / 2
        x += w
    ax.annotate("", xy=(centres["x[i+1]"], 4.7), xytext=(68, 2.2), arrowprops=dict(arrowstyle="-|>", color=POS_C, lw=1.6))
    ax.annotate("", xy=(centres["x[i−1]"], 4.7), xytext=(72, 2.2), arrowprops=dict(arrowstyle="-|>", color=POS_C, lw=1.6))
    ax.text(70, 1.4, "positions from BOTH sides (north and south): true\nthis is what \"it looks at both sides\" means",
            ha="center", va="center", fontsize=12, color=POS_C, fontweight="bold")
    ax.annotate("", xy=(centres["r[i]"], 4.7), xytext=(19, 2.2), arrowprops=dict(arrowstyle="-|>", color=RATE_C, lw=1.6))
    ax.text(19, 1.4, "ONE rate, from the angle of the NORTH boundary only\n(brie.py 820-826), used on both sides",
            ha="center", va="center", fontsize=12, color=RATE_C, fontweight="bold")

    def boundary(ax, rs, rn, title, verdict, col):
        blank(ax)
        lo, hi = 12, 17
        coast(ax, s, lo=lo, hi=hi, base=-150)
        ax.set_ylim(-150, 330)
        ax.set_title(title, fontsize=13, pad=8)
        d = s[15] - s[14]
        south, north = 2 * rs * d, -2 * rn * d + 0.0
        for i, rr in ((14, rs), (15, rn)):
            ax.add_patch(Rectangle((i - 0.4, -125), 0.8, 36, fc=RATE_C if rr > 1e-9 else "white", ec=INK_MUTED, lw=0.8))
            ax.text(i, -107, f"rate {1000 * rr:.0f}", ha="center", va="center", fontsize=10.5,
                    color="white" if rr > 1e-9 else INK_MUTED, fontweight="bold")
        ax.plot([14.5, 14.5], [-150, 300], color=INK, lw=1.4, ls=(0, (4, 3)))
        ax.text(14.5, 310, "the boundary 14 | 15", ha="center", fontsize=10.5, color=INK_MUTED)
        ax.text(14.5, 160, "sand flows south ←", ha="center", fontsize=10, color=INK_MUTED)
        for i, val, lab, x0 in ((14, south, "14 writes", 13.4), (15, -north, "15 writes", 15.6)):
            if abs(val) > 1e-6:                    # 15's arrow leaves it toward the boundary; 14's arrives from it
                a0, a1 = ((15.45, 200), (14.6, 200)) if i == 15 else ((14.4, 200), (13.55, 200))
                ax.add_patch(FancyArrowPatch(a0, a1, arrowstyle="-|>", mutation_scale=22, lw=4.5, color=col))
            ax.text(x0, 255, f"{lab}:\n{abs(val):.1f} m {'arrives' if i == 14 else 'leaves'}" if abs(val) > 1e-6
                    else f"{lab}:\nnothing leaves", ha="center", fontsize=11,
                    bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=col, lw=1.2))
        ax.text(14.5, 60, verdict, ha="center", fontsize=14, fontweight="bold", color=col,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=col, lw=1.4))
        ax.text(lo - 0.4, -146, "north →", fontsize=10, color=INK_MUTED, va="bottom")

    boundary(fig.add_subplot(gs[1, 0]), r[14], r[14], "If both sections used the same rate",
             "books agree: 0 m", GOOD)
    boundary(fig.add_subplot(gs[1, 1]), r[14], r[15], "What BRIE does: each uses its own north boundary's rate",
             f"+{edges[14][0]:.1f} m appears", BAD)
    ax = fig.add_subplot(gs[2, :])
    blank(ax)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 1)
    ax.text(5, 0.62, "Add up every section's change and the positions cancel, but each boundary leaves behind\n"
            "(rate of the section to its south − rate of the section to its north) × (the step across the boundary)",
            ha="center", va="center", fontsize=13, linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.5", fc="#f7f7f7", ec=INK_MUTED, lw=0.8))
    ax.text(5, 0.08, "Neighbours' rates match (a straight coast, or one curving smoothly): about 0.   "
            "Rates jump between neighbours (a bend or a kink): sand made or lost.", ha="center", fontsize=12,
            color=INK)
    fig.suptitle("\"BRIE looks at both sides\" — for the positions, yes; for the rate, no", fontsize=17,
                 fontweight="bold")
    png = workspace_twin(FIGS) / "brie_sand_both_sides_one_rate.png"
    save(fig, png, dpi=200, close=True)
    record_caption(png, (
        "**Tests:** nothing new; this figure answers the objection that BRIE cannot make or lose sand because "
        "its step looks at both sides of each section. **How to read:** top, BRIE's update for one section "
        "(brie.py 1311-1324): the shoreline positions come from both neighbours (blue), the rate r from the angle "
        "of the section's north boundary only (gold; brie.py 820-826, 1294-1309) and is applied to both sides. "
        "Middle, the boundary at the south edge of a 100 m headland in one year of an unmodified `brie.Brie` "
        "(alongshore step only, option A waves): with one rate for both sections their books agree; with "
        "BRIE's own rates, section 14 (whose north boundary is the headland's edge) books 3.7 m arriving and "
        "section 15 (whose north boundary is flat, rate 0 under these waves) books nothing leaving. Bottom, the "
        "sum over the island. **Shows:** the objection is right about the positions and, for a single rate, "
        "about conservation; the error comes from each section having its own rate, so it appears only where "
        "neighbouring rates differ, at a bend."))


def brie_numbers():
    b = Brie(name="walk", barrier_model=False, inlet_model=False, ast_model=True, **WAVES,
             alongshore_section_length=500.0, alongshore_section_count=N, time_step=1.0, time_step_count=3,
             sea_level_rise_rate=0.0)
    b._x_s = X0.copy()
    x = X0
    theta = 180 * np.arctan2(x[NEXT] - x, b._dy) / np.pi                       # 820-828
    idx = np.maximum(1, np.minimum(b._wave_climl, np.round(90 - theta).astype(int)))
    D = b._coast_diff[idx]
    r = np.maximum(0, D * b._dt / 2 / b._dy ** 2)                                # 1294-1309
    north = r * (x[NEXT] - x)                                                    # 1315-1324, split in two
    south = r * (x[PREV] - x)
    ang = np.arange(-30, 31)
    curve = b._coast_diff[np.maximum(1, np.minimum(b._wave_climl, np.round(90 - ang).astype(int)))]
    b.update()
    full = b._x_s - X0
    return dict(theta=theta, D=D, r=r, north=north, south=south, ang=ang, curve=curve, full=full, dy=b._dy, dt=b._dt)


def fmt(v, d=2):
    s = f"{v:+.{d}f}" if abs(v) >= 10 ** -d / 2 else "0"
    return s.replace("-", "−")


def code_panel(ax, active):
    blank(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.0, 0.97, "BRIE's code (brie.py, condensed)", fontsize=11, color=INK_MUTED, va="top")
    for k, (lines, src, key) in enumerate(CODE):
        y = 0.7 - k * 0.27
        on = key == active
        ax.text(0.0, y, f"{lines:>9}", fontsize=10, family=MONO, color=INK_MUTED, va="center")
        ax.text(0.12, y, src, fontsize=11.2, family=MONO, color=INK if on else "#9a9a9a", va="center",
                fontweight="bold" if on else "normal",
                bbox=dict(boxstyle="round,pad=0.3", fc="#fff4d6" if on else "none", ec=RATE_C if on else "none",
                          lw=1.2))


def coast_panel(ax, hi_sec=(), hi_col=POS_C, arrows=(), cross=None, show_x=True, rate=None, theta=None):
    blank(ax)
    sea = 1000 - X0                                  # seaward height of each section
    for i in range(N):
        ax.add_patch(Rectangle((i, -40), 1, sea[i] + 40, fc=SAND, ec="none"))
        ax.add_patch(Rectangle((i, sea[i]), 1, 240 - sea[i], fc=SEA, ec="none"))
    xs = np.repeat(np.r_[np.arange(N), N], 2)[1:-1]
    ax.plot(xs, np.repeat(sea, 2), color=INK, lw=2.2)
    for i in range(1, N):
        ax.plot([i, i], [-40, max(sea[i - 1], sea[i])], color="#b8913f", lw=0.8, ls=(0, (2, 2)))
    for i in range(N):
        on = i in hi_sec
        ax.add_patch(Rectangle((i + 0.04, -38), 0.92, sea[i] + 36, fc="none", ec=hi_col if on else "none", lw=3))
        ax.text(i + 0.5, -62, f"section {i}", ha="center", fontsize=10, fontweight="bold" if on else "normal",
                color=hi_col if on else INK)
        if show_x:
            ax.text(i + 0.5, -12, f"x = {X0[i]:.0f}", ha="center", fontsize=9.5, family=MONO,
                    color=POS_C if on else INK)
        if theta is not None and not np.isnan(theta[i]):
            ax.text(i + 0.5, 205, f"θ = {fmt(theta[i], 1)}°", ha="center", fontsize=9.5, color=ANG_C)
        if rate is not None and not np.isnan(rate[i]):
            ax.text(i + 0.5, 175, f"r = {rate[i]:.3f}", ha="center", fontsize=9.5, color=RATE_C, fontweight="bold")
    for (a, bb_, col, lab) in arrows:
        ya = 150
        ax.add_patch(FancyArrowPatch((a + 0.5, ya), (bb_ + 0.5, ya), arrowstyle="-|>", mutation_scale=18, lw=2.4,
                                     color=col, connectionstyle="arc3,rad=-0.35"))
        if lab:
            ax.text((a + bb_) / 2 + 0.5, ya + 32, lab, ha="center", fontsize=9.5, color=col)
    if cross is not None:
        ax.text(cross + 0.5, 120, "(not read)", ha="center", fontsize=10, color=INK_MUTED)
    ax.text(N - 0.05, 228, "north →", ha="right", fontsize=10, color=INK_MUTED, va="top")
    ax.text(0.05, 228, "sea", fontsize=10, color="#5b86a8", va="top")
    ax.set_xlim(-0.05, N + 0.05)
    ax.set_ylim(-75, 240)


def note(ax, text, fs=11.5, color=INK):
    blank(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.0, 1.0, text, fontsize=fs, va="top", color=color, linespacing=1.45, wrap=True)


def frame_layout():
    fig = plt.figure(figsize=(12.5, 8.2), dpi=84)
    ax_title = fig.add_axes([0.02, 0.9, 0.96, 0.09])
    ax_code = fig.add_axes([0.02, 0.66, 0.96, 0.23])
    ax_coast = fig.add_axes([0.02, 0.27, 0.56, 0.38])
    ax_side = fig.add_axes([0.62, 0.05, 0.36, 0.58])
    ax_note = fig.add_axes([0.02, 0.02, 0.56, 0.22])
    for a in (ax_title,):
        blank(a)
    return fig, ax_title, ax_code, ax_coast, ax_side, ax_note


def title(ax, step, text):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.0, 0.75, step, fontsize=11, color=INK_MUTED, va="center")
    ax.text(0.0, 0.25, text, fontsize=16, fontweight="bold", va="center")


def table(ax, rows, header, hi_row=None, col_colors=None, fs=10.5):
    blank(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ncol = len(header)
    xs = np.linspace(0.02, 0.98, ncol + 1)[:-1] + 0.48 / ncol
    for j, h in enumerate(header):
        ax.text(xs[j], 0.96, h, ha="center", va="top", fontsize=fs, fontweight="bold",
                color=(col_colors or [INK] * ncol)[j])
    for k, row in enumerate(rows):
        y = 0.86 - k * 0.12
        if k == hi_row:
            ax.add_patch(Rectangle((0.0, y - 0.055), 1.0, 0.11, fc="#fff4d6", ec="none"))
        for j, v in enumerate(row):
            ax.text(xs[j], y, v, ha="center", va="center", fontsize=fs, family=MONO,
                    color=(col_colors or [INK] * ncol)[j])


def gif_equation(nums):
    from PIL import Image
    th, r, north, south = nums["theta"], nums["r"], nums["north"], nums["south"]
    frames, durs = [], []

    def grab(fig, ms):
        fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert("RGB"))
        durs.append(ms)
        plt.close(fig)

    # 0. the coast
    fig, at, ac, ao, asd, an = frame_layout()
    title(at, "BRIE's alongshore step, line by line", "A coast of six 500 m sections with a 100 m headland")
    code_panel(ac, None)
    coast_panel(ao)
    note(an, "x is BRIE's shoreline position for each section (x_s), measured landward:\n"
         "the headland sections (2, 3, 4) sit 100 m further seaward, so their x is 900.\n"
         "BRIE does three things each year, the three lines above. We follow them in order.")
    table(asd, [(str(i), f"{X0[i]:.0f}") for i in range(N)], ["section i", "x[i] (m)"], col_colors=[INK, POS_C])
    grab(fig, 6500)
    # 1. what np.r_[1:ny, 0] is
    fig, at, ac, ao, asd, an = frame_layout()
    title(at, "Line 1 of 3: the angle", "Which neighbour does the angle line read?")
    code_panel(ac, "angle")
    coast_panel(ao)
    note(an, "x_s[np.r_[1:ny, 0]] is the same list of positions shifted by one place:\n"
         "for section i it gives x[i+1], the NEXT section to the north (the last one wraps to 0).\n"
         "So the angle uses x[i+1] − x[i]: section i and its NORTH neighbour only.")
    table(asd, [(str(i), str(NEXT[i]), str(PREV[i])) for i in range(N)],
          ["i", "np.r_[1:ny, 0]\n= i+1 (north)", "i−1 (south)\nnot in this line"],
          col_colors=[INK, POS_C, INK_MUTED])
    grab(fig, 9000)
    # 2. each section's angle, one at a time
    for i in range(N):
        fig, at, ac, ao, asd, an = frame_layout()
        title(at, "Line 1 of 3: the angle", f"Section {i}: its angle comes from x[{NEXT[i]}] − x[{i}]")
        code_panel(ac, "angle")
        shown = np.where(np.arange(N) <= i, th, np.nan)
        coast_panel(ao, hi_sec=(i, NEXT[i]), arrows=[(i, NEXT[i], ANG_C, "reads north")] if NEXT[i] != 0 else [],
                    cross=PREV[i], theta=shown)
        d = X0[NEXT[i]] - X0[i]
        note(an, f"θ[{i}] = atan( (x[{NEXT[i]}] − x[{i}]) / 500 m ) = atan( ({X0[NEXT[i]]:.0f} − {X0[i]:.0f}) / 500 )"
             f" = {fmt(th[i], 1)}°\n"
             + ("Section i−1 = " + str(PREV[i]) + " is not used for this angle."
                if d == 0 else "This boundary is the headland's edge, so the angle is not zero."))
        table(asd, [(str(k), fmt(X0[NEXT[k]] - X0[k], 0) if k <= i else "", fmt(th[k], 1) + "°" if k <= i else "")
                    for k in range(N)], ["i", "x[i+1] − x[i]", "θ[i]"], hi_row=i, col_colors=[INK, POS_C, ANG_C])
        grab(fig, 3000 if d == 0 else 4500)
    # 3. angle -> rate
    fig, at, ac, ao, asd, an = frame_layout()
    title(at, "Line 2 of 3: the rate", "Each section's ONE rate comes from its ONE angle")
    code_panel(ac, "rate")
    coast_panel(ao, theta=th, rate=r)
    note(an, "BRIE looks up the wave climate's diffusivity D at the angle (the curve on the right), then\n"
         f"r = max(0, D × dt / 2 / dy²) = max(0, D × 1 / 2 / 500²).  Facing straight into these waves\n"
         f"(θ = 0) D is {nums['D'][0]:.0f}, negative, so BRIE sets r = 0. At ±11.3° D is positive and r is not 0.\n"
         "Section 2 is part of the headland, but its north side (2 | 3) is flat, so its rate is 0.")
    asd.clear()
    asd.set_position([0.69, 0.08, 0.29, 0.52])
    asd.plot(nums["ang"], nums["curve"], color=INK, lw=1.8)
    asd.set_ylim(-40000, 165000)
    asd.axhline(0, color=INK_MUTED, lw=0.7)
    for k in range(N):
        asd.scatter(th[k], nums["D"][k], s=70, color=RATE_C if r[k] > 0 else INK_MUTED, zorder=5)
    for a_, lab in ((-11.31, f"sections 1: D = {nums['D'][1]:.0f}\n→ r = {r[1]:.3f}"),
                    (11.31, f"section 4: D = {nums['D'][4]:.0f}\n→ r = {r[4]:.3f}"), (0, "sections 0, 2, 3, 5: D < 0 → r = 0")):
        asd.annotate(lab, xy=(a_, np.interp(a_, nums["ang"], nums["curve"])),
                     xytext=(a_ - 12 if a_ > 5 else a_ + 2 if a_ < -5 else -29, 42000 if a_ > 5 else 22000 if a_ < -5 else -33000),
                     fontsize=9, color=RATE_C if a_ else INK_MUTED, arrowprops=dict(arrowstyle="-", color=INK_MUTED, lw=0.6))
    asd.set_xlabel("angle θ of the section's north boundary (°)", fontsize=9.5)
    asd.set_ylabel("BRIE's diffusivity D (m²/yr)", fontsize=9.5)
    asd.set_title("BRIE's table for these waves (coast_diff)", fontsize=10)
    asd.tick_params(labelsize=8.5)
    grab(fig, 12000)
    # 4. the update, section by section
    for i in (1, 2, 4, 5):
        fig, at, ac, ao, asd, an = frame_layout()
        title(at, "Line 3 of 3: the update", f"Section {i}: positions from BOTH sides, times its ONE rate")
        code_panel(ac, "update")
        coast_panel(ao, hi_sec=(PREV[i], i, NEXT[i]), rate=np.where(np.arange(N) == i, r, np.nan),
                    arrows=[(i, NEXT[i], POS_C, ""), (i, PREV[i], POS_C, "")])
        note(an, f"change[{i}] = r[{i}] × ( x[{NEXT[i]}] − 2·x[{i}] + x[{PREV[i]}] )\n"
             f"          = r[{i}] × ( x[{NEXT[i]}] − x[{i}] )  +  r[{i}] × ( x[{PREV[i]}] − x[{i}] )\n"
             f"          =   {r[i]:.3f} × {fmt(X0[NEXT[i]] - X0[i], 0)}          +  {r[i]:.3f} × {fmt(X0[PREV[i]] - X0[i], 0)}\n"
             f"          =   {fmt(north[i])} (with north)  +  {fmt(south[i])} (with south)  =  {fmt(north[i] + south[i])} m",
             fs=11, color=INK)
        an.texts[0].set_family(MONO)
        table(asd, [(str(k), f"{r[k]:.3f}", fmt(north[k]), fmt(south[k])) for k in range(N)],
              ["i", "r[i]", "with north\nr[i]·(x[i+1]−x[i])", "with south\nr[i]·(x[i−1]−x[i])"], hi_row=i,
              col_colors=[INK, RATE_C, POS_C, POS_C], fs=10)
        grab(fig, 7000)
    # 5. the books at each boundary
    fig, at, ac, ao, asd, an = frame_layout()
    title(at, "Putting it together", "The same boundary, booked by two sections with two different rates")
    code_panel(ac, "update")
    coast_panel(ao, hi_sec=(1, 2, 4, 5), hi_col=BAD, rate=r)
    rows = []
    for k in range(N):
        a_, b_ = k, NEXT[k]
        rows.append((f"{a_} | {b_}", fmt(north[a_]), fmt(south[b_]), fmt(north[a_] + south[b_])))
    total = (north + south).sum()
    note(an, "At the boundary between sections 1 and 2:\n"
         f"  section 1 books  r[1] × (x[2] − x[1]) = {r[1]:.3f} × (−100) = {fmt(north[1])} m\n"
         f"  section 2 books  r[2] × (x[1] − x[2]) = {r[2]:.3f} × (+100) = {fmt(south[2])} m\n"
         "They should cancel. They don't, because r[1] ≠ r[2].", fs=11)
    an.texts[0].set_family(MONO)
    table(asd, rows + [("all", "", "", fmt(total))], ["boundary", "south section\nbooks", "north section\nbooks",
                                                      "left over"], col_colors=[INK, POS_C, POS_C, BAD], fs=10)
    grab(fig, 11000)
    # 6. one shared rate
    fig, at, ac, ao, asd, an = frame_layout()
    title(at, "The comparison", "If both sections used the same rate for their shared boundary")
    code_panel(ac, None)
    coast_panel(ao, hi_sec=(1, 2), hi_col=GOOD, rate=r)
    note(an, f"Give boundary 1 | 2 one rate, {r[1]:.3f}, used by both sections:\n"
         f"  section 1 books  {r[1]:.3f} × (−100) = {fmt(north[1])} m\n"
         f"  section 2 books  {r[1]:.3f} × (+100) = {fmt(-north[1])} m\n"
         "They cancel: whatever leaves one section arrives in the other. Every boundary like this → 0.",
         fs=11, color=GOOD)
    an.texts[0].set_family(MONO)
    table(asd, [(f"{k} | {NEXT[k]}", fmt(north[k]), fmt(-north[k]), "0") for k in range(N)] + [("all", "", "", "0")],
          ["boundary", "south section\nbooks", "north section\nbooks", "left over"],
          col_colors=[INK, POS_C, POS_C, GOOD], fs=10)
    grab(fig, 9000)
    # 7. the real year
    fig, at, ac, ao, asd, an = frame_layout()
    title(at, "The whole year, from Brie.update() itself",
          f"The island's shoreline moved {fmt(nums['full'].sum(), 1)} m landward in total")
    code_panel(ac, None)
    coast_panel(ao, rate=r)
    note(an, f"The lines above are the explicit half of BRIE's step (left over: {fmt(total)} m). BRIE also solves an\n"
         "implicit half with the same per-section rates, so the year's total is larger. Nothing added or\n"
         f"removed sand, yet the summed x grew by {nums['full'].sum():.1f} m: {nums['full'].sum():.1f} m × 500 m × "
         "19.8 m ≈ 148,000 m³ of sand gone.\nA copy of these lines reproduces Brie.update() to three decimals.", fs=11)
    table(asd, [(str(k), fmt(nums["full"][k])) for k in range(N)] + [("sum", fmt(nums["full"].sum()))],
          ["section", "change in x this year (m)"], col_colors=[INK, BAD], fs=10.5)
    grab(fig, 12000)
    out = workspace_twin(FIGS) / "brie_sand_equation_step_by_step.gif"
    frames[0].save(out, save_all=True, append_images=frames[1:], duration=durs, loop=0, optimize=True)
    print("wrote", out, f"| {len(frames)} frames, {sum(durs) / 1000:.0f} s")


def fig_worked_example(nums):
    th, r, north, south = nums["theta"], nums["r"], nums["north"], nums["south"]
    fig, axes = plt.subplots(1, 4, figsize=(17, 5.6), gridspec_kw=dict(width_ratios=(1, 1.1, 1.3, 1.2)))
    steps = [("1  The angle (brie.py 820-828)", "θ[i] = atan( (x[i+1] − x[i]) / 500 )",
              [(str(i), f"{X0[i]:.0f}", f"{X0[NEXT[i]]:.0f}", fmt(th[i], 1)) for i in range(N)],
              ["i", "x[i]", "x[i+1]", "θ[i] °"], [INK, POS_C, POS_C, ANG_C]),
             ("2  The rate (1294-1309)", "r[i] = max(0, D(θ[i]) × dt / 2 / dy²)",
              [(str(i), fmt(th[i], 1), f"{nums['D'][i]:.0f}", f"{r[i]:.3f}") for i in range(N)],
              ["i", "θ[i] °", "D(θ)", "r[i]"], [INK, ANG_C, INK, RATE_C]),
             ("3  The update (1315-1324)", "change[i] = r[i]·(x[i+1]−x[i]) + r[i]·(x[i−1]−x[i])",
              [(str(i), f"{r[i]:.3f}", fmt(north[i]), fmt(south[i])) for i in range(N)],
              ["i", "r[i]", "with north", "with south"], [INK, RATE_C, POS_C, POS_C]),
             ("4  Each boundary's two books", "south's book + north's book should be 0",
              [(f"{k}|{NEXT[k]}", fmt(north[k]), fmt(south[NEXT[k]]), fmt(north[k] + south[NEXT[k]])) for k in range(N)]
              + [("all", "", "", fmt((north + south).sum()))],
              ["boundary", "south", "north", "left over"], [INK, POS_C, POS_C, BAD])]
    for ax, (t, f, rows, hdr, cols) in zip(axes, steps):
        table(ax, rows, hdr, col_colors=cols, fs=11)
        ax.set_title(f"{t}\n{f}", fontsize=11.5, loc="left", family=MONO if False else None)
    fig.suptitle("BRIE's alongshore step on six sections (x = 1000, 1000, 900, 900, 900, 1000 m): the angle reads only the "
                 "north neighbour, so each section has one rate", fontsize=13, fontweight="bold")
    fig.text(0.5, 0.02, f"Full year from Brie.update(): the summed x grows by {nums['full'].sum():.1f} m "
             "(≈ 148,000 m³ of sand lost) with nothing removing any. Option A waves, 500 m sections.",
             ha="center", fontsize=11.5, color=BAD)
    fig.tight_layout(rect=(0, 0.05, 1, 0.92))
    png = workspace_twin(FIGS) / "brie_sand_equation_worked_example.png"
    save(fig, png, dpi=200, close=True)
    record_caption(png, (
        "**Tests:** nothing new; BRIE's alongshore step worked through by hand on six sections, every number from "
        "a real `brie.Brie` (barrier and inlet models off, option A waves, 500 m sections). **How to read:** "
        "left to right, the three lines of BRIE's step and the bookkeeping they imply. (1) each section's angle "
        "comes from x[i+1] − x[i] only: `x_s[np.r_[1:ny, 0]]` is the list shifted by one, the next section to "
        "the north. (2) the angle picks one diffusivity from BRIE's wave-climate table, and BRIE sets negative "
        "ones to zero. (3) the update uses the positions on both sides but that one rate. (4) at each boundary "
        "the two sections' books do not cancel where their rates differ. **Shows:** where the angle, and so the "
        "rate, comes from, and why the positions on both sides do not stop sand being made or lost."))

def main():
    apply_style()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig_both_sides()
    nums = brie_numbers()
    fig_worked_example(nums)
    gif_equation(nums)


if __name__ == "__main__":
    main()
