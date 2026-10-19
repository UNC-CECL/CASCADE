#!/usr/bin/env python3
"""
BRIE's sand budget over 50 years: the headland told step by step and then run, and the total on three coasts.

    python over_time_figures.py   ->  figures/brie_sand_total_over_time.png
                                      figures/brie_sand_headland_50_years.gif   (git-ignored)

Unmodified BRIE (alongshore step only) against the same coast with one rate per boundary. The straight coast
keeps its sand exactly, a gentle bump loses 10.6 m in 50 years and a sharp 100 m headland 117.9 m.

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


def fig_real_brie(results):
    fig = plt.figure(figsize=figsize("double", height=3.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=(1, 1, 1.15))
    for j, name in enumerate(("gentle bump", "sharp step")):
        ax = fig.add_subplot(gs[0, j])
        real, shared, prof, _ = results[name]
        xs = np.repeat(np.r_[I - 0.5, I[-1] + 0.5], 2)[1:-1]
        ax.plot(xs, np.repeat(prof[0], 2), color=INK_MUTED, lw=1.0, ls=(0, (3, 2)), label="start")
        ax.plot(xs, np.repeat(prof["shared"], 2), color=SHARED_C, lw=1.4, label=f"year {YEARS}, one rate per boundary")
        ax.plot(xs, np.repeat(prof[YEARS], 2), color=SHAPE_C[name], lw=1.6, label=f"year {YEARS}, BRIE as written")
        ax.set_xlim(6, 26)
        ax.set_ylim(-10, 110)
        ax.set_xlabel("domain (500 m each), south → north")
        if j == 0:
            ax.set_ylabel("Shoreline (m), + seaward")
            ax.legend(frameon=False, fontsize=6.3, loc="upper left")
        open_frame(ax)
        _letter_inside(ax, j)
        ax.set_title(f"{name.capitalize()}: start and year {YEARS}", loc="center", fontsize=8.5)
    ax = fig.add_subplot(gs[0, 2])
    yrs = np.arange(YEARS + 1)
    for name in SHAPES:
        real, shared, _, b = results[name]
        ax.plot(yrs, real, color=SHAPE_C[name], lw=1.8, label=f"{name}, BRIE as written")
    ax.plot(yrs, results["sharp step"][1], color=SHARED_C, lw=1.4, ls=(0, (4, 2)), label="sharp step, one rate per boundary")
    b = results["sharp step"][3]
    lost = results["sharp step"][0][-1]
    vol = lost * DY * depth(b)
    ax.annotate(f"{lost:+.0f} m of shoreline\n≈ {vol / 1e6:+.1f} million m³", xy=(YEARS, lost), xytext=(22, lost * 0.55),
                fontsize=7, color=C["GROIN"], arrowprops=dict(arrowstyle="-", color=C["GROIN"], lw=0.6))
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlabel("Year")
    ax.set_ylabel("Sand made (+) or lost (−)\n(m of summed shoreline)")
    ax.legend(frameon=False, fontsize=6.3, loc="lower left")
    open_frame(ax)
    _letter_inside(ax, 2)
    ax.set_title("The island's total sand, nothing added", loc="center", fontsize=8.5)
    fig.suptitle("Unmodified BRIE, alongshore step only: where its sand budget holds and where it does not", fontsize=10)
    png = FIGS / "brie_sand_total_over_time.png"
    save(fig, png, dpi=250, close=True)
    g = results["gentle bump"][0][-1]
    record_caption(png, (
        "**Tests:** whether BRIE's alongshore step keeps the island's sand constant when nothing adds or removes "
        "any, using an unmodified `brie.Brie` with its barrier and inlet models off, so `update()` moves the "
        "shoreline by the alongshore step alone. Option A waves, 30 domains of 500 m on BRIE's periodic coast, "
        f"{YEARS} one-year steps. **How to read:** (a) and (b) the shoreline at the start (dashed), after {YEARS} "
        "years of BRIE as written, and after the same years with each boundary's rate shared by the two domains it "
        "joins (green). (c) the island's summed shoreline against its start: zero means sand is conserved. "
        f"The volume uses BRIE's active depth ({depth(b):.1f} m) and 500 m domains. **Shows:** on a straight coast "
        f"nothing is made or lost (0.000 m); a gentle 100 m bump loses {g:.1f} m in {YEARS} years; a sharp 100 m "
        f"step loses {lost:.1f} m, {abs(lost) / 300:.0%} of the sand in it. With one rate per boundary every case "
        "conserves to round-off, so the loss comes from the two domains at a bend using different rates for the "
        "boundary between them, not from the waves or the equation."))


# 4. The GIF: the story told in stages, then 50 years side by side
def gif_headland():
    from PIL import Image
    b, r, theta, s, edges = headland_numbers()
    real, shared, _, _ = run(HEADLAND)
    hist_real, hist_sh = [s.copy()], [s.copy()]
    bb = brie_on(HEADLAND, YEARS)
    x = bb._x_s.copy()
    for _ in range(YEARS):
        bb.update()
        x = shared_step(bb, x)
        hist_real.append(seaward(bb).copy())
        hist_sh.append(seaward(bb, x).copy())
    (s14, _), (s17, _) = edges[14], edges[17]
    frames, durs = [], []

    def grab(fig, ms):
        fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert("RGB"))
        durs.append(ms)

    def stage(headline, sub, draw, ms):
        fig, ax = plt.subplots(figsize=(10, 6.2), dpi=90)
        fig.subplots_adjust(left=0.02, right=0.98, top=0.83, bottom=0.04)
        blank(ax)
        coast(ax, s, base=-160)
        ax.set_ylim(-160, 330)
        ax.text(HI - 0.6, -150, "north →", fontsize=11, color=INK_MUTED, ha="right", va="bottom")
        fig.text(0.03, 0.95, headline, fontsize=17, fontweight="bold", va="top")
        fig.text(0.03, 0.885, sub, fontsize=12, color=INK_MUTED, va="top")
        draw(ax)
        grab(fig, ms)
        plt.close(fig)

    def rates_drawn(ax, upto=HI):
        for i in range(LO, min(upto, HI)):
            fast = r[i] > 1e-9
            ax.add_patch(Rectangle((i - 0.42, -120), 0.84, 32, fc=GAIN if fast else "white", ec=INK_MUTED, lw=0.8))
            ax.text(i, -104, f"{1000 * r[i]:.0f}" if fast else "0", ha="center", va="center", fontsize=10,
                    color="white" if fast else INK_MUTED, fontweight="bold")
        ax.text(LO - 0.4, -135, "rate per 1,000, from each section's NORTH boundary", fontsize=10, va="top")

    stage("A headland on BRIE's coast", "Nothing adds or removes sand. The island's total should never change.",
          lambda ax: ax.text(16, 125, "100 m headland", ha="center", fontsize=12, color=INK_MUTED), 2600)
    for upto in range(LO + 1, HI + 1):
        i = upto - 1

        def look(ax, i=i, upto=upto):
            rates_drawn(ax, upto)
            arrow(ax, i, 150, i + 1, 150, col=GAIN if r[i] > 1e-9 else INK_MUTED, lw=2.2, ms=16)
            ax.text(i + 0.5, 175, "looks north", ha="center", fontsize=9.5, color=INK_MUTED)
        stage("Step 1: each section picks ONE rate", "from the angle of the boundary on its north side only",
              look, 650 if upto < HI else 1800)
    stage("Step 2: it uses that rate on BOTH its sides", "so the two sections at a boundary can use different "
          "rates for the same handover", lambda ax: rates_drawn(ax), 2600)

    def edge_south(ax):
        rates_drawn(ax)
        bubble(ax, 13.3, 230, f"section 14:\n\"{s14:.1f} m arrives\"", fs=11)
        bubble(ax, 15.6, 230, "section 15:\n\"nothing leaves\"", fs=11)
        ax.text(14.5, 70, f"+{s14:.1f} m\nappears", ha="center", fontsize=14, fontweight="bold", color=GAIN,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=GAIN, lw=1.5))
    stage("Step 3: the books at the south edge", "section 14 sees a steep north side; section 15 sees a flat one",
          edge_south, 3400)

    def edge_north(ax):
        edge_south(ax)
        bubble(ax, 17.0, 230, f"section 17:\n\"{-s17:.1f} m leaves\"", fs=11, ec=BAD)
        bubble(ax, 19.2, 230, "section 18:\n\"nothing arrives\"", fs=11, ec=BAD)
        ax.text(17.5, 70, f"{s17:.1f} m\nvanishes".replace("-", "−"), ha="center", fontsize=14, fontweight="bold",
                color=BAD, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=BAD, lw=1.5))
    stage("Step 4: the books at the north edge", "the same thing, the other way round",
          edge_north, 3800)
    vol1 = abs(real[1]) * DY * depth(b)
    stage(f"Year 1: the island has {abs(real[1]):.0f} m of shoreline less",
          f"about {round(vol1, -3):,.0f} m³ of sand, although nothing removed any", edge_north, 3400)
    # 50 years side by side
    for y in list(range(0, YEARS + 1)):
        fig, axes = plt.subplots(2, 1, figsize=(10, 6.2), dpi=90)
        fig.subplots_adjust(left=0.02, right=0.74, top=0.78, bottom=0.04, hspace=0.32)
        for ax, h, lab, col, tot in ((axes[0], hist_real, "BRIE as written", BAD, real),
                                     (axes[1], hist_sh, "if both sides used the same rate", GOOD, shared)):
            blank(ax)
            coast(ax, h[y], base=-60, outline=col)
            ax.plot(np.repeat(np.r_[np.arange(LO, HI) - 0.5, HI - 0.5], 2)[1:-1], np.repeat(s[LO:HI], 2),
                    color=INK_MUTED, lw=1, ls=(0, (3, 2)))
            ax.set_ylim(-60, 130)
            ax.set_title(lab, fontsize=13, loc="left", color=col, fontweight="bold")
            v = tot[y]
            fig.text(0.86, ax.get_position().y0 + ax.get_position().height * 0.55,
                     ("0 m" if abs(v) < 0.05 else f"{v:+.0f} m".replace("-", "−")), ha="center", fontsize=26,
                     fontweight="bold", color=col)
            fig.text(0.86, ax.get_position().y0 + ax.get_position().height * 0.25, "sand made or lost",
                     ha="center", fontsize=10.5, color=INK_MUTED)
        fig.text(0.03, 0.95, f"Year {y} of {YEARS}", fontsize=17, fontweight="bold", va="top")
        fig.text(0.03, 0.885, "dashed: the start. The same coast, the same waves; only the bookkeeping differs.",
                 fontsize=12, color=INK_MUTED, va="top")
        grab(fig, 1500 if y == 0 else 4500 if y == YEARS else 110)
        plt.close(fig)
    vol = abs(real[-1]) * DY * depth(b)
    fig, ax = plt.subplots(figsize=(10, 6.2), dpi=90)
    blank(ax)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.text(5, 7.4, f"After {YEARS} years: {abs(real[-1]):.0f} m of shoreline lost", ha="center", fontsize=22,
            fontweight="bold", color=BAD)
    ax.text(5, 5.9, f"≈ {vol / 1e6:.1f} million m³ of sand, about {round(vol / TRUCK_M3, -3):,.0f} truckloads,\n"
            f"from {abs(real[-1]) / 3:.0%} of a 300 m headland... with nothing removing any.",
            ha="center", fontsize=14, linespacing=1.5)
    ax.text(5, 3.3, "BRIE's equation conserves sand. Its step does too, on straight and gently curving coasts.\n"
            "At a sharp bend, like a groin's fillet, the two sides of a boundary use different rates,\n"
            "and the difference is sand made or lost.", ha="center", fontsize=12.5, color=INK, linespacing=1.6)
    grab(fig, 6000)
    plt.close(fig)
    out = FIGS / "brie_sand_headland_50_years.gif"
    frames[0].save(out, save_all=True, append_images=frames[1:], duration=durs, loop=0, optimize=True)
    print("wrote", out, f"| {len(frames)} frames")

def main():
    apply_style()
    FIGS.mkdir(parents=True, exist_ok=True)
    results = {name: run(shape) for name, shape in SHAPES.items()}
    fig_real_brie(results)
    gif_headland()


if __name__ == "__main__":
    main()
