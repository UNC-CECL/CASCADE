#!/usr/bin/env python3
"""
How much BRIE's step makes or loses depends on the coast's shape: how much the rate jumps between neighbours.

    python shape_figures.py [--years 40] [--no-gif]
        ->  figures/brie_sand_where_it_matters.png       straight, gentle, sharp: one number each
            figures/brie_sand_how_sharp.png              one year against the turn of a step
            figures/brie_sand_four_shapes.png            fillet, headland, bump, bay: books and totals
            figures/brie_sand_smooth_medium_sharp.gif    one feature three ways, 50 years (+ _end.png)
            figures/brie_sand_four_shapes_step_by_step.gif

Unmodified BRIE (alongshore step only) against one rate per boundary. Under these waves a flat boundary has rate 0,
so a slope meeting flat coast is a jump in rate even when the slope is gentle.

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
def _trapezoid(ramp):                          # 100 m plateau, ramps of `ramp` sections each side, centred on 16
    up = np.clip((I - (14 - ramp)) / (ramp + 1), 0, 1)
    down = np.clip(((18 + ramp) - I) / (ramp + 1), 0, 1)
    return 100.0 * np.minimum(up, down)


COASTS = {                                     # same height, from smooth to sharp
    "very smooth": 100 * np.exp(-0.5 * ((I - 16) / 3.5) ** 2),
    "medium": _trapezoid(2),
    "sharp": np.where((I >= 14) & (I < 19), 100.0, 0.0),
}
COAST_C = {"very smooth": GOOD, "medium": C["LATE"], "sharp": BAD}
FOUR_SHAPES = {                                   # seaward + (m)
    "groin-like fillet": np.select([I == 15, I == 16, I == 17, I == 18, I == 14, I == 13, I == 12],
                                   [100, 60, 30, 10, -60, -30, -10], 0.0),
    "sharp headland": np.where((I >= 15) & (I < 18), 100.0, 0.0),
    "gentle bump": 100 * np.exp(-0.5 * ((I - 15) / 3.0) ** 2),
    "bay": np.where((I >= 14) & (I < 18), -80.0, 0.0),
}
WINDOW = (8, 23)                             # domains drawn
MADE_C, LOST_C, SHARED_C = C["ADDED"], C["GROIN"], C["REF"]
HEIGHTS = np.arange(-200, 201, 10)           # step heights swept in brie_sand_how_sharp (m, seaward +)
WIN_LO, WIN_HI = 6, 27                       # sections drawn in the three-coast GIF
GIF = workspace_twin(FIGS) / "brie_sand_four_shapes_step_by_step.gif"
# -----------------------------------------------------------------------------


# 3. Where it matters: three coasts, one number each
def fig_where():
    shapes = (("straight coast", SHAPES["straight"]), ("gentle bump", SHAPES["gentle bump"]),
              ("sharp headland", HEADLAND))
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.8))
    for ax, (name, sh) in zip(axes, shapes):
        real, *_rest, b = run(sh)
        blank(ax)
        coast(ax, sh, lo=6, hi=25, base=-80)
        ax.set_ylim(-80, 260)
        lost = real[-1]
        col = GOOD if abs(lost) < 1 else (C["LATE"] if abs(lost) < 30 else BAD)
        ax.text(15.5, 205, f"{lost:+.0f} m".replace("-", "−") if abs(lost) >= 0.5 else "0 m", ha="center",
                fontsize=30, fontweight="bold", color=col)
        vol = abs(lost) * DY * depth(b)
        sub = "nothing made or lost" if abs(lost) < 0.5 else f"≈ {vol / 1e6:.2f} million m³ gone" if vol >= 1e6 \
            else f"≈ {vol / 1e3:,.0f},000 m³ gone"
        ax.text(15.5, 160, sub, ha="center", fontsize=11.5, color=col)
        ax.set_title(name, fontsize=14, pad=8)
    fig.suptitle(f"After {YEARS} years of BRIE with nothing adding or removing sand: the sharper the bend, "
                 "the bigger the error", fontsize=14, fontweight="bold", y=0.97)
    fig.text(0.5, 0.02, "Unmodified BRIE, alongshore step only. A straight coast is exact; a gentle bend nearly; "
             "a sharp bend (like a groin's fillet) is where it matters.", ha="center", fontsize=11, color=INK_MUTED)
    fig.tight_layout(rect=(0, 0.06, 1, 0.9))
    png = workspace_twin(FIGS) / "brie_sand_where_it_matters.png"
    save(fig, png, dpi=200, close=True)
    record_caption(png, (
        f"**Tests:** whether unmodified BRIE keeps the island's sand constant over {YEARS} years on three "
        "coasts (alongshore step only, option A waves, 500 m sections, nothing added or removed). **How to read:** "
        "each coast as a sketch, with the change in the island's summed shoreline after "
        f"{YEARS} years and the volume it stands for (× 500 m × BRIE's active depth). **Shows:** BRIE's sand "
        "budget holds on a straight coast and nearly on a gentle bend; at a sharp bend it does not. This is "
        "why colleagues who use BRIE on smooth coasts see it conserve sand, and why a groin, which makes a sharp "
        "bend, exposes the error."))


def fig_how_sharp():
    made = []
    for h in HEIGHTS:
        real, *_ = run(np.where((I >= 15) & (I < 18), float(h), 0.0), years=1)
        made.append(real[-1])
    made = np.array(made)
    angle = np.degrees(np.arctan2(HEIGHTS, DY))
    fig, ax = plt.subplots(figsize=figsize("single", height=3.0), constrained_layout=True)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.plot(angle, made, color=C["GROIN"], lw=1.8, marker="o", ms=2.5)
    ax.axvspan(-5, 5, color=SHARED_C, alpha=0.12, lw=0)
    ax.text(0, made.min() * 0.9, "gentle:\nnear zero", ha="center", fontsize=6.6, color=SHARED_C)
    ax.set_xlabel("Turn of the coast at the step's edges (°)\n(step height / 500 m domain; + seaward step)")
    ax.set_ylabel("Sand made (+) or lost (−)\nin one year (m of shoreline)")
    open_frame(ax)
    ax.set_title("The sharper the bend, the more sand is lost\n(a step on a coast facing the waves)", loc="center",
                 fontsize=8.5)
    png = workspace_twin(FIGS) / "brie_sand_how_sharp.png"
    save(fig, png, dpi=250, close=True)
    record_caption(png, (
        "**Tests:** how the sand BRIE's step makes or loses depends on how sharply the coast bends. **How to "
        "read:** each point is one year of an unmodified `brie.Brie` (alongshore step only, option A waves) "
        "started from a step three domains wide, of height −200 to +200 m (seaward +), plotted against the "
        "angle the coast turns at the step's edges. The shaded band is ±5°. **Shows:** a straight coast "
        "conserves exactly and a gently turning one nearly; the error grows roughly with the square of the turn. "
        "On this coast, facing the waves, a step either way loses sand; on other shapes and orientations the "
        "error can make sand instead (the groin-like fillet in brie_sand_four_shapes gains, and the fixed-shift groin made "
        "+539 m in 50 years at −20°). A groin's fillet or a sharp natural bend is where it matters."))


def turn(shape):
    return float(np.degrees(np.arctan(np.abs(np.diff(shape)).max() / DY)))


# 6. Three coasts, 50 years, BRIE as written beside one rate per boundary
def gif_smooth_medium_sharp():
    from PIL import Image
    runs = {}
    for name, shape in COASTS.items():
        b = brie_on(shape, YEARS)
        x = b._x_s.copy()
        real, sh, dis = [shape.copy()], [shape.copy()], []
        for _ in range(YEARS + 1):
            r, _ = rates(b, b._x_s)
            sea = seaward(b)
            dis.append((r - np.roll(r, -1)) * (np.roll(sea, -1) - sea) * 2)   # boundary i|i+1, this year's books
            if len(real) > YEARS:
                break
            b.update()
            x = shared_step(b, x)
            real.append(seaward(b).copy())
            sh.append(seaward(b, x).copy())
        runs[name] = (np.array(real), np.array(sh), np.array(dis), b)
    dmax = max(np.abs(v[2]).max() for v in runs.values())
    frames, durs = [], []

    def grab(fig, ms):
        fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert("RGB"))
        durs.append(ms)
        plt.close(fig)

    for y in range(YEARS + 1):
        fig, axes = plt.subplots(3, 1, figsize=(10.5, 7.6), dpi=88)
        fig.subplots_adjust(left=0.02, right=0.72, top=0.86, bottom=0.03, hspace=0.42)
        fig.text(0.03, 0.965, "The same 100 m feature, three shapes, unmodified BRIE", fontsize=16, fontweight="bold",
                 va="top")
        fig.text(0.03, 0.915, f"year {y} of {YEARS}   ·   solid: BRIE as written   ·   green: one rate per boundary   "
                 "·   red dots: the two sides disagree", fontsize=10.5, color=INK_MUTED, va="top")
        for ax, (name, (real, sh, dis, b)) in zip(axes, runs.items()):
            blank(ax)
            coast(ax, real[y], lo=WIN_LO, hi=WIN_HI, base=-40, outline=COAST_C[name])
            cells = np.arange(WIN_LO, WIN_HI)
            xs = np.repeat(np.r_[cells - 0.5, cells[-1] + 0.5], 2)[1:-1]
            ax.plot(xs, np.repeat(COASTS[name][cells], 2), color=INK_MUTED, lw=0.9, ls=(0, (3, 2)))
            ax.plot(xs, np.repeat(sh[y][cells], 2), color=GOOD, lw=1.3, alpha=0.9)
            d = dis[min(y, len(dis) - 1)]
            for i in cells[:-1]:
                if abs(d[i]) > 0.02 * dmax:
                    ax.scatter(i + 0.5, max(real[y][i], real[y][i + 1]) + 18, s=20 + 380 * abs(d[i]) / dmax,
                               color=BAD, alpha=0.85, zorder=5, edgecolor="white", lw=0.6)
            ax.set_ylim(-40, 150)
            ax.set_title(f"{name}  (steepest turn {turn(COASTS[name]):.0f}°)", fontsize=12.5, loc="left",
                         color=COAST_C[name], fontweight="bold")
            tot = real[y].sum() - real[0].sum()
            pos = ax.get_position()
            col = GOOD if abs(tot) < 0.5 else COAST_C[name] if name != "very smooth" else C["LATE"]
            fig.text(0.86, pos.y0 + pos.height * 0.62, "0 m" if abs(tot) < 0.05 else f"{tot:+.1f} m".replace("-", "−"),
                     ha="center", fontsize=24, fontweight="bold", color=BAD if tot < -1 else col)
            fig.text(0.86, pos.y0 + pos.height * 0.25, "sand made or lost\n(green line: 0 m)", ha="center",
                     fontsize=9.5, color=INK_MUTED)
        grab(fig, 2200 if y == 0 else 5000 if y == YEARS else 120)
    out = workspace_twin(FIGS) / "brie_sand_smooth_medium_sharp.gif"
    frames[0].save(out, save_all=True, append_images=frames[1:], duration=durs, loop=0, optimize=True)
    frames[-1].save(workspace_twin(FIGS) / "brie_sand_smooth_medium_sharp_end.png")
    print("wrote", out, f"| {len(frames)} frames")
    for name, (real, *_rest) in runs.items():
        print(f"{name:12s} steepest turn {turn(COASTS[name]):4.1f} deg  year {YEARS}: {real[-1].sum() - real[0].sum():+8.2f} m")


# Each boundary's first-order books for one year: south side, north side, and their sum (made +, lost -)
def books(b, x, shared=False):
    r, _ = rates(b, x)
    s = seaward(b, x)
    d = np.roll(s, -1) - s                         # boundary i|i+1
    rn = r if shared else np.roll(r, -1)
    south, north = 2 * r * d, -2 * rn * d
    return south + 0.0, north + 0.0, south + north + 0.0


# Year-by-year states for one shape: BRIE as written (real Brie.update) and one rate per boundary
def history(shape, years):
    b = brie_on(shape, years)
    x = b._x_s.copy()
    out = [(seaward(b).copy(), seaward(b, x).copy(), books(b, b._x_s)[2], books(b, x, True)[2])]
    for _ in range(years):
        b.update()
        x = shared_step(b, x)
        out.append((seaward(b).copy(), seaward(b, x).copy(), books(b, b._x_s)[2], books(b, x, True)[2]))
    return out, b


def steps(ax, s, **kw):
    xs = np.repeat(np.r_[I - 0.5, I[-1] + 0.5], 2)[1:-1]
    return ax.plot(xs, np.repeat(s, 2), **kw)


def coast_fill(ax, s):
    xs = np.repeat(np.r_[I - 0.5, I[-1] + 0.5], 2)[1:-1]
    ys = np.repeat(s, 2)
    ax.fill_between(xs, -400, ys, color=C["ADDED_FILL"], lw=0)
    ax.fill_between(xs, ys, 400, color=C["WATER"], alpha=0.45, lw=0)


def book_bars(ax, made, lim):
    bnd = I + 0.5
    keep = (bnd > WINDOW[0]) & (bnd < WINDOW[1])
    ax.bar(bnd[keep], made[keep], width=0.55, color=[MADE_C if m > 0 else LOST_C for m in made[keep]])
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_ylim(-lim, lim)
    ax.set_xlim(*WINDOW)


def fig_four(hist, years):
    names = list(FOUR_SHAPES)
    fig = plt.figure(figsize=figsize("double", height=7.4), constrained_layout=True)
    gs = fig.add_gridspec(3, len(names), height_ratios=(1, 0.75, 0.8))
    k = 0
    blim = 1.15 * max(np.abs(h[0][2]).max() for h, _ in hist.values())
    tlim = 1.15 * max(abs(h[-1][0].sum() - h[0][0].sum()) for h, _ in hist.values())
    for j, name in enumerate(names):
        h, b = hist[name]
        ax = fig.add_subplot(gs[0, j])
        coast_fill(ax, h[0][0])
        steps(ax, h[0][0], color=INK_MUTED, lw=1.0, ls=(0, (3, 2)), label="start")
        steps(ax, h[-1][1], color=SHARED_C, lw=1.3, label="one rate per boundary")
        steps(ax, h[-1][0], color=INK, lw=1.6, label="BRIE as written")
        ax.set_xlim(*WINDOW)
        ax.set_ylim(-110, 140)
        ax.set_xticks([])
        if j == 0:
            ax.set_ylabel("Shoreline (m), + seaward")
            ax.legend(frameon=False, fontsize=5.8, loc="upper left")
        else:
            ax.set_yticklabels([])
        open_frame(ax)
        _letter_inside(ax, k)
        k += 1
        ax.set_title(f"{name}\nstart and year {years}", loc="center", fontsize=8)
        ax = fig.add_subplot(gs[1, j])
        book_bars(ax, h[0][2], blim)
        if j == 0:
            ax.set_ylabel("Year 1: books disagree\nby (m; + made, − lost)")
        else:
            ax.set_yticklabels([])
        ax.set_xlabel("boundary between domains", fontsize=6.5)
        ax.tick_params(axis="x", labelsize=6)
        open_frame(ax)
        _letter_inside(ax, k)
        k += 1
        ax = fig.add_subplot(gs[2, j])
        yrs = np.arange(years + 1)
        tot = np.array([s[0].sum() for s in h]) - h[0][0].sum()
        tsh = np.array([s[1].sum() for s in h]) - h[0][1].sum()
        ax.plot(yrs, tsh, color=SHARED_C, lw=1.3)
        ax.plot(yrs, tot, color=LOST_C if tot[-1] < 0 else MADE_C, lw=1.8)
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.text(years * 0.97, tot[-1] + (0.12 if tot[-1] >= 0 else -0.12) * tlim,
                f"{tot[-1]:+.1f} m".replace("-", "−"), fontsize=7, ha="right", va="center",
                color=LOST_C if tot[-1] < 0 else MADE_C, fontweight="bold")
        ax.set_ylim(-tlim, tlim)
        ax.set_xlabel("Year")
        if j == 0:
            ax.set_ylabel("Island's sand made (+)\nor lost (−), m of shoreline")
        else:
            ax.set_yticklabels([])
        open_frame(ax)
        _letter_inside(ax, k)
        k += 1
    fig.suptitle("Unmodified BRIE on four coastline shapes: the books disagree at every bend, and the total drifts",
                 fontsize=10)
    png = workspace_twin(FIGS) / "brie_sand_four_shapes.png"
    save(fig, png, dpi=250, close=True)
    tot = {n: hist[n][0][-1][0].sum() - hist[n][0][0][0].sum() for n in names}
    record_caption(png, (
        "**Tests:** whether BRIE's alongshore step keeps the island's sand constant on four coastline shapes, "
        "with an unmodified `brie.Brie` (barrier and inlet models off, so only the alongshore step moves the "
        f"shoreline; option A waves; 30 domains of 500 m; nothing added or removed; {years} years). **How to "
        "read:** top, each shape at the start (dashed) and after "
        f"{years} years, BRIE as written (black) and with one rate per boundary (green). Middle, year 1 at every "
        "boundary: what the domain to the south books across it minus what the domain to the north books "
        "(first order, 2 r Δ each); gold is sand made, red sand lost, and an empty boundary is one where the two "
        "agree. Bottom, the island's summed shoreline against its start (green: one rate per boundary, zero "
        "throughout). **Shows:** wherever the coast bends, the two domains at a boundary use different rates and "
        "the books disagree; along straight stretches they agree. Over the years, "
        + ", ".join(f"the {n} {'gains' if v > 0 else 'loses'} {abs(v):.1f} m" for n, v in tot.items())
        + ". With one rate per boundary all four conserve to round-off."))


def gif_four_shapes(hist, years):
    from PIL import Image
    names = list(FOUR_SHAPES)
    blim = 1.15 * max(np.abs(np.array([s[2] for s in h])).max() for h, _ in hist.values())
    fig = plt.figure(figsize=(9.0, 6.2), dpi=95)
    gs = fig.add_gridspec(2, 2, width_ratios=(1.6, 1), height_ratios=(1, 0.62), left=0.08, right=0.98, top=0.82,
                          bottom=0.09, hspace=0.35, wspace=0.22)
    imgs, durs = [], []
    for name in names:
        h, _ = hist[name]
        tot = np.array([s[0].sum() for s in h]) - h[0][0].sum()
        tsh = np.array([s[1].sum() for s in h]) - h[0][1].sum()
        tl = 1.2 * max(5.0, np.abs(tot).max())
        for y in range(years + 1):
            fig.clear()
            fig.text(0.08, 0.95, f"Unmodified BRIE, alongshore step only: a {name}", fontsize=12, fontweight="bold",
                     va="top")
            fig.text(0.08, 0.905, f"year {y} of {years}   ·   nothing adds or removes sand, so the island's total "
                     "should stay at zero", fontsize=9, color=INK_MUTED, va="top")
            s_real, s_sh, made, _ = h[y]
            ax = fig.add_subplot(gs[0, 0])
            coast_fill(ax, s_real)
            steps(ax, h[0][0], color=INK_MUTED, lw=1.0, ls=(0, (3, 2)))
            steps(ax, s_sh, color=SHARED_C, lw=1.4)
            steps(ax, s_real, color=INK, lw=2.0)
            for bnd in I[(I + 0.5 > WINDOW[0]) & (I + 0.5 < WINDOW[1])] + 0.5:
                ax.axvline(bnd, color=INK_MUTED, lw=0.5, ls=(0, (1, 2)), alpha=0.6)
            ax.set_xlim(*WINDOW)
            ax.set_ylim(-110, 140)
            ax.set_xticks(np.arange(WINDOW[0] + 1, WINDOW[1]), [str(i) for i in np.arange(WINDOW[0] + 1, WINDOW[1])],
                          fontsize=7)
            ax.set_ylabel("Shoreline (m), + seaward")
            ax.set_title("black: BRIE as written · green: one rate per boundary · dashed: start", fontsize=8.5,
                         loc="left")
            open_frame(ax)
            ax = fig.add_subplot(gs[1, 0])
            book_bars(ax, made, blim)
            ax.set_ylabel("books disagree by\n(m this year)", fontsize=8)
            ax.set_xlabel("boundary between domains (dotted lines above), south → north", fontsize=8)
            ax.set_title("at each boundary: south domain's book minus north domain's (gold made, red lost)",
                         fontsize=8.5, loc="left")
            ax.tick_params(labelsize=7)
            open_frame(ax)
            ax = fig.add_subplot(gs[:, 1])
            ax.plot(np.arange(y + 1), tsh[:y + 1], color=SHARED_C, lw=1.6, label="one rate per boundary")
            ax.plot(np.arange(y + 1), tot[:y + 1], color=LOST_C if tot[y] < 0 else MADE_C, lw=2.2,
                    label="BRIE as written")
            ax.axhline(0, color=INK_MUTED, lw=0.7)
            ax.set_xlim(0, years)
            ax.set_ylim(-tl, tl)
            ax.text(0.5, 0.95, f"{tot[y]:+.1f} m".replace("-", "−"), transform=ax.transAxes, ha="center", va="top",
                    fontsize=16, fontweight="bold", color=LOST_C if tot[y] < -0.05 else MADE_C if tot[y] > 0.05 else SHARED_C)
            ax.set_xlabel("Year", fontsize=8)
            ax.set_ylabel("Island's sand made (+) or lost (−)\n(m of summed shoreline)", fontsize=8)
            ax.legend(frameon=False, fontsize=7.5, loc="lower left")
            ax.set_title("the island's total", fontsize=8.5, loc="left")
            ax.tick_params(labelsize=7)
            open_frame(ax)
            fig.canvas.draw()
            imgs.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert("RGB"))
            durs.append(1400 if y == 0 else 2600 if y == years else 160 if y > 3 else 700)
    plt.close(fig)
    GIF.parent.mkdir(parents=True, exist_ok=True)
    imgs[0].save(GIF, save_all=True, append_images=imgs[1:], duration=durs, loop=0, optimize=True)
    print("wrote", GIF, f"| {len(imgs)} frames")


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--years", type=int, default=40, help="years for the four-shapes figure and GIF")
    ap.add_argument("--no-gif", action="store_true")
    a = ap.parse_args()
    apply_style()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig_where()
    fig_how_sharp()
    hist = {name: history(shape, a.years) for name, shape in FOUR_SHAPES.items()}
    fig_four(hist, a.years)
    if not a.no_gif:
        gif_smooth_medium_sharp()
        gif_four_shapes(hist, a.years)


if __name__ == "__main__":
    main()
