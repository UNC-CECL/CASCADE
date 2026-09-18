"""
coastsat_lrr_windows.py
==============================================================================
Observed shoreline change rate along the island, one panel per rate window,
every panel on the SAME y axis.

WHY
    Each window under coastsat_lrr/<start>_<end>/ ships its own
    domain_lrr_bar.png, autoscaled to that window. Put four of them side by
    side and a 2 m/yr swing in 1984-2004 is drawn as tall as a 7 m/yr swing in
    2010-2024, so the eye reads the wrong story. These figures pin one y range
    across the windows (Hannah, 2026-09-15: "the y axis among these plots must
    be held to the same bounds so they are easily comparable").

WHAT IS DRAWN
    The per-domain mean LRR (the model's grading target, `mean_lrr` in
    domain_lrr_summary.csv) as a LINE between domain centres over the 90 GIS
    domains (a step at the bin edges was tried 2026-09-15 and Hannah's
    advisor asked for a line). The line is coloured by sign and filled to
    zero: RdBu blue where the shore accreted, red where it eroded - the same
    reading as the per-window domain_lrr_bar.png, and Hannah's call over a
    BrBG pair (2026-09-15). Two dotted lines are +/- one standard deviation
    across the CoastSat transects inside each domain.
    Village spans are the house bands; the Buxton groin and the two piers are
    hairlines from the site config. Nothing else is on the canvas; the
    captions carry the method.

Y BOUNDS
    Symmetric: the largest |mean| over all windows plus a 1 m pad, rounded up
    to the next metre, mirrored about zero, ticks every 2 m/yr. The std lines
    are NOT in the bound - one wide domain at Cape Point was pushing every
    panel to +/-9 with nothing above 7 (the 2026-09-15 first cut). The value
    used is written to supporting/y_bounds.txt.

OUTPUT   data/hatteras_init/5-scr/coastsat_lrr_windows/
    lrr_<start>_<end>.png          one figure per window
    lrr_four_windows.png           2 x 2: the 1984-start period in the left
                                   column, the 1996-start period in the right
    pdf/<same stems>.pdf           the vector copies, in their own folder
                                   (Hannah, 2026-09-15)
    supporting/
        lrr_windows_wide.csv       the four means and stds side by side
        y_bounds.txt               the bounds every panel uses
        <figure>.pdf, CAPTIONS.md
    CAPTIONS.md                    written through hat_figure_style.caption()

USAGE
    python coastsat_lrr_windows.py                      # the four default windows
    python coastsat_lrr_windows.py --windows 1984_2004 2004_2024
==============================================================================
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

from site_layer.hat_figure_style import (  # noqa: E402
    C, DOMAIN_AXIS_LABEL, INK, INK_MUTED, STRUCTURE_LABEL_PT, apply_style,
    caption, figsize, open_frame, save, structures, support_dir, town_bands,
    _title,
)
from site_layer.hat_observed_rates import COASTSAT_LRR_WINDOWS, domain_csv  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS  # noqa: E402

# The 2 x 2 layout is by model period: each column is one CHAIN of windows,
# the second starting where the first ends (1984-2004 then 2004-2024 is the
# 1984-start period; 1996-2010 then 2010-2024 the 1996-start). Lettered across
# then down. Any window set that is not two chains of two falls back to one
# column.
DEFAULT_WINDOWS = [(1984, 2004), (1996, 2010), (2004, 2024), (2010, 2024)]
OUT_DIR = COASTSAT_LRR_WINDOWS

N_DOMAINS = 90
Y_LABEL = "Shoreline change rate (m/yr)"
Y_PAD_M = 1.0
Y_TICK_M = 2.0

# The RdBu poles and their light fills: red erosion, blue accretion, as the
# per-window bar charts already read. In this figure the pair means SIGN, not
# vintage; no vintage is drawn here, so the two readings never meet.
#
# The fills are the house light pair lightened by a third toward white
# (Hannah, 2026-09-15): at full strength the fill carried more weight than the
# lines over it, and on the comparison figure the black model line is what
# should read first. Both figures take these constants, so they move together.
FILL_LIGHTEN = 1 / 3


def _lighten(hex_colour: str, amount: float) -> str:
    """Mix a hex colour toward white by `amount` (0 = unchanged, 1 = white)."""
    r, g, b = (int(hex_colour[i:i + 2], 16) for i in (1, 3, 5))
    mix = lambda c: round(c + (255 - c) * amount)
    return "#{:02x}{:02x}{:02x}".format(mix(r), mix(g), mix(b))


C_ACCRETE = C["LATE"]
C_ACCRETE_FILL = _lighten(C["LATE_FILL"], FILL_LIGHTEN)
C_ERODE = C["EARLY"]
C_ERODE_FILL = _lighten(C["EARLY_FILL"], FILL_LIGHTEN)


# -----------------------------------------------------------------------------
# data
# -----------------------------------------------------------------------------
def load_window(start: int, end: int) -> pd.DataFrame:
    """domain 1..90 with mean_lrr / std_lrr; a domain with no fit is NaN."""
    df = pd.read_csv(domain_csv(start, end))
    df["domain_number"] = df["domain_number"].astype(int)
    full = pd.DataFrame({"domain_number": np.arange(1, N_DOMAINS + 1)})
    return full.merge(df[["domain_number", "mean_lrr", "std_lrr", "n_valid"]],
                      on="domain_number", how="left")


def shared_bounds(frames: list[pd.DataFrame]) -> float:
    """Half-range: the largest |mean| over every window plus the pad, rounded
    up to the next whole metre per year."""
    extreme = max(float(np.nanmax(df["mean_lrr"].abs())) for df in frames)
    return float(math.ceil(extreme + Y_PAD_M))


def signed_segments(x, y):
    """The line as (segment, colour) pairs. A segment that crosses zero is
    split where it crosses, so each half carries its own sign and the colour
    changes exactly where the fill does."""
    segs, cols = [], []
    for x0, y0, x1, y1 in zip(x[:-1], y[:-1], x[1:], y[1:]):
        if np.isnan(y0) or np.isnan(y1):
            continue
        if (y0 < 0) != (y1 < 0) and y0 != y1:
            xc = x0 + (x1 - x0) * (0.0 - y0) / (y1 - y0)
            segs += [[(x0, y0), (xc, 0.0)], [(xc, 0.0), (x1, y1)]]
            cols += [C_ERODE if y0 < 0 else C_ACCRETE,
                     C_ERODE if y1 < 0 else C_ACCRETE]
        else:
            segs.append([(x0, y0), (x1, y1)])
            cols.append(C_ERODE if (y0 + y1) / 2 < 0 else C_ACCRETE)
    return segs, cols


# -----------------------------------------------------------------------------
# one panel
# -----------------------------------------------------------------------------
STRUCTURE_LABEL_PT_GRID = 4.0  # the 2 x 2, whose panels are half the width


# structures() moved to hat_figure_style 2026-09-15 (shared with
# duneline_vs_coastsat); imported above.


def draw_panel(ax, df: pd.DataFrame, half: float, label: bool = True,
               label_pt: float = STRUCTURE_LABEL_PT, std: bool = True,
               line_lw: float = 1.0, fill_y=None, fill_outline_lw: float = 0.8):
    """The observed panel. `std`, `line_lw` and `fill_y` exist for the
    comparison figure that lays a scoring target over this
    (scripts/analyze_output/compare_runs/HAT_rate_windows.py):
    with `fill_y` given, THAT series takes the fill and a light outline, and
    the per-domain means are only the thin line over it, so the reference is
    the shape and the data the line (Hannah, 2026-09-15, option A)."""
    x = df["domain_number"].to_numpy(dtype=float)
    y = df["mean_lrr"].to_numpy(dtype=float)
    s = df["std_lrr"].fillna(0).to_numpy(dtype=float)

    ax.set_xlim(0.5, N_DOMAINS + 0.5)
    ax.set_ylim(-half, half)
    town_bands(ax, label=label)

    ax.axhline(0, color=INK_MUTED, lw=0.6, zorder=2)
    f = y if fill_y is None else np.asarray(fill_y, dtype=float)
    ok = ~np.isnan(f)
    ax.fill_between(x, 0, f, where=ok & (f >= 0), interpolate=True,
                    color=C_ACCRETE_FILL, lw=0, zorder=3)
    ax.fill_between(x, 0, f, where=ok & (f < 0), interpolate=True,
                    color=C_ERODE_FILL, lw=0, zorder=3)
    if fill_y is not None:
        segs, cols = signed_segments(x, f)
        ax.add_collection(LineCollection(segs, colors=cols,
                                         linewidths=fill_outline_lw,
                                         capstyle="round", zorder=5))
    if std:
        for band in (y - s, y + s):
            ax.plot(x, band, color=INK_MUTED, lw=0.5, ls=(0, (1, 1.6)),
                    zorder=4)
    segs, cols = signed_segments(x, y)
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=line_lw,
                                     capstyle="round", zorder=5))
    structures(ax, label, label_pt)

    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(Y_TICK_M))
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    open_frame(ax)


def caption_text(windows, half: float, grid: bool) -> str:
    wins = ", ".join(f"{a}–{b}" for a, b in windows)
    if grid:
        head = ("Observed shoreline change rate by GIS domain (1 at Cape "
                f"Point, 90 at Pea Island) for {wins}: the 1984-start period "
                "in the left column, the 1996-start period in the right, the "
                "earlier window of each above the later.")
    else:
        head = ("Observed shoreline change rate by GIS domain (1 at Cape "
                f"Point, 90 at Pea Island), {wins}.")
    return (head + " The line is the mean linear regression rate of the "
            "CoastSat transects inside each 500 m domain "
            "(domain_lrr_summary.csv), blue and filled where the shoreline "
            "moved seaward, red where it moved landward; the dotted lines "
            "are ±1 standard deviation across those transects. Village spans "
            "are shaded; the solid hairline is the Buxton groin and the "
            "dotted hairlines are the Avon and Rodanthe piers. The y axis is "
            f"held at ±{half:g} m/yr on every panel, the largest |mean| over "
            "the four windows plus 1 m rounded up, so the panels are directly "
            "comparable.")


# -----------------------------------------------------------------------------
# figures
# -----------------------------------------------------------------------------
def _save(fig, stem):
    """PNG in the folder, PDF under supporting/ -- the house save() does that
    itself since 2026-09-15 (it put the PDF beside the PNG, and this script
    kept a pdf/ folder of its own, before then)."""
    out = save(fig, OUT_DIR / stem)
    plt.close(fig)
    return out


def single_figure(start, end, df, half):
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.36),
                           constrained_layout=True)
    draw_panel(ax, df, half)
    ax.set_title(f"{start}–{end}", loc="center")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel(Y_LABEL)
    caption(fig, caption_text([(start, end)], half, grid=False))
    return _save(fig, f"lrr_{start}_{end}")


def _chains(windows):
    """Windows linked end-to-start, each chain sorted by start, chains by
    their first start: [(1984,2004),(2004,2024)], [(1996,2010),(2010,2024)]."""
    rest = sorted(windows)
    chains = []
    while rest:
        chain = [rest.pop(0)]
        while True:
            nxt = next((w for w in rest if w[0] == chain[-1][1]), None)
            if nxt is None:
                break
            rest.remove(nxt)
            chain.append(nxt)
        chains.append(chain)
    return chains


def grid_figure(windows, frames, half):
    """2 x 2 when the windows form two chains of two (a column per chain,
    the earlier window above); one column otherwise."""
    chains = _chains(windows)
    grid = len(chains) == 2 and all(len(c) == 2 for c in chains)
    if grid:
        nrow, ncol = 2, 2
        lookup = dict(zip(windows, frames))
        cells = [(r, c, (chain[r], lookup[chain[r]])) for r in range(2)
                 for c, chain in enumerate(chains)]
        fig, axes = plt.subplots(nrow, ncol, sharex=True, sharey=True,
                                 figsize=figsize("double", height=4.6),
                                 constrained_layout=True)
    else:
        nrow, ncol = len(windows), 1
        cells = [(i, 0, (w, f)) for i, (w, f) in enumerate(zip(windows, frames))]
        fig, axes = plt.subplots(nrow, ncol, sharex=True, sharey=True,
                                 figsize=figsize("double", height=1.55 * nrow + 0.6),
                                 constrained_layout=True, squeeze=False)
    axes = np.asarray(axes).reshape(nrow, ncol)
    for i, (r, c, ((start, end), df)) in enumerate(cells):
        ax = axes[r, c]
        draw_panel(ax, df, half, label=(i == 0),
                   label_pt=STRUCTURE_LABEL_PT_GRID)
        _title(ax, i, f"{start}–{end}")
        if c > 0:
            ax.tick_params(labelleft=False)
    for ax in axes[-1, :]:
        ax.set_xlabel(DOMAIN_AXIS_LABEL)
    fig.supylabel(Y_LABEL, fontsize=9)
    caption(fig, caption_text([w for _, _, (w, _) in cells], half, grid=grid))
    return _save(fig, "lrr_four_windows")


# -----------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 2)[1])
    ap.add_argument("--windows", nargs="+", metavar="START_END",
                    help="rate windows, e.g. 1984_2004 (default: the four)")
    args = ap.parse_args(argv)

    if args.windows:
        windows = [tuple(int(v) for v in w.split("_")) for w in args.windows]
    else:
        windows = DEFAULT_WINDOWS

    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    frames = [load_window(a, b) for a, b in windows]
    half = shared_bounds(frames)
    (support_dir(OUT_DIR) / "y_bounds.txt").write_text(
        f"y axis on every panel: -{half:g} to +{half:g} m/yr\n"
        f"= ceil(max |mean_lrr| + {Y_PAD_M:g}) over "
        + ", ".join(f"{a}-{b}" for a, b in windows) + "\n"
        "(the std lines are not in the bound)\n", encoding="utf-8")

    wide = pd.DataFrame({"domain_number": np.arange(1, N_DOMAINS + 1)})
    for (a, b), df in zip(windows, frames):
        wide[f"mean_{a}_{b}"] = df["mean_lrr"].round(3)
        wide[f"std_{a}_{b}"] = df["std_lrr"].round(3)
    wide.to_csv(support_dir(OUT_DIR) / "lrr_windows_wide.csv", index=False)

    written = []
    for (a, b), df in zip(windows, frames):
        written += single_figure(a, b, df, half)
    written += grid_figure(windows, frames, half)

    print(f"y bounds  +/-{half:g} m/yr")
    for p in written:
        print("wrote    ", p.relative_to(_REPO))


if __name__ == "__main__":
    main()
