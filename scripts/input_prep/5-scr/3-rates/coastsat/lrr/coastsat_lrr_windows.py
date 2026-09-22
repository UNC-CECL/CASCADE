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

OUTPUT   data/hatteras_init/5-scr/3-rates/coastsat/lrr/  (since 2026-09-19)
    ONLY the 2 x 2 is drawn now. Until 2026-09-19 this wrote
    4-comparisons/coastsat_windows/ with one folder per window and the
    --overlay halves figure too; those duplicated rates_figures.py's window
    and chain figures and were archived (--overlay is retired). The listing
    below is as it was.
    lrr_four_windows.png           2 x 2: the 1984-start period in the left
                                   column, the 1996-start period in the right
    supporting/
        lrr_windows_wide.csv       the four means and stds side by side
        y_bounds.txt               the bounds every panel uses
        lrr_four_windows.pdf, CAPTIONS.md
    <start>_<end>/lrr_<start>_<end>.png    one figure per window, its PDF and
                                           caption under its own supporting/

    --overlay 1996_2024 (2026-09-18) draws ONLY, into 1996_2024/:
    lrr_1996_2024_halves.png       two stacked panels: (a) the long window
                                   filled by sign, the model-input fill
                                   footprints as bars above it; (b) the two
                                   default windows that chain across it
                                   (1996-2010 grey, 2010-2024 black). No std
                                   lines; the y axis is the tightest whole
                                   metre holding every line, NOT the shared
                                   bound. Also published to
                                   output/figures/shoreline/
    supporting/lrr_1996_2024_halves.csv   the three means side by side
    The long window is context, not a grading target, so it is NOT added to
    the default four or to the 2 x 2.

USAGE
    python coastsat_lrr_windows.py                      # the four default windows
    python coastsat_lrr_windows.py --windows 1984_2004 2004_2024
    python coastsat_lrr_windows.py --overlay 1996_2024
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
    caption, figsize, figure_dir, open_frame, save, structures, support_dir,
    town_bands, _title,
)
from site_layer.hat_observed_rates import COASTSAT_LRR_WINDOWS, domain_csv  # noqa: E402
from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_ANNOTATIONS, HATTERAS_NOURISHMENT_PROJECTS,
)

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
# coastsat_vs_duneline); imported above.


def draw_panel(ax, df: pd.DataFrame, half: float, label: bool = True,
               label_pt: float = STRUCTURE_LABEL_PT, std: bool = True,
               line_lw: float = 1.0, fill_y=None, fill_outline_lw: float = 0.8):
    """The observed panel. `std`, `line_lw` and `fill_y` exist for the
    comparison figure that lays a scoring target over this
    (scripts/analyze_output/compare_runs/rate_windows.py):
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
    return _save(fig, f"{start}_{end}/lrr_{start}_{end}")


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
# the long window and its halves (--overlay)
# -----------------------------------------------------------------------------
# TWO STACKED PANELS (Hannah, 2026-09-18). The first cut laid the halves over
# the long window as dotted and dashed ink lines in one panel; they were more
# variable than the long window, crossed it everywhere, and pulled the eye off
# the quantity the figure is about. Now (a) is the long window alone, filled
# by sign, and (b) the two model periods as plain lines.
#
# The halves are a LUMINANCE pair, not the house vintage pair: panel (a)
# already spends red/blue on SIGN, and red "earlier" in (b) directly under red
# "landward" in (a) would read as one thing. Lighter grey is the earlier
# period, ink the later -- it survives greyscale and colour deficiency.
HALF_COLOURS = [C["BASE"], INK]
FILL_BAR_C = INK


def fills_in(start: int, end: int):
    """(year, first_gis, last_gis) of every ENABLED model-input fill placed
    inside the window, read from the site config the hindcast reads -- the
    model-input footprint, not the wider or narrower record span (Hannah,
    2026-09-18)."""
    return sorted((p.year, min(p.gis_domains), max(p.gis_domains))
                  for p in HATTERAS_NOURISHMENT_PROJECTS
                  if p.enabled and start <= p.year <= end)


def draw_fills(ax, fills, half: float, label_pt: float = STRUCTURE_LABEL_PT):
    """A bar just ABOVE the frame over each fill footprint, the year on it.
    A mark, not a shade: the village bands already shade. Above rather than
    along the bottom because the pier labels stand at the bottom, and Avon's
    ran through the 2022 label there (first cut, 2026-09-18).

    Placed in AXES fractions, not data units, so the same call works on the
    m/yr and the metre figures (duneline_windows.py reuses it); `half` is
    kept for the old callers and not used."""
    trans = ax.get_xaxis_transform()
    y = 1.025
    for year, lo, hi in fills:
        ax.plot([lo - 0.45, hi + 0.45], [y, y], color=FILL_BAR_C, lw=2.2,
                solid_capstyle="butt", zorder=6, clip_on=False,
                transform=trans)
        ax.text((lo + hi) / 2, y + 0.02, f"{year} fill", ha="center",
                va="bottom", fontsize=label_pt, color=FILL_BAR_C, zorder=6,
                clip_on=False, transform=trans)


# SHOALS (Hannah, 2026-09-18: "lighter so it doesn't take away from the
# shoreline change, I just want to see where it is"). Full-height BOXES, a
# thin amber outline over a sparse, faint amber hatch, with NO fill: the
# other shoal figures' solid wash would tint the sign fill in (a) and merge
# with the village greys, while a hatch stays readable over both. A thin
# bottom strip was tried first; Hannah asked for hatched boxes instead.
# The house amber (C["ADDED"]) the other alongshore figures use for shoals.
SHOAL_C = C["ADDED"]
SHOAL_HATCH = "///"
SHOAL_HATCH_ALPHA = 0.30
SHOAL_HATCH_LW = 0.5         # matplotlib 3.9: a global rc, read at draw
SHOAL_EDGE_ALPHA = 0.55
SHOAL_TEXT = "#8a620e"       # the label colour of the other shoal figures


def draw_shoals(ax, label: bool = True, label_pt: float = STRUCTURE_LABEL_PT):
    """Each shoal zone of the site config as a hatched, outlined box the
    full height of the panel, behind the data, named at the bottom when
    `label`. Hatch and outline are two patches because matplotlib 3.9 takes
    the hatch colour from the edge colour."""
    matplotlib.rcParams["hatch.linewidth"] = SHOAL_HATCH_LW
    for name, (lo, hi) in HATTERAS_ANNOTATIONS.shoal_zones.items():
        x0, w = lo - 0.5, (hi + 0.5) - (lo - 0.5)
        common = dict(transform=ax.get_xaxis_transform(), facecolor="none",
                      zorder=0.8, clip_on=True)
        ax.add_patch(plt.Rectangle((x0, 0.0), w, 1.0, hatch=SHOAL_HATCH,
                                   edgecolor=SHOAL_C, lw=0,
                                   alpha=SHOAL_HATCH_ALPHA, **common))
        ax.add_patch(plt.Rectangle((x0, 0.0), w, 1.0, edgecolor=SHOAL_C,
                                   lw=0.6, alpha=SHOAL_EDGE_ALPHA, **common))
        if label:
            ax.text((lo + hi) / 2, 0.02, name,
                    transform=ax.get_xaxis_transform(), ha="center",
                    va="bottom", fontsize=label_pt, color=SHOAL_TEXT,
                    zorder=1)


def halves_of(start: int, end: int):
    """The chain of default windows that runs start -> end end-to-start."""
    for chain in _chains(DEFAULT_WINDOWS):
        if chain[0][0] == start and chain[-1][1] == end:
            return chain
    raise SystemExit(f"no chain of default windows runs {start}-{end}; "
                     f"have {DEFAULT_WINDOWS}")


def tight_bound(frames) -> float:
    """The smallest whole metre per year that holds every line drawn, no pad.
    Hannah asked for a tighter axis than the shared +/-8 (2026-09-18) and
    named +/-6; 2010-2024 reaches +6.9 at GIS 1, so a fixed 6 would clip a
    measured value. The bound is computed so it can never clip."""
    return float(math.ceil(max(float(np.nanmax(df["mean_lrr"].abs()))
                               for df in frames)))


def overlay_caption(long_w, halves, half: float, shared: float, fills) -> str:
    (a, b), (h1, h2) = long_w, halves
    fill_txt = "; ".join(f"{y} at GIS {lo}–{hi}" for y, lo, hi in fills)
    return ("Observed shoreline change rate by GIS domain (1 at Cape Point, "
            f"90 at Pea Island). (a) {a}–{b}: the mean linear regression rate "
            "of the CoastSat transects inside each 500 m domain "
            "(domain_lrr_summary.csv), blue and filled where the shoreline "
            "moved seaward, red where it moved landward. (b) The same mean "
            f"over the two model periods, {h1[0]}–{h1[1]} (grey) and "
            f"{h2[0]}–{h2[1]} (black). The spread across transects is in the "
            "domain tables, not drawn. Black bars above (a) mark the beach "
            "fills placed inside the window, at the footprint the hindcast "
            f"uses ({fill_txt}); the rates there include the placed sand, and "
            "a fill is a step that a single slope fits poorly. The hatched "
            "amber boxes in both panels mark the offshore "
            "shoals (" + "; ".join(
                f"{n} GIS {lo}–{hi}" for n, (lo, hi)
                in HATTERAS_ANNOTATIONS.shoal_zones.items()) + "). "
            "Village spans are shaded; the solid hairline is the Buxton groin and the "
            "dotted hairlines are the Avon and Rodanthe piers. Both panels "
            f"share a y axis of ±{half:g} m/yr, the smallest whole metre that "
            "holds every value; the single-window figures use "
            f"±{shared:g}, so this figure is not read panel-for-panel against "
            f"them. The {a}–{b} rate is context: no model run is graded "
            "against it.")


def draw_halves(ax, halves, frames, half):
    """Panel (b): the two periods as solid lines over the same furniture as
    (a), unlabelled -- (a) carries the village and structure names."""
    ax.set_xlim(0.5, N_DOMAINS + 0.5)
    ax.set_ylim(-half, half)
    town_bands(ax, label=False)
    ax.axhline(0, color=INK_MUTED, lw=0.6, zorder=2)
    handles = []
    for (h0, h1), df, col in zip(halves, frames, HALF_COLOURS):
        (ln,) = ax.plot(df["domain_number"], df["mean_lrr"], color=col,
                        lw=1.1, zorder=5, label=f"{h0}–{h1}")
        handles.append(ln)
    structures(ax, label=False)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(Y_TICK_M))
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    open_frame(ax)
    return handles


def overlay_figure(long_w, halves, frames, half, shared):
    """(a) the long window filled by sign; (b) its halves as lines."""
    (a, b), (df_long, *df_halves) = long_w, frames
    fills = fills_in(a, b)
    fig, (ax_a, ax_b) = plt.subplots(
        2, 1, sharex=True, sharey=True, constrained_layout=True,
        figsize=figsize("double", height=4.9))
    draw_panel(ax_a, df_long, half, std=False)
    draw_fills(ax_a, fills, half)
    draw_shoals(ax_a, label=True)
    _title(ax_a, 0, f"{a}–{b}")
    handles = draw_halves(ax_b, halves, df_halves, half)
    draw_shoals(ax_b, label=False)
    _title(ax_b, 1, "The two model periods")
    ax_b.set_xlabel(DOMAIN_AXIS_LABEL)
    fig.supylabel(Y_LABEL, fontsize=9)
    fig.legend(handles=handles, loc="outside lower center", ncol=2,
               frameon=False)
    caption(fig, overlay_caption(long_w, halves, half, shared, fills))
    return _save_both(fig, f"{a}_{b}", f"lrr_{a}_{b}_halves"), fills


def _save_both(fig, window, stem):
    """Into the comparison folder, and published to output/figures/shoreline/
    (finished figures publish by subject, 2026-09-18). Each copy gets its
    CAPTIONS.md entry through the caption() wrapper."""
    out = save(fig, OUT_DIR / window / stem)
    out += save(fig, figure_dir("shoreline") / stem)
    plt.close(fig)
    return out


def run_overlay(spec: str):
    a, b = (int(v) for v in spec.split("_"))
    halves = halves_of(a, b)
    frames = [load_window(a, b)] + [load_window(*w) for w in halves]
    # The single-window figures' bound, named in the caption so a reader
    # knows this figure's axis differs from theirs.
    shared = shared_bounds([load_window(*w) for w in DEFAULT_WINDOWS])
    half = tight_bound(frames)

    table = pd.DataFrame({"domain_number": np.arange(1, N_DOMAINS + 1)})
    for (w0, w1), df in zip([(a, b)] + halves, frames):
        table[f"mean_{w0}_{w1}"] = df["mean_lrr"].round(3)
    table.to_csv(support_dir(OUT_DIR / f"{a}_{b}") / f"lrr_{a}_{b}_halves.csv",
                 index=False)

    written, fills = overlay_figure((a, b), halves, frames, half, shared)
    print(f"y bounds  +/-{half:g} m/yr (single-window figures: +/-{shared:g})")
    print("fills     " + ", ".join(f"{y} GIS {lo}-{hi}" for y, lo, hi in fills))
    for p in written:
        print("wrote    ", p.relative_to(_REPO))


# -----------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 2)[1])
    ap.add_argument("--windows", nargs="+", metavar="START_END",
                    help="rate windows, e.g. 1984_2004 (default: the four)")
    ap.add_argument("--overlay", metavar="START_END",
                    help="a long window drawn over the default windows that "
                         "chain across it, e.g. 1996_2024; draws only that")
    args = ap.parse_args(argv)

    if args.overlay:
        # RETIRED 2026-09-19: the halves overlay duplicated the 3-rates chain
        # figure (lrr/chains/lrr_chain_1996_2010_2024) and was archived.
        sys.exit("--overlay is retired; see 3-rates/coastsat/lrr/chains/")

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

    # Since 2026-09-19 only the 2 x 2 is drawn, into 3-rates/coastsat/lrr/
    # (OUT_DIR). The per-window figures are rates_figures.py's; drawing them
    # here too would overwrite 3-rates/coastsat/lrr/<w>/lrr_<w>.png.
    written = grid_figure(windows, frames, half)

    print(f"y bounds  +/-{half:g} m/yr")
    for p in written:
        print("wrote    ", p.relative_to(_REPO))


if __name__ == "__main__":
    main()
