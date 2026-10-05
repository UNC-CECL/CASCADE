"""
Observed shoreline change rate along the island, one panel per window, every panel on the same y axis.

    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_windows.py
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_windows.py --windows 1984_2004 2004_2024

Writes the multi-window figure and one per window, with captions. Details: scripts/input_prep/5-scr/3-rates/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-03
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

# --- CONFIG ------------------------------------------------------------------
# The 2 x 2 layout is by model period
DEFAULT_WINDOWS = [(1984, 2004), (1996, 2010), (2004, 2024), (2010, 2024)]
OUT_DIR = COASTSAT_LRR_WINDOWS

N_DOMAINS = 90
Y_LABEL = "Shoreline change rate (m/yr)"
Y_PAD_M = 1.0
Y_TICK_M = 2.0

# The RdBu poles and their light fills
FILL_LIGHTEN = 1 / 3
# -----------------------------------------------------------------------------


# Mix a hex colour toward white by `amount` (0 = unchanged, 1 = white)
def _lighten(hex_colour: str, amount: float) -> str:
    r, g, b = (int(hex_colour[i:i + 2], 16) for i in (1, 3, 5))
    mix = lambda c: round(c + (255 - c) * amount)
    return "#{:02x}{:02x}{:02x}".format(mix(r), mix(g), mix(b))


C_ACCRETE = C["LATE"]
C_ACCRETE_FILL = _lighten(C["LATE_FILL"], FILL_LIGHTEN)
C_ERODE = C["EARLY"]
C_ERODE_FILL = _lighten(C["EARLY_FILL"], FILL_LIGHTEN)


# domain 1..90 with mean_lrr / std_lrr
def load_window(start: int, end: int) -> pd.DataFrame:
    df = pd.read_csv(domain_csv(start, end))
    df["domain_number"] = df["domain_number"].astype(int)
    full = pd.DataFrame({"domain_number": np.arange(1, N_DOMAINS + 1)})
    return full.merge(df[["domain_number", "mean_lrr", "std_lrr", "n_valid"]],
                      on="domain_number", how="left")


# Half-range
def shared_bounds(frames: list[pd.DataFrame]) -> float:
    extreme = max(float(np.nanmax(df["mean_lrr"].abs())) for df in frames)
    return float(math.ceil(extreme + Y_PAD_M))


# The line as (segment, colour) pairs
def signed_segments(x, y):
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


# One panel

STRUCTURE_LABEL_PT_GRID = 4.0  # the 2 x 2, whose panels are half the width


# structures() lives in hat_figure_style since 2026-09-15; imported above


# The observed panel
def draw_panel(ax, df: pd.DataFrame, half: float, label: bool = True,
               label_pt: float = STRUCTURE_LABEL_PT, std: bool = True,
               line_lw: float = 1.0, fill_y=None, fill_outline_lw: float = 0.8):
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


# The caption for a set of windows
def caption_text(windows, half: float, grid: bool,
                 bound_over: str = "the four windows") -> str:
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
            f"{bound_over} plus 1 m rounded up, so the panels are directly "
            "comparable.")


# Figures

# PNG in the folder and PDF under supporting/, which the house save() does itself
def _save(fig, stem):
    out = save(fig, OUT_DIR / stem)
    plt.close(fig)
    return out


# One window's panel as its own figure
def single_figure(start, end, df, half, bound_over="the four windows",
                  title=None):
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.36),
                           constrained_layout=True)
    draw_panel(ax, df, half)
    ax.set_title(title or f"{start}–{end}", loc="center")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel(Y_LABEL)
    caption(fig, caption_text([(start, end)], half, grid=False,
                              bound_over=bound_over))
    return _save(fig, f"{start}_{end}/lrr_{start}_{end}")


# Candidate windows over a reference (--candidates, Hannah 2026-10-02)

# The reference is grey and filled; the candidates are lines, so sign colours stay off them
REF_LINE_C = C["BASE"]
# Points between the frame and a title, clearing the fill bars and their labels
TITLE_PAD_FILLS = 22
# Every figure of the candidate study shares one y axis, taken over all of these
CANDIDATE_STUDY_WINDOWS = ((1996, 2015), (2010, 2026), (1996, 2026), (2010, 2020))


# The candidate study's shared y bound, and the windows it was taken over as text
def candidate_half():
    half = shared_bounds([load_window(*w) for w in CANDIDATE_STUDY_WINDOWS])
    over = ", ".join(f"{a}–{b}" for a, b in CANDIDATE_STUDY_WINDOWS)
    return half, over
REF_FILL_C = C["BASE_FILL"]
CANDIDATE_COLOURS = [INK, C["ACCENT"]]


# The last CoastSat date any transect in the window reaches
def record_end(start: int, end: int) -> str:
    t = pd.read_csv(domain_csv(start, end).parent / "transect_lrr_full.csv",
                    usecols=["end_date"])
    return str(t["end_date"].max())


# The candidates' caption
def candidates_caption(ref, cands, half, fills) -> str:
    (a, b) = ref
    lines = " and ".join(f"{s}–{e} ({n})" for (s, e), n
                         in zip(cands, ["black", "purple"]))
    # Name any end year the record covers less than half of
    ends = sorted({(e, d) for s, e in [ref] + cands
                   if (d := record_end(s, e)) < f"{e}-07-01"})
    end_txt = "".join(
        f" The record ends {d}, so a window to {e} holds only that much of {e}."
        for e, d in ends)
    fill_txt = "; ".join(f"{y} at GIS {lo}–{hi}" for y, lo, hi in fills)
    return ("Observed shoreline change rate by GIS domain (1 at Cape Point, "
            f"90 at Pea Island) for the candidate windows {lines}, over the "
            f"{a}–{b} rate (grey, filled). Each line is the mean linear "
            "regression rate of the CoastSat transects inside each 500 m "
            "domain (domain_lrr_summary.csv), fitted over the calendar window "
            "(1 January of the first year to 31 December of the last), "
            "seaward positive; the spread across transects is in the domain "
            "tables, not drawn." + end_txt
            + (f" Black bars above the frame mark the beach fills placed inside "
               f"{a}–{b} at the footprint the hindcast uses ({fill_txt}); the "
               "rates there include the placed sand." if fills else "")
            + " The hatched amber boxes mark the offshore shoals ("
            + "; ".join(f"{n} GIS {lo}–{hi}" for n, (lo, hi)
                        in HATTERAS_ANNOTATIONS.shoal_zones.items())
            + "). Village spans are shaded; the solid hairline is the Buxton "
            "groin and the dotted hairlines are the Avon and Rodanthe piers. "
            f"The y axis is ±{half:g} m/yr, the largest |mean| over "
            f"{candidate_half()[1]} plus 1 m rounded up, shared by every figure "
            "of the candidate windows.")


# One panel: the reference filled grey, each candidate a line on top
def candidates_figure(ref, cands, frames, half):
    df_ref, *df_cands = frames
    fills = fills_in(*ref)
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40),
                           constrained_layout=True)
    ax.set_xlim(0.5, N_DOMAINS + 0.5)
    ax.set_ylim(-half, half)
    town_bands(ax, label=True)
    ax.axhline(0, color=INK_MUTED, lw=0.6, zorder=2)
    x = df_ref["domain_number"].to_numpy(float)
    y = df_ref["mean_lrr"].to_numpy(float)
    ax.fill_between(x, 0, y, where=~np.isnan(y), interpolate=True,
                    color=REF_FILL_C, lw=0, zorder=3)
    (h_ref,) = ax.plot(x, y, color=REF_LINE_C, lw=0.9, zorder=4,
                       label=f"{ref[0]}–{ref[1]}")
    handles = [h_ref]
    for (s, e), df, col in zip(cands, df_cands, CANDIDATE_COLOURS):
        (ln,) = ax.plot(df["domain_number"], df["mean_lrr"], color=col,
                        lw=1.2, zorder=6, label=f"{s}–{e}")
        handles.append(ln)
    structures(ax, True, STRUCTURE_LABEL_PT)
    draw_shoals(ax, label=True)
    if fills:
        draw_fills(ax, fills, half)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(Y_TICK_M))
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    open_frame(ax)
    ax.set_title("Shoreline change rate, "
                 + " and ".join(f"{s}–{e}" for s, e in cands)
                 + f", against the {ref[0]}–{ref[1]} rate", loc="center",
                 pad=TITLE_PAD_FILLS if fills else None)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel(Y_LABEL)
    fig.legend(handles=handles, loc="outside lower center", ncol=len(handles),
               frameon=False)
    caption(fig, candidates_caption(ref, cands, half, fills))
    stem = "_and_".join(f"{s}_{e}" for s, e in cands)
    return _save(fig, f"{ref[0]}_{ref[1]}/lrr_{stem}_over_{ref[0]}_{ref[1]}")


# Each window on its own, then the candidates over the reference, all on one y axis
def run_candidates(cand_specs, ref_spec):
    cands = [tuple(int(v) for v in w.split("_")) for w in cand_specs]
    ref = tuple(int(v) for v in ref_spec.split("_"))
    frames = [load_window(*ref)] + [load_window(*w) for w in cands]
    half, over = candidate_half()
    written = []
    for (a, b), df in zip([ref] + cands, frames):
        written += single_figure(a, b, df, half, bound_over=over,
                                 title=f"Shoreline change rate, {a}–{b}")
    written += candidates_figure(ref, cands, frames, half)

    table = pd.DataFrame({"domain_number": np.arange(1, N_DOMAINS + 1)})
    for (a, b), df in zip([ref] + cands, frames):
        table[f"mean_{a}_{b}"] = df["mean_lrr"].round(3)
    stem = "_and_".join(f"{s}_{e}" for s, e in cands)
    table.to_csv(support_dir(OUT_DIR / f"{ref[0]}_{ref[1]}")
                 / f"lrr_{stem}_over_{ref[0]}_{ref[1]}.csv", index=False)
    print(f"y bounds  +/-{half:g} m/yr")
    for p in written:
        print("wrote    ", p.relative_to(_REPO))


# Windows linked end-to-start, each chain sorted by start, chains by their first start
def _chains(windows):
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


# 2 x 2 when the windows form two chains of two (a column per chain, the earlier window above)
def grid_figure(windows, frames, half):
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


# The long window and its halves (--overlay)

# TWO STACKED PANELS (Hannah, 2026-09-18)
HALF_COLOURS = [C["BASE"], INK]
FILL_BAR_C = INK


# (year, first_gis, last_gis) of every enabled model-input fill inside the window, from the config
def fills_in(start: int, end: int):
    return sorted((p.year, min(p.gis_domains), max(p.gis_domains))
                  for p in HATTERAS_NOURISHMENT_PROJECTS
                  if p.enabled and start <= p.year <= end)


# A bar just ABOVE the frame over each fill footprint, the year(s) on it; a refilled span is one bar
def draw_fills(ax, fills, half: float, label_pt: float = STRUCTURE_LABEL_PT):
    trans = ax.get_xaxis_transform()
    y = 1.025
    spans: dict = {}
    for year, lo, hi in fills:
        spans.setdefault((lo, hi), []).append(year)
    for (lo, hi), years in spans.items():
        ax.plot([lo - 0.45, hi + 0.45], [y, y], color=FILL_BAR_C, lw=2.2,
                solid_capstyle="butt", zorder=6, clip_on=False,
                transform=trans)
        label = f"{years[0]} fill" if len(years) == 1 else \
            f"{', '.join(str(v) for v in years)} fills"
        ax.text((lo + hi) / 2, y + 0.02, label, ha="center",
                va="bottom", fontsize=label_pt, color=FILL_BAR_C, zorder=6,
                clip_on=False, transform=trans)


# Shoals as hatched outline boxes, light enough not to compete (2026-09-18)
SHOAL_C = C["ADDED"]
SHOAL_HATCH = "///"
SHOAL_HATCH_ALPHA = 0.30
SHOAL_HATCH_LW = 0.5         # matplotlib 3.9: a global rc, read at draw
SHOAL_EDGE_ALPHA = 0.55
SHOAL_TEXT = "#8a620e"       # the label colour of the other shoal figures


# Each shoal zone as a hatched, outlined box the panel's full height behind the data, named if asked
def draw_shoals(ax, label: bool = True, label_pt: float = STRUCTURE_LABEL_PT):
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


# The chain of default windows that runs start -> end end-to-start
def halves_of(start: int, end: int):
    for chain in _chains(DEFAULT_WINDOWS):
        if chain[0][0] == start and chain[-1][1] == end:
            return chain
    raise SystemExit(f"no chain of default windows runs {start}-{end}; "
                     f"have {DEFAULT_WINDOWS}")


# The smallest whole metre per year that holds every line drawn, no pad
def tight_bound(frames) -> float:
    return float(math.ceil(max(float(np.nanmax(df["mean_lrr"].abs()))
                               for df in frames)))


# The retired overlay's caption
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


# Panel (b)
def draw_halves(ax, halves, frames, half):
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


# (a) the long window filled by sign
def overlay_figure(long_w, halves, frames, half, shared):
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


# Into the comparison folder and published to output/figures/2-observations/shoreline/, captioned
def _save_both(fig, window, stem):
    out = save(fig, OUT_DIR / window / stem)
    out += save(fig, figure_dir("observations", "shoreline") / stem)
    plt.close(fig)
    return out


# The retired halves overlay
def run_overlay(spec: str):
    a, b = (int(v) for v in spec.split("_"))
    halves = halves_of(a, b)
    frames = [load_window(a, b)] + [load_window(*w) for w in halves]
    # The single-window figures' bound, named in the caption
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


# Run: the stacked figure and each window's
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 2)[1])
    ap.add_argument("--windows", nargs="+", metavar="START_END",
                    help="rate windows, e.g. 1984_2004 (default: the four)")
    ap.add_argument("--overlay", metavar="START_END",
                    help="a long window drawn over the default windows that "
                         "chain across it, e.g. 1996_2024; draws only that")
    ap.add_argument("--candidates", nargs="+", metavar="START_END",
                    help="candidate windows, each drawn alone and then over "
                         "--reference, e.g. 1996_2015 2010_2026")
    ap.add_argument("--reference", metavar="START_END",
                    help="the long window the candidates are drawn over")
    args = ap.parse_args(argv)

    if args.candidates:
        if not args.reference:
            sys.exit("--candidates needs --reference")
        apply_style()
        run_candidates(args.candidates, args.reference)
        return

    if args.overlay:
        # --overlay retired 2026-09-19 (see lrr/chains/)
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

    # Since 2026-09-19 only the 2 x 2 is drawn, into 3-rates/coastsat/lrr/ (OUT_DIR)
    written = grid_figure(windows, frames, half)

    print(f"y bounds  +/-{half:g} m/yr")
    for p in written:
        print("wrote    ", p.relative_to(_REPO))


if __name__ == "__main__":
    main()
