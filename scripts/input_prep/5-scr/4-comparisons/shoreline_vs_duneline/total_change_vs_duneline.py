"""
total_change_vs_duneline.py
==============================================================================
TOTAL shoreline change against the dune line's MEASURED net change, per GIS
domain, over 1996-2010, 2010-2024 and the whole 1996-2024. Built 2026-09-21
(Hannah, by interview: "the difference between the measured dune line change
between 1996-2010 and 2010-2024, and compare that with the LRR net position
change over those same periods"); renamed from lrr_net_change.py and merged
with projected_vs_duneline.py the same day, after a second interview settled
the vocabulary.

WHY "TOTAL" AND NOT "PROJECTED".  A rate turned into a distance is named by
the window it was FITTED on, never by the arithmetic:

    TOTAL      the rate is evaluated over the SAME window it was fitted on.
               Every window here is total: LRR(1996-2010) x 14 yr,
               LRR(2010-2024) x 14 yr, LRR(1996-2024) x 28 yr. Nothing is
               extrapolated.
    PROJECTED  the rate is carried onto a window it was NOT fitted on. None
               of that happens here; it lives in
               3-rates/coastsat/projected/ and, against the model, in
               output/comparisons/target_comparison/projected/.

THE MERGE (2026-09-21, Hannah, by interview).  This folder absorbed
`shoreline_vs_duneline/projected/1996_2024`, which was the 1996-2024 rate x
the 25.72 yr DUNE-LINE INTERVAL rather than the 28 calendar years. That is
the same fit window with a shorter span, so it was never a projection either.
Its three numbers were already carried here as the `*_dune_interval_m`
columns -- checked identical to 0 m before the merge -- so the headline
figure is the 28 yr CALENDAR span and the dune interval is a column and a
caption line. Hannah chose the calendar span so every window in the tree is
read the same way. The old folder is in superseded_20260921/.

    Hannah, 2026-09-21: "Dont fit on the exact dune dates, the year is most
    important" -- so the fit window stays the stored CALENDAR window,
    1 January to 31 December.

WHAT IS COMPARED, per 500 m GIS domain, SEAWARD POSITIVE, in metres
    shoreline   3-rates/coastsat/lrr/<window>/transect_lrr_full.csv,
                lrr_m_yr x (end_year - start_year), mean over the ~10 CoastSat
                transects of the domain. 14 yr in each half, 28 yr over the
                whole -- the CALENDAR interval (Hannah's choice, 2026-09-21),
                matching the year labels and 3-rates/coastsat/total_change.
    dune line   3-rates/duneline/endpoint/<window>/ as stored, end line minus
                start line, mean over the ~5 dune transects. Read, never
                recomputed.
    beach width shoreline change - dune-line change; positive = the beach
                widened (the waterline gained on the dune).

THE INTERVAL MISMATCH, reported and not corrected.  The dune line is two
photographs, 11.63 yr apart in the first half (1997-10-12 to 2009-05-30),
14.09 yr in the second (to 2023-07-01, assumed) and 25.72 yr over the whole,
while the shoreline is carried over the full calendar span in each. So in the
first half the gap holds about 2.4 yr of shoreline drift on top of beach-width
change. Every table carries the interval-matched alternative beside the
headline value (`*_dune_interval_m`, the same rate x the dune line's own
span) and every PROVENANCE.md reports both, so the size of the term is
visible rather than argued about.

FIGURE TITLES carry quantity, window and method (Hannah, 2026-09-21), e.g.
"Total shoreline change vs dune line, 1996-2010 (LRR x 14 yr)".

WHAT IS DRAWN
    <window>/    two_panel, shaded_gap, overlay -- the three presentations,
                 one folder per window.
    chains/      1996-2024 stacked over its two halves on one y axis, the
                 shape of endpoint_net_change/chains.
    difference/  second half minus first half on both sides: where the
                 shoreline trend sped up or slowed, did the dune line follow?

OUTPUT   data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/total_change/
    <window>/domain_comparison.csv, the three PNGs, PROVENANCE.md,
             supporting/ (PDFs, CAPTIONS.md)
    chains/total_change_chain_1996_2010_2024.png,
             supporting/{domain_comparison,island_summary}.csv
    difference/total_change_difference.png, supporting/domain_difference.csv

USAGE
    python scripts/input_prep/5-scr/4-comparisons/shoreline_vs_duneline/total_change_vs_duneline.py
    python ... --start-year 1996 --end-year 2010
==============================================================================
"""

from __future__ import annotations

import argparse
import datetime as dt
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)

import rates_figures as rf  # noqa: E402  (the 3-rates drawing helpers)
from rates_figures import cw, plt  # noqa: E402
from duneline_vs_coastsat import beach_width_handles, shade_beach_width  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK, INK_MUTED, _title, apply_style, caption, figsize,
    compare_header, mark_offaxis, offaxis_clause, open_frame, save, structures,
    support_dir,
    town_bands,
)
from site_layer.hat_observed_rates import (  # noqa: E402
    COASTSAT_LRR_ROOT, DUNELINE_ENDPOINT_ROOT, ENDPOINT_DOMAIN_FILE,
    ENDPOINT_TRANSECT_FILE, PROJECTED_VS_DUNELINE_ENDPOINT,
    TOTAL_CHANGE_VS_DUNELINE,
)

N = cw.N_DOMAINS
# How the SHORELINE side is read. The dune side never varies -- it is the
# measured endpoint between the two digitized lines bounding the window -- so
# this is the only choice the script makes, and it is what the output folder
# is named for. Mirrors 3-rates/coastsat/{total_change,projected}.
PRODUCTS = {
    "total": {
        "root": TOTAL_CHANGE_VS_DUNELINE,
        "noun": "Total shoreline change",
        # each window's own rate, evaluated over that window: nothing extrapolated
        "rate_window": lambda s, e: (s, e),
        "windows": [(1996, 2010), (2010, 2024), (1996, 2024)],
        "stem": "coastsat_total_change_vs_duneline",
        # the halves AND the whole exist, so the stack is the 3-panel chain
        "stacked": "chain",
        "difference": True,
    },
    "projected": {
        "root": PROJECTED_VS_DUNELINE_ENDPOINT,
        "noun": "Projected shoreline change",
        "rate_window": lambda s, e: (1996, 2024),
        # no 1996_2024: there the rate window IS the change window, so the
        # answer is the TOTAL product, not a second copy under another name.
        "windows": [(1996, 2010), (2010, 2024)],
        "stem": "coastsat_projected_vs_duneline",
        # BOTH halves carry the same 1996-2024 rate over the same 14 yr, so
        # the shoreline side is IDENTICAL in them. That makes
        # change_between_periods zero by construction on that side -- not
        # drawn -- but it is exactly what makes the STACK worth drawing: one
        # prediction against two different dune-line outcomes, which is the
        # question this product exists to ask (Hannah, 2026-09-22).
        "stacked": "halves",
        "difference": False,
    },
}
PROD = PRODUCTS["total"]
WHOLE = (1996, 2024)
HALVES = [(1996, 2010), (2010, 2024)]
WINDOWS = HALVES + [WHOLE]          # the two the interview asked for, then context
CHAIN_WINDOWS = [WHOLE] + HALVES    # whole on top, as in net_change/chains
CHAIN_STEM = "coastsat_total_change_vs_duneline_1996_2010_2024_stacked"

C_SHORE = "#2166ac"     # the house shoreline blue of net_change/
C_DUNE = "#b2182b"      # and its dune red
C_GAP_TOWN = 0.055      # village bands as a strip: the gap already shades grey
C_DIFF_GAP = "0.88"     # the difference figure's neutral gap (not a beach width)
LW = 1.1
Y_LABEL = "Net change in position (m)"
# The shared fixed metre axis (Hannah, 2026-09-22); see
# coastsat_total_change.Y_HALF_M for why it is fixed and not a floor.
Y_HALF_M = 100.0
Y_TICK_M = 20.0


# --------------------------------------------------------------------------
# load
# --------------------------------------------------------------------------

def load(window):
    """Per-domain change from both sides, plus the per-transect shoreline.

    The shoreline is the stored calendar-window LRR carried over the CALENDAR
    interval; `*_dune_interval_m` is the same rate over the dune line's own
    span, kept beside it as the diagnostic for the mismatch.
    """
    s, e = window
    w = f"{s}_{e}"
    years = float(e - s)
    # `w` indexes the DUNE side and the change window; the shoreline rate may
    # come from a different window (see PRODUCTS).

    dune_t = pd.read_csv(DUNELINE_ENDPOINT_ROOT / w / ENDPOINT_TRANSECT_FILE)
    dune_d = pd.read_csv(DUNELINE_ENDPOINT_ROOT / w / ENDPOINT_DOMAIN_FILE)
    meta = dune_t.iloc[0]
    dune_years = float(meta["interval_yr"])

    rs, re_ = PROD["rate_window"](s, e)
    cs = pd.read_csv(COASTSAT_LRR_ROOT / f"{rs}_{re_}" / "transect_lrr_full.csv")
    cs = cs[cs["domain_number"].between(1, N)].copy()
    cs["domain_number"] = cs["domain_number"].astype(int)
    cs["shoreline_change_m"] = cs["lrr_m_yr"] * years
    cs["shoreline_change_dune_interval_m"] = cs["lrr_m_yr"] * dune_years

    g = cs.groupby("domain_number")
    dom = pd.DataFrame({
        "n_coastsat_transects": g.size(),
        "mean_lrr_m_yr": g["lrr_m_yr"].mean(),
        "shoreline_change_m": g["shoreline_change_m"].mean(),
        "shoreline_change_dune_interval_m": g["shoreline_change_dune_interval_m"].mean(),
    }).reindex(range(1, N + 1))
    d = dune_d.set_index(dune_d["domain_number"].astype(int))
    dom["n_dune_transects"] = d["n_transects"].reindex(dom.index)
    dom["dune_change_m"] = d["mean_change_m"].reindex(dom.index)
    dom["dune_rate_m_yr"] = dom["dune_change_m"] / dune_years
    dom["beach_width_change_m"] = dom["shoreline_change_m"] - dom["dune_change_m"]
    dom["beach_width_change_dune_interval_m"] = (
        dom["shoreline_change_dune_interval_m"] - dom["dune_change_m"])
    dom = dom.rename_axis("domain_number").reset_index()
    dom.insert(1, "lrr_interval_yr", round(years, 4))
    dom.insert(2, "dune_interval_yr", round(dune_years, 4))
    return dom.round(3), cs, meta, years, dune_years


def stats(dom, years, dune_years):
    ok = dom.dropna(subset=["shoreline_change_m", "dune_change_m"])
    x, y = ok["shoreline_change_m"], ok["dune_change_m"]
    bw, bwm = ok["beach_width_change_m"], ok["beach_width_change_dune_interval_m"]
    slope, icpt = np.polyfit(x, y, 1)
    return {
        "n_domains": len(ok),
        "lrr_interval_yr": round(years, 4),
        "dune_interval_yr": round(dune_years, 4),
        "mean_shoreline_change_m": round(x.mean(), 2),
        "mean_dune_change_m": round(y.mean(), 2),
        "mean_beach_width_change_m": round(bw.mean(), 2),
        "mean_beach_width_dune_interval_m": round(bwm.mean(), 2),
        "beach_width_min_m": round(bw.min(), 2),
        "beach_width_max_m": round(bw.max(), 2),
        "domains_beach_narrowed": int((bw < 0).sum()),
        "r_dune_vs_shoreline": round(float(np.corrcoef(x, y)[0, 1]), 3),
        "slope_dune_on_shoreline": round(float(slope), 3),
        "intercept_m": round(float(icpt), 2),
        "rmse_dune_minus_shoreline_m": round(float(np.sqrt(((y - x) ** 2).mean())), 2),
        "domains_same_sign": int((np.sign(x) == np.sign(y)).sum()),
        "domains_shoreline_landward": int((x < 0).sum()),
        "domains_dune_landward": int((y < 0).sum()),
    }


# --------------------------------------------------------------------------
# per-window figures (three presentations, from the absorbed projected_vs_duneline)
# --------------------------------------------------------------------------

def _empty(dom):
    return rf._frame(dom.assign(_nan=np.nan), "_nan")


def _marks(ax, half, label, window):
    cw.draw_shoals(ax, label=label)
    fills = cw.fills_in(*window)
    if fills and label:
        cw.draw_fills(ax, fills, half)


def _shoreline_panel(ax, dom, cs, half, window, label=True):
    cw.draw_panel(ax, rf._frame(dom, "shoreline_change_m"), half, label=label, std=False)
    tt, x = rf._along(cs)
    n_out = rf._dots(ax, x, tt["shoreline_change_m"].to_numpy(float), half)
    ax.plot(dom["domain_number"], dom["dune_change_m"], color=INK, lw=LW,
            zorder=12, solid_capstyle="round")
    _marks(ax, half, label, window)
    # The axis is fixed at +/-100 m, so a domain-mean line can leave it. Mark
    # it at the edge and hand the values back for the caption; nothing that
    # walks off the panel should do so silently.
    off = [("the shoreline change",
            mark_offaxis(ax, dom["domain_number"], dom["shoreline_change_m"],
                         half, color=C_SHORE)),
           ("the dune-line change",
            mark_offaxis(ax, dom["domain_number"], dom["dune_change_m"],
                         half, color=INK))]
    return n_out, off


def _gap_panel(ax, dom, half, window, label=True):
    cw.draw_panel(ax, _empty(dom), half, label=label, std=False)
    x = dom["domain_number"].to_numpy(float)
    ys = dom["shoreline_change_m"].to_numpy(float)
    yd = dom["dune_change_m"].to_numpy(float)
    shade_beach_width(ax, x, ys, yd)
    ax.plot(x, yd, color=C_DUNE, lw=LW, zorder=12)
    ax.plot(x, ys, color=C_SHORE, lw=LW, zorder=12)
    _marks(ax, half, label, window)
    return [("the shoreline change", mark_offaxis(ax, x, ys, half, color=C_SHORE)),
            ("the dune-line change", mark_offaxis(ax, x, yd, half, color=C_DUNE))]


def _width_panel(ax, dom, half, window, label=False):
    cw.draw_panel(ax, _empty(dom), half, label=label, std=False)
    x = dom["domain_number"].to_numpy(float)
    y = dom["beach_width_change_m"].to_numpy(float)
    shade_beach_width(ax, x, y, np.zeros_like(y))
    ax.plot(x, y, color=INK, lw=LW, zorder=12)
    _marks(ax, half, label, window)
    return [("the beach-width change", mark_offaxis(ax, x, y, half, color=INK))]


def _house_handles():
    return [(Line2D([], [], color=cw.C_ACCRETE, lw=1.0), Line2D([], [], color=cw.C_ERODE, lw=1.0)),
            (Line2D([], [], color=cw.C_ACCRETE, marker="o", ms=2.2, lw=0),
             Line2D([], [], color=cw.C_ERODE, marker="o", ms=2.2, lw=0)),
            Line2D([], [], color=INK, lw=LW)]


def _legend(fig, handles, labels, ncol):
    fig.legend(handles, labels, loc="outside lower center", ncol=ncol, frameon=False,
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})


def _dune_side(meta, dune_years):
    """"dune line: 1997-10-12 -> 2023-07-01 (25.7 yr), measured" for the header."""
    return (f"dune line: {meta['start_date']} → {meta['end_date']}"
            + (" (assumed)" if bool(meta["end_date_assumed"]) else "")
            + f",  {dune_years:.1f} yr,  measured")


def _header(window, meta, years, dune_years):
    """The two sides named on the canvas (Hannah, 2026-09-22).

    The interval mismatch is the reason this is worth the space: the shoreline
    is carried over the full CALENDAR span and the dune line spans whatever
    its two photographs do, so the gap between them holds that difference as
    well as beach-width change. Stated here, it cannot be missed.
    """
    s, e = window
    rs, re_ = PROD["rate_window"](s, e)
    note = ("" if (rs, re_) == (s, e) else
            "  (fitted on the FULL period, carried onto this half)")
    return [f"{s}–{e}   ·   shoreline: CoastSat LRR {rs}–{re_} × {years:.0f} yr{note}",
            _dune_side(meta, dune_years)]


def _base_caption(window, meta, years, dune_years, st):
    s, e = window
    return (
        f"The shoreline change is the CoastSat trend's: for each transect the "
        f"{s}–{e} linear regression rate (the ordinary-least-squares slope through "
        f"every satellite position from 1 January {s} to 31 December {e}) "
        f"multiplied by {years:.0f} yr, the calendar interval, averaged over the "
        "~10 transects in each 500 m domain. The rate is fitted and evaluated "
        "inside the same window, so nothing is extrapolated. The dune-line change "
        f"is measured: the {int(meta['end_vintage'])} digitized dune line minus the "
        f"{int(meta['start_vintage'])} line along the 100 m transects, averaged over "
        f"the ~5 in each domain, spanning {meta['start_date']} to {meta['end_date']}"
        + (" (assumed)" if bool(meta["end_date_assumed"]) else "")
        + f", {dune_years:.2f} yr. Both are seaward positive, in metres. Beach-width "
        "change is shoreline change minus dune-line change: positive where the beach "
        "widened, negative where it narrowed"
        + (f"; because the two spans differ by {years - dune_years:+.1f} yr, it also "
           f"holds that much shoreline drift (over the dune line's own "
           f"{dune_years:.2f} yr the mean would be "
           f"{st['mean_beach_width_dune_interval_m']:+.1f} m rather than "
           f"{st['mean_beach_width_change_m']:+.1f} m)"
           if abs(years - dune_years) > 0.25 else "")
        + ". " + rf._marks_clause(s, e))


def _stats_clause(st, half):
    return (f" Over the {st['n_domains']} domains the shoreline changed "
            f"{st['mean_shoreline_change_m']:+.1f} m on average and the dune line "
            f"{st['mean_dune_change_m']:+.1f} m; the beach widened by "
            f"{st['mean_beach_width_change_m']:+.1f} m (range "
            f"{st['beach_width_min_m']:+.1f} to {st['beach_width_max_m']:+.1f} m) and "
            f"narrowed in {st['domains_beach_narrowed']} of them. r(shoreline, dune "
            f"line) = {st['r_dune_vs_shoreline']:.2f}. The y axis is ±{half:g} m on "
            "every panel.")


def _pad_title(ax, window):
    """Lift the centred title clear of the fill bars.

    draw_fills puts its bars at 1.025 in axes fractions with the year label
    above them, so a title at the default pad lands on top of "2022 fill".
    _title() has already set the bold letter at the left; re-setting only the
    centred string keeps it and moves both (pad is per-axes in matplotlib).
    """
    if cw.fills_in(*window):
        ax.set_title(ax.get_title(loc="center"), loc="center", pad=20)


def halves_overlay_figure(loaded, out_dir):
    """The two halves' overlay panels on one sheet, 1996-2010 above 2010-2024.

    TITLE SPACE (Hannah, 2026-09-22: "be strategic ... concise yet organized").
    The method is the SAME in both panels -- that is the whole point of the
    projected product -- so it is stated ONCE in the header and never repeated.
    Each panel title then carries only what actually differs between them: the
    window, and the two dune-line dates with their real interval. Nothing is
    said twice, and nothing a reader needs is only in the caption.
    """
    windows = [w for w in PROD["windows"] if w[1] - w[0] == 14]
    half, tick = Y_HALF_M, Y_TICK_M
    fig, axes = plt.subplots(len(windows), 1, sharex=True, sharey=True,
                             constrained_layout=True,
                             figsize=figsize("double", height=5.6))
    stats, off, n_out_total = [], [], 0
    for i, (ax, w) in enumerate(zip(axes, windows)):
        dom, cs, meta, years, dune_years = loaded[w]
        n_out, o = _shoreline_panel(ax, dom, cs, half, w, label=(i == 0))
        n_out_total += n_out
        off += [(f"{w[0]}–{w[1]} {lab}", pts) for lab, pts in o]
        ax.yaxis.set_major_locator(MultipleLocator(tick))
        # Only what differs: the window and the dune line's own dates.
        if i and cw.fills_in(*w):
            cw.draw_fills(ax, cw.fills_in(*w), half)
        _title(ax, i, f"{w[0]}–{w[1]}   ·   dune line {meta['start_date']} → "
                      f"{meta['end_date']}"
                      + (" (assumed)" if bool(meta["end_date_assumed"]) else "")
                      + f",  {dune_years:.1f} yr")
        _pad_title(ax, w)
        stats.append((w, dom, dune_years))
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    fig.supylabel(Y_LABEL, fontsize=9)

    shared_rate = len({PROD["rate_window"](*w) for w in windows}) == 1
    rs, re_ = PROD["rate_window"](*windows[0])
    _legend(fig, _house_handles(),
            [(f"{PROD['noun']} (CoastSat LRR {rs}–{re_} × 14 yr)" if shared_rate
              else f"{PROD['noun']} (own CoastSat LRR × 14 yr)"),
             "Individual CoastSat transects",
             "Total dune line change (measured)"], ncol=3)
    # Said once above the panels, because it is what the two panels have in
    # common -- and for `projected` it is also the whole point of the figure.
    compare_header(fig, f"{PROD['noun']} vs dune line   ·   shoreline is "
                        + (f"the SAME CoastSat LRR {rs}–{re_} × 14 yr in both panels"
                           if shared_rate else
                           "each panel's OWN CoastSat LRR × its own 14 yr"))

    per = "; ".join(
        f"{w[0]}–{w[1]} shoreline {d['shoreline_change_m'].mean():+.1f} m, dune line "
        f"{d['dune_change_m'].mean():+.1f} m, beach width "
        f"{d['beach_width_change_m'].mean():+.1f} m, r = "
        f"{d[['shoreline_change_m', 'dune_change_m']].corr().iloc[0, 1]:.2f}"
        for w, d, _ in stats)
    caption(fig, (
        (f"The {rs}–{re_} CoastSat linear regression rate × 14 yr — one long-term "
         "trend, PROJECTED onto each half — " if shared_rate else
         "Each half's OWN CoastSat linear regression rate × its own 14 yr — "
         "fitted and evaluated inside the same window, so nothing is "
         "extrapolated — ")
        + "against the dune line's measured net "
        "change over that same half, by GIS domain (1 at Cape Point, 90 at Pea "
        "Island). The coloured line and fill are the shoreline (blue seaward, red "
        "landward) with its single transects as dots"
        + (f" ({n_out_total} beyond the axis in total, open circles at its edge)"
           if n_out_total else "")
        + "; the black line is the dune line. "
        + ("**The shoreline side is identical in the two panels** — same rate, same "
           "14 yr — so every difference between them is the dune line's, which is "
           "what this figure is for. " if shared_rate else
           "Both sides change between the panels here: the shoreline because each "
           "half has its own fitted rate, the dune line because each half has its "
           "own pair of photographs. Its companion, "
           "`coastsat_projected_vs_duneline_endpoint/all_windows_stacked/`, holds "
           "the shoreline side FIXED at the 1996–2024 rate, so the two sheets "
           "together separate what the trend predicts from what it was fitted on. ")
        + "Where the black "
        "line sits below the fill's edge the beach widened. The dune line spans its "
        "own two photographs, 11.6 yr in the first half and 14.1 yr in the second "
        "(the 2023 date assumed), against the shoreline's 14 calendar years in both, "
        f"so each gap also holds that much shoreline drift. Per window: {per}. "
        "Seaward positive, in metres; the y axis is ±"
        f"{half:g} m on both panels, the fixed metre axis. "
        + _marks_note(windows) + offaxis_clause(off, half)))
    written = save(fig, out_dir / f"{PROD['stem']}_{windows[0][0]}_{windows[0][1]}_"
                                  f"{windows[-1][1]}_halves_overlay")
    plt.close(fig)
    return written


def _marks_note(windows):
    """The shoal/village/fill sentence, said once for the stacked figure."""
    return rf._marks_clause(*windows[0]).replace(
        "Black bars above the panel mark the beach fills placed inside the window",
        "Black bars above the TOP panel mark the beach fills placed inside it")


def window_figures(window, dom, cs, meta, years, dune_years, st, out_dir):
    half, tick = Y_HALF_M, Y_TICK_M
    s, e = window
    base = _base_caption(window, meta, years, dune_years, st) + _stats_clause(st, half)
    # The fit window is IN the label, not implied by the folder (Hannah,
    # 2026-09-21): "CoastSat LRR 1996–2010 × 14 yr" says where the rate came
    # from, and the reader can see it matches the window in the title -- which
    # is the whole difference between total change and a projection.
    rs, re_ = PROD["rate_window"](s, e)
    lab_shore = (f"{PROD['noun']} (CoastSat LRR {rs}–{re_} × {years:.0f} yr)")
    lab_dune = "Total dune line change (measured)"
    written = []

    # 1. two panels
    fig, axes = plt.subplots(2, 1, sharex=True, constrained_layout=True,
                             figsize=figsize("double", height=5.6))
    n_out, off = _shoreline_panel(axes[0], dom, cs, half, window, label=True)
    off = off + _width_panel(axes[1], dom, half, window)
    for ax in axes:
        ax.yaxis.set_major_locator(MultipleLocator(tick))
    axes[0].set_ylabel(Y_LABEL)
    axes[1].set_ylabel("Beach-width change (m)")
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    # Quantity, window, method (Hannah, 2026-09-21). "Total" because the rate
    # is fitted on the window it is evaluated over; see the module docstring.
    _title(axes[0], 0, f"{PROD['noun']} vs dune line, {s}–{e} "
                       f"(CoastSat LRR {rs}–{re_} × {years:.0f} yr)")
    _title(axes[1], 1, f"Beach-width change, {s}–{e}")
    _pad_title(axes[0], window)
    _legend(fig, _house_handles() + beach_width_handles(),
            [lab_shore, "Individual CoastSat transects", lab_dune,
             "Beach widened", "Beach narrowed"], ncol=3)
    compare_header(fig, _header(window, meta, years, dune_years))
    caption(fig, (
        "Trend-implied shoreline change against measured dune-line change by GIS "
        f"domain (1 at Cape Point, 90 at Pea Island), {s}–{e}. (a) The shoreline as "
        "the coloured line and fill (blue seaward, red landward) with its single "
        "transects as dots"
        + (f" ({n_out} beyond the axis, open circles at its edge)" if n_out else "")
        + "; the dune line as the black line. (b) Beach-width change, the shoreline "
        "change minus the dune-line change: solid grey where the beach widened, "
        "hatched where it narrowed. " + base + offaxis_clause(off, half)))
    written += save(fig, out_dir / f"{PROD['stem']}_{s}_{e}_two_panel")
    plt.close(fig)

    # 2. shaded gap
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40), constrained_layout=True)
    off = _gap_panel(ax, dom, half, window)
    ax.yaxis.set_major_locator(MultipleLocator(tick))
    ax.set_ylabel(Y_LABEL)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    _legend(fig, [Line2D([], [], color=C_SHORE, lw=LW), Line2D([], [], color=C_DUNE, lw=LW)]
            + beach_width_handles(),
            [lab_shore, lab_dune, "Beach widened", "Beach narrowed"], ncol=2)
    compare_header(fig, _header(window, meta, years, dune_years))
    caption(fig, (
        "Trend-implied shoreline change (blue) and measured dune-line change (red) by "
        f"GIS domain (1 at Cape Point, 90 at Pea Island), {s}–{e}, domain means. The "
        "space between them is the beach-width change: solid grey where the beach "
        "widened (the blue line above the red), hatched where it narrowed. " + base + offaxis_clause(off, half)))
    written += save(fig, out_dir / f"{PROD['stem']}_{s}_{e}_shaded_gap")
    plt.close(fig)

    # 3. overlay
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40), constrained_layout=True)
    n_out, off = _shoreline_panel(ax, dom, cs, half, window, label=True)
    ax.yaxis.set_major_locator(MultipleLocator(tick))
    ax.set_ylabel(Y_LABEL)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    _legend(fig, _house_handles(),
            [lab_shore, "Individual CoastSat transects", lab_dune], ncol=2)
    compare_header(fig, _header(window, meta, years, dune_years))
    caption(fig, (
        "Trend-implied shoreline change against measured dune-line change by GIS "
        f"domain (1 at Cape Point, 90 at Pea Island), {s}–{e}: the shoreline as the "
        "coloured line and fill (blue seaward, red landward) with its single "
        "transects as dots"
        + (f" ({n_out} beyond the axis, open circles at its edge)" if n_out else "")
        + "; the dune line as the black line. Where the black line is below the "
        "fill's edge, the beach widened. " + base + offaxis_clause(off, half)))
    written += save(fig, out_dir / f"{PROD['stem']}_{s}_{e}_overlay")
    plt.close(fig)
    return half, written


def write_provenance(window, out_dir, meta, years, dune_years, st, half):
    s, e = window
    drift = years - dune_years
    (out_dir / "PROVENANCE.md").write_text("\n".join([
        f"# 4-comparisons/shoreline_vs_duneline/{PROD['root'].name}/{s}_{e}"
        " - provenance",
        "",
        f"Written {dt.datetime.now():%Y-%m-%d %H:%M} by "
        "scripts/input_prep/5-scr/4-comparisons/shoreline_vs_duneline/"
        "total_change_vs_duneline.py.",
        "",
        f"- Shoreline: `3-rates/coastsat/lrr/{s}_{e}/transect_lrr_full.csv`, lrr_m_yr "
        f"x {years:.0f} yr (the CALENDAR interval; the fit runs 1 January {s} to 31 "
        f"December {e}). Fitted and evaluated in the same window, so nothing is "
        "extrapolated - this is the fitted trend's net change, not a projection.",
        f"- Dune line: `3-rates/duneline/endpoint/{s}_{e}/` as stored "
        f"({int(meta['start_vintage'])} and {int(meta['end_vintage'])} lines), "
        f"{meta['start_date']} to {meta['end_date']}"
        + (" (ASSUMED)" if bool(meta["end_date_assumed"]) else "")
        + f", {dune_years:.2f} yr.",
        "- Beach-width change = shoreline change - dune-line change; positive = widened.",
        "- Domain means only: the two use different transects "
        "(~10 CoastSat, ~5 dune per 500 m domain).",
        "",
        "## The interval mismatch",
        "",
        f"The shoreline is carried over {years:.0f} calendar years, the dune line "
        f"measured over {dune_years:.2f}, a difference of {drift:+.2f} yr (Hannah "
        "chose the calendar interval, 2026-09-21: the year is what the period means). "
        "So the beach-width gap holds that much shoreline drift on top of true width "
        "change. Nothing is corrected for it; the interval-matched value is carried "
        "beside the headline one in every table "
        "(`*_dune_interval_m` = the same rate x the dune line's own span).",
        "",
        f"Mean beach-width change: **{st['mean_beach_width_change_m']:+.1f} m** at "
        f"{years:.0f} yr, {st['mean_beach_width_dune_interval_m']:+.1f} m at "
        f"{dune_years:.2f} yr - a {abs(st['mean_beach_width_change_m'] - st['mean_beach_width_dune_interval_m']):.1f} m term.",
        "",
        "## Island summary",
        "",
        f"Domain mean shoreline {st['mean_shoreline_change_m']:+.1f} m, dune line "
        f"{st['mean_dune_change_m']:+.1f} m, beach width "
        f"{st['mean_beach_width_change_m']:+.1f} m (range {st['beach_width_min_m']:+.1f} "
        f"to {st['beach_width_max_m']:+.1f}); beach narrowed in "
        f"{st['domains_beach_narrowed']} of {st['n_domains']} domains. "
        f"r(shoreline, dune line) = {st['r_dune_vs_shoreline']:.2f}, slope "
        f"{st['slope_dune_on_shoreline']:.2f}, RMSE "
        f"{st['rmse_dune_minus_shoreline_m']:.1f} m; they agree in sign in "
        f"{st['domains_same_sign']} of {st['n_domains']}. The shoreline is landward in "
        f"{st['domains_shoreline_landward']}, the dune line in "
        f"{st['domains_dune_landward']}. y axis ±{half:g} m.",
        "",
    ]), encoding="utf-8")


# --------------------------------------------------------------------------
# the chain: the whole over its two halves
# --------------------------------------------------------------------------

def _chain_panel(ax, dom, half, label):
    ax.set_xlim(0.5, N + 0.5)
    ax.set_ylim(-half, half)
    town_bands(ax, label=label, strip=C_GAP_TOWN)
    cw.draw_shoals(ax, label=False)
    ax.axhline(0, color=INK_MUTED, lw=0.6, zorder=2)
    x = dom["domain_number"].to_numpy(float)
    ys = dom["shoreline_change_m"].to_numpy(float)
    yd = dom["dune_change_m"].to_numpy(float)
    shade_beach_width(ax, x, ys, yd)
    ax.plot(x, yd, color=C_DUNE, lw=LW, zorder=5)
    ax.plot(x, ys, color=C_SHORE, lw=LW, zorder=5)
    structures(ax, label)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(20.0))
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    open_frame(ax)


def chain_figure(loaded, summaries, out):
    ext = max(float(np.nanmax(np.abs(loaded[w][0][c])))
              for w in CHAIN_WINDOWS
              for c in ("shoreline_change_m", "dune_change_m"))
    half = float(math.ceil(ext / rf.Y_STEP_M) * rf.Y_STEP_M)
    by = {f"{s}_{e}": st for (s, e), st in zip(WINDOWS, summaries)}
    whole = by[f"{WHOLE[0]}_{WHOLE[1]}"]
    h1, h2 = (by[f"{s}_{e}"] for s, e in HALVES)
    # the shoreline halves need NOT add up to the whole: each is its own fit
    add = h1["mean_shoreline_change_m"] + h2["mean_shoreline_change_m"]
    dune_add = h1["mean_dune_change_m"] + h2["mean_dune_change_m"]

    apply_style()
    fig, axes = plt.subplots(3, 1, sharex=True, sharey=True, constrained_layout=True,
                             figsize=figsize("double", height=6.6))
    for i, (ax, w) in enumerate(zip(axes, CHAIN_WINDOWS)):
        dom, _, meta, years, _ = loaded[w]
        _chain_panel(ax, dom, half, label=(i == 0))
        _title(ax, i, f"Total shoreline change, {w[0]}–{w[1]} "
                      f"(CoastSat LRR {w[0]}–{w[1]} × {years:.0f} yr)")
        if i == 0:
            _pad_title(ax, WHOLE)   # only axes[0] carries the fill bars here
    fills = cw.fills_in(*WHOLE)
    cw.draw_fills(axes[0], fills, half)
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    fig.supylabel(Y_LABEL, fontsize=9)
    compare_header(fig, [
        "shoreline: each panel's own CoastSat LRR × that window's calendar span",
        "dune line: measured between the two digitized lines bounding it  ·  "
        + ";  ".join(f"{w[0]}–{w[1]} {loaded[w][4]:.1f} yr" for w in CHAIN_WINDOWS)])
    fig.legend(handles=[Line2D([], [], color=C_SHORE, lw=1.2,
                               label="Total shoreline change (each panel's own "
                                     "CoastSat LRR × that window)"),
                        Line2D([], [], color=C_DUNE, lw=1.2,
                               label="Total dune line change (measured)")]
               + beach_width_handles(),
               loc="outside lower center", ncol=4, frameon=False)
    caption(fig, (
        "Trend-implied shoreline change and measured dune-line change by GIS domain "
        "(1 at Cape Point, 90 at Pea Island), seaward positive, in metres: "
        "(a) 1996–2024, (b) 1996–2010, (c) 2010–2024. Blue: the CoastSat linear "
        "regression rate fitted to every satellite position inside that calendar "
        "window, multiplied by the window's length (28, 14 and 14 yr), averaged over "
        "the ~10 transects of each 500 m domain; each rate is fitted and evaluated "
        "over the same years, so none is extrapolated. Red: the dune line, the end "
        "digitized line minus the start line along the 100 m transects, averaged over "
        "the ~5 of each domain (1997-10-12, 2009-05-30 and 2023-07-01, the last "
        "assumed). The gap between them is the change in beach width: solid grey "
        "where the beach widened (blue above red), hatched where it narrowed. The "
        "dune line's two halves add up to its whole exactly "
        f"({h1['mean_dune_change_m']:+.1f} and {h2['mean_dune_change_m']:+.1f} give "
        f"{dune_add:+.1f} m against {whole['mean_dune_change_m']:+.1f}); the "
        "shoreline's need not, because each half is its own fit — "
        f"{h1['mean_shoreline_change_m']:+.1f} and "
        f"{h2['mean_shoreline_change_m']:+.1f} give {add:+.1f} m against the "
        f"{whole['mean_shoreline_change_m']:+.1f} m of the single 28-yr trend, and "
        "the difference is the trend break the halves resolve. The dune line is "
        "measured over 11.63 and 14.09 yr against the shoreline's 14, so each gap "
        "also holds that much drift. Black bars above (a) mark the beach fills inside "
        "the window at the footprint the hindcast uses; hatched amber boxes mark the "
        "offshore shoals. Village spans are the grey strip along the top of each "
        "panel; the solid hairline is the Buxton groin and the dotted hairlines the "
        "Avon and Rodanthe piers. All panels share a y axis of "
        f"±{half:g} m."))
    written = save(fig, out / CHAIN_STEM)
    plt.close(fig)
    return half, written


# --------------------------------------------------------------------------
# the difference: second half minus first half, both sides
# --------------------------------------------------------------------------

def difference_table(loaded):
    a = loaded[HALVES[0]][0].set_index("domain_number")
    b = loaded[HALVES[1]][0].set_index("domain_number")
    d = pd.DataFrame(index=a.index)
    d["shoreline_first_m"] = a["shoreline_change_m"]
    d["shoreline_second_m"] = b["shoreline_change_m"]
    d["d_shoreline_m"] = b["shoreline_change_m"] - a["shoreline_change_m"]
    d["dune_first_m"] = a["dune_change_m"]
    d["dune_second_m"] = b["dune_change_m"]
    d["d_dune_m"] = b["dune_change_m"] - a["dune_change_m"]
    # rate form: the interval mismatch cancels out of a difference of rates
    d["d_shoreline_m_yr"] = b["mean_lrr_m_yr"] - a["mean_lrr_m_yr"]
    d["d_dune_m_yr"] = b["dune_rate_m_yr"] - a["dune_rate_m_yr"]
    d["same_sign"] = np.sign(d["d_shoreline_m"]) == np.sign(d["d_dune_m"])
    return d.round(3)


def difference_figure(d, out):
    ok = d.dropna(subset=["d_shoreline_m", "d_dune_m"])
    r_m = float(np.corrcoef(ok["d_shoreline_m"], ok["d_dune_m"])[0, 1])
    ok_r = d.dropna(subset=["d_shoreline_m_yr", "d_dune_m_yr"])
    r_rate = float(np.corrcoef(ok_r["d_shoreline_m_yr"], ok_r["d_dune_m_yr"])[0, 1])
    ext = float(np.nanmax(np.abs(d[["d_shoreline_m", "d_dune_m"]].to_numpy(float))))
    half = float(math.ceil((ext + 5) / rf.Y_STEP_M) * rf.Y_STEP_M)

    apply_style()
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40), constrained_layout=True)
    ax.set_xlim(0.5, N + 0.5)
    ax.set_ylim(-half, half)
    town_bands(ax, label=True, strip=C_GAP_TOWN)
    cw.draw_shoals(ax, label=True)
    ax.axhline(0, color=INK_MUTED, lw=0.8, zorder=2)
    x = d.index.to_numpy(float)
    ys = d["d_shoreline_m"].to_numpy(float)
    yd = d["d_dune_m"].to_numpy(float)
    # NOT shade_beach_width: both lines are differences, so the gap between
    # them is disagreement, not a width. One neutral fill, no hatch, so the
    # figure does not borrow a vocabulary that would read as beach width.
    ax.fill_between(x, ys, yd, color=C_DIFF_GAP, lw=0, zorder=3)
    ax.plot(x, yd, color=C_DUNE, lw=LW, zorder=5)
    ax.plot(x, ys, color=C_SHORE, lw=LW, zorder=5)
    structures(ax, True)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(20.0 if half > 60 else 10.0))
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    open_frame(ax)
    ax.set_ylabel("Second half − first half (m)")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    # Each half is fitted on itself, so there is no single fit window to name
    # in the title; the legend carries it per line.
    ax.set_title("Total shoreline change and dune line, second half − first half "
                 "(each half's own CoastSat LRR × 14 yr)")
    compare_header(fig, [
        "shoreline: CoastSat LRR 2010–2024 × 14 yr  −  CoastSat LRR 1996–2010 × 14 yr",
        "dune line: measured 2009–2023  −  measured 1997–2009"])
    fig.legend(handles=[Line2D([], [], color=C_SHORE, lw=1.2,
                               label="Shoreline (CoastSat LRR 2010–2024 × 14 yr "
                                     "− CoastSat LRR 1996–2010 × 14 yr)"),
                        Line2D([], [], color=C_DUNE, lw=1.2,
                               label="Total dune line change (measured), "
                                     "2009–2023 − 1997–2009"),
                        Patch(facecolor=C_DIFF_GAP, edgecolor="none",
                              label="Gap between them")],
               loc="outside lower center", ncol=2, frameon=False)
    caption(fig, (
        "How much each side changed between the two periods, by GIS domain (1 at Cape "
        "Point, 90 at Pea Island): the second half's net change minus the first "
        "half's, seaward positive. Above zero the feature did better in 2010–2024 "
        "than in 1996–2010 (it retreated less, or advanced more); below zero it did "
        "worse. Blue: the CoastSat trend, each half's linear regression rate times "
        "its 14 calendar years. Red: the dune line, each half measured between its "
        "own two digitized lines. The grey fill is the gap between them, where the "
        "two features changed pace by different amounts; it is not a beach width, "
        "since both lines are differences, so it carries no widened/narrowed hatching. "
        f"Over {len(ok)} domains r = {r_m:.2f} in metres and "
        f"{r_rate:.2f} as rates (m/yr, which removes the dune line's unequal spans of "
        "11.63 and 14.09 yr); the two moved the same way in "
        f"{int(ok['same_sign'].sum())} of {len(ok)} domains. Island means: shoreline "
        f"{ok['d_shoreline_m'].mean():+.1f} m, dune line {ok['d_dune_m'].mean():+.1f} "
        "m. Village spans are the grey strip along the top; the solid hairline is the "
        "Buxton groin and the dotted hairlines the Avon and Rodanthe piers; hatched "
        f"amber boxes mark the offshore shoals. y axis ±{half:g} m."))
    written = save(fig, out / "coastsat_total_change_vs_duneline_change_between_periods")
    plt.close(fig)
    return r_m, r_rate, half, written


# --------------------------------------------------------------------------

def main() -> int:
    global PROD
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--product", choices=sorted(PRODUCTS), default="total",
                    help="how the SHORELINE side is read. total (default): each "
                         "window's own LRR over its own years. projected: the "
                         "1996-2024 LRR carried onto each 14-yr half. The dune "
                         "side is the measured endpoint either way.")
    ap.add_argument("--start-year", type=int)
    ap.add_argument("--end-year", type=int)
    args = ap.parse_args()
    PROD = PRODUCTS[args.product]
    one = (args.start_year, args.end_year) if args.start_year else None
    todo = [one] if one else PROD["windows"]
    if one and one not in PROD["windows"]:
        raise SystemExit(
            f"{one[0]}_{one[1]} is not a --product {args.product} window "
            f"({', '.join(f'{a}_{b}' for a, b in PROD['windows'])}). "
            "For the full period the rate window IS the change window, so the "
            "answer is --product total.")

    apply_style()
    loaded, summaries, written = {}, [], []
    for window in todo:
        dom, cs, meta, years, dune_years = load(window)
        st = stats(dom, years, dune_years)
        loaded[window] = (dom, cs, meta, years, dune_years)
        summaries.append(st)
        out_dir = PROD["root"] / f"{window[0]}_{window[1]}"
        out_dir.mkdir(parents=True, exist_ok=True)
        dom.to_csv(out_dir / "domain_comparison.csv", index=False)
        half, w = window_figures(window, dom, cs, meta, years, dune_years, st, out_dir)
        write_provenance(window, out_dir, meta, years, dune_years, st, half)
        written += w
        print(f"{window[0]}-{window[1]}  x{years:.0f} yr (dune {dune_years:.2f})  "
              f"shoreline {st['mean_shoreline_change_m']:+6.1f} m  "
              f"dune {st['mean_dune_change_m']:+6.1f} m  "
              f"beach {st['mean_beach_width_change_m']:+6.1f} m "
              f"({st['mean_beach_width_dune_interval_m']:+.1f} m interval-matched)  "
              f"r={st['r_dune_vs_shoreline']:.2f}")

    if one:
        for p in written:
            print("wrote    ", p.relative_to(_REPO))
        return 0

    if PROD["stacked"] == "halves":
        stack_out = PROD["root"] / "all_windows_stacked"
        stack_out.mkdir(parents=True, exist_ok=True)
        written += halves_overlay_figure(loaded, stack_out)
        print(f"\n(no change_between_periods for --product {args.product}: both "
              f"halves carry the SAME {PROD['rate_window'](1996, 2010)[0]}-"
              f"{PROD['rate_window'](1996, 2010)[1]} rate over the same 14 yr, so "
              "the difference between them is zero by construction on the "
              "shoreline side)")
        for p in written:
            print("wrote    ", p.relative_to(_REPO))
        return 0

    # chain
    chain_out = PROD["root"] / "all_windows_stacked"
    chain_out.mkdir(parents=True, exist_ok=True)
    sup = support_dir(chain_out)
    rows = [loaded[w][0].assign(window=f"{w[0]}_{w[1]}") for w in CHAIN_WINDOWS]
    cols = ["window", "domain_number", "lrr_interval_yr", "dune_interval_yr",
            "mean_lrr_m_yr", "shoreline_change_m", "shoreline_change_dune_interval_m",
            "dune_change_m", "dune_rate_m_yr", "beach_width_change_m",
            "beach_width_change_dune_interval_m", "n_coastsat_transects",
            "n_dune_transects"]
    pd.concat(rows)[cols].to_csv(sup / "domain_comparison.csv", index=False)
    summ = pd.DataFrame(summaries)
    summ.insert(0, "window", [f"{s}_{e}" for s, e in WINDOWS])
    summ.to_csv(sup / "island_summary.csv", index=False)
    half, w = chain_figure(loaded, summaries, chain_out)
    written += w
    # The two-panel overlay of the halves ALONGSIDE the three-panel chain, not
    # instead of it (Hannah, 2026-09-22). The chain answers "how does the whole
    # period relate to its halves"; this one is the direct counterpart of the
    # projected product's sheet, so the two can be laid side by side to see what
    # changes when the rate is fitted per half instead of over the full record.
    written += halves_overlay_figure(loaded, chain_out)

    # difference
    diff_out = PROD["root"] / "change_between_periods"
    diff_out.mkdir(parents=True, exist_ok=True)
    d = difference_table(loaded)
    d.to_csv(support_dir(diff_out) / "domain_difference.csv")
    r_m, r_rate, dhalf, w = difference_figure(d, diff_out)
    written += w

    print(f"\nchain y bounds +/-{half:g} m; difference +/-{dhalf:g} m, "
          f"r={r_m:.2f} (m) / {r_rate:.2f} (m/yr)")
    print(summ[["window", "mean_shoreline_change_m", "mean_dune_change_m",
                "mean_beach_width_change_m", "mean_beach_width_dune_interval_m",
                "r_dune_vs_shoreline", "slope_dune_on_shoreline",
                "domains_same_sign"]].to_string(index=False))
    for p in written:
        print("wrote    ", p.relative_to(_REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
