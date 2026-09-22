"""
projected_vs_duneline.py
==============================================================================
The shoreline change the long-term rate projects against the dune line's net
change, per GIS domain, 1996-2024. Built 2026-09-19 (Hannah, by interview),
the dune-line counterpart of 3-rates/coastsat/lrr_projected.

SHORELINE   per CoastSat transect, the 1996-2024 LRR (3-rates/coastsat/lrr)
            times the DUNE-LINE interval, 1997-10-12 to 2023-07-01
            (25.72 yr), not the 28 yr of lrr_projected, so both changes cover
            the same dates. Per domain the mean of its transects.
DUNE LINE   3-rates/duneline/endpoint/1996_2024: the 2023 digitized line minus
            the 1997 line along the 100 m transects, per domain the mean
            (read, not recomputed). The 2023 flight date is assumed (1 July).
BEACH WIDTH shoreline change - dune-line change. The shoreline is the seaward
            edge of the beach and the dune line its landward edge, so a
            positive value is a beach that WIDENED (the shoreline moved
            seaward of where the dune line went), negative one that narrowed.

The two observations use different transects (CoastSat's ~10 per domain,
the dune line's ~5), so they meet only as domain means. SEAWARD POSITIVE.

THREE FIGURES (Hannah asked for all three presentations)
    two_panel    top: the shoreline (house-style fill, transect dots) with the
                 dune line as a black line; bottom: beach-width change.
    shaded_gap   both as lines on one panel, shoreline blue and dune line red
                 as in net_change_1996_2024, the gap between them grey.
    overlay      the top panel of two_panel on its own.
    One y range for all three.

OUTPUT   data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/projected/1996_2024/
    domain_comparison.csv       per domain: mean LRR, projected shoreline
                                change, dune-line change, beach-width change,
                                the transect counts
    projected_vs_duneline_{two_panel,shaded_gap,overlay}.png  (PDF and
                                CAPTIONS.md under supporting/)
    PROVENANCE.md

USAGE
    python scripts/input_prep/5-scr/projected_vs_duneline/projected_vs_duneline.py
==============================================================================
"""

from __future__ import annotations

import datetime as dt
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
for _sub in ("", "CoastSat", "duneline_vs_coastsat"):
    sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / _sub))

import rates_figures as rf  # noqa: E402  (the 3-rates drawing helpers)
from rates_figures import cw, plt  # noqa: E402
from duneline_vs_coastsat import beach_width_handles, shade_beach_width  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK, apply_style, caption, figsize, save,
)
from site_layer.hat_observed_rates import (  # noqa: E402
    COASTSAT_LRR_ROOT, DUNELINE_ENDPOINT_ROOT, ENDPOINT_DOMAIN_FILE,
    ENDPOINT_TRANSECT_FILE, SHORELINE_VS_DUNELINE,
)

N = cw.N_DOMAINS
START, END = 1996, 2024
C_SHORE = "#2166ac"     # net_change_1996_2024's shoreline blue
C_DUNE = "#b2182b"      # and its dune red
LW = 1.1
Y_LABEL = "Net change in position (m)"


def load():
    w = f"{START}_{END}"
    dune_t = pd.read_csv(DUNELINE_ENDPOINT_ROOT / w / ENDPOINT_TRANSECT_FILE)
    dune_d = pd.read_csv(DUNELINE_ENDPOINT_ROOT / w / ENDPOINT_DOMAIN_FILE)
    meta = dune_t.iloc[0]
    years = float(meta["interval_yr"])
    cs = pd.read_csv(COASTSAT_LRR_ROOT / w / "transect_lrr_full.csv")
    cs = cs[cs["domain_number"].between(1, N)].copy()
    cs["domain_number"] = cs["domain_number"].astype(int)
    cs["projected_change_m"] = cs["lrr_m_yr"] * years
    g = cs.groupby("domain_number")
    dom = pd.DataFrame({
        "n_coastsat_transects": g.size(),
        "mean_lrr_m_yr": g["lrr_m_yr"].mean(),
        "shoreline_projected_change_m": g["projected_change_m"].mean(),
    }).reindex(range(1, N + 1))
    d = dune_d.set_index(dune_d["domain_number"].astype(int))
    dom["n_dune_transects"] = d["n_transects"].reindex(dom.index)
    dom["duneline_change_m"] = d["mean_change_m"].reindex(dom.index)
    dom["beach_width_change_m"] = dom["shoreline_projected_change_m"] - dom["duneline_change_m"]
    dom = dom.round(3).rename_axis("domain_number").reset_index()
    dom.insert(1, "interval_yr", round(years, 4))
    return dom, cs, meta, years


def _frame(dom, col):
    return rf._frame(dom, col)


def _empty(dom):
    return rf._frame(dom.assign(_nan=np.nan), "_nan")


def _marks(ax, half, label):
    cw.draw_shoals(ax, label=label)
    fills = cw.fills_in(START, END)
    if fills and label:
        cw.draw_fills(ax, fills, half)


def _shoreline_panel(ax, dom, cs, half, label=True):
    """House-style shoreline fill and dots, the dune line as a black line."""
    cw.draw_panel(ax, _frame(dom, "shoreline_projected_change_m"), half,
                  label=label, std=False)
    tt, x = rf._along(cs)
    n_out = rf._dots(ax, x, tt["projected_change_m"].to_numpy(float), half)
    ax.plot(dom["domain_number"], dom["duneline_change_m"], color=INK, lw=LW,
            zorder=12, solid_capstyle="round")
    _marks(ax, half, label)
    return n_out


def _gap_panel(ax, dom, half, label=True):
    cw.draw_panel(ax, _empty(dom), half, label=label, std=False)
    x = dom["domain_number"].to_numpy(float)
    ys = dom["shoreline_projected_change_m"].to_numpy(float)
    yd = dom["duneline_change_m"].to_numpy(float)
    shade_beach_width(ax, x, ys, yd)
    ax.plot(x, yd, color=C_DUNE, lw=LW, zorder=12)
    ax.plot(x, ys, color=C_SHORE, lw=LW, zorder=12)
    _marks(ax, half, label)


def _width_panel(ax, dom, half, label=False):
    cw.draw_panel(ax, _empty(dom), half, label=label, std=False)
    x = dom["domain_number"].to_numpy(float)
    y = dom["beach_width_change_m"].to_numpy(float)
    shade_beach_width(ax, x, y, np.zeros_like(y))
    ax.plot(x, y, color=INK, lw=LW, zorder=12)
    _marks(ax, half, label)


def _house_handles():
    return [(Line2D([], [], color=cw.C_ACCRETE, lw=1.0), Line2D([], [], color=cw.C_ERODE, lw=1.0)),
            (Line2D([], [], color=cw.C_ACCRETE, marker="o", ms=2.2, lw=0),
             Line2D([], [], color=cw.C_ERODE, marker="o", ms=2.2, lw=0)),
            Line2D([], [], color=INK, lw=LW)]


def _legend(fig, handles, labels, ncol):
    fig.legend(handles, labels, loc="outside lower center", ncol=ncol, frameon=False,
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})


def _base_caption(meta, years):
    return (
        f"The shoreline change is projected: for each CoastSat transect the "
        f"{START}–{END} linear regression rate (the ordinary-least-squares slope "
        f"through every satellite position from 1 January {START} to 31 December "
        f"{END}) multiplied by {years:.1f} yr, the interval between the two dune "
        f"lines ({meta['start_date']} to {meta['end_date']}, the end date assumed), "
        "averaged over the ~10 transects in each 500 m domain. The dune-line change "
        f"is observed: the {int(meta['end_vintage'])} digitized dune line minus the "
        f"{int(meta['start_vintage'])} line along the 100 m transects, averaged over "
        "the ~5 in each domain. Both cover the same dates and are seaward positive, "
        "in metres. Beach-width change is shoreline change minus dune-line change: "
        "positive where the beach widened, negative where it narrowed. "
        + rf._marks_clause(START, END))


def figures(dom, cs, meta, years):
    ext = np.nanmax(np.abs(dom[["shoreline_projected_change_m", "duneline_change_m",
                                "beach_width_change_m"]].to_numpy(float)))
    half = float(math.ceil((ext + 5) / rf.Y_STEP_M) * rf.Y_STEP_M)
    tick = 20.0 if half > 60 else 10.0
    out_dir = (SHORELINE_VS_DUNELINE / "superseded_20260921"
               / f"projected_dune_interval_{START}_{END}")
    out_dir.mkdir(parents=True, exist_ok=True)
    ok = dom.dropna(subset=["shoreline_projected_change_m", "duneline_change_m"])
    r = np.corrcoef(ok["shoreline_projected_change_m"], ok["duneline_change_m"])[0, 1]
    bw = ok["beach_width_change_m"]
    stats = (f" Over the {len(ok)} domains, r(shoreline, dune line) = {r:.2f}; the "
             f"beach widened by {bw.mean():+.1f} m on average (range {bw.min():+.1f} "
             f"to {bw.max():+.1f} m). The y axis is ±{half:g} m on every panel.")
    written = []

    # 1. two panels
    fig, axes = plt.subplots(2, 1, sharex=True, constrained_layout=True,
                             figsize=figsize("double", height=5.6))
    n_out = _shoreline_panel(axes[0], dom, cs, half, label=True)
    _width_panel(axes[1], dom, half)
    for ax in axes:
        ax.yaxis.set_major_locator(MultipleLocator(tick))
    axes[0].set_ylabel(Y_LABEL)
    axes[1].set_ylabel("Beach-width change (m)")
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    from site_layer.hat_figure_style import _title
    _title(axes[0], 0, "Shoreline and dune-line change")
    _title(axes[1], 1, "Beach-width change")
    _legend(fig, _house_handles() + beach_width_handles(),
            [f"Projected shoreline change (LRR × {years:.1f} yr)",
             "Individual CoastSat transects", "Dune-line change (endpoint)",
             "Beach widened", "Beach narrowed"], ncol=3)
    caption(fig, (
        "Projected shoreline change against observed dune-line change by GIS domain "
        f"(1 at Cape Point, 90 at Pea Island), {meta['start_date'][:4]}–"
        f"{meta['end_date'][:4]}. (a) The shoreline as the coloured line and fill "
        "(blue seaward, red landward) with its single transects as dots"
        + (f" ({n_out} beyond the axis, open circles at its edge)" if n_out else "")
        + "; the dune line as the black line. (b) Beach-width change, the shoreline "
        "change minus the dune-line change: solid grey where the beach widened, hatched where it narrowed. " + _base_caption(meta, years)
        + stats))
    written += save(fig, out_dir / "projected_vs_duneline_two_panel")
    plt.close(fig)

    # 2. shaded gap
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40), constrained_layout=True)
    _gap_panel(ax, dom, half)
    ax.yaxis.set_major_locator(MultipleLocator(tick))
    ax.set_ylabel(Y_LABEL)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    _legend(fig, [Line2D([], [], color=C_SHORE, lw=LW), Line2D([], [], color=C_DUNE, lw=LW)]
            + beach_width_handles(),
            [f"Projected shoreline change (LRR × {years:.1f} yr)",
             "Dune-line change (endpoint)", "Beach widened", "Beach narrowed"], ncol=2)
    caption(fig, (
        "Projected shoreline change (blue) and observed dune-line change (red) by GIS "
        "domain (1 at Cape Point, 90 at Pea Island), domain means. The space between "
        "them is the beach-width change: solid grey where the beach widened (the blue "
        "line above the red), hatched where it narrowed. "
        + _base_caption(meta, years) + stats))
    written += save(fig, out_dir / "projected_vs_duneline_shaded_gap")
    plt.close(fig)

    # 3. overlay
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40), constrained_layout=True)
    n_out = _shoreline_panel(ax, dom, cs, half, label=True)
    ax.yaxis.set_major_locator(MultipleLocator(tick))
    ax.set_ylabel(Y_LABEL)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    _legend(fig, _house_handles(),
            [f"Projected shoreline change (LRR × {years:.1f} yr)",
             "Individual CoastSat transects", "Dune-line change (endpoint)"], ncol=2)
    caption(fig, (
        "Projected shoreline change against observed dune-line change by GIS domain "
        "(1 at Cape Point, 90 at Pea Island): the shoreline as the coloured line and "
        "fill (blue seaward, red landward) with its single transects as dots"
        + (f" ({n_out} beyond the axis, open circles at its edge)" if n_out else "")
        + "; the dune line as the black line. Where the black line is below the "
        "fill's edge, the beach widened. " + _base_caption(meta, years) + stats))
    written += save(fig, out_dir / "projected_vs_duneline_overlay")
    plt.close(fig)
    return out_dir, r, half, written


def main() -> int:
    apply_style()
    dom, cs, meta, years = load()
    out_dir = (SHORELINE_VS_DUNELINE / "superseded_20260921"
               / f"projected_dune_interval_{START}_{END}")
    out_dir.mkdir(parents=True, exist_ok=True)
    dom.to_csv(out_dir / "domain_comparison.csv", index=False)
    out_dir, r, half, _ = figures(dom, cs, meta, years)
    ok = dom.dropna(subset=["shoreline_projected_change_m", "duneline_change_m"])
    bw = ok["beach_width_change_m"]
    (out_dir / "PROVENANCE.md").write_text("\n".join([
        f"# 4-comparisons/shoreline_vs_duneline/projected/{START}_{END} - provenance",
        "",
        f"Written {dt.datetime.now():%Y-%m-%d %H:%M} by "
        "scripts/input_prep/5-scr/projected_vs_duneline/projected_vs_duneline.py.",
        "",
        f"- Shoreline: `3-rates/coastsat/lrr/{START}_{END}/transect_lrr_full.csv`, "
        f"lrr_m_yr x {years:.4f} yr (the dune-line interval, so both cover "
        f"{meta['start_date']} to {meta['end_date']}; the end date is ASSUMED). "
        "Not the 28-yr `3-rates/coastsat/lrr_projected`.",
        f"- Dune line: `3-rates/duneline/endpoint/{START}_{END}/` as stored "
        f"({int(meta['start_vintage'])} and {int(meta['end_vintage'])} lines).",
        "- Beach-width change = shoreline change - dune-line change; positive = widened.",
        "- Domain means only: the two use different transects.",
        "",
        "## Island summary",
        "",
        f"Domain mean shoreline {ok['shoreline_projected_change_m'].mean():+.1f} m, "
        f"dune line {ok['duneline_change_m'].mean():+.1f} m, beach width "
        f"{bw.mean():+.1f} m (range {bw.min():+.1f} to {bw.max():+.1f}); beach "
        f"narrowed in {int((bw < 0).sum())} of {len(ok)} domains. "
        f"r(shoreline, dune line) = {r:.2f}. y axis ±{half:g} m.",
        "",
    ]), encoding="utf-8")
    print(f"x{years:.2f} yr  r={r:.2f}  shoreline {ok['shoreline_projected_change_m'].mean():+.1f} m  "
          f"dune {ok['duneline_change_m'].mean():+.1f} m  beach {bw.mean():+.1f} m "
          f"(narrowed in {int((bw < 0).sum())})  -> {out_dir.relative_to(_REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
