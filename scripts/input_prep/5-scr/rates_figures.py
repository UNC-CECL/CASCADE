"""
rates_figures.py
==============================================================================
One house-style figure per window for every product under 5-scr/3-rates/,
written beside its tables (Hannah, 2026-09-18: "I wanted all of these to have
figures"). This replaces the autoscaled quick-looks the LRR fit used to draw
(archived in 5-scr/archive/coastsat_lrr_quicklooks_20260918/).

WHAT EACH FIGURE SHOWS
    The per-domain value as the sign-coloured line and fill of the
    coastsat_lrr_windows panel (blue seaward, red landward; the drawing is
    imported), the individual transects behind it as small dots coloured by
    their OWN sign (the same blue / red; Hannah, 2026-09-18), and the
    village bands, groin and piers, the offshore shoals as faint hatched boxes,
    and the model-input beach fills inside the window as bars above the panel.

        coastsat/lrr/<w>/lrr_<w>.png                    m/yr, +/-1 std dotted
        coastsat/endpoint/<w>/coastsat_endpoint_<w>.png m
        duneline/endpoint/<w>/duneline_endpoint_<w>.png m
        coastsat/5yr_bins/<w>/lrr_5yr_bins_<w>.png      m/yr, one panel per bin
                                                        (coastsat_5yr_bins_figure.py)

    and per MODEL CHAIN (1984-2004-2024, 1996-2010-2024) the chain's two
    windows stacked, earlier above, on the same axis (2026-09-18):
        <product root>/chains/<stem>_chain_<y0>_<y1>_<y2>.png
    for lrr, coastsat endpoint and duneline endpoint. 5yr_bins has no chain
    figure: it covers only the 1996 chain, and its 1996_2024 figure is it.

Y AXES
    lrr        the bound the window figures use: the largest |domain mean|
               over every window plus 1 m, rounded up (+/-8 m/yr today)
    endpoint   ONE bound for both endpoint products, the smallest multiple of
               10 m holding every domain mean of every window, so a shoreline
               figure reads against its dune-line figure directly
    The transect dots are NOT in the bound; a dot outside it is drawn at the
    edge as an open marker and counted in the caption.

USAGE
    python scripts/input_prep/5-scr/rates_figures.py
==============================================================================
"""

from __future__ import annotations

import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "CoastSat"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

import coastsat_lrr_windows as cw  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, apply_style, caption, figsize, save,
)
from site_layer.hat_observed_rates import (  # noqa: E402
    COASTSAT_ENDPOINT_ROOT, COASTSAT_LRR_ROOT, DUNELINE_ENDPOINT_ROOT,
    ENDPOINT_DOMAIN_FILE, ENDPOINT_TRANSECT_FILE, windows,
)
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS  # noqa: E402

N = cw.N_DOMAINS
# Each transect dot takes the colour of its own sign, the line's blue / red
# (the model_vs_observed per-domain dots use the same pair). Grey until 2026-09-18.
DOT_ALPHA = 0.4
DOT_S = 3.0
Y_STEP_M = 10.0
BINS_SCRIPT = (_REPO / "scripts" / "input_prep" / "5-scr" / "CoastSat_timeseries"
               / "coastsat_5yr_bins_figure.py")


def _windows(root):
    return sorted(p.name for p in root.iterdir()
                  if p.is_dir() and p.name[:4].isdigit() and "_" in p.name)


def _along(t):
    """x for each transect: its domain plus an even spread inside it."""
    t = t.sort_values(["domain_number"]).reset_index(drop=True)
    rank = t.groupby("domain_number").cumcount()
    n = t.groupby("domain_number")["domain_number"].transform("count")
    return t, (t["domain_number"] - 0.5 + (rank + 0.5) / n).to_numpy(float)


def _dots(ax, x, y, half):
    ok = np.isfinite(y)
    inside = ok & (np.abs(y) <= half)
    col = np.where(y < 0, cw.C_ERODE, cw.C_ACCRETE)
    ax.scatter(x[inside], y[inside], s=DOT_S, c=col[inside], alpha=DOT_ALPHA,
               linewidths=0, zorder=3.5)
    out = ok & ~inside
    if out.any():
        ax.scatter(x[out], np.clip(y[out], -half * 0.985, half * 0.985), s=9,
                   facecolors="none", edgecolors=col[out], linewidths=0.6, zorder=3.5)
    return int(out.sum())


def _frame(domain_df, value, std=None):
    df = pd.DataFrame({"domain_number": np.arange(1, N + 1)})
    d = domain_df.set_index("domain_number")
    df["mean_lrr"] = [d[value].get(g, np.nan) for g in range(1, N + 1)]
    df["std_lrr"] = ([d[std].get(g, 0.0) for g in range(1, N + 1)] if std else 0.0)
    return df


def _marks_clause(start, end):
    fills = cw.fills_in(start, end)
    shoals = "; ".join(f"{n} GIS {lo}–{hi}" for n, (lo, hi)
                       in HATTERAS_ANNOTATIONS.shoal_zones.items())
    fill_txt = ("Black bars above the panel mark the beach fills placed inside the "
                "window, at the footprint the hindcast uses ("
                + "; ".join(f"{y} at GIS {lo}–{hi}" for y, lo, hi in fills) + "). "
                if fills else "")
    return (fill_txt + f"Hatched amber boxes mark the offshore shoals ({shoals}). "
            "Village spans are shaded; the solid hairline is the Buxton groin and "
            "the dotted hairlines are the Avon and Rodanthe piers.")


def _draw(frame, x, y, half, y_label, tick, fills, std):
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40),
                           constrained_layout=True)
    cw.draw_panel(ax, frame, half, std=std)
    n_out = _dots(ax, x, y, half)
    cw.draw_shoals(ax, label=True)
    if fills:
        cw.draw_fills(ax, fills, half)
    ax.yaxis.set_major_locator(MultipleLocator(tick))
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel(y_label)
    return fig, ax, n_out


def _legend(fig, what, std):
    h = [(Line2D([], [], color=cw.C_ACCRETE, lw=1.0), Line2D([], [], color=cw.C_ERODE, lw=1.0)),
         (Line2D([], [], color=cw.C_ACCRETE, marker="o", ms=2.2, lw=0),
          Line2D([], [], color=cw.C_ERODE, marker="o", ms=2.2, lw=0))]
    # Sentence case, the quantity named by its method; "domain mean" and the
    # seaward / landward colours are stated in each caption (Hannah, 2026-09-19).
    labels = [what, "Individual transects"]
    if std:
        h.append(Line2D([], [], color=cw.INK_MUTED, lw=0.5, ls=(0, (1, 1.6))))
        labels.append("±1 standard deviation")
    from matplotlib.legend_handler import HandlerTuple
    fig.legend(h, labels, loc="outside lower center", ncol=len(h), frameon=False,
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})


# -----------------------------------------------------------------------------
def lrr_figures():
    wins = windows()
    frames = {w: pd.read_csv(COASTSAT_LRR_ROOT / f"{w[0]}_{w[1]}" / "domain_lrr_summary.csv")
              for w in wins}
    half = cw.shared_bounds([_frame(f.assign(domain_number=f.domain_number.astype(int)),
                                    "mean_lrr") for f in frames.values()])
    out = []
    for (s, e), dom in frames.items():
        dom["domain_number"] = dom["domain_number"].astype(int)
        t = pd.read_csv(COASTSAT_LRR_ROOT / f"{s}_{e}" / "transect_lrr_full.csv")
        t = t[t["domain_number"].between(1, N)]
        t, x = _along(t)
        fig, ax, n_out = _draw(_frame(dom, "mean_lrr", "std_lrr"), x,
                               t["lrr_m_yr"].to_numpy(float), half,
                               cw.Y_LABEL, cw.Y_TICK_M, cw.fills_in(s, e), std=True)
        _legend(fig, "Shoreline change rate (LRR)", std=True)
        caption(fig, (
            f"Observed shoreline change rate by GIS domain (1 at Cape Point, 90 at "
            f"Pea Island), {s}–{e}: for each CoastSat transect the ordinary-least-"
            "squares slope of shoreline position against date over the calendar "
            f"window (1 January {s} to 31 December {e}, at least 3 positions), "
            "seaward positive. The coloured line and fill are the mean of the ~10 "
            "transects in each 500 m domain, blue where the shoreline moved seaward "
            "and red where it moved landward; the dotted lines are ±1 standard "
            "deviation across them; the dots are the single transects, blue or "
            "red by their own sign"
            + (f" ({n_out} beyond the axis, drawn as open circles at its edge)" if n_out else "")
            + ". " + _marks_clause(s, e)
            + f" The y axis is ±{half:g} m/yr, the bound of every window figure "
            "(the largest |domain mean| over all windows plus 1 m, rounded up)."
            + (" This window is context: no model run is graded against it."
               if (s, e) == (1996, 2024) else "")))
        out += save(fig, COASTSAT_LRR_ROOT / f"{s}_{e}" / f"lrr_{s}_{e}")
        plt.close(fig)
    return out, half


def endpoint_figures():
    products = [("coastsat", COASTSAT_ENDPOINT_ROOT, "shoreline"),
                ("duneline", DUNELINE_ENDPOINT_ROOT, "dune-line")]
    data = {}
    for key, root, _ in products:
        for w in _windows(root):
            data[(key, w)] = (pd.read_csv(root / w / ENDPOINT_DOMAIN_FILE),
                              pd.read_csv(root / w / ENDPOINT_TRANSECT_FILE))
    extreme = max(float(np.nanmax(np.abs(d["mean_change_m"]))) for d, _ in data.values())
    half = float(math.ceil(extreme / Y_STEP_M) * Y_STEP_M)
    tick = 10.0 if half <= 60 else 20.0
    out = []
    for key, root, what in products:
        for w in _windows(root):
            dom, tr = data[(key, w)]
            s, e = (int(v) for v in w.split("_"))
            meta = tr.iloc[0]
            v0, v1 = int(meta["start_vintage"]), int(meta["end_vintage"])
            tr = tr[tr["domain_number"].between(1, N)]
            tr, x = _along(tr)
            fig, ax, n_out = _draw(_frame(dom, "mean_change_m"), x,
                                   tr["change_m"].to_numpy(float), half,
                                   f"Net change in {what} position (m)", tick,
                                   cw.fills_in(v0, v1), std=False)
            _legend(fig, f"Net {what} change (endpoint)", std=False)
            dates = (f"{meta['start_date']} to {meta['end_date']}"
                     + (" (the end date assumed; no flight date is known)"
                        if bool(meta["end_date_assumed"]) else ""))
            if key == "coastsat":
                how = ("for each CoastSat transect, the mean satellite shoreline position "
                       "within six months of the end dune-line image date minus the "
                       f"same about the start date ({dates}); the ~10 transects of each "
                       "500 m domain are averaged")
            else:
                how = (f"the {v1} digitized dune line minus the {v0} line, each "
                       "measured along the 100 m transects from a fixed offshore "
                       f"datum ({dates}); the ~5 transects of each 500 m domain are "
                       "averaged")
            caption(fig, (
                f"Net change in {what} position by GIS domain (1 at Cape Point, 90 at "
                f"Pea Island), {v0}–{v1}, the lines standing in for the model years "
                f"{s} and {e}: {how}. Seaward positive, in metres; this is net "
                "displacement, not a rate. The coloured line and fill are the domain "
                "means, blue seaward and red landward; the dots are the single "
                "transects, blue or red by their own sign"
                + (f" ({n_out} beyond the axis, drawn as open circles at its edge)" if n_out else "")
                + ". " + _marks_clause(v0, v1)
                + f" The y axis is ±{half:g} m, shared by every CoastSat and dune-line "
                "endpoint figure so the two can be read against each other."))
            out += save(fig, root / w / f"{key}_endpoint_{w}")
            plt.close(fig)
    return out, half


# -----------------------------------------------------------------------------
# the two model chains, one figure each per product (Hannah, 2026-09-18)
# -----------------------------------------------------------------------------
CHAINS = [((1984, 2004), (2004, 2024)), ((1996, 2010), (2010, 2024))]


def _load(product, s, e):
    """(frame, x, y, fills, vintages) for one product and window."""
    if product == "lrr":
        root = COASTSAT_LRR_ROOT / f"{s}_{e}"
        dom = pd.read_csv(root / "domain_lrr_summary.csv")
        dom["domain_number"] = dom["domain_number"].astype(int)
        t = pd.read_csv(root / "transect_lrr_full.csv")
        t = t[t["domain_number"].between(1, N)]
        t, x = _along(t)
        return (_frame(dom, "mean_lrr", "std_lrr"), x, t["lrr_m_yr"].to_numpy(float),
                cw.fills_in(s, e), (s, e))
    root = (COASTSAT_ENDPOINT_ROOT if product == "coastsat" else DUNELINE_ENDPOINT_ROOT) / f"{s}_{e}"
    dom = pd.read_csv(root / ENDPOINT_DOMAIN_FILE)
    t = pd.read_csv(root / ENDPOINT_TRANSECT_FILE)
    v0, v1 = int(t.iloc[0]["start_vintage"]), int(t.iloc[0]["end_vintage"])
    t = t[t["domain_number"].between(1, N)]
    t, x = _along(t)
    return (_frame(dom, "mean_change_m"), x, t["change_m"].to_numpy(float),
            cw.fills_in(v0, v1), (v0, v1))


def chain_figures(half_lrr, half_end):
    """One figure per chain per product: the chain's two windows stacked,
    earlier above, on the product's shared axis (the same bound as its
    single-window figures). Written to <product root>/chains/."""
    from site_layer.hat_figure_style import _title
    products = [
        ("lrr", COASTSAT_LRR_ROOT, half_lrr, cw.Y_LABEL, cw.Y_TICK_M, True,
         "Shoreline change rate (LRR)", "lrr"),
        ("coastsat", COASTSAT_ENDPOINT_ROOT, half_end, "Net change in shoreline position (m)",
         10.0 if half_end <= 60 else 20.0, False, "Net shoreline change (endpoint)", "coastsat_endpoint"),
        ("duneline", DUNELINE_ENDPOINT_ROOT, half_end, "Net change in dune-line position (m)",
         10.0 if half_end <= 60 else 20.0, False, "Net dune-line change (endpoint)", "duneline_endpoint"),
    ]
    out = []
    for product, root, half, ylab, tick, std, what, stem in products:
        for chain in CHAINS:
            s0, e1 = chain[0][0], chain[-1][1]
            fig, axes = plt.subplots(len(chain), 1, sharex=True, sharey=True,
                                     constrained_layout=True,
                                     figsize=figsize("double", height=5.2))
            titles, n_out = [], 0
            for i, (ax, (s, e)) in enumerate(zip(axes, chain)):
                frame, x, y, fills, (v0, v1) = _load(product, s, e)
                cw.draw_panel(ax, frame, half, label=(i == 0), std=std)
                n_out += _dots(ax, x, y, half)
                cw.draw_shoals(ax, label=(i == 0))
                if fills:
                    cw.draw_fills(ax, fills, half)
                ax.yaxis.set_major_locator(MultipleLocator(tick))
                label = (f"{s}–{e}" if product == "lrr"
                         else f"{v0}–{v1}" + ("" if (v0, v1) == (s, e) else f" (for {s}–{e})"))
                _title(ax, i, label)
                titles.append(label)
            axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
            fig.supylabel(ylab, fontsize=9)
            _legend(fig, what, std=std)
            if product == "lrr":
                body = ("the OLS slope of each CoastSat transect's shoreline position "
                        "against date over the calendar window, averaged per 500 m "
                        "domain (the line and fill; dotted ±1 standard deviation)")
                unit = f"±{half:g} m/yr, the bound of every window figure"
            elif product == "coastsat":
                body = ("the mean CoastSat position within six months of the end "
                        "dune-line image date minus the same about the start date, "
                        "averaged per 500 m domain (the line and fill), in metres")
                unit = f"±{half:g} m, shared with every endpoint figure"
            else:
                body = ("the end digitized dune line minus the start line along the "
                        "100 m transects, averaged per 500 m domain (the line and "
                        "fill), in metres")
                unit = f"±{half:g} m, shared with every endpoint figure"
            caption(fig, (
                f"The {s0} → {chain[0][1]} → {e1} model chain, "
                + " above ".join(f"({chr(97 + i)}) {t}" for i, t in enumerate(titles))
                + f", by GIS domain (1 at Cape Point, 90 at Pea Island): {body}. "
                "Seaward positive; blue and filled where the feature moved seaward, "
                "red where it moved landward; the dots are the single transects, "
                "blue or red by their own sign"
                + (f" ({n_out} beyond the axis, drawn as open circles at its edge)"
                   if n_out else "")
                + ". Black bars above a panel mark the beach fills placed in that "
                "window at the footprint the hindcast uses; hatched amber boxes mark "
                "the offshore shoals; village spans are shaded; the solid hairline is "
                "the Buxton groin and the dotted hairlines are the Avon and Rodanthe "
                f"piers. Both panels share a y axis of {unit}."
                + (" The 2023 dune-line flight date is not known and is assumed to "
                   "be 1 July." if product != "lrr" and e1 == 2024 else "")))
            out += save(fig, root / "chains" / f"{stem}_chain_{s0}_{chain[0][1]}_{e1}")
            plt.close(fig)
    return out


def main() -> int:
    apply_style()
    written, half = lrr_figures()
    print(f"coastsat/lrr        {len(written) // 2} figures, y +/-{half:g} m/yr")
    w2, half2 = endpoint_figures()
    print(f"*/endpoint          {len(w2) // 2} figures, y +/-{half2:g} m")
    w3 = chain_figures(half, half2)
    print(f"*/chains            {len(w3) // 2} figures (1984-2004-2024, 1996-2010-2024)")
    r = subprocess.run([sys.executable, str(BINS_SCRIPT)], capture_output=True,
                       text=True, encoding="utf-8",
                       env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    print(r.stdout.strip() or r.stderr.strip())
    return r.returncode


if __name__ == "__main__":
    sys.exit(main())
