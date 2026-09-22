"""
coastsat_5yr_bins_figure.py
==============================================================================
When, inside a window, did the shoreline change? The CoastSat LRR in
successive 5-year bins, one panel per bin, per GIS domain, in the house style.
Written 2026-09-18 when the bins were rebuilt on the 1996 -> 2010 -> 2024
chain; it writes beside the table, as every 3-rates product does.

READS    3-rates/coastsat/5yr_bins/<window>/lrr_bins_5yr.csv
         (coastsat_lrr_5year_bins.py: per transect an OLS over the bin's
         positions, domain MEAN, |rate| > 50 m/yr dropped, bins under 3.75 yr
         dropped)
DRAWS    one stacked panel per bin, filled blue where the shoreline moved
         seaward and red where it moved landward (the coastsat_lrr_windows
         panel drawing, imported); ONE y axis for all three figures, the largest
         |rate| over GIS 2-90 plus 1 m rounded up, with GIS 1 clipped and
         labelled where it runs past (Hannah, 2026-09-18). Village bands, the groin and
         piers, the offshore shoals as faint hatched boxes; a model-input
         beach fill is marked above the panel of the bin it falls in.
WRITES   3-rates/coastsat/5yr_bins/<window>/lrr_5yr_bins_<window>.png
         (+ PDF and CAPTIONS.md under supporting/)

USAGE
    python scripts/input_prep/5-scr/3-rates/coastsat/5yr_bins/coastsat_5yr_bins_figure.py
==============================================================================
"""

from __future__ import annotations

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

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

import coastsat_lrr_windows as cw  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, _title, apply_style, caption, figsize, save,
)
from site_layer.hat_observed_rates import TIMESERIES_LRR  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS  # noqa: E402

WINDOWS = ["1996_2010", "2010_2024", "1996_2024"]
# Beside its table since 2026-09-18 (every 3-rates product carries its own
# figure; Hannah). It sat in 4-comparisons/coastsat_5yr_bins/ for an hour.
OUT_ROOT = TIMESERIES_LRR
N = cw.N_DOMAINS


def load(window):
    """[(label, first_year, last_year, frame), ...] oldest first."""
    df = pd.read_csv(TIMESERIES_LRR / window / "lrr_bins_5yr.csv", index_col=0)
    out = []
    for label, row in df.iterrows():
        a, b = (int(x) for x in str(label).replace("–", "-").split("-")[:2])
        vals = pd.to_numeric(row, errors="coerce")
        frame = pd.DataFrame({"domain_number": np.arange(1, N + 1)})
        frame["mean_lrr"] = [vals.get(str(g), np.nan) for g in range(1, N + 1)]
        frame["std_lrr"] = 0.0
        out.append((f"{a}–{b}", a, b, frame))
    return out


def tick_step(half):
    """A tick every 2, 5 or 10 m/yr, so a panel carries at most ~8 labels."""
    return 2.0 if half <= 8 else 5.0 if half <= 20 else 10.0


# THE AXIS IS CAPPED (Hannah, 2026-09-18). Cape Point (GIS 1) reaches about
# +35 m/yr in the 2020s bins -- the shoal attaching -- and at full scale that
# one domain set a +/-36 axis that flattened every other bin, all of which stay
# within +/-13. The bound is taken over GIS 2-90 of EVERY window, so the three
# figures share it; a value beyond it is clipped at the edge, marked with a
# triangle and labelled with its value, never dropped.
CAP_EXCLUDES = (1,)


def shared_cap():
    vals = [np.nanmax(np.abs(f.loc[~f["domain_number"].isin(CAP_EXCLUDES), "mean_lrr"]))
            for w in WINDOWS for *_, f in load(w)]
    return float(math.ceil(max(vals) + cw.Y_PAD_M))


def mark_clipped(ax, frame, half):
    """A triangle at the axis edge for every domain beyond it, labelled with
    the domain and its value. Returns [(gis, value), ...]."""
    out = []
    for g, v in zip(frame["domain_number"], frame["mean_lrr"]):
        if np.isfinite(v) and abs(v) > half:
            y = half if v > 0 else -half
            ax.plot([g], [y * 0.97], marker="^" if v > 0 else "v", ms=4.5,
                    color=cw.C_ACCRETE if v > 0 else cw.C_ERODE, zorder=7,
                    clip_on=False)
            ax.text(g + 0.9, y * 0.93, f"GIS {int(g)}: {v:+.1f}", fontsize=6.5,
                    ha="left", va="top" if v > 0 else "bottom", color=cw.INK,
                    zorder=7)
            out.append((int(g), float(v)))
    return out


def figure(window, half):
    bins = load(window)
    clipped = []
    n = len(bins)
    fig, axes = plt.subplots(n, 1, sharex=True, sharey=True, squeeze=False,
                             constrained_layout=True,
                             figsize=figsize("double", height=1.45 * n + 0.8))
    axes = axes[:, 0]
    for i, (ax, (label, a, b, frame)) in enumerate(zip(axes, bins)):
        cw.draw_panel(ax, frame, half, label=(i == 0), std=False)
        clipped += [(label, g, v) for g, v in mark_clipped(ax, frame, half)]
        cw.draw_shoals(ax, label=(i == 0))
        fills = cw.fills_in(a, b)
        if fills:
            cw.draw_fills(ax, fills, half)
        _title(ax, i, label)
        # draw_panel ticks every 2 m/yr, which is unreadable once one bin's
        # spike (Cape Point, GIS 1, 2020-2024) sets a +/-36 axis
        ax.yaxis.set_major_locator(MultipleLocator(tick_step(half)))
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    fig.supylabel(cw.Y_LABEL, fontsize=9)

    s, e = window.split("_")
    shoal_txt = "; ".join(f"{nm} GIS {lo}–{hi}" for nm, (lo, hi)
                          in HATTERAS_ANNOTATIONS.shoal_zones.items())
    caption(fig, (
        f"Observed shoreline change rate by GIS domain (1 at Cape Point, 90 at "
        f"Pea Island) in successive 5-year bins of {s}–{e}, oldest at the top. "
        "Each panel is the mean linear regression rate of the CoastSat transects "
        "inside each 500 m domain over that bin's positions only (calendar years "
        "inclusive; a transect needs at least 3 positions, rates beyond ±50 m/yr "
        "are dropped, and a bin shorter than 3.75 years is not drawn), blue and "
        "filled where the shoreline moved seaward, red where it moved landward. "
        "A bin's rate is a short fit through noisy positions, so single-domain "
        "spikes are weak evidence; the pattern across bins is the point. Black "
        "bars mark a beach fill at the footprint the hindcast uses, above the "
        "bin it was placed in; hatched amber boxes mark the offshore shoals "
        f"({shoal_txt}). Village spans are shaded; the solid hairline is the "
        "Buxton groin and the dotted hairlines are the Avon and Rodanthe piers. "
        f"Every panel of all three 5-year-bin figures shares a y axis of "
        f"±{half:g} m/yr: the largest |rate| over GIS 2–90 of every bin plus "
        "1 m, rounded up. Cape Point (GIS 1) is left out of that bound; where "
        "a value runs past the axis it is clipped at the edge and marked with "
        "a triangle and its value"
        + (" (" + "; ".join(f"{lab}: GIS {g} at {v:+.1f} m/yr" for lab, g, v in clipped)
           + ")" if clipped else "")
        + ". The window-length figures in 3-rates/coastsat/lrr use ±8, so these are "
        "not read against them panel for panel."))
    out = save(fig, OUT_ROOT / window / f"lrr_5yr_bins_{window}")
    plt.close(fig)
    return out, half, [b[0] for b in bins]


def main() -> int:
    apply_style()
    cap = shared_cap()
    for w in WINDOWS:
        out, half, labels = figure(w, cap)
        print(f"{w}  bins {', '.join(labels)}  y +/-{half:g}  -> "
              f"{out[0].relative_to(_REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
