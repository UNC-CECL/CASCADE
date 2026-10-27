"""
How Barrier3D's logistic dune growth raises a 1, 5 and 10 m dune over 30 years, at each period's rmin-rmax.

    python scripts/figure_making/model/dune_growth_examples.py

Writes output/figures/4-model-mechanics/barrier3d/barrier3d_dune_growth_examples.png. The growth is
Barrier3D's own equation (barrier3d.py, update_dune_domain), h += r * h * (1 - h / Dmax), on the front dune
row with no storms; r is drawn once per 10 m dune cell between rmin and rmax and kept for the run. Rates
and Dmax are read from the runner and the parameter file at the values of 2026-10-22.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-22
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402

# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(_p for _p in Path(__file__).resolve().parents  if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
sys.path.insert(0, str(_HAT_DIR))
from site_layer.hat_figure_style import (INK_MUTED, SMOOTH_RAMP, _title, apply_style, figsize,  # noqa: E402
                                         figure_dir, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
OUT = figure_dir("mechanics", "barrier3d") / "barrier3d_dune_growth_examples.png"
PARAMS = PROJECT_ROOT / "data" / "hatteras_init" / "Hatteras-CASCADE-parameters.yaml"
# DUNE_GROWTH_RATES in scripts/real-world-cases/Hatteras/hatteras_ms/HAT_hindcast.py
RATES = {"1996 start (1996–2009)": (0.01, 0.03), "2009 start (2009–2025)": (0.02, 0.04)}
START_HEIGHTS_M = (1.0, 5.0, 10.0)      # dune height above the berm
YEARS = 30
COLOURS = SMOOTH_RAMP[1:]               # an ordered variable: light to dark blue
# -----------------------------------------------------------------------------


def grow(h0, r, dmax, years):
    """Barrier3D's front-row dune growth, year by year, no storms."""
    h = [h0]
    for _ in range(years):
        h.append(h[-1] + r * h[-1] * (1 - h[-1] / dmax))
    return np.array(h)


def main():
    apply_style()
    p = yaml.safe_load(PARAMS.read_text())
    berm, dmax = p["BermEl"], p["Dmaxel"] - p["BermEl"]
    t = np.arange(YEARS + 1)

    fig, axes = plt.subplots(1, 2, figsize=figsize("double", 0.42), sharey=True, constrained_layout=True)
    lines, gains = [], {}
    for i, (ax, (label, (rlo, rhi))) in enumerate(zip(axes, RATES.items())):
        ax.axhline(dmax, color=INK_MUTED, lw=0.8, ls="--")
        ax.text(0.2, dmax + 0.15, f"Dmax {dmax:.1f} m (Dmaxel {p['Dmaxel']:.1f} m NAVD88)",
                fontsize=7, color=INK_MUTED, va="bottom")
        for h0, c in zip(START_HEIGHTS_M, COLOURS):
            lo, mid, hi = (grow(h0, r, dmax, YEARS) for r in (rlo, (rlo + rhi) / 2, rhi))
            ax.fill_between(t, lo, hi, color=c, alpha=0.3, lw=0)
            ln, = ax.plot(t, mid, color=c, lw=1.4, label=f"{h0:g} m start")
            for end, r in ((lo[-1], rlo), (hi[-1], rhi)):
                gains[(label, h0, r)] = end - h0
            if i == 0:
                lines.append(ln)
        _title(ax, i, f"{label}: r {rlo}–{rhi} /yr")
        ax.set_xlabel("Years")
        ax.set_xlim(0, YEARS)
        ax.set_ylim(0, dmax + 1.0)
        ax.grid(axis="y")
        open_frame(ax)
    axes[0].set_ylabel("Dune height above berm (m)")
    fig.legend(handles=lines, loc="outside lower center", ncol=len(lines), frameon=False)
    save(fig, OUT, close=True)

    gain_txt = "; ".join(f"{label.split(' (')[0]}, {h0:g} m: +{gains[(label, h0, lo)]:.2f} to "
                         f"+{gains[(label, h0, hi)]:.2f} m"
                         for label, (lo, hi) in RATES.items() for h0 in START_HEIGHTS_M)
    record_caption(OUT, (
        f"Barrier3D's dune growth applied to a dune starting {', '.join(f'{h:g}' for h in START_HEIGHTS_M)} m "
        f"above the berm ({berm} m NAVD88), for {YEARS} years with no storms. Each year h grows by r·h·(1 − h/Dmax), "
        f"Dmax = {dmax:.1f} m above the berm; r is drawn once per 10 m dune cell between rmin and rmax and kept "
        f"for the run. Line: the mid-range r; band: rmin to rmax. Front dune row; the second row grows at a third "
        f"of this. Gain after {YEARS} years: {gain_txt}. Written by "
        f"scripts/figure_making/model/dune_growth_examples.py."))
    print(OUT)
    print(gain_txt.replace("; ", "\n"))


if __name__ == "__main__":
    main()
