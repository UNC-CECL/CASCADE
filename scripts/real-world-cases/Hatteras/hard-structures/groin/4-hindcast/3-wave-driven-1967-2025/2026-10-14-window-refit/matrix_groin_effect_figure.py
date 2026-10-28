#!/usr/bin/env python3
"""
Does the groin bring the hindcast closer to the observations? Full management, edgeBE, with and without it.

    python matrix_groin_effect_figure.py   # -> figures/matrix_groin_effect_full_management_edgeBE.png

Net shoreline change over the calibration (1996-2009) and test (2009-2025) windows, seaward positive, at
every GIS domain: the CoastSat net-change target, the groin run (0.16 and 0.10, off from 2017) and its
no-groin twin. Smoothed as the matrix comparison figures are (7-domain LOWESS, GIS 1-10 raw). Left: the
whole reach, scored over GIS 2-89; right: the domains around the groin, scored over GIS 2-12.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-14
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

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in HERE.parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as REPO  # noqa: E402
for p in (_HAT_DIR, _HAT_DIR / "hatteras_ms", _HAT_DIR / "hatteras_ms" / "experiments",
          _HAT_DIR / "analyze_output" / "compare_runs" / "matrix_vs_observed"):
    sys.path.insert(0, str(p))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import site_layer.hat_scoring as common  # noqa: E402
from matrix_vs_observed import model_end_change  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hat_observed_rates import net_change_domain_csv  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
MATRIX = REPO / "output" / "raw_runs" / "matrix"
RUNS = {  # (start, end): (label, run-name stem without the groin token, groin strength)
    (1996, 2009): ("Calibration 1996-2009", "1996_2009/edgeBE/HAT_1996_2009_edgeBE_offsetmetres_road_bdm", 0.16),
    (2009, 2025): ("Test 2009-2025", "2009_2025/edgeBE/HAT_2009_2025_edgeBE_offsetmetres_road_bdm_nourish", 0.10),
}
ZOOM = (1, 15)                 # GIS domains drawn around the groin (GIS 5|6)
LOCAL = (2, 12)                # scored near the groin; GIS 1 is Cape Point, outside the interior
GROIN_FACE = 5.5
OBSERVED = dict(color="black", lw=1.6)
NOGROIN = dict(color="#7f7f7f", lw=1.4)
GROIN = dict(color="#7b3294", lw=1.4)
OUT = workspace_twin(HERE) / "figures" / "matrix_groin_effect_full_management_edgeBE"
# -----------------------------------------------------------------------------


def run_change(stem, token):
    name = Path(stem).name + "_" + token
    return common.smooth_like_target(model_end_change(MATRIX / Path(stem).parent / name))


def observed(key):
    d = pd.read_csv(net_change_domain_csv(*key))
    return d.set_index(d.columns[0])["net_change_lowess7_m"]


def score(model, obs, lo, hi):
    r = (model - obs).loc[lo:hi].dropna()
    return float(r.mean()), float(np.sqrt((r ** 2).mean()))


def main():
    fs.apply_style()
    fig, axes = plt.subplots(2, 2, figsize=fs.figsize("double", height=5.6),
                             gridspec_kw=dict(width_ratios=(3, 1.2)), constrained_layout=True)
    for row, (key, (label, stem, s)) in enumerate(RUNS.items()):
        obs, on, off = observed(key), run_change(stem, "groin"), run_change(stem, "nogroin")
        for col, ax in enumerate(axes[row]):
            lo, hi = (2, 89) if col == 0 else LOCAL
            rmse = {}
            for name, series, style in (("observed", obs, OBSERVED), ("no groin", off, NOGROIN),
                                        ("groin", on, GROIN)):
                if series is not obs:
                    bias, rmse[name] = score(series, obs, lo, hi)
                    print(f"{label}, GIS {lo}-{hi}, {name}: RMSE {rmse[name]:.1f} m, bias {bias:+.1f} m")
                ax.plot(series.index, series.values, zorder=4 if series is obs else 3, **style)
            scores = f"RMSE GIS {lo}-{hi}: no groin {rmse['no groin']:.1f} m, groin {rmse['groin']:.1f} m"
            ax.axhline(0, color=fs.INK_MUTED, lw=0.6, zorder=1)
            if col == 1:
                a, b = ZOOM
                for series, style in ((obs, OBSERVED), (off, NOGROIN), (on, GROIN)):
                    ax.plot(series.index, series.values, "o", ms=2.5, color=style["color"], zorder=4)
                ax.set_xlim(a - 0.5, b + 0.5)
                span = pd.concat([x.loc[a:b] for x in (obs, on, off)])
                pad = 0.08 * (span.max() - span.min())
                ax.set_ylim(span.min() - pad, span.max() + pad)
                ax.axvline(GROIN_FACE, color=fs.C["GROIN"], lw=0.8)
                ax.text(GROIN_FACE + 0.2, 0.97, "groin", transform=ax.get_xaxis_transform(),
                        color=fs.C["GROIN"], fontsize=6.5, va="top")
                ax.set_title(f"Around the groin, RMSE GIS {lo}-{hi}\nno groin {rmse['no groin']:.1f} m, "
                             f"groin {rmse['groin']:.1f} m", fontsize=7, loc="left")
            else:
                ax.set_xlim(0.5, 90.5)
                # GIS 1 (Cape Point) runs off the scale in the test window; the axis follows GIS 2-89
                inner = pd.concat([x.loc[2:89] for x in (obs, on, off)])
                pad = 0.08 * (inner.max() - inner.min())
                ax.set_ylim(inner.min() - pad, inner.max() + pad)
                fs.town_bands(ax)
                fs.structures(ax)
                ax.set_ylabel("Net shoreline change (m)\nseaward +")
                ax.set_title(f"{label} (groin strength {s:.2f})\n{scores}", fontsize=8, loc="left")
            if row == len(RUNS) - 1:
                ax.set_xlabel(fs.DOMAIN_AXIS_LABEL)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], label="observed (CoastSat net-change target)", **OBSERVED),
                        Line2D([], [], label="model, no groin", **NOGROIN),
                        Line2D([], [], label="model, groin", **GROIN)],
               loc="outside lower center", ncol=3, frameon=False, fontsize=7.5)
    fig.suptitle("Full management, edgeBE: does the groin bring the model closer to the observed change?",
                 fontsize=9)
    paths = fs.save(fig, OUT, dpi=300, close=True)
    print(*paths, sep="\n")


if __name__ == "__main__":
    main()
