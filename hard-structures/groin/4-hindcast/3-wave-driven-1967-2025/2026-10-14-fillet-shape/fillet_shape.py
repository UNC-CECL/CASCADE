#!/usr/bin/env python3
"""
Is the model's fillet the right shape? The 1967-1978 build-up, observed against the corrected groin.

    python fillet_shape.py   # -> fillet_shape.csv, figures/fillet_shape_1967_1978.png

The wet/dry photo lines (1974, 1978) against the 1967-start runs of ../2026-10-13-pregroin-planform/
(the "full life" runs: 1967 wet/dry shoreline at GIS 2-12, groin from 1969, storms 1968-1983 replaying
1984-1999), at no groin, 0.16 and the whole-life 0.25. Change since 1967 at GIS 2-12, seaward positive,
with each line's straight line between GIS 2 and the last covered domain taken off, so only the shape
the groins made is compared (pregroin_planform.shape_change does the same).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-14
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OBSERVED = (REPO / "hard-structures" / "groin" / "1-observations" / "analysis" / "1-aerial_photo_gap"
            / "Change_from_wetdry_1967_D2_D12.csv")      # landward +, from the 1967 wet/dry line
RUNS = REPO / "output" / "raw_runs" / "experiments" / "groin-wave-driven" / "2026-10-13-pregroin-full-life"
START = 1967                   # row k of a run's shoreline matrix is the end of 1967 + k
YEARS = (1974, 1978)
GIS = np.arange(2, 13)
GROIN_FACE = 5.5
MODELS = {  # label: (strength, style)
    "model, no groin": (0.0, dict(color="#7f7f7f", lw=1.4)),
    "model, groin 0.16": (0.16, dict(color="#c2a5cf", lw=1.4)),
    "model, groin 0.25 (whole-life)": (0.25, dict(color="#7b3294", lw=1.6)),
}
OBS_STYLE = dict(color="black", lw=1.8)
OUT = HERE / "figures" / "fillet_shape_1967_1978"
# -----------------------------------------------------------------------------


# Seaward + change since 1967 with its end-to-end line (GIS 2 to the last covered domain) taken off
def shape(landward):
    v = -pd.Series(landward, index=GIS, dtype=float)
    ok = v.dropna()
    return v - np.interp(GIS, [ok.index[0], ok.index[-1]], [ok.iloc[0], ok.iloc[-1]])


def observed(year):
    t = pd.read_csv(OBSERVED).set_index("Domain_ID")
    return shape(t[f"change_from_wetdry_1967_wetdry_{year}_m"].reindex(GIS).to_numpy())


def modelled(strength, year):
    x = np.load(next(RUNS.glob(f"periods_1969-{strength:g}/*/*/*/*/*_shoreline_matrix.npy")))
    pads = HATTERAS_DOMAINS.gis_to_pad(GIS)
    return shape(x[year - START, pads] - x[0, pads])


def main():
    rows = []
    fs.apply_style()
    fig, axes = plt.subplots(1, len(YEARS), figsize=fs.figsize("double", height=3.4), sharey=True,
                             constrained_layout=True)
    for ax, year in zip(axes, YEARS):
        obs = observed(year)
        rows.append(pd.DataFrame({"year": year, "gis": GIS, "series": "observed", "change_m": obs.values}))
        ax.axhline(0, color=fs.INK_MUTED, lw=0.6)
        ax.axvline(GROIN_FACE, color=fs.C["GROIN"], lw=0.8)
        ax.text(GROIN_FACE + 0.1, 0.03, "groin", transform=ax.get_xaxis_transform(), color=fs.C["GROIN"],
                fontsize=6.5)
        ax.text(3.5, 0.97, "downdrift (south)", transform=ax.get_xaxis_transform(), ha="center",
                fontsize=6.5, color=fs.INK_MUTED, va="top")
        ax.text(9, 0.97, "updrift (north)", transform=ax.get_xaxis_transform(), ha="center",
                fontsize=6.5, color=fs.INK_MUTED, va="top")
        for label, (s, style) in MODELS.items():
            m = modelled(s, year)
            rows.append(pd.DataFrame({"year": year, "gis": GIS, "series": label, "change_m": m.values}))
            ax.plot(GIS, m.values, "o-", ms=2.5, label=label, **style)
        ax.plot(GIS, obs.values, "o-", ms=3, label="observed (wet/dry photo)", zorder=5, **OBS_STYLE)
        ax.set_title(f"1967 to {year}", fontsize=8.5, loc="left")
        ax.set_xlim(1.5, 12.5)
        ax.set_xlabel(fs.DOMAIN_AXIS_LABEL)
    axes[0].set_ylabel("Shoreline change since 1967 (m), seaward +\n(end-to-end line removed)")
    axes[0].legend(loc="lower right", frameon=False, fontsize=6.5)
    fig.suptitle("The fillet the groins built, 1967-1978: observed against the corrected groin", fontsize=9)
    paths = fs.save(fig, OUT, dpi=300, close=True)
    table = pd.concat(rows).pivot_table(index=["year", "series"], columns="gis", values="change_m")
    table.round(1).to_csv(HERE / "fillet_shape.csv")
    pd.set_option("display.width", 200)
    print(table.round(1).to_string())
    print(*paths, sep="\n")


if __name__ == "__main__":
    main()
