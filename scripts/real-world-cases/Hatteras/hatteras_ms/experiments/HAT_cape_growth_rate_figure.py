"""
Did the faster Cape Point dune growth change anything? Cape rate 0.10-0.20 against the calibrated rate, both windows.

    python scripts/hatteras_ms/experiments/HAT_cape_growth_rate_figure.py

Reads the cape_calibrated and cape_r010_020 runs of HAT_cape_growth_rate.py
(zeroBE, full management, groin on, quintic 30-buffer closure; only GIS 1-5's
rmin/rmax differ). One row per window. Left: shoreline position change, end
minus start, both runs smoothed with the 7-domain LOWESS exactly as CoastSat is
(smooth_model), against CoastSat net change. Right: the median dune crest over
GIS 1-5 each model year, both runs, against the surveyed median crest there.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-21
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

_HERE = Path(__file__).resolve()
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(_p for _p in _HERE.parents  if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
sys.path[:0] = [str(_HAT_DIR), str(_HERE.parent)]
import HAT_cape_growth_rate as X  # noqa: E402
import HAT_dune_ceiling_uniform as U  # noqa: E402
import HAT_dune_crest_vs_surveys as V  # noqa: E402
from cascade_pipeline.coastsat_lowess import CoastSatDataset, build_coastsat_series  # noqa: E402
from cascade_pipeline.plotting.net_change_comparison import LOWESS_CONFIG, smooth_model  # noqa: E402
from cascade_pipeline.shoreline import compute_change_rate  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hat_observed_rates import net_change_dir  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS as A, HATTERAS_DOMAINS as D  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = X.EXP_DIR / "figures" / "cape_growth_rate_effect_position_change_and_crest.png"
RUNS = {"cape_calibrated": ("CASCADE, cape at the island rate", fs.INK_MUTED, (0, (3, 1.6))),
        "cape_r010_020": ("CASCADE, cape r 0.10–0.20", A.model_color, "-")}
GIS_ZOOM = (1, 30)
ISLAND_RATE = {1996: "0.01–0.03", 2009: "0.02–0.04"}
BERM_NAVD88 = 1.7
OBS_COLOUR = "#6baed6"
# -----------------------------------------------------------------------------


def observed(start, end):
    csv = net_change_dir(start, end) / f"transect_net_change_{start}_{end}.csv"
    t = pd.read_csv(csv)
    s = build_coastsat_series([CoastSatDataset(label="CoastSat", period_start=start, csv_path=str(csv),
                                               rate_col="net_change_m")],
                              active_period_start=start, lowess_config=LOWESS_CONFIG, domains=D)[0]
    w = s["windows"][0]
    return t.groupby("domain_number").net_change_m.mean(), np.asarray(w["gis_x"]), np.asarray(w["smoothed"])


def model_position(rd):
    m = np.load(next(rd.glob("*_shoreline_matrix.npy")))
    return smooth_model(compute_change_rate(m, span_years=1, flip_sign=True))[D.start_real_index:D.end_real_index]


def cape_crest(rd):
    c = U.load(rd)
    cells = np.concatenate([np.asarray(c.barrier3d[D.gis_to_pad(g)].DuneDomain, float).max(axis=2)
                            for g in X.CAPE], axis=1) * 10 + BERM_NAVD88
    return np.median(cells, axis=1)


def main():
    fs.apply_style()
    fig, axes = plt.subplots(2, 2, figsize=fs.figsize(height=6.2), gridspec_kw={"width_ratios": (1.7, 1)})
    gis = np.arange(1, 91)
    rows = []
    for row, start in zip(axes, X.STARTS):
        rd = {v: X.run_dir(v, start) for v in RUNS}
        end = int(rd["cape_r010_020"].parent.parent.name.split("_")[1])
        dmean, sx, sy = observed(start, end)

        ax = row[0]
        south = dmean.index <= LOWESS_CONFIG.skip_southern_domains
        ax.plot(dmean.index[south], dmean[south], color=OBS_COLOUR,
                lw=1.6, ls="--", label="CoastSat, domain mean (unsmoothed, GIS 1–10)")
        keep = sx > LOWESS_CONFIG.skip_southern_domains
        ax.plot(sx[keep], sy[keep], color=OBS_COLOUR, lw=2.0, label="CoastSat (7-domain LOWESS)")
        pos = {}
        for v, (lab, col, ls) in RUNS.items():
            pos[v] = model_position(rd[v])
            ax.plot(gis, pos[v], color=col, ls=ls, lw=2.0, label=f"{lab} (7-domain LOWESS)")
        ax.axhline(0, color=fs.INK_MUTED, lw=0.8, ls="--")
        ax.axvline(5.5, color="#b2182b", lw=1.0, ls=":")
        ax.text(5.6, 0.97, "groin", transform=ax.get_xaxis_transform(), color="#b2182b", fontsize=7, va="top")
        ax.axvspan(0.5, 5.5, color="0.93", zorder=0)
        ax.set_xlim(GIS_ZOOM[0] - 0.5, GIS_ZOOM[1] + 0.5)
        ax.set_ylabel(f"Position change {start}–{end} (m)")
        ax.text(0.01, 0.03, f"{start}–{end}; island r {ISLAND_RATE[start]}", transform=ax.transAxes, fontsize=8)

        ax = row[1]
        for v, (lab, col, ls) in RUNS.items():
            cr = cape_crest(rd[v])
            ax.plot(start + np.arange(cr.size), cr, color=col, ls=ls, lw=2.0)
        sv = V.surveys(start)
        sv = sv[sv.gis.isin(X.CAPE)].dropna(subset=["crest"]).groupby("survey").crest.median()
        ax.plot(sv.index, sv.values, "o", color=OBS_COLOUR, ms=5, label="Surveyed median crest")
        ax.set_ylabel("Median dune crest, GIS 1–5 (m NAVD88)")
        ax.set_xlabel("Year")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
        ax.legend(loc="upper left", fontsize=7, frameon=False)

        for v in RUNS:
            rows.append(dict(window=f"{start}-{end}", run=v,
                             cape_change_m=" ".join(f"{x:+.0f}" for x in pos[v][:5]),
                             gis6_30_mean_m=round(float(pos[v][5:30].mean()), 1),
                             gis31_90_mean_m=round(float(pos[v][30:].mean()), 1)))
        rows.append(dict(window=f"{start}-{end}", run="CoastSat domain mean",
                         cape_change_m=" ".join(f"{dmean.loc[g]:+.0f}" for g in X.CAPE)))
    axes[1, 0].set_xlabel(fs.DOMAIN_AXIS_LABEL)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, fontsize=7.5, frameon=False, bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=250, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    t = pd.DataFrame(rows)
    t.to_csv(X.EXP_DIR / "tables" / "cape_growth_rate_effect.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_colwidth", 60):
        print(t.to_string(index=False))
    print(OUT)


if __name__ == "__main__":
    main()
