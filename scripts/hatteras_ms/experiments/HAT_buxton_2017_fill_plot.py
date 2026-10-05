"""
Before and after the 2017 Buxton fill: model rates with and without it against CoastSat, 2010-2026.

    python scripts/hatteras_ms/experiments/HAT_buxton_2017_fill_plot.py [--smoothed | --rates-only]

Reads the two runs and tables written by HAT_buxton_2017_fill.py; draws the
alongshore rates and the model-minus-observed residual over the south end.
--smoothed passes the model rates through the target's own smoothing first.
--rates-only draws the rates panel alone: CoastSat smoothed at 7 domains everywhere, the model
smoothed the same but raw at GIS 1-10.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-03
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(_HERE.parent))
import HAT_buxton_2017_fill as X  # noqa: E402
from site_layer.hat_figure_style import (C, C_1984, C_1997, DOMAIN_AXIS_LABEL, INK, INK_MUTED, _title, apply_style,  # noqa: E402
                                         figsize, open_frame, record_caption, save, structures, town_bands)
from site_layer.hat_observed_rates import lrr_csv  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS, SCORE_INTERIOR_GIS  # noqa: E402
from cascade_pipeline.coastsat_lowess import spliced_lowess_series  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
XLIM = (0.5, 30.5)
OUT = X.EXP_DIR / "figures" / "buxton_2017_fill_before_after.png"
OUT_SMOOTHED = X.EXP_DIR / "figures" / "buxton_2017_fill_before_after_smoothed.png"
OUT_RATES = X.EXP_DIR / "figures" / "buxton_2017_fill_before_after_rates_lowess7.png"
SPLICE = 10                     # GIS 1..SPLICE keep raw values, as in the target
# -----------------------------------------------------------------------------


# Model rates through the target's LOWESS: one point per domain centre, same window and splice
def smooth(series, skip=SPLICE):
    ids = series.index.to_numpy(int)
    centres = (ids - 0.5) * HATTERAS_DOMAINS.domain_spacing_m
    sm, _ = spliced_lowess_series(ids, centres, series.to_numpy(float), X.TARGET_WINDOW,
                                  skip=skip, domains=HATTERAS_DOMAINS)
    return sm.reindex(series.index)


def rmse(m, t, lo, hi):
    r = (m - t).loc[lo:hi].dropna()
    return float(np.sqrt((r ** 2).mean()))


# CoastSat through a 7-domain LOWESS everywhere; the model the same but raw at GIS 1..SPLICE
def rates_only():
    from cascade_pipeline.coastsat_lowess import CoastSatDataset, LowessConfig, build_coastsat_series
    from cascade_pipeline.hindcast import build_target_table
    apply_style()
    per = pd.read_csv(X.EXP_DIR / "tables" / "per_domain_lrr.csv", index_col=0)
    cfg = LowessConfig(window_domains=(X.TARGET_WINDOW,), skip_southern_domains=0)
    series = build_coastsat_series([CoastSatDataset(label=f"CoastSat {X.WINDOW[0]}", period_start=X.WINDOW[0],
                                                    csv_path=str(lrr_csv(*X.WINDOW)))], X.WINDOW[0], cfg,
                                   domains=HATTERAS_DOMAINS)
    obs = build_target_table(series[0], cfg, HATTERAS_DOMAINS, X.TARGET_WINDOW).set_index("gis_domain")
    obs = obs["target_lrr_m_yr"].reindex(per.index)
    mod = {v: smooth(per[f"model_{v}"]) for v in X.VARIANTS}
    lo, hi = min(X.FOOTPRINT_GIS), max(X.FOOTPRINT_GIS)
    fp_rmse = {v: rmse(mod[v], obs, lo, hi) for v in X.VARIANTS}
    int_rmse = {v: rmse(mod[v], obs, *SCORE_INTERIOR_GIS) for v in X.VARIANTS}
    print({v: (round(fp_rmse[v], 3), round(int_rmse[v], 3)) for v in X.VARIANTS})

    tr = pd.read_csv(lrr_csv(*X.WINDOW)).dropna(subset=["domain_number"])
    tr = tr[tr.domain_number.between(*XLIM)]
    g = per.index.values
    fig, a = plt.subplots(figsize=figsize("double", height=4.2))
    fig.subplots_adjust(left=0.07, right=0.99, top=0.93, bottom=0.21)
    jit = np.random.default_rng(0).uniform(-0.3, 0.3, len(tr))
    a.scatter(tr.domain_number + jit, tr.lrr_m_yr, s=3, color=C["BASE_FILL"], lw=0, zorder=1)
    a.axvspan(lo - 0.5, hi + 0.5, color=C["BASE_FILL"], alpha=0.35, lw=0, zorder=0)
    a.plot(g, obs, color="black", lw=1.8, zorder=4)
    a.plot(g, mod["without2017"], color=C_1984, lw=1.5, zorder=3)
    a.plot(g, mod["with2017"], color=C_1997, lw=1.5, zorder=5)
    a.axhline(0, color=INK_MUTED, lw=0.5, zorder=2)
    a.set_xlim(*XLIM)
    a.set_ylim(-5, 11)
    a.set_xticks(range(1, 31, 2))
    a.set_ylabel("Shoreline change rate (m/yr)")
    a.set_xlabel(DOMAIN_AXIS_LABEL)
    a.grid(axis="y")
    open_frame(a)
    a.set_title("Shoreline change rate (LRR), 2010–2026")
    town_bands(a, strip=0.06)
    structures(a)
    handles = [Line2D([], [], color="black", lw=1.8, label="Observed (CoastSat)"),
               Line2D([], [], color=C["BASE_FILL"], marker="o", ms=3, ls="none", label="Individual transects"),
               Line2D([], [], color=C_1984, lw=1.5,
                      label=f"Model without the 2017 fill (RMSE {fp_rmse['without2017']:.2f})"),
               Line2D([], [], color=C_1997, lw=1.5,
                      label=f"Model with the 2017 fill (RMSE {fp_rmse['with2017']:.2f})"),
               Patch(color=C["BASE_FILL"], alpha=0.7, label="2017 fill, GIS 6–15")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
    save(fig, OUT_RATES, close=True)
    record_caption(OUT_RATES, (
        "2010–2026 edgeBE full-management hindcast with and without the 2017–18 Buxton fill "
        "(2.6 million cubic yards, placed over GIS 6–15 and fired in 2017), southern 30 domains. "
        "Linear-regression shoreline change rate per domain. CoastSat passes through a 7-domain LOWESS over "
        "every domain, GIS 1–10 included. Both model runs pass through the same LOWESS but keep their raw "
        "values at GIS 1–10, where one value per domain is too few for the smoother to hold the steep end. "
        "Grey points are individual transects. Positive is seaward. Legend RMSE is over GIS 6–15, m/yr. "
        f"Over the island interior the RMSE is {int_rmse['with2017']:.2f} with the fill and "
        f"{int_rmse['without2017']:.2f} without."))
    print(OUT_RATES)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoothed", action="store_true")
    ap.add_argument("--rates-only", action="store_true")
    args = ap.parse_args()
    if args.rates_only:
        rates_only()
        return
    smoothed = args.smoothed
    apply_style()
    per = pd.read_csv(X.EXP_DIR / "tables" / "per_domain_lrr.csv", index_col=0)
    sc = pd.read_csv(X.EXP_DIR / "tables" / "scores.csv").set_index("variant")
    out = OUT_SMOOTHED if smoothed else OUT
    if smoothed:
        for v in X.VARIANTS:
            per[f"model_{v}"] = smooth(per[f"model_{v}"])
        rows = {v: dict(gis6_15_rmse_vs_target=rmse(per[f"model_{v}"], per.target_lrr_m_yr,
                                                    min(X.FOOTPRINT_GIS), max(X.FOOTPRINT_GIS)),
                        interior_rmse_m_yr=rmse(per[f"model_{v}"], per.target_lrr_m_yr, *SCORE_INTERIOR_GIS))
                for v in X.VARIANTS}
        sc = pd.DataFrame(rows).T
        sc.index.name = "variant"
        sc.to_csv(X.EXP_DIR / "tables" / "scores_smoothed.csv")
        print(sc.round(3))
    tr = pd.read_csv(lrr_csv(*X.WINDOW)).dropna(subset=["domain_number"])
    tr = tr[tr.domain_number.between(*XLIM)]
    g = per.index.values
    fp = (min(X.FOOTPRINT_GIS) - 0.5, max(X.FOOTPRINT_GIS) + 0.5)

    fig, (a, b) = plt.subplots(2, 1, figsize=figsize("double", height=5.2), sharex=True,
                               constrained_layout=True, gridspec_kw=dict(height_ratios=(1.4, 1)))

    # (a) rates
    jit = np.random.default_rng(0).uniform(-0.3, 0.3, len(tr))
    a.scatter(tr.domain_number + jit, tr.lrr_m_yr, s=3, color=C["BASE_FILL"], lw=0, zorder=1)
    a.plot(g, per.target_lrr_m_yr, color=INK, lw=1.6, zorder=4)
    mk = dict(marker=None) if smoothed else dict(marker="o", ms=2.5)
    a.plot(g, per.model_without2017, color=C["BASE"], lw=1.3, ls=(0, (4, 2)), zorder=3, **mk)
    a.plot(g, per.model_with2017, color=C["ACCENT"], lw=1.4, zorder=5, **mk)
    a.axhline(0, color=INK_MUTED, lw=0.5, zorder=2)
    a.set_xlim(*XLIM)
    a.set_ylim(-5, 11)
    a.axvspan(*fp, color=C["ACCENT_FILL"], alpha=0.25, lw=0, zorder=0)
    a.set_ylabel("Shoreline change rate (m/yr)")
    a.grid(axis="y")
    open_frame(a)
    _title(a, 0, "Shoreline change rate (LRR), 2010–2026")

    # (b) residual
    w = 0.38
    for off, col, key, hatch in ((-w / 2, C["BASE"], "without2017", None), (w / 2, C["ACCENT"], "with2017", None)):
        res = per[f"model_{key}"] - per.target_lrr_m_yr
        b.bar(g + off, res, width=w, color=col, lw=0, zorder=3)
    b.axhline(0, color=INK, lw=0.6, zorder=4)
    b.axvspan(*fp, color=C["ACCENT_FILL"], alpha=0.25, lw=0, zorder=0)
    b.set_ylim(-4.8, 5.2)
    b.set_ylabel("Model − observed (m/yr)")
    b.set_xlabel(DOMAIN_AXIS_LABEL)
    b.set_xticks(range(1, 31, 2))
    b.grid(axis="y")
    open_frame(b)
    _title(b, 1, "Residual against the observed rate")
    for ax in (a, b):
        town_bands(ax, strip=0.06)

    fp_rmse = {k: sc.loc[k, "gis6_15_rmse_vs_target"] for k in X.VARIANTS}
    handles = [Line2D([], [], color=INK, lw=1.6, label="Observed (CoastSat, smoothed)"),
               Line2D([], [], color=C["BASE_FILL"], marker="o", ms=3, ls="none", label="Individual transects"),
               Line2D([], [], color=C["BASE"], lw=1.3, ls=(0, (4, 2)),
                      label=f"Model without the 2017 fill (RMSE {fp_rmse['without2017']:.2f})", **mk),
               Line2D([], [], color=C["ACCENT"], lw=1.4,
                      label=f"Model with the 2017 fill (RMSE {fp_rmse['with2017']:.2f})", **mk),
               Patch(color=C["ACCENT_FILL"], alpha=0.5, label="2017 fill, GIS 6–15")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    structures(a)
    save(fig, out, close=True)

    s = sc.round(2)
    how = ("Both model curves pass through the same smoothing as the target (7-domain LOWESS on the domain "
           "values, raw at GIS 1–10). " if smoothed else "")
    record_caption(out, (how +
        "2010–2026 edgeBE full-management hindcast with and without the 2017–18 Buxton fill "
        "(2.6 million cubic yards, placed over GIS 6–15 and fired in 2017), southern 30 domains. "
        "(a) Model linear-regression shoreline change rate per domain against the CoastSat 2010–2026 target "
        "(7-domain LOWESS; raw domain means at GIS 1–10, where smoothing is suppressed); grey points are "
        "individual transects. Positive is seaward. (b) Model minus target per domain. Legend RMSE is over "
        f"GIS 6–15, m/yr. Over the island interior the RMSE is {s.loc['with2017', 'interior_rmse_m_yr']} "
        f"with the fill and {s.loc['without2017', 'interior_rmse_m_yr']} without. Domains north of GIS 17 "
        "are unchanged between the two runs."))
    print(out)


if __name__ == "__main__":
    main()
