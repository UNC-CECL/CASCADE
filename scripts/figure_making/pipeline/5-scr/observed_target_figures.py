"""
observed_target_figures.py
==============================================================================
How the observed shoreline-change target the hindcast is scored against is
built from CoastSat, for the two current windows (1996-2010, 2010-2024).

    python scripts/figure_making/pipeline/5-scr/observed_target_figures.py

Writes output/figures/3-model-inputs/5-observed-target/observed_target_<window>.png.

THE STEPS DRAWN (the producers' own functions, not re-implemented)
    1. One transect's CoastSat shoreline positions (chainage, + seaward) in the
       calendar window, and the OLS slope through them: the LRR
       (5-scr/lib/coastsat_lrr.compute_lrr; the calendar-window filter of
       coastsat_domain_lrr.py, START y0-01-01, END y1-12-31). The slope drawn
       is checked against transect_lrr_full.csv.
    2. Transects grouped into their GIS domain (transect_domain_lookup.csv)
       and averaged: the raw domain mean.
    3. LOWESS at transect resolution over a 7-domain (3.5 km) window, averaged
       back to domains, with GIS 1-10 kept as raw domain means
       (cascade_pipeline.coastsat_lowess.spliced_lowess_series, the same two
       steps hindcast.build_target_table applies to the scoring target).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "5-scr" / "lib"))

from site_layer import hat_observed_rates as obs  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, C, C_1997, INK, INK_MUTED, DOMAIN_AXIS_LABEL, figsize, figure_dir, save,
    record_caption, _title, open_frame, town_bands,
)
from cascade_pipeline.coastsat_lowess import (  # noqa: E402
    CoastSatDataset, load_transect_data, lowess_transect_values, spliced_lowess_series,
    DEFAULT_LOWESS,
)
from cascade_pipeline.domains import DEFAULT_DOMAINS  # noqa: E402
from coastsat_lrr import load_timeseries, filter_dates, compute_lrr  # noqa: E402

OUT = figure_dir("inputs", "5-observed-target")
WINDOWS = ((1996, 2010), (2010, 2024))
EXAMPLE_GIS = 45
ZOOM = (38, 52)
TARGET_WINDOW = 7   # 10 until 2026-09-28, with the runner


def timeseries_file(transect_id):
    site = "_".join(transect_id.split("_")[:3])
    return obs.COASTSAT_TIMESERIES / f"{site}_timeseries" / f"{transect_id}.csv"


def fig_target(window):
    y0, y1 = window
    tag = f"{y0}_{y1}"
    csv = obs.COASTSAT_LRR_ROOT / tag / "transect_lrr_full.csv"
    full = pd.read_csv(csv)
    ds = CoastSatDataset(label=tag, period_start=y0, csv_path=str(csv))
    dom_ids, lrr, along = load_transect_data(ds)
    target, frac = spliced_lowess_series(dom_ids, along, lrr, TARGET_WINDOW)
    smooth_t, _ = lowess_transect_values(along, lrr, TARGET_WINDOW)
    raw_means = pd.Series(lrr).groupby(dom_ids).mean()
    skip = DEFAULT_LOWESS.skip_southern_domains

    # one transect, the middle of the example domain
    ex = full[full.domain_number == EXAMPLE_GIS].sort_values("transect_id")
    row = ex.iloc[len(ex) // 2]
    ts = load_timeseries(str(timeseries_file(row.transect_id)))
    inwin = filter_dates(ts, f"{y0}-01-01", f"{y1}-12-31")
    fit = compute_lrr(inwin)
    if abs(fit["lrr_m_yr"] - row.lrr_m_yr) > 1e-3:
        raise RuntimeError(f"{row.transect_id}: recomputed LRR {fit['lrr_m_yr']} != table {row.lrr_m_yr}")
    t_all = ts.date.dt.year + ts.date.dt.dayofyear / 365.25
    t_in = inwin.date.dt.year + inwin.date.dt.dayofyear / 365.25
    x_yrs = (inwin.date - inwin.date.min()).dt.total_seconds() / 86400 / 365.25
    b, a = np.polyfit(x_yrs, inwin.chainage_m, 1)

    fig = plt.figure(figsize=figsize("double", height=7.6), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, width_ratios=[1, 1.35], height_ratios=[1, 1, 1.15])
    ax_t = fig.add_subplot(gs[0, 0])
    ax_t.plot(t_all, ts.chainage_m, ".", color="0.82", ms=2.5, label="outside the window")
    ax_t.plot(t_in, inwin.chainage_m, ".", color=C_1997, ms=3, label="in the window")
    xx = np.linspace(0, x_yrs.max(), 50)
    ax_t.plot(t_in.min() + xx, a + b * xx, color=INK, lw=1.4, label=f"OLS: {b:+.2f} m/yr")
    ax_t.axvspan(y0, y1 + 1, color="0.95", lw=0, zorder=0)
    ax_t.set_ylabel("shoreline position (m, + seaward)")
    ax_t.set_xlabel("year")
    open_frame(ax_t)
    ax_t.legend(frameon=False, fontsize=7, loc="best")
    _title(ax_t, 0, f"One transect, GIS {EXAMPLE_GIS}")

    ax_z = fig.add_subplot(gs[0, 1])
    zsel = (dom_ids >= ZOOM[0]) & (dom_ids <= ZOOM[1])
    x_dom = along / DEFAULT_DOMAINS.domain_spacing_m + DEFAULT_DOMAINS.first_gis_id - 0.5
    ax_z.plot(x_dom[zsel], lrr[zsel], "o", color="0.6", ms=3, label="transect LRR")
    for g in range(ZOOM[0], ZOOM[1] + 1):
        if g in raw_means.index:
            ax_z.hlines(raw_means[g], g - 0.5, g + 0.5, color=INK, lw=1.6)
        if g % 2:
            ax_z.axvspan(g - 0.5, g + 0.5, color="0.95", lw=0, zorder=0)
    ex_x = x_dom[(dom_ids == EXAMPLE_GIS)]
    k = int(np.argmin(np.abs(np.sort(ex_x) - np.sort(ex_x)[len(ex_x) // 2])))
    ax_z.plot(np.sort(ex_x)[k], row.lrr_m_yr, "o", mfc="none", mec=C_1997, ms=7, mew=1.2)
    ax_z.axhline(0, color=INK_MUTED, lw=0.5)
    ax_z.set_xlim(ZOOM[0] - 0.5, ZOOM[1] + 0.5)
    ax_z.set_xlabel(DOMAIN_AXIS_LABEL)
    ax_z.set_ylabel("LRR (m/yr, + accretion)")
    open_frame(ax_z)
    ax_z.legend(handles=[plt.Line2D([], [], color="0.6", marker="o", ls="", ms=3, label="transect LRR"),
                         plt.Line2D([], [], color=INK, lw=1.6, label="domain mean")],
                frameon=False, fontsize=7, loc="best")
    _title(ax_z, 1, f"Transects into domains, GIS {ZOOM[0]}-{ZOOM[1]}")

    ax = fig.add_subplot(gs[1:, :])
    ax.axhline(0, color=INK_MUTED, lw=0.5)
    ax.plot(x_dom, lrr, ".", color="0.8", ms=2.5, label="transect LRR", zorder=1)
    ax.step(raw_means.index, raw_means.values, where="mid", color="0.45", lw=0.9, label="domain mean", zorder=2)
    order = np.argsort(x_dom)
    ax.plot(x_dom[order], smooth_t[order], color=C_1997, lw=1.0, alpha=0.8,
            label=f"LOWESS, {TARGET_WINDOW} domains, at transect resolution", zorder=3)
    tn = target[target.index > skip]
    ts_ = target[target.index <= skip]
    ax.plot(tn.index, tn.values, "o-", color=INK, lw=1.6, ms=2.5, label="target: LOWESS averaged to domains",
            zorder=4)
    ax.plot(ts_.index, ts_.values, "s", color=C["ACCENT"], ms=4, label=f"target: GIS 1-{skip} raw domain means",
            zorder=5)
    ax.axvspan(0.5, skip + 0.5, color=C["ACCENT_FILL"], alpha=0.25, lw=0, zorder=0)
    lo, hi = np.nanpercentile(lrr, [0.5, 99.5])
    ax.set_ylim(min(lo, np.nanmin(target)) - 0.5, max(hi, np.nanmax(target)) + 0.5)
    ax.set_xlim(0.5, 90.5)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel(f"LRR {y0}-{y1} (m/yr, + accretion)")
    town_bands(ax)
    open_frame(ax)
    ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="lower right")
    _title(ax, 2, "The scoring target")

    out = save(fig, OUT / f"observed_target_{tag}.png")
    plt.close(fig)
    n_t = int(np.isfinite(lrr).sum())
    record_caption(out[0],
        f"How the observed shoreline-change target for {y0}-{y1} is built from CoastSat. (a) One transect "
        f"({row.transect_id}, GIS {EXAMPLE_GIS}): its satellite-derived shoreline positions, + seaward, and the "
        f"ordinary least-squares slope through those in the calendar window {y0}-01-01 to {y1}-12-31 (blue; "
        f"{fit['n_obs']} positions): the linear regression rate, LRR = {row.lrr_m_yr:+.2f} m/yr (the value in "
        "transect_lrr_full.csv, recomputed here). (b) The transects of GIS "
        f"{ZOOM[0]}-{ZOOM[1]}, each placed in its domain by transect_domain_lookup.csv, and the domain mean "
        "(black); the circled point is the transect of panel a. (c) The whole reach: every transect's LRR "
        f"(grey, {n_t} transects), the raw domain means (grey steps), a LOWESS through the transects over a "
        f"{TARGET_WINDOW}-domain (5 km) window at transect resolution (blue, lowess frac {frac:.3f}), and the "
        f"target the hindcast is scored against: that LOWESS averaged back to domains for GIS {skip + 1}-90 "
        f"(black) and the raw domain means for GIS 1-{skip} (purple squares, shaded), where boundary effects "
        "near Cape Point dominate the smoother (cascade_pipeline.coastsat_lowess.spliced_lowess_series, the "
        "steps hindcast.build_target_table applies). GIS 1 is Cape Point, GIS 90 Pea Island.")
    return out


def main():
    apply_style()
    for w in WINDOWS:
        print(fig_target(w)[0].relative_to(REPO))


if __name__ == "__main__":
    main()
