"""
Why 0.02-0.04 is the dune growth rate for the 2009-2025 test period, one point per figure.

    python scripts/hatteras_ms/experiments/HAT_test_period_growth_figures.py

Figures for the test-period growth study (HAT_test_period_growth_rate.py, run it
first), the 2009 counterparts of the stable-dune fit figures. The 2009 rate was
fitted island-wide, so these are too: every 10 m dune cell against the 2014-2019
lidar (corrected to the 2009 lidar on NC-12). Written to
data/.../2009-start/dune-growth-rate/figures/model-fit/:
  test_growth_fit_trajectory     the adopted rate against every survey: median and low dunes
  test_growth_fit_error_by_rate  every rate tested: 2019 error and RMSE over the surveys
  test_growth_fit_spread         crest change per cell, 2009-2019: observed against three rates
The error figure recomputes the scored table's numbers from the runs and stops if they differ.
Runs are read with their own buffer count (15 a side), not today's site config.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-21
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_test_period_growth_rate as T  # noqa: E402
import HAT_dune_crest_vs_surveys as V  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402

G, U = T.G, T.U

# --- CONFIG ------------------------------------------------------------------
OUT = T.OUT_DIR / "figures" / "model-fit"
START = T.START
ADOPTED = (0.02, 0.04)
CALIBRATION = T.ADOPTED                    # 0.01-0.03, the 1996-2009 rate, run as the reference
RATES = (CALIBRATION,) + T.RATES
SPREAD_RATES = {CALIBRATION: ("Model, r 0.01-0.03 (calibration rate)", "#92c5de"),
                ADOPTED: ("Model, r 0.02-0.04 (adopted)", fs.INK),
                (0.03, 0.07): ("Model, r 0.03-0.07", "#f4a582")}
C_OBS = fs.C_1997                          # house colour of the test period
# -----------------------------------------------------------------------------


def save(fig, name, caption):
    path = OUT / f"{name}.png"
    fs.save(fig, path, close=True)
    fs.record_caption(path, caption + " Producer: HAT_test_period_growth_figures.py.")
    print(path, flush=True)


def rate_label(r):
    return f"r = {r[0]:g}" if r[0] == r[1] else f"r {r[0]:g}-{r[1]:g}"


_CRESTS = {}


# Model crest per saved state, GIS 1-90 x 50 cells, m NAVD88, with the run's own buffer count
def model_crests(r):
    if r not in _CRESTS:
        T.point(*(T.REFERENCE if r == CALIBRATION else (T.TAG, T.EXP_DIR)))
        c = U.load(G.run_dir(START, G.name(T.CEILING, r, groin=True)))
        pad = (len(c.barrier3d) - 90) // 2
        cc = [np.asarray(c.barrier3d[pad + g - 1].DuneDomain, float).max(axis=2) * 10 for g in range(1, 91)]
        _CRESTS[r] = np.stack(cc, axis=1).reshape(cc[0].shape[0], -1) + V.BERM_NAVD88
    return _CRESTS[r]


def state_index(state):
    return int(state[:4]) - V.PERIODS[START][1]


# The surveys as a cell x survey table, on the cells with a crest in every survey; and each survey's model state
def survey_wide():
    s = V.surveys(START)
    w = s.pivot_table(index=["gis", "cell"], columns="survey", values="crest").dropna()
    states = s.groupby("survey").model_state.first()
    idx = np.array([(g - 1) * 50 + c for g, c in w.index])
    return w, states, idx


def fig_trajectory():
    w, states, idx = survey_wide()
    first = V.PERIODS[START][1]
    fig, axes = plt.subplots(1, 2, figsize=fs.figsize("double", height=3.0), layout="constrained")
    for i, (ax, q, name) in enumerate(zip(axes, (50, 10), ("Median dune crest", "Low dunes: 10th-percentile crest"))):
        m = model_crests(ADOPTED)[:, idx]
        yrs = first + np.arange(m.shape[0])
        if q == 50:
            ax.fill_between(yrs, *np.percentile(m, [25, 75], axis=1), color=fs.INK, alpha=0.12, lw=0,
                            label="Model, r 0.02-0.04: middle half of cells")
        ax.plot(yrs, np.percentile(m, q, axis=1), color=fs.INK, lw=2.0, label="Model, r 0.02-0.04 (adopted)")
        for r, ls in ((CALIBRATION, ":"), ((0.03, 0.07), "--")):
            line = np.percentile(model_crests(r)[:, idx], q, axis=1)
            ax.plot(yrs, line, color=fs.INK_MUTED, lw=0.9, ls=ls)
            ax.text(yrs[-1] + 0.2, line[-1], rate_label(r), fontsize=6, va="center", color=fs.INK_MUTED)
        x = np.array([state_index(states[y]) for y in w.columns]) + first + 0.08
        v = w.to_numpy()
        med = np.percentile(v, q, axis=0)
        if q == 50:
            lo, hi = np.percentile(v, 25, axis=0), np.percentile(v, 75, axis=0)
            ax.errorbar(x, med, yerr=[med - lo, hi - med], fmt="o", ms=5, mfc="white", mec=C_OBS, mew=1.4,
                        ecolor=C_OBS, elinewidth=1, capsize=0, zorder=6,
                        label="Observed: lidar survey, median and middle half")
        else:
            ax.plot(x, med, "o", ms=5, mfc="white", mec=C_OBS, mew=1.4, zorder=6, label="Observed: lidar survey")
        ax.plot(first, np.percentile(m[0], q), "s", ms=4.5, color=C_OBS, zorder=6, label="Observed: 2009 lidar (the start)")
        ax.text(0.97, 0.05, name, transform=ax.transAxes, ha="right", fontsize=7, fontweight="bold")
        ax.set(xlabel="Year (1 January)", xlim=(first - 0.4, yrs[-1] + 1.8), xticks=range(first, yrs[-1] + 1, 4))
        ax.grid(axis="y", color=fs.GRID_C, lw=0.5)
        fs.open_frame(ax)
        fs._letter_inside(ax, i)
    axes[0].set_ylabel("Dune crest (m NAVD88)")
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, frameon=False, fontsize=6.5, loc="outside lower center", ncol=4)
    save(fig, "test_growth_fit_trajectory",
         f"The adopted test-period growth rate against the 2014-2019 lidar, island-wide ({len(idx)} dune cells with a "
         f"crest in every survey). (a) Median crest; line and band: the model with r 0.02-0.04 (median, middle half of "
         f"the cells); circles and bars: the lidar. (b) The 10th-percentile crest, the low dunes. Square: the 2009 lidar, "
         f"which is the model's start. Dotted: the calibration rate, 0.01-0.03; dashed: 0.03-0.07. 2009 start, Dmaxel "
         f"12.2 m, groin on, ends re-solved 2026-10-20. The rate is fitted to these surveys, so this window tests "
         f"overwash and shoreline, not dune growth.")


# Model minus lidar, island-wide median and 10th percentile, each survey on the cells it covers (the table's rule)
def survey_errors(r):
    mod = model_crests(r)
    rows = []
    for year, g in V.surveys(START).groupby("survey"):
        g = g.sort_values(["gis", "cell"])
        obs, m = g.crest.to_numpy(), mod[state_index(g.model_state.iloc[0])]
        ok = np.isfinite(obs)
        rows.append(dict(survey=int(year), median=np.median(m[ok]) - np.median(obs[ok]),
                         p10=np.percentile(m[ok], 10) - np.percentile(obs[ok], 10)))
    return pd.DataFrame(rows)


def fig_errors():
    table = pd.read_csv(T.OUT_DIR / "tables" / "test_period_growth_scores.csv")
    rows = []
    for r in RATES:
        e = survey_errors(r)
        key = f"{r[0]:.2f}-{r[1]:.2f}"
        t = table[table.rate.str.startswith(key)].iloc[0]
        for col, tcol in (("median", "median_err_2014_2019"), ("p10", "p10_err_2014_2019")):
            ref = np.array([float(x) for x in t[tcol].split(" / ")])
            if np.abs(e[col].round(2).to_numpy() - ref).max() > 0.011:
                raise SystemExit(f"{key} {col}: runs give {e[col].round(2).tolist()}, table has {ref.tolist()}")
        rows.append(dict(rate=r, end_median=e["median"].iloc[-1], end_p10=e.p10.iloc[-1],
                         rmse_median=np.sqrt((e["median"] ** 2).mean()), rmse_p10=np.sqrt((e.p10 ** 2).mean())))
    d = pd.DataFrame(rows)
    order = sorted(RATES, key=lambda r: (np.mean(r), r[1] - r[0]))
    y = {r: k for k, r in enumerate(order)}
    fig, axes = plt.subplots(1, 2, figsize=fs.figsize("double", height=2.8), layout="constrained", sharey=True,
                             gridspec_kw=dict(width_ratios=[1.4, 1]))
    for i, (ax, kind, xl) in enumerate(zip(axes, ("end", "rmse"), ("Model minus lidar, island-wide,\n2019 survey (m)",
                                                                   "RMSE over the 2014-2019\nsurveys (m)"))):
        ax.axhspan(y[ADOPTED] - 0.45, y[ADOPTED] + 0.45, color="0.9", zorder=0)
        for col, mk, lab in (("median", "o", "Median crest"), ("p10", "^", "Low dunes (10th percentile)")):
            ax.plot(d[f"{kind}_{col}"], d.rate.map(y), ls="none", marker=mk, ms=5,
                    mfc=fs.INK if col == "median" else "white", mec=fs.INK, mew=1.1, label=lab)
        if kind == "end":
            ax.axvline(0, color=fs.INK, lw=0.8, ls="--")
        ax.set(xlabel=xl)
        ax.grid(axis="x", color=fs.GRID_C, lw=0.5)
        fs.open_frame(ax)
        fs._letter_inside(ax, i)
    names = {ADOPTED: "  (adopted)", CALIBRATION: "  (calibration rate)"}
    axes[0].set_yticks(range(len(order)), [rate_label(r) + names.get(r, "") for r in order])
    axes[0].set_ylim(-0.6, len(order) - 0.4)
    axes[0].set_ylabel("Growth rate (rmin-rmax)")
    fig.legend(*axes[0].get_legend_handles_labels(), frameon=False, fontsize=6.5, loc="outside lower center", ncol=2)
    save(fig, "test_growth_fit_error_by_rate",
         "Every growth rate tested on the 2009 start against the 2014-2019 lidar, island-wide, ordered by mean rate. "
         "(a) Model minus lidar at the 2019 survey: median crest (filled) and the low dunes' 10th-percentile crest "
         "(open). (b) RMSE over the five surveys. 0.02-0.04 has the median closest to the lidar and the lowest median "
         "RMSE; the calibration rate falls behind, faster rates overshoot. The low dunes fit the calibration rate "
         "better (2019 error -0.13 vs +0.51 m, RMSE 0.36 vs 0.45 m): the adopted rate builds them too high. Same numbers as tables/test_period_growth_scores.csv; each survey scored on the cells it covers.")


def fig_spread():
    w, states, idx = survey_wide()
    last = w.columns.max()
    obs = w[last].to_numpy() - model_crests(ADOPTED)[0][idx]          # model year 0 is the 2009 lidar
    series = [(f"Observed: lidar, 2009 to {last}", obs, C_OBS)]
    for r, (lab, col) in SPREAD_RATES.items():
        m = model_crests(r)
        series.append((lab, m[state_index(states[last])][idx] - m[0][idx], col))
    fig, ax = plt.subplots(figsize=fs.figsize("double", height=3.0), layout="constrained")
    for k, (lab, v, col) in enumerate(series):
        yk = len(series) - 1 - k
        parts = ax.violinplot(v, positions=[yk], vert=False, widths=0.8, showextrema=False)
        for b in parts["bodies"]:
            b.set_facecolor(col)
            b.set_alpha(0.35)
            b.set_edgecolor("none")
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        ax.plot([q1, q3], [yk, yk], color=col, lw=3, solid_capstyle="butt")
        ax.plot(med, yk, marker="o", ms=4.5, mfc="white", mec=col, mew=1.3)
    ax.axvline(np.median(obs), color=C_OBS, lw=0.6, ls="--")
    ticks = [f"{lab}\nmedian {np.median(v):+.2f} m, middle half {np.subtract(*np.percentile(v, [75, 25])):.2f} m wide"
             for lab, v, _ in series]
    ax.set_yticks(range(len(series)), ticks[::-1], fontsize=6.5)
    ax.set(xlabel=f"Crest change per 10 m cell, 2009-{last} (m)", xlim=(-2.0, 3.0))
    ax.grid(axis="x", color=fs.GRID_C, lw=0.5)
    fs.open_frame(ax)
    save(fig, "test_growth_fit_spread",
         f"Crest change per 10 m cell, island-wide, from the 2009 lidar (the model's start) to the {last} survey "
         f"({len(obs)} cells with a crest in every survey). Shape: distribution; bar: middle half; circle: median; "
         f"dashed: the observed median. 0.02-0.04 comes closest to the observed median; every rate gives less spread than "
         f"the lidar. Cells that lost their dune to overwash sit in the left tail. Barrier3D draws "
         f"each cell's r at random within the range, so the model's spread is not placed where the real fast and slow "
         f"dunes are.")


def main():
    fs.apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    fig_trajectory()
    fig_errors()
    fig_spread()


if __name__ == "__main__":
    main()
