"""
The observed net shoreline change per model period: CoastSat end-window mean minus start-window mean.

    python scripts/input_prep/5-scr/3-rates/coastsat/net_change/coastsat_net_change.py
    python scripts/input_prep/5-scr/3-rates/coastsat/net_change/coastsat_net_change.py --periods 1996_2009

The calibration and test target of the DEM-to-DEM plan. Windows come from
hat_observed_rates.NET_CHANGE_WINDOWS; the means from mean_shoreline/<window>/.
Per transect, per domain (raw), and per domain at 7-domain LOWESS with GIS 1-10
raw, the smoothing the model side gets. Seaward positive. Writes
5-scr/3-rates/coastsat/net_change/<start>_<end>/.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-05
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from statsmodels.nonparametric.smoothers_lowess import lowess  # noqa: E402

_REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))

from cascade_pipeline.coastsat_lowess import DEFAULT_LOWESS  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK_MUTED, _title, apply_style, compare_header, figsize,
    mark_offaxis, offaxis_clause, open_frame, record_caption, save, structures, town_bands)
from site_layer.hat_observed_rates import (  # noqa: E402
    NET_CHANGE_CENTRES, NET_CHANGE_ROOT, NET_CHANGE_WINDOWS, mean_shoreline_csv,
    mean_shoreline_label, net_change_dir, net_change_domain_csv)
from site_layer.hatteras_site_config import HATTERAS_PERIODS, run_years  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
N_DOMAINS = 90
SMOOTH_DOMAINS = max(DEFAULT_LOWESS.window_domains)   # 7, the group's range
RAW_SOUTH = DEFAULT_LOWESS.skip_southern_domains      # GIS 1-10 stay raw
POS_HALF = 100.0                                      # m, as target_comparison
C_RAW = "#92c5de"                                     # domain means pale and thick, LOWESS dark and thin
C_SMOOTH = "#2166ac"
LW_RAW, LW_SMOOTH = 2.8, 1.2
ROLE = {(1996, 2009): "Calibration", (2009, 2025): "Test"}
# -----------------------------------------------------------------------------


# Per-domain series indexed by GIS 1-90
def _gis(values):
    return pd.Series(values, index=pd.RangeIndex(1, N_DOMAINS + 1, name="gis_domain"), dtype=float)


# The same smoothing matrix_vs_observed.smoothed gives the model: LOWESS over domain means, GIS 1-10 raw
def smooth_like_model(series):
    x = series.index.to_numpy(dtype=float)
    y = series.to_numpy(dtype=float)
    ok = np.isfinite(y)
    out = pd.Series(np.nan, index=series.index)
    out[ok] = lowess(y[ok], x[ok], frac=SMOOTH_DOMAINS / len(x), return_sorted=False)
    raw = series.index <= RAW_SOUTH
    out[raw] = series[raw]
    return out


# One window's per-transect means, the included ones only
def window_means(window):
    t = pd.read_csv(mean_shoreline_csv(*window))
    t = t[t["included"].astype(bool) & t["domain_number"].between(1, N_DOMAINS)]
    return t.set_index("transect_id")


# Per transect: end mean minus start mean, for transects in both windows
def transect_net_change(period):
    start_w, end_w = NET_CHANGE_WINDOWS[period]
    s, e = window_means(start_w), window_means(end_w)
    both = s.index.intersection(e.index)
    out = pd.DataFrame({
        "domain_number": s.loc[both, "domain_number"].astype(int),
        "start_mean_chainage_m": s.loc[both, "mean_chainage_m"],
        "end_mean_chainage_m": e.loc[both, "mean_chainage_m"],
        "start_n_obs": s.loc[both, "n_obs"].astype(int),
        "end_n_obs": e.loc[both, "n_obs"].astype(int),
    })
    out["net_change_m"] = out["end_mean_chainage_m"] - out["start_mean_chainage_m"]
    out["se_net_change_m"] = np.hypot(s.loc[both, "se_chainage_m"], e.loc[both, "se_chainage_m"])
    dropped = sorted(set(s.index) ^ set(e.index))
    return out.sort_values(["domain_number"]).rename_axis("transect_id"), dropped


# Per domain: transect count, raw mean, the standard error of that mean, and the smoothed curve
def domain_net_change(tr):
    g = tr.groupby("domain_number")
    raw = _gis(g["net_change_m"].mean().reindex(range(1, N_DOMAINS + 1)).values)
    n = g.size().reindex(range(1, N_DOMAINS + 1)).fillna(0).astype(int).values
    se = _gis((g["se_net_change_m"].apply(lambda v: np.sqrt((v ** 2).sum())) / g.size())
              .reindex(range(1, N_DOMAINS + 1)).values)
    return pd.DataFrame({"n_transects": n, "net_change_m": raw, "se_net_change_m": se,
                         "net_change_lowess7_m": smooth_like_model(raw)}, index=raw.index)


# Years between the two window centres
def centre_interval(period):
    a, b = (pd.Timestamp(d) for d in NET_CHANGE_CENTRES[period])
    return (b - a).days / 365.25


# Raw domain means and the smoothed curve for one period, on a given axis
def draw(ax, i, period, dom):
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(1, N_DOMAINS)
    ax.set_ylim(-POS_HALF, POS_HALF)
    ax.grid(axis="y")
    open_frame(ax)
    town_bands(ax, label=(i == 0))
    x = dom.index
    ax.plot(x, dom["net_change_m"], color=C_RAW, lw=LW_RAW, zorder=5)
    ax.plot(x, dom["net_change_lowess7_m"], color=C_SMOOTH, lw=LW_SMOOTH, zorder=6)
    off = [("domain means", mark_offaxis(ax, x, dom["net_change_m"].values, POS_HALF, color=C_RAW)),
           ("LOWESS", mark_offaxis(ax, x, dom["net_change_lowess7_m"].values, POS_HALF,
                                   color=C_SMOOTH))]
    s, e = NET_CHANGE_WINDOWS[period]
    _title(ax, i, f"{ROLE[period]} period, {period[0]}-{period[1]}: "
                  f"mean {e[0]} to {e[1]} minus mean {s[0]} to {s[1]}")
    ax.set_ylabel("Net change (m)")
    ax.legend(handles=[
        Line2D([], [], color=C_RAW, lw=LW_RAW,
               label=f"domain mean, island mean {dom['net_change_m'].mean():+.1f} m"),
        Line2D([], [], color=C_SMOOTH, lw=LW_SMOOTH,
               label=f"7-domain LOWESS, GIS 1-{RAW_SOUTH} raw")],
        loc="lower center", fontsize=6.5, frameon=False)
    return offaxis_clause(off, POS_HALF, "m")


# The caption sentence that says what a period's two windows are
def window_sentence(period):
    s, e = NET_CHANGE_WINDOWS[period]
    c0, c1 = NET_CHANGE_CENTRES[period]
    last = HATTERAS_PERIODS[period[0]]["last_model_year"]
    tail = (" CoastSat stops on 2026-01-13, so the end mean actually spans 2025-02-17 to "
            "2026-01-13." if period == (2009, 2025) else "")
    return (f"{ROLE[period]} period {period[0]}-{period[1]}: the CoastSat mean over {e[0]} to "
            f"{e[1]} (centred on {c1}) minus the mean over {s[0]} to {s[1]} (centred on {c0}), "
            f"{centre_interval(period):.2f} yr apart; the model runs {run_years(period[0])} "
            f"years, 1 Jan {period[0]} to 1 Jan {last + 1}.{tail}")


def figure(period, dom, out):
    fig, ax = plt.subplots(figsize=figsize("double", height=3.4), constrained_layout=True)
    clause = draw(ax, 0, period, dom)
    structures(ax, label=True)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    tag = f"{period[0]}_{period[1]}"
    png = out / f"coastsat_net_change_{tag}.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        "Observed net shoreline change per GIS domain, the target the model's end-minus-start "
        f"shoreline is scored against. {window_sentence(period)} Pale line: the mean over each "
        "domain's transects (transects that pass the 10-position minimum in both windows). Dark "
        f"line: the same at {SMOOTH_DOMAINS}-domain LOWESS, GIS 1-{RAW_SOUTH} raw, which is "
        "how the model side is smoothed. Seaward positive." + clause))
    return png


def both_figure(results):
    periods = [p for p in NET_CHANGE_WINDOWS if p in results]
    fig, axes = plt.subplots(len(periods), 1, figsize=figsize("double", height=2.9 * len(periods)),
                             sharex=True, constrained_layout=True)
    axes = np.atleast_1d(axes)
    clauses = [draw(ax, i, p, results[p]) for i, (ax, p) in enumerate(zip(axes, periods))]
    for i, ax in enumerate(axes):
        structures(ax, label=(i == len(axes) - 1))
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    compare_header(fig, "CoastSat net shoreline change between DEM-centred window means")
    png = NET_CHANGE_ROOT / "coastsat_net_change_both_periods.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        "Observed net shoreline change per GIS domain for the two model periods. "
        + " ".join(f"({chr(97 + i)}) {window_sentence(p)}" for i, p in enumerate(periods))
        + f" The calibration end window and the test start window are the same mean. Pale: domain "
        f"means; dark: {SMOOTH_DOMAINS}-domain LOWESS, GIS 1-{RAW_SOUTH} raw. Seaward positive."
        + "".join(clauses)))
    return png


def write_provenance(period, tr, dom, dropped, out):
    s, e = NET_CHANGE_WINDOWS[period]
    c0, c1 = NET_CHANGE_CENTRES[period]
    built = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    interior = dom.loc[2:89]
    text = f"""# coastsat/net_change/{period[0]}_{period[1]} -- the observed net change, {ROLE[period].lower()} period

Written by `scripts/input_prep/5-scr/3-rates/coastsat/net_change/coastsat_net_change.py` on {built}.

## What this is

The **target** of the DEM-to-DEM plan. Each period gets sources and sinks fitted so that the
model's end-minus-start shoreline matches this, not the LRR. The LRR also carries storms, fills
and the groin.

| | |
|---|---|
| start window | `mean_shoreline/{mean_shoreline_label(*s)}/`, centred on {c0} |
| end window | `mean_shoreline/{mean_shoreline_label(*e)}/`, centred on {c1} |
| between centres | {centre_interval(period):.2f} yr |
| model run | {run_years(period[0])} years, 1 Jan {period[0]} to 1 Jan {HATTERAS_PERIODS[period[0]]['last_model_year'] + 1} |
| transects in both windows | {len(tr)} (dropped, in one window only: {len(dropped)}) |
| domain mean net change, GIS 1-90 | {dom['net_change_m'].mean():+.1f} m (range {dom['net_change_m'].min():+.1f} to {dom['net_change_m'].max():+.1f}) |
| interior GIS 2-89, smoothed | mean {interior['net_change_lowess7_m'].mean():+.1f} m, sd {interior['net_change_lowess7_m'].std():.1f} m |
| median domain standard error | {dom['se_net_change_m'].median():.1f} m |

## Conventions

- Seaward positive. Net change = end mean chainage minus start mean chainage, per transect,
  then the mean of a domain's transects.
- `net_change_lowess7_m` is LOWESS over the 90 domain means at frac 7/90, with GIS 1-{RAW_SOUTH}
  left raw. This is what `matrix_vs_observed.smoothed` does to the model series, so both sides
  are smoothed alike.
- `se_net_change_m` combines the two window means' standard errors (in quadrature, then over
  the domain's transects). It is sampling noise only, not tide or datum error.
- The model's interval differs from the window-centre interval: the run starts on 1 Jan of the
  DEM year and steps whole years. This is reported here, not corrected.
"""
    if period == (2009, 2025):
        text += """- The end window is 2025-08-17 +/-6 months, not +/-1 yr, because CoastSat stops on
  2026-01-13. Its data covers 2025-02-17 to 2026-01-13.
"""
    (out / "PROVENANCE.md").write_text(text, encoding="utf-8")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--periods", nargs="+", default=None,
                    help="start_end tokens, e.g. 1996_2009 (default: every NET_CHANGE_WINDOWS period)")
    a = ap.parse_args(argv)
    periods = ([tuple(int(v) for v in p.split("_")) for p in a.periods] if a.periods
               else list(NET_CHANGE_WINDOWS))
    apply_style()
    results = {}
    for period in periods:
        out = net_change_dir(*period)
        out.mkdir(parents=True, exist_ok=True)
        tr, dropped = transect_net_change(period)
        dom = domain_net_change(tr)
        tag = f"{period[0]}_{period[1]}"
        tr.round(4).to_csv(out / f"transect_net_change_{tag}.csv")
        dom.round(4).to_csv(net_change_domain_csv(*period))
        png = figure(period, dom, out)
        write_provenance(period, tr, dom, dropped, out)
        results[period] = dom
        print(f"{tag}: {len(tr)} transects, island mean {dom['net_change_m'].mean():+.1f} m, "
              f"{centre_interval(period):.2f} yr between centres; wrote {png.relative_to(_REPO)}")
    if len(results) == len(NET_CHANGE_WINDOWS):
        print(f"wrote {both_figure(results).relative_to(_REPO)}")


if __name__ == "__main__":
    main()
