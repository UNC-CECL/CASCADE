"""
The observed net shoreline change per model period: CoastSat end-window mean minus start-window mean.

    python scripts/input_prep/5-scr/3-rates/coastsat/net_change/coastsat_net_change.py
    python scripts/input_prep/5-scr/3-rates/coastsat/net_change/coastsat_net_change.py --periods 1996_2009

The calibration and test target of the DEM-to-DEM plan. Windows come from
hat_observed_rates.NET_CHANGE_WINDOWS; the means from mean_shoreline/<window>/.
A third, full period 1996-2025 (1996 start window to 2025 end window) is context
for the full-window run, not a target, so it stays out of NET_CHANGE_WINDOWS.
Per transect, per domain (raw), and per domain at 7-domain LOWESS with GIS 1-10
raw, the smoothing the model side gets. Seaward positive. Writes
5-scr/3-rates/coastsat/net_change/<start>_<end>/.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-12
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

# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in Path(__file__).resolve().parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as _REPO  # noqa: E402
sys.path.insert(0, str(_HAT_DIR))

from cascade_pipeline.coastsat_lowess import DEFAULT_LOWESS  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK_MUTED, _letter, apply_style, figsize, mark_offaxis, offaxis_clause,
    open_frame, record_caption, save, structures, town_bands)
from site_layer.hat_observed_rates import (  # noqa: E402
    NET_CHANGE_CENTRES, NET_CHANGE_ROOT, NET_CHANGE_WINDOWS, mean_shoreline_csv,
    mean_shoreline_label, net_change_dir, net_change_domain_csv)
from site_layer.hatteras_site_config import HATTERAS_PERIODS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
N_DOMAINS = 90
SMOOTH_DOMAINS = max(DEFAULT_LOWESS.window_domains)   # 7, the group's range
RAW_SOUTH = DEFAULT_LOWESS.skip_southern_domains      # GIS 1-10 stay raw
POS_HALF = 100.0                                      # m, as target_comparison
C_RAW = "#92c5de"                                     # domain means pale and thick, LOWESS dark and thin
C_SMOOTH = "#2166ac"
LW_RAW, LW_SMOOTH = 2.8, 1.2
ROLE = {(1996, 2009): "Calibration", (2009, 2025): "Test", (1996, 2025): "Full"}
# The full period: calibration start window to test end window; context, not graded
FULL = (1996, 2025)
FULL_LAST_MODEL_YEAR = 2025                           # HAT_full_window_1996_2025.py: 30 years, to 1 Jan 2026
WINDOWS = {**NET_CHANGE_WINDOWS,
           FULL: (NET_CHANGE_WINDOWS[(1996, 2009)][0], NET_CHANGE_WINDOWS[(2009, 2025)][1])}
CENTRES = {**NET_CHANGE_CENTRES,
           FULL: (NET_CHANGE_CENTRES[(1996, 2009)][0], NET_CHANGE_CENTRES[(2009, 2025)][1])}
DATA_END = "2026-01-13"                               # the last CoastSat position
# -----------------------------------------------------------------------------

# What every caption says the plotted quantity is
WHAT = ("Net change in shoreline position per GIS domain: for each transect, the mean CoastSat "
        "shoreline in a window around the end date minus the mean in a window around the start "
        "date, then averaged over the domain's transects. Start windows are the DEM flight "
        "midpoint ±6 months (the island offset's window); end windows are the last modelled "
        "year ±6 months (1 Jul to 30 Jun).")


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
    start_w, end_w = WINDOWS[period]
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
    a, b = (pd.Timestamp(d) for d in CENTRES[period])
    return (b - a).days / 365.25


# A window's half-width as text, from its length: +/-6 mo (1 yr long) or +/-1 yr (2 yr long)
def half_width(window):
    days = (pd.Timestamp(window[1]) - pd.Timestamp(window[0])).days
    return "±6 mo" if days < 548 else "±1 yr"


# The model run that goes with a period: (years, last model year)
def model_run(period):
    last = (FULL_LAST_MODEL_YEAR if period == FULL
            else HATTERAS_PERIODS[period[0]]["last_model_year"])
    return last - period[0] + 1, last


# Raw domain means and the smoothed curve for one period, on a given axis
def draw(ax, i, period, dom, letter=True):
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
    if letter:
        ax.set_title(_letter(i), loc="left", fontweight="bold")
    ax.set_ylabel("Net shoreline change (m)")
    for y, va, word in ((0.97, "top", "seaward ↑"), (0.03, "bottom", "landward ↓")):
        ax.text(0.995, y, word, transform=ax.transAxes, ha="right", va=va,
                fontsize=7, color=INK_MUTED, zorder=7)
    s, e = WINDOWS[period]
    c0, c1 = CENTRES[period]
    island = f"{dom['net_change_m'].mean():+.1f}".replace("-", "−")
    ax.legend(handles=[
        Line2D([], [], color=C_RAW, lw=LW_RAW, label=f"Domain mean (island mean {island} m)"),
        Line2D([], [], color=C_SMOOTH, lw=LW_SMOOTH, label=f"{SMOOTH_DOMAINS}-domain LOWESS")],
        title=(f"Mean shoreline {c1[:4]} ({half_width(e)}) − "
               f"mean shoreline {c0[:4]} ({half_width(s)})"),
        loc="lower center", fontsize=6.5, title_fontsize=6.5, frameon=False)
    return offaxis_clause(off, POS_HALF, "m")


# The caption sentence that says what a period's two windows are
def window_sentence(period):
    s, e = WINDOWS[period]
    c0, c1 = CENTRES[period]
    years, last = model_run(period)
    run = "the full-window run" if period == FULL else "the model"
    tail = (f" CoastSat stops on {DATA_END}, so the end mean actually spans {e[0]} to "
            f"{DATA_END}." if e[1] > DATA_END else "")
    return (f"{ROLE[period]} period, {period[0]}-{period[1]}: the mean CoastSat shoreline over "
            f"{e[0]} to {e[1]} ({half_width(e)} of {c1}) minus the mean over {s[0]} to {s[1]} "
            f"({half_width(s)} of {c0}), {centre_interval(period):.2f} yr apart; {run} runs "
            f"{years} years, 1 Jan {period[0]} to 1 Jan {last + 1}.{tail}")


def figure(period, dom, out):
    fig, ax = plt.subplots(figsize=figsize("double", height=3.4), constrained_layout=True)
    clause = draw(ax, 0, period, dom, letter=False)
    structures(ax, label=True)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    tag = f"{period[0]}_{period[1]}"
    png = out / f"coastsat_net_change_{tag}.png"
    save(fig, png, dpi=300, close=True)
    role = ("Context for the full-window run, not a calibration or test target. Because the 2009 "
            "mean cancels between the two model periods, a transect's 1996-2025 change is its 1996-2009 change plus its "
            "2009-2025 change." if period == FULL else
            "The target the model's end-minus-start shoreline is scored against.")
    record_caption(png, (
        f"{WHAT} {role} {window_sentence(period)} Pale line: the mean over each "
        "domain's transects (transects that pass the 10-position minimum in both windows). Dark "
        f"line: the same at {SMOOTH_DOMAINS}-domain LOWESS, GIS 1-{RAW_SOUTH} raw, which is "
        "how the model side is smoothed. The legend title names the two means differenced. "
        "Seaward positive." + clause))
    return png


def both_figure(results):
    periods = list(NET_CHANGE_WINDOWS)
    fig, axes = plt.subplots(len(periods), 1, figsize=figsize("double", height=2.9 * len(periods)),
                             sharex=True, constrained_layout=True)
    axes = np.atleast_1d(axes)
    clauses = [draw(ax, i, p, results[p]) for i, (ax, p) in enumerate(zip(axes, periods))]
    for i, ax in enumerate(axes):
        structures(ax, label=(i == len(axes) - 1))
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    png = NET_CHANGE_ROOT / "coastsat_net_change_both_periods.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        f"{WHAT} The two model periods. "
        + " ".join(f"({chr(97 + i)}) {window_sentence(p)}" for i, p in enumerate(periods))
        + f" The calibration end window and the test start window are the same mean. Pale: domain "
        f"means; dark: {SMOOTH_DOMAINS}-domain LOWESS, GIS 1-{RAW_SOUTH} raw. The legend titles "
        "name the two means differenced. Seaward positive."
        + "".join(clauses)))
    return png


def write_provenance(period, tr, dom, dropped, out):
    s, e = WINDOWS[period]
    c0, c1 = CENTRES[period]
    years, last = model_run(period)
    built = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    interior = dom.loc[2:89]
    if period == FULL:
        what = """Context, not a target: the full record from the calibration start window to the test end
window, for reading the full-window run (`scripts/real-world-cases/Hatteras/hatteras_ms/experiments/windows/HAT_full_window_1996_2025.py`)
against. Because the 2009 mean cancels between the two model periods, a transect's change here is its 1996-2009 change plus
its 2009-2025 change."""
    else:
        what = """The **target** of the DEM-to-DEM plan. Each period gets sources and sinks fitted so that the
model's end-minus-start shoreline matches this, not the LRR. The LRR also carries storms, fills
and the groin."""
    text = f"""# coastsat/net_change/{period[0]}_{period[1]} -- the observed net change, {ROLE[period].lower()} period

Written by `scripts/input_prep/5-scr/3-rates/coastsat/net_change/coastsat_net_change.py` on {built}.

## What this is

{what}

| | |
|---|---|
| start window | `mean_shoreline/{mean_shoreline_label(*s)}/`, centred on {c0} |
| end window | `mean_shoreline/{mean_shoreline_label(*e)}/`, centred on {c1} |
| between centres | {centre_interval(period):.2f} yr |
| model run | {years} years, 1 Jan {period[0]} to 1 Jan {last + 1} |
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
    if e[1] > DATA_END:
        text += f"""- The end window runs past the CoastSat record, which stops on {DATA_END}. Its data
  covers {e[0]} to {DATA_END}.
"""
    (out / "PROVENANCE.md").write_text(text, encoding="utf-8")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--periods", nargs="+", default=None,
                    help="start_end tokens, e.g. 1996_2009 (default: the two model periods and 1996_2025)")
    a = ap.parse_args(argv)
    periods = ([tuple(int(v) for v in p.split("_")) for p in a.periods] if a.periods
               else list(WINDOWS))
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
    if all(p in results for p in NET_CHANGE_WINDOWS):
        print(f"wrote {both_figure(results).relative_to(_REPO)}")


if __name__ == "__main__":
    main()
