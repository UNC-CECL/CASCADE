"""
duck_rslr_analysis.py
=====================
Relative sea level rise (RSLR) rate for the Duck, NC gauge
(NOAA CO-OPS station 8651370), fitted over each hindcast window.

Reads the NOAA monthly mean sea level file (seasonal cycle removed), fits an
OLS trend inside each window, and writes everything to
data/hatteras_init/3-env-forcings/2-rslr/ (2026-09-15 layout; the folder was rslr/ until 2026-09-18):

    record/   duck_8651370_meantrend.csv       the NOAA download, untouched
    fits/     duck_rslr_rates.csv              ONE ROW PER WINDOW: slope, CI,
                                               n, and the value the site
                                               config carries
              duck_rslr_timeseries_<w>.csv     the monthly record inside the
                                               window with its fitted trend
                                               and residual
    figures/  duck_rslr_full_record.png/.pdf   the record with the windows
              duck_rslr_windows.png/.pdf       one panel per window
              duck_rslr_residuals.png/.pdf     residuals per window
              CAPTIONS.md                      written by caption()

THE RATES FILE IS NEW (2026-09-15). Until then the fitted slopes existed only
as annotations burned onto the figures and as hand-typed literals in
scripts/site_layer/hatteras_site_config.py (HATTERAS_PERIODS[...]["sea_level_rise_rate"],
rounded to 0.001 m/yr). The config still carries those literals -- this
script does NOT feed the model -- but the file is the record they were read
from, and its `config_m_yr` column is the rounded value so the two can be
diffed.

STYLE. Drawn under the house standard (scripts/site_layer/hat_figure_style.py) since
2026-09-15: printed width, Arial, panel letters, nothing on the canvas that
belongs in a caption. Windows are drawn in the vintage pair -- the earlier of
two windows in red, the later in blue -- and the two pairs (1984-2004 with
2004-2024; 1996-2010 with 2010-2024) never share a panel, so the pair rule
holds everywhere. The NOAA full-record trend is the reference green.

Units: all computed rates are in metres per year [m/yr]; the CSVs carry
mm/yr beside them.

Author: Hannah Henry
Date:   5/4/2026; restyled and split into record/fits/figures 2026-09-15
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# =============================================================================
# CONFIGURATION
# =============================================================================

# The gauge record and every product of this script live in the DATA tree;
# only the script lives here (2026-09-12). Anchored on this file so it follows
# the checkout.
_SCRIPTS = Path(__file__).resolve().parents[3]          # .../scripts
import sys as _envsys
from pathlib import Path as _EnvP
_envsys.path.insert(0, str(next(_q for _q in _EnvP(__file__).resolve().parents
                                if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_env_forcings as _env  # noqa: E402
_RSLR_DATA = _env.RSLR_ROOT              # 2-rslr/ since 2026-09-18
RECORD_DIR = _env.RSLR_RECORD_DIR
FITS_DIR = _env.RSLR_FITS_DIR
FIGURES_DIR = _env.RSLR_FIGURES_DIR

DATA_FILE = RECORD_DIR / "duck_8651370_meantrend.csv"
OUTPUT_PREFIX = "duck_rslr"

# --- Windows ---
# (start_year, end_year) inclusive of both calendar years, so 1984-2004 fits
# on 21 years of monthly values. The end year is the survey that closes the
# window, which is why the windows overlap at their joints.
# Windows 3 and 4 added 2026-09-11. They OVERLAP windows 1 and 2 on purpose --
# these are four hindcast windows over one record, not a partition of it.
#
# Each window has a PAIR INDEX (0 = earlier of its pair, 1 = later) that picks
# its colour, and a PAIR that decides which panel of the full-record figure
# it is drawn on.
WINDOWS = [
    # start, end,  pair, pair_index
    (1984, 2004, 0, 0),
    (2004, 2024, 0, 1),
    (1996, 2010, 1, 0),
    (2010, 2024, 1, 1),
]

# --- Fit settings ---
MIN_MONTHS = 24
CONFIG_DECIMALS = 3     # the precision hatteras_site_config.py stores rates at

# =============================================================================
# END CONFIGURATION
# =============================================================================

# `scripts/` holds hat_figure_style; this file is three levels below it.
sys.path.insert(0, str(_SCRIPTS))
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, C_1984, C_1997, C_1984_FILL, C_1997_FILL, INK_MUTED,
    _title, apply_style, caption, figsize, open_frame, save,
)

PAIR_COLOUR = {0: C_1984, 1: C_1997}
PAIR_FILL = {0: C_1984_FILL, 1: C_1997_FILL}


def window_label(start: int, end: int) -> str:
    return f"{start}–{end}"


def window_token(start: int, end: int) -> str:
    return f"{start}_{end}"


def _fmt(rate_m_yr: float, ci_m_yr: float) -> str:
    """'+0.0040 ± 0.0006 m/yr' — four places, so nothing rounds to zero."""
    sign = "+" if rate_m_yr >= 0 else ""
    return f"{sign}{rate_m_yr:.4f} ± {ci_m_yr:.4f} m/yr"


# ---------------------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------------------

def load_noaa_meantrend(filepath: Path) -> pd.DataFrame:
    """
    Load a NOAA CO-OPS monthly mean trend CSV file.

    These files have 4 metadata header lines, a blank line, then a column
    header line, then data. Values are metres relative to the station's MSL
    datum (the file's own header says so; an earlier version of this script
    labelled the axis MLLW, which was wrong).
    """
    df = pd.read_csv(
        filepath,
        skiprows=6,                   # 4 metadata lines + blank line + header
        names=["Year", "Month", "Monthly_MSL", "Linear_Trend",
               "High_Conf", "Low_Conf"],
        usecols=range(6),
    )
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["Monthly_MSL", "Year", "Month"]).reset_index(drop=True)

    # Decimal year: year + (month - 0.5) / 12 centres each month
    df["decimal_year"] = df["Year"] + (df["Month"] - 0.5) / 12.0

    print(f"Loaded {len(df)} monthly records: "
          f"{int(df['Year'].min())} – {int(df['Year'].max())}")
    return df


# ---------------------------------------------------------------------------
# 2. FIT LINEAR TREND (OLS + confidence intervals)
# ---------------------------------------------------------------------------

def fit_linear_trend(df: pd.DataFrame, start_year: int, end_year: int) -> dict:
    """
    OLS trend on Monthly_MSL within [start_year, end_year].

    Returns slope, intercept, R², p-value, the 95% CI half-width on the slope,
    dense predicted arrays with the CI band on the MEAN, and the subset.
    """
    mask = (df["Year"] >= start_year) & (df["Year"] <= end_year)
    df_fit = df[mask].copy()
    n = len(df_fit)
    if n < MIN_MONTHS:
        raise ValueError(
            f"Only {n} monthly observations in {start_year}–{end_year}; "
            f"need at least {MIN_MONTHS} for a trend.")

    t = df_fit["decimal_year"].values
    y = df_fit["Monthly_MSL"].values
    slope, intercept, r_value, p_value, se_slope = stats.linregress(t, y)

    t_crit = stats.t.ppf(0.975, df=n - 2)
    ci95 = t_crit * se_slope

    t_pred = np.linspace(t.min(), t.max(), 500)
    y_pred = slope * t_pred + intercept
    t_mean = t.mean()
    ss_t = np.sum((t - t_mean) ** 2)
    se_mean = se_slope * np.sqrt(1 / n + (t_pred - t_mean) ** 2 / ss_t)

    return {
        "slope_m_yr": slope,
        "intercept": intercept,
        "r_squared": r_value ** 2,
        "p_value": p_value,
        "ci95_m_yr": ci95,
        "t_predicted": t_pred,
        "y_predicted": y_pred,
        "y_ci_upper": y_pred + t_crit * se_mean,
        "y_ci_lower": y_pred - t_crit * se_mean,
        "df_fit": df_fit,
        "n": n,
    }


def print_summary(result: dict, start: int, end: int) -> None:
    print(f"  {window_label(start, end)}:  "
          f"{_fmt(result['slope_m_yr'], result['ci95_m_yr'])}   "
          f"n={result['n']}  R²={result['r_squared']:.3f}  "
          f"p={result['p_value']:.1e}  "
          f"config={round(result['slope_m_yr'], CONFIG_DECIMALS):.3f}")


# ---------------------------------------------------------------------------
# 3. TABLES
# ---------------------------------------------------------------------------

def export_rates(windows, results) -> Path:
    """
    One row per window. `config_m_yr` is the slope at the precision the site
    config stores, so a diff against HATTERAS_PERIODS is one column.
    """
    rows = []
    for (start, end, _pair, _idx), r in zip(windows, results):
        rows.append({
            "window": window_token(start, end),
            "start_year": start,
            "end_year": end,
            "n_months": r["n"],
            "slope_m_yr": r["slope_m_yr"],
            "ci95_m_yr": r["ci95_m_yr"],
            "slope_mm_yr": r["slope_m_yr"] * 1000,
            "ci95_mm_yr": r["ci95_m_yr"] * 1000,
            "intercept_m": r["intercept"],
            "r_squared": r["r_squared"],
            "p_value": r["p_value"],
            "config_m_yr": round(r["slope_m_yr"], CONFIG_DECIMALS),
        })
    out = pd.DataFrame(rows)
    fname = FITS_DIR / f"{OUTPUT_PREFIX}_rates.csv"
    out.to_csv(fname, index=False, float_format="%.6g")
    print(f"Saved: {fname.relative_to(_RSLR_DATA)}")
    return fname


def export_timeseries(result: dict, start_year: int, end_year: int) -> Path:
    """The monthly record inside the window, its fitted trend, the residual."""
    df_fit = result["df_fit"].copy()
    t = df_fit["decimal_year"].values
    y_trend = result["slope_m_yr"] * t + result["intercept"]
    residuals = df_fit["Monthly_MSL"].values - y_trend

    out = pd.DataFrame({
        "Year": df_fit["Year"].values.astype(int),
        "Month": df_fit["Month"].values.astype(int),
        "Decimal_Year": t,
        "Monthly_MSL_m": df_fit["Monthly_MSL"].values,
        "Monthly_MSL_mm": df_fit["Monthly_MSL"].values * 1000,
        "Trend_m": y_trend,
        "Trend_mm": y_trend * 1000,
        "Residual_m": residuals,
        "Residual_mm": residuals * 1000,
    })
    fname = FITS_DIR / f"{OUTPUT_PREFIX}_timeseries_{window_token(start_year, end_year)}.csv"
    out.to_csv(fname, index=False, float_format="%.5f")
    print(f"Saved: {fname.relative_to(_RSLR_DATA)}")
    return fname


# ---------------------------------------------------------------------------
# 4. FIGURES
# ---------------------------------------------------------------------------

def _draw_window(ax, result, colour, fill, label):
    """Trend line and its 95% band for one window, on top of whatever record
    the caller has already drawn."""
    ax.fill_between(result["t_predicted"], result["y_ci_lower"],
                    result["y_ci_upper"], color=fill, alpha=0.6, lw=0,
                    zorder=3)
    ax.plot(result["t_predicted"], result["y_predicted"], color=colour,
            lw=1.2, zorder=4, label=label)


def plot_full_record(df, windows, results) -> Path:
    """
    The whole record, twice: (a) with the 1984-2004 / 2004-2024 pair,
    (b) with the 1996-2010 / 2010-2024 pair. The two pairs overlap in time,
    so on one panel the four bands would sit on top of each other.
    """
    n_pairs = 1 + max(w[2] for w in windows)
    fig, axes = plt.subplots(n_pairs, 1, figsize=figsize("double", aspect=0.36 * n_pairs),
                             sharex=True, sharey=True, constrained_layout=True)
    axes = np.atleast_1d(axes)

    noaa_mask = df["Linear_Trend"].notna()
    for p, ax in enumerate(axes):
        ax.plot(df["decimal_year"], df["Monthly_MSL"], color=C["BASE"],
                lw=0.6, alpha=0.7, zorder=1, label="Monthly mean sea level, seasonal cycle removed")
        if noaa_mask.sum() > 10:
            ax.plot(df.loc[noaa_mask, "decimal_year"],
                    df.loc[noaa_mask, "Linear_Trend"], color=C["REF"], lw=1.0,
                    ls="--", zorder=2, label="NOAA trend, full record")
        pair_windows = [(w, r) for w, r in zip(windows, results) if w[2] == p]
        for (start, end, _pair, idx), r in pair_windows:
            ax.axvspan(start, end + 1, color=PAIR_FILL[idx], alpha=0.15, lw=0, zorder=0)
            _draw_window(ax, r, PAIR_COLOUR[idx], PAIR_FILL[idx],
                         f"Trend {window_label(start, end)}")
        _title(ax, p, " and ".join(window_label(w[0], w[1]) for w, _ in pair_windows))
        ax.set_ylabel("Mean sea level (m, MSL datum)")
        open_frame(ax)
        ax.grid(True, axis="y")
    axes[-1].set_xlabel("Year")
    axes[-1].set_xlim(df["decimal_year"].min() - 0.5, df["decimal_year"].max() + 0.5)

    # One legend for the figure: the record and reference once, then the
    # window trends of each panel in order.
    handles, labels = [], []
    for ax in axes:
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l not in labels:
                handles.append(h)
                labels.append(l)
    fig.legend(handles, labels, loc="outside lower center", ncol=3, frameon=False)

    rate_text = "; ".join(
        f"{window_label(w[0], w[1])} {_fmt(r['slope_m_yr'], r['ci95_m_yr'])}"
        for w, r in zip(windows, results))
    caption(fig,
            "Monthly mean sea level at Duck, NC (NOAA CO-OPS station 8651370), "
            "seasonal cycle removed, in metres relative to the station's MSL "
            "datum, with the OLS trend fitted inside each hindcast window and its "
            "95% confidence band. The earlier window of each pair is red, the "
            "later blue; the shaded spans mark the windows. The dashed green "
            f"line is NOAA's trend over the full record. Fitted rates (± 95% CI): {rate_text}. "
            "Table: fits/duck_rslr_rates.csv.")
    out = save(fig, FIGURES_DIR / f"{OUTPUT_PREFIX}_full_record", close=True)
    print(f"Saved: {out[0].relative_to(_RSLR_DATA)}")
    return out[0]


def plot_windows(windows, results) -> Path:
    """One panel per window: the monthly record inside it, the trend and its
    band. Shared y so the slopes compare by eye."""
    n = len(windows)
    ncol = 2 if n > 1 else 1
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize("double", aspect=0.33 * nrow),
                             sharey=True, constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()

    all_y = np.concatenate([r["df_fit"]["Monthly_MSL"].values for r in results])
    pad = (all_y.max() - all_y.min()) * 0.06
    for i, (ax, (start, end, _pair, idx), r) in enumerate(zip(axes, windows, results)):
        df_fit = r["df_fit"]
        ax.plot(df_fit["decimal_year"], df_fit["Monthly_MSL"], color=C["BASE"],
                lw=0.7, zorder=1, label="Monthly mean sea level")
        _draw_window(ax, r, PAIR_COLOUR[idx], PAIR_FILL[idx], "Fitted trend")
        ax.axhline(0, color=INK_MUTED, lw=0.5, ls=":", zorder=0)
        _title(ax, i, window_label(start, end))
        ax.set_xlim(start - 0.5, end + 1.0)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=8))
        ax.set_ylim(all_y.min() - pad, all_y.max() + pad)
        open_frame(ax)
        ax.grid(True, axis="y")
        if i % ncol == 0:
            ax.set_ylabel("Mean sea level (m, MSL datum)")
        if i >= n - ncol:
            ax.set_xlabel("Year")
    for ax in axes[n:]:
        ax.set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=len(labels), frameon=False)

    stats_text = "; ".join(
        f"({chr(ord('a') + i)}) {_fmt(r['slope_m_yr'], r['ci95_m_yr'])}, "
        f"n = {r['n']}, R² = {r['r_squared']:.2f}"
        for i, r in enumerate(results))
    caption(fig,
            "Each hindcast window of the Duck record on its own axes, with the "
            "monthly mean sea level (seasonal cycle removed) and the OLS trend; "
            "the 95% confidence band on the trend is drawn but is narrower than "
            "the line over most of each window. The y axis is shared. Red is the earlier "
            f"window of a pair, blue the later. {stats_text}. "
            "Table: fits/duck_rslr_rates.csv.")
    out = save(fig, FIGURES_DIR / f"{OUTPUT_PREFIX}_windows", close=True)
    print(f"Saved: {out[0].relative_to(_RSLR_DATA)}")
    return out[0]


def plot_residuals(windows, results) -> Path:
    """Observed minus fitted trend, per window, with a 12-month running mean."""
    n = len(windows)
    ncol = 2 if n > 1 else 1
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize("double", aspect=0.28 * nrow),
                             sharey=True, constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()

    for i, (ax, (start, end, _pair, idx), r) in enumerate(zip(axes, windows, results)):
        df_fit = r["df_fit"]
        t = df_fit["decimal_year"].values
        resid = df_fit["Monthly_MSL"].values - (r["slope_m_yr"] * t + r["intercept"])
        ax.bar(t, resid, width=1 / 12, color=C["BASE_FILL"], edgecolor="none",
               zorder=1, label="Monthly residual")
        ax.axhline(0, color=INK_MUTED, lw=0.5, zorder=0)
        running = pd.Series(resid).rolling(12, center=True).mean().values
        ax.plot(t, running, color=PAIR_COLOUR[idx], lw=1.3, zorder=3,
                label="12-month running mean")
        _title(ax, i, window_label(start, end))
        ax.set_xlim(start - 0.5, end + 1.0)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=8))
        open_frame(ax)
        if i % ncol == 0:
            ax.set_ylabel("Residual (m)")
        if i >= n - ncol:
            ax.set_xlabel("Year")
    for ax in axes[n:]:
        ax.set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=len(labels), frameon=False)

    caption(fig,
            "Residuals of the monthly mean sea level from the fitted trend in "
            "each hindcast window (bars), with a 12-month centred running mean "
            "(line, coloured as the window's trend). The running mean shows the "
            "interannual variability the linear trend does not carry; the model "
            "takes the trend only.")
    out = save(fig, FIGURES_DIR / f"{OUTPUT_PREFIX}_residuals", close=True)
    print(f"Saved: {out[0].relative_to(_RSLR_DATA)}")
    return out[0]


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--no-figures", action="store_true",
                    help="fit and write the tables only")
    args = ap.parse_args(argv)

    for d in (FITS_DIR, FIGURES_DIR):
        d.mkdir(parents=True, exist_ok=True)

    df = load_noaa_meantrend(DATA_FILE)

    print("\nRSLR rate, Duck NC (station 8651370), OLS on monthly MSL:")
    results = []
    for (start, end, _pair, _idx) in WINDOWS:
        r = fit_linear_trend(df, start, end)
        results.append(r)
        print_summary(r, start, end)
    print()

    export_rates(WINDOWS, results)
    for (start, end, _pair, _idx), r in zip(WINDOWS, results):
        export_timeseries(r, start, end)

    if not args.no_figures:
        apply_style()
        plot_full_record(df, WINDOWS, results)
        plot_windows(WINDOWS, results)
        plot_residuals(WINDOWS, results)


if __name__ == "__main__":
    main()
