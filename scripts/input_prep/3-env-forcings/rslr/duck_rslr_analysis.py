"""
duck_rslr_analysis.py
=====================
Relative Sea Level Rise (RSLR) rate analysis for Duck, NC
NOAA CO-OPS Tide Gauge Station 8651370

Reads the NOAA monthly mean sea level (seasonal cycle removed) file,
fits a linear trend to one or more user-specified time periods, and
produces publication/presentation-quality figures for each period plus
a side-by-side comparison figure when multiple periods are provided.

Units: All computed rates are in meters per year [m/yr].
       Set DISPLAY_MM = True to show results in mm/yr instead.

Author: Hannah Henry
Date:   5/4/2026
"""

# =============================================================================
# CONFIGURATION — edit these values before running
# =============================================================================

# Path to the NOAA CO-OPS mean trend CSV file
DATA_FILE = r"/scripts/input_prep/3-env-forcings/rslr/duck_8651370_meantrend.csv"

# --- Periods of interest ---
# Add as many (start_year, end_year, label) tuples as you like.
# The label is used in figure titles, legends, and filenames.
# Years are inclusive of the full calendar year.
PERIODS = [
    (1984, 2004, "Period 1: 1984–2004"),
    (2004, 2024, "Period 2: 2004–2024"),
]

# --- Output settings ---
SAVE_FIGURES    = True         # False = display only, do not save
OUTPUT_PREFIX   = "duck_rslr"  # Prefix for all saved filenames
FIGURE_DPI      = 150          # DPI for saved figures (150 = good for slides)
FIGURE_FORMAT   = "png"        # "png", "pdf", or "svg"

# --- Units ---
DISPLAY_MM      = False         # True = mm/yr in plots; False = m/yr

# --- Color scheme (presentation-friendly) ---
# Colors for each period — extend these lists if you add more periods
PERIOD_COLORS = [
    "#C94040",   # Period 1 trend: warm red
    "#2E6FA3",   # Period 2 trend: steel blue
    "#5E9E5E",   # Period 3 trend (if added): green
    "#9B59B6",   # Period 4 trend (if added): purple
]
PERIOD_CI_COLORS = [
    "#F0A0A0",   # Period 1 CI shading
    "#A0C4E8",   # Period 2 CI shading
    "#A8D5A8",
    "#D5A8E8",
]
PERIOD_HIGHLIGHT_COLORS = [
    "#FFF0B2",   # Period 1 background highlight
    "#D6EAF8",   # Period 2 background highlight
    "#D5F5E3",
    "#E8D5F5",
]

COLOR_MSL  = "#5B8DB8"   # Monthly MSL observations
COLOR_NOAA = "#7BAE7F"   # NOAA precomputed trend (all-record)

# =============================================================================
# END CONFIGURATION
# =============================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import warnings
warnings.filterwarnings("ignore")


def _fmt(rate, ci, unit):
    """
    Format a rate ± CI string with enough decimal places to avoid
    rounding to zero.  mm/yr uses 2 dp; m/yr uses 4 dp.
    """
    dp   = 2 if unit == "mm/yr" else 4
    sign = "+" if rate >= 0 else ""
    return f"{sign}{rate:.{dp}f} ± {ci:.{dp}f} {unit}"


# ---------------------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------------------

def load_noaa_meantrend(filepath):
    """
    Load a NOAA CO-OPS monthly mean trend CSV file.

    These files have 4 metadata header lines, a blank line, then a
    column-header line, then data. Trailing commas are stripped by pandas.

    Returns a DataFrame with a decimal_year column added.
    """
    df = pd.read_csv(
        filepath,
        skiprows=6,                   # 4 metadata lines + blank line + header
        names=["Year", "Month", "Monthly_MSL", "Linear_Trend",
               "High_Conf", "Low_Conf"],
        usecols=range(6),
    )

    # Coerce all columns to numeric (drops any stray text rows)
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["Monthly_MSL", "Year", "Month"]).reset_index(drop=True)

    # Decimal year: year + (month - 0.5) / 12  — centers each month
    df["decimal_year"] = df["Year"] + (df["Month"] - 0.5) / 12.0

    print(f"Loaded {len(df)} monthly records: "
          f"{int(df['Year'].min())} – {int(df['Year'].max())}")
    return df


# ---------------------------------------------------------------------------
# 2. FIT LINEAR TREND (OLS + confidence intervals)
# ---------------------------------------------------------------------------

def fit_linear_trend(df, start_year, end_year):
    """
    Fit an OLS linear trend to Monthly_MSL within [start_year, end_year].

    Returns a dict with slope, intercept, R², p-value, 95% CI half-width,
    dense predicted arrays for plotting, and the subset DataFrame.
    """
    mask   = (df["Year"] >= start_year) & (df["Year"] <= end_year)
    df_fit = df[mask].copy()
    n      = len(df_fit)

    if n < 24:
        raise ValueError(
            f"Only {n} monthly observations in {start_year}–{end_year}. "
            "Need at least 24 for a meaningful trend estimate."
        )

    t = df_fit["decimal_year"].values
    y = df_fit["Monthly_MSL"].values

    slope, intercept, r_value, p_value, se_slope = stats.linregress(t, y)

    # 95% CI on slope: t_{0.975, n-2} × SE(slope)
    t_crit = stats.t.ppf(0.975, df=n - 2)
    ci95   = t_crit * se_slope

    # Dense arrays for smooth trend line and CI band
    t_pred     = np.linspace(t.min(), t.max(), 500)
    y_pred     = slope * t_pred + intercept
    t_mean     = t.mean()
    ss_t       = np.sum((t - t_mean) ** 2)
    se_mean    = se_slope * np.sqrt(1 / n + (t_pred - t_mean) ** 2 / ss_t)
    y_ci_upper = y_pred + t_crit * se_mean
    y_ci_lower = y_pred - t_crit * se_mean

    return {
        "slope_m_yr":   slope,
        "slope_mm_yr":  slope * 1000,
        "intercept":    intercept,
        "r_squared":    r_value ** 2,
        "p_value":      p_value,
        "ci95_m_yr":    ci95,
        "ci95_mm_yr":   ci95 * 1000,
        "t_predicted":  t_pred,
        "y_predicted":  y_pred,
        "y_ci_upper":   y_ci_upper,
        "y_ci_lower":   y_ci_lower,
        "df_fit":       df_fit,
        "n":            n,
    }


# ---------------------------------------------------------------------------
# 3. PRINT RESULTS SUMMARY
# ---------------------------------------------------------------------------

def print_summary(result, label, display_mm=True):
    rate, ci, unit = (
        (result["slope_mm_yr"], result["ci95_mm_yr"], "mm/yr")
        if display_mm else
        (result["slope_m_yr"],  result["ci95_m_yr"],  "m/yr")
    )
    print(f"\n{'='*60}")
    print(f"  RSLR RATE SUMMARY — Duck, NC (Station 8651370)")
    print(f"{'='*60}")
    print(f"  Period          :  {label}")
    print(f"  Observations (n):  {result['n']} monthly values")
    print(f"  RSLR rate       :  {rate:+.3f} ± {ci:.3f} {unit}  (95% CI)")
    print(f"  R²              :  {result['r_squared']:.4f}")
    print(f"  p-value         :  {result['p_value']:.2e}")
    print(f"{'='*60}")


# ---------------------------------------------------------------------------
# 4. FIGURE A — Full record with all periods highlighted
# ---------------------------------------------------------------------------

def plot_full_record(df, periods, results,
                     display_mm=True, save=True, prefix="duck_rslr",
                     dpi=150, fmt="png"):
    """
    Full-record time series with each analysis period highlighted and
    its fitted trend overlaid. The NOAA all-record trend is shown for
    reference.
    """
    scale  = 1000 if display_mm else 1
    unit   = "mm/yr" if display_mm else "m/yr"
    ylabel = ("Monthly Mean Sea Level (mm, relative to MLLW datum)"
              if display_mm else
              "Monthly Mean Sea Level (m, relative to MLLW datum)")

    fig, ax = plt.subplots(figsize=(14, 5))

    # All monthly observations (muted)
    ax.plot(df["decimal_year"], df["Monthly_MSL"] * scale,
            color=COLOR_MSL, alpha=0.30, linewidth=0.7, zorder=1,
            label="Monthly MSL") #(seasonal cycle removed)

    # NOAA precomputed all-record trend
    noaa_mask = df["Linear_Trend"].notna()
    if noaa_mask.sum() > 10:
        ax.plot(df.loc[noaa_mask, "decimal_year"],
                df.loc[noaa_mask, "Linear_Trend"] * scale,
                color=COLOR_NOAA, linewidth=1.5, linestyle="--", alpha=0.75,
                zorder=2, label="NOAA all-record trend")

    # Each period: highlight, brighter MSL line, CI band, trend line
    for i, ((start, end, label), result) in enumerate(zip(periods, results)):
        color_hi = PERIOD_HIGHLIGHT_COLORS[i % len(PERIOD_HIGHLIGHT_COLORS)]
        color_tr = PERIOD_COLORS[i % len(PERIOD_COLORS)]
        color_ci = PERIOD_CI_COLORS[i % len(PERIOD_CI_COLORS)]
        rate     = result["slope_mm_yr"] if display_mm else result["slope_m_yr"]
        ci_val   = result["ci95_mm_yr"]  if display_mm else result["ci95_m_yr"]

        ax.axvspan(start, end + 1, color=color_hi, alpha=0.50, zorder=0)

        mask_in = (df["Year"] >= start) & (df["Year"] <= end)
        ax.plot(df.loc[mask_in, "decimal_year"],
                df.loc[mask_in, "Monthly_MSL"] * scale,
                color=COLOR_MSL, alpha=0.85, linewidth=0.9, zorder=3)

        ax.fill_between(result["t_predicted"],
                        result["y_ci_upper"] * scale,
                        result["y_ci_lower"] * scale,
                        color=color_ci, alpha=0.50, zorder=4)

        ax.plot(result["t_predicted"], result["y_predicted"] * scale,
                color=color_tr, linewidth=2.2, zorder=5,
                label=f"{label}:  {_fmt(rate, ci_val, unit)}")

    ax.set_xlabel("Year", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title("Relative Sea Level — Duck, NC  (NOAA Station 8651370)",
                 fontsize=13, fontweight="bold")
    ax.legend(loc="upper left", fontsize=9, framealpha=0.92)
    ax.grid(True, linestyle=":", alpha=0.4)
    plt.tight_layout()

    if save:
        fname = f"{prefix}_full_record.{fmt}"
        plt.savefig(fname, dpi=dpi, bbox_inches="tight")
        print(f"Saved: {fname}")

    return fig


# ---------------------------------------------------------------------------
# 5. FIGURE B — Side-by-side comparison of all periods (shared y-axis)
# ---------------------------------------------------------------------------

def plot_comparison(periods, results,
                    display_mm=True, save=True, prefix="duck_rslr",
                    dpi=150, fmt="png"):
    """
    One subplot per period, side by side, sharing the same y-axis so
    trends are directly visually comparable across periods.
    """
    n_periods = len(periods)
    scale     = 1000 if display_mm else 1
    unit      = "mm/yr" if display_mm else "m/yr"
    ylabel    = ("Monthly Mean Sea Level (mm)"
                 if display_mm else "Monthly Mean Sea Level (m)")

    fig, axes = plt.subplots(1, n_periods,
                             figsize=(6.5 * n_periods, 5),
                             sharey=True)
    if n_periods == 1:
        axes = [axes]

    # Shared y limits across all periods
    all_y  = np.concatenate([r["df_fit"]["Monthly_MSL"].values * scale
                              for r in results])
    y_pad  = (all_y.max() - all_y.min()) * 0.08
    y_lims = (all_y.min() - y_pad, all_y.max() + y_pad)

    for i, (ax, (start, end, label), result) in enumerate(
            zip(axes, periods, results)):

        df_fit   = result["df_fit"]
        color_tr = PERIOD_COLORS[i % len(PERIOD_COLORS)]
        color_ci = PERIOD_CI_COLORS[i % len(PERIOD_CI_COLORS)]
        rate     = result["slope_mm_yr"] if display_mm else result["slope_m_yr"]
        ci_val   = result["ci95_mm_yr"]  if display_mm else result["ci95_m_yr"]

        # CI band
        ax.fill_between(result["t_predicted"],
                        result["y_ci_upper"] * scale,
                        result["y_ci_lower"] * scale,
                        color=color_ci, alpha=0.55,
                        label="95% confidence interval", zorder=2)

        # Monthly bars + line
        ax.bar(df_fit["decimal_year"], df_fit["Monthly_MSL"] * scale,
               width=0.07, color=COLOR_MSL, alpha=0.45, zorder=1)
        ax.plot(df_fit["decimal_year"], df_fit["Monthly_MSL"] * scale,
                color=COLOR_MSL, linewidth=1.0, alpha=0.8, zorder=3,
                label="Monthly MSL")

        # Trend line
        ax.plot(result["t_predicted"], result["y_predicted"] * scale,
                color=color_tr, linewidth=2.5, zorder=4,
                label=f"Trend: {_fmt(rate, ci_val, unit)}")

        # Zero reference
        ax.axhline(0, color="gray", linewidth=0.8, linestyle=":", zorder=0)

        # Stats annotation box
        p_str = (f"{result['p_value']:.2e}"
                 if result["p_value"] > 1e-10 else "< 1e-10")
        annotation = (
            f"RSLR: {_fmt(rate, ci_val, unit)}\n"
            f"R² = {result['r_squared']:.3f}\n"
            f"p = {p_str}\n"
            f"n = {result['n']}"
        )
        ax.text(0.03, 0.97, annotation,
                transform=ax.transAxes, fontsize=9.5,
                verticalalignment="top",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                          edgecolor="gray", alpha=0.88))

        ax.set_title(label, fontsize=11, fontweight="bold", pad=6)
        ax.set_xlabel("Year", fontsize=11)
        ax.set_xlim(start - 0.5, end + 1.0)
        ax.set_ylim(y_lims)
        ax.grid(True, linestyle=":", alpha=0.4)
        ax.legend(loc="lower right", fontsize=8.5, framealpha=0.9)

        if i == 0:
            ax.set_ylabel(ylabel, fontsize=11)

    fig.suptitle("Relative Sea Level Rise — Duck, NC  (Station 8651370)",
                 fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()

    if save:
        fname = f"{prefix}_comparison.{fmt}"
        plt.savefig(fname, dpi=dpi, bbox_inches="tight")
        print(f"Saved: {fname}")

    return fig


# ---------------------------------------------------------------------------
# 6. FIGURE C — Individual period close-up
# ---------------------------------------------------------------------------

def plot_period_closeup(result, start_year, end_year, label,
                        period_idx=0, display_mm=True, save=True,
                        prefix="duck_rslr", dpi=150, fmt="png"):
    """
    Close-up of a single analysis period: monthly bars, fitted trend,
    95% CI band, and annotated rate box. Designed for individual slides.
    """
    df_fit   = result["df_fit"]
    scale    = 1000 if display_mm else 1
    unit     = "mm/yr" if display_mm else "m/yr"
    ylabel   = ("Monthly Mean Sea Level (mm)"
                if display_mm else "Monthly Mean Sea Level (m)")
    rate     = result["slope_mm_yr"] if display_mm else result["slope_m_yr"]
    ci_val   = result["ci95_mm_yr"]  if display_mm else result["ci95_m_yr"]
    color_tr = PERIOD_COLORS[period_idx % len(PERIOD_COLORS)]
    color_ci = PERIOD_CI_COLORS[period_idx % len(PERIOD_CI_COLORS)]

    fig, ax = plt.subplots(figsize=(11, 5))

    ax.fill_between(result["t_predicted"],
                    result["y_ci_upper"] * scale,
                    result["y_ci_lower"] * scale,
                    color=color_ci, alpha=0.55,
                    label="95% confidence interval", zorder=2)

    ax.bar(df_fit["decimal_year"], df_fit["Monthly_MSL"] * scale,
           width=0.07, color=COLOR_MSL, alpha=0.45, zorder=1)

    ax.plot(df_fit["decimal_year"], df_fit["Monthly_MSL"] * scale,
            color=COLOR_MSL, linewidth=1.0, alpha=0.8, zorder=3,
            label="Monthly MSL") #(seasonal cycle removed)

    ax.plot(result["t_predicted"], result["y_predicted"] * scale,
            color=color_tr, linewidth=2.5, zorder=4,
            label=f"Trend: {_fmt(rate, ci_val, unit)}")

    ax.axhline(0, color="gray", linewidth=0.8, linestyle=":", zorder=0)

    p_str = (f"{result['p_value']:.2e}"
             if result["p_value"] > 1e-10 else "< 1e-10")
    annotation = (
        f"RSLR rate: {_fmt(rate, ci_val, unit)}\n"
        f"R² = {result['r_squared']:.3f}   p = {p_str}\n"
        f"n = {result['n']} monthly observations"
    )
    ax.text(0.02, 0.97, annotation,
            transform=ax.transAxes, fontsize=10,
            verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="gray", alpha=0.88))

    ax.set_xlabel("Year", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(f"RSLR — {label}  |  Duck, NC (Station 8651370)",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="lower right", fontsize=9, framealpha=0.9)
    ax.set_xlim(start_year - 0.5, end_year + 1.0)
    ax.grid(True, linestyle=":", alpha=0.4)
    plt.tight_layout()

    if save:
        fname = f"{prefix}_period_{start_year}_{end_year}.{fmt}"
        plt.savefig(fname, dpi=dpi, bbox_inches="tight")
        print(f"Saved: {fname}")

    return fig


# ---------------------------------------------------------------------------
# 7. FIGURE D — Residuals from the trend (per period)
# ---------------------------------------------------------------------------

def plot_residuals(result, start_year, end_year, label,
                   display_mm=True, save=True, prefix="duck_rslr",
                   dpi=150, fmt="png"):
    """
    Residuals (observed minus fitted trend) with a 12-month running mean.
    Useful for showing interannual variability and checking fit quality.
    """
    df_fit    = result["df_fit"].copy()
    t         = df_fit["decimal_year"].values
    y_trend   = result["slope_m_yr"] * t + result["intercept"]
    scale     = 1000 if display_mm else 1
    residuals = (df_fit["Monthly_MSL"].values - y_trend) * scale
    ylabel    = "Residual (mm)" if display_mm else "Residual (m)"

    fig, ax = plt.subplots(figsize=(11, 4))

    ax.bar(df_fit["decimal_year"], residuals,
           width=0.07, color=COLOR_MSL, alpha=0.65, zorder=2)
    ax.axhline(0, color="gray", linewidth=1.0, linestyle=":", zorder=1)

    resid_series = pd.Series(residuals, index=df_fit["decimal_year"].values)
    ax.plot(df_fit["decimal_year"],
            resid_series.rolling(12, center=True).mean(),
            color="darkorange", linewidth=1.8, zorder=3,
            label="12-month running mean")

    ax.set_xlabel("Year", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(f"MSL Residuals from Trend — {label}  |  Duck, NC",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.set_xlim(start_year - 0.5, end_year + 1.0)
    ax.grid(True, linestyle=":", alpha=0.4)
    plt.tight_layout()

    if save:
        fname = f"{prefix}_residuals_{start_year}_{end_year}.{fmt}"
        plt.savefig(fname, dpi=dpi, bbox_inches="tight")
        print(f"Saved: {fname}")

    return fig


# ---------------------------------------------------------------------------
# 8. EXPORT TIME SERIES TO CSV (per period)
# ---------------------------------------------------------------------------

def export_timeseries(result, start_year, end_year, prefix="duck_rslr"):
    """
    Export the analysis-period time series and fitted trend to CSV.
    Useful for ArcGIS, supplementary dissertation data, or further analysis.
    """
    df_fit    = result["df_fit"].copy()
    t         = df_fit["decimal_year"].values
    y_trend   = result["slope_m_yr"] * t + result["intercept"]
    residuals = df_fit["Monthly_MSL"].values - y_trend

    out = pd.DataFrame({
        "Year":           df_fit["Year"].values.astype(int),
        "Month":          df_fit["Month"].values.astype(int),
        "Decimal_Year":   t,
        "Monthly_MSL_m":  df_fit["Monthly_MSL"].values,
        "Monthly_MSL_mm": df_fit["Monthly_MSL"].values * 1000,
        "Trend_m":        y_trend,
        "Trend_mm":       y_trend * 1000,
        "Residual_m":     residuals,
        "Residual_mm":    residuals * 1000,
    })

    fname = f"{prefix}_timeseries_{start_year}_{end_year}.csv"
    out.to_csv(fname, index=False, float_format="%.5f")
    print(f"Saved: {fname}")
    return out


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    # --- Load data ---
    df = load_noaa_meantrend(DATA_FILE)

    # --- Fit trend for each period and print results ---
    results = []
    for (start, end, label) in PERIODS:
        result = fit_linear_trend(df, start, end)
        results.append(result)
        print_summary(result, label, display_mm=DISPLAY_MM)

    # --- Figure A: Full record with all periods highlighted ---
    plot_full_record(
        df, PERIODS, results,
        display_mm=DISPLAY_MM, save=SAVE_FIGURES,
        prefix=OUTPUT_PREFIX, dpi=FIGURE_DPI, fmt=FIGURE_FORMAT
    )

    # --- Figure B: Side-by-side comparison (only generated if >1 period) ---
    if len(PERIODS) > 1:
        plot_comparison(
            PERIODS, results,
            display_mm=DISPLAY_MM, save=SAVE_FIGURES,
            prefix=OUTPUT_PREFIX, dpi=FIGURE_DPI, fmt=FIGURE_FORMAT
        )

    # --- Figures C & D + CSV export: one set per period ---
    for i, ((start, end, label), result) in enumerate(zip(PERIODS, results)):

        plot_period_closeup(
            result, start, end, label,
            period_idx=i,
            display_mm=DISPLAY_MM, save=SAVE_FIGURES,
            prefix=OUTPUT_PREFIX, dpi=FIGURE_DPI, fmt=FIGURE_FORMAT
        )

        plot_residuals(
            result, start, end, label,
            display_mm=DISPLAY_MM, save=SAVE_FIGURES,
            prefix=OUTPUT_PREFIX, dpi=FIGURE_DPI, fmt=FIGURE_FORMAT
        )

        export_timeseries(result, start, end, prefix=OUTPUT_PREFIX)

    plt.show()
