#!/usr/bin/env python3
"""
Sensitivity Results vs CoastSat — Raw + Smoothing Window Comparison
====================================================================
Post-processing script for HAT_waveSensitivity_1984_2004.py comparison.

Reads the shoreline change rate CSVs saved during a sensitivity session
(no CASCADE re-run required) and plots each parameter sweep against:
  - Raw CoastSat LRR (faded period colour, ±1 std envelope)
  - LOESS-smoothed CoastSat for each window size in COMPARE_WINDOWS_DOMAINS

Outputs (written to OUTPUT_DIR)
--------------------------------
  wave_height_vs_coastsat.png
  wave_period_vs_coastsat.png
  wave_asymmetry_vs_coastsat.png
  wave_angle_high_fraction_vs_coastsat.png
  sensitivity_overview_2x2.png      <- all four parameters in one grid

Usage
-----
1. Point SESSION_DIR to the timestamped comparison folder from the sensitivity run.
2. Point COASTSAT_CSV_1984_2004 / _2004_2024 to your CoastSat CSVs.
3. python plot_sensitivity_vs_coastsat.py
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory
from statsmodels.nonparametric.smoothers_lowess import lowess
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# SECTION 1: SESSION FOLDER  <- edit this each time
# =============================================================================
# Point to the timestamped folder created by HAT_waveSensitivity_1984_2004.py.
# Structure expected:
#   SESSION_DIR/
#     wave_height/
#       HAT_1984_2004_wvSens_wave_height_1p0/
#         HAT_1984_2004_wvSens_wave_height_1p0_shoreline_change_rate.csv
#       ...
#     wave_period/ ...
#     wave_asymmetry/ ...
#     wave_angle_high_fraction/ ...

SESSION_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\output\sensitivity_analysis\HAT_1984_2004_waveSensitivity_20260508_154602"

# Period being plotted — must match the sensitivity run
START_YEAR = 1984
END_YEAR   = 2004

# =============================================================================
# SECTION 2: COASTSAT DATA
# =============================================================================

PROJECT_BASE_DIR  = r"C:\Users\hanna\PycharmProjects\CASCADE"
COASTSAT_BASE_DIR = os.path.join(
    PROJECT_BASE_DIR, "scripts", "input_prep", "CoastSat"
)

COASTSAT_CSV_1984_2004 = os.path.join(
    COASTSAT_BASE_DIR, "1984_2004", "domain_lrr_summary.csv"
)
COASTSAT_CSV_2004_2024 = os.path.join(
    COASTSAT_BASE_DIR, "2004_2024_specific_dates", "domain_lrr_summary.csv"
)

CS_DOMAIN_COL = "domain_number"
CS_LRR_COL    = "mean_lrr"
CS_STD_COL    = "std_lrr"     # set to None if your CSV has no std column

# Which CoastSat period to use as the primary reference
# (should match the sensitivity run period)
CS_ACTIVE_PERIOD = 1984   # <- 1984 or 2004

# =============================================================================
# SECTION 3: SMOOTHING WINDOWS
# =============================================================================
# Window sizes in CASCADE domains (500 m each).
# Must have the same number of entries as C_WINDOWS below.
#   5  domains  ->  2.5 km  ->  frac = 0.056
#  10  domains  ->  5.0 km  ->  frac = 0.111   <- recommended calibration window
#  15  domains  ->  7.5 km  ->  frac = 0.167

COMPARE_WINDOWS_DOMAINS = [5, 10, 15]

DOMAIN_MIN = 1
DOMAIN_MAX = 90

# =============================================================================
# SECTION 4: COLOURS AND STYLE
# =============================================================================

# CoastSat period colours (matches coastsat_smoothed_final.py)
C_CS_1984 = "#1F4E79"   # dark blue
C_CS_2004 = "#833C00"   # dark brown-red

# Smoothed window line colours (same order as COMPARE_WINDOWS_DOMAINS)
C_WINDOWS = ["#2ca02c", "#ff7f0e", "#9467bd"]   # green, orange, purple

# Model sweep lines — viridis colormap applied per parameter value
MODEL_LW    = 1.8   # model line width
MODEL_ALPHA = 0.82

# =============================================================================
# SECTION 5: GEOGRAPHIC ANNOTATION CONSTANTS
# =============================================================================

TOWN_SPANS = {
    "Buxton":      (7,  8),
    "Avon":        (21, 31),
    "Tri-Village": (68, 83),
}
VILLAGE_LINES = {"Salvo": 69, "Waves": 74, "Rodanthe": 80}
PIERS         = {"Avon Pier": 26, "Rodanthe Pier": 79}
GROINS        = {"Buxton Groin": 6}
WIMBLE_SHOALS = (60, 74)

C_TOWN_SPAN    = "#90AFC5"
C_WIMBLE       = "#E0A800"
C_VILLAGE_LINE = "0.40"
C_PIER         = "#1565C0"
C_GROIN        = "#B71C1C"

# =============================================================================
# SECTION 6: OUTPUT
# =============================================================================

OUTPUT_DIR = os.path.join(SESSION_DIR, "coastsat_comparison_figures")

# Parameter metadata — controls figure labels and value ordering
PARAM_META = {
    "wave_height": {
        "label": "Wave Height",
        "units": "m",
        "axis_label": "Hs (m)",
    },
    "wave_period": {
        "label": "Wave Period",
        "units": "s",
        "axis_label": "Tp (s)",
    },
    "wave_asymmetry": {
        "label": "Wave Asymmetry",
        "units": "",
        "axis_label": "Asymmetry",
    },
    "wave_angle_high_fraction": {
        "label": "Wave Angle High Fraction",
        "units": "",
        "axis_label": "Angle high frac",
    },
}

# =============================================================================
# HELPER: DISCOVER SENSITIVITY RUNS
# =============================================================================

def discover_param_runs(session_dir, param_name):
    """
    Scan session_dir/{param_name}/ for per-value run subfolders and return
    {param_value (float) -> rate_csv_path (str)}, sorted by value.

    Folder naming convention (from HAT_waveSensitivity_1984_2004.py):
      HAT_{YYYY}_{YYYY}_wvSens_{param_name}_{value_str}/
    where value_str uses 'p' in place of '.' (e.g. 1p2 -> 1.2).
    """
    param_dir = os.path.join(session_dir, param_name)
    if not os.path.isdir(param_dir):
        return {}

    marker = f"_wvSens_{param_name}_"
    runs   = {}

    for folder_name in sorted(os.listdir(param_dir)):
        full_path = os.path.join(param_dir, folder_name)
        if not os.path.isdir(full_path):
            continue
        if marker not in folder_name:
            continue

        val_str = folder_name.split(marker)[-1]
        try:
            val = float(val_str.replace("p", "."))
        except ValueError:
            continue

        csvs = glob.glob(os.path.join(full_path, "*_shoreline_change_rate.csv"))
        if csvs:
            runs[val] = csvs[0]

    return dict(sorted(runs.items()))


def load_rate_csv(csv_path, start_real_idx=15, end_real_idx=105):
    """
    Load a *_shoreline_change_rate.csv and return (gis_ids, rates) for
    real domains only.
    """
    df = pd.read_csv(csv_path)
    real = df.dropna(subset=["gis_domain_id"]).copy()
    real["gis_domain_id"] = real["gis_domain_id"].astype(int)
    real = real.sort_values("gis_domain_id").reset_index(drop=True)
    return real["gis_domain_id"].values, real["model_rate_m_per_yr"].values

# =============================================================================
# HELPER: COASTSAT LOADING + SMOOTHING
# =============================================================================

def load_coastsat(path, period_label):
    """
    Load a CoastSat domain_lrr_summary CSV.
    Returns a DataFrame with columns: domain, cs_lrr, cs_std
    or None if the file is missing.
    """
    if path is None or not os.path.exists(path):
        print(f"  CoastSat ({period_label}): not found — {path}")
        return None

    df = pd.read_csv(path)
    # Strip any filename prefix that got concatenated onto the header
    df.columns = [c.split("csv")[-1] if "csv" in c else c for c in df.columns]

    dom_col = CS_DOMAIN_COL
    lrr_col = CS_LRR_COL

    if dom_col not in df.columns or lrr_col not in df.columns:
        print(f"  CoastSat ({period_label}): missing columns. Found: {list(df.columns)}")
        return None

    df[dom_col] = pd.to_numeric(df[dom_col], errors="coerce")
    df[lrr_col] = pd.to_numeric(df[lrr_col], errors="coerce")
    df = df.dropna(subset=[dom_col, lrr_col])
    df[dom_col] = df[dom_col].astype(int)
    df = df[(df[dom_col] >= DOMAIN_MIN) & (df[dom_col] <= DOMAIN_MAX)]

    # std column (optional)
    if CS_STD_COL and CS_STD_COL in df.columns:
        df[CS_STD_COL] = pd.to_numeric(df[CS_STD_COL], errors="coerce").fillna(0)
        std_vals = df[CS_STD_COL].values
    else:
        std_vals = np.zeros(len(df))

    result = pd.DataFrame({
        "domain": df[dom_col].values,
        "cs_lrr": df[lrr_col].values,
        "cs_std": std_vals,
    }).sort_values("domain").reset_index(drop=True)

    print(f"  CoastSat ({period_label}): {len(result)} domains  "
          f"range {result['cs_lrr'].min():+.2f} to {result['cs_lrr'].max():+.2f} m/yr")
    return result


def domains_to_frac(n_domains):
    """Convert a window size in CASCADE domains to a LOESS frac value."""
    return n_domains / (DOMAIN_MAX - DOMAIN_MIN + 1)


def apply_loess(domains, values, frac):
    """Apply LOESS smoothing; returns smoothed values at the same positions."""
    valid = ~np.isnan(values)
    if valid.sum() < 5:
        return values.copy()
    result   = lowess(values[valid], domains[valid], frac=frac, return_sorted=True)
    smoothed = np.full_like(values, np.nan)
    smoothed[valid] = np.interp(domains[valid], result[:, 0], result[:, 1])
    return smoothed


def compute_all_smoothed(cs_df):
    """
    Compute LOESS-smoothed CoastSat LRR for every window in
    COMPARE_WINDOWS_DOMAINS.  Returns {n_domains: smoothed_array}.
    """
    d = cs_df["domain"].values.astype(float)
    v = cs_df["cs_lrr"].values
    return {
        n: apply_loess(d, v, frac=domains_to_frac(n))
        for n in COMPARE_WINDOWS_DOMAINS
    }

# =============================================================================
# GEOGRAPHIC ANNOTATION HELPERS
# =============================================================================

def add_geographic_annotations(ax):
    trans = blended_transform_factory(ax.transData, ax.transAxes)

    # 1. Wimble Shoals
    wlo, whi = WIMBLE_SHOALS
    ax.axvspan(wlo - 0.5, whi + 0.5,
               color=C_WIMBLE, alpha=0.10, zorder=0,
               hatch="///", edgecolor=C_WIMBLE, linewidth=0)
    ax.text((wlo + whi) / 2.0, 0.04, "Wimble Shoals\ninfluence",
            transform=trans, ha="center", va="bottom", fontsize=7,
            color="#7A5800", style="italic",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.80))

    # 2. Community spans
    for span_label, (d_lo, d_hi) in TOWN_SPANS.items():
        ax.axvspan(d_lo - 0.5, d_hi + 0.5,
                   color=C_TOWN_SPAN, alpha=0.14, zorder=0)
        ax.text((d_lo + d_hi) / 2.0, 0.90, span_label, transform=trans,
                ha="center", va="top", fontsize=8, color="0.25", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

    # 3. Village center lines
    for vname, dom in VILLAGE_LINES.items():
        ax.axvline(dom, color=C_VILLAGE_LINE, lw=0.9, ls="--", alpha=0.65, zorder=1)
        ax.text(dom, 0.84, vname, transform=trans,
                ha="center", va="top", fontsize=7.5, color="0.30",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.80))

    # 4. Piers
    for pname, dom in PIERS.items():
        ax.axvline(dom, color=C_PIER, lw=1.0, ls="-.", alpha=0.80, zorder=2)
        ax.text(dom, 0.76, pname, transform=trans,
                ha="center", va="top", fontsize=7, color=C_PIER, rotation=90,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.80))

    # 5. Groins
    for gname, dom in GROINS.items():
        ax.axvline(dom, color=C_GROIN, lw=1.1, ls=":", alpha=0.85, zorder=2)
        ax.text(dom, 0.76, gname, transform=trans,
                ha="center", va="top", fontsize=7, color=C_GROIN, rotation=90,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.80))


def annotation_legend_handles():
    return [
        Patch(fc=C_TOWN_SPAN, alpha=0.30, label="Community"),
        Patch(fc=C_WIMBLE, alpha=0.25, hatch="///",
              edgecolor=C_WIMBLE, linewidth=0, label="Wimble Shoals"),
        Line2D([0], [0], color=C_VILLAGE_LINE, lw=0.9, ls="--",
               label="Village center"),
        Line2D([0], [0], color=C_PIER, lw=1.0, ls="-.", label="Pier"),
        Line2D([0], [0], color=C_GROIN, lw=1.1, ls=":", label="Groin"),
    ]

# =============================================================================
# AXIS STYLING
# =============================================================================

def style_ax(ax):
    ax.set_xlim(DOMAIN_MIN - 0.5, DOMAIN_MAX + 0.5)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis="both", which="major", direction="in", length=5, labelsize=10)
    ax.tick_params(axis="both", which="minor", direction="in", length=3)
    ax.grid(True, which="major", ls=":", lw=0.6, alpha=0.4, color="gray")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(1.1)
    ax.axhline(0, color="#2c2c2c", lw=1.1, ls="--", alpha=0.55)


def add_accretion_erosion_labels(ax):
    ybot, ytop = ax.get_ylim()
    zero_frac  = (0 - ybot) / (ytop - ybot)
    ax.text(1.0, zero_frac + (1 - zero_frac) / 2, "Accretion \u25b2",
            transform=ax.transAxes, fontsize=9, color="#555555",
            ha="right", va="center", style="italic")
    ax.text(1.0, zero_frac / 2, "Erosion \u25bc",
            transform=ax.transAxes, fontsize=9, color="#555555",
            ha="right", va="center", style="italic")

# =============================================================================
# MAIN FIGURE: one parameter, all model lines + raw + all smoothing windows
# =============================================================================

def plot_param_vs_coastsat(
    param_name, runs, cs_active, cs_period_color, smoothed_dict,
    out_path, cs_ref=None, cs_ref_color=None,
):
    """
    Full annotated figure for one wave parameter sweep.

    Layout (all on one panel):
      - Raw CoastSat (faded period colour, ±1 std envelope)
      - LOESS-smoothed CoastSat for each window in COMPARE_WINDOWS_DOMAINS
      - Cascade model lines (viridis, one per parameter value)
      - Geographic annotation layer
      - Grouped legend

    Parameters
    ----------
    param_name      : str
    runs            : dict {param_value -> (gis_ids, rates)}
    cs_active       : pd.DataFrame — CoastSat for the active period
    cs_period_color : str — colour for the active CoastSat period
    smoothed_dict   : dict {n_domains -> smoothed_array} for the active period
    out_path        : str — comparison PNG path
    cs_ref          : pd.DataFrame or None — secondary CoastSat period (reference)
    cs_ref_color    : str or None
    """
    meta   = PARAM_META.get(param_name, {"label": param_name, "units": "", "axis_label": param_name})
    values = sorted(runs.keys())
    n_vals = len(values)
    cmap   = plt.get_cmap("viridis", n_vals)

    fig, ax = plt.subplots(figsize=(14, 6.0), constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # -- Geographic annotations (drawn first, underneath data) ----------------
    add_geographic_annotations(ax)

    # -- CoastSat raw (active period) -----------------------------------------
    d_cs  = cs_active["domain"].values
    lrr   = cs_active["cs_lrr"].values
    std   = cs_active["cs_std"].values

    ax.fill_between(d_cs, lrr - std, lrr + std,
                    color=cs_period_color, alpha=0.07, zorder=1)
    ax.plot(d_cs, lrr,
            color=cs_period_color, lw=1.0, alpha=0.28,
            marker="o", ms=2.0, zorder=2,
            label=f"CoastSat raw ({START_YEAR}\u2013{END_YEAR})")

    # -- CoastSat secondary period (reference, if provided) -------------------
    if cs_ref is not None and cs_ref_color is not None:
        d_ref = cs_ref["domain"].values
        ax.plot(d_ref, cs_ref["cs_lrr"].values,
                color=cs_ref_color, lw=1.0, alpha=0.20,
                ls="--", zorder=2,
                label=f"CoastSat raw (ref period, faded)")

    # -- LOESS-smoothed CoastSat (one line per window) ------------------------
    for n_dom, win_col in zip(COMPARE_WINDOWS_DOMAINS, C_WINDOWS):
        smoothed = smoothed_dict[n_dom]
        km       = n_dom * 0.5
        frac     = domains_to_frac(n_dom)
        ax.plot(d_cs, smoothed,
                color=win_col, lw=2.2, alpha=0.90, zorder=4,
                label=f"LOESS {n_dom}-domain ({km:.1f} km,  frac={frac:.3f})")

    # -- Cascade model lines (viridis, one per parameter value) ---------------
    for k, val in enumerate(values):
        gis_ids, rates = runs[val]
        col = cmap(k / max(n_vals - 1, 1))
        lbl = f"{meta['label']} = {val}{meta['units']}"
        ax.plot(gis_ids, rates, color=col,
                lw=MODEL_LW, alpha=MODEL_ALPHA, zorder=5, label=lbl)

    # -- Axis formatting ------------------------------------------------------
    style_ax(ax)

    # Set ylim from data before placing accretion/erosion labels
    all_data = np.concatenate(
        [lrr]
        + [smoothed_dict[n] for n in COMPARE_WINDOWS_DOMAINS]
        + [r for _, r in runs.values()]
    )
    finite = all_data[np.isfinite(all_data)]
    if finite.size:
        ypad = (finite.max() - finite.min()) * 0.08
        ax.set_ylim(finite.min() - ypad, finite.max() + ypad)

    add_accretion_erosion_labels(ax)

    # -- Axis labels & title --------------------------------------------------
    ax.set_xlabel("CASCADE Model Domain (500 m alongshore)",
                  fontsize=11, fontweight="bold", labelpad=6)
    ax.set_ylabel("Shoreline Change Rate (m/yr)",
                  fontsize=11, fontweight="bold", labelpad=6)
    ax.text(0.0, 1.01, "\u2190 S  |  Cape Hatteras",
            transform=ax.transAxes, fontsize=9, color="#444444",
            ha="left", va="bottom", style="italic", clip_on=False)
    ax.text(1.0, 1.01, "Pea Island  |  N \u2192",
            transform=ax.transAxes, fontsize=9, color="#444444",
            ha="right", va="bottom", style="italic", clip_on=False)
    ax.set_title(
        f"Sensitivity to {meta['label']}  |  Hatteras Island, NC  |  "
        f"{START_YEAR}\u2013{END_YEAR}",
        fontsize=12, fontweight="bold", pad=10, color="#1a2a3a",
    )

    # -- Legend: grouped into three labelled sections -------------------------
    # Group 1: Observed (CoastSat raw + smoothed windows)
    # Group 2: Model runs (one per parameter value)
    # Group 3: Geographic annotation proxies
    #
    # A dummy handle with an empty label creates visual spacing in the legend.

    def spacer(label=""):
        return Line2D([0], [0], color="none", label=label)

    observed_handles = [
        Line2D([0], [0], color=cs_period_color, lw=1.0, alpha=0.50,
               marker="o", ms=3, label=f"CoastSat raw ({START_YEAR}\u2013{END_YEAR})"),
    ]
    for n_dom, win_col in zip(COMPARE_WINDOWS_DOMAINS, C_WINDOWS):
        km   = n_dom * 0.5
        frac = domains_to_frac(n_dom)
        observed_handles.append(
            Line2D([0], [0], color=win_col, lw=2.2,
                   label=f"LOESS {n_dom}-domain ({km:.1f} km)")
        )

    model_handles = [
        Line2D([0], [0], color=cmap(k / max(n_vals - 1, 1)), lw=MODEL_LW,
               label=f"{val}{meta['units']}")
        for k, val in enumerate(values)
    ]

    all_handles = (
        [spacer("  Observed CoastSat")]
        + observed_handles
        + [spacer("  Model  \u2014  " + meta["label"])]
        + model_handles
        + [spacer("  Geographic reference")]
        + annotation_legend_handles()
    )

    ax.legend(
        handles=all_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.01),
        fontsize=8.5,
        framealpha=0.95,
        edgecolor="#cccccc",
        ncol=5,
        handlelength=1.8,
    )

    # -- Caption --------------------------------------------------------------
    window_str = ", ".join(
        f"{n} domains ({n * 0.5:.1f} km)" for n in COMPARE_WINDOWS_DOMAINS
    )
    fig.text(
        0.012, 0.002,
        f"LOESS windows: {window_str}.  "
        f"Raw CoastSat shown faded for reference.  "
        f"Model: CASCADE  |  Observed: CoastSat LRR per 500-m domain.",
        fontsize=7.5, color="#666666", ha="left", va="bottom", style="italic",
    )

    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {os.path.basename(out_path)}")
    plt.close(fig)

# =============================================================================
# 2×2 OVERVIEW FIGURE
# =============================================================================

def plot_overview_2x2(all_runs, cs_active, cs_period_color, all_smoothed, out_path):
    """
    2x2 grid: one panel per wave parameter, compact version.
    Each panel shows: raw CoastSat (faded), all smoothing windows, all model lines.
    No geographic annotations (too cluttered at this scale) — just data.
    """
    params = [p for p in PARAM_META if p in all_runs and all_runs[p]]
    n_params = len(params)
    ncols = 2
    nrows = (n_params + 1) // 2

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(13 * ncols / 2, 5.2 * nrows),
        constrained_layout=True,
    )
    fig.patch.set_facecolor("white")
    axes_flat = np.array(axes).flatten()

    d_cs  = cs_active["domain"].values
    lrr   = cs_active["cs_lrr"].values
    std   = cs_active["cs_std"].values

    for ax_idx, param_name in enumerate(params):
        ax   = axes_flat[ax_idx]
        runs = all_runs[param_name]
        meta = PARAM_META[param_name]
        values = sorted(runs.keys())
        n_vals = len(values)
        cmap   = plt.get_cmap("viridis", n_vals)

        ax.set_facecolor("white")
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_linewidth(1.1)

        # Raw CoastSat
        ax.fill_between(d_cs, lrr - std, lrr + std,
                        color=cs_period_color, alpha=0.07, zorder=1)
        ax.plot(d_cs, lrr, color=cs_period_color, lw=0.9, alpha=0.25,
                marker="o", ms=1.5, zorder=2)

        # Smoothed windows
        for n_dom, win_col in zip(COMPARE_WINDOWS_DOMAINS, C_WINDOWS):
            km   = n_dom * 0.5
            frac = domains_to_frac(n_dom)
            ax.plot(d_cs, all_smoothed[n_dom],
                    color=win_col, lw=1.8, alpha=0.90, zorder=4,
                    label=f"LOESS {n_dom}-dom ({km:.1f} km)")

        # Model lines
        for k, val in enumerate(values):
            gis_ids, rates = runs[val]
            col = cmap(k / max(n_vals - 1, 1))
            ax.plot(gis_ids, rates,
                    color=col, lw=1.5, alpha=0.82, zorder=5,
                    label=f"{val}{meta['units']}")

        ax.axhline(0, color="#2c2c2c", lw=0.9, ls="--", alpha=0.5)

        ax.set_xlim(DOMAIN_MIN - 0.5, DOMAIN_MAX + 0.5)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
        ax.tick_params(axis="both", which="major", direction="in",
                       length=4, labelsize=9)
        ax.grid(True, which="major", ls=":", lw=0.5, alpha=0.35, color="gray")

        ax.set_xlabel("GIS Domain ID", fontsize=9, fontweight="bold")
        ax.set_ylabel("Rate (m/yr)", fontsize=9, fontweight="bold")
        ax.set_title(f"Sensitivity to {meta['label']}",
                     fontsize=10, fontweight="bold", color="#1a2a3a")

        # Compact legend: smoothing windows first, then model values
        handles = []
        # Raw proxy
        handles.append(
            Line2D([0], [0], color=cs_period_color, lw=0.9, alpha=0.5,
                   marker="o", ms=2, label="CoastSat raw")
        )
        # Smoothed window proxies
        for n_dom, win_col in zip(COMPARE_WINDOWS_DOMAINS, C_WINDOWS):
            handles.append(
                Line2D([0], [0], color=win_col, lw=1.8,
                       label=f"LOESS {n_dom}-dom")
            )
        # Model proxies
        for k, val in enumerate(values):
            col = cmap(k / max(n_vals - 1, 1))
            handles.append(
                Line2D([0], [0], color=col, lw=1.5,
                       label=f"{val}{meta['units']}")
            )
        ax.legend(handles=handles, fontsize=7, framealpha=0.9,
                  loc="lower right", ncol=2)

    # Hide unused panels
    for ax in axes_flat[n_params:]:
        ax.set_visible(False)

    fig.suptitle(
        f"Wave Climate Sensitivity vs CoastSat  |  Hatteras Island, NC  |  "
        f"{START_YEAR}\u2013{END_YEAR}",
        fontsize=13, fontweight="bold", color="#1a2a3a",
    )

    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {os.path.basename(out_path)}")
    plt.close(fig)

# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("Sensitivity vs CoastSat — Raw + Smoothing Window Comparison")
    print("=" * 70)
    print(f"Session: {SESSION_DIR}")
    print(f"Output : {OUTPUT_DIR}\n")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ── Load CoastSat ─────────────────────────────────────────────────────────
    print("Loading CoastSat data...")
    cs_1984 = load_coastsat(COASTSAT_CSV_1984_2004, "1984_2004")
    cs_2004 = load_coastsat(COASTSAT_CSV_2004_2024, "2004-2024")

    # Identify active vs reference period
    if CS_ACTIVE_PERIOD == 1984:
        cs_active     = cs_1984
        cs_ref        = cs_2004
        active_color  = C_CS_1984
        ref_color     = C_CS_2004
    else:
        cs_active     = cs_2004
        cs_ref        = cs_1984
        active_color  = C_CS_2004
        ref_color     = C_CS_1984

    if cs_active is None:
        print(f"ERROR: Active CoastSat period ({CS_ACTIVE_PERIOD}) could not be loaded.")
        return

    # ── Compute smoothed CoastSat (once, reused for every parameter figure) ──
    print("\nComputing LOESS smoothed CoastSat...")
    all_smoothed = compute_all_smoothed(cs_active)
    for n_dom in COMPARE_WINDOWS_DOMAINS:
        km   = n_dom * 0.5
        frac = domains_to_frac(n_dom)
        print(f"  {n_dom}-domain window  ({km:.1f} km,  frac={frac:.3f})  ... OK")

    # ── Discover and load sensitivity runs ───────────────────────────────────
    print(f"\nScanning session folder for sensitivity runs...")
    all_runs = {}

    for param_name in PARAM_META:
        discovered = discover_param_runs(SESSION_DIR, param_name)
        if not discovered:
            print(f"  {param_name}: no runs found — skipping")
            continue

        runs = {}
        for val, csv_path in discovered.items():
            try:
                gis_ids, rates = load_rate_csv(csv_path)
                runs[val] = (gis_ids, rates)
            except Exception as e:
                print(f"    Could not load {os.path.basename(csv_path)}: {e}")

        if runs:
            all_runs[param_name] = runs
            print(f"  {PARAM_META[param_name]['label']:35s}  "
                  f"{len(runs)} values: {sorted(runs.keys())}")

    if not all_runs:
        print("\nNo sensitivity runs found. Check SESSION_DIR.")
        return

    # ── Generate one figure per parameter ────────────────────────────────────
    print(f"\nGenerating figures...")

    for param_name, runs in all_runs.items():
        meta     = PARAM_META[param_name]
        out_path = os.path.join(OUTPUT_DIR, f"{param_name}_vs_coastsat.png")
        print(f"  {meta['label']} ...")
        plot_param_vs_coastsat(
            param_name      = param_name,
            runs            = runs,
            cs_active       = cs_active,
            cs_period_color = active_color,
            smoothed_dict   = all_smoothed,
            out_path        = out_path,
            cs_ref          = cs_ref,
            cs_ref_color    = ref_color,
        )

    # ── 2×2 overview ─────────────────────────────────────────────────────────
    if len(all_runs) > 1:
        print("  2x2 overview ...")
        plot_overview_2x2(
            all_runs        = all_runs,
            cs_active       = cs_active,
            cs_period_color = active_color,
            all_smoothed    = all_smoothed,
            out_path        = os.path.join(OUTPUT_DIR, "sensitivity_overview_2x2.png"),
        )

    print(f"\nDone.  All figures in:\n  {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
