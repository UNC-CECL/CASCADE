"""
CoastSat LRR Smoothing — Hatteras Island
=========================================
Always runs both transect-based and domain-averaged smoothing.

  "domain"    Smooth pre-averaged domain LRR values (original approach).
              Input : domain_lrr_summary.csv — one row per CASCADE domain.

  "transect"  Smooth individual transect LRR values first, then aggregate
              the smoothed signal back to domain resolution for CASCADE.
              Input : transect_lrr_full.csv  — one row per CoastSat transect.

Within "transect" mode, TRANSECT_X_AXIS controls what the smoother uses as x:
  "transect_id"   sequential integer derived from sort order
  "along_coast_m" cumulative along-coast distance derived from domain position

along_coast_m is derived automatically from domain number if not present in
the CSV: each domain's transects are spread evenly across its 500 m band.
Physical spacing for the LOESS frac is always estimated from along_coast_m.

Outputs — domain-space figures (both modes)
-------------------------------------------
  overview_smoothed.png              raw + LOESS overlay, both periods
  smoothed_only_comparison.png       clean version for presentations
  combined_periods.png               both periods on one panel
  smoothing_sensitivity_*.png        3-panel bandwidth sensitivity
  window_comparison.png              all window sizes overlaid
  coastsat_<mode>_smoothed_table.csv domain-level raw + smoothed values

Additional outputs — transect mode only
---------------------------------------
  transect_smoothed_overview.png     raw transect scatter + LOESS in transect space
  transect_window_comparison.png     window sensitivity in transect space
  coastsat_transect_lrr_*.csv        full transect-level table with lrr_smooth
"""

# pathlib must be imported before the CONFIG block because every path below is
# built from PROJECT_BASE_DIR at module level.
import pathlib

# ANCHORED, NOT TYPED. Every path below used to be an absolute literal: the
# output one had lost its drive ("/scripts/input_prep/...") and so wrote its
# figures to C:\scripts\ instead of into the repository, and the input ones
# still spelled the folder "input_preperation" and pointed at a CoastSat tree
# that has since moved under 5-scr. Anchoring on the pyproject.toml at the repo
# root makes all of them follow the checkout and survive this file changing
# depth.
PROJECT_BASE_DIR = next(
    q for q in pathlib.Path(__file__).resolve().parents
    if (q / "pyproject.toml").exists()
)


# ============================================================
# CONFIG
# ============================================================

# ── Smoothing x-axis ─────────────────────────────────────────
# "along_coast_m" is recommended — keeps the physical window consistent.
# "transect_id" is available but causes non-uniform spacing artefacts.
TRANSECT_X_AXIS = "along_coast_m"   # "along_coast_m" | "transect_id"

# ── Domain-mode inputs ───────────────────────────────────────
DOMAIN_CSV_1984_2004 = str(PROJECT_BASE_DIR / "scripts" / "input_prep" / "5-scr" / "CoastSat"
                        / "1984_2004" / "domain_lrr_summary.csv")
DOMAIN_CSV_2004_2024 = str(PROJECT_BASE_DIR / "scripts" / "input_prep" / "5-scr" / "CoastSat"
                        / "2004_2024" / "domain_lrr_summary.csv")
CS_DOMAIN_COL = "domain_number"
CS_LRR_COL    = "mean_lrr"
CS_STD_COL    = "std_lrr"

# ── Transect-mode inputs ─────────────────────────────────────
# Point to your transect_lrr_full.csv files for each period
TRANSECT_CSV_1984_2004 = str(PROJECT_BASE_DIR / "scripts" / "input_prep" / "5-scr" / "CoastSat"
                          / "1984_2004" / "transect_lrr_full.csv")
TRANSECT_CSV_2004_2024 = str(PROJECT_BASE_DIR / "scripts" / "input_prep" / "5-scr" / "CoastSat"
                          / "2004_2024" / "transect_lrr_full.csv")

# Column names in your transect CSV (transect_lrr_full.csv)
T_TRANSECT_ID_COL = "transect_id"    # string IDs — converted to sequential int internally
T_ALONG_COAST_COL = None             # not in CSV — derived from domain position automatically
T_DOMAIN_COL      = "domain_number"
T_LRR_COL         = "lrr_m_yr"
T_STD_COL         = "unc_m_yr"       # uncertainty column; set to None to skip

# Optional: filter to only point_in_polygon domain matches
FILTER_POINT_IN_POLYGON = False
T_MATCH_METHOD_COL      = "match_method"

# ── Domain geometry ──────────────────────────────────────────
DOMAIN_MIN       = 1
DOMAIN_MAX       = 90
DOMAIN_SPACING_M = 500   # metres per CASCADE domain

# ── LOESS window ─────────────────────────────────────────────
# Physical window width in km — applies to both modes.
# Converted to a frac automatically based on data resolution.
#   2.5 km = 5 domains | 3.5 km = 7 domains | 4.0 km = 8 domains
LOESS_WINDOW_KM = 3.5   # primary smoothing window (7 domains)

# Window sizes (km) tested in sensitivity / comparison figures
COMPARE_WINDOWS_KM = [2.5, 3.5, 5.0]   # 5, 7, 10 domains

# ── Geographic annotations — domain space ────────────────────
TOWN_SPANS = {
    "Buxton":      ( 7,  8),
    "Avon":        (21, 31),
    "Tri-Village": (68, 83),
}
VILLAGE_LINES = {
    "Salvo":    69,
    "Waves":    74,
    "Rodanthe": 80,
}
PIERS  = {"Avon Pier": 26, "Rodanthe Pier": 79}
GROINS = {"Buxton Groin": 6}
WIMBLE_SHOALS = (60, 74)

# ── Geographic annotations — transect / along-coast space ────
# The same features as above, in metres, so the transect-space figures carry
# the reference marks their legend already advertises. Left as None these were
# silently skipped and transect_overview.png shipped a legend for annotations
# it never drew.
#
# Derived from the domain-space values, not measured independently, using the
# two conventions this script already uses:
#   - domain d occupies [(d-1)*500, d*500)      (load_transect_csv, ~line 306)
#   - a domain plots at its band centre (d-0.5)*500   (~line 671)
# So an inclusive span lo..hi runs ((lo-1)*500, hi*500), and a point feature at
# domain d sits at (d-0.5)*500. Re-derive these if DOMAIN_SPACING_M changes.
T_TOWN_SPANS = {
    "Buxton":      ( 3000,  4000),   # domains  7-8
    "Avon":        (10000, 15500),   # domains 21-31
    "Tri-Village": (33500, 41500),   # domains 68-83
}
T_VILLAGE_LINES = {"Salvo": 34250, "Waves": 36750, "Rodanthe": 39750}
T_PIERS         = {"Avon Pier": 12750, "Rodanthe Pier": 39250}
T_GROINS        = {"Buxton Groin": 2750}   # domain 6; see note below
T_WIMBLE_SHOALS = (29500, 37000)           # domains 60-74

# NOTE: GROINS above puts the Buxton groin at domain 6, so its band centre is
# 2750 m. The hindcast places it at GIS 5.5 -- the boundary between domains 5
# and 6, i.e. 2500 m. Translated as-configured rather than silently re-sited;
# the 250 m gap is a question for the domain-space value, not for this one.

# ── Annotation colors ────────────────────────────────────────
C_TOWN_SPAN    = "#90AFC5"
C_WIMBLE       = "#E0A800"
C_VILLAGE_LINE = "0.40"
C_PIER         = "#1565C0"
C_GROIN        = "#B71C1C"

# ── Period colors ────────────────────────────────────────────
C_CS_1984 = "#1F4E79"
C_CS_2004 = "#833C00"
C_WINDOWS  = ["#0072B2", "#E69F00", "#CC79A7"]  # blue=5dom, amber=7dom, pink=10dom

# ── Output ───────────────────────────────────────────────────
OUTPUT_DIR = str(PROJECT_BASE_DIR / "scripts" / "input_prep" / "6-scr-smooth"
                 / "HAT_loess_method_comparison_output")

# ============================================================
# IMPORTS
# ============================================================
import os
import sys

# Windows consoles default to cp1252, which cannot encode the arrows and
# en-dashes in the progress output -- the script died on its first status
# line. UTF-8 here so it runs the same from PyCharm, a terminal or a
# scheduled call.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
from statsmodels.nonparametric.smoothers_lowess import lowess
import warnings
warnings.filterwarnings("ignore")

# Subfolders are created automatically in main()

# ============================================================
# STYLE
# ============================================================
plt.rcParams.update({
    "font.family": "Arial",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linestyle": ":",
})

# ============================================================
# FRAC / WINDOW HELPERS
# ============================================================

def km_to_frac(window_km, n_points, spacing_m):
    """Convert a physical window width (km) to a LOESS frac for n_points."""
    k = (window_km * 1000.0) / spacing_m
    return float(np.clip(k / n_points, 0.02, 1.0))


def estimate_spacing(x_values):
    """Median spacing between consecutive sorted x values (positive diffs only)."""
    arr   = np.sort(np.asarray(x_values, dtype=float))
    diffs = np.diff(arr)
    pos   = diffs[diffs > 0]
    return float(np.median(pos)) if len(pos) else 1.0


def domain_frac(window_km=LOESS_WINDOW_KM):
    """LOESS frac for domain-space smoothing at a given physical window width."""
    n = DOMAIN_MAX - DOMAIN_MIN + 1
    return km_to_frac(window_km, n, DOMAIN_SPACING_M)


def transect_frac(n_transects, spacing_m, window_km=LOESS_WINDOW_KM):
    """LOESS frac for transect-space smoothing at a given physical window width."""
    return km_to_frac(window_km, n_transects, spacing_m)

# ============================================================
# DATA LOADING
# ============================================================

def load_domain_csv(path, period_label):
    """Load domain-averaged LRR summary CSV. Returns standardised DataFrame or None.
    Tolerates common column name variations and strips filename-prefix corruption
    (e.g. column named "domain_lrr_summary.csvdomain_number" instead of "domain_number").
    """
    if path is None or not os.path.exists(path):
        print(f"  Domain CSV ({period_label}): SKIPPED — not found: {path}")
        return None
    df = pd.read_csv(path)

    # Strip any filename prefix accidentally prepended to column names
    # e.g. "domain_lrr_summary.csvdomain_number" -> "domain_number"
    df.columns = [c.split(".csv")[-1] if ".csv" in c else c for c in df.columns]

    # Resolve domain column — try configured name then common alternatives
    domain_col = CS_DOMAIN_COL
    if domain_col not in df.columns:
        # Also check for columns that contain the target name as a substring
        matches = [c for c in df.columns if CS_DOMAIN_COL in c or c in CS_DOMAIN_COL]
        fallbacks = matches + ["domain", "Domain", "domain_id", "DOMAIN"]
        for fb in fallbacks:
            if fb in df.columns:
                domain_col = fb
                print(f"  Note: '{CS_DOMAIN_COL}' not found, using '{domain_col}' instead")
                break
        else:
            print(f"  ERROR: cannot find domain column. Available: {list(df.columns)}")
            return None

    for col in [domain_col, CS_LRR_COL, CS_STD_COL]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=[domain_col, CS_LRR_COL])
    df[domain_col] = df[domain_col].astype(int)
    df = df[(df[domain_col] >= DOMAIN_MIN) & (df[domain_col] <= DOMAIN_MAX)]
    df = df[[domain_col, CS_LRR_COL, CS_STD_COL]].rename(columns={
        domain_col:  "domain",
        CS_LRR_COL:  "cs_lrr",
        CS_STD_COL:  "cs_std",
    }).sort_values("domain").reset_index(drop=True)
    print(f"  Domain CSV ({period_label}): {len(df)} domains  "
          f"LRR range {df['cs_lrr'].min():+.2f}–{df['cs_lrr'].max():+.2f} m/yr")
    return df

def load_transect_csv(path, period_label):
    """
    Load transect-level LRR CSV. Returns standardised DataFrame or None.

    Handles:
      - String transect IDs (e.g. 'usa_NC_0032_0021') — sorted by domain then
        ID string, then replaced with a sequential integer (1, 2, 3 …)
      - Missing along_coast_m — derived by spreading each domain's transects
        evenly across its 500 m band (domain 1 → 0–500 m, domain 2 → 500–1000 m …)
      - Physical spacing for the LOESS frac always estimated from along_coast_m
    """
    if path is None or not os.path.exists(path):
        print(f"  Transect CSV ({period_label}): SKIPPED — not found: {path}")
        return None

    df = pd.read_csv(path)

    # Coerce required numeric columns
    for col in [T_DOMAIN_COL, T_LRR_COL]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if T_STD_COL and T_STD_COL in df.columns:
        df[T_STD_COL] = pd.to_numeric(df[T_STD_COL], errors="coerce")

    # Drop rows missing domain or LRR
    df = df.dropna(subset=[T_TRANSECT_ID_COL, T_LRR_COL, T_DOMAIN_COL])
    df[T_DOMAIN_COL] = df[T_DOMAIN_COL].astype(int)
    df = df[(df[T_DOMAIN_COL] >= DOMAIN_MIN) & (df[T_DOMAIN_COL] <= DOMAIN_MAX)]

    # Optional match-method filter
    if FILTER_POINT_IN_POLYGON and T_MATCH_METHOD_COL in df.columns:
        before = len(df)
        df = df[df[T_MATCH_METHOD_COL] == "point_in_polygon"]
        print(f"  point_in_polygon filter: {before} → {len(df)} transects")

    # Sort by domain then transect ID string (zero-padded suffix sorts correctly)
    df = df.sort_values([T_DOMAIN_COL, T_TRANSECT_ID_COL]).reset_index(drop=True)

    # Replace string transect IDs with sequential integers based on sort order
    df["transect_id"] = np.arange(1, len(df) + 1)

    # Derive along_coast_m from domain position if not present in CSV.
    # Each domain's transects are spread evenly across its 500 m band so that
    # physical spacing can be estimated for the LOESS frac calculation.
    if T_ALONG_COAST_COL is None or T_ALONG_COAST_COL not in df.columns:
        def _spread_within_domain(grp):
            n         = len(grp)
            dom_start = (grp[T_DOMAIN_COL].iloc[0] - 1) * DOMAIN_SPACING_M
            offsets   = (np.arange(n) + 0.5) * (DOMAIN_SPACING_M / n)
            grp       = grp.copy()
            grp["along_coast_m"] = dom_start + offsets
            return grp
        df = df.groupby(T_DOMAIN_COL, group_keys=False).apply(_spread_within_domain)
        print(f"  along_coast_m: derived from domain position (not in CSV)")
    else:
        df["along_coast_m"] = pd.to_numeric(df[T_ALONG_COAST_COL], errors="coerce")

    # Standardise remaining column names
    rename = {T_DOMAIN_COL: "domain", T_LRR_COL: "lrr"}
    if T_STD_COL and T_STD_COL in df.columns:
        rename[T_STD_COL] = "lrr_std"
    df = df.rename(columns=rename).reset_index(drop=True)

    spacing = estimate_spacing(df["along_coast_m"].values)
    print(f"  Transect CSV ({period_label}): {len(df)} transects  "
          f"est. spacing {spacing:.1f} m  "
          f"LRR range {df['lrr'].min():+.2f}–{df['lrr'].max():+.2f} m/yr")
    return df

# ============================================================
# SMOOTHING
# ============================================================

def apply_loess(x, values, frac):
    """LOESS smoother. Returns smoothed array at the same x positions."""
    valid = ~np.isnan(values)
    if valid.sum() < 5:
        return values.copy()
    smoothed = np.full_like(values, np.nan, dtype=float)
    result   = lowess(values[valid], x[valid], frac=frac, return_sorted=True)
    smoothed[valid] = np.interp(x[valid], result[:, 0], result[:, 1])
    return smoothed


def smooth_domain_df(df, window_km=LOESS_WINDOW_KM):
    """Apply LOESS in domain space. Returns copy of df with cs_lrr_smooth column."""
    df   = df.copy()
    frac = domain_frac(window_km)
    df["cs_lrr_smooth"] = apply_loess(
        df["domain"].values.astype(float),
        df["cs_lrr"].values,
        frac,
    )
    return df


def smooth_transect_df(df, window_km=LOESS_WINDOW_KM):
    """
    Apply LOESS in transect space. Returns copy of df with lrr_smooth column.

    Physical spacing for the frac is always estimated from along_coast_m
    (derived or real), so the window is correct in physical kilometres
    regardless of whether transect_id or along_coast_m is the plot x-axis.
    """
    df      = df.copy()
    spacing = estimate_spacing(df["along_coast_m"].values)
    frac    = transect_frac(len(df), spacing, window_km)

    x = (df["along_coast_m"].values.astype(float)
         if TRANSECT_X_AXIS == "along_coast_m"
         else df["transect_id"].values.astype(float))

    df["lrr_smooth"] = apply_loess(x, df["lrr"].values, frac)
    df["_x_smooth"]  = x   # stored so plot functions don't recompute
    return df


def aggregate_to_domains(t_df):
    """
    Average smoothed (and raw) transect values within each CASCADE domain.

    Returns a domain-level DataFrame matching the domain CSV schema so all
    domain-space plot functions work unchanged:
      domain        — CASCADE domain number
      cs_lrr        — mean of raw transect LRRs within the domain
      cs_std        — std  of raw transect LRRs within the domain
      cs_lrr_smooth — mean of smoothed transect LRRs within the domain
    """
    grp = t_df.groupby("domain")
    domain_df = pd.DataFrame({
        "domain":        grp["lrr"].mean().index,
        "cs_lrr":        grp["lrr"].mean().values,
        "cs_std":        grp["lrr"].std(ddof=1).values,
        "cs_lrr_smooth": grp["lrr_smooth"].mean().values,
    }).reset_index(drop=True)
    return domain_df.sort_values("domain").reset_index(drop=True)

# ============================================================
# ANNOTATION HELPERS
# ============================================================

def _add_annotations(ax, town_spans, wimble, village_lines, piers, groins):
    """Draw geographic reference annotations. Any position set to None is skipped."""
    trans = blended_transform_factory(ax.transData, ax.transAxes)

    lo, hi = wimble
    if lo is not None and hi is not None:
        ax.axvspan(lo, hi, color=C_WIMBLE, alpha=0.12, zorder=0)
        ax.text((lo + hi) / 2, 0.03, "Wimble\nShoals",
                transform=trans, fontsize=7, ha="center",
                color=C_WIMBLE, style="italic")

    for name, (lo, hi) in town_spans.items():
        if lo is None or hi is None:
            continue
        ax.axvspan(lo, hi, color=C_TOWN_SPAN, alpha=0.18, zorder=0)
        ax.text((lo + hi) / 2, 0.94, name,
                transform=trans, fontsize=7.5, ha="center",
                color="#2C5F80", fontweight="bold")

    for name, pos in village_lines.items():
        if pos is None:
            continue
        ax.axvline(pos, color=C_VILLAGE_LINE, lw=0.8, ls="--", zorder=2)
        ax.text(pos, 0.88, name, transform=trans, fontsize=6.5,
                ha="center", color=C_VILLAGE_LINE, rotation=90, va="top")

    for name, pos in piers.items():
        if pos is None:
            continue
        ax.axvline(pos, color=C_PIER, lw=1.0, ls="-.", zorder=2)
        ax.text(pos, 0.76, name, transform=trans, fontsize=6.5,
                ha="center", color=C_PIER, rotation=90, va="top")

    for name, pos in groins.items():
        if pos is None:
            continue
        ax.axvline(pos, color=C_GROIN, lw=1.0, ls=":", zorder=2)
        ax.text(pos, 0.76, name, transform=trans, fontsize=6.5,
                ha="center", color=C_GROIN, rotation=90, va="top")


def add_domain_annotations(ax):
    _add_annotations(ax, TOWN_SPANS, WIMBLE_SHOALS, VILLAGE_LINES, PIERS, GROINS)


def add_transect_annotations(ax):
    _add_annotations(ax, T_TOWN_SPANS, T_WIMBLE_SHOALS,
                     T_VILLAGE_LINES, T_PIERS, T_GROINS)


def annotation_legend_handles():
    return [
        Patch(color=C_TOWN_SPAN,  alpha=0.4,  label="Community span"),
        Patch(color=C_WIMBLE,     alpha=0.25, label="Wimble Shoals influence"),
        Line2D([0], [0], color=C_VILLAGE_LINE, lw=1, ls="--", label="Village center"),
        Line2D([0], [0], color=C_PIER,         lw=1, ls="-.", label="Pier"),
        Line2D([0], [0], color=C_GROIN,        lw=1, ls=":",  label="Groin"),
    ]


def style_domain_axis(ax, is_bottom=True):
    ax.set_xlim(DOMAIN_MIN - 0.5, DOMAIN_MAX + 0.5)
    ax.axhline(0, color="k", lw=0.6, ls="--", alpha=0.4)
    if is_bottom:
        ax.set_xlabel("CASCADE Domain Number", fontsize=11, fontweight="bold")


def style_transect_axis(ax, x_values, is_bottom=True):
    ax.set_xlim(x_values.min() - 1, x_values.max() + 1)
    ax.axhline(0, color="k", lw=0.6, ls="--", alpha=0.4)
    if is_bottom:
        label = ("Along-coast Distance (m)" if TRANSECT_X_AXIS == "along_coast_m"
                 else "Transect ID (sequential)")
        ax.set_xlabel(label, fontsize=11, fontweight="bold")

# ============================================================
# DOMAIN-SPACE FIGURES
# Works identically for both modes — receives a domain-level DataFrame
# regardless of whether it came from load_domain_csv or aggregate_to_domains.
# ============================================================

def _domain_two_panel(d1984, d2004, show_raw, out_path, suptitle):
    configs   = [(d1984, "1984–2004", C_CS_1984),
                 (d2004, "2004–2024", C_CS_2004)]
    fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True)
    fig.suptitle(suptitle, fontsize=14, fontweight="bold", y=1.01)

    for i, (ax, (df, period, color)) in enumerate(zip(axes, configs)):
        is_bottom = (i == len(configs) - 1)
        if df is None:
            ax.text(0.5, 0.5, f"No data — {period}",
                    transform=ax.transAxes, ha="center", fontsize=12)
            continue
        if show_raw:
            ax.plot(df["domain"], df["cs_lrr"],
                    color=color, lw=0.8, alpha=0.30, marker="o", ms=2,
                    zorder=1, label="Raw domain LRR")
            ax.fill_between(df["domain"],
                            df["cs_lrr"] - df["cs_std"],
                            df["cs_lrr"] + df["cs_std"],
                            color=color, alpha=0.07, zorder=0)
        ax.plot(df["domain"], df["cs_lrr_smooth"],
                color=color, lw=2.5, zorder=3, label="LOESS smoothed")
        ax.set_ylabel("Rate (m/yr)", fontsize=10, fontweight="bold")
        ax.set_title(period, fontsize=11, fontweight="bold", loc="left", pad=5)
        style_domain_axis(ax, is_bottom)
        add_domain_annotations(ax)
        handles, _ = ax.get_legend_handles_labels()
        if is_bottom:
            handles += annotation_legend_handles()
        ax.legend(handles=handles, fontsize=8.5, framealpha=0.95,
                  loc="upper center", ncol=2 if is_bottom else 1)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


def plot_domain_overview(d1984, d2004, out_path, suffix=""):
    _domain_two_panel(
        d1984, d2004, show_raw=True, out_path=out_path,
        suptitle=(f"Shoreline Change Rate — Transect-Based Smoothing\nHatteras Island, NC{suffix}"),
    )


def plot_domain_smoothed_only(d1984, d2004, out_path, suffix=""):
    _domain_two_panel(
        d1984, d2004, show_raw=False, out_path=out_path,
        suptitle=(f"Shoreline Change Rate — Transect-Based, Smoothed Curve\nHatteras Island, NC{suffix}"),
    )


def plot_domain_combined(d1984, d2004, out_path, suffix=""):
    fig, ax = plt.subplots(figsize=(16, 6))
    fig.suptitle(
        f"Shoreline Change Rate — Both Periods, Transect-Based\nHatteras Island, NC{suffix}",
        fontsize=14, fontweight="bold",
    )
    for df, period, color in [(d1984, "1984–2004", C_CS_1984),
                               (d2004, "2004–2024", C_CS_2004)]:
        if df is None:
            continue
        ax.fill_between(df["domain"],
                        df["cs_lrr_smooth"] - df["cs_std"],
                        df["cs_lrr_smooth"] + df["cs_std"],
                        color=color, alpha=0.10)
        ax.plot(df["domain"], df["cs_lrr_smooth"],
                color=color, lw=2.5, label=period)
    ax.set_ylabel("Rate (m/yr)", fontsize=11, fontweight="bold")
    style_domain_axis(ax, is_bottom=True)
    add_domain_annotations(ax)
    handles, _ = ax.get_legend_handles_labels()
    handles += annotation_legend_handles()
    ax.legend(handles=handles, fontsize=9, framealpha=0.95, loc="upper center", ncol=2)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


def plot_domain_sensitivity(df, period_label, color, out_path, suffix=""):
    """3-panel bandwidth sensitivity figure for a single period (domain space)."""
    fracs  = [domain_frac(w) for w in COMPARE_WINDOWS_KM]
    labels = [f"{w:.1f} km  ({int(round(w * 1000 / DOMAIN_SPACING_M))} domains,"
              f"  frac={f:.3f})"
              for w, f in zip(COMPARE_WINDOWS_KM, fracs)]

    fig, axes = plt.subplots(len(COMPARE_WINDOWS_KM), 1,
                             figsize=(16, 4 * len(COMPARE_WINDOWS_KM)), sharex=True)
    fig.suptitle(
        f"Window Sensitivity — Domain-Averaged Smoothing: {period_label}\nHatteras Island, NC{suffix}",
        fontsize=13, fontweight="bold", y=1.01,
    )
    for j, (ax, frac, lbl, km) in enumerate(zip(axes, fracs, labels, COMPARE_WINDOWS_KM)):
        is_bottom = (j == len(fracs) - 1)
        m = smooth_domain_df(df, window_km=km)
        ax.plot(df["domain"], df["cs_lrr"],
                color=color, lw=0.8, alpha=0.25, marker="o", ms=2, label="Raw")
        ax.plot(m["domain"], m["cs_lrr_smooth"],
                color=color, lw=2.5, label="LOESS")
        ax.text(0.01, 0.97, lbl, transform=ax.transAxes, fontsize=8.5, va="top",
                bbox=dict(boxstyle="round", fc="white", alpha=0.88, ec="0.7"))
        ax.set_ylabel("Rate (m/yr)", fontsize=10, fontweight="bold")
        style_domain_axis(ax, is_bottom)
        add_domain_annotations(ax)
        handles, _ = ax.get_legend_handles_labels()
        if is_bottom:
            handles += annotation_legend_handles()
        ax.legend(handles=handles, fontsize=8.5, framealpha=0.95,
                  loc="upper center", ncol=2 if is_bottom else 1)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


def plot_domain_window_comparison(d1984, d2004, out_path, suffix=""):
    """All window sizes overlaid in domain space — both periods."""
    configs   = [(d1984, "1984–2004", C_CS_1984),
                 (d2004, "2004–2024", C_CS_2004)]
    fig, axes = plt.subplots(2, 1, figsize=(16, 11), sharex=True)
    fig.suptitle(
        f"Window Comparison — Domain-Averaged Smoothing\nHatteras Island, NC{suffix}",
        fontsize=14, fontweight="bold", y=1.01,
    )
    for i, (ax, (df, period, pcol)) in enumerate(zip(axes, configs)):
        is_bottom = (i == len(configs) - 1)
        if df is None:
            ax.text(0.5, 0.5, f"No data — {period}",
                    transform=ax.transAxes, ha="center", fontsize=12)
            continue
        ax.plot(df["domain"], df["cs_lrr"],
                color=pcol, lw=1.0, alpha=0.25, marker="o", ms=2.5, label="Raw")
        for km, wc in zip(COMPARE_WINDOWS_KM, C_WINDOWS):
            m    = smooth_domain_df(df, window_km=km)
            ndom = int(round(km * 1000 / DOMAIN_SPACING_M))
            ax.plot(m["domain"], m["cs_lrr_smooth"],
                    color=wc, lw=2.5, label=f"{km:.1f} km  ({ndom} domains)")
        ax.set_ylabel("Rate (m/yr)", fontsize=11, fontweight="bold")
        ax.set_title(period, fontsize=12, fontweight="bold", loc="left", pad=5)
        style_domain_axis(ax, is_bottom)
        add_domain_annotations(ax)
        handles, _ = ax.get_legend_handles_labels()
        if is_bottom:
            handles += annotation_legend_handles()
        ax.legend(handles=handles, fontsize=8.5, framealpha=0.95,
                  loc="upper center", ncol=2 if is_bottom else 1)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")

# ============================================================
# TRANSECT-SPACE FIGURES  (produced only in "transect" mode)
# ============================================================

def plot_transect_overview(t1984, t2004, d1984, d2004, out_path):
    """
    Raw transect scatter + LOESS smoothed curve in along-coast space,
    with domain-averaged LRR overlaid as a dashed line with open markers.
    Requested by advisor to show how much variability domain averaging collapses.
    """
    configs   = [(t1984, d1984, "1984–2004", C_CS_1984),
                 (t2004, d2004, "2004–2024", C_CS_2004)]
    fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True)
    fig.suptitle(
        "Shoreline Change Rate — Individual Transects and Domain Averages\nHatteras Island, NC",
        fontsize=14, fontweight="bold", y=1.01,
    )
    for i, (ax, (t_df, d_df, period, color)) in enumerate(zip(axes, configs)):
        is_bottom = (i == len(configs) - 1)
        if t_df is None:
            ax.text(0.5, 0.5, f"No data — {period}",
                    transform=ax.transAxes, ha="center", fontsize=12)
            continue
        x = t_df["_x_smooth"].values

        # Domain boundary lines at every 500 m — drawn first so they sit behind data
        for domain_n in range(DOMAIN_MIN, DOMAIN_MAX + 2):
            boundary_m = (domain_n - 1) * DOMAIN_SPACING_M
            ax.axvline(boundary_m, color="0.82", lw=0.4, ls="-", zorder=0, alpha=0.8)

        # Raw transect LRRs as scatter
        ax.scatter(x, t_df["lrr"], color=color, s=4, alpha=0.20,
                   zorder=1, label="Individual transect LRR")
        # LOESS smoothed transect curve
        ax.plot(x, t_df["lrr_smooth"], color=color, lw=2.5,
                zorder=3, label="LOESS smoothed (transect-based)")
        # Domain averages: x = centre of each 500 m domain in along-coast metres
        if d_df is not None:
            x_dom = (d_df["domain"].values - 0.5) * DOMAIN_SPACING_M
            ax.plot(x_dom, d_df["cs_lrr"].values,
                    color="black", lw=1.5, ls="--",
                    marker="o", ms=5, markerfacecolor="white", markeredgewidth=1.5,
                    zorder=4, alpha=0.75, label="Domain average LRR")
        ax.set_ylabel("Rate (m/yr)", fontsize=10, fontweight="bold")
        ax.set_title(period, fontsize=11, fontweight="bold", loc="left", pad=5)
        style_transect_axis(ax, x, is_bottom)
        add_transect_annotations(ax)
        handles, _ = ax.get_legend_handles_labels()
        if is_bottom:
            handles += annotation_legend_handles()
        ncols = 4 if is_bottom else 1
        ax.legend(handles=handles, fontsize=8.5, framealpha=0.95,
                  loc="lower center", ncol=ncols)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")

def plot_transect_window_comparison(t1984, t2004, out_path):
    """Window sensitivity in transect space — all km windows overlaid, both periods."""
    configs = [(t1984, "1984–2004", C_CS_1984),
               (t2004, "2004–2024", C_CS_2004)]
    xlabel = ("Along-coast Distance (m)" if TRANSECT_X_AXIS == "along_coast_m"
             else "Transect ID (sequential)")
    fig, axes = plt.subplots(2, 1, figsize=(16, 11), sharex=True)
    fig.suptitle(
        f"Window Comparison — Transect Space\nHatteras Island, NC",
        fontsize=14, fontweight="bold", y=1.01,
    )
    for i, (ax, (df, period, pcol)) in enumerate(zip(axes, configs)):
        is_bottom = (i == len(configs) - 1)
        if df is None:
            ax.text(0.5, 0.5, f"No data — {period}",
                    transform=ax.transAxes, ha="center", fontsize=12)
            continue
        x       = df["_x_smooth"].values
        spacing = estimate_spacing(df["along_coast_m"].values)
        ax.scatter(x, df["lrr"], color=pcol, s=3, alpha=0.18, label="Raw")
        for km, wc in zip(COMPARE_WINDOWS_KM, C_WINDOWS):
            frac     = transect_frac(len(df), spacing, km)
            smoothed = apply_loess(x, df["lrr"].values, frac)
            ndom     = int(round(km * 1000 / DOMAIN_SPACING_M))
            ax.plot(x, smoothed, color=wc, lw=2.5,
                    label=f"{km:.1f} km  ({ndom} domains,  frac={frac:.3f})")
        ax.set_ylabel("Rate (m/yr)", fontsize=11, fontweight="bold")
        ax.set_title(period, fontsize=12, fontweight="bold", loc="left", pad=5)
        style_transect_axis(ax, x, is_bottom)
        add_transect_annotations(ax)
        handles, _ = ax.get_legend_handles_labels()
        if is_bottom:
            handles += annotation_legend_handles()
        ax.legend(handles=handles, fontsize=8.5, framealpha=0.95,
                  loc="upper center", ncol=2 if is_bottom else 1)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


# ============================================================
# TRANSECT-SMOOTHED WINDOWS IN DOMAIN SPACE
# Smooths at transect level for each window, aggregates to domain
# resolution, then plots with domain number and geographic annotations.
# ============================================================

def plot_transect_windows_domain_space(t1984, t2004, out_path):
    """
    For each window in COMPARE_WINDOWS_KM:
      1. Apply LOESS at transect resolution
      2. Aggregate smoothed values to domain means
      3. Plot against CASCADE domain number with geographic annotations

    Replaces plot_domain_window_comparison in transect mode so all curves
    shown are transect-based — no domain-averaged smoothing is mixed in.
    """
    configs   = [(t1984, "1984\u20132004", C_CS_1984),
                 (t2004, "2004\u20132024", C_CS_2004)]
    fig, axes = plt.subplots(2, 1, figsize=(16, 11), sharex=True)
    fig.suptitle(
        "Window Comparison — Transect-Based Smoothing\nHatteras Island, NC",
        fontsize=14, fontweight="bold", y=1.01,
    )

    for i, (ax, (t_df, period, pcol)) in enumerate(zip(axes, configs)):
        is_bottom = (i == len(configs) - 1)
        if t_df is None:
            ax.text(0.5, 0.5, f"No data — {period}",
                    transform=ax.transAxes, ha="center", fontsize=12)
            continue

        # Raw domain means — mean of raw transect LRRs, same for all windows
        d_raw = aggregate_to_domains(t_df)
        ax.plot(d_raw["domain"], d_raw["cs_lrr"],
                color=pcol, lw=1.0, alpha=0.25, marker="o", ms=2.5, label="Raw")

        spacing = estimate_spacing(t_df["along_coast_m"].values)

        # One smoothed curve per window: smooth transects → aggregate to domains
        for km, wc in zip(COMPARE_WINDOWS_KM, C_WINDOWS):
            frac     = transect_frac(len(t_df), spacing, km)
            t_smooth = smooth_transect_df(t_df, window_km=km)
            d_agg    = aggregate_to_domains(t_smooth)
            ndom     = int(round(km * 1000 / DOMAIN_SPACING_M))
            ax.plot(d_agg["domain"], d_agg["cs_lrr_smooth"],
                    color=wc, lw=2.5,
                    label=f"{km:.1f} km  ({ndom} domains,  frac={frac:.3f})")

        ax.set_ylabel("Rate (m/yr)", fontsize=11, fontweight="bold")
        ax.set_title(period, fontsize=12, fontweight="bold", loc="left", pad=5)
        style_domain_axis(ax, is_bottom)
        add_domain_annotations(ax)
        handles, _ = ax.get_legend_handles_labels()
        if is_bottom:
            handles += annotation_legend_handles()
        ax.legend(handles=handles, fontsize=8.5, framealpha=0.95,
                  loc="upper center", ncol=2 if is_bottom else 1)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")

# ============================================================

# ============================================================
# TRANSECT SENSITIVITY (matches plot_transect_windows_domain_space)
# ============================================================

def plot_transect_sensitivity(t_df, period_label, color, out_path, suffix=""):
    """
    3-panel bandwidth sensitivity — transect mode.
    Smooths at transect level for each window then aggregates to domains,
    matching plot_transect_windows_domain_space exactly.
    """
    fig, axes = plt.subplots(len(COMPARE_WINDOWS_KM), 1,
                             figsize=(16, 4 * len(COMPARE_WINDOWS_KM)), sharex=True)
    fig.suptitle(
        f"Window Sensitivity — Transect-Based Smoothing: {period_label}\nHatteras Island, NC{suffix}",
        fontsize=13, fontweight="bold", y=1.01,
    )
    spacing = estimate_spacing(t_df["along_coast_m"].values)
    d_raw   = aggregate_to_domains(t_df)

    for j, (ax, km) in enumerate(zip(axes, COMPARE_WINDOWS_KM)):
        is_bottom = (j == len(COMPARE_WINDOWS_KM) - 1)
        frac  = transect_frac(len(t_df), spacing, km)
        ndom  = int(round(km * 1000 / DOMAIN_SPACING_M))
        label = f"{km:.1f} km  ({ndom} domains,  frac={frac:.3f})"
        t_smooth = smooth_transect_df(t_df, window_km=km)
        d_agg    = aggregate_to_domains(t_smooth)
        ax.plot(d_raw["domain"], d_raw["cs_lrr"],
                color=color, lw=0.8, alpha=0.25, marker="o", ms=2, label="Raw")
        ax.plot(d_agg["domain"], d_agg["cs_lrr_smooth"],
                color=color, lw=2.5, label="LOESS")
        ax.text(0.01, 0.97, label, transform=ax.transAxes, fontsize=8.5, va="top",
                bbox=dict(boxstyle="round", fc="white", alpha=0.88, ec="0.7"))
        ax.set_ylabel("Rate (m/yr)", fontsize=10, fontweight="bold")
        style_domain_axis(ax, is_bottom)
        add_domain_annotations(ax)
        handles, _ = ax.get_legend_handles_labels()
        if is_bottom:
            handles += annotation_legend_handles()
        ax.legend(handles=handles, fontsize=8.5, framealpha=0.95,
                  loc="center", ncol=2 if is_bottom else 1)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")



# ============================================================
# METHOD COMPARISON — transect-based vs domain-averaged LOESS
# Both smoothing approaches overlaid for each window size.
# This is the figure Laura requested to directly compare methods.
# ============================================================

def plot_method_comparison(t1984, t2004, da1984, da2004, out_path):
    """
    For each window in COMPARE_WINDOWS_KM, plots both:
      — Solid line  : transect-based LOESS (smooth transects → aggregate to domains)
      — Dashed line : domain-averaged LOESS (smooth domain means directly)

    Both are shown in domain space with geographic annotations.
    Window sizes are labelled in km and equivalent domain count.

    This directly answers Laura's question: are the two approaches different?
    """
    configs   = [(t1984, da1984, "1984–2004", C_CS_1984),
                 (t2004, da2004, "2004–2024", C_CS_2004)]
    fig, axes = plt.subplots(2, 1, figsize=(16, 11), sharex=True)
    fig.suptitle(
        "LOESS Smoothing — Transect-Based vs Domain-Averaged\n"
        "Hatteras Island, NC  (solid = transect-based,  dashed = domain-averaged)",
        fontsize=14, fontweight="bold", y=1.01,
    )

    for i, (ax, (t_df, da_df, period, pcol)) in enumerate(zip(axes, configs)):
        is_bottom = (i == len(configs) - 1)
        if t_df is None and da_df is None:
            ax.text(0.5, 0.5, f"No data — {period}",
                    transform=ax.transAxes, ha="center", fontsize=12)
            continue

        # Raw domain means as faint background markers
        if t_df is not None:
            d_raw = aggregate_to_domains(t_df)
            ax.plot(d_raw["domain"], d_raw["cs_lrr"],
                    color=pcol, lw=0, marker="o", ms=3, alpha=0.20,
                    zorder=1, label="Raw domain LRR")

        spacing = estimate_spacing(t_df["along_coast_m"].values) if t_df is not None else None

        handles_extra = []
        for km, wc in zip(COMPARE_WINDOWS_KM, C_WINDOWS):
            ndom  = int(round(km * 1000 / DOMAIN_SPACING_M))
            label = f"{km:.1f} km  ({ndom} domains)"

            # — Transect-based: smooth transects then aggregate
            if t_df is not None:
                frac_t   = transect_frac(len(t_df), spacing, km)
                t_smooth = smooth_transect_df(t_df, window_km=km)
                d_agg    = aggregate_to_domains(t_smooth)
                line_t, = ax.plot(d_agg["domain"], d_agg["cs_lrr_smooth"],
                                  color=wc, lw=2.5, ls="-", zorder=4,
                                  label=f"{label} — transect")

            # — Domain-averaged: smooth the 90 domain means directly
            if da_df is not None:
                m = smooth_domain_df(da_df, window_km=km)
                line_d, = ax.plot(m["domain"], m["cs_lrr_smooth"],
                                  color=wc, lw=2.5, ls="--", zorder=3, alpha=0.80,
                                  label=f"{label} — domain avg")

        ax.set_ylabel("Rate (m/yr)", fontsize=11, fontweight="bold")
        ax.set_title(period, fontsize=12, fontweight="bold", loc="left", pad=5)
        style_domain_axis(ax, is_bottom)
        add_domain_annotations(ax)

        # Build legend: raw + window pairs + annotations (bottom panel only)
        handles, labels = ax.get_legend_handles_labels()
        if is_bottom:
            handles += annotation_legend_handles()
        ax.legend(handles=handles, fontsize=8.0, framealpha=0.95,
                  loc="lower center", ncol=3 if is_bottom else 2)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


def plot_method_comparison_single(t1984, t2004, da1984, da2004,
                                   window_km, out_path):
    """
    Single-window method comparison: transect-based vs domain-averaged LOESS.
    Shows one window size only so the two curves can be read clearly.
    """
    ndom  = int(round(window_km * 1000 / DOMAIN_SPACING_M))
    label = f"{window_km:.1f} km  ({ndom} domains)"

    configs   = [(t1984, da1984, "1984–2004", C_CS_1984),
                 (t2004, da2004, "2004–2024", C_CS_2004)]
    fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True)
    fig.suptitle(
        f"LOESS Smoothing Comparison — {label}\n"
        "Hatteras Island, NC  "
        "(solid = transect-based,  dashed = domain-averaged)",
        fontsize=14, fontweight="bold", y=1.01,
    )

    for i, (ax, (t_df, da_df, period, pcol)) in enumerate(zip(axes, configs)):
        is_bottom = (i == len(configs) - 1)
        if t_df is None and da_df is None:
            ax.text(0.5, 0.5, f"No data — {period}",
                    transform=ax.transAxes, ha="center", fontsize=12)
            continue

        # Raw domain means as faint background dots
        if t_df is not None:
            d_raw = aggregate_to_domains(t_df)
            ax.plot(d_raw["domain"], d_raw["cs_lrr"],
                    color=pcol, lw=0, marker="o", ms=3, alpha=0.20,
                    zorder=1, label="Raw domain LRR")

        # Transect-based LOESS (solid)
        if t_df is not None:
            spacing = estimate_spacing(t_df["along_coast_m"].values)
            t_smooth = smooth_transect_df(t_df, window_km=window_km)
            d_agg    = aggregate_to_domains(t_smooth)
            ax.plot(d_agg["domain"], d_agg["cs_lrr_smooth"],
                    color=pcol, lw=2.8, ls="-", zorder=4,
                    label=f"Transect-based LOESS  ({label})")

        # Domain-averaged LOESS (dashed)
        if da_df is not None:
            m = smooth_domain_df(da_df, window_km=window_km)
            ax.plot(m["domain"], m["cs_lrr_smooth"],
                    color=pcol, lw=2.8, ls="--", zorder=3, alpha=0.75,
                    label=f"Domain-averaged LOESS  ({label})")

        ax.set_ylabel("Rate (m/yr)", fontsize=11, fontweight="bold")
        ax.set_title(period, fontsize=12, fontweight="bold", loc="left", pad=5)
        style_domain_axis(ax, is_bottom)
        add_domain_annotations(ax)
        handles, _ = ax.get_legend_handles_labels()
        if is_bottom:
            handles += annotation_legend_handles()
        ax.legend(handles=handles, fontsize=8.5, framealpha=0.95,
                  loc="lower center", ncol=3 if is_bottom else 2)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")

# ============================================================
# MAIN
# Always runs both transect-based and domain-averaged smoothing.
# Outputs are organised into clearly labelled subfolders.
# ============================================================

def main():
    print("=" * 65)
    print("CoastSat LRR Smoothing — Hatteras Island")
    print(f"  X-axis      : {TRANSECT_X_AXIS}")
    print(f"  Window      : {LOESS_WINDOW_KM} km  "
          f"({int(round(LOESS_WINDOW_KM * 1000 / DOMAIN_SPACING_M))} domains)")
    print(f"  Compare     : {COMPARE_WINDOWS_KM} km  "
          f"({[int(round(w * 1000 / DOMAIN_SPACING_M)) for w in COMPARE_WINDOWS_KM]} domains)")
    print("  Both methods always run — outputs saved to subfolders.")
    print("=" * 65)

    # ── Create subfolders ──────────────────────────────────────────
    DIR_T       = os.path.join(OUTPUT_DIR, "01_transect_based")
    DIR_D       = os.path.join(OUTPUT_DIR, "02_domain_averaged")
    DIR_C       = os.path.join(OUTPUT_DIR, "03_cascade_inputs")
    DIR_COMPARE = os.path.join(OUTPUT_DIR, "04_method_comparison")
    for d in [DIR_T, DIR_D, DIR_C, DIR_COMPARE]:
        os.makedirs(d, exist_ok=True)

    # ── Load transect data ──────────────────────────────────────
    print("\nLoading transect data...")
    t1984_raw = load_transect_csv(TRANSECT_CSV_1984_2004, "1984–2004")
    t2004_raw = load_transect_csv(TRANSECT_CSV_2004_2024, "2004–2024")
    t1984 = smooth_transect_df(t1984_raw) if t1984_raw is not None else None
    t2004 = smooth_transect_df(t2004_raw) if t2004_raw is not None else None
    td1984 = aggregate_to_domains(t1984) if t1984 is not None else None
    td2004 = aggregate_to_domains(t2004) if t2004 is not None else None

    # ── Load domain-averaged data ────────────────────────────────
    print("\nLoading domain-averaged data...")
    da1984 = load_domain_csv(DOMAIN_CSV_1984_2004, "1984–2004")
    da2004 = load_domain_csv(DOMAIN_CSV_2004_2024, "2004–2024")
    if da1984 is not None: da1984 = smooth_domain_df(da1984)
    if da2004 is not None: da2004 = smooth_domain_df(da2004)

    # ── 01: Transect-based figures ───────────────────────────────
    print("\n[01] Transect-based figures → 01_transect_based/")

    # Transect overview: raw scatter + LOESS + domain averages overlaid
    plot_transect_overview(t1984, t2004, td1984, td2004,
        os.path.join(DIR_T, "transect_overview.png"))

    # Domain-space overview using transect-smoothed values
    plot_domain_overview(td1984, td2004,
        os.path.join(DIR_T, "overview_smoothed.png"), "  [transect-based]")
    plot_domain_smoothed_only(td1984, td2004,
        os.path.join(DIR_T, "smoothed_only.png"), "  [transect-based]")
    plot_domain_combined(td1984, td2004,
        os.path.join(DIR_T, "combined_periods.png"), "  [transect-based]")

    # Window comparison in domain space (transect-smoothed)
    plot_transect_windows_domain_space(t1984, t2004,
        os.path.join(DIR_T, "window_comparison.png"))

    # Sensitivity — transect-based
    if t1984 is not None:
        plot_transect_sensitivity(t1984, "1984–2004", C_CS_1984,
            os.path.join(DIR_T, "sensitivity_1984_2004.png"))
    if t2004 is not None:
        plot_transect_sensitivity(t2004, "2004–2024", C_CS_2004,
            os.path.join(DIR_T, "sensitivity_2004_2024.png"))

    # ── 02: Domain-averaged figures ──────────────────────────────
    print("\n[02] Domain-averaged figures → 02_domain_averaged/")

    plot_domain_window_comparison(da1984, da2004,
        os.path.join(DIR_D, "window_comparison.png"), "  [domain-averaged]")

    if da1984 is not None:
        plot_domain_sensitivity(da1984, "1984–2004", C_CS_1984,
            os.path.join(DIR_D, "sensitivity_1984_2004.png"))
    if da2004 is not None:
        plot_domain_sensitivity(da2004, "2004–2024", C_CS_2004,
            os.path.join(DIR_D, "sensitivity_2004_2024.png"))

    # ── 04: Method comparison (Laura's request) ────────────────────
    print("\n[04] Method comparison figures → 04_method_comparison/")

    plot_method_comparison(t1984, t2004, da1984, da2004,
        os.path.join(DIR_COMPARE, "transect_vs_domain_smoothing.png"))

    # Individual window figures — one per window size
    for km in COMPARE_WINDOWS_KM:
        ndom  = int(round(km * 1000 / DOMAIN_SPACING_M))
        fname = f"transect_vs_domain_{ndom}domains_{km:.1f}km.png"
        plot_method_comparison_single(t1984, t2004, da1984, da2004,
            km, os.path.join(DIR_COMPARE, fname))

    # ── 03: CASCADE inputs (transect-based) ───────────────────────
    print("\n[03] Exporting CASCADE inputs → 03_cascade_inputs/")

    parts = []
    for df, lbl in [(td1984, "1984_2004"), (td2004, "2004_2024")]:
        if df is None:
            continue
        t = df[["domain", "cs_lrr", "cs_std", "cs_lrr_smooth"]].copy()
        t.columns = ["domain"] + [f"{c}_{lbl}" for c in ["cs_lrr", "cs_std", "cs_lrr_smooth"]]
        parts.append(t)
    if parts:
        from functools import reduce
        table = reduce(lambda a, b: a.merge(b, on="domain", how="outer"), parts)
        fname = "cascade_lrr_inputs_transect_based.csv"
        table.sort_values("domain").to_csv(os.path.join(DIR_C, fname), index=False)
        print(f"  Saved: {fname}")

    # Also save full transect-level tables for reference
    for df, lbl in [(t1984, "1984_2004"), (t2004, "2004_2024")]:
        if df is None:
            continue
        cols = ["transect_id", "along_coast_m", "domain", "lrr", "lrr_smooth"]
        if "lrr_std" in df.columns:
            cols.insert(4, "lrr_std")
        fname = f"transect_lrr_{lbl}.csv"
        df[cols].to_csv(os.path.join(DIR_C, fname), index=False)
        print(f"  Saved: {fname}")

    print("\n" + "=" * 65)
    print("Done!  Output structure:")
    print(f"  01_transect_based/   — transect-smoothed figures")
    print(f"  02_domain_averaged/  — domain-averaged figures")
    print(f"  03_cascade_inputs/   — CASCADE-ready CSV + transect tables")
    print("=" * 65)


if __name__ == "__main__":
    main()
