"""
CoastSat LRR Smoothing — Hatteras Island
=========================================
Supports two smoothing modes controlled by SMOOTHING_MODE in CONFIG:

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

# ============================================================
# CONFIG
# ============================================================

# ── Smoothing mode ───────────────────────────────────────────
SMOOTHING_MODE  = "domain"      # "domain" | "transect"
TRANSECT_X_AXIS = "along_coast_m"   # "transect_id" | "along_coast_m"
#   transect_id   : LOESS x = sequential integer (1, 2, 3 …)
#   along_coast_m : LOESS x = derived along-coast distance in metres

# ── Domain-mode inputs ───────────────────────────────────────
DOMAIN_CSV_1984_2004 = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\1984_2004\domain_lrr_summary.csv"
DOMAIN_CSV_2004_2024 = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\2004_2024\domain_lrr_summary.csv"
CS_DOMAIN_COL = "domain_number"
CS_LRR_COL    = "mean_lrr"
CS_STD_COL    = "std_lrr"

# ── Transect-mode inputs ─────────────────────────────────────
# Point to your transect_lrr_full.csv files for each period
TRANSECT_CSV_1984_2004 = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\1984_2004\transect_lrr_full.csv"
TRANSECT_CSV_2004_2024 = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\2004_2024\transect_lrr_full.csv"

# Column names in your transect CSV (transect_lrr_full.csv)
T_TRANSECT_ID_COL = "transect_id"    # string IDs — converted to sequential int internally
T_ALONG_COAST_COL = None             # not in CSV — derived from domain position automatically
T_DOMAIN_COL      = "domain_number"
T_LRR_COL         = "lrr_m_yr"
T_STD_COL         = "unc_m_yr"       # uncertainty column; set to None to skip

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
COMPARE_WINDOWS_KM = [2.5, 3.0, 3.5, 4.0]   # 5, 6, 7, 8 domains

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
# Fill in values matching your transect IDs or along-coast distances (m)
# once you know them. Any entry left as None is silently skipped.
T_TOWN_SPANS = {
    "Buxton":      (None, None),
    "Avon":        (None, None),
    "Tri-Village": (None, None),
}
T_VILLAGE_LINES = {"Salvo": None, "Waves": None, "Rodanthe": None}
T_PIERS         = {"Avon Pier": None, "Rodanthe Pier": None}
T_GROINS        = {"Buxton Groin": None}
T_WIMBLE_SHOALS = (None, None)

# ── Annotation colors ────────────────────────────────────────
C_TOWN_SPAN    = "#90AFC5"
C_WIMBLE       = "#E0A800"
C_VILLAGE_LINE = "0.40"
C_PIER         = "#1565C0"
C_GROIN        = "#B71C1C"

# ── Period colors ────────────────────────────────────────────
C_CS_1984 = "#1F4E79"
C_CS_2004 = "#833C00"
C_WINDOWS  = ["#2ca02c", "#ff7f0e", "#9467bd", "#e377c2"]

# ── Output ───────────────────────────────────────────────────
OUTPUT_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\shoreline_change_smoothing\coastsat_smoothing_domainmethod"

# ============================================================
# IMPORTS
# ============================================================
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
from statsmodels.nonparametric.smoothers_lowess import lowess
import warnings
warnings.filterwarnings("ignore")

os.makedirs(OUTPUT_DIR, exist_ok=True)

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
    """Load domain-averaged LRR summary CSV. Returns standardised DataFrame or None."""
    if path is None or not os.path.exists(path):
        print(f"  Domain CSV ({period_label}): SKIPPED — not found: {path}")
        return None
    df = pd.read_csv(path)
    for col in [CS_DOMAIN_COL, CS_LRR_COL, CS_STD_COL]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=[CS_DOMAIN_COL, CS_LRR_COL])
    df[CS_DOMAIN_COL] = df[CS_DOMAIN_COL].astype(int)
    df = df[(df[CS_DOMAIN_COL] >= DOMAIN_MIN) & (df[CS_DOMAIN_COL] <= DOMAIN_MAX)]
    df = df[[CS_DOMAIN_COL, CS_LRR_COL, CS_STD_COL]].rename(columns={
        CS_DOMAIN_COL: "domain",
        CS_LRR_COL:    "cs_lrr",
        CS_STD_COL:    "cs_std",
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
                  loc="best", ncol=2 if is_bottom else 1)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


def plot_domain_overview(d1984, d2004, out_path, suffix=""):
    _domain_two_panel(
        d1984, d2004, show_raw=True, out_path=out_path,
        suptitle=(f"CoastSat Shoreline Change Rate — LOESS Smoothed\n"
                  f"Hatteras Island, NC{suffix}"),
    )


def plot_domain_smoothed_only(d1984, d2004, out_path, suffix=""):
    _domain_two_panel(
        d1984, d2004, show_raw=False, out_path=out_path,
        suptitle=(f"CoastSat Shoreline Change Rate — Smoothed Only\n"
                  f"Hatteras Island, NC{suffix}"),
    )


def plot_domain_combined(d1984, d2004, out_path, suffix=""):
    fig, ax = plt.subplots(figsize=(16, 6))
    fig.suptitle(
        f"CoastSat Shoreline Change Rate — Both Periods\nHatteras Island, NC{suffix}",
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
    ax.legend(handles=handles, fontsize=9, framealpha=0.95, loc="best", ncol=2)
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
        f"LOESS Bandwidth Sensitivity — {period_label}\nHatteras Island, NC{suffix}",
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
                  loc="best", ncol=2 if is_bottom else 1)
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
        f"LOESS Window Comparison — Domain Space\nHatteras Island, NC{suffix}",
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
                  loc="lower right", ncol=2 if is_bottom else 1)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")

# ============================================================
# TRANSECT-SPACE FIGURES  (produced only in "transect" mode)
# ============================================================

def plot_transect_overview(t1984, t2004, out_path):
    """Raw transect scatter + LOESS curve in transect / along-coast space."""
    configs = [(t1984, "1984–2004", C_CS_1984),
               (t2004, "2004–2024", C_CS_2004)]
    xlabel  = ("Along-coast Distance (m)" if TRANSECT_X_AXIS == "along_coast_m"
               else "Transect ID (sequential)")
    fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True)
    fig.suptitle(
        f"CoastSat Transect LRR — LOESS Smoothed\n"
        f"Hatteras Island, NC  (x = {xlabel})",
        fontsize=14, fontweight="bold", y=1.01,
    )
    for i, (ax, (df, period, color)) in enumerate(zip(axes, configs)):
        is_bottom = (i == len(configs) - 1)
        if df is None:
            ax.text(0.5, 0.5, f"No data — {period}",
                    transform=ax.transAxes, ha="center", fontsize=12)
            continue
        x = df["_x_smooth"].values
        ax.scatter(x, df["lrr"], color=color, s=4, alpha=0.20,
                   zorder=1, label="Raw transect LRR")
        ax.plot(x, df["lrr_smooth"], color=color, lw=2.5,
                zorder=3, label="LOESS smoothed")
        ax.set_ylabel("Rate (m/yr)", fontsize=10, fontweight="bold")
        ax.set_title(period, fontsize=11, fontweight="bold", loc="left", pad=5)
        style_transect_axis(ax, x, is_bottom)
        add_transect_annotations(ax)
        handles, _ = ax.get_legend_handles_labels()
        if is_bottom:
            handles += annotation_legend_handles()
        ax.legend(handles=handles, fontsize=8.5, framealpha=0.95,
                  loc="best", ncol=2 if is_bottom else 1)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


def plot_transect_window_comparison(t1984, t2004, out_path):
    """Window sensitivity in transect space — all km windows overlaid, both periods."""
    configs = [(t1984, "1984–2004", C_CS_1984),
               (t2004, "2004–2024", C_CS_2004)]
    xlabel  = ("Along-coast Distance (m)" if TRANSECT_X_AXIS == "along_coast_m"
               else "Transect ID (sequential)")
    fig, axes = plt.subplots(2, 1, figsize=(16, 11), sharex=True)
    fig.suptitle(
        f"LOESS Window Comparison — Transect Space\nHatteras Island, NC  (x = {xlabel})",
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
                  loc="lower right", ncol=2 if is_bottom else 1)
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
        "LOESS Window Comparison — Transect-Smoothed, Domain Space\n"
        "Hatteras Island, NC",
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
                  loc="lower right", ncol=2 if is_bottom else 1)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


# ============================================================
# TRANSECT VALUES vs DOMAIN AVERAGES
# Shows raw transect LRRs as a dotted line alongside domain-averaged
# LRRs as a solid line — both periods, domain space.
# ============================================================

def plot_transect_vs_domain_averages(t1984, t2004, d1984, d2004, out_path):
    """
    Overlays two series per period:
      - Transect LRR values as a dotted line (connected in along-coast order,
        x mapped to domain-equivalent position = along_coast_m / DOMAIN_SPACING_M)
      - Domain-averaged LRR as a solid line with circle markers at domain centroids

    Allows the advisor to see how much within-domain variability the
    domain-averaging step collapses before smoothing.
    """
    configs   = [(t1984, d1984, "1984\u20132004", C_CS_1984),
                 (t2004, d2004, "2004\u20132024", C_CS_2004)]
    fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True)
    fig.suptitle(
        "CoastSat LRR \u2014 Transect Values vs Domain Averages\n"
        "Hatteras Island, NC",
        fontsize=14, fontweight="bold", y=1.01,
    )

    for i, (ax, (t_df, d_df, period, color)) in enumerate(zip(axes, configs)):
        is_bottom = (i == len(configs) - 1)
        if t_df is None or d_df is None:
            ax.text(0.5, 0.5, f"No data \u2014 {period}",
                    transform=ax.transAxes, ha="center", fontsize=12)
            continue

        # Transect LRRs as dotted line — x = along_coast_m / domain_spacing
        # gives a continuous domain-equivalent x so both series share the axis
        x_transect = t_df["along_coast_m"].values / DOMAIN_SPACING_M
        ax.plot(x_transect, t_df["lrr"].values,
                color=color, lw=0.9, ls=":", alpha=0.65,
                label="Individual transect LRR")

        # Domain averages as solid line with circle markers
        ax.plot(d_df["domain"], d_df["cs_lrr"],
                color=color, lw=2.0, ls="-",
                marker="o", ms=4, markerfacecolor="white", markeredgewidth=1.5,
                label="Domain average LRR")

        ax.set_ylabel("Rate (m/yr)", fontsize=10, fontweight="bold")
        ax.set_title(period, fontsize=11, fontweight="bold", loc="left", pad=5)
        style_domain_axis(ax, is_bottom)
        add_domain_annotations(ax)
        handles, _ = ax.get_legend_handles_labels()
        if is_bottom:
            handles += annotation_legend_handles()
        ax.legend(handles=handles, fontsize=8.5, framealpha=0.95,
                  loc="upper center", ncol=3 if is_bottom else 2)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")

# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 65)
    print("CoastSat LRR Smoothing — Hatteras Island")
    print(f"  Mode        : {SMOOTHING_MODE}")
    if SMOOTHING_MODE == "transect":
        print(f"  X-axis      : {TRANSECT_X_AXIS}")
    print(f"  Window      : {LOESS_WINDOW_KM} km  "
          f"({int(round(LOESS_WINDOW_KM * 1000 / DOMAIN_SPACING_M))} domains)")
    print(f"  Compare     : {COMPARE_WINDOWS_KM} km  "
          f"({[int(round(w * 1000 / DOMAIN_SPACING_M)) for w in COMPARE_WINDOWS_KM]} domains)")
    print("=" * 65)

    # ── Load & smooth ──────────────────────────────────────
    print("\nLoading and smoothing data...")

    if SMOOTHING_MODE == "domain":
        d1984 = load_domain_csv(DOMAIN_CSV_1984_2004, "1984–2004")
        d2004 = load_domain_csv(DOMAIN_CSV_2004_2024, "2004–2024")
        if d1984 is not None:
            d1984 = smooth_domain_df(d1984)
        if d2004 is not None:
            d2004 = smooth_domain_df(d2004)
        t1984 = t2004 = None

    else:  # transect mode
        t1984_raw = load_transect_csv(TRANSECT_CSV_1984_2004, "1984–2004")
        t2004_raw = load_transect_csv(TRANSECT_CSV_2004_2024, "2004–2024")
        t1984 = smooth_transect_df(t1984_raw) if t1984_raw is not None else None
        t2004 = smooth_transect_df(t2004_raw) if t2004_raw is not None else None
        d1984 = aggregate_to_domains(t1984) if t1984 is not None else None
        d2004 = aggregate_to_domains(t2004) if t2004 is not None else None

    mode_suffix = (f"  [{SMOOTHING_MODE} mode, "
                   f"{LOESS_WINDOW_KM} km / "
                   f"{int(round(LOESS_WINDOW_KM * 1000 / DOMAIN_SPACING_M))} domain window]")

    # ── Domain-space figures ───────────────────────────────
    print("\nGenerating domain-space figures...")

    plot_domain_overview(d1984, d2004,
        os.path.join(OUTPUT_DIR, "overview_smoothed.png"), mode_suffix)

    plot_domain_smoothed_only(d1984, d2004,
        os.path.join(OUTPUT_DIR, "smoothed_only_comparison.png"), mode_suffix)

    plot_domain_combined(d1984, d2004,
        os.path.join(OUTPUT_DIR, "combined_periods.png"), mode_suffix)

    if d1984 is not None:
        plot_domain_sensitivity(d1984, "1984–2004", C_CS_1984,
            os.path.join(OUTPUT_DIR, "smoothing_sensitivity_1984_2004.png"), mode_suffix)

    if d2004 is not None:
        plot_domain_sensitivity(d2004, "2004–2024", C_CS_2004,
            os.path.join(OUTPUT_DIR, "smoothing_sensitivity_2004_2024.png"), mode_suffix)

    # In transect mode window_comparison uses transect-smoothed curves
    # aggregated to domain space. In domain mode smooths domain means directly.
    if SMOOTHING_MODE == "transect":
        plot_transect_windows_domain_space(t1984, t2004,
            os.path.join(OUTPUT_DIR, "window_comparison.png"))
    else:
        plot_domain_window_comparison(d1984, d2004,
            os.path.join(OUTPUT_DIR, "window_comparison.png"), mode_suffix)

    # ── Transect-space figures ─────────────────────────────
    if SMOOTHING_MODE == "transect":
        print("\nGenerating transect-space figures...")
        plot_transect_overview(t1984, t2004,
            os.path.join(OUTPUT_DIR, "transect_smoothed_overview.png"))
        plot_transect_window_comparison(t1984, t2004,
            os.path.join(OUTPUT_DIR, "transect_window_comparison.png"))
        # Domain-averaged smoothing in same format for direct comparison
        plot_domain_window_comparison(d1984, d2004,
            os.path.join(OUTPUT_DIR, "window_comparison_domain_averaged.png"),
            "  [domain-averaged smoothing]")
        # Transect values vs domain averages overlay
        plot_transect_vs_domain_averages(t1984, t2004, d1984, d2004,
            os.path.join(OUTPUT_DIR, "transect_vs_domain_averages.png"))

    # ── Export tables ──────────────────────────────────────
    print("\nExporting tables...")

    # Domain-level table (both modes)
    table_parts = []
    for df, label in [(d1984, "1984_2004"), (d2004, "2004_2024")]:
        if df is None:
            continue
        t = df[["domain", "cs_lrr", "cs_std", "cs_lrr_smooth"]].copy()
        t.columns = ["domain"] + [f"{c}_{label}"
                                   for c in ["cs_lrr", "cs_std", "cs_lrr_smooth"]]
        table_parts.append(t)
    if table_parts:
        from functools import reduce
        table = reduce(lambda a, b: a.merge(b, on="domain", how="outer"), table_parts)
        fname = f"coastsat_{SMOOTHING_MODE}_smoothed_table.csv"
        table.sort_values("domain").to_csv(os.path.join(OUTPUT_DIR, fname), index=False)
        print(f"  Saved: {fname}")

    # Transect-level tables (transect mode only)
    if SMOOTHING_MODE == "transect":
        for df, label in [(t1984, "1984_2004"), (t2004, "2004_2024")]:
            if df is None:
                continue
            cols = ["transect_id", "along_coast_m", "domain", "lrr", "lrr_smooth"]
            if "lrr_std" in df.columns:
                cols.insert(4, "lrr_std")
            fname = f"coastsat_transect_lrr_{label}.csv"
            df[cols].to_csv(os.path.join(OUTPUT_DIR, fname), index=False)
            print(f"  Saved: {fname}")

    print("\n" + "=" * 65)
    print("Done!")
    print(f"  Mode: {SMOOTHING_MODE}  |  "
          f"Window: {LOESS_WINDOW_KM} km  "
          f"({int(round(LOESS_WINDOW_KM * 1000 / DOMAIN_SPACING_M))} domains)")
    print("  Domain-space figures:")
    print("    overview_smoothed.png")
    print("    smoothed_only_comparison.png")
    print("    combined_periods.png")
    print("    smoothing_sensitivity_1984_2004.png")
    print("    smoothing_sensitivity_2004_2024.png")
    print("    window_comparison.png")
    if SMOOTHING_MODE == "transect":
        print("  Transect-space figures:")
        print("    transect_smoothed_overview.png")
        print("    transect_window_comparison.png")
    print("=" * 65)


if __name__ == "__main__":
    main()
