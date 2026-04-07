"""
CoastSat LRR Smoothing — Hatteras Island
=========================================
Applies LOESS smoothing to CoastSat shoreline change rates across two
time periods and produces publication-ready figures.

Outputs (saved to OUTPUT_DIR)
------------------------------
  overview_smoothed.png          – 2-panel both periods, raw + LOESS overlay
  smoothed_only_comparison.png   – 2-panel both periods, smoothed lines only
  combined_periods.png           – single panel combining both periods (smoothed)
  smoothing_sensitivity_*.png    – 3-panel bandwidth sensitivity per period
  window_comparison.png          – NEW: raw + 3 smoothing windows overlaid,
                                    both periods, for window selection

Smoothing method: LOESS (locally weighted scatterplot smoothing)
  - Applied independently to each period's CoastSat series
  - Preserves large-scale spatial patterns while removing per-domain noise
"""

# ============================================================
# CONFIG
# ============================================================

# --- CoastSat CSVs ---
COASTSAT_CSV_1984_2004 = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\1984_2004\domain_lrr_summary.csv"
COASTSAT_CSV_2004_2024 = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\2004_2024\domain_lrr_summary.csv"

CS_DOMAIN_COL = "domain_number"
CS_LRR_COL    = "mean_lrr"
CS_STD_COL    = "std_lrr"

# --- Domain range ---
DOMAIN_MIN = 1
DOMAIN_MAX = 90

# --- LOESS bandwidth (fraction of data used per local fit) ---
# 0.10 = ~9 domains  → more local, preserves more variation
# 0.167 = ~15 domains → recommended default (1.5km smoothing window)
# 0.20 = ~18 domains → smoother, loses finer spatial patterns
LOESS_FRAC = 0.167

# --- Window comparison: domain counts to test in window_comparison.png ---
# Each value is the number of CASCADE domains (~500 m each) used in the
# local LOESS fit.  Fracs are computed as n / (DOMAIN_MAX - DOMAIN_MIN + 1).
#   5  domains → frac ≈ 0.056  (2.5 km window)
#  10  domains → frac ≈ 0.111  (5.0 km window)
#  15  domains → frac ≈ 0.167  (7.5 km window) ← matches LOESS_FRAC default
COMPARE_WINDOWS_DOMAINS = [5, 10, 15]   # ← edit here (in domains)

# --- Geographic annotations ---

# Community shaded spans: label → (domain_lo, domain_hi), inclusive
TOWN_SPANS = {
    "Buxton":      (7,  8),
    "Avon":        (21, 31),
    "Tri-Village": (68, 83),   # Salvo / Waves / Rodanthe combined area
}

# Individual village center lines within Tri-Village
VILLAGE_LINES = {
    "Salvo":    69,
    "Waves":    74,
    "Rodanthe": 80,
}

# Pier reference lines
PIERS = {
    "Avon Pier":      26,
    "Rodanthe Pier":  79,
}

# Groin reference lines
GROINS = {
    "Buxton Groin": 6,
}

# Wimble Shoals offshore shoal influence zone: (domain_lo, domain_hi)
# Spans domains 60–74; overlaps with the southern portion of Tri-Village.
WIMBLE_SHOALS = (60, 74)

# --- Annotation colors ---
C_TOWN_SPAN    = "#90AFC5"   # steel blue  — community spans
C_WIMBLE       = "#E0A800"   # amber       — Wimble Shoals fill
C_VILLAGE_LINE = "0.40"      # dark gray   — Salvo / Waves / Rodanthe lines
C_PIER         = "#1565C0"   # medium blue — pier lines
C_GROIN        = "#B71C1C"   # dark red    — groin lines

# --- Output ---
OUTPUT_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\shoreline_change_smoothing\coastsat_smoothing_1984_2024_TEST 5 labels"

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

# Color palette
C_CS_1984 = "#1F4E79"   # dark blue  — 1984–2004
C_CS_2004 = "#833C00"   # dark brown — 2004–2024

# Colors for the three smoothing windows in the comparison figure.
# Raw CoastSat inherits the period color (faded); windows are distinct.
C_WINDOWS = ["#2ca02c", "#ff7f0e", "#9467bd"]   # green, orange, purple

# ============================================================
# DATA LOADING
# ============================================================

def load_coastsat(path, domain_col, lrr_col, std_col, period_label):
    if path is None or not os.path.exists(path):
        print(f"  CoastSat ({period_label}): SKIPPED — not found: {path}")
        return None
    df = pd.read_csv(path)
    df[domain_col] = pd.to_numeric(df[domain_col], errors="coerce")
    df[lrr_col]    = pd.to_numeric(df[lrr_col],    errors="coerce")
    df[std_col]    = pd.to_numeric(df[std_col],     errors="coerce")
    df = df.dropna(subset=[domain_col, lrr_col])
    df[domain_col] = df[domain_col].astype(int)
    df = df[(df[domain_col] >= DOMAIN_MIN) & (df[domain_col] <= DOMAIN_MAX)]
    df = df[[domain_col, lrr_col, std_col]].rename(columns={
        domain_col: "domain",
        lrr_col:    "cs_lrr",
        std_col:    "cs_std",
    }).sort_values("domain").reset_index(drop=True)
    print(f"  CoastSat ({period_label}): {len(df)} domains  "
          f"range {df['cs_lrr'].min():+.2f} to {df['cs_lrr'].max():+.2f} m/yr")
    return df


# ============================================================
# SMOOTHING
# ============================================================

def apply_loess(domains, values, frac=LOESS_FRAC):
    """Apply LOESS smoothing. Returns smoothed values at same domain positions."""
    valid = ~np.isnan(values)
    if valid.sum() < 5:
        return values.copy()
    smoothed = np.full_like(values, np.nan)
    result = lowess(values[valid], domains[valid], frac=frac, return_sorted=True)
    smoothed[valid] = np.interp(domains[valid], result[:, 0], result[:, 1])
    return smoothed


def add_smoothed_columns(df, frac=LOESS_FRAC):
    """Add LOESS-smoothed LRR column to a CoastSat dataframe."""
    df = df.copy()
    d = df["domain"].values.astype(float)
    df["cs_lrr_smooth"] = apply_loess(d, df["cs_lrr"].values, frac)
    return df


def domains_to_frac(n_domains):
    """
    Convert a window size in CASCADE domains to a LOESS frac value.

    Always divides by the total island domain range (DOMAIN_MAX - DOMAIN_MIN + 1)
    so that fracs are consistent regardless of how many domains have valid data
    in a given CSV.  This ensures:
      5  domains → frac ≈ 0.056
      10 domains → frac ≈ 0.111
      15 domains → frac ≈ 0.167  (matches the LOESS_FRAC default)
    """
    return n_domains / (DOMAIN_MAX - DOMAIN_MIN + 1)


# ============================================================
# SHARED AXIS HELPERS
# ============================================================

def add_annotations(ax):
    """
    Add all geographic reference annotations to an axis.

    Layer order (bottom → top):
      1. Wimble Shoals influence zone  (hatched amber fill, bottom label)
      2. Community shaded spans        (steel-blue fill, top labels)
      3. Village center lines          (dashed gray,  y=0.88)
      4. Pier lines                    (dash-dot blue, y=0.76, rotated)
      5. Groin lines                   (dotted red,    y=0.76, rotated)

    All label y-positions use blended axes-fraction coordinates so they
    stay fixed relative to the panel height regardless of data range.
    """
    # Blended transform: data x, axes-fraction y
    trans = blended_transform_factory(ax.transData, ax.transAxes)

    # ------------------------------------------------------------------
    # 1. Wimble Shoals influence zone
    #    Placed first so community spans render on top of it.
    #    Label at the bottom of the panel to avoid crowding the top tier.
    # ------------------------------------------------------------------
    wlo, whi = WIMBLE_SHOALS
    ax.axvspan(wlo - 0.5, whi + 0.5,
               color=C_WIMBLE, alpha=0.10, zorder=0,
               hatch="///", edgecolor=C_WIMBLE, linewidth=0)
    ax.text((wlo + whi) / 2.0, 0.04,
            "Wimble Shoals\ninfluence", transform=trans,
            ha="center", va="bottom", fontsize=7, color="#7A5800",
            style="italic",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.80))

    # ------------------------------------------------------------------
    # 2. Community / town spans
    #    Labels at y=0.90 (one step below the panel frac-label tier at 0.97).
    # ------------------------------------------------------------------
    for span_label, (d_lo, d_hi) in TOWN_SPANS.items():
        ax.axvspan(d_lo - 0.5, d_hi + 0.5,
                   color=C_TOWN_SPAN, alpha=0.14, zorder=0)
        ax.text((d_lo + d_hi) / 2.0, 0.90,
                span_label, transform=trans,
                ha="center", va="top", fontsize=8, color="0.25",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

    # ------------------------------------------------------------------
    # 3. Village center lines (within Tri-Village span)
    #    Dashed gray; labels at y=0.79 so they sit clearly below span labels.
    # ------------------------------------------------------------------
    for vname, dom in VILLAGE_LINES.items():
        ax.axvline(dom, color=C_VILLAGE_LINE, lw=0.9, ls="--",
                   alpha=0.65, zorder=1)
        ax.text(dom, 0.79, vname, transform=trans,
                ha="center", va="top", fontsize=7.5, color="0.30",
                bbox=dict(boxstyle="round,pad=0.15", fc="white",
                          ec="none", alpha=0.80))

    # ------------------------------------------------------------------
    # 4. Pier lines
    #    Dash-dot blue; rotated labels at y=0.74.
    #    Avon Pier (26) sits inside the Avon span — offset label downward
    #    so it clears the span label above.
    # ------------------------------------------------------------------
    for pname, dom in PIERS.items():
        ax.axvline(dom, color=C_PIER, lw=1.0, ls="-.", alpha=0.80, zorder=2)
        ax.text(dom, 0.74, pname, transform=trans,
                ha="center", va="top", fontsize=7, color=C_PIER,
                rotation=90,
                bbox=dict(boxstyle="round,pad=0.15", fc="white",
                          ec="none", alpha=0.80))

    # ------------------------------------------------------------------
    # 5. Groin lines
    #    Dotted dark-red; rotated labels at y=0.74.
    #    Buxton Groin (6) is just left of the Buxton span (7–8).
    # ------------------------------------------------------------------
    for gname, dom in GROINS.items():
        ax.axvline(dom, color=C_GROIN, lw=1.1, ls=":", alpha=0.85, zorder=2)
        ax.text(dom, 0.74, gname, transform=trans,
                ha="center", va="top", fontsize=7, color=C_GROIN,
                rotation=90,
                bbox=dict(boxstyle="round,pad=0.15", fc="white",
                          ec="none", alpha=0.80))


def annotation_legend_handles():
    """
    Return proxy artists explaining the annotation layer types.
    Append these to a plot's legend handle list so readers can decode
    all reference marks without scanning every label individually.
    """
    return [
        Patch(fc=C_TOWN_SPAN,  alpha=0.30, label="Community"),
        Patch(fc=C_WIMBLE,     alpha=0.25, hatch="///",
              edgecolor=C_WIMBLE, linewidth=0, label="Wimble Shoals influence"),
        Line2D([0], [0], color=C_VILLAGE_LINE, lw=0.9, ls="--",
               label="Village center"),
        Line2D([0], [0], color=C_PIER,  lw=1.0, ls="-.", label="Pier"),
        Line2D([0], [0], color=C_GROIN, lw=1.1, ls=":",  label="Groin"),
    ]


def style_domain_axis(ax):
    ax.set_xlim(DOMAIN_MIN - 0.5, DOMAIN_MAX + 0.5)
    ax.set_xticks(range(DOMAIN_MIN, DOMAIN_MAX + 1, 10))
    ax.set_xticks(range(DOMAIN_MIN, DOMAIN_MAX + 1, 5), minor=True)
    ax.set_xlabel("CASCADE Model Domain (500 m alongshore)",
                  fontsize=11, fontweight="bold")
    ax.axhline(0, color="black", lw=1.1, ls="--", alpha=0.55)
    ax.annotate("Accretion ▲", xy=(DOMAIN_MAX + 0.5, 0.15),
                xycoords=("data", "axes fraction"),
                ha="right", va="bottom", fontsize=8, color="0.4")
    ax.annotate("Erosion ▼", xy=(DOMAIN_MAX + 0.5, 0.10),
                xycoords=("data", "axes fraction"),
                ha="right", va="top", fontsize=8, color="0.4")


# ============================================================
# FIGURE 1 — PRIMARY: Raw + LOESS overlay, both periods
# ============================================================

def plot_overview_smoothed(cs_1984, cs_2004, out_path):
    """
    2-panel figure: raw CoastSat points (faded) + LOESS overlay (bold).
    """
    fig, axes = plt.subplots(2, 1, figsize=(16, 11), sharex=True)
    fig.suptitle("CoastSat Shoreline Change Rates — Hatteras Island, NC\n"
                 "Raw (faded) + LOESS smoothed (bold)",
                 fontsize=14, fontweight="bold", y=1.01)

    configs = [
        (axes[0], cs_1984, "1984–2004", C_CS_1984),
        (axes[1], cs_2004, "2004–2024", C_CS_2004),
    ]

    for i, (ax, df, label, color) in enumerate(configs):
        is_bottom = (i == len(configs) - 1)

        if df is None:
            ax.text(0.5, 0.5, f"No data — {label}",
                    transform=ax.transAxes, ha="center", fontsize=12)
            continue

        m = add_smoothed_columns(df)
        d = m["domain"]

        # Raw (faded, thin)
        ax.plot(d, m["cs_lrr"], color=color, lw=1.0, alpha=0.30,
                marker="o", ms=2.5, zorder=2, label=f"Raw {label}")

        # ±1 std shading (very faint)
        ax.fill_between(d,
                        m["cs_lrr"] - m["cs_std"],
                        m["cs_lrr"] + m["cs_std"],
                        color=color, alpha=0.08, zorder=1)

        # LOESS smoothed (bold)
        ax.plot(d, m["cs_lrr_smooth"], color=color, lw=2.8,
                label=f"LOESS smoothed (frac={LOESS_FRAC})", zorder=4)

        ax.set_ylabel("Shoreline Change Rate (m/yr)", fontsize=11, fontweight="bold")
        ax.set_title(label, fontsize=12, fontweight="bold", loc="left", pad=6)
        style_domain_axis(ax)
        add_annotations(ax)

        # Annotation legend proxies in the bottom panel only
        handles, lbls = ax.get_legend_handles_labels()
        if is_bottom:
            handles = handles + annotation_legend_handles()
        ax.legend(handles=handles, fontsize=8.5, framealpha=0.95,
                  loc="lower right", ncol=2 if is_bottom else 1)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


# ============================================================
# FIGURE 2 — SMOOTHED ONLY (clean, for presentations)
# ============================================================

def plot_smoothed_only(cs_1984, cs_2004, out_path):
    """
    2-panel: LOESS smoothed lines only, no raw data.
    Cleanest version for presentations or dissertation figures.
    """
    fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True)
    fig.suptitle("CoastSat Shoreline Change Rates — Hatteras Island, NC\n"
                 f"LOESS smoothed (frac={LOESS_FRAC})",
                 fontsize=14, fontweight="bold", y=1.01)

    configs = [
        (axes[0], cs_1984, "1984–2004", C_CS_1984),
        (axes[1], cs_2004, "2004–2024", C_CS_2004),
    ]

    for i, (ax, df, label, color) in enumerate(configs):
        is_bottom = (i == len(configs) - 1)

        if df is None:
            continue
        m = add_smoothed_columns(df)
        d = m["domain"]

        ax.plot(d, m["cs_lrr_smooth"], color=color, lw=3.0,
                label=f"CoastSat {label} (smoothed)")

        # ±1 std shading around smoothed line
        ax.fill_between(d,
                        m["cs_lrr_smooth"] - m["cs_std"],
                        m["cs_lrr_smooth"] + m["cs_std"],
                        color=color, alpha=0.12, label="±1 std")

        ax.set_ylabel("Shoreline Change Rate (m/yr)", fontsize=11, fontweight="bold")
        ax.set_title(label, fontsize=12, fontweight="bold", loc="left", pad=6)
        style_domain_axis(ax)
        add_annotations(ax)

        handles, lbls = ax.get_legend_handles_labels()
        if is_bottom:
            handles = handles + annotation_legend_handles()
        ax.legend(handles=handles, fontsize=8.5, framealpha=0.95,
                  loc="lower right", ncol=2 if is_bottom else 1)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


# ============================================================
# FIGURE 3 — COMBINED: Both periods on one panel (smoothed)
# ============================================================

def plot_combined_periods(cs_1984, cs_2004, out_path):
    """
    Single panel: both periods of smoothed CoastSat on one axis.
    Good for directly comparing the two periods.
    """
    fig, ax = plt.subplots(figsize=(16, 6))

    for df, label, color in [
        (cs_1984, "1984–2004", C_CS_1984),
        (cs_2004, "2004–2024", C_CS_2004),
    ]:
        if df is None:
            continue
        m = add_smoothed_columns(df)
        d = m["domain"]
        ax.plot(d, m["cs_lrr_smooth"], color=color, lw=2.5,
                label=f"CoastSat {label}")
        ax.fill_between(d,
                        m["cs_lrr_smooth"] - m["cs_std"] * 0.5,
                        m["cs_lrr_smooth"] + m["cs_std"] * 0.5,
                        color=color, alpha=0.10)

    ax.set_ylabel("Shoreline Change Rate (m/yr)", fontsize=12, fontweight="bold")
    ax.set_title("CoastSat Shoreline Change Rates — Both Periods (LOESS smoothed)",
                 fontsize=13, fontweight="bold", pad=12)
    style_domain_axis(ax)
    add_annotations(ax)

    handles, lbls = ax.get_legend_handles_labels()
    ax.legend(handles=handles + annotation_legend_handles(),
              fontsize=9, framealpha=0.95, ncol=2)

    fig.text(0.5, -0.04,
             f"Domain-averaged LRR smoothed with LOESS (frac={LOESS_FRAC}). "
             f"Shading = ±0.5 std of CoastSat transects per domain.",
             ha="center", fontsize=8, color="0.4", style="italic")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


# ============================================================
# FIGURE 4 — SMOOTHING SENSITIVITY (frac comparison)
# ============================================================

def plot_smoothing_sensitivity(df, period_label, out_path,
                               window_domains=COMPARE_WINDOWS_DOMAINS):
    """
    3-panel showing the effect of each LOESS window on CoastSat data.
    One smoothed line per panel so individual window behavior is clear.
    Fracs are derived from COMPARE_WINDOWS_DOMAINS / 90 domains, matching
    exactly what plot_window_comparison uses.
    """
    n_total = DOMAIN_MAX - DOMAIN_MIN + 1   # 90 domains
    fracs = [domains_to_frac(n) for n in window_domains]
    panel_labels = [
        f"frac={f:.3f}  ({n} domains, {n*0.5:.1f} km window)"
        for n, f in zip(window_domains, fracs)
    ]
    # Mark whichever window matches LOESS_FRAC as the current default
    panel_labels = [
        lbl + "  ← current default" if abs(f - LOESS_FRAC) < 1e-4 else lbl
        for lbl, f in zip(panel_labels, fracs)
    ]

    color = C_CS_1984 if "1984" in period_label else C_CS_2004

    fig, axes = plt.subplots(3, 1, figsize=(16, 13), sharex=True, sharey=True)
    fig.suptitle(f"LOESS Smoothing Sensitivity — CoastSat {period_label}\n"
                 f"Effect of bandwidth choice",
                 fontsize=13, fontweight="bold", y=1.01)

    for j, (ax, frac, plabel) in enumerate(zip(axes, fracs, panel_labels)):
        is_bottom = (j == len(fracs) - 1)

        m = add_smoothed_columns(df, frac=frac)
        d = m["domain"]

        # Raw faded
        ax.plot(d, m["cs_lrr"], color=color, lw=0.8, alpha=0.25,
                marker="o", ms=2, zorder=1, label="Raw")

        # Smoothed
        ax.plot(d, m["cs_lrr_smooth"], color=color, lw=2.5,
                label="LOESS smoothed", zorder=3)

        ax.text(0.01, 0.97, plabel,
                transform=ax.transAxes, fontsize=8.5, va="top",
                bbox=dict(boxstyle="round", fc="white", alpha=0.88, ec="0.7"))

        ax.set_ylabel("Rate (m/yr)", fontsize=10, fontweight="bold")
        style_domain_axis(ax)
        add_annotations(ax)

        handles, lbls = ax.get_legend_handles_labels()
        if is_bottom:
            handles = handles + annotation_legend_handles()
        # Use loc="best" so matplotlib picks whitespace automatically —
        # the data range differs between periods and panel corners vary.
        ax.legend(handles=handles, fontsize=8.5, framealpha=0.95,
                  loc="best", ncol=2 if is_bottom else 1)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


# ============================================================
# FIGURE 5 — WINDOW COMPARISON: All 3 windows overlaid (NEW)
# ============================================================

def plot_window_comparison(cs_1984, cs_2004, out_path,
                           window_domains=COMPARE_WINDOWS_DOMAINS):
    """
    2-panel figure (one per period) showing the raw CoastSat LRR (faded)
    overlaid by LOESS-smoothed lines for each window size in
    COMPARE_WINDOWS_DOMAINS.  All smoothed lines share one panel so spatial
    patterns and differences between window choices are directly visible.

    Window sizes are specified in CASCADE domains and converted to LOESS fracs
    using the actual number of valid domains in each dataset.

    Parameters
    ----------
    cs_1984, cs_2004 : pd.DataFrame or None
        Loaded CoastSat data for each period.
    out_path : str
        Full path for the output PNG.
    window_domains : list of int
        Number of CASCADE domains for each smoothing window to compare.
        Defaults to COMPARE_WINDOWS_DOMAINS from CONFIG.
    """
    configs = [
        (cs_1984, "1984–2004", C_CS_1984),
        (cs_2004, "2004–2024", C_CS_2004),
    ]

    # Build figure — 2 rows (periods), 1 column; shared x-axis
    fig, axes = plt.subplots(2, 1, figsize=(16, 11), sharex=True)
    fig.suptitle(
        "CoastSat Shoreline Change Rate — LOESS Window Comparison\n"
        "Hatteras Island, NC",
        fontsize=14, fontweight="bold", y=1.01,
    )

    for i, (ax, df, period_label, period_color) in enumerate(configs):
        is_bottom = (i == len(configs) - 1)

        if df is None:
            ax.text(0.5, 0.5, f"No data — {period_label}",
                    transform=ax.transAxes, ha="center", fontsize=12)
            ax.set_title(period_label, fontsize=12, fontweight="bold",
                         loc="left", pad=6)
            continue

        d = df["domain"]
        ax.plot(d, df["cs_lrr"],
                color=period_color, lw=1.0, alpha=0.25,
                marker="o", ms=2.5, zorder=2,
                label="Raw CoastSat LRR")

        # ±1 std envelope (very faint)
        ax.fill_between(d,
                        df["cs_lrr"] - df["cs_std"],
                        df["cs_lrr"] + df["cs_std"],
                        color=period_color, alpha=0.07, zorder=1)

        # --- One smoothed line per window size ---
        for n_dom, win_color in zip(window_domains, C_WINDOWS):
            frac = domains_to_frac(n_dom)
            smoothed = apply_loess(d.values.astype(float),
                                   df["cs_lrr"].values, frac=frac)
            km = n_dom * 0.5   # 500 m domains → km
            ax.plot(d, smoothed,
                    color=win_color, lw=2.5, zorder=4,
                    label=f"{n_dom}-domain window  ({km:.1f} km,  frac={frac:.3f})")

        # --- Panel formatting ---
        ax.set_ylabel("Shoreline Change Rate (m/yr)",
                      fontsize=11, fontweight="bold")
        ax.set_title(period_label, fontsize=12, fontweight="bold",
                     loc="left", pad=6)
        style_domain_axis(ax)
        add_annotations(ax)

        handles, lbls = ax.get_legend_handles_labels()
        if is_bottom:
            handles = handles + annotation_legend_handles()
        ax.legend(handles=handles, fontsize=8.5, framealpha=0.95,
                  loc="lower right", ncol=2 if is_bottom else 1)

    # Shared caption
    window_str = ", ".join(
        f"{n} domains ({n*0.5:.1f} km)" for n in window_domains
    )
    fig.text(
        0.5, -0.03,
        f"LOESS windows tested: {window_str}.  "
        "Raw CoastSat LRR shown faded for reference.",
        ha="center", fontsize=8, color="0.4", style="italic",
    )

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
    print(f"LOESS bandwidth: frac={LOESS_FRAC}")
    print(f"Window comparison: {COMPARE_WINDOWS_DOMAINS} domains")
    print("=" * 65)

    # Load
    print("\nLoading data...")
    cs_1984 = load_coastsat(COASTSAT_CSV_1984_2004, CS_DOMAIN_COL,
                            CS_LRR_COL, CS_STD_COL, "1984–2004")
    cs_2004 = load_coastsat(COASTSAT_CSV_2004_2024, CS_DOMAIN_COL,
                            CS_LRR_COL, CS_STD_COL, "2004–2024")

    # Generate figures
    print("\nGenerating figures...")

    plot_overview_smoothed(cs_1984, cs_2004,
        os.path.join(OUTPUT_DIR, "overview_smoothed.png"))

    plot_smoothed_only(cs_1984, cs_2004,
        os.path.join(OUTPUT_DIR, "smoothed_only_comparison.png"))

    plot_combined_periods(cs_1984, cs_2004,
        os.path.join(OUTPUT_DIR, "combined_periods.png"))

    if cs_1984 is not None:
        plot_smoothing_sensitivity(cs_1984, "1984–2004",
            os.path.join(OUTPUT_DIR, "smoothing_sensitivity_1984_2004.png"))

    if cs_2004 is not None:
        plot_smoothing_sensitivity(cs_2004, "2004–2024",
            os.path.join(OUTPUT_DIR, "smoothing_sensitivity_2004_2024.png"))

    # NEW: overlaid window comparison
    plot_window_comparison(cs_1984, cs_2004,
        os.path.join(OUTPUT_DIR, "window_comparison.png"))

    # Save smoothed data table
    table_parts = []
    for df, label in [(cs_1984, "1984_2004"), (cs_2004, "2004_2024")]:
        if df is None:
            continue
        m = add_smoothed_columns(df)
        cols = ["domain", "cs_lrr", "cs_std", "cs_lrr_smooth"]
        t = m[cols].copy()
        t.columns = ["domain"] + [f"{c}_{label}" for c in cols if c != "domain"]
        table_parts.append(t)

    if table_parts:
        from functools import reduce
        table = reduce(lambda a, b: a.merge(b, on="domain", how="outer"), table_parts)
        table = table.sort_values("domain")
        out_csv = os.path.join(OUTPUT_DIR, "coastsat_smoothed_table.csv")
        table.to_csv(out_csv, index=False)
        print(f"  Saved: coastsat_smoothed_table.csv")

    print("\n" + "=" * 65)
    print("Done! Key figures:")
    print("  overview_smoothed.png        ← raw + smoothed overlay, both periods")
    print("  smoothed_only_comparison.png ← clean version for presentations")
    print("  combined_periods.png         ← both periods on one panel")
    print("  smoothing_sensitivity_*.png  ← effect of bandwidth choice")
    print("  window_comparison.png        ← all 3 windows overlaid ← NEW")
    print("=" * 65)


if __name__ == "__main__":
    main()
