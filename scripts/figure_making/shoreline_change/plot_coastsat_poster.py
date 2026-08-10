"""
CoastSat Poster Figure — Hatteras Island
=========================================
Produces a single publication-quality poster figure showing LOESS-smoothed
CoastSat shoreline change rates for both calibration periods:
  - Period 1 (calibration): 1984–2004
  - Period 2 (test):        2004–2024

Smoothing method: LOESS with a fixed 10-domain (5.0 km) window.
Output: poster_two_periods_10dom.png  (saved to OUTPUT_DIR)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
from statsmodels.nonparametric.smoothers_lowess import lowess
import warnings
warnings.filterwarnings("ignore")

# ============================================================
# CONFIG  — edit paths and parameters here
# ============================================================

# Input CSVs — pre-aggregated domain-level LRR summaries
COASTSAT_CSV_1984_2004 = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\1984_2004\domain_lrr_summary.csv"
COASTSAT_CSV_2004_2024 = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\2004_2024\domain_lrr_summary.csv"

CS_DOMAIN_COL = "domain_number"
CS_LRR_COL    = "mean_lrr"
CS_STD_COL    = "std_lrr"

# Domain range
DOMAIN_MIN = 1
DOMAIN_MAX = 90

# LOESS smoothing window in CASCADE domains (500 m each)
#   10 domains → 5.0 km window → frac ≈ 0.111
POSTER_WINDOW_DOMAINS = 10

# Output
OUTPUT_DIR  = r"/scripts/figure_making/shoreline_change/poster_plot"
OUTPUT_FILE = "poster_two_periods_10dom.png"

# --- Line colors ---
C_CS_1984 = "#1F4E79"   # dark blue       — 1984–2004
C_CS_2004 = "#833C00"   # dark brown-red  — 2004–2024

# --- Annotation colors ---
C_TOWN_SPAN    = "#90AFC5"   # steel blue  — community spans
C_WIMBLE       = "#E0A800"   # amber       — Wimble Shoals fill
C_VILLAGE_LINE = "0.40"      # dark gray   — Salvo / Waves / Rodanthe lines
C_PIER         = "#1565C0"   # medium blue — pier lines
C_GROIN        = "#B71C1C"   # dark red    — groin lines

# --- Geographic annotations ---
TOWN_SPANS = {
    "Buxton":      (7,  8),
    "Avon":        (21, 31),
    "Tri-Village": (68, 83),
}
VILLAGE_LINES = {
    "Salvo":    69,
    "Waves":    74,
    "Rodanthe": 80,
}
PIERS = {
    "Avon Pier":     26,
    "Rodanthe Pier": 79,
}
GROINS = {
    "Buxton Groin": 6,
}
WIMBLE_SHOALS = (60, 74)

# ============================================================
# STYLE
# ============================================================
plt.rcParams.update({
    "font.family":        "Arial",
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.grid":          True,
    "grid.alpha":         0.25,
    "grid.linestyle":     ":",
})

# ============================================================
# DATA LOADING
# ============================================================

def load_coastsat(path, period_label):
    """Load a domain-level CoastSat LRR summary CSV."""
    if path is None or not os.path.exists(path):
        print(f"  CoastSat ({period_label}): SKIPPED — not found: {path}")
        return None
    df = pd.read_csv(path)
    # Fix corrupted header: filename sometimes gets prepended to the first
    # column name on export (e.g. "domain_lrr_summary.csvdomain_number").
    df.columns = [c.split("csv")[-1] if "csv" in c else c for c in df.columns]
    df[CS_DOMAIN_COL] = pd.to_numeric(df[CS_DOMAIN_COL], errors="coerce")
    df[CS_LRR_COL]    = pd.to_numeric(df[CS_LRR_COL],    errors="coerce")
    df[CS_STD_COL]    = pd.to_numeric(df[CS_STD_COL],    errors="coerce")
    df = df.dropna(subset=[CS_DOMAIN_COL, CS_LRR_COL])
    df[CS_DOMAIN_COL] = df[CS_DOMAIN_COL].astype(int)
    df = df[(df[CS_DOMAIN_COL] >= DOMAIN_MIN) & (df[CS_DOMAIN_COL] <= DOMAIN_MAX)]
    df = df[[CS_DOMAIN_COL, CS_LRR_COL, CS_STD_COL]].rename(columns={
        CS_DOMAIN_COL: "domain",
        CS_LRR_COL:    "cs_lrr",
        CS_STD_COL:    "cs_std",
    }).sort_values("domain").reset_index(drop=True)
    print(f"  CoastSat ({period_label}): {len(df)} domains  "
          f"range {df['cs_lrr'].min():+.2f} to {df['cs_lrr'].max():+.2f} m/yr")
    return df


def apply_loess(domains, values, frac):
    """Apply LOESS smoothing. Returns smoothed values at the same domain positions."""
    valid = ~np.isnan(values)
    if valid.sum() < 5:
        return values.copy()
    smoothed = np.full_like(values, np.nan)
    result = lowess(values[valid], domains[valid], frac=frac, return_sorted=True)
    smoothed[valid] = np.interp(domains[valid], result[:, 0], result[:, 1])
    return smoothed

# ============================================================
# ANNOTATION HELPERS
# ============================================================

def add_annotations(ax, village_label_y=0.43):
    """
    Add all geographic reference annotations to an axis.

    Layer order (bottom → top):
      1. Wimble Shoals influence zone  (hatched amber fill, bottom label)
      2. Community shaded spans        (steel-blue fill, top labels)
      3. Village center lines          (dashed gray, village_label_y)
      4. Pier lines                    (dash-dot blue, y=0.74, rotated)
      5. Groin lines                   (dotted red,    y=0.74, rotated)

    Parameters
    ----------
    village_label_y : float
        Axes-fraction y for Salvo/Waves/Rodanthe labels.  Default 0.25 keeps
        them below the data curves in this figure's wide y-range (~−6 to +4).
        Increase toward 0.79 for figures with a narrower y-range.
    """
    trans = blended_transform_factory(ax.transData, ax.transAxes)

    # 1. Wimble Shoals
    wlo, whi = WIMBLE_SHOALS
    ax.axvspan(wlo - 0.5, whi + 0.5,
               color=C_WIMBLE, alpha=0.10, zorder=0,
               hatch="///", edgecolor=C_WIMBLE, linewidth=0)
    ax.text((wlo + whi) / 2.0, 0.04,
            "Wimble Shoals\ninfluence", transform=trans,
            ha="center", va="bottom", fontsize=7, color="#7A5800",
            style="italic",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.80))

    # 2. Community / town spans
    for span_label, (d_lo, d_hi) in TOWN_SPANS.items():
        ax.axvspan(d_lo - 0.5, d_hi + 0.5,
                   color=C_TOWN_SPAN, alpha=0.14, zorder=0)
        ax.text((d_lo + d_hi) / 2.0, 0.90,
                span_label, transform=trans,
                ha="center", va="top", fontsize=8, color="0.25",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

    # 3. Village center lines
    # x_nudge shifts only the label (not the line) for villages near a pier.
    village_label_nudge = {
        "Rodanthe": +0,
    }
    for vname, dom in VILLAGE_LINES.items():
        ax.axvline(dom, color=C_VILLAGE_LINE, lw=0.9, ls="--",
                   alpha=0.80, zorder=1)
        label_x = dom + village_label_nudge.get(vname, 0.0)
        ax.text(label_x, village_label_y, vname, transform=trans,
                ha="center", va="top", fontsize=7.5, color="0.30",
                bbox=dict(boxstyle="round,pad=0.15", fc="white",
                          ec="none", alpha=0.80))

    # 4. Pier lines
    for pname, dom in PIERS.items():
        ax.axvline(dom, color=C_PIER, lw=1.0, ls="-.", alpha=0.80, zorder=2)
        ax.text(dom, 0.80, pname, transform=trans,
                ha="center", va="top", fontsize=7, color=C_PIER,
                rotation=90,
                bbox=dict(boxstyle="round,pad=0.15", fc="white",
                          ec="none", alpha=0.80))

    # 5. Groin lines
    for gname, dom in GROINS.items():
        ax.axvline(dom, color=C_GROIN, lw=1.1, ls=":", alpha=0.85, zorder=2)
        ax.text(dom, 0.55, gname, transform=trans,
                ha="center", va="top", fontsize=7, color=C_GROIN,
                rotation=90,
                bbox=dict(boxstyle="round,pad=0.15", fc="white",
                          ec="none", alpha=0.80))


def annotation_legend_handles():
    """Proxy artists for the annotation legend entries."""
    return [
        Patch(fc=C_TOWN_SPAN, alpha=0.30, label="Community"),
        Patch(fc=C_WIMBLE, alpha=0.25, hatch="///",
              edgecolor=C_WIMBLE, linewidth=0, label="Wimble Shoals influence"),
        Line2D([0], [0], color=C_VILLAGE_LINE, lw=0.9, ls="--",
               label="Village center"),
        Line2D([0], [0], color=C_PIER,  lw=1.0, ls="-.", label="Pier"),
        Line2D([0], [0], color=C_GROIN, lw=1.1, ls=":",  label="Groin"),
    ]

# ============================================================
# PLOT
# ============================================================

def plot_poster(cs_1984, cs_2004, out_path, window_domains=POSTER_WINDOW_DOMAINS):
    """
    Single-panel poster figure: both CoastSat periods smoothed with a fixed
    LOESS window, fill-toward-zero shading, and full geographic annotations.

    Parameters
    ----------
    cs_1984, cs_2004 : pd.DataFrame or None
    out_path : str
    window_domains : int
        LOESS window in CASCADE domains.  Default: POSTER_WINDOW_DOMAINS (10).
    """
    frac = window_domains / (DOMAIN_MAX - DOMAIN_MIN + 1)
    km   = window_domains * 0.5

    fig, ax = plt.subplots(figsize=(13, 6.0))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Annotations drawn first so data lines render on top
    add_annotations(ax)

    configs = [
        (cs_1984, "1984–2004", C_CS_1984, "calibration"),
        (cs_2004, "2004–2024", C_CS_2004, "test"),
    ]

    all_smoothed = []
    for df, label, color, role in configs:
        if df is None:
            continue
        d        = df["domain"].values.astype(float)
        smoothed = apply_loess(d, df["cs_lrr"].values, frac=frac)
        all_smoothed.append(smoothed)

        # Fill toward zero
        ax.fill_between(d, smoothed, 0,
                        where=(smoothed < 0),  color=color, alpha=0.12,
                        interpolate=True)
        ax.fill_between(d, smoothed, 0,
                        where=(smoothed >= 0), color=color, alpha=0.10,
                        interpolate=True)

        # Main smoothed line
        ax.plot(d, smoothed, color=color, lw=2.8, zorder=5,
                label=f"Period {'1' if '1984' in label else '2'}: {label}  ({role})")

    # Zero reference line
    ax.axhline(0, color="#2c2c2c", lw=1.2, ls="--", alpha=0.65, zorder=3)

    # Axes formatting
    ax.set_xlim(DOMAIN_MIN - 0.5, DOMAIN_MAX + 0.5)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis="both", which="major", labelsize=11, direction="in", length=5)
    ax.tick_params(axis="both", which="minor", direction="in", length=3)
    ax.grid(True, which="major", linestyle=":", linewidth=0.6, alpha=0.4, color="gray")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(1.2)

    # Lock ylim to smoothed data range
    if all_smoothed:
        combined = np.concatenate(all_smoothed)
        pad = (combined.max() - combined.min()) * 0.06
        ax.set_ylim(combined.min() - pad, combined.max() + pad)

    # Accretion / Erosion side labels
    # Adjust ACCRETION_Y and EROSION_Y (0.0 = panel bottom, 1.0 = panel top)
    ACCRETION_Y = 0.85
    EROSION_Y = 0.15

    ax.text(1.0, ACCRETION_Y, "Accretion ▲",
            transform=ax.transAxes, fontsize=9.5, color="#555555",
            ha="right", va="center", style="italic")
    ax.text(1.0, EROSION_Y, "Erosion ▼",
            transform=ax.transAxes, fontsize=9.5, color="#555555",
            ha="right", va="center", style="italic")

    # Axis labels
    ax.set_xlabel("CASCADE Model Domain (500 m alongshore)",
                  fontsize=12, fontweight="bold", labelpad=8)
    ax.set_ylabel("Shoreline Change Rate (m/yr)",
                  fontsize=12, fontweight="bold", labelpad=8)

    # N/S orientation
    ax.text(0.0, 1.01, "← S  |  Cape Hatteras",
            transform=ax.transAxes, fontsize=9, color="#444444",
            ha="left", va="bottom", style="italic", clip_on=False)
    ax.text(1.0, 1.01, "Pea Island  |  N →",
            transform=ax.transAxes, fontsize=9, color="#444444",
            ha="right", va="bottom", style="italic", clip_on=False)

    # Title
    ax.set_title(
        "CoastSat Shoreline Change Rates — Hatteras Island, North Carolina, USA\n"
        f"LOESS smoothed  ({window_domains}-domain window,  {km:.1f} km,  frac={frac:.3f})",
        fontsize=13, fontweight="bold", pad=12, color="#1a2a3a",
    )

    # Legend
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles=handles + annotation_legend_handles(),
              loc="lower center", bbox_to_anchor=(0.5, 0.02),
              fontsize=9.5, framealpha=0.95, edgecolor="#cccccc",
              frameon=True, ncol=4)

    # Caption
    fig.text(
        0.012, 0.005,
        f"Rates calculated as Linear Regression Rate (LRR) averaged across all CoastSat transects "
        f"within each 500-m CASCADE model domain, smoothed with LOESS "
        f"({window_domains}-domain / {km:.1f} km window).  "
        f"Shoreline data: CoastSat satellite-derived shorelines (Vos et al., 2019).",
        fontsize=8, color="#666666", ha="left", va="bottom", style="italic",
    )

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")

# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 55)
    print("CoastSat Poster Figure — Hatteras Island")
    print(f"LOESS window: {POSTER_WINDOW_DOMAINS} domains  "
          f"({POSTER_WINDOW_DOMAINS * 0.5:.1f} km)")
    print("=" * 55)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("\nLoading data...")
    cs_1984 = load_coastsat(COASTSAT_CSV_1984_2004, "1984–2004")
    cs_2004 = load_coastsat(COASTSAT_CSV_2004_2024, "2004–2024")

    print("\nGenerating poster figure...")
    plot_poster(cs_1984, cs_2004,
                os.path.join(OUTPUT_DIR, OUTPUT_FILE))

    print("\n" + "=" * 55)
    print(f"Done!  →  {OUTPUT_FILE}")
    print("=" * 55)


if __name__ == "__main__":
    main()
