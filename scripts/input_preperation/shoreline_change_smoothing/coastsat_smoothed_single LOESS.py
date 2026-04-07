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
LOESS_FRAC = 0.111

# --- Town locations for reference lines ---
TOWNS = {
    "Buxton": 8,
    "Avon": 26,
    "Salvo": 69,
    "Waves": 74,
    "Rodanthe": 80,
}

# --- Output ---
OUTPUT_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\shoreline_change_smoothing\coastsat_smoothing_10domains"

# ============================================================
# IMPORTS
# ============================================================
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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


# ============================================================
# SHARED AXIS HELPERS
# ============================================================

def add_town_lines(ax, ymax=None):
    for town, dom in TOWNS.items():
        ax.axvline(dom, color="0.55", lw=0.9, ls="--", alpha=0.7, zorder=0)
        y_pos = ax.get_ylim()[1] * 0.93 if ymax is None else ymax * 0.93
        ax.text(dom, y_pos, town, ha="center", va="top",
                fontsize=8, color="0.35",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))


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

    for ax, df, label, color in configs:
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
        ax.legend(fontsize=9.5, framealpha=0.95, loc="lower right")
        style_domain_axis(ax)
        add_town_lines(ax)

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

    for ax, df, label, color in configs:
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
        ax.legend(fontsize=9.5, framealpha=0.95, loc="lower right")
        style_domain_axis(ax)
        add_town_lines(ax)

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
    ax.legend(fontsize=10, framealpha=0.95)
    style_domain_axis(ax)
    add_town_lines(ax)

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

def plot_smoothing_sensitivity(df, period_label, out_path):
    """
    3-panel showing effect of different LOESS bandwidths on CoastSat data.
    """
    fracs = [0.10, 0.15, 0.20]
    labels_frac = [
        "frac=0.10  (~9 domains)",
        "frac=0.15  (~13 domains) ← default",
        "frac=0.20  (~18 domains)",
    ]

    color = C_CS_1984 if "1984" in period_label else C_CS_2004

    fig, axes = plt.subplots(3, 1, figsize=(16, 13), sharex=True, sharey=True)
    fig.suptitle(f"LOESS Smoothing Sensitivity — CoastSat {period_label}\n"
                 f"Effect of bandwidth choice",
                 fontsize=13, fontweight="bold", y=1.01)

    for ax, frac, flabel in zip(axes, fracs, labels_frac):
        m = add_smoothed_columns(df, frac=frac)
        d = m["domain"]

        # Raw faded
        ax.plot(d, m["cs_lrr"], color=color, lw=0.8, alpha=0.25,
                marker="o", ms=2, zorder=1, label="Raw")

        # Smoothed
        ax.plot(d, m["cs_lrr_smooth"], color=color, lw=2.5,
                label="LOESS smoothed", zorder=3)

        ax.text(0.01, 0.97, flabel,
                transform=ax.transAxes, fontsize=8.5, va="top",
                bbox=dict(boxstyle="round", fc="white", alpha=0.88, ec="0.7"))

        ax.set_ylabel("Rate (m/yr)", fontsize=10, fontweight="bold")
        ax.legend(fontsize=9, framealpha=0.95, loc="lower right")
        style_domain_axis(ax)
        add_town_lines(ax)

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
    print("=" * 65)


if __name__ == "__main__":
    main()
