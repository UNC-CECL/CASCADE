"""
CoastSat LRR Comparison — Two Periods
======================================
Compares domain-level Linear Regression Rates (LRR) between two
CoastSat periods:
    Period 1 : 1984–2004
    Period 2 : 2004–2024

Both loaded from domain_lrr_summary.csv files produced by
coastsat_domain_lrr.py.

Outputs (saved to OUTPUT_DIR)
-------
  line_comparison_two_periods.png   – along-island line plot
  scatter_two_periods.png           – per-domain scatter + regression
  difference_two_periods.png        – Period 2 minus Period 1 bar chart
  comparison_table.csv              – all values and stats per domain
"""

# ============================================================
# CONFIG  –  edit paths before running
# ============================================================

COASTSAT_CSV_P1 = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\1984_2004\domain_lrr_summary.csv"
COASTSAT_CSV_P2 = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\2004_2024\domain_lrr_summary.csv"

PERIOD_1_LABEL = "1984–2004"
PERIOD_2_LABEL = "2004–2024"

CS_DOMAIN_COL = "domain_number"
CS_LRR_COL    = "mean_lrr"
CS_STD_COL    = "std_lrr"

DOMAIN_MIN = 1
DOMAIN_MAX = 90

OUTPUT_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\two_period_comparison"

# ============================================================
# IMPORTS
# ============================================================
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ============================================================
# HELPERS
# ============================================================

def load_coastsat(path, domain_col, lrr_col, std_col, period_label):
    if path is None or not os.path.exists(path):
        raise FileNotFoundError(
            f"CoastSat ({period_label}): file not found: {path}")
    df = pd.read_csv(path)
    df[domain_col] = pd.to_numeric(df[domain_col], errors="coerce")
    df[lrr_col]    = pd.to_numeric(df[lrr_col],    errors="coerce")
    df[std_col]    = pd.to_numeric(df[std_col],     errors="coerce")
    df = df.dropna(subset=[domain_col, lrr_col])
    df[domain_col] = df[domain_col].astype(int)
    df = df[(df[domain_col] >= DOMAIN_MIN) & (df[domain_col] <= DOMAIN_MAX)]
    df = df[[domain_col, lrr_col, std_col]].rename(columns={
        domain_col: "domain",
        lrr_col   : "lrr",
        std_col   : "std",
    }).sort_values("domain").reset_index(drop=True)
    print(f"CoastSat ({period_label}): {len(df)} domains  "
          f"| LRR range: {df['lrr'].min():+.2f} to {df['lrr'].max():+.2f} m/yr")
    return df


def merge_periods(p1, p2):
    """Inner join on domain — only domains present in BOTH periods."""
    merged = p1.merge(p2, on="domain", suffixes=("_p1", "_p2"), how="inner")
    merged["difference"] = merged["lrr_p2"] - merged["lrr_p1"]
    return merged.sort_values("domain").reset_index(drop=True)


def regression_stats(merged):
    x = merged["lrr_p1"].values
    y = merged["lrr_p2"].values
    slope, intercept, r, p, _ = stats.linregress(x, y)
    rmse = np.sqrt(np.mean((y - x) ** 2))
    bias = np.mean(y - x)
    return dict(slope=slope, intercept=intercept, r2=r**2,
                rmse=rmse, bias=bias, p=p, n=len(x))


# ============================================================
# PLOT FUNCTIONS
# ============================================================

def plot_line_comparison(merged, out_path):
    """Along-island line plot: Period 1 vs Period 2 with ±1 std shading."""
    fig, ax = plt.subplots(figsize=(16, 6))

    ax.plot(merged["domain"], merged["lrr_p1"],
            color="#2166ac", lw=2.5, marker="o", ms=4,
            label=f"CoastSat {PERIOD_1_LABEL}")
    ax.fill_between(merged["domain"],
                    merged["lrr_p1"] - merged["std_p1"],
                    merged["lrr_p1"] + merged["std_p1"],
                    color="#2166ac", alpha=0.12)

    ax.plot(merged["domain"], merged["lrr_p2"],
            color="#b2182b", lw=2.5, marker="s", ms=4,
            label=f"CoastSat {PERIOD_2_LABEL}")
    ax.fill_between(merged["domain"],
                    merged["lrr_p2"] - merged["std_p2"],
                    merged["lrr_p2"] + merged["std_p2"],
                    color="#b2182b", alpha=0.12)

    ax.axhline(0, color="black", lw=1.2, ls="--", alpha=0.6)
    ax.set_xlabel("Domain Number", fontsize=13, fontweight="bold")
    ax.set_ylabel("LRR (m/yr)\n← Erosion | Accretion →", fontsize=12, fontweight="bold")
    ax.set_title(f"CoastSat LRR: {PERIOD_1_LABEL} vs {PERIOD_2_LABEL}  "
                 f"(shading = ±1 std)",
                 fontsize=14, fontweight="bold", pad=14)
    ax.legend(fontsize=10, framealpha=0.95)
    ax.grid(True, alpha=0.25, ls=":")
    ax.set_xlim(DOMAIN_MIN - 1, DOMAIN_MAX + 1)
    ax.set_xticks(range(DOMAIN_MIN, DOMAIN_MAX + 1, 10))
    ax.set_xticks(range(DOMAIN_MIN, DOMAIN_MAX + 1, 5), minor=True)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


def plot_scatter(merged, out_path):
    """Per-domain scatter: Period 1 (x) vs Period 2 (y), 1:1 + regression."""
    s = regression_stats(merged)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(merged["lrr_p1"], merged["lrr_p2"],
               color="steelblue", s=45, alpha=0.75, zorder=3,
               label="Domains")

    all_vals = pd.concat([merged["lrr_p1"], merged["lrr_p2"]])
    lim = [all_vals.min() - 0.5, all_vals.max() + 0.5]

    ax.plot(lim, lim, "k--", lw=1.5, alpha=0.55, label="1:1 line", zorder=1)

    x_fit = np.linspace(lim[0], lim[1], 100)
    y_fit = s["slope"] * x_fit + s["intercept"]
    ax.plot(x_fit, y_fit, color="#b2182b", lw=2.2, alpha=0.85,
            label=f"Best fit  R²={s['r2']:.2f}", zorder=2)

    txt = (f"R²   = {s['r2']:.3f}\n"
           f"RMSE = {s['rmse']:.3f} m/yr\n"
           f"Bias = {s['bias']:+.3f} m/yr\n"
           f"n    = {s['n']} domains")
    ax.text(0.04, 0.97, txt, transform=ax.transAxes, fontsize=10,
            va="top", family="monospace",
            bbox=dict(boxstyle="round", facecolor="white",
                      alpha=0.88, edgecolor="grey"))

    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel(f"CoastSat LRR {PERIOD_1_LABEL} (m/yr)", fontsize=13, fontweight="bold")
    ax.set_ylabel(f"CoastSat LRR {PERIOD_2_LABEL} (m/yr)", fontsize=13, fontweight="bold")
    ax.set_title(f"CoastSat LRR: {PERIOD_1_LABEL} vs {PERIOD_2_LABEL}  (per domain)",
                 fontsize=13, fontweight="bold", pad=12)
    ax.legend(fontsize=10, framealpha=0.95)
    ax.grid(True, alpha=0.25, ls=":")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


def plot_difference(merged, out_path):
    """
    Bar chart of Period 2 − Period 1 by domain.
    Blue = more accretional in P2, Red = more erosional in P2.
    """
    diff   = merged["difference"]
    colors = ["#2166ac" if v >= 0 else "#b2182b" for v in diff]
    mean_d = diff.mean()

    fig, ax = plt.subplots(figsize=(16, 5))
    ax.bar(merged["domain"], diff, color=colors, edgecolor="none", width=0.85)
    ax.axhline(0, color="black", lw=1.2)
    ax.axhline(mean_d, color="grey", lw=1.8, ls="--",
               label=f"Mean difference ({mean_d:+.2f} m/yr)")

    ax.set_xlabel("Domain Number", fontsize=13, fontweight="bold")
    ax.set_ylabel(f"{PERIOD_2_LABEL} − {PERIOD_1_LABEL} (m/yr)",
                  fontsize=12, fontweight="bold")
    ax.set_title(f"LRR Change Between Periods: {PERIOD_2_LABEL} minus {PERIOD_1_LABEL}\n"
                 f"Blue = more accretional in {PERIOD_2_LABEL}  |  "
                 f"Red = more erosional in {PERIOD_2_LABEL}",
                 fontsize=13, fontweight="bold", pad=12)
    ax.legend(fontsize=10, framealpha=0.95)
    ax.grid(True, axis="y", alpha=0.25, ls=":")
    ax.set_xlim(DOMAIN_MIN - 1, DOMAIN_MAX + 1)
    ax.set_xticks(range(DOMAIN_MIN, DOMAIN_MAX + 1, 10))
    ax.set_xticks(range(DOMAIN_MIN, DOMAIN_MAX + 1, 5), minor=True)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 65)
    print("CoastSat LRR Comparison — Two Periods")
    print("=" * 65)

    # ---- Load ----
    print("\nLoading data...")
    cs_p1 = load_coastsat(COASTSAT_CSV_P1, CS_DOMAIN_COL,
                          CS_LRR_COL, CS_STD_COL, PERIOD_1_LABEL)
    cs_p2 = load_coastsat(COASTSAT_CSV_P2, CS_DOMAIN_COL,
                          CS_LRR_COL, CS_STD_COL, PERIOD_2_LABEL)

    # ---- Merge ----
    merged = merge_periods(cs_p1, cs_p2)
    print(f"\nDomains in both periods: {len(merged)}")

    # ---- Stats summary ----
    s = regression_stats(merged)
    print(f"\n--- Period-to-Period Statistics ---")
    print(f"  Domains compared        : {s['n']}")
    print(f"  R²                      : {s['r2']:.3f}")
    print(f"  RMSE                    : {s['rmse']:.3f} m/yr")
    print(f"  Bias (P2 − P1)          : {s['bias']:+.3f} m/yr")
    print(f"  Best-fit slope          : {s['slope']:.3f}  (1.0 = no change)")

    # ---- Plots ----
    print("\nGenerating plots...")
    plot_line_comparison(merged,
        os.path.join(OUTPUT_DIR, "line_comparison_two_periods.png"))
    plot_scatter(merged,
        os.path.join(OUTPUT_DIR, "scatter_two_periods.png"))
    plot_difference(merged,
        os.path.join(OUTPUT_DIR, "difference_two_periods.png"))

    # ---- Save comparison table ----
    out_csv = os.path.join(OUTPUT_DIR, "comparison_table.csv")
    merged.to_csv(out_csv, index=False)
    print(f"  Saved: comparison_table.csv")

    print("\n" + "=" * 65)
    print("Done!")
    print("=" * 65)


if __name__ == "__main__":
    main()
