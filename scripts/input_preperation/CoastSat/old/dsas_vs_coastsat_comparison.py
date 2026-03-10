"""
DSAS vs CoastSat LRR Comparison — Hatteras Island
===================================================
Compares domain-level Linear Regression Rates (LRR) from:

  DSAS     : LRR aggregated from DSAS transect-level output across
             digitized shorelines (1978, 1997, 2019).
             Column used: MEAN_LRR (mean across transects per domain)
             Periods: 1978–1997 and 1997–2019

  CoastSat : LRR from satellite imagery (100s of observations per transect)
             Column used: mean_lrr from domain_lrr_summary.csv
             Periods: 1978–1997 and 1997–2019

This gives a direct apples-to-apples comparison:
  same periods, same domain aggregation, different data sources.

Key methodological difference:
  DSAS LRR uses ~5 transects per domain from 3 digitized shorelines.
  CoastSat LRR uses ~10 transects per domain from 100s of satellite images.
  Differences reflect both data source uncertainty and shoreline definition
  (DSAS = MHW digitized, CoastSat = wet/dry line from Landsat).

Outputs (saved to OUTPUT_DIR)
-------
  line_comparison_1978_1997.png   – along-island line plot, period 1
  line_comparison_1997_2019.png   – along-island line plot, period 2
  scatter_1978_1997.png           – per-domain scatter + regression, period 1
  scatter_1997_2019.png           – per-domain scatter + regression, period 2
  difference_1978_1997.png        – CoastSat minus DSAS bar chart, period 1
  difference_1997_2019.png        – CoastSat minus DSAS bar chart, period 2
  both_periods_overview.png       – 2-panel overview of both periods
  comparison_table.csv            – all values and statistics per domain
"""

# ============================================================
# CONFIG  –  edit paths before running
# ============================================================

# --- DSAS CSVs (uploaded files) ---
DSAS_CSV_1978_1997 = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\shoreline_change\dsas_1978_1997_rates.csv"
DSAS_CSV_1997_2019 = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\shoreline_change\dsas_1997_2019_rates.csv"

DSAS_DOMAIN_COL = "domain_id"   # domain number column
DSAS_LRR_COL    = "MEAN_LRR"    # mean LRR across transects per domain
DSAS_STD_COL    = "STD_LRR"     # std across transects per domain

# --- CoastSat CSVs (from coastsat_domain_lrr.py) ---
COASTSAT_CSV_1978_1997 = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\1978_1997\transect_lrr_full.csv"
COASTSAT_CSV_1997_2019 = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\1997_2019\transect_lrr_full.csv"

CS_DOMAIN_COL = "domain_number"
CS_LRR_COL    = "mean_lrr"      # "mean_lrr" or "median_lrr"
CS_STD_COL    = "std_lrr"

# --- Domain range: your 90 real Hatteras domains ---
DOMAIN_MIN = 1
DOMAIN_MAX = 90

# --- Rodanthe highlight ---
RODANTHE_MIN = 77
RODANTHE_MAX = 83

# --- Output ---
OUTPUT_DIR = r"C:/path/to/comparison_outputs"

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

def load_dsas(path, domain_col, lrr_col, std_col, period_label):
    df = pd.read_csv(path)
    df[domain_col] = pd.to_numeric(df[domain_col], errors="coerce")
    df[lrr_col]    = pd.to_numeric(df[lrr_col],    errors="coerce")
    df[std_col]    = pd.to_numeric(df[std_col],     errors="coerce")
    df = df.dropna(subset=[domain_col, lrr_col])
    df[domain_col] = df[domain_col].astype(int)
    df = df[(df[domain_col] >= DOMAIN_MIN) & (df[domain_col] <= DOMAIN_MAX)]
    df = df[[domain_col, lrr_col, std_col]].rename(columns={
        domain_col: "domain",
        lrr_col   : "dsas_lrr",
        std_col   : "dsas_std",
    }).sort_values("domain").reset_index(drop=True)
    print(f"DSAS ({period_label}): {len(df)} domains with valid LRR  "
          f"| range: {df['dsas_lrr'].min():+.2f} to {df['dsas_lrr'].max():+.2f} m/yr")
    return df


def load_coastsat(path, domain_col, lrr_col, std_col, period_label):
    if path is None or not os.path.exists(path):
        print(f"CoastSat ({period_label}): SKIPPED — file not found: {path}")
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
        lrr_col   : "cs_lrr",
        std_col   : "cs_std",
    }).sort_values("domain").reset_index(drop=True)
    print(f"CoastSat ({period_label}): {len(df)} domains with valid LRR  "
          f"| range: {df['cs_lrr'].min():+.2f} to {df['cs_lrr'].max():+.2f} m/yr")
    return df


def merge_datasets(dsas, coastsat):
    """Inner join on domain — only domains present in BOTH datasets."""
    merged = dsas.merge(coastsat, on="domain", how="inner")
    merged["difference"] = merged["cs_lrr"] - merged["dsas_lrr"]
    return merged.sort_values("domain").reset_index(drop=True)


def regression_stats(merged):
    x = merged["dsas_lrr"].values
    y = merged["cs_lrr"].values
    slope, intercept, r, p, _ = stats.linregress(x, y)
    rmse = np.sqrt(np.mean((y - x) ** 2))
    bias = np.mean(y - x)
    return dict(slope=slope, intercept=intercept, r2=r**2,
                rmse=rmse, bias=bias, p=p, n=len(x))


def rodanthe_stats(merged):
    rod = merged[(merged["domain"] >= RODANTHE_MIN) &
                 (merged["domain"] <= RODANTHE_MAX)]
    return rod


# ============================================================
# PLOT FUNCTIONS
# ============================================================

def plot_line_comparison(merged, period_label, out_path):
    """Along-island line plot: DSAS vs CoastSat with ±1 std shading."""
    fig, ax = plt.subplots(figsize=(16, 6))

    # DSAS
    ax.plot(merged["domain"], merged["dsas_lrr"],
            color="#2166ac", lw=2.5, marker="o", ms=4,
            label=f"DSAS LRR {period_label}")
    ax.fill_between(merged["domain"],
                    merged["dsas_lrr"] - merged["dsas_std"],
                    merged["dsas_lrr"] + merged["dsas_std"],
                    color="#2166ac", alpha=0.12)

    # CoastSat
    ax.plot(merged["domain"], merged["cs_lrr"],
            color="#b2182b", lw=2.5, marker="s", ms=4,
            label=f"CoastSat LRR {period_label}")
    ax.fill_between(merged["domain"],
                    merged["cs_lrr"] - merged["cs_std"],
                    merged["cs_lrr"] + merged["cs_std"],
                    color="#b2182b", alpha=0.12)

    ax.axvspan(RODANTHE_MIN, RODANTHE_MAX, alpha=0.13,
               color="orange", label="Rodanthe", zorder=0)
    ax.axhline(0, color="black", lw=1.2, ls="--", alpha=0.6)

    ax.set_xlabel("Domain Number", fontsize=13, fontweight="bold")
    ax.set_ylabel("LRR (m/yr)\n← Erosion | Accretion →", fontsize=12, fontweight="bold")
    ax.set_title(f"DSAS vs CoastSat LRR — {period_label}  (shading = ±1 std)",
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


def plot_scatter(merged, period_label, out_path):
    """
    Per-domain scatter: DSAS (x) vs CoastSat (y).
    1:1 line + regression + stats annotation.
    Rodanthe domains highlighted.
    """
    s = regression_stats(merged)
    rod_mask = ((merged["domain"] >= RODANTHE_MIN) &
                (merged["domain"] <= RODANTHE_MAX))

    fig, ax = plt.subplots(figsize=(8, 8))

    ax.scatter(merged.loc[~rod_mask, "dsas_lrr"],
               merged.loc[~rod_mask, "cs_lrr"],
               color="steelblue", s=45, alpha=0.75, zorder=3,
               label="Other domains")
    ax.scatter(merged.loc[rod_mask, "dsas_lrr"],
               merged.loc[rod_mask, "cs_lrr"],
               color="darkorange", s=65, marker="D", zorder=4,
               label=f"Rodanthe (D{RODANTHE_MIN}–{RODANTHE_MAX})")

    # Label Rodanthe points
    for _, row in merged[rod_mask].iterrows():
        ax.annotate(str(int(row["domain"])),
                    xy=(row["dsas_lrr"], row["cs_lrr"]),
                    xytext=(5, 4), textcoords="offset points",
                    fontsize=8, color="darkorange")

    # Axis limits
    all_vals = pd.concat([merged["dsas_lrr"], merged["cs_lrr"]])
    lim = [all_vals.min() - 0.5, all_vals.max() + 0.5]

    # 1:1 line
    ax.plot(lim, lim, "k--", lw=1.5, alpha=0.55, label="1:1 line", zorder=1)

    # Best-fit line
    x_fit = np.linspace(lim[0], lim[1], 100)
    y_fit = s["slope"] * x_fit + s["intercept"]
    ax.plot(x_fit, y_fit, color="#b2182b", lw=2.2, alpha=0.85,
            label=f"Best fit  R²={s['r2']:.2f}", zorder=2)

    # Stats box
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
    ax.set_xlabel("DSAS LRR (m/yr)", fontsize=13, fontweight="bold")
    ax.set_ylabel("CoastSat LRR (m/yr)", fontsize=13, fontweight="bold")
    ax.set_title(f"DSAS vs CoastSat — {period_label}  (per domain)",
                 fontsize=13, fontweight="bold", pad=12)
    ax.legend(fontsize=10, framealpha=0.95)
    ax.grid(True, alpha=0.25, ls=":")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


def plot_difference(merged, period_label, out_path):
    """
    Bar chart of CoastSat − DSAS by domain.
    Blue = CoastSat more accretional, Red = CoastSat more erosional.
    """
    diff   = merged["difference"]
    colors = ["#2166ac" if v >= 0 else "#b2182b" for v in diff]
    mean_d = diff.mean()

    fig, ax = plt.subplots(figsize=(16, 5))
    ax.bar(merged["domain"], diff, color=colors, edgecolor="none", width=0.85)
    ax.axhline(0, color="black", lw=1.2)
    ax.axhline(mean_d, color="grey", lw=1.8, ls="--",
               label=f"Mean difference ({mean_d:+.2f} m/yr)")
    ax.axvspan(RODANTHE_MIN, RODANTHE_MAX, alpha=0.13,
               color="orange", label="Rodanthe", zorder=0)

    ax.set_xlabel("Domain Number", fontsize=13, fontweight="bold")
    ax.set_ylabel("CoastSat − DSAS (m/yr)", fontsize=12, fontweight="bold")
    ax.set_title(f"Dataset Difference: CoastSat minus DSAS — {period_label}\n"
                 f"Blue = CoastSat more accretional  |  Red = CoastSat more erosional",
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


def plot_both_periods_overview(merged_1978, merged_1997, out_path):
    """
    Two-panel overview: one panel per calibration period.
    Both DSAS and CoastSat shown on each panel.
    """
    fig, axes = plt.subplots(2, 1, figsize=(16, 11), sharex=True)

    for ax, merged, period_label, dsas_color, cs_color in [
        (axes[0], merged_1978, "1978–1997", "#2166ac", "#b2182b"),
        (axes[1], merged_1997, "1997–2019", "#1a9641", "#d73027"),
    ]:
        if merged is None:
            ax.text(0.5, 0.5, f"No data for {period_label}",
                    transform=ax.transAxes, ha="center", fontsize=13)
            continue

        s = regression_stats(merged)

        ax.plot(merged["domain"], merged["dsas_lrr"],
                color=dsas_color, lw=2.5, marker="o", ms=3.5,
                label=f"DSAS {period_label}")
        ax.fill_between(merged["domain"],
                        merged["dsas_lrr"] - merged["dsas_std"],
                        merged["dsas_lrr"] + merged["dsas_std"],
                        color=dsas_color, alpha=0.10)

        ax.plot(merged["domain"], merged["cs_lrr"],
                color=cs_color, lw=2.5, marker="s", ms=3.5,
                label=f"CoastSat {period_label}")
        ax.fill_between(merged["domain"],
                        merged["cs_lrr"] - merged["cs_std"],
                        merged["cs_lrr"] + merged["cs_std"],
                        color=cs_color, alpha=0.10)

        ax.axvspan(RODANTHE_MIN, RODANTHE_MAX, alpha=0.13,
                   color="orange", label="Rodanthe", zorder=0)
        ax.axhline(0, color="black", lw=1.1, ls="--", alpha=0.6)

        # Stats annotation
        txt = f"R²={s['r2']:.2f}  RMSE={s['rmse']:.2f} m/yr  Bias={s['bias']:+.2f} m/yr"
        ax.text(0.02, 0.97, txt, transform=ax.transAxes, fontsize=9,
                va="top", bbox=dict(boxstyle="round", facecolor="white",
                                    alpha=0.85, edgecolor="grey"))

        ax.set_ylabel("LRR (m/yr)\n← Erosion | Accretion →",
                      fontsize=11, fontweight="bold")
        ax.legend(fontsize=10, framealpha=0.95, loc="upper right")
        ax.grid(True, alpha=0.25, ls=":")

    axes[-1].set_xlabel("Domain Number", fontsize=13, fontweight="bold")
    axes[-1].set_xlim(DOMAIN_MIN - 1, DOMAIN_MAX + 1)
    axes[-1].set_xticks(range(DOMAIN_MIN, DOMAIN_MAX + 1, 10))
    axes[-1].set_xticks(range(DOMAIN_MIN, DOMAIN_MAX + 1, 5), minor=True)
    fig.suptitle("DSAS vs CoastSat LRR — Both Calibration Periods",
                 fontsize=15, fontweight="bold", y=1.01)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 65)
    print("DSAS vs CoastSat LRR Comparison — Hatteras Island")
    print("=" * 65)

    # ---- Load ----
    print("\nLoading data...")
    dsas_1978 = load_dsas(DSAS_CSV_1978_1997, DSAS_DOMAIN_COL,
                          DSAS_LRR_COL, DSAS_STD_COL, "1978–1997")
    dsas_1997 = load_dsas(DSAS_CSV_1997_2019, DSAS_DOMAIN_COL,
                          DSAS_LRR_COL, DSAS_STD_COL, "1997–2019")
    cs_1978   = load_coastsat(COASTSAT_CSV_1978_1997, CS_DOMAIN_COL,
                              CS_LRR_COL, CS_STD_COL, "1978–1997")
    cs_1997   = load_coastsat(COASTSAT_CSV_1997_2019, CS_DOMAIN_COL,
                              CS_LRR_COL, CS_STD_COL, "1997–2019")

    # ---- Merge ----
    merged_1978 = merge_datasets(dsas_1978, cs_1978) if cs_1978 is not None else None
    merged_1997 = merge_datasets(dsas_1997, cs_1997) if cs_1997 is not None else None

    # ---- Stats summary ----
    for merged, label in [(merged_1978, "1978–1997"), (merged_1997, "1997–2019")]:
        if merged is None:
            continue
        s   = regression_stats(merged)
        rod = rodanthe_stats(merged)
        print(f"\n--- {label} ---")
        print(f"  Domains compared   : {s['n']}")
        print(f"  R²                 : {s['r2']:.3f}")
        print(f"  RMSE               : {s['rmse']:.3f} m/yr")
        print(f"  Bias (CS − DSAS)   : {s['bias']:+.3f} m/yr")
        print(f"  Best-fit slope     : {s['slope']:.3f}  (1.0 = perfect)")
        print(f"  Rodanthe DSAS mean : {rod['dsas_lrr'].mean():+.3f} m/yr")
        print(f"  Rodanthe CS mean   : {rod['cs_lrr'].mean():+.3f} m/yr")
        print(f"  Rodanthe diff mean : {rod['difference'].mean():+.3f} m/yr")

    # ---- Plots ----
    print("\nGenerating plots...")
    if merged_1978 is not None:
        plot_line_comparison(merged_1978, "1978–1997",
            os.path.join(OUTPUT_DIR, "line_comparison_1978_1997.png"))
        plot_scatter(merged_1978, "1978–1997",
            os.path.join(OUTPUT_DIR, "scatter_1978_1997.png"))
        plot_difference(merged_1978, "1978–1997",
            os.path.join(OUTPUT_DIR, "difference_1978_1997.png"))

    if merged_1997 is not None:
        plot_line_comparison(merged_1997, "1997–2019",
            os.path.join(OUTPUT_DIR, "line_comparison_1997_2019.png"))
        plot_scatter(merged_1997, "1997–2019",
            os.path.join(OUTPUT_DIR, "scatter_1997_2019.png"))
        plot_difference(merged_1997, "1997–2019",
            os.path.join(OUTPUT_DIR, "difference_1997_2019.png"))

    plot_both_periods_overview(merged_1978, merged_1997,
        os.path.join(OUTPUT_DIR, "both_periods_overview.png"))

    # ---- Save comparison table ----
    table_parts = []
    for merged, label in [(merged_1978, "1978_1997"), (merged_1997, "1997_2019")]:
        if merged is None:
            continue
        t = merged.copy()
        t.columns = ["domain"] + [f"{c}_{label}" for c in t.columns if c != "domain"]
        table_parts.append(t)

    if table_parts:
        from functools import reduce
        table = reduce(lambda a, b: a.merge(b, on="domain", how="outer"), table_parts)
        table = table.sort_values("domain")
        out_csv = os.path.join(OUTPUT_DIR, "comparison_table.csv")
        table.to_csv(out_csv, index=False)
        print(f"  Saved: comparison_table.csv")

    print("\n" + "=" * 65)
    print("Done!")
    print("=" * 65)


if __name__ == "__main__":
    main()
