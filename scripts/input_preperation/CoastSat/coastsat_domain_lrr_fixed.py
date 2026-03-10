"""
CoastSat Domain-Level LRR Summary
===================================
Step 2 of 2 in the domain-level LRR workflow.

Requires:
  1. transect_domain_lookup.csv   – from coastsat_domain_mapping.py
  2. CoastSat time-series CSVs    – one per transect
  3. coastsat_lrr_analysis.py     – in the same directory (or on PYTHONPATH)

Outputs:
  domain_lrr_summary.csv  –  one row per domain with aggregated LRR stats
  domain_lrr_plot.png     –  bar chart coloured by erosion/accretion
  transect_lrr_full.csv   –  full transect-level results with domain assignments

Usage
-----
Edit the CONFIG section below, then run:
    python coastsat_domain_lrr.py
"""

# ============================================================
# CONFIG
# ============================================================

# Path to the lookup table produced by coastsat_domain_mapping.py
LOOKUP_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\transect_domain_lookup.csv"

# Root folder containing all site subfolders (e.g. usa_NC_0032_timeseries, usa_NC_0033_timeseries, ...)
# The script will automatically find every CSV in every subfolder one level down.
# Example: r"C:/Users/hahenry/Downloads"
ROOT_DATA_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\coastsat_timeseries"

# Optional: only include subfolders whose names contain this string.
# Set to "" to include ALL subfolders under ROOT_DATA_DIR.
SITE_FILTER = "usa_NC"    # e.g. "usa_NC" to match usa_NC_0032_timeseries, usa_NC_0033_timeseries, etc.

# Date range for LRR calculation
START_DATE = "1997-01-01"
END_DATE   = "2019-12-31"

# Minimum observations per transect to include it in domain summaries
MIN_OBS = 3

# Output directory
OUTPUT_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\1997_2019"

# CASCADE buffer domains to EXCLUDE from summaries
# (e.g., the 15 buffer domains on each end of your 90+30 setup)
# Set to an empty list [] to include all domains
BUFFER_DOMAINS = []   # all 90 real domains included; CoastSat transects won't match buffer domains anyway

# ============================================================
# IMPORTS
# ============================================================
import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator
import warnings
warnings.filterwarnings("ignore")

# Import LRR functions from the companion script
# (assumes both scripts are in the same directory)
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)
from coastsat_lrr_analysis import (
    load_timeseries, filter_dates, compute_lrr, _empty_lrr
)


# ============================================================
# FUNCTIONS
# ============================================================

def collect_csv_map(root_dir: str, site_filter: str = "") -> dict:
    """
    Auto-discover all time-series CSVs under root_dir.

    Walks one level of subfolders (e.g. usa_NC_0032_timeseries/) and
    collects every CSV inside them. Optionally filters to subfolders
    whose names contain site_filter.

    Returns:
        { transect_id_stem : full_filepath }
        e.g. { 'usa_NC_0032_0011' : 'C:/Downloads/usa_NC_0032_timeseries/usa_NC_0032_0011.csv' }
    """
    csv_map = {}
    if not os.path.isdir(root_dir):
        print(f"⚠️  ROOT_DATA_DIR not found: {root_dir}")
        return csv_map

    subfolders = [
        os.path.join(root_dir, d)
        for d in sorted(os.listdir(root_dir))
        if os.path.isdir(os.path.join(root_dir, d))
        and (site_filter == "" or site_filter in d)
    ]

    print(f"Found {len(subfolders)} site folder(s) under {root_dir}:")
    for sf in subfolders:
        csvs = glob.glob(os.path.join(sf, "*.csv"))
        for fpath in csvs:
            stem = os.path.splitext(os.path.basename(fpath))[0]
            csv_map[stem] = fpath
        print(f"  {os.path.basename(sf):40s}  →  {len(csvs)} CSVs")

    return csv_map


def compute_all_lrr(lookup: pd.DataFrame, csv_map: dict,
                    start: str, end: str, min_obs: int) -> pd.DataFrame:
    """
    For every transect in the lookup table, find its CSV, compute LRR,
    and return a merged DataFrame with domain assignments.
    """
    records = []
    missing = []

    for _, row in lookup.iterrows():
        tid = str(row["transect_id"])

        if tid not in csv_map:
            missing.append(tid)
            result = _empty_lrr(0)
        else:
            try:
                df_raw  = load_timeseries(csv_map[tid])
                df_filt = filter_dates(df_raw, start, end)
                if len(df_filt) < min_obs:
                    result = _empty_lrr(len(df_filt))
                else:
                    result = compute_lrr(df_filt)
            except Exception as e:
                print(f"  ERROR {tid}: {e}")
                result = _empty_lrr(0)

        records.append({
            "transect_id"  : tid,
            "domain_number": row["domain_number"],
            "match_method" : row.get("match_method", ""),
            **result
        })

    if missing:
        print(f"\n  ⚠️  {len(missing)} transects in lookup not found in CSVs.")
        print(f"     First few: {missing[:5]}")

    return pd.DataFrame(records)


def domain_summary(transect_df: pd.DataFrame,
                   buffer_domains: list) -> pd.DataFrame:
    """
    Aggregate transect-level LRR results to domain level.

    Returns one row per domain with:
      n_transects    – total transects in domain
      n_valid        – transects with valid LRR
      mean_lrr       – mean LRR across transects (m/yr)
      median_lrr     – median LRR
      std_lrr        – standard deviation
      min_lrr        – most erosional transect
      max_lrr        – most accretionary transect
      pct_eroding    – % of transects with negative LRR
    """
    df = transect_df.copy()
    if buffer_domains:
        df = df[~df["domain_number"].isin(buffer_domains)]

    df_valid = df[df["lrr_m_yr"].notna()]

    summary = (
        df_valid.groupby("domain_number")["lrr_m_yr"]
        .agg(
            n_valid     = "count",
            mean_lrr    = "mean",
            median_lrr  = "median",
            std_lrr     = "std",
            min_lrr     = "min",
            max_lrr     = "max",
        )
        .reset_index()
    )

    # Add pct_eroding
    def pct_eroding(group):
        return (group < 0).sum() / len(group) * 100 if len(group) > 0 else np.nan

    eroding = df_valid.groupby("domain_number")["lrr_m_yr"].apply(pct_eroding).reset_index()
    eroding.columns = ["domain_number", "pct_eroding"]
    summary = summary.merge(eroding, on="domain_number", how="left")

    # Add total transect count (including those with too few obs)
    total = df.groupby("domain_number")["transect_id"].count().reset_index()
    total.columns = ["domain_number", "n_transects"]
    summary = summary.merge(total, on="domain_number", how="left")

    # Round for readability
    for col in ["mean_lrr", "median_lrr", "std_lrr", "min_lrr", "max_lrr", "pct_eroding"]:
        summary[col] = summary[col].round(3)

    return summary.sort_values("domain_number").reset_index(drop=True)


def plot_domain_lrr(summary: pd.DataFrame, start: str, end: str,
                    out_path: str, metric: str = "mean_lrr"):
    """
    Bar chart of domain-level LRR coloured by magnitude.
    metric: 'mean_lrr' or 'median_lrr'
    """
    df = summary.dropna(subset=[metric]).sort_values("domain_number")

    # Diverging colourmap centred on zero
    norm   = mcolors.TwoSlopeNorm(vmin=df[metric].min(),
                                   vcenter=0,
                                   vmax=max(df[metric].max(), 0.01))
    cmap   = cm.RdBu
    colors = [cmap(norm(v)) for v in df[metric]]

    fig, axes = plt.subplots(2, 1, figsize=(14, 9),
                              gridspec_kw={"height_ratios": [3, 1]})

    # --- Top: bar chart ---
    ax = axes[0]
    bars = ax.bar(df["domain_number"].astype(str), df[metric],
                  color=colors, edgecolor="none", width=0.8)
    ax.axhline(0, color="black", lw=1.0, zorder=5)

    # Error bars (±1 std) if showing mean
    if metric == "mean_lrr" and "std_lrr" in df.columns:
        ax.errorbar(range(len(df)), df["mean_lrr"],
                    yerr=df["std_lrr"],
                    fmt="none", color="black", capsize=3, lw=0.8, zorder=6)

    ax.set_xlabel("CASCADE Domain Number", fontsize=11)
    ax.set_ylabel("LRR (m/yr)", fontsize=11)
    ax.set_title(f"Domain-Level Linear Regression Rates  [{start} – {end}]\n"
                 f"({'Mean' if metric == 'mean_lrr' else 'Median'} across CoastSat transects per domain)",
                 fontsize=12)
    ax.tick_params(axis="x", rotation=90, labelsize=8)
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_xlim(-0.5, len(df) - 0.5)

    # Colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation="vertical", fraction=0.02, pad=0.02)
    cbar.set_label("LRR (m/yr)", fontsize=9)

    # --- Bottom: n_valid per domain ---
    ax2 = axes[1]
    ax2.bar(df["domain_number"].astype(str), df["n_valid"],
            color="steelblue", edgecolor="none", width=0.8, alpha=0.7)
    ax2.set_xlabel("CASCADE Domain Number", fontsize=10)
    ax2.set_ylabel("N transects", fontsize=10)
    ax2.tick_params(axis="x", rotation=90, labelsize=8)
    ax2.grid(True, axis="y", alpha=0.3)
    ax2.set_xlim(-0.5, len(df) - 0.5)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.show()
    print(f"Domain LRR plot saved: {out_path}")


def plot_transect_scatter(transect_df: pd.DataFrame,
                          buffer_domains: list, out_path: str):
    """
    Scatter of individual transect LRRs coloured by domain,
    sorted by domain number.  Good for seeing within-domain spread.
    """
    df = transect_df.copy()
    if buffer_domains:
        df = df[~df["domain_number"].isin(buffer_domains)]
    df = df[df["lrr_m_yr"].notna()].sort_values("domain_number")

    cmap   = cm.tab20
    domains = df["domain_number"].unique()
    d_norm  = {d: i / max(len(domains) - 1, 1) for i, d in enumerate(sorted(domains))}
    colors  = [cmap(d_norm[d]) for d in df["domain_number"]]

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.scatter(range(len(df)), df["lrr_m_yr"],
               c=colors, s=25, zorder=3, alpha=0.8)
    ax.axhline(0, color="black", lw=1.0)

    # Mark domain boundaries
    prev_d = None
    for i, (_, row) in enumerate(df.iterrows()):
        if row["domain_number"] != prev_d:
            if prev_d is not None:
                ax.axvline(i - 0.5, color="grey", lw=0.5, ls="--", zorder=1)
            prev_d = row["domain_number"]

    ax.set_xlabel("Transect (ordered by domain)", fontsize=11)
    ax.set_ylabel("LRR (m/yr)", fontsize=11)
    ax.set_title("Individual Transect LRRs by Domain", fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.show()


# ============================================================
# MAIN
# ============================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ---- Load lookup table ----
    print(f"Loading lookup table: {LOOKUP_CSV}")
    lookup = pd.read_csv(LOOKUP_CSV)
    print(f"  {len(lookup)} transects, {lookup['domain_number'].nunique()} unique domains\n")

    # ---- Map CSVs ----
    csv_map = collect_csv_map(ROOT_DATA_DIR, SITE_FILTER)
    print(f"Found {len(csv_map)} time-series CSVs in data folders.\n")

    # ---- Compute LRR for all transects ----
    print(f"Computing LRR for period: {START_DATE} → {END_DATE}")
    print("-" * 55)
    transect_df = compute_all_lrr(lookup, csv_map, START_DATE, END_DATE, MIN_OBS)

    valid = transect_df["lrr_m_yr"].notna().sum()
    print(f"\nTransect LRR complete: {valid}/{len(transect_df)} with valid results.")

    # ---- Domain summary ----
    summary = domain_summary(transect_df, BUFFER_DOMAINS)

    print(f"\n{'='*55}")
    print(f"DOMAIN-LEVEL SUMMARY  ({START_DATE} – {END_DATE})")
    print(f"{'='*55}")
    print(summary.to_string(index=False))
    print(f"{'='*55}")

    # ---- Overall stats ----
    valid_s = summary.dropna(subset=["mean_lrr"])
    print(f"\nStudy-area mean LRR : {valid_s['mean_lrr'].mean():+.3f} m/yr")
    print(f"Eroding domains     : {(valid_s['mean_lrr'] < 0).sum()} / {len(valid_s)}")
    print(f"Accreting domains   : {(valid_s['mean_lrr'] > 0).sum()} / {len(valid_s)}")

    # ---- Save outputs ----
    transect_out = os.path.join(OUTPUT_DIR, "transect_lrr_full.csv")
    summary_out  = os.path.join(OUTPUT_DIR, "domain_lrr_summary.csv")
    transect_df.to_csv(transect_out, index=False)
    summary.to_csv(summary_out, index=False)
    print(f"\nSaved: {transect_out}")
    print(f"Saved: {summary_out}")

    # ---- Plots ----
    plot_domain_lrr(summary, START_DATE, END_DATE,
                    os.path.join(OUTPUT_DIR, "domain_lrr_bar.png"),
                    metric="mean_lrr")

    plot_transect_scatter(transect_df, BUFFER_DOMAINS,
                          os.path.join(OUTPUT_DIR, "transect_lrr_scatter.png"))

    return summary, transect_df


# ============================================================
if __name__ == "__main__":
    summary, transect_df = main()
