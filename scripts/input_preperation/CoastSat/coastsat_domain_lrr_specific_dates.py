"""
CoastSat Domain-Level LRR Summary
===================================
Step 2 of 2 in the domain-level LRR workflow.

Requires:
  1. transect_domain_lookup.csv   – from coastsat_domain_mapping.py
  2. CoastSat time-series CSVs    – one per transect
  3. coastsat_lrr_analysis_old.py     – in the same directory (or on PYTHONPATH)

Outputs:
  domain_lrr_summary.csv  –  one row per domain with aggregated LRR stats
  domain_lrr_plot.png     –  bar chart coloured by erosion/accretion
  transect_lrr_full.csv   –  full transect-level results with domain assignments

Usage
-----
Edit the CONFIG section below, then run:
    python coastsat_domain_lrr_fixed.py

FILTER MODE GUIDE
-----------------
This script supports two mutually exclusive filtering modes.

  MODE A — MATCH_YEARS (recommended for DSAS comparison)
    Set MATCH_YEARS to the calendar years used in your DSAS analysis.
    CoastSat will only use observations from those years, giving a
    like-for-like comparison with DSAS (which uses a small set of
    digitized USGS shorelines from fixed survey years).

    Example — Period 1997–2019, USGS shorelines from 1997 and 2019:
        MATCH_YEARS = [1997, 2019]

    Example — Period 1978–1997, USGS shorelines from 1978 and 1997:
        MATCH_YEARS = [1978, 1997]

    If a transect has no CoastSat imagery in one of the target years it
    will still use whatever observations exist in the remaining years.
    The field 'n_obs' in the output tells you how many points were used.

    Optional: set YEAR_WINDOW_DAYS > 0 to expand each target year into a
    ±N-day window around Jan 1. For example YEAR_WINDOW_DAYS = 182
    captures the full calendar year (±182 ≈ ±6 months). Leave at 0 to
    use whole calendar years (Jan 1 – Dec 31).

  MODE B — DATE RANGE (original behaviour, all imagery in range)
    Leave MATCH_YEARS = [] and set START_DATE / END_DATE.
    CoastSat uses every available observation within the date range.
    Good for producing dense-observation LRRs, but not directly
    comparable to DSAS when DSAS is anchored to specific survey years.
"""

# ============================================================
# CONFIG
# ============================================================

# Path to the lookup table produced by coastsat_domain_mapping.py
LOOKUP_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\transect_domain_lookup.csv"

# Root folder containing all site subfolders
ROOT_DATA_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\coastsat_timeseries"

# Only include subfolders whose names contain this string. "" = all.
SITE_FILTER = "usa_NC"

# ------------------------------------------------------------------
# FILTER MODE — choose A or B (see module docstring above)
# ------------------------------------------------------------------

# MODE A: Survey-date-matched filter (recommended for DSAS comparison)
#
# List the EXACT acquisition dates of the USGS shorelines you digitized
# in DSAS — one string per shoreline, format "YYYY-MM-DD".
# CoastSat will collect imagery within ±WINDOW_DAYS of each date.
#
# How to find these dates:
#   Open your DSAS shoreline feature class / shapefile and check the
#   date attribute (often "Date_", "ShorDate", or similar) for each
#   shoreline.  Use the actual image/survey acquisition date, not the
#   year of the DSAS analysis period.
#
# Examples:
#   Period 1997–2019, 3 USGS shorelines:
#     SURVEY_DATES = ["1997-09-15", "2008-04-22", "2019-10-03"]
#
#   Period 1978–1997, 3 USGS shorelines:
#     SURVEY_DATES = ["1978-06-10", "1986-08-01", "1997-09-15"]
#
# Set SURVEY_DATES = [] to fall back to MODE B (continuous date range).
SURVEY_DATES = [] #["1997-09-27", "2009-08-19", "2019-12-08"]   # <-- replace with your actual dates

# Half-width of the search window around each survey date, in days.
# ±30 days (~1 month) is a reasonable default.
# Increase (e.g. 60 or 90) if CoastSat has sparse coverage due to
# cloud cover or satellite availability in your region.
WINDOW_DAYS = 90

# MODE B: Continuous date range — only used when SURVEY_DATES = []
START_DATE = "1997-01-01"
END_DATE   = "2019-12-31"

# ------------------------------------------------------------------

# Minimum observations per transect to include it in domain summaries.
# Note: with MATCH_YEARS filtering you may get as few as 1–2 obs per
# target year. Recommended minimum for a meaningful LRR is 2; lower
# values will produce unreliable slopes. Adjust to your data density.
MIN_OBS = 2

# Output directory
OUTPUT_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\1997_2019_specific_dates"

# CASCADE buffer domains to EXCLUDE from summaries
BUFFER_DOMAINS = []

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
import warnings
warnings.filterwarnings("ignore")

# Import LRR functions from the companion script
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)
from coastsat_lrr_analysis import (
    load_timeseries, filter_dates, filter_to_dates, filter_to_years, compute_lrr, _empty_lrr
)


# ============================================================
# FUNCTIONS
# ============================================================

def collect_csv_map(root_dir: str, site_filter: str = "") -> dict:
    """
    Auto-discover all time-series CSVs under root_dir.

    Returns:
        { transect_id_stem : full_filepath }
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


def apply_filter(df_raw: pd.DataFrame,
                 survey_dates: list,
                 window_days: int,
                 start: str,
                 end: str) -> pd.DataFrame:
    """
    Apply the appropriate filter based on config.

    If survey_dates is non-empty → exact-date window filter (MODE A).
    Otherwise → continuous date range filter (MODE B).
    """
    if survey_dates:
        return filter_to_dates(df_raw, survey_dates, window_days)
    else:
        return filter_dates(df_raw, start, end)


def compute_all_lrr(lookup: pd.DataFrame,
                    csv_map: dict,
                    survey_dates: list,
                    window_days: int,
                    start: str,
                    end: str,
                    min_obs: int) -> pd.DataFrame:
    """
    For every transect in the lookup table, find its CSV, apply the
    configured filter, compute LRR, and return a merged DataFrame.
    """
    records = []
    missing = []

    if survey_dates:
        filter_label = f"±{window_days} days around {survey_dates}"
    else:
        filter_label = f"{start} → {end}"

    for _, row in lookup.iterrows():
        tid = str(row["transect_id"])
        date_counts = {}   # obs count per survey-date anchor for this transect

        if tid not in csv_map:
            missing.append(tid)
            result = _empty_lrr(0)
        else:
            try:
                df_raw  = load_timeseries(csv_map[tid])
                df_filt = apply_filter(df_raw, survey_dates, window_days, start, end)

                # Count images found per survey-date anchor (only when date-matching)
                if survey_dates and "survey_date" in df_filt.columns:
                    date_counts = df_filt.groupby("survey_date").size().to_dict()
                    # Ensure all anchors appear even if 0 images found
                    for sd in survey_dates:
                        date_counts.setdefault(sd, 0)

                if len(df_filt) < min_obs:
                    result = _empty_lrr(len(df_filt))
                else:
                    result = compute_lrr(df_filt)
            except Exception as e:
                print(f"  ERROR {tid}: {e}")
                result = _empty_lrr(0)

        # Prefix per-date counts so they're easy to identify in the CSV
        date_count_cols = {f"n_obs_{sd}": date_counts.get(sd, 0)
                           for sd in survey_dates} if survey_dates else {}

        records.append({
            "transect_id"  : tid,
            "domain_number": row["domain_number"],
            "match_method" : row.get("match_method", ""),
            **result,
            **date_count_cols,
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
      n_transects  – total transects in domain
      n_valid      – transects with valid LRR
      mean_lrr     – mean LRR across transects (m/yr)
      median_lrr   – median LRR
      std_lrr      – standard deviation
      min_lrr      – most erosional transect
      max_lrr      – most accretionary transect
      pct_eroding  – % of transects with negative LRR
    """
    df = transect_df.copy()
    if buffer_domains:
        df = df[~df["domain_number"].isin(buffer_domains)]

    df_valid = df[df["lrr_m_yr"].notna()]

    summary = (
        df_valid.groupby("domain_number")["lrr_m_yr"]
        .agg(
            n_valid    = "count",
            mean_lrr   = "mean",
            median_lrr = "median",
            std_lrr    = "std",
            min_lrr    = "min",
            max_lrr    = "max",
        )
        .reset_index()
    )

    def pct_eroding(group):
        return (group < 0).sum() / len(group) * 100 if len(group) > 0 else np.nan

    eroding = df_valid.groupby("domain_number")["lrr_m_yr"].apply(pct_eroding).reset_index()
    eroding.columns = ["domain_number", "pct_eroding"]
    summary = summary.merge(eroding, on="domain_number", how="left")

    total = df.groupby("domain_number")["transect_id"].count().reset_index()
    total.columns = ["domain_number", "n_transects"]
    summary = summary.merge(total, on="domain_number", how="left")

    # Mean images per transect — this is what actually changes between window
    # sizes, unlike n_valid which only reflects how many transects have any data
    mean_nobs = df_valid.groupby("domain_number")["n_obs"].mean().reset_index()
    mean_nobs.columns = ["domain_number", "mean_n_obs"]
    mean_nobs["mean_n_obs"] = mean_nobs["mean_n_obs"].round(1)
    summary = summary.merge(mean_nobs, on="domain_number", how="left")

    # Per-survey-date image counts (columns named n_obs_YYYY-MM-DD)
    date_cols = [c for c in df.columns if c.startswith("n_obs_")]
    for dc in date_cols:
        agg = df_valid.groupby("domain_number")[dc].mean().reset_index()
        agg.columns = ["domain_number", f"mean_{dc}"]
        agg[f"mean_{dc}"] = agg[f"mean_{dc}"].round(1)
        summary = summary.merge(agg, on="domain_number", how="left")

    for col in ["mean_lrr", "median_lrr", "std_lrr", "min_lrr", "max_lrr", "pct_eroding"]:
        summary[col] = summary[col].round(3)

    return summary.sort_values("domain_number").reset_index(drop=True)


def plot_domain_lrr(summary: pd.DataFrame,
                    period_label: str,
                    out_path: str,
                    metric: str = "mean_lrr"):
    """Bar chart of domain-level LRR coloured by magnitude."""
    df = summary.dropna(subset=[metric]).sort_values("domain_number")

    norm   = mcolors.TwoSlopeNorm(vmin=df[metric].min(),
                                   vcenter=0,
                                   vmax=max(df[metric].max(), 0.01))
    cmap   = cm.RdBu
    colors = [cmap(norm(v)) for v in df[metric]]

    fig, axes = plt.subplots(2, 1, figsize=(14, 9),
                              gridspec_kw={"height_ratios": [3, 1]})

    ax = axes[0]
    ax.bar(df["domain_number"].astype(str), df[metric],
           color=colors, edgecolor="none", width=0.8)
    ax.axhline(0, color="black", lw=1.0, zorder=5)

    if metric == "mean_lrr" and "std_lrr" in df.columns:
        ax.errorbar(range(len(df)), df["mean_lrr"],
                    yerr=df["std_lrr"],
                    fmt="none", color="black", capsize=3, lw=0.8, zorder=6)

    ax.set_xlabel("CASCADE Domain Number", fontsize=11)
    ax.set_ylabel("LRR (m/yr)", fontsize=11)
    ax.set_title(f"Domain-Level Linear Regression Rates  [{period_label}]\n"
                 f"({'Mean' if metric == 'mean_lrr' else 'Median'} across CoastSat transects per domain)",
                 fontsize=12)
    ax.tick_params(axis="x", rotation=90, labelsize=8)
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_xlim(-0.5, len(df) - 0.5)

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation="vertical", fraction=0.02, pad=0.02)
    cbar.set_label("LRR (m/yr)", fontsize=9)

    ax2 = axes[1]
    date_cols = sorted([c for c in df.columns if c.startswith("mean_n_obs_")])

    if date_cols:
        # Stacked bars: one colour per survey-date anchor
        palette = ["#4393c3", "#f4a582", "#92c46c", "#d6604d"]
        bottom  = np.zeros(len(df))
        x_labels = df["domain_number"].astype(str)
        for i, dc in enumerate(date_cols):
            label = dc.replace("mean_n_obs_", "")   # e.g. "1997-09-27"
            vals  = df[dc].fillna(0).values
            ax2.bar(x_labels, vals, bottom=bottom,
                    color=palette[i % len(palette)],
                    edgecolor="none", width=0.8, alpha=0.85,
                    label=label)
            bottom += vals
        ax2.legend(fontsize=8, title="Survey date", title_fontsize=8,
                   loc="upper right")
        ax2.set_ylabel("Mean images\nper transect", fontsize=10)
    else:
        # Fallback when not using date-matching: show total n_obs
        ax2.bar(df["domain_number"].astype(str),
                df["mean_n_obs"].fillna(0),
                color="steelblue", edgecolor="none", width=0.8, alpha=0.7)
        ax2.set_ylabel("Mean images\nper transect", fontsize=10)

    ax2.set_xlabel("CASCADE Domain Number", fontsize=10)
    ax2.tick_params(axis="x", rotation=90, labelsize=8)
    ax2.grid(True, axis="y", alpha=0.3)
    ax2.set_xlim(-0.5, len(df) - 0.5)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.show()
    print(f"Domain LRR plot saved: {out_path}")


def plot_transect_scatter(transect_df: pd.DataFrame,
                          buffer_domains: list, out_path: str):
    """Scatter of individual transect LRRs by domain."""
    df = transect_df.copy()
    if buffer_domains:
        df = df[~df["domain_number"].isin(buffer_domains)]
    df = df[df["lrr_m_yr"].notna()].sort_values("domain_number")

    cmap    = cm.tab20
    domains = df["domain_number"].unique()
    d_norm  = {d: i / max(len(domains) - 1, 1) for i, d in enumerate(sorted(domains))}
    colors  = [cmap(d_norm[d]) for d in df["domain_number"]]

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.scatter(range(len(df)), df["lrr_m_yr"],
               c=colors, s=25, zorder=3, alpha=0.8)
    ax.axhline(0, color="black", lw=1.0)

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

    # ---- Determine filter mode and build period label ----
    if SURVEY_DATES:
        dates_str    = ", ".join(SURVEY_DATES)
        period_label = f"±{WINDOW_DAYS}-day windows around [{dates_str}]"
        filter_mode  = "MODE A — survey-date-matched"
    else:
        period_label = f"{START_DATE} – {END_DATE}"
        filter_mode  = "MODE B — continuous date range"

    print("=" * 65)
    print(f"CoastSat Domain LRR  |  {filter_mode}")
    print(f"Filter               :  {period_label}")
    print(f"MIN_OBS              :  {MIN_OBS}")
    print("=" * 65)

    # ---- Load lookup table ----
    print(f"\nLoading lookup table: {LOOKUP_CSV}")
    lookup = pd.read_csv(LOOKUP_CSV)
    print(f"  {len(lookup)} transects, {lookup['domain_number'].nunique()} unique domains\n")

    # ---- Map CSVs ----
    csv_map = collect_csv_map(ROOT_DATA_DIR, SITE_FILTER)
    print(f"\nFound {len(csv_map)} time-series CSVs.\n")

    # ---- Compute LRR for all transects ----
    print(f"Computing LRR — {period_label}")
    print("-" * 65)
    transect_df = compute_all_lrr(
        lookup, csv_map,
        survey_dates = SURVEY_DATES,
        window_days  = WINDOW_DAYS,
        start        = START_DATE,
        end          = END_DATE,
        min_obs      = MIN_OBS,
    )

    valid = transect_df["lrr_m_yr"].notna().sum()
    print(f"\nTransect LRR complete: {valid}/{len(transect_df)} with valid results.")

    # ---- Domain summary ----
    summary = domain_summary(transect_df, BUFFER_DOMAINS)

    print(f"\n{'='*65}")
    print(f"DOMAIN-LEVEL SUMMARY  ({period_label})")
    print(f"{'='*65}")
    print(summary.to_string(index=False))
    print(f"{'='*60}")

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
    plot_domain_lrr(summary, period_label,
                    os.path.join(OUTPUT_DIR, "domain_lrr_bar.png"),
                    metric="mean_lrr")

    plot_transect_scatter(transect_df, BUFFER_DOMAINS,
                          os.path.join(OUTPUT_DIR, "transect_lrr_scatter.png"))

    return summary, transect_df


# ============================================================
if __name__ == "__main__":
    summary, transect_df = main()
