"""
CoastSat Transect-Level LRR Plot
=================================
Compute and plot Linear Regression Rates (LRR) for a user-specified
list of CoastSat transect IDs over any date range.

No lookup table or domain mapping required — just point it at your
timeseries folder, list the transects you care about, and run.

Supports two date-filtering modes (same as coastsat_domain_lrr_specific_dates.py):
  MODE A — SURVEY_DATES : collect imagery within ±WINDOW_DAYS of specific dates
            (best for DSAS-comparable analysis anchored to survey dates)
  MODE B — START/END    : all imagery in a continuous date window
            (best for general trend analysis)

Outputs
-------
  transect_lrr_results.csv   – LRR, R², p-value, uncertainty, n_obs per transect
  transect_lrr_bar.png       – bar chart coloured by erosion/accretion magnitude
  transect_lrr_timeseries/   – (optional) one time-series PNG per transect

Usage
-----
Edit the CONFIG section below, then run:
    python coastsat_transect_lrr_plot.py
"""

# ============================================================
# CONFIG  –  edit these values before running
# ============================================================

# ---- Transect IDs to process ----
#
# TWO WAYS to specify transects — use whichever suits your case:
#
# OPTION 1 — RANGE (recommended for a contiguous stretch of coast)
#   Give the southernmost and northernmost transect IDs.
#   Both must share the same site prefix (e.g. usa_NC_0032).
#   The script expands the numeric suffix into the full range inclusively.
#
#   Example: IDs 85–112 on site usa_NC_0032
TRANSECT_RANGE_START = "usa_NC_0032_0085"   # southernmost transect
TRANSECT_RANGE_END   = "usa_NC_0032_0112"   # northernmost transect
#
# OPTION 2 — EXPLICIT LIST (for non-contiguous or multi-site selections)
#   Populate TRANSECT_IDS_EXPLICIT with exact filename stems.
#   Set both TRANSECT_RANGE_START and TRANSECT_RANGE_END to "" to use this.
TRANSECT_IDS_EXPLICIT = [
    # "usa_NC_0033_0002",
    # "usa_NC_0034_0005",
]
#
# If TRANSECT_RANGE_START is set, OPTION 1 takes precedence over OPTION 2.

# ---- Data location ----
# Root folder containing site subfolders (e.g. usa_NC_0032_timeseries/).
# The script searches all subfolders one level deep.
ROOT_DATA_DIR = r"/scripts/input_prep/CoastSat/coastsat_timeseries"

# Optional: only search subfolders whose names contain this string.
# Set to "" to search ALL subfolders under ROOT_DATA_DIR.
SITE_FILTER = "usa_NC"

# ---- Date filter (choose MODE A or MODE B) ----

# MODE A — exact survey dates (recommended for DSAS comparison)
# Set to [] to fall back to MODE B.
SURVEY_DATES = []   # e.g. ["1997-09-15", "2008-04-22", "2019-10-03"]
WINDOW_DAYS  = 90   # ±days around each survey date

# MODE B — continuous date range (used when SURVEY_DATES = [])
START_DATE = "2000-01-01"
END_DATE   = "2025-12-31"

# ---- Quality control ----
# Minimum observations required to compute a valid LRR.
MIN_OBS = 3

# ---- Plot options ----
# Label shown on the x-axis under each bar. Options:
#   "full"   → full transect ID (e.g. usa_NC_0033_0002)
#   "short"  → last two underscore-separated parts (e.g. 0033_0002)
#   "index"  → numeric index (1, 2, 3...)
X_LABEL_STYLE = "short"

# If True, save a time-series PNG for every transect in a subfolder.
PLOT_TIMESERIES = True

# Show each bar annotated with its LRR value (m/yr)?
ANNOTATE_BARS = True

# ---- Output ----
OUTPUT_DIR = r"/scripts/input_prep/5-scr/CoastSat\custom\buxton_2000_2025"

# ============================================================
# IMPORTS
# ============================================================
import os
import sys
import glob
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm

warnings.filterwarnings("ignore")

# Import core LRR functions from the companion analysis script.
# Assumes coastsat_lrr_analysis.py is in the same directory.
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)
from coastsat_lrr_analysis import (
    load_timeseries,
    filter_dates,
    filter_to_dates,
    compute_lrr,
    _empty_lrr,
)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def resolve_transect_ids() -> list:
    """
    Build the final list of transect IDs from the CONFIG above.

    If TRANSECT_RANGE_START is set, expands the numeric suffix range
    (inclusive) into a full list. Both endpoints must share the same
    site prefix (e.g. 'usa_NC_0032').

    Falls back to TRANSECT_IDS_EXPLICIT if no range is specified.
    """
    if TRANSECT_RANGE_START and TRANSECT_RANGE_END:
        # Split into prefix and zero-padded numeric suffix
        # e.g. "usa_NC_0032_0085"  →  prefix="usa_NC_0032", start_num=85
        def split_id(tid):
            parts = tid.rsplit("_", 1)
            if len(parts) != 2:
                raise ValueError(f"Cannot parse transect ID for range: '{tid}'. "
                                 f"Expected format like 'usa_NC_0032_0085'.")
            prefix, suffix = parts
            return prefix, int(suffix), len(suffix)  # keep zero-pad width

        prefix_s, start_num, width_s = split_id(TRANSECT_RANGE_START)
        prefix_e, end_num,   width_e = split_id(TRANSECT_RANGE_END)

        if prefix_s != prefix_e:
            raise ValueError(
                f"TRANSECT_RANGE_START and TRANSECT_RANGE_END must share the same "
                f"site prefix.\n  Start: '{TRANSECT_RANGE_START}' → prefix '{prefix_s}'\n"
                f"  End:   '{TRANSECT_RANGE_END}' → prefix '{prefix_e}'"
            )

        width = max(width_s, width_e)
        lo, hi = min(start_num, end_num), max(start_num, end_num)
        ids = [f"{prefix_s}_{str(n).zfill(width)}" for n in range(lo, hi + 1)]
        print(f"Range resolved: {TRANSECT_RANGE_START} → {TRANSECT_RANGE_END}  "
              f"({len(ids)} transects)")
        return ids

    elif TRANSECT_IDS_EXPLICIT:
        print(f"Using explicit list: {len(TRANSECT_IDS_EXPLICIT)} transect(s)")
        return list(TRANSECT_IDS_EXPLICIT)

    else:
        raise ValueError("No transects specified. Set TRANSECT_RANGE_START/END "
                         "or populate TRANSECT_IDS_EXPLICIT in the CONFIG.")


def discover_csvs(root_dir: str, site_filter: str = "") -> dict:
    """
    Walk one level of subfolders under root_dir and collect every CSV.

    Returns
    -------
    dict  { transect_id_stem : full_filepath }
    e.g.  { 'usa_NC_0033_0002' : 'C:/…/usa_NC_0033_timeseries/usa_NC_0033_0002.csv' }
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

    for sf in subfolders:
        for fpath in glob.glob(os.path.join(sf, "*.csv")):
            stem = os.path.splitext(os.path.basename(fpath))[0]
            csv_map[stem] = fpath

    print(f"Discovered {len(csv_map)} CSV(s) across {len(subfolders)} subfolder(s).")
    return csv_map


def apply_filter(df_raw: pd.DataFrame) -> pd.DataFrame:
    """Apply the date filter configured above (MODE A or MODE B)."""
    if SURVEY_DATES:
        return filter_to_dates(df_raw, SURVEY_DATES, WINDOW_DAYS)
    else:
        return filter_dates(df_raw, START_DATE, END_DATE)


def make_x_label(transect_id: str, idx: int) -> str:
    """Format the x-axis tick label for a transect based on X_LABEL_STYLE."""
    if X_LABEL_STYLE == "short":
        parts = transect_id.rsplit("_", 2)
        return "_".join(parts[-2:]) if len(parts) >= 2 else transect_id
    elif X_LABEL_STYLE == "index":
        return str(idx + 1)
    else:  # "full"
        return transect_id


def period_label() -> str:
    """Human-readable description of the active date filter."""
    if SURVEY_DATES:
        dates_str = ", ".join(SURVEY_DATES)
        return f"±{WINDOW_DAYS}-day windows around [{dates_str}]"
    return f"{START_DATE} – {END_DATE}"


# ============================================================
# COMPUTATION
# ============================================================

def compute_transect_lrr(transect_ids: list, csv_map: dict) -> pd.DataFrame:
    """
    Compute LRR for each requested transect ID.

    Returns a DataFrame with one row per transect:
        transect_id, lrr_m_yr, unc_m_yr, r_squared, p_value,
        n_obs, start_date, end_date, status
    """
    records = []

    for i, tid in enumerate(transect_ids):
        if tid not in csv_map:
            print(f"  MISSING  {tid}  (no CSV found)")
            records.append({"transect_id": tid, "status": "missing", **_empty_lrr(0)})
            continue

        try:
            df_raw  = load_timeseries(csv_map[tid])
            df_filt = apply_filter(df_raw)

            if len(df_filt) < MIN_OBS:
                print(f"  SKIP     {tid}  ({len(df_filt)} obs < MIN_OBS={MIN_OBS})")
                records.append({"transect_id": tid, "status": f"too_few_obs ({len(df_filt)})",
                                 **_empty_lrr(len(df_filt))})
            else:
                result = compute_lrr(df_filt)
                status = "ok"
                print(f"  OK       {tid}  "
                      f"n={result['n_obs']:4d}  "
                      f"LRR={result['lrr_m_yr']:+7.3f} m/yr  "
                      f"R²={result['r_squared']:.3f}  "
                      f"p={result['p_value']:.4f}")
                records.append({"transect_id": tid, "status": status, **result})

        except Exception as e:
            print(f"  ERROR    {tid}: {e}")
            records.append({"transect_id": tid, "status": f"error: {e}", **_empty_lrr(0)})

    return pd.DataFrame(records)


# ============================================================
# PLOTS
# ============================================================

def plot_lrr_bar(results: pd.DataFrame, out_path: str):
    """
    Horizontal bar chart of LRR per transect, coloured by magnitude.
    Error bars show 95% confidence interval (±unc_m_yr).
    """
    df = results.dropna(subset=["lrr_m_yr"]).copy()
    if df.empty:
        print("⚠️  No valid LRR values to plot.")
        return

    # Build x-axis labels
    df["x_label"] = [
        make_x_label(tid, i)
        for i, tid in enumerate(df["transect_id"])
    ]

    # Diverging colourmap centred on zero
    vmin = df["lrr_m_yr"].min()
    vmax = df["lrr_m_yr"].max()
    # Guard against all-same-sign results
    vcenter = 0.0
    if vmin >= 0:
        vmin = -0.001
    if vmax <= 0:
        vmax = 0.001
    norm   = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    cmap   = cm.RdBu
    colors = [cmap(norm(v)) for v in df["lrr_m_yr"]]

    fig, ax = plt.subplots(figsize=(max(8, len(df) * 0.55 + 2), 6))

    x = np.arange(len(df))
    bars = ax.bar(x, df["lrr_m_yr"], color=colors, edgecolor="none",
                  width=0.7, zorder=3)

    # Error bars (±unc from 95% CI)
    if "unc_m_yr" in df.columns:
        ax.errorbar(x, df["lrr_m_yr"],
                    yerr=df["unc_m_yr"].fillna(0),
                    fmt="none", color="black", capsize=4, lw=1.0, zorder=4)

    # Zero line
    ax.axhline(0, color="black", lw=1.2, zorder=5)

    # Optional bar annotations
    if ANNOTATE_BARS:
        for xi, (_, row) in zip(x, df.iterrows()):
            v = row["lrr_m_yr"]
            offset = (vmax - vmin) * 0.025
            va = "bottom" if v >= 0 else "top"
            y_pos = v + (offset if v >= 0 else -offset)
            ax.text(xi, y_pos, f"{v:+.2f}", ha="center", va=va,
                    fontsize=7.5, color="black")

    # Axes formatting
    ax.set_xticks(x)
    ax.set_xticklabels(df["x_label"], rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("LRR (m/yr)", fontsize=11)
    ax.set_xlabel("Transect", fontsize=11)
    ax.set_title(f"Shoreline Change Rate (LRR)\n[{period_label()}]", fontsize=12)
    ax.grid(True, axis="y", alpha=0.3, zorder=1)
    ax.set_xlim(-0.6, len(df) - 0.4)

    # Colourbar
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation="vertical", fraction=0.02, pad=0.03)
    cbar.set_label("LRR (m/yr)", fontsize=9)

    # Summary stats in a text box
    mean_lrr   = df["lrr_m_yr"].mean()
    n_eroding  = (df["lrr_m_yr"] < 0).sum()
    n_accreting = (df["lrr_m_yr"] >= 0).sum()
    stats_text = (f"n={len(df)}  |  mean={mean_lrr:+.2f} m/yr\n"
                  f"eroding: {n_eroding}  |  accreting: {n_accreting}")
    ax.text(0.01, 0.99, stats_text, transform=ax.transAxes,
            va="top", ha="left", fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Bar chart saved: {out_path}")


def plot_all_timeseries(results: pd.DataFrame, csv_map: dict, out_dir: str):
    """
    Save one time-series PNG per transect showing raw observations,
    the filtered window, and the fitted LRR trend line.
    """
    os.makedirs(out_dir, exist_ok=True)
    saved = 0

    for _, row in results.iterrows():
        tid = row["transect_id"]
        if tid not in csv_map:
            continue

        try:
            df_raw  = load_timeseries(csv_map[tid])
            df_filt = apply_filter(df_raw)

            fig, ax = plt.subplots(figsize=(10, 4))

            # All observations (background)
            ax.scatter(df_raw["date"], df_raw["chainage_m"],
                       color="lightgray", s=15, zorder=1, label="All obs.")

            # Filtered observations
            ax.scatter(df_filt["date"], df_filt["chainage_m"],
                       color="steelblue", s=25, zorder=2, label="Selected obs.")

            # Trend line
            if not np.isnan(row["lrr_m_yr"]) and len(df_filt) >= 2:
                t0   = df_filt["date"].min()
                yrs  = (df_filt["date"] - t0).dt.total_seconds() / (86400.0 * 365.25)
                coeffs = np.polyfit(yrs.values, df_filt["chainage_m"].values, 1)
                y_hat  = np.polyval(coeffs, yrs)
                label  = (f"LRR = {row['lrr_m_yr']:+.2f} ± {row['unc_m_yr']:.2f} m/yr  "
                          f"(R²={row['r_squared']:.2f}, n={int(row['n_obs'])})")
                ax.plot(df_filt["date"], y_hat, color="crimson", lw=2.0,
                        zorder=3, label=label)

            ax.set_title(f"Transect: {tid}", fontsize=11)
            ax.set_xlabel("Date")
            ax.set_ylabel("Cross-shore distance (m)")
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)

            # Shade the active filter window
            if not SURVEY_DATES:
                ax.axvspan(
                    pd.Timestamp(START_DATE, tz="UTC"),
                    pd.Timestamp(END_DATE,   tz="UTC"),
                    alpha=0.08, color="steelblue", zorder=0, label="_nolegend_"
                )

            plt.tight_layout()
            fname = os.path.join(out_dir, f"{tid}_timeseries.png")
            plt.savefig(fname, dpi=120, bbox_inches="tight")
            plt.close()
            saved += 1

        except Exception as e:
            print(f"  Time-series plot failed for {tid}: {e}")

    print(f"Time-series plots saved: {saved} → {out_dir}")


# ============================================================
# MAIN
# ============================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ---- Resolve transect list from config ----
    transect_ids = resolve_transect_ids()

    print("=" * 60)
    print("CoastSat Transect LRR Plot")
    print(f"Period    : {period_label()}")
    print(f"Transects : {len(transect_ids)} requested")
    print(f"MIN_OBS   : {MIN_OBS}")
    print("=" * 60)

    # ---- Discover CSVs ----
    csv_map = discover_csvs(ROOT_DATA_DIR, SITE_FILTER)

    # Check which requested transects are available before running
    found    = [t for t in transect_ids if t in csv_map]
    missing  = [t for t in transect_ids if t not in csv_map]
    if missing:
        print(f"\n⚠️  {len(missing)} transect(s) not found in data folders:")
        for m in missing:
            print(f"     {m}")
    print(f"\nProcessing {len(found)}/{len(transect_ids)} transect(s)...\n")

    # ---- Compute LRR ----
    print("-" * 60)
    results = compute_transect_lrr(transect_ids, csv_map)
    print("-" * 60)

    # ---- Summary ----
    valid = results.dropna(subset=["lrr_m_yr"])
    print(f"\n{'='*60}")
    print(f"RESULTS  ({period_label()})")
    print(f"{'='*60}")
    display_cols = ["transect_id", "lrr_m_yr", "unc_m_yr",
                    "r_squared", "p_value", "n_obs", "status"]
    print(results[[c for c in display_cols if c in results.columns]].to_string(index=False))
    print(f"{'='*60}")

    if len(valid):
        print(f"\nMean LRR   : {valid['lrr_m_yr'].mean():+.3f} m/yr")
        print(f"Median LRR : {valid['lrr_m_yr'].median():+.3f} m/yr")
        print(f"Eroding    : {(valid['lrr_m_yr'] < 0).sum()} / {len(valid)}")
        print(f"Accreting  : {(valid['lrr_m_yr'] >= 0).sum()} / {len(valid)}")

    # ---- Save CSV ----
    csv_out = os.path.join(OUTPUT_DIR, "transect_lrr_results.csv")
    results.to_csv(csv_out, index=False)
    print(f"\nSaved: {csv_out}")

    # ---- Bar chart ----
    bar_out = os.path.join(OUTPUT_DIR, "transect_lrr_bar.png")
    plot_lrr_bar(results, bar_out)

    # ---- Time-series plots (optional) ----
    if PLOT_TIMESERIES:
        ts_dir = os.path.join(OUTPUT_DIR, "transect_lrr_timeseries")
        plot_all_timeseries(results, csv_map, ts_dir)

    return results


# ============================================================
if __name__ == "__main__":
    results = main()
