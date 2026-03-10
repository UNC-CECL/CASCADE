"""
CoastSat Transect LRR (Linear Regression Rate) Analysis
========================================================
Loads CoastSat time-series CSVs from one or more folders,
filters to a user-defined date range, and computes the
Linear Regression Rate (LRR) for every transect.

Expected file naming convention:
    <site>_<transect_id>.csv   e.g.  usa_NC_0033_0002.csv

Expected CSV format (CoastSat standard):
    Column 1 – "dates UTC"     : ISO-8601 datetime string
    Column 2 – "chainage (m)"  : cross-shore distance in metres

Usage
-----
Edit the CONFIG section below, then run:
    python coastsat_lrr_analysis.py
"""

# ============================================================
# CONFIG  –  edit these values before running
# ============================================================

# List every folder that contains time-series CSVs.
# All CSVs in these folders will be processed.
DATA_FOLDERS = [
    r"C:/path/to/your/coastsat_data/site1",
    r"C:/path/to/your/coastsat_data/site2",
    # add more as needed…
]

# Date range filter (inclusive).  Use None to include all dates.
START_DATE = "1984-01-01"   # e.g. "1997-01-01"  or  None
END_DATE   = "2026-01-01"   # e.g. "2019-12-31"  or  None

# Minimum number of observations required to compute an LRR.
MIN_OBS = 10

# Output CSV path (set to None to skip saving).
OUTPUT_CSV = "lrr_results.csv"

# ============================================================
# IMPORTS
# ============================================================
import os
import glob
import warnings
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.cm as cm

warnings.filterwarnings("ignore", category=FutureWarning)


# ============================================================
# FUNCTIONS
# ============================================================

def parse_transect_id(filepath: str) -> str:
    """Extract transect ID from filename (everything after last '_')."""
    basename = os.path.splitext(os.path.basename(filepath))[0]
    return basename  # keep full name; change if you need only suffix


def load_timeseries(filepath: str) -> pd.DataFrame:
    """
    Load a CoastSat time-series CSV.
    Returns a DataFrame with columns ['date', 'chainage_m'].
    """
    df = pd.read_csv(filepath, header=0)
    df.columns = [c.strip() for c in df.columns]          # strip whitespace
    date_col    = df.columns[0]
    chainage_col = df.columns[1]

    df = df.rename(columns={date_col: "date", chainage_col: "chainage_m"})
    df["date"] = pd.to_datetime(df["date"], utc=True)
    df["chainage_m"] = pd.to_numeric(df["chainage_m"], errors="coerce")
    df = df.dropna(subset=["date", "chainage_m"]).sort_values("date").reset_index(drop=True)
    return df


def filter_dates(df: pd.DataFrame, start: str | None, end: str | None) -> pd.DataFrame:
    """Clip DataFrame to [start, end] date range (inclusive)."""
    if start:
        df = df[df["date"] >= pd.Timestamp(start, tz="UTC")]
    if end:
        df = df[df["date"] <= pd.Timestamp(end, tz="UTC")]
    return df.reset_index(drop=True)


def compute_lrr(df: pd.DataFrame) -> dict:
    """
    Compute Linear Regression Rate (LRR) from a time-series DataFrame.

    Returns a dict with:
        lrr_m_yr   – slope in m/yr
        r_squared  – R² of the regression
        p_value    – p-value of the slope
        unc_m_yr   – 95 % confidence interval half-width (m/yr)
        n_obs      – number of observations used
        start_date – earliest date in filtered series
        end_date   – latest date in filtered series
    """
    if len(df) < 2:
        return _empty_lrr(len(df))

    # Convert dates to decimal years for regression
    t0   = df["date"].min()
    days = (df["date"] - t0).dt.total_seconds() / 86400.0
    yrs  = days / 365.25

    x = yrs.values
    y = df["chainage_m"].values

    slope, intercept, r, p, se = stats.linregress(x, y)

    # 95 % CI half-width  (t * SE)
    from scipy.stats import t as t_dist
    dof = len(x) - 2
    t95 = t_dist.ppf(0.975, dof) if dof > 0 else np.nan
    unc = t95 * se

    return {
        "lrr_m_yr"  : round(slope, 4),
        "r_squared" : round(r**2, 4),
        "p_value"   : round(p, 6),
        "unc_m_yr"  : round(unc, 4),
        "n_obs"     : len(x),
        "start_date": df["date"].min().strftime("%Y-%m-%d"),
        "end_date"  : df["date"].max().strftime("%Y-%m-%d"),
    }


def _empty_lrr(n: int) -> dict:
    return dict(lrr_m_yr=np.nan, r_squared=np.nan, p_value=np.nan,
                unc_m_yr=np.nan, n_obs=n, start_date=None, end_date=None)


def collect_csv_files(folders: list[str]) -> list[str]:
    """Gather all CSV files from the provided folder list."""
    files = []
    for folder in folders:
        found = glob.glob(os.path.join(folder, "*.csv"))
        files.extend(found)
    if not files:
        print("⚠️  No CSV files found. Check your DATA_FOLDERS paths.")
    return sorted(files)


# ============================================================
# PLOTTING HELPERS
# ============================================================

def plot_timeseries(df_raw: pd.DataFrame, df_filtered: pd.DataFrame,
                    lrr_result: dict, transect_id: str):
    """Plot raw + filtered observations with the LRR trend line."""
    fig, ax = plt.subplots(figsize=(10, 4))

    # All observations (faded)
    ax.scatter(df_raw["date"], df_raw["chainage_m"],
               color="lightgray", s=15, zorder=1, label="All obs.")

    # Filtered observations
    ax.scatter(df_filtered["date"], df_filtered["chainage_m"],
               color="steelblue", s=20, zorder=2, label="Selected period")

    # LRR trend line
    if not np.isnan(lrr_result["lrr_m_yr"]) and len(df_filtered) >= 2:
        t0   = df_filtered["date"].min()
        days = (df_filtered["date"] - t0).dt.total_seconds() / 86400.0
        yrs  = days / 365.25
        y_hat = (np.polyfit(yrs.values, df_filtered["chainage_m"].values, 1)[0] * yrs +
                 np.polyfit(yrs.values, df_filtered["chainage_m"].values, 1)[1])
        ax.plot(df_filtered["date"], y_hat, color="crimson", lw=2,
                label=f"LRR = {lrr_result['lrr_m_yr']:+.2f} m/yr  (R²={lrr_result['r_squared']:.2f})")

    ax.set_title(f"Transect {transect_id}", fontsize=12)
    ax.set_xlabel("Date")
    ax.set_ylabel("Cross-shore distance (m)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_lrr_summary(results_df: pd.DataFrame):
    """Bar chart of LRR values coloured by sign."""
    df = results_df.dropna(subset=["lrr_m_yr"]).sort_values("lrr_m_yr")
    colors = ["steelblue" if v >= 0 else "crimson" for v in df["lrr_m_yr"]]

    fig, ax = plt.subplots(figsize=(max(8, len(df) * 0.25), 5))
    ax.bar(df["transect_id"], df["lrr_m_yr"], color=colors, edgecolor="none")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("Transect")
    ax.set_ylabel("LRR (m/yr)")
    ax.set_title(f"Linear Regression Rates  [{START_DATE} – {END_DATE}]")
    ax.tick_params(axis="x", rotation=90, labelsize=7)
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.show()


# ============================================================
# MAIN
# ============================================================

def main():
    csv_files = collect_csv_files(DATA_FOLDERS)
    print(f"Found {len(csv_files)} CSV file(s).\n")

    records = []

    for fpath in csv_files:
        transect_id = parse_transect_id(fpath)

        try:
            df_raw      = load_timeseries(fpath)
            df_filtered = filter_dates(df_raw, START_DATE, END_DATE)

            if len(df_filtered) < MIN_OBS:
                print(f"  SKIP  {transect_id}  ({len(df_filtered)} obs < MIN_OBS={MIN_OBS})")
                result = _empty_lrr(len(df_filtered))
            else:
                result = compute_lrr(df_filtered)
                print(f"  OK    {transect_id}  "
                      f"n={result['n_obs']:4d}  "
                      f"LRR={result['lrr_m_yr']:+7.3f} m/yr  "
                      f"R²={result['r_squared']:.3f}")

        except Exception as e:
            print(f"  ERROR {transect_id}: {e}")
            result = _empty_lrr(0)

        records.append({"transect_id": transect_id, "filepath": fpath, **result})

    results_df = pd.DataFrame(records)

    # ---- Summary stats ----
    valid = results_df.dropna(subset=["lrr_m_yr"])
    print(f"\n{'='*55}")
    print(f"  Transects processed : {len(results_df)}")
    print(f"  With valid LRR      : {len(valid)}")
    if len(valid):
        print(f"  Mean LRR            : {valid['lrr_m_yr'].mean():+.3f} m/yr")
        print(f"  Median LRR          : {valid['lrr_m_yr'].median():+.3f} m/yr")
        print(f"  Min / Max LRR       : {valid['lrr_m_yr'].min():+.3f} / {valid['lrr_m_yr'].max():+.3f} m/yr")
    print(f"{'='*55}\n")

    # ---- Save results ----
    if OUTPUT_CSV:
        results_df.to_csv(OUTPUT_CSV, index=False)
        print(f"Results saved to: {OUTPUT_CSV}")

    # ---- Plots ----
    if len(valid) > 1:
        plot_lrr_summary(results_df)

    return results_df


# ============================================================
# SINGLE-TRANSECT INTERACTIVE MODE
# (call this from a notebook or REPL to inspect one transect)
# ============================================================

def inspect_transect(filepath: str,
                     start: str | None = START_DATE,
                     end: str | None   = END_DATE):
    """
    Load, filter, compute LRR, and plot a single transect CSV.

    Example:
        from coastsat_lrr_analysis import inspect_transect
        inspect_transect("data/usa_NC_0033_0002.csv", "1997-01-01", "2019-12-31")
    """
    tid        = parse_transect_id(filepath)
    df_raw     = load_timeseries(filepath)
    df_filt    = filter_dates(df_raw, start, end)
    result     = compute_lrr(df_filt)

    print(f"\nTransect : {tid}")
    print(f"Period   : {start} → {end}")
    print(f"N obs    : {result['n_obs']}")
    print(f"LRR      : {result['lrr_m_yr']:+.4f} m/yr")
    print(f"95% CI   : ±{result['unc_m_yr']:.4f} m/yr")
    print(f"R²       : {result['r_squared']:.4f}")
    print(f"p-value  : {result['p_value']:.6f}")

    plot_timeseries(df_raw, df_filt, result, tid)
    return result


# ============================================================
if __name__ == "__main__":
    results = main()
