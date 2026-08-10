"""
CoastSat Transect LRR (Linear Regression Rate) Analysis
========================================================
Loads CoastSat time-series CSVs from one or more folders,
filters to a user-defined date range OR a set of specific years,
and computes the Linear Regression Rate (LRR) for every transect.

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
DATA_FOLDERS = [
    r"C:/path/to/your/coastsat_data/site1",
    r"C:/path/to/your/coastsat_data/site2",
]

# Date range filter (inclusive). Used when MATCH_YEARS is empty.
# Set either to None to include all dates.
START_DATE = "1984-01-01"
END_DATE   = "2004-01-01"

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
    """Extract transect ID from filename (stem, no extension)."""
    basename = os.path.splitext(os.path.basename(filepath))[0]
    return basename


def load_timeseries(filepath: str) -> pd.DataFrame:
    """
    Load a CoastSat time-series CSV.
    Returns a DataFrame with columns ['date', 'chainage_m'].
    """
    df = pd.read_csv(filepath, header=0)
    df.columns = [c.strip() for c in df.columns]
    date_col     = df.columns[0]
    chainage_col = df.columns[1]

    df = df.rename(columns={date_col: "date", chainage_col: "chainage_m"})
    df["date"] = pd.to_datetime(df["date"], utc=True)
    df["chainage_m"] = pd.to_numeric(df["chainage_m"], errors="coerce")
    df = df.dropna(subset=["date", "chainage_m"]).sort_values("date").reset_index(drop=True)
    return df


def filter_dates(df: pd.DataFrame, start: str | None, end: str | None) -> pd.DataFrame:
    """
    Clip DataFrame to a continuous [start, end] date range (inclusive).

    Use this when you want all CoastSat observations between two dates,
    e.g. all imagery from 1997-01-01 to 2019-12-31.

    For matching specific USGS shoreline years instead, use filter_to_years().
    """
    if start:
        df = df[df["date"] >= pd.Timestamp(start, tz="UTC")]
    if end:
        df = df[df["date"] <= pd.Timestamp(end, tz="UTC")]
    return df.reset_index(drop=True)


def filter_to_dates(df: pd.DataFrame,
                    survey_dates: list[str],
                    window_days: int = 30) -> pd.DataFrame:
    """
    Retain only CoastSat observations that fall within ±window_days of
    one or more specific USGS shoreline survey dates.

    This is the preferred filter for direct DSAS comparisons because it
    anchors the CoastSat window to the actual survey date of each
    digitized shoreline, not just the calendar year.  For example, if
    your DSAS Period 1997–2019 uses shorelines digitized on
    1997-09-15, 2008-04-22, and 2019-10-03, this function will collect
    CoastSat imagery within ±30 days of each of those dates, ensuring
    seasonal and tidal conditions are as comparable as possible.

    Parameters
    ----------
    df : pd.DataFrame
        Output of load_timeseries() — must have a tz-aware 'date' column.
    survey_dates : list[str]
        ISO-8601 date strings for each USGS shoreline used in DSAS,
        e.g. ["1997-09-15", "2008-04-22", "2019-10-03"].
    window_days : int, optional
        Half-width of the search window in days around each survey date.
        Default 30 (i.e. ±1 month).  Increase if CoastSat has sparse
        coverage in your area (cloud cover, satellite gaps).

    Returns
    -------
    pd.DataFrame
        Filtered copy sorted by date, index reset.
        Column 'survey_date' is added to show which anchor date each
        observation was matched to (nearest anchor wins if windows overlap).

    Notes
    -----
    If the windows for two survey dates overlap, an observation is
    attributed to the nearest anchor date.  This avoids double-counting
    a single image in the regression.

    Examples
    --------
    # Period 1997–2019 with 3 USGS shorelines, ±30-day window
    df_matched = filter_to_dates(
        df_raw,
        survey_dates = ["1997-09-15", "2008-04-22", "2019-10-03"],
        window_days  = 30,
    )

    # Period 1978–1997 with 3 USGS shorelines, ±30-day window
    df_matched = filter_to_dates(
        df_raw,
        survey_dates = ["1978-06-10", "1986-08-01", "1997-09-15"],
        window_days  = 30,
    )
    """
    if not survey_dates:
        return df.copy()

    anchors = [pd.Timestamp(d, tz="UTC") for d in survey_dates]

    # For each observation find the nearest anchor and its distance in days
    def nearest_anchor(obs_date):
        diffs = [(abs((obs_date - a).total_seconds()) / 86400.0, a) for a in anchors]
        return min(diffs, key=lambda x: x[0])

    rows = []
    for _, row in df.iterrows():
        dist_days, anchor = nearest_anchor(row["date"])
        if dist_days <= window_days:
            rows.append({**row.to_dict(), "survey_date": anchor.strftime("%Y-%m-%d")})

    if not rows:
        empty = df.iloc[:0].copy()
        empty["survey_date"] = pd.Series(dtype="str")
        return empty

    out = pd.DataFrame(rows)
    return out.sort_values("date").reset_index(drop=True)


# kept for backwards compatibility — wraps filter_to_dates using Jan 1 of each year
def filter_to_years(df: pd.DataFrame,
                    years: list[int],
                    window_days: int = 0) -> pd.DataFrame:
    """
    Backwards-compatible wrapper.  Prefer filter_to_dates() for new runs.
    If window_days > 0, anchors windows on Jan 1 of each year.
    If window_days == 0, keeps any observation whose calendar year is in the list.
    """
    if not years:
        return df.copy()
    if window_days > 0:
        survey_dates = [f"{yr}-01-01" for yr in years]
        return filter_to_dates(df, survey_dates, window_days)
    else:
        filtered = df[df["date"].dt.year.isin(years)].copy()
        return filtered.sort_values("date").reset_index(drop=True)


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

    t0   = df["date"].min()
    days = (df["date"] - t0).dt.total_seconds() / 86400.0
    yrs  = days / 365.25

    x = yrs.values
    y = df["chainage_m"].values

    slope, intercept, r, p, se = stats.linregress(x, y)

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

    ax.scatter(df_raw["date"], df_raw["chainage_m"],
               color="lightgray", s=15, zorder=1, label="All obs.")
    ax.scatter(df_filtered["date"], df_filtered["chainage_m"],
               color="steelblue", s=20, zorder=2, label="Selected obs.")

    if not np.isnan(lrr_result["lrr_m_yr"]) and len(df_filtered) >= 2:
        t0   = df_filtered["date"].min()
        days = (df_filtered["date"] - t0).dt.total_seconds() / 86400.0
        yrs  = days / 365.25
        coeffs = np.polyfit(yrs.values, df_filtered["chainage_m"].values, 1)
        y_hat  = np.polyval(coeffs, yrs)
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

    valid = results_df.dropna(subset=["lrr_m_yr"])
    print(f"\n{'='*55}")
    print(f"  Transects processed : {len(results_df)}")
    print(f"  With valid LRR      : {len(valid)}")
    if len(valid):
        print(f"  Mean LRR            : {valid['lrr_m_yr'].mean():+.3f} m/yr")
        print(f"  Median LRR          : {valid['lrr_m_yr'].median():+.3f} m/yr")
        print(f"  Min / Max LRR       : {valid['lrr_m_yr'].min():+.3f} / {valid['lrr_m_yr'].max():+.3f} m/yr")
    print(f"{'='*55}\n")

    if OUTPUT_CSV:
        results_df.to_csv(OUTPUT_CSV, index=False)
        print(f"Results saved to: {OUTPUT_CSV}")

    if len(valid) > 1:
        plot_lrr_summary(results_df)

    return results_df


# ============================================================
# SINGLE-TRANSECT INTERACTIVE MODE
# ============================================================

def inspect_transect(filepath: str,
                     start: str | None = START_DATE,
                     end: str | None   = END_DATE,
                     match_years: list[int] | None = None,
                     window_days: int = 0):
    """
    Load, filter, compute LRR, and plot a single transect CSV.

    Pass match_years to use year-matching instead of a continuous range.

    Examples:
        # Full date range
        inspect_transect("data/usa_NC_0033_0002.csv", "1997-01-01", "2019-12-31")

        # Year-matched (DSAS-comparable)
        inspect_transect("data/usa_NC_0033_0002.csv", match_years=[1997, 2019])
    """
    tid     = parse_transect_id(filepath)
    df_raw  = load_timeseries(filepath)

    if match_years:
        df_filt = filter_to_years(df_raw, match_years, window_days)
        period  = f"years={match_years}"
    else:
        df_filt = filter_dates(df_raw, start, end)
        period  = f"{start} → {end}"

    result = compute_lrr(df_filt)

    print(f"\nTransect : {tid}")
    print(f"Period   : {period}")
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
