"""
Step 2: shoreline change rate (OLS slope, LRR) per transect and per zone.

    python shoreline_rates_template.py --start-year 1996 --end-year 2024
    python shoreline_rates_template.py --from 2003-09-19 --to 2024-12-31

Window = 1 Jan start year through 31 Dec end year. --from/--to instead give
the first and last day, both included whole. Input: one CSV per
transect (col 1 date, col 2 cross-shore position in m), filename = transect
id; the lookup from step 1. Needs pandas, numpy, scipy, matplotlib.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# --- CONFIG ------------------------------------------------------------------
TIMESERIES_DIR = Path("data/timeseries")
LOOKUP_CSV     = Path("data/transect_zones.csv")   # optional
OUTPUT_DIR     = Path("output")
START_YEAR, END_YEAR = 1996, 2024

MIN_OBS      = 10
MAX_ABS_RATE = 50.0      # m/yr
MAX_P_VALUE  = 1.0       # off; screening on p drops stable transects (README)
MIN_R2       = 0.0       # off, same reason
SEAWARD_POSITIVE = True
# -----------------------------------------------------------------------------


# Read one transect CSV: date + position, by column order
def load_timeseries(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(columns={df.columns[0]: "date", df.columns[1]: "position_m"})
    df["date"] = pd.to_datetime(df["date"], utc=True, errors="coerce")
    df["position_m"] = pd.to_numeric(df["position_m"], errors="coerce")
    return df.dropna(subset=["date", "position_m"]).sort_values("date")


# transect_id -> zone_id from step 1
def load_lookup(path: Path) -> dict[str, str]:
    if not path.exists():
        print(f"  no lookup at {path}; one zone 'all'")
        return {}
    t = pd.read_csv(path, dtype=str)
    return dict(zip(t.iloc[:, 0].str.strip(), t.iloc[:, 1].str.strip()))


# OLS slope of position vs time, with 95% CI
def compute_lrr(df: pd.DataFrame) -> dict:
    n = len(df)
    if n < 3:
        return {"rate_m_yr": np.nan, "r_squared": np.nan, "p_value": np.nan,
                "unc_m_yr": np.nan, "n_obs": n}
    years = (df["date"] - df["date"].min()).dt.total_seconds() / (86400 * 365.25)
    fit = stats.linregress(years, df["position_m"])
    sign = 1 if SEAWARD_POSITIVE else -1
    return {"rate_m_yr": sign * fit.slope, "r_squared": fit.rvalue ** 2,
            "p_value": fit.pvalue,
            "unc_m_yr": stats.t.ppf(0.975, n - 2) * fit.stderr,   # 95% CI half-width
            "n_obs": n,
            "first_date": df["date"].min().date().isoformat(),
            "last_date": df["date"].max().date().isoformat()}


# Screen out bad measurements, not small signals
def keep(r: pd.Series) -> bool:
    return (np.isfinite(r["rate_m_yr"]) and r["n_obs"] >= MIN_OBS
            and r["p_value"] <= MAX_P_VALUE and r["r_squared"] >= MIN_R2
            and abs(r["rate_m_yr"]) <= MAX_ABS_RATE)


# Window from the command line: whole calendar years, or first and last day.
# Returns (lo, hi, tag), lo inclusive, hi exclusive (the day after the last day)
def choose_window(a: argparse.Namespace):
    if bool(a.date_from) != bool(a.date_to):
        raise SystemExit("give both --from and --to, or neither")
    if a.date_from:
        if a.start_year is not None or a.end_year is not None:
            raise SystemExit("give --from/--to or --start-year/--end-year, not both")
        lo = pd.Timestamp(a.date_from, tz="UTC")
        last = pd.Timestamp(a.date_to, tz="UTC")
        if last < lo:
            raise SystemExit("--to is before --from")
        return lo, last + pd.Timedelta(days=1), f"{lo.date()}_{last.date()}"
    y0 = START_YEAR if a.start_year is None else a.start_year
    y1 = END_YEAR if a.end_year is None else a.end_year
    if y1 < y0:
        raise SystemExit("--end-year is before --start-year")
    return (pd.Timestamp(f"{y0}-01-01", tz="UTC"),
            pd.Timestamp(f"{y1 + 1}-01-01", tz="UTC"), f"{y0}_{y1}")


# Bar chart of zone rates, blue seaward
def draw_bars(per_zone: pd.DataFrame, tag: str, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar(per_zone["zone_id"].astype(str), per_zone["mean_rate_m_yr"],
           yerr=per_zone["mean_unc_m_yr"], edgecolor="none",
           color=["#2166ac" if v > 0 else "#b2182b" for v in per_zone["mean_rate_m_yr"]],
           error_kw={"ecolor": "0.4", "lw": 0.8})
    ax.axhline(0, color="0.2", lw=0.8)
    ax.set(xlabel="zone", ylabel="shoreline change (m/yr)",
           title=f"Shoreline change rate, {tag.replace('_', ' to ')} (blue seaward)")
    if len(per_zone) > 25:
        ax.set_xticks([])
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


# Run: window, fit, screen, average by zone, write
def main() -> None:
    ap = argparse.ArgumentParser(description="Shoreline change rate per zone.")
    ap.add_argument("--start-year", type=int)
    ap.add_argument("--end-year", type=int)
    ap.add_argument("--from", dest="date_from", metavar="YYYY-MM-DD")
    ap.add_argument("--to", dest="date_to", metavar="YYYY-MM-DD")
    ap.add_argument("--timeseries", type=Path, default=TIMESERIES_DIR)
    ap.add_argument("--out", type=Path, default=OUTPUT_DIR)
    a = ap.parse_args()
    lo, hi, tag = choose_window(a)

    # Find the transect files and the lookup
    files = sorted(a.timeseries.rglob("*.csv"))
    if not files:
        raise SystemExit(f"no CSVs in {a.timeseries.resolve()}")
    lookup = load_lookup(LOOKUP_CSV)

    # Fit every transect inside the window
    rows = []
    for f in files:
        df = load_timeseries(f)
        df = df[(df["date"] >= lo) & (df["date"] < hi)]
        rows.append({"transect_id": f.stem,
                     "zone_id": lookup.get(f.stem) if lookup else "all",
                     **compute_lrr(df)})
    per_transect = pd.DataFrame(rows)
    per_transect["kept"] = per_transect.apply(keep, axis=1)
    print(f"{len(files)} transects, {lo.date()} to {(hi - pd.Timedelta(days=1)).date()}: "
          f"{per_transect['kept'].sum()} pass the screen, "
          f"{per_transect['zone_id'].isna().sum()} not in the lookup")

    # Average the kept transects by zone
    good = per_transect[per_transect["kept"]]
    if good.empty:
        raise SystemExit("nothing passed the screen")
    per_zone = (good.groupby("zone_id")
                    .agg(mean_rate_m_yr=("rate_m_yr", "mean"),
                         std_rate_m_yr=("rate_m_yr", "std"),
                         mean_unc_m_yr=("unc_m_yr", "mean"),
                         n_transects=("rate_m_yr", "size"))
                    .reset_index())

    # Write both tables
    a.out.mkdir(parents=True, exist_ok=True)
    per_transect.to_csv(a.out / f"rates_per_transect_{tag}.csv", index=False)
    per_zone.to_csv(a.out / f"rates_per_zone_{tag}.csv", index=False)

    # Bar chart of zone rates
    draw_bars(per_zone, tag, a.out / f"rates_per_zone_{tag}.png")

    print(per_zone.to_string(index=False, float_format=lambda v: f"{v:8.2f}"))


if __name__ == "__main__":
    main()
