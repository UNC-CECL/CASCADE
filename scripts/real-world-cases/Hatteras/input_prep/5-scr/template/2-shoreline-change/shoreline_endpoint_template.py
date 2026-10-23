"""
Step 2: net shoreline change (end mean - start mean) per transect and per zone.

    python shoreline_endpoint_template.py --start-year 1996 --end-year 2010
    python shoreline_endpoint_template.py --start-year 1996 --end-year 2024 --end-window-years 3
    python shoreline_endpoint_template.py --start-date 1997-10-12 --end-date 2009-05-30

Default ends: the whole start year and the whole end year. Rate = change /
interval between the two window centres. --start-date/--end-date centres
+/-6-month windows on survey dates instead. Input as shoreline_rates_template.py.
Needs pandas, numpy, matplotlib.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- CONFIG ------------------------------------------------------------------
TIMESERIES_DIR = Path("data/timeseries")
LOOKUP_CSV     = Path("data/transect_zones.csv")   # optional
OUTPUT_DIR     = Path("output")
START_YEAR, END_YEAR = 1996, 2024
END_WINDOW_YEARS = 1
HALF_WINDOW_DAYS = 182.5     # survey-date mode only

MIN_OBS_PER_END = 3
MAX_ABS_RATE    = 50.0       # m/yr
SEAWARD_POSITIVE = True
DAYS_PER_YEAR   = 365.25
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


# End windows as whole calendar years: (lo, hi, centre) per end, lo inclusive, hi exclusive
def calendar_ends(start_year: int, end_year: int, n: int):
    if 2 * n > end_year - start_year + 1:
        raise SystemExit(f"{n}-year end windows overlap")
    def win(y):
        lo, hi = pd.Timestamp(f"{y}-01-01", tz="UTC"), pd.Timestamp(f"{y + n}-01-01", tz="UTC")
        return lo, hi, lo + (hi - lo) / 2
    return win(start_year), win(end_year - n + 1)


# End windows centred on survey dates
def survey_ends(start_date: str, end_date: str, half_days: float):
    half = pd.Timedelta(days=half_days)
    def win(d):
        c = pd.Timestamp(d, tz="UTC")
        return c - half, c + half, c
    a, b = win(start_date), win(end_date)
    if a[1] > b[0]:
        raise SystemExit("survey windows overlap")
    return a, b


# Mean position at each end, and the change between them
def endpoint(df: pd.DataFrame, ends, interval_yr: float) -> dict:
    out, mean_dates = {}, []
    for tag, (lo, hi, _) in zip(("start", "end"), ends):
        sel = df[(df["date"] >= lo) & (df["date"] < hi)]
        out[f"position_{tag}_m"] = sel["position_m"].mean()
        out[f"n_{tag}"] = len(sel)
        mean_dates.append(sel["date"].mean() if len(sel) else pd.NaT)
    change = (out["position_end_m"] - out["position_start_m"]) * (1 if SEAWARD_POSITIVE else -1)
    gap = mean_dates[1] - mean_dates[0]
    return {**out, "change_m": change, "rate_m_yr": change / interval_yr,
            "interval_yr": interval_yr,
            "interval_obs_yr": gap.total_seconds() / 86400 / DAYS_PER_YEAR
                               if pd.notna(gap) else np.nan}   # reported, not used


# Screen out ends with too few passes
def keep(r: pd.Series) -> bool:
    return (np.isfinite(r["change_m"]) and r["n_start"] >= MIN_OBS_PER_END
            and r["n_end"] >= MIN_OBS_PER_END and abs(r["rate_m_yr"]) <= MAX_ABS_RATE)


# End windows, output tag and interval from the command line: survey dates or calendar years
def choose_ends(a: argparse.Namespace):
    if bool(a.start_date) != bool(a.end_date):
        raise SystemExit("give both --start-date and --end-date, or neither")
    if a.start_date:
        ends = survey_ends(a.start_date, a.end_date, HALF_WINDOW_DAYS)
        tag = f"{a.start_date}_{a.end_date}"
        interval_yr = (ends[1][2] - ends[0][2]).days / DAYS_PER_YEAR
        return ends, tag, interval_yr
    if a.end_year <= a.start_year:
        raise SystemExit("--end-year must be after --start-year")
    ends = calendar_ends(a.start_year, a.end_year, a.end_window_years)
    tag = f"{a.start_year}_{a.end_year}"
    if a.end_window_years != 1:
        tag += f"_ends{a.end_window_years}yr"
    interval_yr = float(round((ends[1][2] - ends[0][2]).days / DAYS_PER_YEAR))
    return ends, tag, interval_yr


# Bar chart of zone change, blue seaward
def draw_bars(per_zone: pd.DataFrame, tag: str, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar(per_zone["zone_id"].astype(str), per_zone["mean_change_m"],
           yerr=per_zone["std_change_m"], edgecolor="none",
           color=["#2166ac" if v > 0 else "#b2182b" for v in per_zone["mean_change_m"]],
           error_kw={"ecolor": "0.4", "lw": 0.8})
    ax.axhline(0, color="0.2", lw=0.8)
    ax.set(xlabel="zone", ylabel="net shoreline change (m)",
           title=f"Net shoreline change, {tag.replace('_', ' to ', 1).replace('_', ', ')} (blue seaward)")
    if len(per_zone) > 25:
        ax.set_xticks([])
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


# Run: set ends, difference, screen, average by zone, write
def main() -> None:
    ap = argparse.ArgumentParser(description="Net shoreline change per zone.")
    ap.add_argument("--start-year", type=int, default=START_YEAR)
    ap.add_argument("--end-year", type=int, default=END_YEAR)
    ap.add_argument("--end-window-years", type=int, default=END_WINDOW_YEARS)
    ap.add_argument("--start-date")
    ap.add_argument("--end-date")
    ap.add_argument("--timeseries", type=Path, default=TIMESERIES_DIR)
    ap.add_argument("--out", type=Path, default=OUTPUT_DIR)
    a = ap.parse_args()

    # Choose the end windows and the interval
    ends, tag, interval_yr = choose_ends(a)

    # Find the transect files and the lookup
    files = sorted(a.timeseries.rglob("*.csv"))
    if not files:
        raise SystemExit(f"no CSVs in {a.timeseries.resolve()}")
    lookup = load_lookup(LOOKUP_CSV)

    # Difference every transect
    rows = [{"transect_id": f.stem,
             "zone_id": lookup.get(f.stem) if lookup else "all",
             **endpoint(load_timeseries(f), ends, interval_yr)} for f in files]
    per_transect = pd.DataFrame(rows)
    per_transect["kept"] = per_transect.apply(keep, axis=1)
    for name, (lo, hi, _) in zip(("start", "end"), ends):
        print(f"  {name} window {lo.date()} to {(hi - pd.Timedelta(seconds=1)).date()}")
    print(f"{len(files)} transects, interval {interval_yr:.2f} yr: "
          f"{per_transect['kept'].sum()} pass the screen, "
          f"{per_transect['zone_id'].isna().sum()} not in the lookup")

    # Average the kept transects by zone
    good = per_transect[per_transect["kept"]]
    if good.empty:
        raise SystemExit("nothing passed the screen")
    per_zone = (good.groupby("zone_id")
                    .agg(mean_change_m=("change_m", "mean"),
                         std_change_m=("change_m", "std"),
                         mean_rate_m_yr=("rate_m_yr", "mean"),
                         pct_landward=("change_m", lambda s: 100.0 * (s < 0).mean()),
                         n_transects=("change_m", "size"))
                    .reset_index()
                    .assign(interval_yr=interval_yr))

    # Write both tables
    a.out.mkdir(parents=True, exist_ok=True)
    per_transect.to_csv(a.out / f"endpoint_per_transect_{tag}.csv", index=False)
    per_zone.to_csv(a.out / f"endpoint_per_zone_{tag}.csv", index=False)

    # Bar chart of zone change
    draw_bars(per_zone, tag, a.out / f"endpoint_per_zone_{tag}.png")

    print(per_zone.to_string(index=False, float_format=lambda v: f"{v:8.2f}"))


if __name__ == "__main__":
    main()
