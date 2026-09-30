"""
shoreline_endpoint_template.py
==============================================================================
Net shoreline change, from satellite time series to a table and a figure,
in one standalone file: the mean position at the end of a window minus the
mean position at its start.

This is a TEMPLATE, not a product. It is the endpoint reading the 5-scr
folder makes (3-rates/coastsat/endpoint/, and the observed side of
total_change/), stripped of everything specific to Hatteras Island. It is the
companion of shoreline_rates_template.py: same input, same lookup, same
window, a different question.

It imports nothing from this repository. It needs pandas, numpy, matplotlib.

RATE OR ENDPOINT -- they are not two ways of getting the same number
--------------------------------------------------------------------
    LRR (the rates template)   a straight line through EVERY position in the
                               window. Asks: what is the trend?
    ENDPOINT (this file)       where the shoreline sat at the end minus where
                               it sat at the start. Asks: how far did it go?

They agree on a shoreline that moves steadily. They disagree when something
happened late or early in the window -- a nourishment, a storm, a jump -- and
that disagreement is information, not error. Compute both.

THE PROCESS
-----------
    1  LOAD      one time series per transect, as in the rates template.
    2  ENDS      the two end windows. By default the WHOLE start year and the
                 WHOLE end year: a window 1996-2010 compares the mean over
                 1996-01-01..1996-12-31 with the mean over 2010-01-01..
                 2010-12-31, the same calendar years the rate template fits.
                 A single satellite pass carries metres of tide, wave set-up
                 and cloud-edge noise; a year of passes averages one full
                 seasonal cycle, so seasonal width cancels.
    3  DIFFERENCE  end mean minus start mean, metres. Divided by the interval
                 between the two window centres for a rate.
    4  SCREEN    drop a transect with too few passes at EITHER end. Not on
                 the size of the change -- the same argument as the rates
                 template's note on significance screening.
    5  GROUP     average surviving transects into zones.
    6  WRITE     per-zone and per-transect tables, and a figure.

THE INTERVAL
------------
The rate divides by the calendar interval between the two window centres --
14.0 yr for 1996-2010 with one-year ends -- not by the gap between the mean
dates of the passes that happened to fall in each window. The period is the
calendar years; the pass dates are an accident of cloud cover. The pass-date
gap is written beside it (interval_obs_yr) so you can see how far apart they
are, rather than silently adjusting either.

ANCHORING ON A SURVEY INSTEAD (optional)
----------------------------------------
To difference like for like against something measured on a day -- a
digitized dune line, an aerial photo, a lidar flight -- give --start-date and
--end-date. Each end window is then +/- HALF_WINDOW_DAYS about that date, and
the rate divides by the interval between the two dates. Use this only for that
comparison; for a period's change, keep the calendar years.

USAGE
    python shoreline_endpoint_template.py --start-year 1996 --end-year 2010
    python shoreline_endpoint_template.py --start-year 1996 --end-year 2024 --end-window-years 3
    python shoreline_endpoint_template.py --start-date 1997-10-12 --end-date 2009-05-30

INPUT this expects -- identical to shoreline_rates_template.py
    TIMESERIES_DIR/   one CSV per transect, searched recursively; the filename
                      without .csv is the transect id; col 1 = date,
                      col 2 = cross-shore position (m)
    LOOKUP_CSV        transect_id, zone_id (transect_zone_join_template.py)
==============================================================================
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# =============================================================================
# CONFIG  -- the only part you should need to edit
# =============================================================================

TIMESERIES_DIR = Path("data/timeseries")
LOOKUP_CSV     = Path("data/transect_zones.csv")   # optional
OUTPUT_DIR     = Path("output")

# The window, in whole calendar years, BOTH INCLUDED -- the same convention
# as shoreline_rates_template.py, so a 1996-2010 endpoint and a 1996-2010
# rate describe the same span.
START_YEAR = 1996
END_YEAR   = 2024

# How many calendar years to average at each end. 1 = the start year alone
# and the end year alone. More years is less noise and a blurrier "when": a
# 3-year end window over 2022-2024 is centred on mid-2023, not on 2024.
END_WINDOW_YEARS = 1

# Survey anchoring (see the docstring). Half-width of each end window in days.
HALF_WINDOW_DAYS = 182.5   # +/- six months: one full seasonal cycle

# Step 4, the screen.
MIN_OBS_PER_END = 3        # passes needed in EACH end window
MAX_ABS_RATE    = 50.0     # m/yr; bigger is a georeferencing error

SEAWARD_POSITIVE = True    # False if your chainage grows LANDWARD

DAYS_PER_YEAR = 365.25

# =============================================================================
# 1  LOAD  -- the same two functions as the rates template
# =============================================================================

def transect_id_from_name(path: Path) -> str:
    """The whole filename without .csv (CoastSat's transect id)."""
    return path.stem


def load_timeseries(path: Path) -> pd.DataFrame:
    """One transect's series -> DataFrame['date', 'position_m'], read by position."""
    df = pd.read_csv(path, header=0)
    df = df.rename(columns={df.columns[0]: "date", df.columns[1]: "position_m"})
    df["date"] = pd.to_datetime(df["date"], utc=True, errors="coerce")
    df["position_m"] = pd.to_numeric(df["position_m"], errors="coerce")
    df = df.dropna(subset=["date", "position_m"])
    return df.sort_values("date").reset_index(drop=True)


def load_lookup(path: Path) -> dict[str, str]:
    """transect_id -> zone_id. Missing file means one zone called 'all'."""
    if not path.exists():
        print(f"  no lookup at {path}; every transect goes in one zone")
        return {}
    t = pd.read_csv(path, dtype=str)
    return dict(zip(t.iloc[:, 0].str.strip(), t.iloc[:, 1].str.strip()))


# =============================================================================
# 2  ENDS
# =============================================================================

def calendar_ends(start_year: int, end_year: int, n_years: int):
    """Two (lo, hi, centre) windows of whole calendar years, lo inclusive,
    hi exclusive. The centre is what the rate's interval is measured between."""
    if 2 * n_years > end_year - start_year + 1:
        raise SystemExit(f"{n_years}-year end windows overlap inside "
                         f"{start_year}-{end_year}; use fewer years")
    def win(first_year):
        lo = pd.Timestamp(f"{first_year}-01-01", tz="UTC")
        hi = pd.Timestamp(f"{first_year + n_years}-01-01", tz="UTC")
        return lo, hi, lo + (hi - lo) / 2
    return win(start_year), win(end_year - n_years + 1)


def survey_ends(start_date: str, end_date: str, half_days: float):
    """Two (lo, hi, centre) windows of +/- half_days about survey dates."""
    half = pd.Timedelta(days=half_days)
    def win(d):
        c = pd.Timestamp(d, tz="UTC")
        return c - half, c + half, c
    a, b = win(start_date), win(end_date)
    if a[1] > b[0]:
        raise SystemExit("the two survey windows overlap; shorten HALF_WINDOW_DAYS")
    return a, b


# =============================================================================
# 3  DIFFERENCE
# =============================================================================

def endpoint(df: pd.DataFrame, ends, interval_yr: float) -> dict:
    out = {}
    for tag, (lo, hi, _c) in zip(("start", "end"), ends):
        sel = df[(df["date"] >= lo) & (df["date"] < hi)]
        out[f"position_{tag}_m"] = sel["position_m"].mean() if len(sel) else np.nan
        out[f"n_{tag}"] = len(sel)
        out[f"first_{tag}"] = sel["date"].min().date().isoformat() if len(sel) else None
        out[f"last_{tag}"] = sel["date"].max().date().isoformat() if len(sel) else None
        out[f"_mean_date_{tag}"] = sel["date"].mean() if len(sel) else pd.NaT

    change = out["position_end_m"] - out["position_start_m"]
    if not SEAWARD_POSITIVE:
        change = -change
    out["change_m"] = change
    out["rate_m_yr"] = change / interval_yr
    gap = out.pop("_mean_date_end") - out.pop("_mean_date_start")
    out["interval_obs_yr"] = (gap.total_seconds() / 86400.0 / DAYS_PER_YEAR
                              if pd.notna(gap) else np.nan)
    return out


# =============================================================================
# 4  SCREEN
# =============================================================================

def keep(row: pd.Series) -> bool:
    return (np.isfinite(row["change_m"])
            and row["n_start"] >= MIN_OBS_PER_END
            and row["n_end"] >= MIN_OBS_PER_END
            and abs(row["rate_m_yr"]) <= MAX_ABS_RATE)


# =============================================================================
# the run
# =============================================================================

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[2])
    ap.add_argument("--start-year", type=int, default=START_YEAR)
    ap.add_argument("--end-year", type=int, default=END_YEAR)
    ap.add_argument("--end-window-years", type=int, default=END_WINDOW_YEARS)
    ap.add_argument("--start-date", help="survey anchor (with --end-date)")
    ap.add_argument("--end-date", help="survey anchor (with --start-date)")
    ap.add_argument("--timeseries", type=Path, default=TIMESERIES_DIR)
    ap.add_argument("--out", type=Path, default=OUTPUT_DIR)
    args = ap.parse_args()

    if bool(args.start_date) != bool(args.end_date):
        raise SystemExit("give both --start-date and --end-date, or neither")
    if args.start_date:
        ends = survey_ends(args.start_date, args.end_date, HALF_WINDOW_DAYS)
        tag = f"{args.start_date}_{args.end_date}"
        how = f"+/-{HALF_WINDOW_DAYS:g} d about the survey dates"
    else:
        if args.end_year <= args.start_year:
            raise SystemExit("--end-year must be after --start-year")
        ends = calendar_ends(args.start_year, args.end_year, args.end_window_years)
        tag = f"{args.start_year}_{args.end_year}"
        if args.end_window_years != 1:
            tag += f"_ends{args.end_window_years}yr"
        how = (f"mean of {args.end_window_years} calendar year(s) at each end")
    interval_yr = (ends[1][2] - ends[0][2]).total_seconds() / 86400.0 / DAYS_PER_YEAR
    # Whole calendar years give a whole number of years between centres, up
    # to leap days; round so 1996-2010 reads 14.0, not 13.999.
    if not args.start_date:
        interval_yr = round(interval_yr)

    files = sorted(args.timeseries.rglob("*.csv"))
    if not files:
        raise SystemExit(f"no CSVs in {args.timeseries.resolve()}")
    for label, (lo, hi, _c) in zip(("start", "end"), ends):
        print(f"  {label} window {lo.date()} .. {(hi - pd.Timedelta(seconds=1)).date()}")
    print(f"{len(files)} transect files, {how}, interval {interval_yr:.2f} yr")

    lookup = load_lookup(LOOKUP_CSV)
    rows = []
    for f in files:
        tid = transect_id_from_name(f)
        rows.append({"transect_id": tid, "zone_id": lookup.get(tid, "all"),
                     **endpoint(load_timeseries(f), ends, interval_yr)})
    per_transect = pd.DataFrame(rows)
    per_transect["interval_yr"] = interval_yr
    per_transect["kept"] = per_transect.apply(keep, axis=1)

    n_kept = int(per_transect["kept"].sum())
    print(f"  {n_kept} of {len(per_transect)} transects pass the screen")
    if n_kept == 0:
        raise SystemExit("nothing survived the screen -- check the end windows "
                         "overlap your data, or lower MIN_OBS_PER_END")

    good = per_transect[per_transect["kept"]]
    per_zone = (good.groupby("zone_id")
                    .agg(mean_change_m=("change_m", "mean"),
                         std_change_m=("change_m", "std"),
                         mean_rate_m_yr=("rate_m_yr", "mean"),
                         pct_landward=("change_m", lambda s: 100.0 * (s < 0).mean()),
                         median_n_start=("n_start", "median"),
                         median_n_end=("n_end", "median"),
                         n_transects=("change_m", "size"))
                    .reset_index())
    per_zone["interval_yr"] = interval_yr

    args.out.mkdir(parents=True, exist_ok=True)
    per_transect.to_csv(args.out / f"endpoint_per_transect_{tag}.csv", index=False)
    per_zone.to_csv(args.out / f"endpoint_per_zone_{tag}.csv", index=False)

    fig, ax = plt.subplots(figsize=(10, 4))
    colours = ["#2166ac" if v > 0 else "#b2182b" for v in per_zone["mean_change_m"]]
    ax.bar(per_zone["zone_id"].astype(str), per_zone["mean_change_m"],
           yerr=per_zone["std_change_m"], color=colours, edgecolor="none",
           error_kw={"ecolor": "0.4", "lw": 0.8})
    ax.axhline(0, color="0.2", lw=0.8)
    ax.set_ylabel("net shoreline change (m)")
    ax.set_xlabel("zone")
    ax.set_title(f"Net shoreline change, {tag.replace('_', ' to ')}  "
                 f"(blue seaward, red landward; bars = std across transects)")
    if len(per_zone) > 25:
        ax.set_xticks([])
    fig.tight_layout()
    fig.savefig(args.out / f"endpoint_per_zone_{tag}.png", dpi=200)
    plt.close(fig)

    print(f"  wrote 2 tables and 1 figure to {args.out.resolve()}")
    print(per_zone[["zone_id", "mean_change_m", "mean_rate_m_yr", "n_transects"]]
          .to_string(index=False, float_format=lambda v: f"{v:8.2f}"))


if __name__ == "__main__":
    main()
