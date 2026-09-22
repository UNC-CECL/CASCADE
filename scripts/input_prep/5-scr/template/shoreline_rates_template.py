"""
shoreline_rates_template.py
==============================================================================
A shoreline-change rate, from satellite time series to a table and a figure,
in one standalone file.

This is a TEMPLATE, not a product. It is the process the 5-scr folder runs,
stripped of everything specific to Hatteras Island: no repository layout, no
site config, no house figure style, no path resolvers. Copy it, point the
CONFIG block at your own data, and edit from there.

It imports nothing from this repository, so it runs anywhere pandas, numpy,
scipy and matplotlib are installed.

THE PROCESS, which is the part worth keeping
--------------------------------------------
    1  LOAD      one time series per transect: a date and a cross-shore
                 position. That is all a rate needs.
    2  WINDOW    cut to the period you are fitting. A rate is meaningless
                 without the interval it was fitted over, so the interval
                 travels with the number from here on.
    3  FIT       ordinary least squares through position against time, per
                 transect. The slope is the rate.
    4  SCREEN    drop BAD MEASUREMENTS, before averaging. A transect with
                 four observations and a 40 m/yr slope will otherwise set the
                 tone for its whole zone. Do not drop small signals: see the
                 note on significance screening in CONFIG, which is the one
                 piece of statistics in this file that will bite you.
    5  GROUP     average the surviving transects into the zones you actually
                 reason about (model cells, management reaches, beaches).
    6  WRITE     one table per zone and one per transect, then a figure.
                 Keep the transect table: it is the only way to tell a real
                 signal from one noisy transect.

WHAT TO CHANGE, and what not to
-------------------------------
Change the CONFIG block. Change `load_timeseries` if your CSVs are shaped
differently. Everything else is the method, and changing it changes what the
number means -- so if you do, say so in writing beside the output.

The one thing worth keeping exactly: **screen before you average** (step 4).
It is the difference between a rate and an average of noise.

USAGE
    python shoreline_rates_template.py
    python shoreline_rates_template.py --start 2000-01-01 --end 2020-01-01

INPUT this expects
    TIMESERIES_DIR/
        <anything>_0001.csv      one CSV per transect, the id in the filename
        <anything>_0002.csv      col 1 = date, col 2 = cross-shore position (m)
        ...
    LOOKUP_CSV                   two columns: transect_id, zone_id
                                 Omit it and every transect lands in one zone.
==============================================================================
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# =============================================================================
# CONFIG  -- the only part you should need to edit
# =============================================================================

TIMESERIES_DIR = Path("data/timeseries")   # folder of per-transect CSVs
LOOKUP_CSV     = Path("data/transect_zones.csv")   # transect_id,zone_id (optional)
OUTPUT_DIR     = Path("output")

START_DATE = "1996-01-01"     # inclusive; None for "from the beginning"
END_DATE   = "2024-01-01"     # exclusive; None for "to the end"

# Step 4, the screen. Judgement calls, not universal constants.
MIN_OBS      = 10      # a slope through fewer points than this is noise
MAX_ABS_RATE = 50.0    # m/yr; anything larger is a georeferencing error,
                       # not a shoreline

# SIGNIFICANCE SCREENING IS OFF BY DEFAULT, AND THAT IS DELIBERATE.
# Setting MAX_P_VALUE = 0.05 looks like good practice and quietly biases your
# result. A transect whose shoreline is genuinely STABLE has a true slope near
# zero, so its trend cannot be distinguished from flat, so it fails the test
# and is dropped -- while its eroding neighbours pass. The zone mean is then
# taken over the survivors, which are the transects that moved. You get a
# number biased away from zero, and nothing in the output says so.
#
# This is not hypothetical: run this file against the synthetic set in
# README.md with MAX_P_VALUE = 0.05 and the one transect it discards is the
# one whose true rate is exactly 0.0 m/yr.
#
# Screen on things that indicate a BAD MEASUREMENT (too few observations, a
# physically impossible rate). Do not screen on the size of the signal you
# are trying to measure.
MAX_P_VALUE  = 1.0     # 1.0 = off. See above before changing.
MIN_R2       = 0.0     # 0.0 = off. Same argument: a stable transect has a
                       # low R2 precisely because there is little to fit.

SEAWARD_POSITIVE = True   # False if your chainage grows LANDWARD, which
                          # silently flips every sign in the output

# =============================================================================
# 1  LOAD
# =============================================================================

def transect_id_from_name(path: Path) -> str:
    """The trailing number in the filename, or the whole stem if there is none."""
    m = re.search(r"(\d+)(?=\D*$)", path.stem)
    return m.group(1) if m else path.stem


def load_timeseries(path: Path) -> pd.DataFrame:
    """One transect's series -> DataFrame['date', 'position_m'].

    Reads by POSITION, not by column name, because every provider names these
    two columns differently. If yours has a header row that is not a header,
    or extra columns in between, this is the function to edit.
    """
    df = pd.read_csv(path, header=0)
    df = df.rename(columns={df.columns[0]: "date", df.columns[1]: "position_m"})
    df["date"] = pd.to_datetime(df["date"], utc=True, errors="coerce")
    df["position_m"] = pd.to_numeric(df["position_m"], errors="coerce")
    df = df.dropna(subset=["date", "position_m"])
    return df.sort_values("date").reset_index(drop=True)


# =============================================================================
# 2  WINDOW
# =============================================================================

def filter_dates(df: pd.DataFrame, start: str | None, end: str | None) -> pd.DataFrame:
    if start:
        df = df[df["date"] >= pd.Timestamp(start, tz="UTC")]
    if end:
        df = df[df["date"] < pd.Timestamp(end, tz="UTC")]
    return df


# =============================================================================
# 3  FIT
# =============================================================================

def compute_lrr(df: pd.DataFrame) -> dict:
    """Ordinary least squares: position against time, in m/yr.

    Time is measured in years from the FIRST observation in the window, not
    from year zero, so the intercept stays a sensible number.

    `unc_m_yr` is the half-width of the 95% confidence interval on the slope.
    Report it. A rate of -1.2 +/- 0.3 m/yr and one of -1.2 +/- 4.0 m/yr are
    not the same finding, and only the uncertainty tells them apart.
    """
    n = len(df)
    empty = {"rate_m_yr": np.nan, "r_squared": np.nan, "p_value": np.nan,
             "unc_m_yr": np.nan, "n_obs": n, "start_date": None, "end_date": None}
    if n < 2:
        return empty

    years = (df["date"] - df["date"].min()).dt.total_seconds() / (86400.0 * 365.25)
    slope, _intercept, r, p, stderr = stats.linregress(years.values,
                                                       df["position_m"].values)

    dof = n - 2
    t95 = stats.t.ppf(0.975, dof) if dof > 0 else np.nan

    if not SEAWARD_POSITIVE:
        slope = -slope

    return {"rate_m_yr": slope, "r_squared": r ** 2, "p_value": p,
            "unc_m_yr": t95 * stderr, "n_obs": n,
            "start_date": df["date"].min().date().isoformat(),
            "end_date": df["date"].max().date().isoformat()}


# =============================================================================
# 4  SCREEN  -- before averaging, never after
# =============================================================================

def keep(row: pd.Series) -> bool:
    if not np.isfinite(row["rate_m_yr"]):
        return False
    return (row["n_obs"] >= MIN_OBS
            and row["p_value"] <= MAX_P_VALUE
            and row["r_squared"] >= MIN_R2
            and abs(row["rate_m_yr"]) <= MAX_ABS_RATE)


# =============================================================================
# 5  GROUP
# =============================================================================

def load_lookup(path: Path) -> dict[str, str]:
    """transect_id -> zone_id. Missing file means one zone called 'all'."""
    if not path.exists():
        print(f"  no lookup at {path}; every transect goes in one zone")
        return {}
    t = pd.read_csv(path, dtype=str)
    return dict(zip(t.iloc[:, 0].str.strip(), t.iloc[:, 1].str.strip()))


# =============================================================================
# the run
# =============================================================================

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[2])
    ap.add_argument("--start", default=START_DATE)
    ap.add_argument("--end", default=END_DATE)
    ap.add_argument("--timeseries", type=Path, default=TIMESERIES_DIR)
    ap.add_argument("--out", type=Path, default=OUTPUT_DIR)
    args = ap.parse_args()

    files = sorted(args.timeseries.glob("*.csv"))
    if not files:
        raise SystemExit(f"no CSVs in {args.timeseries.resolve()}")
    print(f"{len(files)} transect files, window {args.start} .. {args.end}")

    lookup = load_lookup(LOOKUP_CSV)

    rows = []
    for f in files:
        tid = transect_id_from_name(f)
        series = filter_dates(load_timeseries(f), args.start, args.end)
        rows.append({"transect_id": tid, "zone_id": lookup.get(tid, "all"),
                     **compute_lrr(series)})

    per_transect = pd.DataFrame(rows)
    per_transect["kept"] = per_transect.apply(keep, axis=1)

    n_kept = int(per_transect["kept"].sum())
    print(f"  {n_kept} of {len(per_transect)} transects pass the screen")
    if n_kept == 0:
        raise SystemExit("nothing survived the screen -- loosen CONFIG, or "
                         "check that the date window overlaps your data")

    good = per_transect[per_transect["kept"]]
    per_zone = (good.groupby("zone_id")
                    .agg(mean_rate_m_yr=("rate_m_yr", "mean"),
                         std_rate_m_yr=("rate_m_yr", "std"),
                         mean_unc_m_yr=("unc_m_yr", "mean"),
                         n_transects=("rate_m_yr", "size"))
                    .reset_index())

    # 6  WRITE. Both tables, always: the per-transect one is what lets you
    # check whether a zone's mean is a signal or one loud transect.
    args.out.mkdir(parents=True, exist_ok=True)
    tag = f"{args.start[:4]}_{args.end[:4]}"
    per_transect.to_csv(args.out / f"rates_per_transect_{tag}.csv", index=False)
    per_zone.to_csv(args.out / f"rates_per_zone_{tag}.csv", index=False)

    fig, ax = plt.subplots(figsize=(10, 4))
    colours = ["#2166ac" if v > 0 else "#b2182b" for v in per_zone["mean_rate_m_yr"]]
    ax.bar(per_zone["zone_id"].astype(str), per_zone["mean_rate_m_yr"],
           yerr=per_zone["mean_unc_m_yr"], color=colours, edgecolor="none",
           error_kw={"ecolor": "0.4", "lw": 0.8})
    ax.axhline(0, color="0.2", lw=0.8)
    ax.set_ylabel("shoreline change (m/yr)")
    ax.set_xlabel("zone")
    ax.set_title(f"Shoreline change rate, {args.start[:4]}-{args.end[:4]}  "
                 f"(blue seaward, red landward)")
    if len(per_zone) > 25:
        ax.set_xticks([])
    fig.tight_layout()
    fig.savefig(args.out / f"rates_per_zone_{tag}.png", dpi=200)
    plt.close(fig)

    print(f"  wrote 2 tables and 1 figure to {args.out.resolve()}")
    print(per_zone.to_string(index=False, float_format=lambda v: f"{v:8.2f}"))


if __name__ == "__main__":
    main()
