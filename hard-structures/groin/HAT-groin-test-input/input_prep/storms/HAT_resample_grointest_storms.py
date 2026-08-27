"""
HAT_resample_grointest_storms.py
============================
Build the storm file for the 1967-1997 groin TEST run by RESAMPLING an existing
CASCADE storm series into a new N-year window. No new storm physics -- this is a
reformat/retime of real events so the record covers the run length.

Why resample instead of create new storms
------------------------------------------
The groin acts through BRIE's wave-driven alongshore transport, which is set by
the WAVE CLIMATE parameters (wave_height, period, asymmetry), NOT by the storm
file. The storm file drives Barrier3D overwash/dune response (cross-shore), which
is very nearly orthogonal to the groin test. So we do not need a 1967-1997 storm
climate -- we need a record of the right LENGTH with realistic event statistics.
Resampling the real 1984-2004 events preserves their Rhigh/Rlow/period/duration
distributions and per-year storm counts while retiming them into N years.

Storm timing is therefore NOT historical for 1967-1997. State this in the run
docstring; it is correct for a function test, not a forced hindcast.

Format (matches historical_storm_creation_v3_HAT.py comparison)
-----------------------------------------------------------
Columns: time, Rhigh, Rlow, period, duration
  time     : 1-based model year (1 .. N)
  Rhigh    : max run-up  [decameters MHW]  (NAVD88 - MHW)/10, as in the creator
  Rlow     : min run-up  [decameters MHW]
  period   : wave period at peak TWL [s]
  duration : storm duration [hrs]
Saved as both .npy (CASCADE input) and .csv (inspection).
"""

import os
import numpy as np
import pandas as pd

# ============================== CONFIG ==============================

# Source storm series to resample (your existing, validated 1984-2004 file).
SOURCE_NPY = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\storms\hindcast_storms\1984_2004\1984_2004_storms_v3_72.npy"
# Target run length. 1967->1997 = 30 model years; CASCADE indexes time 1..N.
TARGET_YEARS = 30

# How to build each target year's storms:
#   "bootstrap_years" : for each target year, copy ALL storms from a randomly
#                       chosen source year (preserves within-year clustering and
#                       realistic annual counts). Recommended.
#   "bootstrap_events": draw individual storms i.i.d. to hit a target annual count
#                       (breaks within-year correlation; use only if you want a
#                       smoother, less clustered series).
RESAMPLE_MODE = "bootstrap_years"

# For "bootstrap_events" only: mean storms/year (default = source mean).
EVENTS_PER_YEAR = None

RANDOM_SEED = 1967   # reproducible

# Duration cap (match your creator's max_storm_dur convention: TRUNCATE, do not
# drop). Set None to leave source durations untouched.
MAX_STORM_DUR = 72
MIN_STORM_DUR = 8

SAVE_DIR  = r"/scripts/groin/HAT-hindcast-groin-test/groin_init/storms/1967_1997"
SAVE_NAME = "1967_1997_grointest_storms"

# ============================ BUILD ============================

def _load_source(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"source storm file not found: {path}")
    arr = np.load(path, allow_pickle=True).astype(float)
    if arr.ndim != 2 or arr.shape[1] != 5:
        raise ValueError(f"expected (N,5) storm array, got {arr.shape}")
    return arr


def _by_year(arr):
    """Split source storms into {source_year: rows[:,1:] (Rhigh,Rlow,period,dur)}."""
    years = arr[:, 0].astype(int)
    out = {}
    for y in np.unique(years):
        out[y] = arr[years == y, 1:]   # drop the time column; we re-stamp it
    return out


def build_storms():
    rng = np.random.default_rng(RANDOM_SEED)
    src = _load_source(SOURCE_NPY)
    by_year = _by_year(src)
    src_years = sorted(by_year.keys())

    rows = []
    if RESAMPLE_MODE == "bootstrap_years":
        for target_year in range(1, TARGET_YEARS + 1):
            chosen = rng.choice(src_years)
            block = by_year[chosen]
            for r in block:
                rows.append([target_year, *r])

    elif RESAMPLE_MODE == "bootstrap_events":
        all_events = src[:, 1:]
        mean_count = (EVENTS_PER_YEAR
                      if EVENTS_PER_YEAR is not None
                      else len(src) / len(src_years))
        for target_year in range(1, TARGET_YEARS + 1):
            n = rng.poisson(mean_count)
            if n == 0:
                continue
            idx = rng.integers(0, len(all_events), size=n)
            for r in all_events[idx]:
                rows.append([target_year, *r])
    else:
        raise ValueError(f"unknown RESAMPLE_MODE: {RESAMPLE_MODE!r}")

    df = pd.DataFrame(rows, columns=["time", "Rhigh", "Rlow", "period", "duration"])

    # Duration cap/floor — truncate to match the creator's convention.
    if MAX_STORM_DUR is not None:
        df["duration"] = df["duration"].clip(upper=MAX_STORM_DUR)
    if MIN_STORM_DUR is not None:
        df = df[df["duration"] >= MIN_STORM_DUR].reset_index(drop=True)

    df = df.sort_values("time").reset_index(drop=True)
    return df, src


def main():
    df, src = build_storms()
    os.makedirs(SAVE_DIR, exist_ok=True)

    npy_path = os.path.join(SAVE_DIR, f"{SAVE_NAME}.npy")
    csv_path = os.path.join(SAVE_DIR, f"{SAVE_NAME}.csv")
    np.save(npy_path, df.to_numpy())
    df.to_csv(csv_path, index=False)

    # ── Report: does the resampled series resemble the source? ───────────────
    def stats(a_rhigh, a_dur, label, nyears, n):
        print(f"  {label:<10} storms={n:4d}  storms/yr={n/nyears:5.2f}  "
              f"Rhigh[{a_rhigh.min():.3f},{a_rhigh.max():.3f}]  "
              f"dur[{a_dur.min():.0f},{a_dur.max():.0f}]")

    src_years = len(np.unique(src[:, 0]))
    print("\nStorm resample complete.")
    print(f"  mode = {RESAMPLE_MODE}, seed = {RANDOM_SEED}")
    stats(src[:, 1], src[:, 4], "SOURCE", src_years, len(src))
    stats(df["Rhigh"].to_numpy(), df["duration"].to_numpy(),
          "RESAMPLED", TARGET_YEARS, len(df))
    print(f"\n  time range: {int(df['time'].min())}..{int(df['time'].max())} "
          f"(need 1..{TARGET_YEARS})")
    empty = set(range(1, TARGET_YEARS + 1)) - set(df["time"].astype(int))
    if empty:
        print(f"  NOTE: {len(empty)} year(s) have no storms: {sorted(empty)} "
              f"(fine -- quiet years are realistic)")
    print(f"\n  Saved: {npy_path}")
    print(f"  Saved: {csv_path}")


if __name__ == "__main__":
    main()
