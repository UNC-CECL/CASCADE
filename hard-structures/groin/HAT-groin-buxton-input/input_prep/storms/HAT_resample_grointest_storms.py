"""
Storm file for the 1967-1997 groin test: real 1984-2004 storms resampled into a 30-year window.

    python HAT_resample_grointest_storms.py

Reads the 1984-2004 storm series (hat_env_forcings); writes SAVE_NAME as .npy
(CASCADE input) and .csv (inspection) into SAVE_DIR. Storm timing is not
historical. Needs numpy and pandas. Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""
from pathlib import Path

import os
import numpy as np
import pandas as pd

# --- CONFIG ------------------------------------------------------------------
_PATH_REPO = next(_p for _p in Path(__file__).resolve().parents
                  if (_p / "pyproject.toml").exists())

# Source storm series to resample: the validated 1984-2004 file
import sys as _envsys
from pathlib import Path as _EnvP
_envsys.path.insert(0, str(next(_q for _q in _EnvP(__file__).resolve().parents
                                if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_env_forcings as _env  # noqa: E402
SOURCE_NPY = str(_env.storm_series_file(1984, 2004))
TARGET_YEARS = 30   # 1967->1997; CASCADE indexes time 1..N

RESAMPLE_MODE = "bootstrap_years"   # or "bootstrap_events" (README)

EVENTS_PER_YEAR = None   # "bootstrap_events" only: mean storms/yr; None = source mean

RANDOM_SEED = 1967   # reproducible

# Duration cap and floor, truncated as the creator does; None leaves durations alone
MAX_STORM_DUR = 72
MIN_STORM_DUR = 8

SAVE_DIR  = str(_PATH_REPO / "hard-structures" / "groin" / "HAT-groin-buxton-input" / "groin_init" / "storms" / "1967_1997")
SAVE_NAME = "1967_1997_grointest_storms"
# -----------------------------------------------------------------------------


# Load the source storm array, checking it is (N, 5)
def _load_source(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"source storm file not found: {path}")
    arr = np.load(path, allow_pickle=True).astype(float)
    if arr.ndim != 2 or arr.shape[1] != 5:
        raise ValueError(f"expected (N,5) storm array, got {arr.shape}")
    return arr


# Source storms split by year: {year: rows of Rhigh, Rlow, period, duration}
def _by_year(arr):
    years = arr[:, 0].astype(int)
    out = {}
    for y in np.unique(years):
        out[y] = arr[years == y, 1:]   # drop the time column; we re-stamp it
    return out


# Resample the source into TARGET_YEARS model years, then cap and floor durations
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

    # Duration cap and floor, truncating as the creator does
    if MAX_STORM_DUR is not None:
        df["duration"] = df["duration"].clip(upper=MAX_STORM_DUR)
    if MIN_STORM_DUR is not None:
        df = df[df["duration"] >= MIN_STORM_DUR].reset_index(drop=True)

    df = df.sort_values("time").reset_index(drop=True)
    return df, src


# Run: resample, save .npy and .csv, report against the source
def main():
    df, src = build_storms()
    os.makedirs(SAVE_DIR, exist_ok=True)

    npy_path = os.path.join(SAVE_DIR, f"{SAVE_NAME}.npy")
    csv_path = os.path.join(SAVE_DIR, f"{SAVE_NAME}.csv")
    np.save(npy_path, df.to_numpy())
    df.to_csv(csv_path, index=False)

    # Report: does the resampled series resemble the source?
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
