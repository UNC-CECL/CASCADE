"""
Stitch the 1967-2017 storm file for the groin hindcast from three sources, on one continuous model-year index.

    python HAT_build_1967_2017_storms.py

Reads the 1967-1997 resample and the 1984-2004 and 2004-2024 real storm series
from INPUT_DIR; writes OUTPUT_PATH (1967 = year 1 .. 2017 = year 51). Needs
numpy. Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""

import os
from pathlib import Path
import numpy as np

# --- CONFIG ------------------------------------------------------------------
REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
STORMS_DIR = REPO / "hard-structures" / "groin" / "3-hindcast" / "1-dipole-1967-2017" / "inputs" / "groin_init" / "storms"
INPUT_DIR = str(STORMS_DIR / "input_storms")  # folder containing the three source .npy files
OUTPUT_PATH = str(STORMS_DIR / "1967_2017" / "1967_2017_groin_storms.npy")

RESAMPLED_FILE = os.path.join(INPUT_DIR, "1967_1997_grointest_storms.npy")
PERIOD1_FILE   = os.path.join(INPUT_DIR, "1984_2004_storms_v3_72.npy")
PERIOD2_FILE   = os.path.join(INPUT_DIR, "2004_2024_storms_v3_72.npy")

# Segments kept, in each source file's own 1-based time index (README)
RESAMPLED_YEARS_KEPT = (1, 17)    # 1967-1983 (drop the file's own 18-30)
PERIOD1_YEARS_KEPT   = (1, 21)    # 1984-2004, all of it
PERIOD2_YEARS_KEPT   = (2, 14)    # 2005-2017 (drop time==1, duplicate of 2004)

# Shifts that make the combined time column one continuous 1..51 index
SHIFT_PERIOD1 = 17   # time 1-21 -> 18-38
SHIFT_PERIOD2 = 37   # time 2-14 -> 39-51

EXPECTED_TOTAL_YEARS = 51   # 1967 (=1) through 2017 (=51) inclusive
# -----------------------------------------------------------------------------


# Stop if any of the three source storm files is missing
def check_inputs_exist():
    for label, path in [("RESAMPLED_FILE", RESAMPLED_FILE),
                        ("PERIOD1_FILE", PERIOD1_FILE),
                        ("PERIOD2_FILE", PERIOD2_FILE)]:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Missing {label}: {os.path.abspath(path)}")
    print("  All three source storm files found.")


# Confirm 2004 is the same storm record in both periods before dropping one copy
def verify_no_duplicate_storms(period1, period2):
    p1_2004 = period1[period1[:, 0] == PERIOD1_YEARS_KEPT[1]][:, 1:]
    p2_2004 = period2[period2[:, 0] == 1][:, 1:]

    if p1_2004.shape != p2_2004.shape or not np.array_equal(p1_2004, p2_2004):
        raise ValueError(
            "Period 1's final year and Period 2's first year do NOT match "
            "-- they are not simply duplicated 2004 records. Stop and check "
            "the assumption behind PERIOD2_YEARS_KEPT before proceeding."
        )
    print(f"  Verified: Period 1 (time=21) and Period 2 (time=1) 2004 storms "
          f"are identical ({p1_2004.shape[0]} storms) -- safe to drop one copy.")


# Run: load, check the 2004 overlap, cut and shift three segments, validate, save
def build_combined_storms():
    check_inputs_exist()

    resampled = np.load(RESAMPLED_FILE)
    period1   = np.load(PERIOD1_FILE)
    period2   = np.load(PERIOD2_FILE)
    print(f"  Loaded resampled: {resampled.shape}, Period 1: {period1.shape}, "
          f"Period 2: {period2.shape}")

    verify_no_duplicate_storms(period1, period2)

    # Segment A: resampled pre-1984 (1967-1983), no shift needed
    lo, hi = RESAMPLED_YEARS_KEPT
    segA = resampled[(resampled[:, 0] >= lo) & (resampled[:, 0] <= hi)].copy()

    # Segment B: real Period 1 (1984-2004), shifted to follow segment A
    lo, hi = PERIOD1_YEARS_KEPT
    segB = period1[(period1[:, 0] >= lo) & (period1[:, 0] <= hi)].copy()
    segB[:, 0] += SHIFT_PERIOD1

    # Segment C: real Period 2 (2005-2017), shifted to follow segment B
    lo, hi = PERIOD2_YEARS_KEPT
    segC = period2[(period2[:, 0] >= lo) & (period2[:, 0] <= hi)].copy()
    segC[:, 0] += SHIFT_PERIOD2

    combined = np.vstack([segA, segB, segC])
    combined = combined[np.argsort(combined[:, 0], kind="stable")]

    # Validate before saving: every model year 1..51 present, nothing extra
    years_present = set(combined[:, 0].astype(int).tolist())
    expected_years = set(range(1, EXPECTED_TOTAL_YEARS + 1))
    missing = expected_years - years_present
    extra   = years_present - expected_years
    if missing:
        raise ValueError(f"Combined storm file is missing model year(s): {sorted(missing)}")
    if extra:
        raise ValueError(f"Combined storm file has unexpected model year(s): {sorted(extra)}")

    print(f"\n  Segment A (1967-1983): {segA.shape[0]:4d} storms, years "
          f"{int(segA[:,0].min())}-{int(segA[:,0].max())}")
    print(f"  Segment B (1984-2004): {segB.shape[0]:4d} storms, years "
          f"{int(segB[:,0].min())}-{int(segB[:,0].max())}")
    print(f"  Segment C (2005-2017): {segC.shape[0]:4d} storms, years "
          f"{int(segC[:,0].min())}-{int(segC[:,0].max())}")
    print(f"\n  Combined: {combined.shape[0]} storms total, "
          f"{EXPECTED_TOTAL_YEARS} model years (1967=1 .. 2017={EXPECTED_TOTAL_YEARS}), "
          f"no gaps or duplicates.")

    np.save(OUTPUT_PATH, combined)
    print(f"\n  Saved: {os.path.abspath(OUTPUT_PATH)}  shape={combined.shape}")

    print("\n  REMINDER: HAT_groin_hindcast_1967_1997.py computes RUN_YEARS = "
          "END_YEAR - START_YEAR (exclusive of END_YEAR). To simulate all 51 "
          "years in this file (through 2017 inclusive), set END_YEAR = 2018, "
          "not 2017.")

    return combined


if __name__ == "__main__":
    build_combined_storms()
