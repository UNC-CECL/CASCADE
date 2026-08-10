"""
HAT_build_1967_2017_storms.py
==============================
Builds the combined 1967-2017 storm file for the extended groin-deterioration
mini hindcast, by stitching together three source files that each use their
own independent 1-based model-time-step index (confirmed against
barrier3d.py: storms are matched via `self._StormSeries[:, 0] == self._time_index`,
where `_time_index` starts at 1 at each run's own first modeled year -- NOT a
calendar year, and NOT 0-based).

SOURCES (see HAT_groin_hindcast_1967_1997.py Section 2c for the nourishment
citation note; storm sourcing is separate from that):
  1967_1997_grointest_storms.npy   - 30-yr bootstrap resample of the 1984-2004
                                      real storm distribution (7.03 storms/yr),
                                      built for the original 1967-1997 test.
                                      Only years 1967-1983 (its own time 1-17)
                                      are reused here -- real storms take over
                                      from 1984 onward, so the resampled years
                                      18-30 (1984-1996) in that file are no
                                      longer needed.
  1984_2004_storms_v3_72.npy       - real storms, Stockdon (2006) runup applied
                                      to WIS wave data + NOAA Duck tide gauge
                                      (station 8651370). Period 1, 21 model
                                      years (1984-2004 inclusive), own time 1-21.
  2004_2024_storms_v3_72.npy       - same methodology, Period 2, 21 model years
                                      (2004-2024 inclusive), own time 1-21.

THE 2004 OVERLAP: Period 1 (time=21) and Period 2 (time=1) both represent the
same real calendar year 2004, confirmed byte-for-byte identical (same 9 storms,
same Rhigh/Rlow/period/duration -- verified below in verify_no_duplicate_storms).
Period 1's copy is kept; Period 2 contributes only its 2005-2017 years (its own
time 2-14) to avoid double-counting 2004's storms.

RESULT: one continuous 1-based time index, 1967=1 through 2017=51 (51 model
years total: 17 resampled + 21 real Period-1 + 13 real Period-2).

IMPORTANT -- END_YEAR convention: HAT_groin_hindcast_1967_1997.py computes
RUN_YEARS = END_YEAR - START_YEAR, which is EXCLUSIVE of END_YEAR (the
original 1967-1997 run, with END_YEAR=1997, only ever simulated through 1996).
To actually simulate through 2017 inclusive (all 51 years in this file), set
END_YEAR = 2018 in that script's Section 3, not 2017.
"""

import os
import numpy as np

# =============================================================================
# CONFIG
# =============================================================================
INPUT_DIR = r"/scripts/groin_module/hindcast_groin_test/groin_init/storms/input_storms"  # folder containing the three source .npy files
OUTPUT_PATH = r"/scripts/groin_module/hindcast_groin_test/groin_init/storms/1967_2017/1967_2017_groin_storms.npy"

RESAMPLED_FILE = os.path.join(INPUT_DIR, "1967_1997_grointest_storms.npy")
PERIOD1_FILE   = os.path.join(INPUT_DIR, "1984_2004_storms_v3_72.npy")
PERIOD2_FILE   = os.path.join(INPUT_DIR, "2004_2024_storms_v3_72.npy")

# Segment boundaries, in each SOURCE file's own 1-based time index.
RESAMPLED_YEARS_KEPT = (1, 17)    # 1967-1983 (drop the file's own 18-30)
PERIOD1_YEARS_KEPT   = (1, 21)    # 1984-2004, all of it
PERIOD2_YEARS_KEPT   = (2, 14)    # 2005-2017 (drop time==1, duplicate of 2004)

# Shifts applied to each segment's time column so the combined file is one
# continuous 1..51 index with no gaps or overlaps.
SHIFT_PERIOD1 = 17   # time 1-21 -> 18-38
SHIFT_PERIOD2 = 37   # time 2-14 -> 39-51

EXPECTED_TOTAL_YEARS = 51   # 1967 (=1) through 2017 (=51) inclusive


# =============================================================================
# HELPERS
# =============================================================================
def check_inputs_exist():
    for label, path in [("RESAMPLED_FILE", RESAMPLED_FILE),
                        ("PERIOD1_FILE", PERIOD1_FILE),
                        ("PERIOD2_FILE", PERIOD2_FILE)]:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Missing {label}: {os.path.abspath(path)}")
    print("  All three source storm files found.")


def verify_no_duplicate_storms(period1, period2):
    """Confirm 2004 (Period 1's time==21, Period 2's time==1) really is the
    same real storm record before we rely on dropping one copy of it."""
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


def build_combined_storms():
    check_inputs_exist()

    resampled = np.load(RESAMPLED_FILE)
    period1   = np.load(PERIOD1_FILE)
    period2   = np.load(PERIOD2_FILE)
    print(f"  Loaded resampled: {resampled.shape}, Period 1: {period1.shape}, "
          f"Period 2: {period2.shape}")

    verify_no_duplicate_storms(period1, period2)

    # Segment A: resampled pre-1984 (1967-1983), no shift needed.
    lo, hi = RESAMPLED_YEARS_KEPT
    segA = resampled[(resampled[:, 0] >= lo) & (resampled[:, 0] <= hi)].copy()

    # Segment B: real Period 1 (1984-2004), shifted to continue after segment A.
    lo, hi = PERIOD1_YEARS_KEPT
    segB = period1[(period1[:, 0] >= lo) & (period1[:, 0] <= hi)].copy()
    segB[:, 0] += SHIFT_PERIOD1

    # Segment C: real Period 2 (2005-2017), shifted to continue after segment B.
    lo, hi = PERIOD2_YEARS_KEPT
    segC = period2[(period2[:, 0] >= lo) & (period2[:, 0] <= hi)].copy()
    segC[:, 0] += SHIFT_PERIOD2

    combined = np.vstack([segA, segB, segC])
    combined = combined[np.argsort(combined[:, 0], kind="stable")]

    # --- Validate before saving -----------------------------------------
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
