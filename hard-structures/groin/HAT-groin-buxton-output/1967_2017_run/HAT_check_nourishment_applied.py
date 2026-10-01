"""
Did the 1971 and 1973 nourishments reach x_s in a saved rig run, or were they skipped inside CASCADE?

    python HAT_check_nourishment_applied.py

Reads the run's saved .npz (RUN_DIR, RUN_NAME) and compares the volume the
runner requested with nourishment_volume_TS, which only the nourish branch
writes. Prints the table; writes nothing.
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""

from pathlib import Path
import os

import numpy as np


# --- CONFIG ------------------------------------------------------------------
# Repo root, found by searching upward (ORGANIZATION.md rule 5)
_PATH_REPO = next(_p for _p in Path(__file__).resolve().parents
                  if (_p / "pyproject.toml").exists())

RUN_DIR = str(_PATH_REPO / "output" / "raw_runs" / "HAT_1967_2018_edge_calibrated_no_groin")
RUN_NAME = "HAT_1967_2018_edge_calibrated_no_groin"

START_YEAR = 1967
NUM_BUFFER_DOMAINS = 15
FIRST_FILE_NUMBER = 2

# Domains and years to check (the historical nourishment schedule)
CHECK_DOMAINS_GIS = [6, 7, 8, 9, 10]
CHECK_YEARS = [1971, 1973]

# Requested volumes (m^3/m) from the run's console log; printed beside, never used in the check
EXPECTED_VOLUME_M3_PER_M = {
    1971: {6: 0.0, 7: 0.0, 8: 305.8, 9: 0.0, 10: 0.0},
    1973: {6: 397.6, 7: 397.6, 8: 397.6, 9: 397.6, 10: 397.6},
}
# -----------------------------------------------------------------------------


# GIS domain id -> padded index
def _gis_to_pad(gis_id):
    return NUM_BUFFER_DOMAINS + (gis_id - FIRST_FILE_NUMBER)


# Run: load the saved Cascade, print requested vs applied per domain and year
def main():
    npz_path = os.path.join(RUN_DIR, f"{RUN_NAME}.npz")
    if not os.path.isfile(npz_path):
        raise FileNotFoundError(f"Saved cascade object not found:\n  {npz_path}")

    data = np.load(npz_path, allow_pickle=True)
    cascade = data["cascade"][0]

    print("=" * 78)
    print("NOURISHMENT APPLICATION CHECK")
    print("=" * 78)

    for year in CHECK_YEARS:
        time_index = year - START_YEAR  # matches barrier3d's time_index convention
        print(f"\nYear {year} (time_index {time_index}):")
        print(f"  {'Domain':<8}{'Requested (m3/m)':<20}{'Actually applied (m3/m)':<26}"
              f"{'narrow_break?':<15}{'Blocked?'}")

        for gis_id in CHECK_DOMAINS_GIS:
            pad = _gis_to_pad(gis_id)
            nourishment_obj = cascade.nourishments[pad]

            requested = EXPECTED_VOLUME_M3_PER_M.get(year, {}).get(gis_id, None)

            applied_ts = getattr(nourishment_obj, "_nourishment_volume_TS", None)
            applied = None
            if applied_ts is not None and 0 <= time_index - 1 < len(applied_ts):
                applied = applied_ts[time_index - 1]

            narrow_break = getattr(nourishment_obj, "_narrow_break", None)

            blocked = "YES" if (requested and requested > 0 and (not applied or applied == 0)) else "no"

            print(f"  D{gis_id:<7}{str(requested):<20}{str(applied):<26}"
                  f"{str(narrow_break):<15}{blocked}")

    print("\nIf 'Actually applied' is 0 while 'Requested' is nonzero, the "
          "nourishment was blocked internally -- narrow_break=1 is the most "
          "likely cause given the code path (see module docstring).")


if __name__ == "__main__":
    main()
