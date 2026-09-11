"""
HAT_check_nourishment_applied.py
==================================
Checks whether historical beach nourishment (1971, 1973) actually reached
x_s inside CASCADE, or was silently skipped by NourishmentManager's
narrow_break early-return (beach_dune_manager.py line ~650) -- which exits
BEFORE the nourish block, regardless of nourish_now/volume having been set
correctly by the runner script.

Console-log prints from the hindcast runner only confirm what the RUNNER
requested (set before cascade.update() is called that year) -- they cannot
see whether Cascade's internal update() actually executed the nourish
branch. This script checks the one signal that can: nourishment_volume_TS,
which is only written on the exact line inside the nourish branch, so a
zero there when a nonzero volume was requested means it was blocked.

Usage: point RUN_DIR at the saved run folder (the one containing the
{run_name}.npz file), then run.
"""

import os
import numpy as np

RUN_DIR = r"/output/raw_runs/HAT_1967_2018_edge_calibrated_no_groin"
RUN_NAME = "HAT_1967_2018_edge_calibrated_no_groin"

START_YEAR = 1967
NUM_BUFFER_DOMAINS = 15
FIRST_FILE_NUMBER = 2

# Domains and years to check (matches the historical nourishment schedule)
CHECK_DOMAINS_GIS = [6, 7, 8, 9, 10]
CHECK_YEARS = [1971, 1973]

# Expected requested volumes (m^3/m), from the actual console log -- used
# only to print a side-by-side comparison, not to change the check itself.
EXPECTED_VOLUME_M3_PER_M = {
    1971: {6: 0.0, 7: 0.0, 8: 305.8, 9: 0.0, 10: 0.0},
    1973: {6: 397.6, 7: 397.6, 8: 397.6, 9: 397.6, 10: 397.6},
}


def _gis_to_pad(gis_id):
    return NUM_BUFFER_DOMAINS + (gis_id - FIRST_FILE_NUMBER)


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
