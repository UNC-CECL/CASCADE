"""
One cell of the groin sweep: a single (M, fraction) rig run in its own process, scored against the 2018 wet/dry change.

    python HAT_groin_sweep_single_combo.py <M> <fraction>

Launched by HAT_groin_sensitivity_sweep.py, not by hand. Runs the "groin"
key of HAT_groin_hindcast_1967_2017.py, saves the modelled D2-D12 profile to
HAT-buxton-hindcast-groin-test/sensitivity_sweep/profiles/ and prints one
line, RESULT_RMSE=<value>, for the sweep to parse.
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import HAT_groin_hindcast_1967_2017 as hc


# --- CONFIG ------------------------------------------------------------------
# Must match HAT_groin_sensitivity_sweep.py
MODEL_FIT_YEAR    = 2017
OBSERVED_FIT_YEAR = 2018
FIT_DOMAINS_GIS = list(range(2, 13))   # D2-D12, full range

WETDRY_CHANGE_TABLE = os.path.join(
    hc.PROJECT_BASE_DIR, "hard-structures", "groin", "HAT-groin-buxton-output",
    "shoreline_position_output",
    "Change_from_wetdry_1967_D2_D12.csv",
)
# -----------------------------------------------------------------------------


# Observed change 1967 -> OBSERVED_FIT_YEAR per fit domain, landward-positive
def load_observed_target():
    df = pd.read_csv(WETDRY_CHANGE_TABLE).set_index("Domain_ID")
    col = f"change_from_wetdry_1967_wetdry_{OBSERVED_FIT_YEAR}_m"
    return np.array([df[col].get(d, np.nan) for d in FIT_DOMAINS_GIS])


# RMSE over domains where both sides have a value
def rmse(modeled, observed):
    mask = np.isfinite(modeled) & np.isfinite(observed)
    if not np.any(mask):
        return np.nan
    return float(np.sqrt(np.mean((modeled[mask] - observed[mask]) ** 2)))


# Run: set M and fraction on the runner, run the groin key, score, save the profile, print the result
def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <M> <fraction>")
    M = float(sys.argv[1])
    fraction = float(sys.argv[2])

    hc.MAKE_FIGURES = False
    hc.MAKE_RUN_GIF = False
    hc.GROIN_TRAPPING_RATE_M_YR = M
    hc.GROIN_DETERIORATION_FRACTION = fraction

    hc.check_inputs_exist()
    island_offset_dam = hc.load_island_offset_dam()
    elevation_files, dune_files = hc.build_file_lists()
    hist_nourish_on, hist_nourish_vol = hc.build_nourishment_arrays_from_manual_inputs()

    run_name = hc.run_one("groin", island_offset_dam, elevation_files, dune_files,
                           hist_nourish_on, hist_nourish_vol)

    m = np.load(os.path.join(hc.OUTPUT_BASE_DIR, run_name,
                              f"{run_name}_shoreline_matrix.npy"))
    row = MODEL_FIT_YEAR - hc.START_YEAR
    if not (0 <= row < m.shape[0]):
        sys.exit(f"MODEL_FIT_YEAR={MODEL_FIT_YEAR} (row {row}) outside "
                 f"this run's {m.shape[0]} modeled years.")

    pos0 = m[0][hc.START_REAL_INDEX:hc.END_REAL_INDEX]
    posN = m[row][hc.START_REAL_INDEX:hc.END_REAL_INDEX]
    change_raw = posN - pos0

    gis_axis = list(range(hc.FIRST_FILE_NUMBER, hc.LAST_FILE_NUMBER + 1))
    idx = [gis_axis.index(d) for d in FIT_DOMAINS_GIS]
    modeled = change_raw[idx]

    observed = load_observed_target()
    err = rmse(modeled, observed)

    # The profiles folder the sweep reads (README)
    profile_dir = os.path.join(hc.PROJECT_BASE_DIR, "hard-structures", "groin",
                                "HAT-buxton-hindcast-groin-test", "sensitivity_sweep", "profiles")
    os.makedirs(profile_dir, exist_ok=True)
    profile_path = os.path.join(profile_dir, f"M{M:g}_frac{fraction:g}.npy")
    np.save(profile_path, modeled)

    print(f"RESULT_RMSE={err}")


if __name__ == "__main__":
    main()
