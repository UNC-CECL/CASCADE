"""
HAT_groin_sweep_single_combo.py
=================================
Runs ONE (M, fraction) combination of the groin sensitivity sweep, in its
own fresh Python process, then exits. Meant to be launched as a subprocess
by HAT_groin_sensitivity_sweep_v2.py, NOT run directly in a loop -- the
whole point is that each combination gets a clean process, so accumulated
state (Cascade/Barrier3D objects, joblib worker pools, whatever else builds
up over 30 sequential simulations in one process) can never carry over from
one combination to the next. That's the fix for the 0xC0000005 (Windows
access violation) crash that showed up partway through an earlier
in-process sweep -- process isolation guarantees a full OS-level cleanup
between every single run, regardless of what the underlying cause was.

Prints exactly one line to stdout on success, in an easily-parsed format:
    RESULT_RMSE=<value>
Everything else printed is normal run_one() logging, which the orchestrator
ignores. A non-zero exit code (including an access violation) is treated by
the orchestrator as "this combination failed" -- it moves on rather than
losing the whole sweep.

Usage: HAT_groin_sweep_single_combo.py <M> <fraction>
"""

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import HAT_groin_hindcast_1967_2017 as hc

# =============================================================================
# CONFIG -- must match HAT_groin_sensitivity_sweep_v2.py exactly
# =============================================================================
MODEL_FIT_YEAR    = 2017
OBSERVED_FIT_YEAR = 2018
FIT_DOMAINS_GIS = list(range(2, 13))   # D2-D12, full range

WETDRY_CHANGE_TABLE = os.path.join(
    hc.PROJECT_BASE_DIR, "scripts", "groin", "HAT-hindcast-groin-test",
    "input_prep", "shoreline_position", "output",
    "Change_from_wetdry_1967_D2_D12.csv",
)


def load_observed_target():
    df = pd.read_csv(WETDRY_CHANGE_TABLE).set_index("Domain_ID")
    col = f"change_from_wetdry_1967_wetdry_{OBSERVED_FIT_YEAR}_m"
    return np.array([df[col].get(d, np.nan) for d in FIT_DOMAINS_GIS])


def rmse(modeled, observed):
    mask = np.isfinite(modeled) & np.isfinite(observed)
    if not np.any(mask):
        return np.nan
    return float(np.sqrt(np.mean((modeled[mask] - observed[mask]) ** 2)))


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

    profile_dir = os.path.join(hc.PROJECT_BASE_DIR, "scripts", "groin",
                                "HAT-hindcast-groin-test", "sensitivity_sweep", "profiles")
    os.makedirs(profile_dir, exist_ok=True)
    profile_path = os.path.join(profile_dir, f"M{M:g}_frac{fraction:g}.npy")
    np.save(profile_path, modeled)

    print(f"RESULT_RMSE={err}")


if __name__ == "__main__":
    main()
