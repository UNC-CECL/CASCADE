import numpy as np
import os
import sys

# =========================================================================
# 1. --- USER DEFINED INPUTS ---
# =========================================================================

# The base path to your data (e.g., '/Users/you/Documents/CASCADE_Data')
# This should be the folder containing the 'dunes' and 'topography_dunes' subfolders.
# IMPORTANT: Adjust this path to match your actual setup.
HATTERAS_BASE = os.path.join(r'/data/hatteras_init')

# The list of cross-shore profile domain IDs (e.g., 30 to 79)
# We will analyze both the DUNE and TOPOGRAPHY files for each number in this list.
FILE_NUMBERS = list(range(58, 61))

# The year of the data you are checking (used in the filename pattern)
DATA_YEAR = '2009'


# =========================================================================
# 2. --- CORE ANALYSIS FUNCTION ---
# =========================================================================

def analyze_profile_pair(base_dir, file_num, data_year):
    """
    Loads and analyzes the dune file and the full topography_dunes file for a single domain,
    reporting both min and max elevations for each file.
    """

    # Construct the file paths using your template logic:
    # 1. DUNE FILE (The alongshore array used for management module)
    dune_path = os.path.join(
        base_dir, 'dunes', data_year,
        f'domain_{file_num}_resampled_dune_{data_year}.npy'
    )

    # 2. TOPOGRAPHY FILE (The full cross-shore profile array for backbarrier)
    elev_path = os.path.join(
        base_dir, 'topography_dunes', data_year,
        f'domain_{file_num}_resampled_topography_{data_year}.npy'
    )

    print(f"--- Domain {file_num} ({data_year}) ---")

    # --- Process DUNE File ---
    if os.path.exists(dune_path):
        try:
            dune_data = np.load(dune_path)
            dune_min = np.min(dune_data)
            dune_max = np.max(dune_data)

            print(f"  DUNE Profile Data:")
            print(f"    MIN Elevation (Dune Crest): {dune_min:.4f} m  <--- CHECK THIS FOR DROWNING")
            print(f"    MAX Elevation (Dune Crest): {dune_max:.4f} m")

        except Exception as e:
            print(f"  DUNE FILE ERROR: Could not load or process {os.path.basename(dune_path)}. Error: {e}")
    else:
        print(f"  DUNE FILE NOT FOUND: {dune_path}")

    # --- Process TOPOGRAPHY File ---
    if os.path.exists(elev_path):
        try:
            elev_data = np.load(elev_path)
            elev_min = np.min(elev_data)
            elev_max = np.max(elev_data)

            print(f"  TOPOGRAPHY Data:")
            print(f"    MIN Elevation (Full Profile): {elev_min:.4f} m")
            print(f"    MAX Elevation (Full Profile): {elev_max:.4f} m  <--- Highest Point in Domain")

        except Exception as e:
            print(f"  TOPOGRAPHY FILE ERROR: Could not load or process {os.path.basename(elev_path)}. Error: {e}")
    else:
        print(f"  TOPOGRAPHY FILE NOT FOUND: {elev_path}")

    print("-" * 35)


# =========================================================================
# 3. --- EXECUTION ---
# =========================================================================
if __name__ == '__main__':
    # 1. Basic Path Check
    if not os.path.isdir(HATTERAS_BASE):
        print(f"FATAL ERROR: Base directory not found at '{HATTERAS_BASE}'")
        print("Please edit the HATTERAS_BASE variable (line 12) to the correct path.")
        sys.exit(1)

    print(f"Starting elevation analysis for {len(FILE_NUMBERS)} domains from base: {HATTERAS_BASE}")
    print("----------------------------------------------------------------")

    # 2. Run Analysis for all defined domains
    for file_num in FILE_NUMBERS:
        analyze_profile_pair(HATTERAS_BASE, file_num, DATA_YEAR)

    print("Analysis Complete.")