# =============================================================================
# Storm Period Splitter
# -----------------------------------------------------------------------------
# Description:
#   Loads the master hindcast storm file (1978-2022), translates the relative
#   year index to the absolute calendar year, and splits the data into two
#   subsets (1978-1997 and 1997-2019) for subsequent analysis.
# =============================================================================

import numpy as np
import pandas as pd
import os

# === USER INPUTS AND CONFIGURATION ===========================================

# Master file path (update this if necessary)
MASTER_STORM_PATH = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\storms\hindcast_storms\HAT_1978_2022_Final_Hindcast_Storms.npy"

# Configuration based on file name and user request
START_YEAR = 1978
SPLIT_YEAR = 1997
END_YEAR_PERIOD_B = 2019

# Column names confirmed by your labmate:
# Note: The 'Year_Index' column is relative to START_YEAR
COLS = ["Year_Index", "Rhigh", "Rlow", "Wave Period", "Duration"]

# Output file paths - FYI this saves to input_preperation>storms>
OUTPUT_FILE_A = f"storms_{START_YEAR}_{SPLIT_YEAR}.npy"
OUTPUT_FILE_B = f"storms_{SPLIT_YEAR}_{END_YEAR_PERIOD_B}.npy"

# === 1. LOAD AND PREPARE DATA ================================================
try:
    storm_arr = np.load(MASTER_STORM_PATH)
    print(f"Loaded master storm file: {MASTER_STORM_PATH}")
    print(f"Array shape: {storm_arr.shape}")

    # Convert to DataFrame for easy filtering and column manipulation
    df = pd.DataFrame(storm_arr, columns=COLS)

    # --- CRITICAL: Calculate Absolute Calendar Year ---
    df['Calendar_Year'] = df['Year_Index'].astype(int) + START_YEAR

except FileNotFoundError:
    print(f"Error: Master file not found at {MASTER_STORM_PATH}")
    exit()

# === 2. DEFINE AND FILTER PERIODS (Inclusive of Boundary Year 1997) =========

# Period A: 1978 to 1997 (inclusive)
filter_A = (df['Calendar_Year'] >= START_YEAR) & (df['Calendar_Year'] <= SPLIT_YEAR)
df_A = df[filter_A]

# Period B: 1997 to 2019 (inclusive)
filter_B = (df['Calendar_Year'] >= SPLIT_YEAR) & (df['Calendar_Year'] <= END_YEAR_PERIOD_B)
df_B = df[filter_B]

print(f"\n--- Period A ({START_YEAR}-{SPLIT_YEAR}) ---")
print(f"Storms found: {len(df_A)}")
if not df_A.empty:
    print(f"Earliest year: {df_A['Calendar_Year'].min()}, Latest year: {df_A['Calendar_Year'].max()}")

print(f"\n--- Period B ({SPLIT_YEAR}-{END_YEAR_PERIOD_B}) ---")
print(f"Storms found: {len(df_B)}")
if not df_B.empty:
    print(f"Earliest year: {df_B['Calendar_Year'].min()}, Latest year: {df_B['Calendar_Year'].max()}")


# === 3. CONVERT BACK TO NUMPY AND SAVE =======================================

# We must ensure the output array only contains the original 5 columns (COLS),
# without the temporary 'Calendar_Year' column.

# Period A Save
if not df_A.empty:
    # Convert back to NumPy array using only the original columns
    array_A = df_A[COLS].values
    np.save(OUTPUT_FILE_A, array_A)
    print(f"\nSuccessfully saved {len(array_A)} storms to: {OUTPUT_FILE_A}")
else:
    print(f"\nWarning: No storms found for Period A. File {OUTPUT_FILE_A} not saved.")

# Period B Save
if not df_B.empty:
    # Convert back to NumPy array using only the original columns
    array_B = df_B[COLS].values
    np.save(OUTPUT_FILE_B, array_B)
    print(f"Successfully saved {len(array_B)} storms to: {OUTPUT_FILE_B}")
else:
    print(f"Warning: No storms found for Period B. File {OUTPUT_FILE_B} not saved.")

print("\n--- Process Complete ---")