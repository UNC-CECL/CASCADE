# =============================================================================
# Storm Period Extractor
# -----------------------------------------------------------------------------
# Description:
#   Loads the master hindcast storm file (1978-2022), translates the relative
#   year index to the absolute calendar year, and extracts the 1984-2004
#   subset for subsequent analysis.
# =============================================================================

import numpy as np
import pandas as pd

# === USER INPUTS AND CONFIGURATION ===========================================

# Master file path (update this if necessary)
MASTER_STORM_PATH = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\storms\hindcast_storms\HAT_1978_2022_Final_Hindcast_Storms.npy"

# The year that Year_Index=0 corresponds to in the master file (do not change)
FILE_ORIGIN_YEAR = 1978

# Analysis period
START_YEAR = 1984
END_YEAR = 2004

# Column names confirmed by your labmate:
# Note: The 'Year_Index' column is relative to START_YEAR
COLS = ["Year_Index", "Rhigh", "Rlow", "Wave Period", "Duration"]

# Output file path - saves to input_preparation > storms >
OUTPUT_FILE = f"storms_{START_YEAR}_{END_YEAR}.npy"

# === 1. LOAD AND PREPARE DATA ================================================
try:
    storm_arr = np.load(MASTER_STORM_PATH)
    print(f"Loaded master storm file: {MASTER_STORM_PATH}")
    print(f"Array shape: {storm_arr.shape}")

    # Convert to DataFrame for easy filtering and column manipulation
    df = pd.DataFrame(storm_arr, columns=COLS)

    # --- CRITICAL: Calculate Absolute Calendar Year ---
    df['Calendar_Year'] = df['Year_Index'].astype(int) + FILE_ORIGIN_YEAR

except FileNotFoundError:
    print(f"Error: Master file not found at {MASTER_STORM_PATH}")
    exit()

# === 2. FILTER TO 1984-2004 (inclusive) =====================================

filter_period = (df['Calendar_Year'] >= START_YEAR) & (df['Calendar_Year'] <= END_YEAR)
df_period = df[filter_period]

print(f"\n--- Period ({START_YEAR}-{END_YEAR}) ---")
print(f"Storms found: {len(df_period)}")
if not df_period.empty:
    print(f"Earliest year: {df_period['Calendar_Year'].min()}, Latest year: {df_period['Calendar_Year'].max()}")

# === 3. CONVERT BACK TO NUMPY AND SAVE =======================================

# Output array uses only the original 5 columns (no temporary Calendar_Year column)
if not df_period.empty:
    array_period = df_period[COLS].values
    np.save(OUTPUT_FILE, array_period)
    print(f"\nSuccessfully saved {len(array_period)} storms to: {OUTPUT_FILE}")
else:
    print(f"\nWarning: No storms found for {START_YEAR}-{END_YEAR}. File {OUTPUT_FILE} not saved.")

print("\n--- Process Complete ---")
