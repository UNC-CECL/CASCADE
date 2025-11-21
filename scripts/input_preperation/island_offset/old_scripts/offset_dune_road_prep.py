import pandas as pd
import numpy as np
from pathlib import Path
import os

"""
Script to clean and align the Dune Offset and Road Setback data files.
This ensures both files contain exactly 104 rows corresponding to the 
CASCADE model's target domains (ID 30 through 133), which resolves 
the "Road setback length != REAL" assertion error.
"""

BASE = Path(r"/").resolve()
os.chdir(str(BASE))
DATA = BASE / "data" / "hatteras_init"

# -------------------------------------------------------------------
# 1. DEFINE THE TARGET DOMAINS
# -------------------------------------------------------------------

# Target domain IDs for the model (30 through 133 inclusive = 104 domains)
TARGET_DOMAIN_IDS = list(range(30, 134))
TARGET_DF = pd.DataFrame(TARGET_DOMAIN_IDS, columns=['domain_id'])
print(f"Targeting {len(TARGET_DOMAIN_IDS)} domains: {TARGET_DOMAIN_IDS[0]} to {TARGET_DOMAIN_IDS[-1]}")

# Define input and output paths based on your file structure
OFFSET_IN_PATH = DATA / "island_offset" / "1978" / "Relative_Dune_1978_Offset.csv"
SETBACK_IN_PATH = DATA / "roads" / "offset" / "1978_road_setback.csv"

# New, cleaned output files that your CASCADE script will use
OFFSET_OUT_PATH = DATA / "island_offset" / "1978" / "Relative_Dune_1978_Offset_CLEAN.csv"
SETBACK_OUT_PATH = DATA / "roads" / "offset" / "1978_road_setback_CLEAN.csv"

# -------------------------------------------------------------------
# 2. PROCESS DUNE OFFSET FILE (Relative_Dune_1978_Offset.csv)
# -------------------------------------------------------------------
print("\nProcessing Dune Offset file...")
try:
    # Load the file you uploaded
    offset_df = pd.read_csv(OFFSET_IN_PATH)
except FileNotFoundError:
    print(f"ERROR: Dune Offset file not found at {OFFSET_IN_PATH}. Please check path.")
    exit()

# Rename the domain ID column (from 'ID' to 'domain_id') for merging
offset_df.rename(columns={'ID': 'domain_id'}, inplace=True)

# Merge: This ensures all 104 TARGET_DOMAIN_IDS are kept, maintaining order.
aligned_offset_df = TARGET_DF.merge(
    offset_df,
    on='domain_id',
    how='left'
)

# Fill any missing values with 0.
if '1978' in aligned_offset_df.columns:
    if aligned_offset_df['1978'].isnull().any():
        print("WARNING: Missing values found in Dune Offset! Filling with 0.")
        aligned_offset_df['1978'].fillna(0, inplace=True)

    # Keep only the value column ('1978') as expected by your final CASCADE script's simple loading logic
    aligned_offset_df = aligned_offset_df[['1978']]
    aligned_offset_df.to_csv(OFFSET_OUT_PATH, index=False)
    print(f"Dune Offset data saved to: {OFFSET_OUT_PATH.name} (Length: {len(aligned_offset_df)})")
else:
    print("ERROR: Could not find '1978' column in Dune Offset file. Check header name.")

# -------------------------------------------------------------------
# 3. PROCESS ROAD SETBACK FILE (1978_road_setback.csv)
# -------------------------------------------------------------------
print("\nProcessing Road Setback file...")
try:
    # Load the file you uploaded
    setback_df = pd.read_csv(SETBACK_IN_PATH)
except FileNotFoundError:
    print(f"ERROR: Road Setback file not found at {SETBACK_IN_PATH}. Please check path.")
    exit()

# Merge: This aligns the existing setback data to the full 104 domain list.
merged_setbacks = TARGET_DF.merge(
    setback_df,
    on='domain_id',
    how='left'
)

# Fill missing values (NaN) with 0.0 meters. This is crucial for domains not in the original setback file.
if merged_setbacks['Road_Setback_m'].isnull().sum() > 0:
    missing_count = merged_setbacks['Road_Setback_m'].isnull().sum()
    print(f"INFO: Missing {missing_count} road setback values. Filling with 0.0m for alignment.")
    merged_setbacks['Road_Setback_m'].fillna(0.0, inplace=True)

# Save the final clean file (retaining the domain_id column for clarity in the file)
merged_setbacks.to_csv(SETBACK_OUT_PATH, index=False)
print(f"Road Setback data saved to: {SETBACK_OUT_PATH.name} (Length: {len(merged_setbacks)})")

print("\nData preparation complete. You can now run your CASCADE script using the '_CLEAN.csv' files.")