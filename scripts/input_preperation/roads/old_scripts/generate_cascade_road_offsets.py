"""
Script: generate_cascade_road_offsets.py
Author: Adapted for Hannah from Benton's scripts
Purpose: Compute road setback values for CASCADE model by subtracting dune offsets from road distances.

Inputs:
- Road offset CSV: distances from offshore datum to road (in meters)
- Dune offset CSV: distances from offshore datum to dune toe (in meters)

Output:
- A CSV with the road setback from dune toe (in meters) for each domain
- Ready to be multiplied by 10 for CASCADE (which uses cm)
"""

import pandas as pd
import numpy as np
import unicodedata
import os

# === USER INPUTS ===

# Full paths to your input CSVs
road_csv = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\roads\offset\1978_road_offset_raw.csv'
dune_csv = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\1978\1978_duneline_offset_raw.csv'

# Output path for final CSV (will auto-create folders if missing)
output_path = r'/data/hatteras_init/roads/offset/1978_road_setback.csv'

# Column names
id_column = 'domain_id'
road_column = 'ORIG_LEN'
dune_column = 'ORIG_LEN'

# === LOAD AND CLEAN DATA ===

# Load CSVs
road_df = pd.read_csv(road_csv)
dune_df = pd.read_csv(dune_csv)

# Normalize and strip column names
road_df.columns = [unicodedata.normalize("NFKD", col).strip() for col in road_df.columns]
dune_df.columns = [unicodedata.normalize("NFKD", col).strip() for col in dune_df.columns]

# Drop any unnamed index columns
road_df = road_df.loc[:, ~road_df.columns.str.contains('^Unnamed')]
dune_df = dune_df.loc[:, ~dune_df.columns.str.contains('^Unnamed')]

# =========================================================
# === FIX: Use 'domain_id_1' for Road Data Grouping ===
# =========================================================

# Group road distances by the CORRECT ID column: 'domain_id_1'
road_grouped = road_df.groupby('domain_id_1')[road_column].min().reset_index()
# Rename offset column AND rename the ID column to 'domain_id' for merging
road_grouped.rename(columns={road_column: 'road_offset', 'domain_id_1': id_column}, inplace=True)

# Group dune distances by the correct ID column: 'domain_id'
dune_grouped = dune_df.groupby(id_column)['ORIG_LEN'].min().reset_index()
dune_grouped.rename(columns={'ORIG_LEN': 'dune_offset'}, inplace=True)

# Merge both grouped datasets
merged = pd.merge(road_grouped, dune_grouped, on=id_column, how='inner')

# Calculate setback in meters
merged['Road_Setback_m'] = merged['road_offset'] - merged['dune_offset']

# Sort by domain ID
merged = merged.sort_values(by=id_column)

# Ensure output folder exists
os.makedirs(os.path.dirname(output_path), exist_ok=True)

# Save final result with domain ID and setback in meters
merged[['domain_id', 'Road_Setback_m']].to_csv(output_path, index=False)

# Debug / confirmation
print(f"Saved CASCADE-ready road setback array to: {output_path}")
print("Merged columns:", merged.columns.tolist())