"""
Hatteras CASCADE Dune Offset Pipeline
====================================

This script:
1. Reads raw dune–baseline intersection CSVs (e.g., 1978, 1997).
2. Calculates the relative dune offset per domain (meters, baseline = minimum).
3. Merges multiple years into a single table (columns = years).
4. Pads the merged table with zeros at the beginning and end for CASCADE
   (15 buffer domains on each side → 135 total rows).

Inputs:
    - One raw CSV per year exported from ArcGIS, with at least:
        'domain_id' : integer domain ID
        'ORIG_LEN'  : distance from baseline to dune (m)
        'LineID'    : transect ID within domain

Outputs (in OUTPUT_DIR):
    - <OUTPUT_BASENAME>_CASCADE_Input_unpadded.csv
        Domain_ID + one column per year (relative offsets, m)
    - <OUTPUT_BASENAME>_CASCADE_Input.csv
        Year columns only (no Domain_ID), unpadded
    - <OUTPUT_BASENAME>_PADDED_135.csv
        Year columns only, padded to 135 rows for CASCADE

Author: Hannah A. Henry (streamlined)
"""

import os
import pandas as pd
import numpy as np

# =============================================================================
# 1. USER CONFIGURATION - LAST RAN FOR WHOLE ISLAND, UPDATE BASED ON STUDY AREA
# =============================================================================

# Mapping: year → raw CSV path
RAW_FILES = {
    1978: r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\1978\1978_duneline_offset_raw.csv",
    1997: r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\1997\1997_duneline_offset_raw.csv",
}

# Output directory
OUTPUT_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\1978_whole_island"

# Base name for output files
OUTPUT_BASENAME = "Whole_Island_Dune_Offsets_1978_1997"

# -------------------------------------------------------------------------
# Domain range for Hatteras (NEW NUMBERING)
# -------------------------------------------------------------------------
# CASCADE expects 105 real domains. For Hatteras these are
# Domain_1 through Domain_105. If you ever shift the window, update
# START_DOMAIN / END_DOMAIN below.
START_DOMAIN = 1
END_DOMAIN = 118
B3D_GRIDS = list(range(START_DOMAIN, END_DOMAIN + 1))

# Padding settings for CASCADE (15 left + 15 right = 30 buffers)
PADDING_ZEROS = 15
TARGET_LENGTH = 148  # 15 + 118 + 15

# Column mapping from raw CSV → standardized names
COL_MAP = {
    "Domain_ID": "domain_id",  # Domain ID
    "Distance": "ORIG_LEN",    # Distance baseline → dune (m)
    "Transect": "LineID",      # Transect ID within each domain
}

# =============================================================================
# 2. FUNCTIONS
# =============================================================================

def calculate_relative_offset(file_path, year, col_map, grids):
    """
    Perform the absolute-to-relative offset calculation for a single year.

    Steps (per year):
        1. Load raw CSV.
        2. Extract domain, distance, transect columns.
        3. For each domain in `grids`:
           - Identify all transects present.
           - For each transect, keep the first record.
           - Compute mean distance across transects.
        4. Compute relative offsets by subtracting the minimum mean distance.
        5. Return DataFrame with ['Domain_ID', '<year>'].

    Parameters
    ----------
    file_path : str
        Path to raw CSV for a single year.
    year : int
        Year label (used as column name).
    col_map : dict
        Mapping from logical names → raw column names.
    grids : list of int
        Domain IDs to process (e.g., 1–105).

    Returns
    -------
    pd.DataFrame or None
    """
    print(f"\n--- Processing {year} ---")
    print(f"Input file: {file_path}")

    try:
        raw_df = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"ERROR: File not found: {file_path}")
        return None

    # Check columns
    for key, col in col_map.items():
        if col not in raw_df.columns:
            print(f"ERROR: Expected column '{col}' for '{key}' not found in {file_path}.")
            return None

    # Standardize columns for internal use
    data_df = pd.DataFrame({
        "B3D_Grid": raw_df[col_map["Domain_ID"]],
        "Distance": raw_df[col_map["Distance"]],
        "Transect": raw_df[col_map["Transect"]],
    })

    mean_distances = []
    seen_domains = []

    # Loop over the target domain IDs for CASCADE
    for grid_id in grids:
        subset = data_df[data_df["B3D_Grid"] == int(grid_id)]

        if subset.empty:
            print(f"  Warning: No data for domain {grid_id} in year {year}.")
            continue

        # Determine transect ID range
        try:
            min_transect = int(subset["Transect"].min())
            max_transect = int(subset["Transect"].max())
        except ValueError:
            print(f"  Warning: Invalid transect data for domain {grid_id} (year {year}). Skipping.")
            continue

        transects = range(min_transect, max_transect + 1)

        # Collect first record per transect
        distances = []
        for t_id in transects:
            t_vals = subset[subset["Transect"] == t_id]
            if not t_vals.empty:
                distances.append(t_vals.iloc[0]["Distance"])

        if not distances:
            print(f"  Warning: No valid transect distances for domain {grid_id} (year {year}). Skipping.")
            continue

        mean_distance = float(np.mean(distances))
        mean_distances.append(mean_distance)
        seen_domains.append(grid_id)

    if not mean_distances:
        print(f"FATAL: No valid domains processed for year {year}.")
        return None

    baseline = min(mean_distances)
    relative_offsets = np.subtract(mean_distances, baseline)

    print(f"Year {year}: {len(mean_distances)} domains processed.")
    print(f"Year {year}: baseline distance = {baseline:.3f} m (min mean).")

    out_df = pd.DataFrame({
        "Domain_ID": seen_domains,
        str(year): relative_offsets,
    })

    return out_df


def merge_years(offset_dfs):
    """
    Merge multiple per-year offset DataFrames on 'Domain_ID'.

    Parameters
    ----------
    offset_dfs : dict
        year → DataFrame with ['Domain_ID', '<year>'].

    Returns
    -------
    pd.DataFrame
        Combined DataFrame with one row per Domain_ID and one column per year.
    """
    merged_df = None

    for year, df in offset_dfs.items():
        if df is None:
            continue

        if merged_df is None:
            merged_df = df.copy()
        else:
            merged_df = pd.merge(
                merged_df,
                df,
                on="Domain_ID",
                how="outer",
            )

    if merged_df is None:
        print("FATAL: No valid year DataFrames to merge.")
        return None

    merged_df = merged_df.sort_values(by="Domain_ID").reset_index(drop=True)
    merged_df = merged_df.fillna(0.0)
    return merged_df


def pad_for_cascade(df, padding_zeros, target_length):
    """
    Pad a DataFrame of year columns with zeros at top and bottom.

    Parameters
    ----------
    df : pd.DataFrame
        Only year columns (e.g., ['1978', '1997']).
    padding_zeros : int
        Number of zero rows to add at top and bottom.
    target_length : int
        Required length after padding (e.g., 135).

    Returns
    -------
    pd.DataFrame or None
    """
    data_columns = list(df.columns)
    current_length = len(df)
    expected_total_padding = target_length - current_length

    if expected_total_padding < 0:
        print(
            f"Warning: current length ({current_length}) > target ({target_length}). "
            f"No padding applied."
        )
        return None

    if expected_total_padding != padding_zeros * 2:
        print(
            f"Warning: To reach length {target_length}, total padding needed is "
            f"{expected_total_padding}, but padding_zeros * 2 = {padding_zeros * 2}."
        )

    print(f"Unpadded length: {current_length} rows.")
    print(f"Padding with {padding_zeros} zeros at top and bottom.")

    zero_block = pd.DataFrame(
        np.zeros((padding_zeros, len(data_columns))),
        columns=data_columns,
    )

    padded = pd.concat(
        [zero_block, df[data_columns], zero_block],
        ignore_index=True,
    )

    final_length = len(padded)
    print(f"Padded length: {final_length} rows.")

    if final_length != target_length:
        print(
            f"ERROR: Final length {final_length} != target {target_length}. "
            f"Check domain range and padding configuration."
        )
        return None

    return padded


# =============================================================================
# 3. MAIN
# =============================================================================

def main():
    # Ensure output directory exists
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR, exist_ok=True)

    # --- 3.1. Compute relative offsets for each year ---
    offset_dfs = {}
    for year, path in RAW_FILES.items():
        df_year = calculate_relative_offset(
            file_path=path,
            year=year,
            col_map=COL_MAP,
            grids=B3D_GRIDS,
        )
        if df_year is None:
            print(f"Year {year}: calculation failed.")
        offset_dfs[year] = df_year

    # --- 3.2. Merge years on Domain_ID ---
    merged = merge_years(offset_dfs)
    if merged is None:
        print("No merged output produced. Exiting.")
        return

    # Save unpadded combined file (Domain_ID + year columns)
    unpadded_path = os.path.join(
        OUTPUT_DIR, f"{OUTPUT_BASENAME}_CASCADE_Input_unpadded.csv"
    )
    merged.to_csv(unpadded_path, index=False)
    print(f"\nUnpadded combined offsets saved to:\n  {unpadded_path}")

    # CASCADE-format (years only, no Domain_ID)
    year_cols = [c for c in merged.columns if c != "Domain_ID"]
    cascade_df = merged[year_cols].copy()

    cascade_unpadded_path = os.path.join(
        OUTPUT_DIR, f"{OUTPUT_BASENAME}_CASCADE_Input.csv"
    )
    cascade_df.to_csv(cascade_unpadded_path, index=False)
    print(f"Unpadded CASCADE-format file saved to:\n  {cascade_unpadded_path}")

    # --- 3.3. Pad for CASCADE (15 left + 15 right = 135 total) ---
    padded_df = pad_for_cascade(
        df=cascade_df,
        padding_zeros=PADDING_ZEROS,
        target_length=TARGET_LENGTH,
    )

    if padded_df is None:
        print("Padded output not created due to errors.")
        return

    padded_path = os.path.join(
        OUTPUT_DIR, f"{OUTPUT_BASENAME}_PADDED_{TARGET_LENGTH}.csv"
    )
    padded_df.to_csv(padded_path, index=False)
    print(f"\nSUCCESS: Padded CASCADE input saved to:\n  {padded_path}")
    print(f"Columns: {', '.join(padded_df.columns)}")
    print(f"Total rows: {len(padded_df)}")


if __name__ == "__main__":
    main()
