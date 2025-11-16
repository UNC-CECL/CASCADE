import pandas as pd
import numpy as np
import os

# --- Configuration ---

# **UNIVERSAL Column Mapping for both 1978 and 1997 files**
# This assumes the Domain ID is now consistently in the 'domain_id' column
# for both raw input files.
COL_MAP = {
    'Domain_ID': 'domain_id',  # The Domain ID column
    'Distance': 'ORIG_LEN',  # The Absolute Offset Distance column (in meters)
    'Transect': 'LineID'  # The Transect ID column
}

# The target domain range (30 through 134, inclusive)
START_DOMAIN = 30
END_DOMAIN = 134
B3D_Grids = list(range(START_DOMAIN, END_DOMAIN + 1))


def calculate_relative_offset(file_path, year, col_map, grids):
    """
    Performs the full absolute-to-relative offset calculation for a single year.
    """

    print(f"\n--- Processing {year} data from {file_path} ---")

    try:
        # 1. Load Data
        Raw_Data = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"Error: File not found: {file_path}. Please ensure the file is in the same directory.")
        return None
    except KeyError as e:
        print(f"Error: Missing expected column in {file_path}. Please check COL_MAP configuration. Missing: {e}")
        return None

    # Extract only the relevant columns and rename for consistency
    Data_DF = pd.DataFrame({
        'B3D_Grid': Raw_Data[col_map['Domain_ID']],
        'Distance': Raw_Data[col_map['Distance']],
        'Transect': Raw_Data[col_map['Transect']]
    })

    # List to store the calculated average distance for each domain ID
    Filtered_Distances = []

    for grid_id in grids:
        # Subset to the current domain
        Subset = Data_DF[Data_DF['B3D_Grid'] == int(grid_id)]

        if Subset.empty:
            print(f"Warning: B3D_Grid ID {grid_id} (required) has no data.")
            Filtered_Distances.append(None)
            continue

        # Find the range of transects used in this domain
        # Use try/except block for robustness against missing transects
        try:
            Min_Transect_Num = int(Subset['Transect'].min())
            Max_Transect_Num = int(Subset['Transect'].max())
            Transects = range(Min_Transect_Num, Max_Transect_Num + 1)
        except ValueError:
            print(f"Warning: B3D_Grid ID {grid_id} has invalid data in the Transect column. Skipping.")
            Filtered_Distances.append(None)
            continue

        Temp_Distance_Array = []
        for transect_id in Transects:
            # Get all entries for this transect
            transect_vals = Subset[Subset['Transect'] == transect_id]

            if not transect_vals.empty:
                # KEEP FIRST ENTRY PER TRANSECT (your original logic)
                distance = transect_vals.iloc[0]['Distance']
                Temp_Distance_Array.append(distance)

        # 3. Calculate Average Dune Distance per Domain
        if Temp_Distance_Array:
            # Calculate the mean distance for the domain and store it
            Mean_Distance = np.mean(Temp_Distance_Array)
            Filtered_Distances.append(Mean_Distance)
        else:
            Filtered_Distances.append(None)

            # --- Step 4: Calculate Relative Offset ---

    # 4a. Separate valid data from missing domains
    Valid_Distances = [d for d in Filtered_Distances if d is not None]
    Valid_Grids = [grids[i] for i, d in enumerate(Filtered_Distances) if d is not None]

    if not Valid_Distances:
        print(f"Fatal Error: No valid domain data found for {year} after filtering. Cannot calculate relative offset.")
        return None

    # 4b. Find the baseline (minimum average offset across all valid domains)
    baseline_distance = min(Valid_Distances)

    # 4c. Calculate relative offset by subtracting the baseline
    Relative_dune_offset = np.subtract(Valid_Distances, baseline_distance)

    print(f"Profiles processed: {len(Valid_Distances)}. Baseline distance: {baseline_distance:.2f} m.")

    # 4d. Create DataFrame for merging (Domain_ID must be preserved)
    relative_dune_df = pd.DataFrame({
        'Domain_ID': Valid_Grids,
        str(year): Relative_dune_offset
    })

    return relative_dune_df


def main():
    """Main function to run the process for both 1978 and 1997 and merge the results."""

    # ----------------------------------------------------------------------
    # Define File Paths and Years
    # ----------------------------------------------------------------------
    data_78 = {
        'path': r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\1978_duneline_offset_raw.csv',
        'year': 1978
    }
    data_97 = {
        'path': r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\1997_duneline_offset_raw.csv',
        'year': 1997
    }

    # ----------------------------------------------------------------------
    # Calculate Offsets for Each Year
    # ----------------------------------------------------------------------
    # Use the single, universal COL_MAP for both years
    df_1978 = calculate_relative_offset(data_78['path'], data_78['year'], COL_MAP, B3D_Grids)
    df_1997 = calculate_relative_offset(data_97['path'], data_97['year'], COL_MAP, B3D_Grids)

    if df_1978 is None or df_1997 is None:
        print("\nProcess aborted because one or both years failed to produce valid data.")
        return

    # ----------------------------------------------------------------------
    # Combine and Export Results
    # ----------------------------------------------------------------------

    df_merged = pd.merge(
        df_1978,
        df_1997,
        on='Domain_ID',
        how='outer'
    ).sort_values(by='Domain_ID').reset_index(drop=True)

    # If a profile is missing in one year, fill it with 0 (zero relative offset).
    df_merged = df_merged.fillna(0)

    # Prepare final CASCADE export format (only offset values, years as header, no Domain_ID)
    output_df = df_merged[['1978', '1997']]

    output_file = '../../../data/hatteras_init/island_offset/Dune_Offsets_1978_1997_CASCADE_Input.csv'

    output_df.to_csv(output_file, index=False)

    print("\n==================================================")
    print("SUCCESS: Combined CASCADE Input File Created")
    print(f"Output File: {output_file}")
    print(f"Total Profiles (Rows): {len(output_df)}")
    print("Columns: 1978, 1997 (Relative Offset in meters)")
    print("==================================================")


if __name__ == '__main__':
    main()