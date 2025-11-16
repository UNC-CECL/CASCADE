import pandas as pd
import numpy as np
import os

# --- Configuration ---
# IMPORTANT: Adjust these paths if your files are in different directories!
INPUT_FILE_NAME = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\roads\offset\1978_road_setback_FINAL.csv'
OUTPUT_FILE_NAME = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\roads\offset\1978_road_setback_2row_FINAL.csv'
ID_COLUMN = 'domain_id'
VALUE_COLUMN = 'Road_Setback_m'


def preprocess_road_setbacks_to_2row():
    """
    Loads the original CSV, extracts Domain ID and Setback Value columns,
    cleans the data, and saves it as a two-row CSV where:
    Row 1: Domain IDs
    Row 2: Road Setback Values
    """
    print(f"--- Processing: {INPUT_FILE_NAME} ---")

    if not os.path.exists(INPUT_FILE_NAME):
        print(f"CRITICAL ERROR: Input file not found: {INPUT_FILE_NAME}")
        return

    try:
        # 1. Load the CSV
        df = pd.read_csv(INPUT_FILE_NAME)

        # 2. Validate columns
        if ID_COLUMN not in df.columns or VALUE_COLUMN not in df.columns:
            print(
                f"ERROR: Required columns ('{ID_COLUMN}' or '{VALUE_COLUMN}') not found. Available columns: {list(df.columns)}. Aborting.")
            return

        # 3. Extract and Clean Data
        # Ensure data is numeric (coercing errors to NaN) and drop any bad rows
        df_clean = df[[ID_COLUMN, VALUE_COLUMN]].apply(pd.to_numeric, errors='coerce').dropna()

        # 4. Create the 2D array for output
        # Get the ID and Value columns as 1D NumPy arrays
        ids = df_clean[ID_COLUMN].to_numpy()
        values = df_clean[VALUE_COLUMN].to_numpy()

        # Stack them vertically (as rows)
        two_row_array = np.vstack([ids, values])

        # 5. Save the Clean Data in the 2-Row Format
        # CRUCIAL: Save without the index (index=False) and without column headers (header=False).
        # We use '%.1f' for float formatting to keep it clean.
        np.savetxt(
            OUTPUT_FILE_NAME,
            two_row_array,
            delimiter=',',
            fmt='%.1f'  # Format output numbers neatly
        )

        print(f"SUCCESS: Cleaned 2-Row data saved to {OUTPUT_FILE_NAME}")
        print(f"         Format: 2 rows ({ID_COLUMN}, {VALUE_COLUMN}), {two_row_array.shape[1]} columns/elements.")

    except Exception as e:
        print(f"FATAL ERROR during processing: {e}")


if __name__ == '__main__':
    preprocess_road_setbacks_to_2row()