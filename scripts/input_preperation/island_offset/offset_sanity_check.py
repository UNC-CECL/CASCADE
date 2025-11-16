import pandas as pd
import os


def check_csv_dimensions(file_name):
    """
    Loads a CSV file and prints its dimensions (shape), column names,
    and data types to help debug array mismatch issues.

    Args:
        file_name (str): The name of the CSV file to load.
    """
    print(f"--- Loading file: {file_name} ---")

    # Simple check to ensure the file exists before attempting to load
    if not os.path.exists(file_name):
        print(f"ERROR: File not found at the specified path: {file_name}")
        print("Please check the file name and path.")
        return

    try:
        # Load the CSV file into a pandas DataFrame.
        # It's good practice to explicitly specify the encoding for robustness.
        df = pd.read_csv(file_name, encoding='utf-8')

        print("\n[1] Overall Data Frame Shape (Rows, Columns):")
        # df.shape returns a tuple (number of rows, number of columns)
        print(df.shape)

        print("\n[2] Detailed Dimensions:")
        num_rows = df.shape[0]
        num_cols = df.shape[1]
        print(f"Number of Rows (Island Cross-Sections/Observations): {num_rows}")
        print(f"Number of Columns (Time periods/Fields): {num_cols}")

        print("\n[3] Column Names and Data Types (dtypes):")
        # This is key for debugging, as mismatched types (e.g., 'object' instead of 'float64')
        # often indicate non-numeric data in a column.
        print(df.dtypes)

        print("\n[4] First 5 Rows of Data:")
        # Displaying the head helps confirm that data was loaded correctly
        print(df.head())

    except Exception as e:
        print(f"\nAn error occurred during file loading or processing: {e}")
        print("This could be due to incorrect delimiters, encoding, or corrupted data.")


# --- Configuration ---
# Change this file name to check a different CSV
# Available files include:
# 'Dune_Offsets_1978_1997_CASCADE_Input.csv'
# 'Buffer_Shoreline_Offsets.csv'
# 'Relative_Dune_1978_Offset_CLEAN.csv'
# 'Road_Dune_Distance_1974.csv'
file_to_check = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\Dune_Offsets_1978_1997_CASCADE_Input.csv'

# Execute the function
check_csv_dimensions(file_to_check)