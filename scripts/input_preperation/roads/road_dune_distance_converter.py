# =============================================================================
# Road-Dune Distance Transformer
# -----------------------------------------------------------------------------
# Description:
#   Converts a long-format CSV (domain_id, Road_Setback_m) into the specific
#   two-line wide format required for CASCADE road-dune distance input.
# =============================================================================

import pandas as pd
import os

# --- Configuration ---
INPUT_FILE = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\roads\offset\1978_road_setback_CLEAN.csv"
OUTPUT_FILE = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\roads\offset\Road_Dune_Distance_1978.csv"

def convert_to_wide_format(input_path, output_path):
    """
    Loads the input CSV, converts it to the two-line wide format, and saves it.
    """
    if not os.path.exists(input_path):
        print(f"Error: Input file not found at {input_path}")
        return

    # 1. Load the data
    try:
        df = pd.read_csv(input_path)
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    # Check if required columns exist
    if 'domain_id' not in df.columns or 'Road_Setback_m' not in df.columns:
        print("Error: Input file must contain 'domain_id' and 'Road_Setback_m' columns.")
        return

    # 2. Sort by domain_id to ensure the order is correct
    # The 'domain_id' is assumed to be the transect identifier, and they must be sequential.
    df = df.sort_values(by='domain_id', ascending=True)

    # 3. Extract and format the Domain IDs (Header Line)
    # The IDs should be integers without decimal points.
    domain_ids = df['domain_id'].astype(int).astype(str).tolist()
    header_line = ",".join(domain_ids)

    # 4. Extract and format the Setback Distances (Data Line)
    # Ensure they are formatted as floats to match the original example (49.0, 24.0, etc.)
    distance_values = df['Road_Setback_m'].apply(lambda x: f"{x:.1f}" if pd.notna(x) else "").tolist()
    data_line = ",".join(distance_values)

    # 5. Write the two lines to the output file
    try:
        with open(output_path, 'w') as f:
            f.write(header_line + "\n")
            f.write(data_line + "\n")
        print(f"\nSuccessfully converted and saved to: {output_path}")
        print(f"Total domains processed: {len(domain_ids)}")
    except Exception as e:
        print(f"Error writing output file: {e}")

# Run the conversion
convert_to_wide_format(INPUT_FILE, OUTPUT_FILE)