import numpy as np
import pandas as pd
import os


def load_road_setbacks(filepath, min_domain_id=30, max_domain_id=128, save_filepath=None):
    """
    Loads, filters (for domain IDs 30-128), and flattens the Road Setbacks
    CSV into a 1D NumPy array. Optionally saves the filtered data to a new CSV.

    Args:
        filepath (str): Path to the input CSV file.
        min_domain_id (int): Starting domain ID (inclusive).
        max_domain_id (int): Ending domain ID (inclusive).
        save_filepath (str, optional): If provided, the filtered data will be
                                       saved to this path. Defaults to None.

    Returns:
        np.ndarray: A 1D array of road setback values for the specified domains.
    """
    if not os.path.exists(filepath):
        print(f"Error: Road Setback file not found at {filepath}")
        return None

    try:
        # 1. Load the entire CSV
        df = pd.read_csv(filepath)

        # 2. Filter the DataFrame for the required domain IDs (30 <= domain_id <= 128)
        filtered_df = df[
            (df['domain_id'] >= min_domain_id) &
            (df['domain_id'] <= max_domain_id)
            ].copy()

        # Check for the expected size
        expected_size = max_domain_id - min_domain_id + 1
        if len(filtered_df) != expected_size:
            print(
                f"Warning: Data size mismatch after filtering. Expected {expected_size} points, but found {len(filtered_df)}.")

        # 3. *** NEW STEP: SAVE THE FILTERED DATA IF A PATH IS PROVIDED ***
        if save_filepath:
            filtered_df.to_csv(save_filepath, index=False)
            print(f"SUCCESS: Filtered data saved to {save_filepath}")

        # 4. Extract the 'Road_Setback_m' column and convert to a flat 1D NumPy array
        road_setback_data = filtered_df['Road_Setback_m'].to_numpy().flatten()

        return road_setback_data

    except KeyError as e:
        print(f"Error: Column {e} not found in the CSV. Check header names.")
        return None
    except Exception as e:
        print(f"Error loading road setback file: {e}")
        return None


# --- SCRIPT EXECUTION BLOCK ---
if __name__ == "__main__":
    # Define your file paths
    INPUT_FILEPATH = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\roads\offset\1978_road_setback_CLEAN.csv'

    # Define the path where you want the NEW filtered file to be saved
    OUTPUT_FILEPATH = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\roads\offset\1978_road_setback_FINAL.csv'

    # Call the function, passing the output path
    road_setback_array = load_road_setbacks(
        filepath=INPUT_FILEPATH,
        save_filepath=OUTPUT_FILEPATH  # Pass the new file path here
    )

    # Print the verification output
    if road_setback_array is not None:
        print("-" * 50)
        print("Road Setback Data Verification (Array in Memory):")
        print(f"Data type: {type(road_setback_array)}")
        print(f"Array shape: {road_setback_array.shape}")
        print(f"First 5 values: {road_setback_array[:5]}")
        print("-" * 50)