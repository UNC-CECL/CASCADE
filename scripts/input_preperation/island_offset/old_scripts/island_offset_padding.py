import pandas as pd
import numpy as np


def pad_data_for_cascade(input_filename, output_filename, padding_zeros=15, target_length=135):
    """
    Loads coastal dune offset data, checks its current length, and pads it with
    zeros at the beginning and end to meet the required length for the CASCADE model.

    Args:
        input_filename (str): The name of the input CSV file.
        output_filename (str): The name for the output padded CSV file.
        padding_zeros (int): The number of zeros to add to the beginning and end.
        target_length (int): The required final length of the array (typically 135).
    """
    try:
        # 1. Load the data
        # Assuming the offsets are in columns named '1978' and '1997'
        df = pd.read_csv(input_filename)

        # Determine the data columns to pad
        data_columns = ['1978', '1997']

        # Check if the expected columns exist
        if not all(col in df.columns for col in data_columns):
            print(f"Error: CSV is missing one or more required columns ({data_columns}).")
            return

        # 2. Check the current length and calculate padding
        current_length = len(df)

        # Calculate the expected total padding needed
        expected_total_padding = target_length - current_length

        if expected_total_padding < 0:
            print(
                f"Warning: Data length ({current_length}) is already greater than the target length ({target_length}). No padding applied.")
            return

        if expected_total_padding != padding_zeros * 2:
            print(
                f"Warning: Current data length is {current_length}. To reach {target_length}, total padding needed is {expected_total_padding}.")
            print(
                f"Your requested padding ({padding_zeros} start + {padding_zeros} end = {padding_zeros * 2}) matches this.")

        print(f"Original data length: {current_length} points.")
        print(f"Applying {padding_zeros} zeros to the beginning and {padding_zeros} to the end.")

        # Create a DataFrame for the padding
        zero_padding = pd.DataFrame(np.zeros((padding_zeros, len(data_columns))), columns=data_columns)

        # 3. Apply padding to the DataFrame
        # Concatenate the start padding, original data, and end padding
        padded_df = pd.concat([zero_padding, df[data_columns], zero_padding], ignore_index=True)

        # 4. Final verification
        final_length = len(padded_df)
        print(f"Final padded data length: {final_length} points.")

        if final_length != target_length:
            print(f"Error: Final length is {final_length}, but target was {target_length}. Check input data structure.")
            return

        # 5. Save the new padded data
        padded_df.to_csv(output_filename, index=False)
        print(f"\nSuccessfully generated {output_filename} with {final_length} values.")

    except FileNotFoundError:
        print(f"Error: The input file '{input_filename}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


# --- Script Execution ---
input_file = r'/data/hatteras_init/island_offset/Dune_Offsets_1978_1997_CASCADE_Input.csv'
output_file = r'/data/hatteras_init/island_offset/Dune_Offsets_1978_1997_PADDED_135.csv'
padding_value = 15
target_value = 135

pad_data_for_cascade(input_file, output_file, padding_value, target_value)