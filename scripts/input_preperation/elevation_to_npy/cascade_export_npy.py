# Script to export clipped and resampled rasters from ArcGIS Pro to .npy NumPy arrays
# Updated to match output format of the CASCADE Clip & Resample Domains tool

import arcpy
import numpy as np
import os

# Set the folder containing clipped/resampled rasters (subfolders named by domain ID)
root_folder = r'C:\Users\hanna\OneDrive - University of North Carolina at Chapel Hill\Ch1_CASCADE_hatteras\Hatteras_CASCADE_Input\domain_elevation\2009_domain_clipresample'

# Set the folder where .npy arrays will be saved
save_path = r'C:\Users\hanna\OneDrive - University of North Carolina at Chapel Hill\Ch1_CASCADE_hatteras\Hatteras_CASCADE_Input\domain_elevation\2009_npy_arrays'
if not os.path.exists(save_path):
    os.makedirs(save_path)

# Walk through each domain subfolder
for domain_id in os.listdir(root_folder):
    domain_folder = os.path.join(root_folder, domain_id)
    if not os.path.isdir(domain_folder):
        continue  # skip if not a folder

    for filename in os.listdir(domain_folder):
        if filename.endswith('.tif') and filename.startswith('resampled_'):
            raster_path = os.path.join(domain_folder, filename)
            raster = arcpy.Raster(raster_path)

            # Convert raster to numpy array
            arr = arcpy.RasterToNumPyArray(raster, nodata_to_value=-10)

            # Save as .npy file using domain ID in name
            output_name = f'{domain_id}_resampled.npy'
            output_path = os.path.join(save_path, output_name)
            np.save(output_path, arr)
            print(f'Saved: {output_path}')
