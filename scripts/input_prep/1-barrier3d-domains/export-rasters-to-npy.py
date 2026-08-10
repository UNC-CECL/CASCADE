# Script designed to be run in ArcGis Pro to export rasters
# and save them as .npy files for future use in CASCADE or for easy use in Python
import arcpy
import numpy as np
save_path_name = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Ocracoke_NPY_Arrays\\Topo_2019\\'
for i in range(11,51):
    name = 'Grid_'+str(i)+'_Resample'
    inRas = arcpy.Raster(name)
    save_name = save_path_name+name
    # Convert Raster to numpy array
    arr = arcpy.RasterToNumPyArray(inRas,nodata_to_value=-10)
    np.save(arr=arr, file=save_name)
    # Calculate percentage of the row for each cell value
    print(save_name+' is saved')