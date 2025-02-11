# Benton Franklin
# 2/11/2025
# Take .npy values exported from ARCGIS pro and identify the locations
# of the dunes and the interior elevation for use in CASCADE
# User will need to input the path locations of where to load and save data
# as well as which edge of the NPY array corresponds to the ocean
# User needs to know the MHW and island berm elevation of study location

import copy
import os
import numpy as np

load_path_name = 'C:\\Users\\frank\\Downloads\\Arrays-Roya\\Arrays\\'
topo_save_path_name = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Processed Topography\\'
dune_save_path_name = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Processed Dune Values\\Dune_2019\\'
MHW = 0.26 # Meters (NAVD 88)
Beach_Berm_Elev = 1.7 # Meters (NAVD 88)

def Process_Raw_Topo_Data(raw_data_path,
                          dune_save_path,
                          topo_save_path,
                          MHW = 0.26,
                          Beach_Berm_Elevation = 1.7,
                          Seaward_Edge = 'Bottom'):
    # Find all Numpy Arrays in raw data folder
    file_names = os.listdir(raw_data_path)
    NPY_Files_Array = []
    for file in file_names:
        if file.endswith('.npy'):
            NPY_Files_Array.append(copy.deepcopy(file))

    # Process the raw elevation data
    for j in range(len(NPY_Files_Array)):
        name = NPY_Files_Array[j]
        # Save the base name before the file extentsion
        final_index = -5
        for chars in range(len(name)):
            if name[chars] == '.':
                final_index = copy.deepcopy(chars)
        base_name = name[:final_index]

        load_name = load_path_name+name

        Raw_Raster = np.load(load_name)
        Raw_Raster -= MHW
        Raw_Raster[Raw_Raster < -1] = -3

        # Transpose if needed
        if Seaward_Edge == 'Right':
            Raw_Raster = Raw_Raster.transpose()
        elif Seaward_Edge == 'Left':
            Raw_Raster_T = copy.deepcopy(Raw_Raster.transpose())
            Raw_Raster = np.flipud(Raw_Raster_T)
        elif Seaward_Edge == 'Top':
            Raw_Raster = np.flipud(Raw_Raster)
        elif Seaward_Edge == 'Bottom':
            Raw_Raster = Raw_Raster

        # New Test
        ProccesedIslandElevationMatrix = np.full((100,50),fill_value=-3.0)
        DuneHeightsVector = np.full((50),fill_value=-3.0)

        for i in range(0,len(ProccesedIslandElevationMatrix[0])):
            R0 = Raw_Raster[:,i]
            R0 = np.flip(R0)

            # Find Begining of beach
            start_beach = next(x for x,val in enumerate(R0) if val > 0.5)

            # Find 'end' of beach
            end_beach = start_beach + 7
            # Find Max dune height on beach
            DuneElevation = max(R0[start_beach:end_beach])
            # Find Max dune values index
            DuneLocation = next(x for x,val in enumerate(R0) if val == DuneElevation)
            # Start Island Immediately after dune line
            StartIsland = (DuneLocation + 1)

            UseElevation = R0[StartIsland:-1]

            # Add values to larger matrix
            ProccesedIslandElevationMatrix[0:len(UseElevation),i] = UseElevation
            DuneElevation = DuneElevation - (Berm_Elev-MHW)
            if DuneElevation < 0:
                DuneElevation = 0.1

            DuneHeightsVector[i] = DuneElevation
          
        # Convert to decameters
        ProccesedIslandElevationMatrix = ProccesedIslandElevationMatrix * .1
        DuneHeightsVector = DuneHeightsVector * .1

        topo_save_name = topo_save_path+base_name+'_topography_2019.npy'
        dune_save_name = dune_save_path+base_name+'_dune_2019.npy'
        np.save(arr=ProccesedIslandElevationMatrix, file=topo_save_name)
        np.save(arr=DuneHeightsVector, file=dune_save_name)
        print(str(j)+' is processed')

Process_Raw_Topo_Data(raw_data_path = load_path_name,
                          dune_save_path = dune_save_path_name,
                          topo_save_path = topo_save_path_name,
                          Beach_Berm_Elev = Beach_Berm_Elev,
                          Seaward_Edge='Right',
                           MHW = MHW
                        )
