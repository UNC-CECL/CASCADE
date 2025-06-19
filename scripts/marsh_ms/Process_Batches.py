import copy

import matplotlib.pyplot as plt
import numpy as np
import os
import copy

import pandas as pd

os.chdir("E:\\Chapter 2\\")
Base_Save_Path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 2\\Updated_Marsh_Accretion_C_Release\\'
#save_path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 2\\Cascade_CSV_Outputs\\Hindcasts\\'

#run_name = 'Geom_1_IL_Baseline_S20.npz' #'Geom_4_IL_10_S49_New_Sink.npz'
#for geos in range(len(Base_Name_List)):
#    temp_run_name = copy.deepcopy(Base_Name_List[geos]+'_Calibrated_Hindcast_2.npz')
#    run_name.append(copy.deepcopy(temp_run_name))

Storm_List = ['Baseline',
              '5',
              '10']

RSLR = [#'IL',
        #'I',
        'IH']

Base_Name_List = ['Geom_1',
                  'Geom_3',
                  'Geom_4',
                  'Geom_5']

Base_Name_List = ['Geom_3']

for Geometries in range(len(Base_Name_List)):
    Base_Name = Base_Name_List[Geometries]

    if Base_Name == 'Geom_1':
        Cascade_Offset = 90
    if Base_Name == 'Geom_3':
        Cascade_Offset = 490
    if Base_Name == 'Geom_4':
        Cascade_Offset = 170
    if Base_Name == 'Geom_5':
        Cascade_Offset = 160
    for RSLR_Rates in range(len(RSLR)):
        for Storm_Mean_Intensities in range(len(Storm_List)):
            Save_Name = Base_Name+'_'+RSLR[RSLR_Rates]+'_'+Storm_List[Storm_Mean_Intensities]
            All_Runs_Shoreline_Location = []
            All_C_Values = []
            for S_Nums in range(0,50):
                # Load a specific Model Run
                if Base_Name == 'Geom_4':
                    Load_Name = Base_Name+'_'+RSLR[RSLR_Rates]+'_'+Storm_List[Storm_Mean_Intensities]+'_S'+str(S_Nums)+'_New_Sink.npz'
                else:
                    Load_Name = Base_Name+'_'+RSLR[RSLR_Rates]+'_'+Storm_List[Storm_Mean_Intensities]+'_S'+str(S_Nums)+'.npz'


                output = np.load(Load_Name, allow_pickle=True)
                cascade = output["cascade"]
                cascade = cascade[0]

                bmft = cascade._bmft_coupler
                bmftc = bmft._bmftc[0]
                O_flux = bmftc.fluxes
                Forest_e = bmftc.Forest_edge
                Marsh_e = bmftc.Marsh_edge

                initial_C = bmftc._marshOM_initial # kg's C accross entire platfrom
                initial_marsh_C_g = initial_C*1000

                # Create initial C layer
                Initial_Marsh_Width = Forest_e[49] - Marsh_e[49]
                Per_m_Marsh_C_g = initial_marsh_C_g/Initial_Marsh_Width
                Initial_Transect_C_Value = np.zeros(len(bmftc.elevation[0]))

                # Add initial C to marsh cells
                Initial_Transect_C_Value[int(Marsh_e[49]):int(Forest_e[49])] = Per_m_Marsh_C_g
                Total_C_Deposited_TS = [Initial_Transect_C_Value]

                #Bmftc._BayOM[yr]

                Shoreline_location_TS = cascade.brie.x_s_save
                All_Reletive_Shoreline = ((-Shoreline_location_TS+1624+Cascade_Offset)+Forest_e[49])[0]


                # Start Year = [50]
                Start_Year = 50
                End_Year = Start_Year+125

                C_autoch = bmftc._organic_dep_autoch# Subtract eroded mass from depositional record
                C_alloch = bmftc._organic_dep_alloch

                Total_Annual_C_Change = np.sum((C_alloch,C_autoch),axis=0)

                for years in range(Start_Year,End_Year):
                    temp_sum = np.sum(Total_Annual_C_Change[Start_Year:years],axis=0)
                    marsh_accumulation = temp_sum[int(Marsh_e[years]):int(Forest_e[years])]
                    new_C_deposition_plus_C0 = np.sum((temp_sum,Initial_Transect_C_Value),axis=0)
                    new_C_deposition_plus_C0[new_C_deposition_plus_C0 < 0] = 0
                    Total_C_Deposited_TS.append(copy.deepcopy(new_C_deposition_plus_C0))

                All_Runs_Shoreline_Location.append(copy.deepcopy(All_Reletive_Shoreline))
                All_C_Values.append(copy.deepcopy(Total_C_Deposited_TS))
                print(Save_Name+'_S'+str(S_Nums)+' loaded')

            Mean_Shoreline_Location = np.nanmean(All_Runs_Shoreline_Location, axis=0)
            Mean_C_Values = np.nanmean(All_C_Values, axis=0)

            Full_Save_Path = Base_Save_Path+Save_Name
            # Save each of the arrays including mean and full values
            np.save(Full_Save_Path+'_Mean_C_Added.npy',Mean_C_Values)
            np.save(Full_Save_Path+'_Mean_Shoreline_Location.npy',Mean_Shoreline_Location)
            np.save(Full_Save_Path+'_All_Shoreline_Location.npy',All_Runs_Shoreline_Location)
            np.save(Full_Save_Path+'_All_C_Added.npy',All_C_Values)

            x = 20