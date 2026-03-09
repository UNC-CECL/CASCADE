import copy

import matplotlib.pyplot as plt
import numpy as np
import os
import copy

import pandas as pd

os.chdir("E:\\Chapter 2\\")
Base_Save_Path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 2\\Updated_Marsh_Accretion_C_Release\\'

Storm_List = ['Baseline',
              '5',
              '10']

RSLR = ['IL',
        'I',
        'IH']

#RSLR = ['IH']
#Storm_List = ['10']


Base_Name_List = ['Geom_1',
                  'Geom_3',
                  'Geom_4',
                  'Geom_5']

Base_Name_List = ['Geom_3']

for Geometries in range(len(Base_Name_List)):
    Base_Name = Base_Name_List[Geometries]

    if Base_Name == 'Geom_1':
        Cascade_Offset = 120
    if Base_Name == 'Geom_3':
        Cascade_Offset = 500
    if Base_Name == 'Geom_4':
        Cascade_Offset = 190
    if Base_Name == 'Geom_5':
        Cascade_Offset = 170
    for RSLR_Rates in range(len(RSLR)):
        for Storm_Mean_Intensities in range(len(Storm_List)):
            Save_Name = Base_Name+'_'+RSLR[RSLR_Rates]+'_'+Storm_List[Storm_Mean_Intensities]
            All_Runs_Shoreline_Location = []
            All_C_Values = []
            for S_Nums in range(0,50):
                # Load a specific Model Run
                if Base_Name == 'Geom_10':
                    Load_Name = Base_Name+'_'+RSLR[RSLR_Rates]+'_'+Storm_List[Storm_Mean_Intensities]+'_S'+str(S_Nums)+'_New_Sink.npz'
                else:
                    Load_Name = Base_Name+'_'+RSLR[RSLR_Rates]+'_'+Storm_List[Storm_Mean_Intensities]+'_S'+str(S_Nums)+'_N_RSLR_60.npz'


                output = np.load(Load_Name, allow_pickle=True)
                cascade = output["cascade"]
                cascade = cascade[0]

                bmft = cascade._bmft_coupler
                bmftc = bmft._bmftc[0]
                O_flux = bmftc.fluxes
                Forest_e = bmftc.Forest_edge
                Marsh_e = bmftc.Marsh_edge
                if Base_Name != 'Geom_4':
                    initial_C = bmftc._marshOM_initial # kg's C accross entire platfrom
                else:
                    initial_C = 9703/2
                initial_marsh_C_g = initial_C*1000

                # Create initial C layer
                Initial_Marsh_Width = Forest_e[49] - Marsh_e[49]
                Per_m_Marsh_C_g = initial_marsh_C_g/Initial_Marsh_Width
                Initial_Transect_C_Value = np.zeros(len(bmftc.elevation[0]))

                # Add initial C to marsh cells
                Initial_Transect_C_Value[int(Marsh_e[49]):int(Forest_e[49])] = Per_m_Marsh_C_g
                Total_C_Deposited_TS = [Initial_Transect_C_Value]
                Bay_C_Deposit_TS = [np.zeros(len(Initial_Transect_C_Value))]

                Shoreline_location_TS = cascade.brie.x_s_save
                All_Reletive_Shoreline = ((-Shoreline_location_TS+1624+Cascade_Offset)+Forest_e[49])[0]


                # Start Year = [50]
                Start_Year = 49
                End_Year = Start_Year+125

                C_autoch = bmftc._organic_dep_autoch# Subtract eroded mass from depositional record
                C_alloch = bmftc._organic_dep_alloch

                Total_Annual_C_Change = np.sum((C_alloch,C_autoch),axis=0)

                All_C_Deposited = []
                for years in range(Start_Year, End_Year):
                    Temp_C = np.sum(Total_Annual_C_Change[Start_Year:years], axis=0)
                    All_C_Deposited.append(copy.deepcopy(Temp_C))
                    bay_accumulation = np.sum(Total_Annual_C_Change[years:years + 1], axis=0)[:int(Marsh_e[years])]
                    temp_total_bay_C_change = np.zeros(len(Total_Annual_C_Change[0]))
                    temp_total_bay_C_change[:int(Marsh_e[years])] = bay_accumulation
                    temp_total_bay_C_change = np.sum((temp_total_bay_C_change, Bay_C_Deposit_TS[years - Start_Year]),
                                                     axis=0)
                    Bay_C_Deposit_TS.append(copy.deepcopy(temp_total_bay_C_change))

                Bay_C_Deposit_TS_Non_Zero = []
                for bay_years in range(len(Bay_C_Deposit_TS)):
                    Temp_Bay = Bay_C_Deposit_TS[bay_years]
                    Temp_Bay[Temp_Bay < 0] = 0
                    Bay_C_Deposit_TS_Non_Zero.append(copy.deepcopy(Temp_Bay))

                Marsh_C_Minus_Bay_Accretion = np.subtract(All_C_Deposited,
                                                          Bay_C_Deposit_TS_Non_Zero[0:len(All_C_Deposited)])
                Total_Marsh_C_Accretion = np.add(Marsh_C_Minus_Bay_Accretion, Initial_Transect_C_Value)
                Total_Marsh_C_Accretion[Total_Marsh_C_Accretion < 0] = 0

                All_Runs_Shoreline_Location.append(copy.deepcopy(All_Reletive_Shoreline))
                All_C_Values.append(copy.deepcopy(Total_Marsh_C_Accretion))
                print(Save_Name+'_S'+str(S_Nums)+' loaded')

            Mean_Shoreline_Location = np.nanmean(All_Runs_Shoreline_Location, axis=0)
            Mean_C_Values = np.nanmean(All_C_Values, axis=0)

            Full_Save_Path = Base_Save_Path+Save_Name
            # Save each of the arrays including mean and full values
            np.save(Full_Save_Path+'_Mean_C_Added_Marsh_N60_F.npy',Mean_C_Values)
            np.save(Full_Save_Path+'_Mean_Shoreline_Location_Marsh_N60_F.npy',Mean_Shoreline_Location)
            np.save(Full_Save_Path+'_All_Shoreline_Location_Marsh_N60_F.npy',All_Runs_Shoreline_Location)
            np.save(Full_Save_Path+'_All_C_Added_Marsh_N60_F.npy',All_C_Values)

            x = 20