import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir("/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/Chapter 2/Cascade_Output_NPZ")

island_names = ['Smith','Hog']


Storm_Intensity = ['_Baseline_','_Two_Percent_','_Five_Percent_','_Ten_Percent_']
Storm_Intensity_Hog = ['_Baseline_','_2_Percent_','_5_Percent_','_10_Percent_']

RSLR = ['Low','Int','High']
for j in range(len(island_names)):
    for k in range(len(RSLR)):
        for z in range(len(Storm_Intensity)):
            bmft_elevation_TS = []
            marsh_edge_TS = []
            forest_edge_TS = []
            shoreline_TS = []
            shoreline_toe_TS = []
            marsh_offset_TS = []
            for i in range(0,50):
                if island_names[j] == 'Hog':
                    run_name = (str(island_names[j]) + str(Storm_Intensity_Hog[z])+ str(
                        RSLR[k]) + '_RSLR_' + str(i) + '.npz')
                elif island_names[j] == 'Metompkin_Bay':
                    if Storm_Intensity[z] == '_Baseline_':
                        run_name = (str(island_names[j]) + str(Storm_Intensity[z])+ str(
                            RSLR[k]) + '_RSLR_' + str(i) + '.npz')
                    else:
                        run_name = (str(island_names[j]) + str(Storm_Intensity[z]) + 'Storms_' + str(
                            RSLR[k]) + '_RSLR_' + str(i) + '.npz')
                else:
                    run_name = (str(island_names[j])+str(Storm_Intensity[z])+'Storms_'+str(RSLR[k])+'_RSLR_'+str(i)+'.npz')
                output = np.load(run_name, allow_pickle=True)
                cascade = output["cascade"]
                cascade = cascade[0]
                b3d = cascade._barrier3d[0]
                shoreline_TS.append(b3d.x_s_TS)
                shoreline_toe_TS.append(b3d.x_t_TS)
                if cascade._marsh_dynamics == True:
                    marsh_edge_TS.append(cascade._bmft_coupler._bmftc[0].Marsh_edge)
                    forest_edge_TS.append(cascade._bmft_coupler._bmftc[0].Forest_edge)
                    bmft_elevation_TS.append(cascade._bmft_coupler._bmftc[0].elevation)
                    marsh_offset_TS.append(cascade._bmft_coupler._x_s_offset_TS[0][0])

            save_name = '/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/Chapter 2/Cascade_TS_Output/'\
                +str(island_names[j])+str(Storm_Intensity[z])+'Storms_'+str(RSLR[k])+'_RSLR'
            if cascade._marsh_dynamics == True:
                Output_info = [
                shoreline_TS,
                shoreline_toe_TS,
                bmft_elevation_TS,
                marsh_edge_TS,
                forest_edge_TS,
                marsh_offset_TS
                ]
            else:
                Output_info = [
                    shoreline_TS,
                    shoreline_toe_TS
                    ]

            np.savez(file=save_name,arr = Output_info)
            print('Saved '+str(island_names[j])+str(Storm_Intensity[z])+'Storms_'+str(RSLR[k])+'_RSLR')