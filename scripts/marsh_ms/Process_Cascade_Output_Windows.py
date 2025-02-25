import numpy as np
import os

#os.chdir("C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 2\\Cascade_Output_NPZ")
os.chdir('E:\\Chapter 2')
#island_names = ['Smith','Hog','Wreck','Metompkin_Marsh','Metompkin_Bay']
island_names = ['Geom_1','Geom_2','Geom_3','Geom_4','Geom_5']
island_names = ['Geom_5']

Storm_Intensity = ['_Baseline_','_5_','_10_']

Save_RSLR = ['IL','I','H']
RSLR = ['IL','I','IH']
for j in range(len(island_names)):
    for k in range(len(RSLR)):
        #for z in range(len(Storm_Intensity)):
        bmft_elevation_TS = []
        marsh_offset_TS = []
        marsh_edge_TS = []
        forest_edge_TS = []
        shoreline_TS = []
        shoreline_toe_TS = []
        for i in range(0,50):
            run_name = (str(island_names[j])+'_'+str(RSLR[k])+str(Storm_Intensity[k])+'S'+str(i)+'.npz')
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
                z = 20

        save_name = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 2\\Cascade_TS_Output\\'+str(island_names[j])+'_'+str(Save_RSLR[k])+str(Storm_Intensity[k])
        if cascade._marsh_dynamics == True:
            Output_info = [
                shoreline_TS,
                shoreline_toe_TS,
                bmft_elevation_TS,
                marsh_edge_TS,
                forest_edge_TS,
                marsh_offset_TS]
        else:
            Output_info = [
                shoreline_TS,
                shoreline_toe_TS
                ]
        output_var_names = ['shoreline_TS',
                'shoreline_toe_TS',
                'bmft_elevation_TS',
                'marsh_edge_TS',
                'forest_edge_TS',
                'BMFT_offset_TS'
                            ]
        for i in range(len(Output_info)):
            np.savez(file=save_name+output_var_names[i],arr = Output_info[i])
        print('Saved '+str(island_names[j])+'_'+str(Save_RSLR[k])+str(Storm_Intensity[k]))