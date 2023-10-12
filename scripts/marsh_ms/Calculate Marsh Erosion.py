# Calculate the amount of marsh erosion for different CASCADE model runs

import numpy as np
import os
# Set WD
os.chdir("/Users/ceclmac/PycharmProjects/CASCADE/Run_output")

Smith_Names = ['Smith_ACC_RSLR1_S1','Smith_ACC_RSLR2_S1','Smith_ACC_RSLR3_S1',
'Smith_ACC_RSLR1_S2','Smith_ACC_RSLR2_S2','Smith_ACC_RSLR3_S2',
'Smith_ACC_RSLR1_S3','Smith_ACC_RSLR2_S3','Smith_ACC_RSLR3_S3']

Wreck_Names = ['Wreck_ACC_RSLR1_S1', 'Wreck_ACC_RSLR2_S1', 'Wreck_ACC_RSLR3_S1',
               'Wreck_ACC_RSLR1_S2', 'Wreck_ACC_RSLR2_S2', 'Wreck_ACC_RSLR3_S2',
               'Wreck_ACC_RSLR1_S3', 'Wreck_ACC_RSLR2_S3', 'Wreck_ACC_RSLR3_S3']

Metompkin_Marsh_Names = ['Metompkin_Marsh_ACC_RSLR1_S1', 'Metompkin_Marsh_ACC_RSLR2_S1', 'Metompkin_Marsh_ACC_RSLR3_S1',
                         'Metompkin_Marsh_ACC_RSLR1_S2', 'Metompkin_Marsh_ACC_RSLR2_S2', 'Metompkin_Marsh_ACC_RSLR3_S2',
                         'Metompkin_Marsh_ACC_RSLR1_S3', 'Metompkin_Marsh_ACC_RSLR2_S3', 'Metompkin_Marsh_ACC_RSLR3_S3']

All_Marsh_Erosion = []
All_Marsh_Erosion_TS = []

name_prefix =Wreck_Names

if name_prefix == Wreck_Names:
    barrier_width = 15
    marsh_size = 25
elif name_prefix == Metompkin_Marsh_Names:
    barrier_width = 15
    marsh_size = 50
elif name_prefix == Smith_Names:
    barrier_width = 10
    marsh_size = 50


for i in range(len(name_prefix)):
    output = np.load(name_prefix[i] + ".npz", allow_pickle=True)
    cascade = output["cascade"]
    cascade = cascade[0]
    b3d = cascade.barrier3d
    ny = np.size(b3d)
    barrierbmft = cascade._bmftc[2]

    marsh_edge_beginning = b3d[1].x_s_TS[1] + barrier_width
    marsh_edge_end_TS = (barrierbmft.Marsh_edge[barrierbmft.startyear: barrierbmft.endyear] / 100) + b3d[0].x_s_TS[
        0] + barrier_width
    marsh_erosion_ts = []
    for j in range(1, len(b3d[1].x_s_TS)):
        if (b3d[1].x_s_TS[j] - marsh_edge_beginning) > 0 and (marsh_size > sum(marsh_erosion_ts)):
        #if (b3d[1].x_s_TS[j] - marsh_edge_beginning) > 0 :
            marsh_erosion = (b3d[1].x_s_TS[j] - marsh_edge_beginning)
            marsh_edge_beginning = b3d[1].x_s_TS[j]
        else:
            marsh_erosion = 0
        marsh_erosion_ts.append(marsh_erosion)



    #marsh_erosion = marsh_edge_beginning - marsh_edge_end_TS

    All_Marsh_Erosion.append(marsh_erosion_ts)
    #All_Marsh_Erosion_TS.append(marsh_edge_end_TS)
#All_Marsh_Erosion = All_Marsh_Erosion*10 # Convert to m from decameters
np.savetxt('/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/VCR Blue Carbon/Distance Traveled/Wreck_Marsh_Erosion_TS_Traveled.csv',
           All_Marsh_Erosion, delimiter=', ', fmt='% s')
