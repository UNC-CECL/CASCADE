import numpy as np
import os

os.chdir("/Users/ceclmac/PycharmProjects/CASCADE/Run_output")

Metompkin_Bay_Names = ['Metompkin_Bay_ACC_RSLR1_S1', 'Metompkin_Bay_ACC_RSLR2_S1', 'Metompkin_Bay_ACC_RSLR3_S1',
                       'Metompkin_Bay_ACC_RSLR1_S2', 'Metompkin_Bay_ACC_RSLR2_S2', 'Metompkin_Bay_ACC_RSLR3_S2',
                       'Metompkin_Bay_ACC_RSLR1_S3', 'Metompkin_Bay_ACC_RSLR2_S3', 'Metompkin_Bay_ACC_RSLR3_S3']

Smith_Names = ['Smith_ACC_RSLR1_S1', 'Smith_ACC_RSLR2_S1', 'Smith_ACC_RSLR3_S1',
               'Smith_ACC_RSLR1_S2', 'Smith_ACC_RSLR2_S2', 'Smith_ACC_RSLR3_S2',
               'Smith_ACC_RSLR1_S3', 'Smith_ACC_RSLR2_S3', 'Smith_ACC_RSLR3_S3']

Wreck_Names = ['Wreck_ACC_RSLR1_S1', 'Wreck_ACC_RSLR2_S1', 'Wreck_ACC_RSLR3_S1',
               'Wreck_ACC_RSLR1_S2', 'Wreck_ACC_RSLR2_S2', 'Wreck_ACC_RSLR3_S2',
               'Wreck_ACC_RSLR1_S3', 'Wreck_ACC_RSLR2_S3', 'Wreck_ACC_RSLR3_S3']

Metompkin_Marsh_Names = ['Metompkin_Marsh_ACC_RSLR1_S1', 'Metompkin_Marsh_ACC_RSLR2_S1', 'Metompkin_Marsh_ACC_RSLR3_S1',
                         'Metompkin_Marsh_ACC_RSLR1_S2', 'Metompkin_Marsh_ACC_RSLR2_S2', 'Metompkin_Marsh_ACC_RSLR3_S2',
                         'Metompkin_Marsh_ACC_RSLR1_S3', 'Metompkin_Marsh_ACC_RSLR2_S3', 'Metompkin_Marsh_ACC_RSLR3_S3']
name_prefix = Metompkin_Bay_Names
#name_prefix =Metompkin_Marsh_Names
#name_prefix = Metompkin_Bay_Names
#name_prefix = Wreck_Names
DistanceTraveled = []
Toe_To_Shoreline = []

for i in range(len(name_prefix)):
    output = np.load(name_prefix[i] + ".npz", allow_pickle=True)
    cascade = output["cascade"]
    cascade = cascade[0]
    b3d = cascade.barrier3d
    ny = np.size(b3d)
    if ny == 4:
        mean_x_s_TS = [b3d[0].x_s_TS,b3d[1].x_s_TS,b3d[2].x_s_TS,b3d[3].x_s_TS]
        arrays = [np.array(x) for x in mean_x_s_TS]
        mean_list =[np.mean(k) for k in zip(*arrays)]
    elif ny == 5:
        mean_x_s_TS = [b3d[0].x_s_TS, b3d[1].x_s_TS, b3d[2].x_s_TS, b3d[3].x_s_TS]
        arrays = [np.array(x) for x in mean_x_s_TS]
        mean_list = [np.mean(k) for k in zip(*arrays)]
    place_holder_distance_traveled = []
    for i2 in range(len(mean_list)):
        place_holder_distance_traveled.append(mean_list[i2]-mean_list[0])

    for i3 in range(1,len(place_holder_distance_traveled)):
        if (place_holder_distance_traveled[i3] - place_holder_distance_traveled[i3-1])<0:
            place_holder_distance_traveled[i3] = place_holder_distance_traveled[i3-1]

    place_holder_distance_traveled.insert(0,name_prefix[i])
    DistanceTraveled.append(place_holder_distance_traveled)



np.savetxt('/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/VCR Blue Carbon/Distance Traveled/Metompkin_Bay_Names_VRSLR_TS_Distance_Traveled.csv',
           DistanceTraveled, delimiter=', ', fmt='% s')