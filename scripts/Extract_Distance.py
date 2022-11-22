import numpy as np
import os

os.chdir("/Users/ceclmac/PycharmProjects/CASCADE/Run_output")


Smith_Run_Name = ['Smith_Marsh_Low_S_1','Smith_Marsh_Low_S_2','Smith_Marsh_Low_S_3',
                  'Smith_S10_1','Smith_S10_2','Smith_S10_3',
                  'Smith_Marsh_S25_1','Smith_Marsh_S25_2','Smith_Marsh_S25_3']
Metompkin_Marsh_Run_Name = ['Metompkin_Marsh_On_Low_S_1','Metompkin_Marsh_On_Low_S_2','Metompkin_Marsh_On_Low_S_3',
                            'Metompkin_Marsh_S10_1','Metompkin_Marsh_S10_2','Metompkin_Marsh_S10_3',
                            'Metompkin_Marsh_On_Low_S25_1','Metompkin_Marsh_On_Low_S25_2','Metompkin_Marsh_On_Low_S25_3'
                            ]
Metompkin_Bay_Run_Name = ['Metompkin_Bay_Low_S_1','Metompkin_Bay_Low_S_2','Metompkin_Bay_Low_S_3',
                            'Metompkin_No_Marsh_S10_1','Metompkin_No_Marsh_S10_2','Metompkin_No_Marsh_S10_3',
                            'Metompkin_No_Marsh_S25_1','Metompkin_No_Marsh_S25_2','Metompkin_No_Marsh_S25_3'
                            ]
Wreck_Run_Name = ['Wreck_Low_S_1','Wreck_Low_S_2','Wreck_Low_S_3',
                  'Wreck_S10_1','Wreck_S10_2','Wreck_S10_3',
                  'Wreck_Low_S25_1','Wreck_Low_S25_2','Wreck_Low_S25_3'
                  ]
name_prefix = Metompkin_Bay_Run_Name
DistanceTraveled = []
Toe_To_Shoreline = []

for i in range(len(name_prefix)):
    output = np.load(name_prefix[i] + ".npz", allow_pickle=True)
    cascade = output["cascade"]
    cascade = cascade[0]
    b3d = cascade.barrier3d
    ny = np.size(b3d)
    place_holder_distance_traveled = [name_prefix[i]]
    place_holder_toe_to_shoreline = [name_prefix[i]]
    for i2 in range(len(b3d)):
        place_holder_distance_traveled.append(b3d[i2].x_s_TS[-1] - b3d[i2].x_s_TS[0])
        place_holder_toe_to_shoreline.append(-b3d[i2].x_s_TS[-1] + b3d[i2].x_b_TS[-1])
    DistanceTraveled.append(place_holder_distance_traveled)
    Toe_To_Shoreline.append(place_holder_toe_to_shoreline)

np.savetxt('/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/VCR Blue Carbon/Distance Traveled/Metompkin_Bay_Distance_Traveled.csv',
           DistanceTraveled, delimiter=', ', fmt='% s')
np.savetxt('/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/VCR Blue Carbon/Distance Traveled/Metompkin_Bay_Toe_Shoreline.csv',
           Toe_To_Shoreline, delimiter=', ', fmt='% s')