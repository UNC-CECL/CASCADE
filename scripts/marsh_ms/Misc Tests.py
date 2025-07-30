import numpy as np

marsh = '/Volumes/BF_Backup/Chapter_2_Runs/Geom_1_IL_10_S0_Test_Coupling_Marsh.npz'
no_marsh = '/Volumes/BF_Backup/Chapter_2_Runs/Geom_1_IL_10_S0_Test_Coupling_No_Marsh.npz'

output_m = np.load(marsh, allow_pickle=True)
output_nm = np.load(no_marsh, allow_pickle=True)

cascade_m = output_m['cascade'][0]
cascade_nm = output_nm['cascade'][0]

distance_m = cascade_m._brie_coupler._brie.x_s_save
distance_nm = cascade_nm._brie_coupler._brie.x_s_save

dif = np.subtract(distance_nm,distance_m)#,axis=0)

t_d_m = distance_m[0][-1] - distance_m[0][0]
t_d_nm = distance_nm[0][-1] - distance_nm[0][0]

import matplotlib.pyplot as plt

plt.plot(cascade_m.barrier3d[0].StormSeries[:,1])
plt.show()
x = 20
