# Load NPZ files and allow you to look at variables

import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir("/Users/ceclmac/PycharmProjects/CASCADE/Run_output")
# run_name='Wreck_ACC_RSLR3_S3' # 5 Length
run_name = "c_level_150"  # 4 length
# run_name='Metompkin_Marsh_S10_3'
# run_name='Smith_S10_3' # 5

name_prefix = run_name
nt_run = 150
number_barrier3d_models = 1

# --------- plot ---------
output = np.load(run_name + ".npz", allow_pickle=True)
cascade = output["cascade"]
cascade = cascade[0]
b3d = cascade.barrier3d
ny = np.size(b3d)
bmft = cascade._bmft_coupler
BMFT_Landscape = bmft._LandscapeTypeWidth_TS[0]
Islandward_Marsh_Edge = bmft._bmftc[0]._Forest_edge[49:]
Bay_Marsh_Edge = bmft._bmftc[0]._Marsh_edge[49:]
Shoreline_Change_TS = cascade.barrier3d[0].ShorelineChangeTS
BB_Transect_TS = []
B3D_Transect_TS = []
Cum_Shoreline_Change = []


for t in range(nt_run - 1):
        # for t in range(int(41)):

        BB_transect = np.flip(
            cascade._bmft_coupler._bmftc[0].elevation[
                cascade._bmft_coupler._bmftc[0].startyear + t - 1,
                int(
                    cascade._bmft_coupler._bmftc[0].Marsh_edge[
                        cascade._bmft_coupler._bmftc[0].startyear + t
                    ]
                ) :,
            ]
        )
        BB_Transect_TS.append(BB_transect)

for i in range(nt_run):
    B3D_Transect_TS.append(np.mean(cascade.barrier3d[0].DomainTS[i], axis=1) * 10)

for i in range(1,150):
    Cum_Shoreline_Change.append(sum(Shoreline_Change_TS[:i]))

Off_TS = []
for k in range(0,149):
    print(k)
    x1 = np.zeros(int(abs(Cum_Shoreline_Change[k])))
    y1 = B3D_Transect_TS[k]
    l1 = x1.tolist()
    l2 = y1.tolist()
    l1.extend(l2)
    Off_TS.append(l1)

plt.plot(Off_TS[0],label =0)
plt.plot(Off_TS[25],label =25)
plt.plot(Off_TS[50],label =50)
plt.plot(Off_TS[75],label =75)
plt.plot(Off_TS[100],label =100)
plt.plot(Off_TS[125],label =125)
plt.legend(loc="upper right")
plt.show()

plt.plot(Off_TS[0],label =0)
plt.plot(Off_TS[1],label =1)
plt.plot(Off_TS[2],label =2)
plt.plot(Off_TS[3],label =3)
plt.plot(Off_TS[4],label =4)
plt.plot(Off_TS[5],label =5)
plt.legend(loc="upper right")
plt.show()

TS_3 = np.array(Off_TS[3])
TS_2 = np.array(Off_TS[2])

TS = TS_3-TS_2

Off_TS_5 = [x1+y1]
Off_TS_5.append(np.zeros(int(abs(Cum_Shoreline_Change[5]))))
Off_TS_5.append(B3D_Transect_TS[5])

plt.plot(B3D_Transect_TS[147])
plt.show()

plt.plot(BB_Transect_TS[0], label ='0')
plt.plot(BB_Transect_TS[25], label ='25')
plt.plot(BB_Transect_TS[50], label ='50')
plt.plot(BB_Transect_TS[75], label = '75')
plt.plot(BB_Transect_TS[100], label = '100')
plt.plot(BB_Transect_TS[125], label = '125')
plt.plot(BB_Transect_TS[149], label = '147')
plt.legend(loc="upper right")
plt.show()

plt.plot(BB_Transect_TS[1], label ='0')
plt.show()

plt.plot(B3D_Transect_TS[1])
plt.show()


plt.plot(BB_Transect_TS[25], label ='0')
plt.show()

plt.plot(B3D_Transect_TS[25])
plt.show()

plt.plot(BB_Transect_TS[50], label ='0')
plt.show()

plt.plot(B3D_Transect_TS[50])
plt.show()

Marsh_1 = 500
Corrected_Location_1 = abs(Islandward_Marsh_Edge[0]-Bay_Marsh_Edge[0])



Marsh_2 = 500 - 8
Corrected_Location_2 = abs(5500 - Islandward_Marsh_Edge[1])


Marsh_100

Corrected_Location_100 = abs(5500 - Islandward_Marsh_Edge[100])

# For 5

Corrected_Location_5 = abs(Islandward_Marsh_Edge[5]-Bay_Marsh_Edge[5])
Location_in_B3D = round(Corrected_Location_5/10)

Len_TS = len(B3D_Transect_TS[5])
Xs = list(range(0,Len_TS,1))

plt.plot(Xs[:-Location_in_B3D],B3D_Transect_TS[5][:-Location_in_B3D])
plt.plot(Xs[-Location_in_B3D:],B3D_Transect_TS[5][-Location_in_B3D:])
plt.show()

# For 100

Corrected_Location_100 = abs(Islandward_Marsh_Edge[100]-Bay_Marsh_Edge[100])
Location_in_B3D = round(Corrected_Location_100/10)

Len_TS = len(B3D_Transect_TS[100])
Xs = list(range(0,Len_TS,1))

plt.plot(Xs[:-Location_in_B3D],B3D_Transect_TS[100][:-Location_in_B3D])
plt.plot(Xs[-Location_in_B3D:],B3D_Transect_TS[100][-Location_in_B3D:])
plt.show()

# For 25

Corrected_Location_25 = abs(Islandward_Marsh_Edge[25]-Bay_Marsh_Edge[25])
Location_in_B3D = round(Corrected_Location_25/10)

Len_TS = len(B3D_Transect_TS[25])
Xs = list(range(0,Len_TS,1))

plt.plot(Xs[:-Location_in_B3D],B3D_Transect_TS[25][:-Location_in_B3D])
plt.plot(Xs[-Location_in_B3D:],B3D_Transect_TS[25][-Location_in_B3D:])
plt.show()

# For 50

Corrected_Location_50 = abs(Islandward_Marsh_Edge[50]-Bay_Marsh_Edge[50])
Location_in_B3D = round(Corrected_Location_50/10)

Len_TS = len(B3D_Transect_TS[50])
Xs = list(range(0,Len_TS,1))

plt.plot(Xs[:-Location_in_B3D],B3D_Transect_TS[50][:-Location_in_B3D])
plt.plot(Xs[-Location_in_B3D:],B3D_Transect_TS[50][-Location_in_B3D:])
plt.show()

# For 75

Corrected_Location_75 = abs(Islandward_Marsh_Edge[75]-Bay_Marsh_Edge[75])
Location_in_B3D = round(Corrected_Location_75/10)

Len_TS = len(B3D_Transect_TS[75])
Xs = list(range(0,Len_TS,1))

plt.plot(Xs[:-Location_in_B3D],B3D_Transect_TS[75][:-Location_in_B3D])
plt.plot(Xs[-Location_in_B3D:],B3D_Transect_TS[75][-Location_in_B3D:])
plt.show()

# For 80

Corrected_Location_80 = abs(Islandward_Marsh_Edge[80]-Bay_Marsh_Edge[80])
Location_in_B3D = round(Corrected_Location_80/10)

Len_TS = len(B3D_Transect_TS[80])
Xs = list(range(0,Len_TS,1))

plt.plot(Xs[:-Location_in_B3D],B3D_Transect_TS[80][:-Location_in_B3D])
plt.plot(Xs[-Location_in_B3D:],B3D_Transect_TS[80][-Location_in_B3D:])
plt.show()

# For 85

Corrected_Location_85 = abs(Islandward_Marsh_Edge[85]-Bay_Marsh_Edge[85])
Location_in_B3D = round(Corrected_Location_85/10)

Len_TS = len(B3D_Transect_TS[85])
Xs = list(range(0,Len_TS,1))

plt.plot(Xs[:-Location_in_B3D],B3D_Transect_TS[85][:-Location_in_B3D])
plt.plot(Xs[-Location_in_B3D:],B3D_Transect_TS[85][-Location_in_B3D:])
plt.show()

# For 95

Corrected_Location_95 = abs(Islandward_Marsh_Edge[95]-Bay_Marsh_Edge[95])
Location_in_B3D = round(Corrected_Location_95/10)

Len_TS = len(B3D_Transect_TS[95])
Xs = list(range(0,Len_TS,1))

plt.plot(Xs[:-Location_in_B3D],B3D_Transect_TS[95][:-Location_in_B3D])
plt.plot(Xs[-Location_in_B3D:],B3D_Transect_TS[95][-Location_in_B3D:])
plt.show()

