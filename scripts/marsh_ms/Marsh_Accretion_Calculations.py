# Load NPZ files and allow you to look at variables

import numpy as np
import matplotlib.pyplot as plt
import os
from sigfig import round

os.chdir("/Users/ceclmac/PycharmProjects/CASCADE/Run_output")
# run_name='Wreck_ACC_RSLR3_S3' # 5 Length
run_name = "Replace 250"  # 4 length
# run_name='Metompkin_Marsh_S10_3'
# run_name='Smith_S10_3' # 5

name_prefix = run_name
nt_run = 250
number_barrier3d_models = 1

# --------- plot ---------
output = np.load(run_name + ".npz", allow_pickle=True)
cascade = output["cascade"]
cascade = cascade[0]
b3d = cascade.barrier3d
ny = np.size(b3d)
bmft = cascade._bmft_coupler
bmft_elev = bmft._bmftc[0].elevation[50:]
BMFT_Landscape = bmft._LandscapeTypeWidth_TS[0]
Islandward_Marsh_Edge = bmft._bmftc[0]._Forest_edge[50:]
Bay_Marsh_Edge = bmft._bmftc[0]._Marsh_edge[50:]
Marsh_Boundary = Islandward_Marsh_Edge - Bay_Marsh_Edge
Shoreline_Change_TS = cascade.barrier3d[0].ShorelineChangeTS
BB_Transect_TS = []
B3D_Transect_TS = []
Cum_Shoreline_Change = []
B3D_Post_Change_TS = bmft._b3d_elev_after_PyBMFT_TS
Post_Change_Transect_TS = []
Shoreface_Movement = b3d[0].x_s_TS

for t in range(nt_run):
        # for t in range(int(41)):

        BB_transect = np.flip(
            cascade._bmft_coupler._bmftc[0].elevation[
                cascade._bmft_coupler._bmftc[0].startyear + t,
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

for i in range(len(B3D_Post_Change_TS)):
    Post_Change_Transect_TS.append(np.mean(B3D_Post_Change_TS[i], axis=1) * 10)

for i in range(nt_run):
    Cum_Shoreline_Change.append(sum(Shoreline_Change_TS[:i]))

Off_TS = []
for k in range(len(B3D_Transect_TS)):
    x1 = np.zeros(int(abs(Cum_Shoreline_Change[k])))
    y1 = B3D_Transect_TS[k]
    l1 = x1.tolist()
    l2 = y1.tolist()
    l1.extend(l2)
    Off_TS.append(l1)

#####
# Interpolate

B3D_Interp_TS = []
Conv_Post_Change_Transect_TS = []
Off_TS_M = []
for l in range(0,len(B3D_Transect_TS)):
    x = np.linspace(
        1, len(B3D_Transect_TS[l]) * 10, num=len(B3D_Transect_TS[l]) * 10
    )
    xp = (
            np.linspace(1, len(B3D_Transect_TS[l]), num=len(B3D_Transect_TS[l]))
            * 10
    )
    xp = xp - 5
    Converted_B3D = np.interp(x, xp, B3D_Transect_TS[l])
    B3D_Interp_TS.append(Converted_B3D)
for l in range(0,len(Post_Change_Transect_TS)):
    x = np.linspace(
        1, len(Post_Change_Transect_TS[l]) * 10, num=len(Post_Change_Transect_TS[l]) * 10
    )
    xp = (
            np.linspace(1, len(Post_Change_Transect_TS[l]), num=len(Post_Change_Transect_TS[l]))
            * 10
    )
    xp = xp - 5
    Converted_B3D = np.interp(x, xp, Post_Change_Transect_TS[l])
    Conv_Post_Change_Transect_TS.append(Converted_B3D)
# Convert offset B3D
for l in range(0,len(Off_TS)):
    x = np.linspace(
        1, len(Off_TS[l]) * 10, num=len(Off_TS[l]) * 10
    )
    xp = (
            np.linspace(1, len(Off_TS[l]), num=len(Off_TS[l]))
            * 10
    )
    xp = xp - 5
    Converted_B3D = np.interp(x, xp, Off_TS[l])
    Off_TS_M.append(Converted_B3D)

Rounded_Offset_TS = []
for i in range(len(bmft._x_s_offset_TS[0])):
    round_temp = round(int(bmft._x_s_offset_TS[0][i])+1)
    Rounded_Offset_TS.append(round_temp)


######
# RSLR Calculations
RSLR = (cascade.barrier3d[0].RSLR[0])*10
RSLR_TS = [0]
for i in range(1,len(cascade.barrier3d[0].RSLR)):
    RSLR_i = (cascade.barrier3d[0].RSLR[i]*i)*10
    RSLR_TS.append(RSLR_i)

max_index_TS = []
for k in range(0,len(BB_Transect_TS)):
    max_index = np.argmax(BB_Transect_TS[k])
    max_index_TS.append(max_index)

x_forest_TS = []
x_marsh_TS = []

for t in range(0,len(BB_Transect_TS)):
    BB_transect = np.flip(
        cascade._bmft_coupler._bmftc[0].elevation[
        cascade._bmft_coupler._bmftc[0].startyear + t - 1,
        int(
            cascade._bmft_coupler._bmftc[0].Marsh_edge[
                cascade._bmft_coupler._bmftc[0].startyear + t
                ]
        ):,
        ]
    )

    x_forest = [cascade._bmft_coupler._bmftc[0].B- int(
        cascade._bmft_coupler._bmftc[0].Forest_edge[
            cascade._bmft_coupler._bmftc[0].startyear + t])]
    x_marsh = [len(BB_transect) - 1]
    x_forest_TS.append(x_forest)
    x_marsh_TS.append(x_marsh)

BMFT_Marsh_TS = []
for k in range(len(BB_Transect_TS)):
    Temp_BMFT_Marsh = BB_Transect_TS[k][int(x_forest_TS[k][0]):int(x_marsh_TS[k][0])]
    BMFT_Marsh_TS.append(Temp_BMFT_Marsh)

B3D_Marsh_TS = []
for l in range(1,len(B3D_Interp_TS)):
    temp_B3D_Marsh = B3D_Interp_TS[l][-int(Marsh_Boundary[l]):]
    B3D_Marsh_TS.append(temp_B3D_Marsh)

Cum_RSLR = (bmft._bmftc[0].msl[bmft._bmftc[0].startyear - 1] + bmft._bmftc[0].amp + (bmft._bmftc[0].RSLRi / 1000)) + RSLR_TS #- RSLR #- 0.008



#################
BMFT_Minus_RSLR = []
for i in range(len(BMFT_Marsh_TS)):
    temp_BMFT_Minus_RSLR = BB_Transect_TS[i] - Cum_RSLR[i] - RSLR
    temp_BMFT_Minus_RSLR[:int(x_forest_TS[i][0])] -=RSLR
    BMFT_Minus_RSLR.append(temp_BMFT_Minus_RSLR)

'''
Rel_RSLR = BB_Transect_TS[0]-Cum_RSLR[1]
start_index= np.argmax(Rel_RSLR)
Forest_Index_TS = x_forest_TS - start_index
Marsh_Index_TS = x_marsh_TS - start_index


BMFT_TS = []
for i in range(len(BB_Transect_TS)):
    temp_BMFT_Minus_RSLR = BB_Transect_TS[i] - Cum_RSLR[i]
    temp_BMFT_Minus_RSLR[:int(Forest_Index_TS[i])] -= RSLR
    BMFT_TS.append(temp_BMFT_Minus_RSLR[np.argmax(temp_BMFT_Minus_RSLR):])
'''

#########
# Test


t1 = BMFT_Minus_RSLR[0][Rounded_Offset_TS[1]:]
t2 = Conv_Post_Change_Transect_TS[1][:len(t1)]
assert_array_almost_equal(x= t1, y=t2, decimal = 3)


t1 = BMFT_Minus_RSLR[150][Rounded_Offset_TS[151]:]
t2 = Conv_Post_Change_Transect_TS[151][:len(t1)]
assert_array_almost_equal(x= t1, y=t2, decimal = 3)

t1 = BMFT_Minus_RSLR[200][Rounded_Offset_TS[201]:]
t2 = Conv_Post_Change_Transect_TS[201][:len(t1)]
assert_array_almost_equal(x= t1, y=t2, decimal = 3)

t1 = BMFT_Minus_RSLR[200][Rounded_Offset_TS[201]:]
t2 = Conv_Post_Change_Transect_TS[201][:len(t1)]
assert_array_almost_equal(x= t1, y=t2, decimal = 3)

t1 = BMFT_Minus_RSLR[240][Rounded_Offset_TS[241]:]
t2 = Conv_Post_Change_Transect_TS[241][:len(t1)]

B3D_Forest_Location_TS = []
for i in range(len(Rounded_Offset_TS)):
    Forest_Temp = int(len(bmft_elev[0]) - Rounded_Offset_TS[i] - bmft._bmftc[0]._Forest_edge[50 + i])
    B3D_Forest_Location_TS.append(Forest_Temp)


#plt.plot(Conv_Post_Change_Transect_TS[0][:B3D_Forest_Location_TS[0]+1], label ='Forest')
plt.plot(Conv_Post_Change_Transect_TS[0][B3D_Forest_Location_TS[0]:], label ='Marsh 0')
plt.plot(Conv_Post_Change_Transect_TS[1][B3D_Forest_Location_TS[1]:], label ='Marsh 1')
#plt.plot(BMFT_TS[0],label='BMFT 0')
plt.legend(loc='upper right')
plt.show()

Marsh_Elev_Change = [list(np.zeros(10000))]
for i in range(len(BB_Transect_TS)-1):
    Dummy_Elev_Change = np.zeros(10000)
    Orginal_Elev = BB_Transect_TS[i]
    New_Elev = BB_Transect_TS[i + 1]
    Dif_Len = len(New_Elev) - len(Orginal_Elev)
    if len(New_Elev) > len(Orginal_Elev):
        add_zeros = np.zeros(Dif_Len)
        Orginal_Elev = np.append(Orginal_Elev, add_zeros)
    elif len(New_Elev) < len(Orginal_Elev):
        Orginal_Elev = Orginal_Elev[:Dif_Len]

    Dif = New_Elev - Orginal_Elev
    Dummy_Elev_Change[x_forest_TS[i][0]:x_marsh_TS[i][0]] = Dif[x_forest_TS[i][0]:x_marsh_TS[i][0]]
    Marsh_Elev_Change.append(list(Dummy_Elev_Change))
Marsh_Elev_Change = np.array(Marsh_Elev_Change)


Total_Marsh_Accretion_TS = [list(Marsh_Elev_Change[0])]
for j in range(1,len(Marsh_Elev_Change)):
    Sum_Marsh = Total_Marsh_Accretion_TS[(j-1)] + Marsh_Elev_Change[j]
    Total_Marsh_Accretion_TS.append(Sum_Marsh)


plt.plot(Total_Marsh_Accretion_TS[0][Rounded_Offset_TS[0]:1350],label='TS 0')
plt.plot(Total_Marsh_Accretion_TS[10][Rounded_Offset_TS[0]:1350],label='TS 10')
plt.plot(Total_Marsh_Accretion_TS[20][Rounded_Offset_TS[0]:1350],label='TS 20')
plt.plot(Total_Marsh_Accretion_TS[40][Rounded_Offset_TS[0]:1350],label='TS 40')
plt.plot(Total_Marsh_Accretion_TS[60][Rounded_Offset_TS[0]:1350],label='TS 60')
plt.plot(Total_Marsh_Accretion_TS[80][Rounded_Offset_TS[0]:1350],label='TS 80')
plt.plot(Total_Marsh_Accretion_TS[100][Rounded_Offset_TS[0]:1350],label='TS 100')
plt.plot(Total_Marsh_Accretion_TS[120][Rounded_Offset_TS[0]:1350],label='TS 120')
plt.plot(Total_Marsh_Accretion_TS[140][Rounded_Offset_TS[0]:1350],label='TS 140')
plt.legend(loc='upper left')
plt.show()

''''''
Dif = BB_Transect_TS[1] - BB_Transect_TS[0]

Orginal_Elev = BB_Transect_TS[i]
New_Elev = BB_Transect_TS[i+1]
Dif_Len = len(New_Elev) - len(Orginal_Elev)
if len(New_Elev) > len(Orginal_Elev):
    add_zeros = np.zeros(Dif_Len)
    Orginal_Elev = np.append(Orginal_Elev,add_zeros)
elif len(New_Elev) < len(Orginal_Elev):
    New_Elev = New_Elev[:Dif_Len]

Dif = New_Elev- Orginal_Elev
Dummy_Elev_Change[x_forest_TS[i][0]:x_marsh_TS[i][0]] = Dif[x_forest_TS[i][0]:x_marsh_TS[i][0]]
