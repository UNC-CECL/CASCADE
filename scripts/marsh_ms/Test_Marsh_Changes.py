import copy

import matplotlib.pyplot as plt
import numpy as np
import os
import copy

import pandas as pd

os.chdir("E:\\Chapter 2\\")
#save_path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 2\\Cascade_CSV_Outputs\\Hindcasts\\'


Base_Name_List = ['Geom_1',
                  'Geom_2',
                  'Geom_3',
                  'Geom_4',
                  'Geom_5']

run_name = 'Geom_5_IH_10_S1.npz' #'Geom_4_IL_10_S49_New_Sink.npz'
#for geos in range(len(Base_Name_List)):
#    temp_run_name = copy.deepcopy(Base_Name_List[geos]+'_Calibrated_Hindcast_2.npz')
#    run_name.append(copy.deepcopy(temp_run_name))

output = np.load(run_name, allow_pickle=True)
cascade = output["cascade"]
cascade = cascade[0]

bmft = cascade._bmft_coupler
bmftc = bmft._bmftc[0]
O_flux = bmftc.fluxes
Forest_e = bmftc.Forest_edge
Marsh_e = bmftc.Marsh_edge
z = bmftc._organic_dep_autoch# Subtract eroded mass from depositional record
x = bmftc._organic_dep_alloch
z1 = z[0]

xyz = np.sum((z,x),axis=0)


plt.plot(z[51][4500:])
plt.plot(z[71][4500:])
plt.plot(z[91][4500:])
plt.plot(z[111][4500:])

plt.show()

plt.plot(z[51][4500:])
plt.plot(z[52][4500:])
plt.plot(z[53][4500:])
plt.plot(z[54][4500:])

plt.show()

Start_Year = 50
End_Year = 175

Cum_C_Deposits = []

for years in range(Start_Year-1,End_Year):
    temp_sum = np.sum(xyz[Start_Year:years],axis=0)
    Cum_C_Deposits.append(copy.deepcopy(temp_sum))

plt.plot(Cum_C_Deposits[0][0:])
plt.plot(Cum_C_Deposits[25][0:])
plt.plot(Cum_C_Deposits[50][0:])
plt.plot(Cum_C_Deposits[75][0:])
plt.plot(Cum_C_Deposits[100][0:])
plt.plot(Cum_C_Deposits[124][0:])

plt.show()

elev = bmftc.elevation

plt.plot(elev[50][4500:])
plt.axvline(x=Forest_e[50]-4500,linestyle ='dashed',color='blue')
plt.axvline(x=Marsh_e[50]-4500,linestyle ='dashed', color='blue')
#plt.show()

plt.plot(elev[70][4500:],color = 'orange')
plt.axvline(x=Forest_e[70]-4500,linestyle ='dashed',color='orange')
plt.axvline(x=Marsh_e[70]-4500,linestyle ='dashed', color='orange')
#plt.show()

plt.plot(elev[90][4500:],color = 'red')
plt.axvline(x=Forest_e[90]-4500,linestyle ='dashed',color='red')
plt.axvline(x=Marsh_e[90]-4500,linestyle ='dashed', color='red')
#plt.show()

plt.plot(elev[110][4500:],color = 'green')
plt.axvline(x=Forest_e[110]-4500,linestyle ='dashed',color='green')
plt.axvline(x=Marsh_e[110]-4500,linestyle ='dashed', color='green')

plt.plot(elev[130][4500:],color = 'black')
plt.axvline(x=Forest_e[130]-4500,linestyle ='dashed',color='black')
plt.axvline(x=Marsh_e[130]-4500,linestyle ='dashed', color='black')
#plt.show()

plt.plot(elev[150][4500:],color = 'grey')
plt.axvline(x=Forest_e[150]-4500,linestyle ='dashed',color='grey')
plt.axvline(x=Marsh_e[150]-4500,linestyle ='dashed', color='grey')
plt.show()


plt.plot(elev[51][4500:])
plt.plot(elev[101][4500:])
plt.plot(elev[151][4500:])
plt.show()

b3d = cascade._barrier3d[0]


z = 20