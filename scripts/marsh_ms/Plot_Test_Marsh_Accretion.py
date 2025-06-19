import copy

import matplotlib.pyplot as plt
import numpy as np
import os
import copy

import pandas as pd

os.chdir("E:\\Chapter 2\\")
#save_path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 2\\Cascade_CSV_Outputs\\Hindcasts\\'




Base_Name_List = ['Geom_1',
                  'Geom_3',
                  'Geom_4',
                  'Geom_5']

Base_Name = Base_Name_List[0]

if Base_Name == 'Geom_1':
    Cascade_Offset = 90
if Base_Name == 'Geom_3':
    Cascade_Offset = 490
if Base_Name == 'Geom_4':
    Cascade_Offset = 170
if Base_Name == 'Geom_5':
    Cascade_Offset = 160

run_name = 'Geom_1_IH_10_S1.npz' #'Geom_4_IL_10_S49_New_Sink.npz'
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

initial_C = bmftc._marshOM_initial # kg's C accross entire platfrom
initial_marsh_C_g = initial_C*1000

# Create initial C layer
Initial_Marsh_Width = Forest_e[49] - Marsh_e[49]
Per_m_Marsh_C_g = initial_marsh_C_g/Initial_Marsh_Width
Initial_Transect_C_Value = np.zeros(len(bmftc.elevation[0]))

# Add initial C to marsh cells
Initial_Transect_C_Value[int(Marsh_e[49]):int(Forest_e[49])] = Per_m_Marsh_C_g
Total_C_Deposited_TS = [Initial_Transect_C_Value]

#Bmftc._BayOM[yr]

Shoreline_location_TS = cascade.brie.x_s_save
All_Reletive_Shoreline = ((-Shoreline_location_TS+1624+90)+Forest_e[49])[0]


# Start Year = [50]
Start_Year = 49
End_Year = Start_Year+125

C_autoch = bmftc._organic_dep_autoch# Subtract eroded mass from depositional record
C_alloch = bmftc._organic_dep_alloch

Total_Annual_C_Change = np.sum((C_alloch,C_autoch),axis=0)

for years in range(Start_Year,End_Year):
    temp_sum = np.sum(Total_Annual_C_Change[Start_Year:years],axis=0)
    marsh_accumulation = temp_sum[int(Marsh_e[years]):]
    temp_total_marsh_C_Change = np.zeros(len(temp_sum))
    temp_total_marsh_C_Change[int(Marsh_e[years]):] = marsh_accumulation
    new_C_deposition_plus_C0 = np.sum((temp_total_marsh_C_Change,Initial_Transect_C_Value),axis=0)
    new_C_deposition_plus_C0[new_C_deposition_plus_C0 < 0] = 0
    Total_C_Deposited_TS.append(copy.deepcopy(new_C_deposition_plus_C0))

plt.plot(Total_C_Deposited_TS[0][0:])
#plt.plot(Total_C_Deposited_TS[10][0:])
plt.plot(Total_C_Deposited_TS[24][0:],alpha=0.8)
#plt.plot(Total_C_Deposited_TS[30][0:])
plt.plot(Total_C_Deposited_TS[49][0:], alpha=0.8)
plt.plot(Total_C_Deposited_TS[74][0:], alpha=0.8)
plt.plot(Total_C_Deposited_TS[99][0:], alpha=0.8)
plt.plot(Total_C_Deposited_TS[124][0:],alpha=0.8)

# plt.axvline(x=All_Reletive_Shoreline[0])
# plt.axvline(x=Marsh_e[49])
# #plt.axvline(x=All_Reletive_Shoreline[24])
# plt.axvline(x=Marsh_e[74])
# #plt.axvline(x=All_Reletive_Shoreline[49])
# plt.axvline(x=Marsh_e[99])
# plt.axvline(x=All_Reletive_Shoreline[74])
# #plt.axvline(x=Marsh_e[124])
# plt.axvline(x=All_Reletive_Shoreline[99])
# #plt.axvline(x=Marsh_e[144],color = 'grey')
# plt.axvline(x=All_Reletive_Shoreline[124])
# #plt.axvline(x=Marsh_e[172],color='black')

plt.xlim(4500,7500)
#plt.ylim(0,3000)
plt.show()

plt.plot(Total_Annual_C_Change[25][0:])
plt.plot(Total_Annual_C_Change[50][0:])
plt.plot(Total_Annual_C_Change[75][0:])
plt.plot(Total_Annual_C_Change[100][0:])
plt.plot(Total_Annual_C_Change[124][0:])
plt.xlim(4000,7000)
plt.ylim(0,3000)
plt.show()

plt.plot(bmftc.elevation[50])
plt.plot(bmftc.elevation[55])
plt.plot(bmftc.elevation[60])
plt.plot(bmftc.elevation[70])
plt.plot(bmftc.elevation[80])

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