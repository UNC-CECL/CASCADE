# import copy
#
# import matplotlib.pyplot as plt
# import numpy as np
# import os
# import copy
#
# import pandas as pd
#
# os.chdir("E:\\Chapter 2\\")
# #save_path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 2\\Cascade_CSV_Outputs\\Hindcasts\\'
#
#
#
#
# Base_Name_List = ['Geom_1',
#                   'Geom_3',
#                   'Geom_4',
#                   'Geom_5']
#
# Base_Name = Base_Name_List[1]
#
# if Base_Name == 'Geom_1':
#     Cascade_Offset = 110
# if Base_Name == 'Geom_3':
#     Cascade_Offset = 490
# if Base_Name == 'Geom_4':
#     Cascade_Offset = 170
# if Base_Name == 'Geom_5':
#     Cascade_Offset = 160
#
# run_name = 'Geom_3_IH_10_S1.npz' #'Geom_4_IL_10_S49_New_Sink.npz'
# #for geos in range(len(Base_Name_List)):
# #    temp_run_name = copy.deepcopy(Base_Name_List[geos]+'_Calibrated_Hindcast_2.npz')
# #    run_name.append(copy.deepcopy(temp_run_name))
#
# output = np.load(run_name, allow_pickle=True)
# cascade = output["cascade"]
# cascade = cascade[0]
#
# bmft = cascade._bmft_coupler
# bmftc = bmft._bmftc[0]
# O_flux = bmftc.fluxes
# Forest_e = bmftc.Forest_edge
# Marsh_e = bmftc.Marsh_edge
#
# initial_C = bmftc._marshOM_initial # kg's C accross entire platfrom
# initial_marsh_C_g = initial_C*1000
#
# # Create initial C layer
# Initial_Marsh_Width = Forest_e[49] - Marsh_e[49]
# Per_m_Marsh_C_g = initial_marsh_C_g/Initial_Marsh_Width
# Initial_Transect_C_Value = np.zeros(len(bmftc.elevation[0]))
#
# # Add initial C to marsh cells
# Initial_Transect_C_Value[int(Marsh_e[49]):int(Forest_e[49])] = Per_m_Marsh_C_g
# Total_C_Deposited_TS = [Initial_Transect_C_Value]
# Total_Bay_C_Deposit = []
# Combined_Test = []
#
# #Bmftc._BayOM[yr]
#
# Shoreline_location_TS = cascade.brie.x_s_save
# All_Reletive_Shoreline = ((-Shoreline_location_TS+1624+Cascade_Offset)+Forest_e[49])[0]
#
#
# # Start Year = [50]
# Start_Year = 50
# End_Year = Start_Year+123
#
# C_autoch = bmftc._organic_dep_autoch# Subtract eroded mass from depositional record
# C_alloch = bmftc._organic_dep_alloch
#
# Total_Annual_C_Change = np.sum((C_alloch,C_autoch),axis=0)
#
# for years in range(Start_Year,End_Year):
#     temp_sum = np.sum(Total_Annual_C_Change[Start_Year:years+1],axis=0)
#     marsh_accumulation = temp_sum[int(Marsh_e[years]):]
#     temp_total_marsh_C_Change = np.zeros(len(temp_sum))
#     temp_total_marsh_C_Change[int(Marsh_e[years]):] = marsh_accumulation
#     bay_accumulation = temp_sum[:int(Marsh_e[years])]
#     temp_total_bay_C_change = np.zeros(len(temp_sum))
#     temp_total_bay_C_change[:int(Marsh_e[years])] = bay_accumulation
#
#
#     new_C_deposition_plus_C0 = np.sum((temp_total_marsh_C_Change,Initial_Transect_C_Value),axis=0)
#     new_C_deposition_plus_C0[new_C_deposition_plus_C0 < 0] = 0
#     temp_total_bay_C_change[temp_total_bay_C_change < 0] = 0
#
#     Total_C_Deposited_TS.append(copy.deepcopy(new_C_deposition_plus_C0))
#     Total_Bay_C_Deposit.append(copy.deepcopy(temp_total_bay_C_change))
#     Combined_Test.append(copy.deepcopy(np.sum((new_C_deposition_plus_C0,temp_total_bay_C_change),axis=0)))
# #
# plt.plot(Total_C_Deposited_TS[0][0:])
# #plt.plot(Total_C_Deposited_TS[10][0:])
# plt.plot(Total_C_Deposited_TS[24][0:],alpha=0.8)
# #plt.plot(Total_C_Deposited_TS[30][0:])
# plt.plot(Total_C_Deposited_TS[49][0:], alpha=0.8)
# plt.plot(Total_C_Deposited_TS[74][0:], alpha=0.8)
# plt.plot(Total_C_Deposited_TS[99][0:], alpha=0.8)
# plt.plot(Total_C_Deposited_TS[123][0:],alpha=0.8)
# plt.axvline(Marsh_e[50+123])
# plt.axvline(Marsh_e[49])
# plt.axvline(Forest_e[49],linestyle='dashed')
# plt.axvline(Forest_e[74],linestyle='dashed')
# plt.axvline(Forest_e[99],linestyle='dashed')
# plt.axvline(Forest_e[124],linestyle='dashed')
# plt.axvline(Forest_e[149],linestyle='dashed')
# plt.axvline(Forest_e[173],linestyle='dashed')
# plt.axvline(All_Reletive_Shoreline[0],color='black')
# plt.axvline(All_Reletive_Shoreline[124],color='black')
# plt.xlim(4200,7500)
# plt.show()
# x = 20
#
# # plt.axvline(x=All_Reletive_Shoreline[0])
# # plt.axvline(x=Marsh_e[49])
# # #plt.axvline(x=All_Reletive_Shoreline[24])
# # plt.axvline(x=Marsh_e[74])
# # #plt.axvline(x=All_Reletive_Shoreline[49])
# # plt.axvline(x=Marsh_e[99])
# # plt.axvline(x=All_Reletive_Shoreline[74])
# # #plt.axvline(x=Marsh_e[124])
# # plt.axvline(x=All_Reletive_Shoreline[99])
# # #plt.axvline(x=Marsh_e[144],color = 'grey')
# # plt.axvline(x=All_Reletive_Shoreline[124])
# # #plt.axvline(x=Marsh_e[172],color='black')
#
# plt.xlim(4500,7500)
# #plt.ylim(0,3000)
# plt.show()
#
#
# plt.plot(Total_Bay_C_Deposit[0])
# plt.plot(Total_Bay_C_Deposit[24])
# plt.plot(Total_Bay_C_Deposit[49])
# plt.plot(Total_Bay_C_Deposit[74])
# plt.plot(Total_Bay_C_Deposit[99])
# #plt.plot(Total_Bay_C_Deposit[124])
# plt.xlim(4500,6000)
# plt.show()
#
# #plt.plot(Total_Bay_C_Deposit[100])
# #plt.plot(Total_C_Deposited_TS[100][0:], alpha=0.8)
# plt.plot(Total_Bay_C_Deposit[20])
# plt.plot(Total_C_Deposited_TS[21][0:], alpha=0.8)
# plt.plot(Combined_Test[20][0:],color='black')
# #plt.plot(Total_Bay_C_Deposit[90])
# #plt.plot(Total_C_Deposited_TS[90][0:], alpha=0.8)
# #plt.plot(Total_Bay_C_Deposit[99])
# #plt.plot(Total_C_Deposited_TS[99][0:], alpha=0.8)
# plt.xlim(4800,5200)
# plt.show()
#
# plt.plot(Total_Bay_C_Deposit[74],color = 'blue', alpha=0.2)
# plt.plot(Total_C_Deposited_TS[75][0:], alpha=0.2,color='red')
# plt.plot(Combined_Test[74][0:],color='black', alpha=0.2,linestyle='dashed')
# plt.xlim(4800,5200)
#
# plt.show()
#
# plt.plot(Total_Bay_C_Deposit[99],color = 'blue', alpha=0.2)
# plt.plot(Total_C_Deposited_TS[100][0:], alpha=0.2,color='red')
# plt.plot(Combined_Test[99][0:],color='black', alpha=0.2,linestyle='dashed')
# plt.xlim(4800,5200)
#
# plt.show()
#
# plt.plot(Total_Bay_C_Deposit[80],color = 'blue', alpha=0.2)
# plt.plot(Total_C_Deposited_TS[81][0:], alpha=0.2,color='red')
# plt.plot(Combined_Test[80][0:],color='black', alpha=0.2,linestyle='dashed')
# plt.axvline(x=Marsh_e[80+50])
# plt.xlim(4400,5200)
# plt.show()
#
# plt.plot(Total_Bay_C_Deposit[85],color = 'blue', alpha=0.2)
# plt.plot(Total_C_Deposited_TS[86][0:], alpha=0.2,color='red')
# plt.plot(Combined_Test[85][0:],color='black', alpha=0.2,linestyle='dashed')
# plt.axvline(x=Marsh_e[85+50])
# plt.xlim(4400,5200)
# plt.show()
#
#
# plt.plot(Total_Bay_C_Deposit[90],color = 'blue', alpha=0.2)
# plt.plot(Total_C_Deposited_TS[91][0:], alpha=0.2,color='red')
# plt.plot(Combined_Test[90][0:],color='black', alpha=0.2,linestyle='dashed')
# plt.axvline(x=Marsh_e[89+50],linestyle='dashed')
# plt.axvline(x=Marsh_e[90+50])
# plt.axvline(x=Marsh_e[91+50],linestyle='dashed')
#
# plt.xlim(4400,5200)
#
# plt.show()
#
# plt.plot(Total_Bay_C_Deposit[90],color = 'red', alpha=0.2)
# plt.plot(Total_Bay_C_Deposit[91],color = 'blue', alpha=0.2)
# plt.show()
#
# plt.plot(Total_C_Deposited_TS[92][0:], alpha=0.2,color='red')
# plt.plot(Combined_Test[91][0:],color='black', alpha=0.2,linestyle='dashed')
# plt.axvline(x=Marsh_e[90+50],linestyle='dashed')
# plt.axvline(x=Marsh_e[91+50])
# plt.axvline(x=Marsh_e[92+50],linestyle='dashed')
# plt.xlim(4400,5200)
# plt.show()
#
#
#
# plt.plot(Total_Bay_C_Deposit[91],color = 'blue', alpha=0.2)
# plt.plot(Total_C_Deposited_TS[92][0:], alpha=0.2,color='red')
# plt.plot(Combined_Test[91][0:],color='black', alpha=0.2,linestyle='dashed')
# plt.axvline(x=Marsh_e[91+50])
# plt.xlim(4400,5200)
# plt.show()
#
# plt.plot(Total_Bay_C_Deposit[93],color = 'blue', alpha=0.2)
# plt.plot(Total_C_Deposited_TS[94][0:], alpha=0.2,color='red')
# plt.plot(Combined_Test[93][0:],color='black', alpha=0.2,linestyle='dashed')
# plt.axvline(x=Marsh_e[93+50])
# plt.axvline(x=Marsh_e[94+50])
#
# plt.xlim(4400,5200)
#
# plt.show()
#
#
# plt.plot(Combined_Test[90][0:],color='black', alpha=0.2,linestyle='dashed')
# plt.plot(Combined_Test[100][0:],color='orange', alpha=0.2,linestyle='dashed')
# plt.plot(Combined_Test[110][0:],color='blue', alpha=0.2,linestyle='dashed')
# plt.plot(Total_C_Deposited_TS[111],color='red',alpha=0.4)
# plt.xlim(4400,5200)
#
# plt.show()
#
#
# plt.plot(Total_Annual_C_Change[25][0:])
# plt.plot(Total_Annual_C_Change[50][0:])
# plt.plot(Total_Annual_C_Change[75][0:])
# plt.plot(Total_Annual_C_Change[100][0:])
# plt.plot(Total_Annual_C_Change[124][0:])
# plt.xlim(4000,7000)
# plt.ylim(0,3000)
# plt.show()
#
# plt.plot(bmftc.elevation[50])
# plt.plot(bmftc.elevation[55])
# plt.plot(bmftc.elevation[60])
# plt.plot(bmftc.elevation[70])
# plt.plot(bmftc.elevation[80])
#
# plt.show()
#
# elev = bmftc.elevation
#
# plt.plot(elev[50])
# plt.axvline(x=Forest_e[50],linestyle ='dashed',color='blue')
# plt.axvline(x=Marsh_e[50],linestyle ='dashed', color='blue')
# #plt.show()
#
# plt.plot(elev[70],color = 'orange')
# plt.axvline(x=Forest_e[70],linestyle ='dashed',color='orange')
# plt.axvline(x=Marsh_e[70],linestyle ='dashed', color='orange')
# #plt.show()
#
# plt.plot(elev[90],color = 'red')
# plt.axvline(x=Forest_e[90],linestyle ='dashed',color='red')
# plt.axvline(x=Marsh_e[90],linestyle ='dashed', color='red')
# #plt.show()
#
# plt.plot(elev[110],color = 'green')
# plt.axvline(x=Forest_e[110],linestyle ='dashed',color='green')
# plt.axvline(x=Marsh_e[110],linestyle ='dashed', color='green')
#
# plt.plot(elev[130],color = 'black')
# plt.axvline(x=Forest_e[130],linestyle ='dashed',color='black')
# plt.axvline(x=Marsh_e[130],linestyle ='dashed', color='black')
# #plt.show()
#
# plt.plot(elev[150],color = 'grey')
# plt.axvline(x=Forest_e[150],linestyle ='dashed',color='grey')
# plt.axvline(x=Marsh_e[150],linestyle ='dashed', color='grey')
#
# plt.plot(elev[173],color = 'grey')
# plt.axvline(x=Forest_e[173],linestyle ='dashed',color='grey')
# plt.axvline(x=Marsh_e[173],linestyle ='dashed', color='grey')
#
# plt.xlim(4500,7500)
# plt.show()

#
# plt.plot(elev[51][4500:])
# plt.plot(elev[101][4500:])
# plt.plot(elev[151][4500:])
# plt.show()
#
# b3d = cascade._barrier3d[0]
#
#
# z = 20
###############################

import copy

import matplotlib.pyplot as plt
import numpy as np
import os
import copy

import pandas as pd

#os.chdir("G:\\Chapter_2_Runs\\")
os.chdir("E:\\Chapter 2\\")

#save_path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 2\\Cascade_CSV_Outputs\\Hindcasts\\'

Base_Name_List = ['Geom_1',
                  'Geom_3',
                  'Geom_4',
                  'Geom_5']

Base_Name_List = ['Geom_3']

Base_Name = Base_Name_List[0]

if Base_Name == 'Geom_1':
    Cascade_Offset = 120
if Base_Name == 'Geom_3':
    Cascade_Offset = 490
if Base_Name == 'Geom_4':
    Cascade_Offset = 170
if Base_Name == 'Geom_5':
    Cascade_Offset = 160

run_name = 'Geom_3_I_5_S9_N_RSLR_60.npz' #'Geom_4_IL_10_S49_New_Sink.npz'
#for geos in range(len(Base_Name_List)):
#    temp_run_name = copy.deepcopy(Base_Name_List[geos]+'_Calibrated_Hindcast_2.npz')
#    run_name.append(copy.deepcopy(temp_run_name))

output = np.load(run_name, allow_pickle=True)
cascade = output["cascade"]
cascade = cascade[0]
#
bmft = cascade._bmft_coupler
bmftc = bmft._bmftc[0]
O_flux = bmftc.fluxes
Forest_e = bmftc.Forest_edge
Marsh_e = bmftc.Marsh_edge
b3d = cascade.barrier3d[0]
DomainTS = b3d.DomainTS
cross_island_TS = []
dune_TS = []
combined_TS = []
combined_TS_NAVD88 = []
Shoreline_Change_TS = []
Total_SLR = []
for t in range (1,len(b3d.ShorelineChangeTS)+1):
    Temp_Shoreline_Change = int(abs(np.sum(b3d.ShorelineChangeTS[0:t])))
    Shoreline_Change_TS.append(copy.deepcopy(Temp_Shoreline_Change))
    Temp_SLR = (np.sum(b3d.RSLR[0:t])*10)+.46
    Total_SLR.append(copy.deepcopy(Temp_SLR))


for years in range(len(DomainTS)):
    Temp_Value = np.mean(DomainTS[years],axis=1)
    cross_island_TS.append(copy.deepcopy(Temp_Value))
    #Temp_Shoreline_Change = int(abs(np.sum(b3d.ShorelineChangeTS[years+1])))
    Shoreline_Buffer = np.zeros(Shoreline_Change_TS[years])
    Temp_Dunes =np.mean(b3d.DuneDomain[years],axis=0)+b3d.BermEl
    dune_TS.append(copy.deepcopy(Temp_Dunes))
    combined_temp_transect = np.concatenate((Shoreline_Buffer,Temp_Dunes,Temp_Value),axis=0)
    combined_TS.append(copy.deepcopy(combined_temp_transect*10))
    combined_TS_NAVD88.append(copy.deepcopy((combined_temp_transect*10)+Total_SLR[years]))

plt.plot(combined_TS[1],color='red')
plt.plot(combined_TS[25])
plt.plot(combined_TS[49])
plt.plot(combined_TS[74],color='blue')
plt.plot(combined_TS[99])
plt.plot(combined_TS[124],color='orange')
plt.show()

fig, ax = plt.subplots()

#ax.plot(combined_TS_NAVD88[0],color='black', alpha=0.8)
ax.plot(combined_TS_NAVD88[1],color='black', alpha=0.8)
ax.plot(combined_TS_NAVD88[24],linestyle='dotted', alpha=0.6, color='#4477AA')
ax.plot(combined_TS_NAVD88[49], linestyle='dashdot', alpha=0.6, color='#EE6677')
ax.plot(combined_TS_NAVD88[74],linestyle='dotted', alpha=0.6, color='#66CCEE')
ax.plot(combined_TS_NAVD88[99], linestyle='dashdot', alpha=0.6, color='#CCBB44')
ax.plot(combined_TS_NAVD88[124],color='#228833', alpha=0.8)
ax.axvline(x=0)
#ax.scatter(x=100,y=0.46)
#ax.scatter(x=50,y=0.46)
plt.xlim(-2,130)

plt.ylim(-1,3.25)
ax.invert_xaxis()
plt.savefig(
    ('C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 2\\Figures\\Revised_RSLR_Figs\\I_5S_S9_Elevation.pdf'),
    format='pdf')
plt.show()


plt.plot(bmftc.elevation[49])
plt.plot(bmftc.elevation[73])
plt.plot(bmftc.elevation[98])
plt.plot(bmftc.elevation[123])
plt.plot(bmftc.elevation[148])
plt.plot(bmftc.elevation[173])
plt.ylim(0,3.5)
plt.xlim(4500,6000)
plt.show()

initial_C = bmftc._marshOM_initial # kg's C accross entire platfrom
initial_marsh_C_g = initial_C*1000

# Create initial C layer
Initial_Marsh_Width = Forest_e[49] - Marsh_e[49]
All_Marsh_Width = np.subtract(Forest_e,Marsh_e)
Per_m_Marsh_C_g = initial_marsh_C_g/Initial_Marsh_Width
Initial_Transect_C_Value = np.zeros(len(bmftc.elevation[0]))

# Add initial C to marsh cells
Initial_Transect_C_Value[int(Marsh_e[49]):int(Forest_e[49])] = Per_m_Marsh_C_g
Total_C_Deposited_TS = [Initial_Transect_C_Value]
Bay_C_Deposit_TS = [np.zeros(len(Initial_Transect_C_Value))]


Shoreline_location_TS = cascade.brie.x_s_save
All_Reletive_Shoreline = ((-Shoreline_location_TS+1624+Cascade_Offset)+Forest_e[49])[0]


# Start Year = [50]
Start_Year = 49
End_Year = Start_Year+125

C_autoch = bmftc._organic_dep_autoch# Subtract eroded mass from depositional record
C_alloch = bmftc._organic_dep_alloch

Total_Annual_C_Change = np.sum((C_alloch,C_autoch),axis=0)
All_C_Deposited = []
for years in range(Start_Year,End_Year):
    Temp_C = np.sum(Total_Annual_C_Change[Start_Year:years],axis=0)
    All_C_Deposited.append(copy.deepcopy(Temp_C))
    bay_accumulation = np.sum(Total_Annual_C_Change[years:years+1],axis=0)[:int(Marsh_e[years])]
    temp_total_bay_C_change = np.zeros(len(Total_Annual_C_Change[0]))
    temp_total_bay_C_change[:int(Marsh_e[years])] = bay_accumulation
    temp_total_bay_C_change = np.sum((temp_total_bay_C_change,Bay_C_Deposit_TS[years-Start_Year]),axis=0)
    Bay_C_Deposit_TS.append(copy.deepcopy(temp_total_bay_C_change))

Bay_C_Deposit_TS_Non_Zero = []
for bay_years in range(len(Bay_C_Deposit_TS)):
    Temp_Bay = Bay_C_Deposit_TS[bay_years]
    Temp_Bay[Temp_Bay<0] = 0
    Bay_C_Deposit_TS_Non_Zero.append(copy.deepcopy(Temp_Bay))


Marsh_C_Minus_Bay_Accretion = np.subtract(All_C_Deposited,Bay_C_Deposit_TS_Non_Zero[0:len(All_C_Deposited)])
Total_Marsh_C_Accretion = np.add(Marsh_C_Minus_Bay_Accretion,Initial_Transect_C_Value)
Total_Marsh_C_Accretion[Total_Marsh_C_Accretion<0] = 0

plt.plot(Total_Marsh_C_Accretion[0])
plt.plot(Total_Marsh_C_Accretion[24])
plt.plot(Total_Marsh_C_Accretion[49])
plt.plot(Total_Marsh_C_Accretion[74])
plt.plot(Total_Marsh_C_Accretion[74])
plt.plot(Total_Marsh_C_Accretion[99])
plt.plot(Total_Marsh_C_Accretion[123])
plt.xlim(4500,6000)
plt.show()

z = 20