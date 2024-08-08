import copy

import numpy as np
import time

import matplotlib.pyplot as plt

import os
import imageio



Change_Rates = np.loadtxt('C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Revised_Offshore_Datum\\All_Shoreline_Change_Rates.csv',skiprows=1,delimiter=',')
Subset_Change_Rates = np.loadtxt('C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Revised_Offshore_Datum\\All_Annual_Change_Rates.csv',skiprows=1,delimiter=',')

# Load in the shoreline change trends looking at 1974 - 2020 and 1988 - 2020
Linear_1974 = Change_Rates[:,0]
Linear_1988 = Change_Rates[:,1]
Endpoint_1974 = Change_Rates[:,2]
Endpoint_1988 = Change_Rates[:,3]

# Load in the specific differences between different years
EP_74_88 = Subset_Change_Rates[:,0]
EP_88_97 = Subset_Change_Rates[:,1]
EP_97_09 = Subset_Change_Rates[:,2]
EP_09_20 = Subset_Change_Rates[:,3]
EP_74_97 = Subset_Change_Rates[:,4]
EP_97_20 = Subset_Change_Rates[:,5]
LRR_74_97 = Subset_Change_Rates[:,6]
LRR_97_20 = Subset_Change_Rates[:,7]


os.chdir('C:\\Users\\frank\\PycharmProjects\\CASCADE\\Run_output')

run_name_batch = []

#run_name_batch.append('Hindcast_High_Complexity_Erosion_Values_1974')
#run_name_batch.append('Hindcast_Complex_Background_Erosion_Values_1988')
#run_name_batch.append('Hindcast_Medium_Complexity_Erosion_Values_1988')
#run_name_batch.append('Hindcast_Only_Endpoint_Erosion_Values_1988')
#run_name_batch.append('Hindcast_High_Complexity_Erosion_Values_from_1988_in_1974')
#run_name_batch.append('Hindcast_Low_Complexity_Erosion_Values_from_1988_in_1974')
#run_name_batch.append('Hindcast_Medium_Complexity_Erosion_Values_from_1988_in_1974')
#run_name_batch.append('Hindcast_High_Complexity_Erosion_Values_from_1974_in_1988')
#run_name_batch.append('Hindcast_High_Complexity_Erosion_Values_from_1974_to_1997')
#run_name_batch.append('Hindcast_High_Complexity_Erosion_Values_from_1974_model_1997_2020')
#run_name_batch.append('Hindcast_High_Complexity_Erosion_Values_from_1988_model_1997_2020')
#run_name_batch.append('Duck_Hindcast_S_0')
#run_name_batch.append('Run_1997_2020_with_1974_1997_ss')
#run_name_batch.append('Run_1997_2020_with_altered_1974_1997_ss')
#run_name_batch.append('Buffer_Revised_Test')
#run_name_batch.append('Ocracoke_Bigger_Storms_5_m_1974_1997') # No sink
#run_name_batch.append('Ocracoke_Bigger_Storms_5_m_1974_1997_middle_sinks') # sinks
#run_name_batch.append('Ocracoke_Bigger_Storms_5_m_1997_2020_Unaltered_Inlet_Edge_Erosion')
#run_name_batch.append('Ocracoke_Bigger_Storms_5_m_1997_2020_altered_Inlet_Edge_Erosion')
#run_name_batch.append('Ocracoke_Bigger_Storms_5_m_1997_2020_unaltered_inlet_edge_middle_sinks')
#run_name_batch.append('Ocracoke_Bigger_Storms_5_m_1997_2020_altered_inlet_edge_middle_sinks')
#run_name_batch.append('Ocracoke_Bigger_Storms_5_m_1997_2020_altered_inlet_edge_middle_sinks_Overwash_Filter_90')
#run_name_batch.append('Ocracoke_Bigger_Storms_5_m_1997_2020_altered_inlet_edge_middle_sinks_Overwash_Filter_0')
#run_name_batch.append('Ocracoke_Bigger_Storms_5_m_1997_2020_altered_inlet_edge_middle_sinks_Overwash_Filter_90')
#run_name_batch.append('Big_Storm_Test_No_Middle_Threshold')
#run_name_batch.append('Ocracoke_Bigger_Storms_5_m_1974_1997_middle_sinks')
run_name_batch.append('Buffer_1997_Test')


'''for i in range(0,10):
   run_name_batch.append('O_Hindcast_S_'+str(i))'''



name_prefix = run_name_batch
#nt_run = 46
#nt_run = 32
nt_run = 23
number_barrier3d_models = 70
buffer_length = 15
All_EP_Change = []

for k in range(0,len(run_name_batch)):
    # --------- plot ---------
    output = np.load(run_name_batch[k] + ".npz", allow_pickle=True)
    cascade = output["cascade"]
    cascade = cascade[0]
    b3d = cascade.barrier3d
    ny = np.size(b3d)

    directory = "C:\\Users\\frank\\PycharmProjects\\CASCADE\\"
    # TMax_MGMT = Needed 0
    # TMAX_Sim = Last simulation year of the model 99
    TMax_Sim = nt_run  # Give length of simulation
    roadway_management_ny = [True] * ny




    """
    NOTE THAT THE BEACH REPRESENTATION IS BASED ON A MODEL SPECIFIED BEACH WIDTH. We set the beach width for the
    remaining time steps after the community has been abandoned to the last managed beach width in order to not have a
    huge jump in the back-barrier position in Barrier3D. OTHERWISE, it is meaningless to the dynamics Barrier3D.
    """
    barrier3d = cascade.barrier3d

    # set up the domain; here we just use the first grid, but that could break in future runs
    BarrierLength = barrier3d[0].BarrierLength
    OriginY = int(barrier3d[0].x_s_TS[0])
    total_shoreline_change = cascade._brie_coupler.brie.x_s_dt
    all_shoreline_change = cascade._brie_coupler.brie.x_s_save

    All_Year_1_Shoreline_Position = all_shoreline_change[:,1]
    All_Year_44_Shoreline_Position = all_shoreline_change[:,-1]

    Year_1_Shoreline_Positions = All_Year_1_Shoreline_Position[buffer_length:-buffer_length]
    Year_1_Shoreline_Positions[0] = 1624
    Year_44_Shoreline_Positions = All_Year_44_Shoreline_Position[buffer_length:-buffer_length]
    EP_Change = ((Year_44_Shoreline_Positions - Year_1_Shoreline_Positions)*-1)/nt_run

    All_EP_Change.append(copy.deepcopy(EP_Change))
    mean_change = np.mean(All_EP_Change, axis=0)

All_OV_Flux = []
for l in range(buffer_length, (number_barrier3d_models - buffer_length-1)):
    total_sum = np.sum(b3d[l]._QowTS)
    All_OV_Flux.append(copy.deepcopy(total_sum))



'''
# For comp with 2nd year
output_2 = np.load('New_Flow_Routing_Test.npz', allow_pickle=True)

cascade_2 = output_2["cascade"]
cascade_2 = cascade_2[0]
b3d_2 = cascade_2.barrier3d

All_OV_Flux_2 = []
for l in range(0,39):
    total_sum_2 = np.sum(b3d_2[l]._QowTS)
    All_OV_Flux_2.append(copy.deepcopy(total_sum_2))

# set up the domain; here we just use the first grid, but that could break in future runs
total_shoreline_change_2 = cascade_2._brie_coupler.brie.x_s_dt
all_shoreline_change_2 = cascade_2._brie_coupler.brie.x_s_save

Year_1_Shoreline_Positions_2 = all_shoreline_change_2[:,1]
Year_1_Shoreline_Positions_2[0] = 1624
Year_44_Shoreline_Positions_2 = all_shoreline_change_2[:,-1]
EP_Change_2 = ((Year_44_Shoreline_Positions_2 - Year_1_Shoreline_Positions_2)*-1)/nt_run
All_Change_2 = np.array(EP_Change_2)
mean_change_2 = np.mean(All_Change_2,axis=0)


# Old vs New

OV_Flux_Dif = np.subtract(All_OV_Flux,All_OV_Flux_2)
[29][30][31][35]

OV_29 = (OV_Flux_Dif[29]/All_OV_Flux[29])*100
OV_30 = (OV_Flux_Dif[30]/All_OV_Flux[30])*100
OV_31 = (OV_Flux_Dif[31]/All_OV_Flux[31])*100
OV_35 = (OV_Flux_Dif[35]/All_OV_Flux[35])*100

Change_OV = [OV_29,OV_30,OV_31,OV_35]


Shoreline_dif = np.subtract(mean_change,All_Change_2)
[28],[30],[31]
sdif_28 =Shoreline_dif[28]/mean_change[28]
sdif_30 =Shoreline_dif[30]/mean_change[30]
sdif_31 =Shoreline_dif[31]/mean_change[31]


# Positive means observed is changing faster than modeled
# Negative means model is overpredicting change on Ocracoke

'''
domain_nums = range(11,50)

#plt.plot(domain_nums, Linear_1988, label = '1988 - 2020 (LRR)')
#plt.plot(domain_nums, Linear_1974, label = '1974 - 2020 (LRR)')
#plt.plot(domain_nums, LRR_74_97, label = '1974 - 1997 (LRR)')
plt.plot(domain_nums, LRR_97_20, label = '1997 - 2020 (LRR)')

#plt.plot(domain_nums, Endpoint_1974, label = '1974 - 2020 (EP)')
plt.axhline(y = 0, color = 'k', linestyle = '--')
plt.plot(domain_nums, mean_change, label= 'Modeled Change Rates')
#plt.plot(domain_nums, All_Change_2, label ='Old FLow Routing')
plt.legend()
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('Model Grid')
plt.show()

plt.plot(domain_nums,All_OV_Flux)#, label = 'New Flow Routing')
#plt.plot(domain_nums,All_OV_Flux_2, label ='Old Flow Routing')
plt.axvline(x = 40, color = 'k', linestyle = '--')
plt.axvline(x = 46, color = 'k', linestyle = '--')
plt.legend()

plt.ylabel('Cum Overwash Amount (m^3/yr)')
plt.xlabel('Model Grid')
plt.show()
'''
EP_74_88
EP_88_97
EP_97_09
EP_09_20
M_EP_Change_74_88
M_EP_Change_88_97
M_EP_Change_97_09
M_EP_Change_09_20


plt.plot(domain_nums, EP_74_88, label = '1974 - 1988 (EP)')

#plt.plot(Endpoint_1988, label = 'Endpoint')
plt.axhline(y = 0, color = 'k', linestyle = '--')
plt.plot(domain_nums, M_EP_Change_74_88, label= 'Modeled (1974 - 1988)')
plt.legend()
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('Grid')
plt.show()

plt.plot(domain_nums, EP_88_97, label = '1988 - 1997 (EP)')
plt.axhline(y = 0, color = 'k', linestyle = '--')
plt.plot(domain_nums, M_EP_Change_88_97, label= 'Modeled (1988 - 1997)')
plt.legend()
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('Grid')
plt.show()

plt.plot(domain_nums, EP_97_09, label = '1997 - 2009 (EP)')
plt.axhline(y = 0, color = 'k', linestyle = '--')
plt.plot(domain_nums, M_EP_Change_97_09, label= 'Modeled (1997 - 2009)')
plt.legend()
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('Grid')
plt.show()

plt.plot(domain_nums, EP_09_20, label = '2009 - 2020 (EP)')
plt.axhline(y = 0, color = 'k', linestyle = '--')
plt.plot(domain_nums, M_EP_Change_09_20, label= 'Modeled (2009 - 2020)')
plt.legend()
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('Grid')
plt.show()


plt.plot(domain_nums, Endpoint_1974, label = '1974 - 2020 (EP)')
plt.plot(domain_nums, Endpoint_1988, label = '1988 - 2020 (EP)')
plt.plot(domain_nums,EP_74_97, label = '1974-1997 (EP)')
plt.axhline(y = 0, color = 'k', linestyle = '--')
#plt.plot(domain_nums, mean_change, label= 'Modeled')
plt.legend()
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('Grid')
plt.show()


plt.plot(domain_nums, Linear_1974, label = '1974 - 2020 (LRR)')
plt.plot(domain_nums, Linear_1988, label = '1988 - 2020 (LRR)')
plt.plot(domain_nums, LRR_74_97, label = '1974 - 1997 (LRR)')

#plt.plot(Endpoint_1988, label = 'Endpoint')
plt.axhline(y = 0, color = 'k', linestyle = '--')
#plt.plot(domain_nums, mean_change, label= 'Modeled')
plt.legend()
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('Grid')
plt.show()'''
'''
#plt.plot(domain_nums, Linear_1988, label = '1988 - 2020 (LRR)')
plt.plot(domain_nums, LRR_74_97, label = '1974 - 1997 (LRR)')
plt.plot(domain_nums, LRR_97_20, label = '1997 - 2020 (LRR)')
plt.plot(domain_nums, Linear_1974, label = '1974 - 2020 (LRR)')
plt.axhline(y = 0, color = 'k', linestyle = '--')
plt.legend()
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('Grid')
plt.show()'''

'''
#plt.plot(domain_nums, Endpoint_1974, label = '1974 - 2020 (EP)')
plt.axhline(y = 0, color = 'k', linestyle = '--')
plt.plot(domain_nums, All_EP_Change[0], label= 'Storm Sequence 0')
plt.plot(domain_nums, All_EP_Change[1], label= 'Storm Sequence 1')
plt.plot(domain_nums, All_EP_Change[2], label= 'Storm Sequence 2')
plt.plot(domain_nums, All_EP_Change[3], label= 'Storm Sequence 3', color = 'r')
plt.plot(domain_nums, All_EP_Change[4], label= 'Storm Sequence 4')
plt.plot(domain_nums, All_EP_Change[5], label= 'Storm Sequence 5')
plt.plot(domain_nums, All_EP_Change[6], label= 'Storm Sequence 6')
plt.plot(domain_nums, All_EP_Change[7], label= 'Storm Sequence 7')
plt.plot(domain_nums, All_EP_Change[8], label= 'Storm Sequence 8')
plt.plot(domain_nums, All_EP_Change[9], label= 'Storm Sequence 9')


plt.legend()
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('Grid')
plt.show()

plt.plot(domain_nums, LRR_74_97)
plt.show()'''

'''plt.plot(b3d[0]._QowTS, label= 'Grid 11')
plt.plot(b3d[1]._QowTS, label= 'Grid 12')
plt.plot(b3d[2]._QowTS, label= 'Grid 13')
plt.plot(b3d[3]._QowTS, label= 'Grid 14')
plt.plot(b3d[4]._QowTS, label= 'Grid 15')
plt.plot(b3d[5]._QowTS, label= 'Grid 16')
plt.plot(b3d[6]._QowTS, label= 'Grid 17')
plt.plot(b3d[7]._QowTS, label= 'Grid 18')
plt.plot(b3d[8]._QowTS, label= 'Grid 19')
plt.plot(b3d[9]._QowTS, label= 'Grid 20')

plt.legend()
plt.ylabel('Qty Overwash')
plt.xlabel('Year')
plt.show()

plt.plot(domain_nums, LRR_74_97)
plt.show()


plt.plot(b3d[10]._QowTS, label= 'Grid 21')
plt.plot(b3d[11]._QowTS, label= 'Grid 22')
plt.plot(b3d[12]._QowTS, label= 'Grid 23')
plt.plot(b3d[13]._QowTS, label= 'Grid 24')
plt.plot(b3d[14]._QowTS, label= 'Grid 25')
plt.plot(b3d[15]._QowTS, label= 'Grid 26')
plt.plot(b3d[16]._QowTS, label= 'Grid 27')
plt.plot(b3d[17]._QowTS, label= 'Grid 28')
plt.plot(b3d[18]._QowTS, label= 'Grid 29')
plt.plot(b3d[19]._QowTS, label= 'Grid 30')

plt.legend()
plt.ylabel('Qty Overwash')
plt.xlabel('Year')
plt.show()

plt.plot(b3d[20]._QowTS, label= 'Grid 31')
plt.plot(b3d[21]._QowTS, label= 'Grid 32')
plt.plot(b3d[22]._QowTS, label= 'Grid 33')
plt.plot(b3d[23]._QowTS, label= 'Grid 34')
plt.plot(b3d[24]._QowTS, label= 'Grid 35')
plt.plot(b3d[25]._QowTS, label= 'Grid 36')
plt.plot(b3d[26]._QowTS, label= 'Grid 37')
plt.plot(b3d[27]._QowTS, label= 'Grid 38')
plt.plot(b3d[28]._QowTS, label= 'Grid 39')
plt.plot(b3d[29]._QowTS, label= 'Grid 40')

plt.legend()
plt.ylabel('Qty Overwash')
plt.xlabel('Year')
plt.show()

plt.plot(b3d[30]._QowTS, label= 'Grid 41')
plt.plot(b3d[31]._QowTS, label= 'Grid 42')
plt.plot(b3d[32]._QowTS, label= 'Grid 43')
plt.plot(b3d[33]._QowTS, label= 'Grid 44')
plt.plot(b3d[34]._QowTS, label= 'Grid 45')
plt.plot(b3d[35]._QowTS, label= 'Grid 46')
plt.plot(b3d[36]._QowTS, label= 'Grid 47')
plt.plot(b3d[37]._QowTS, label= 'Grid 48')
plt.plot(b3d[38]._QowTS, label= 'Grid 49')

plt.legend()
plt.ylabel('Qty Overwash')
plt.xlabel('Year')
plt.show()'''