import copy

import numpy as np
import time

import matplotlib.pyplot as plt

import os
import imageio

data_path = '/Users/ceclmac/PycharmProjects/CASCADE/data/Ocracoke_init_data'

# Set color palette
Change_Rates = np.loadtxt('/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/Chapter 3/Revised_Offshore_Datum/All_Shoreline_Change_Rates.csv',skiprows=1,delimiter=',')
Subset_Change_Rates = np.loadtxt('/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/Chapter 3/Revised_Offshore_Datum/All_Annual_Change_Rates.csv',skiprows=1,delimiter=',')

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

os.chdir('/Users/ceclmac/PycharmProjects/CASCADE/Run_output')

run_name_batch = []

# IL
run_name_batch.append('OCR_IL_Natural_ST0_Sink3')
#run_name_batch.append('OCR_IL_Status_Quo_ST0_Sink3')

# I
#run_name_batch.append('OCR_I_Natural_ST0_Sink3')
#run_name_batch.append('OCR_I_Status_Quo_ST0_Sink3')

# IH
#run_name_batch.append('OCR_IH_Natural_ST0_Sink3')
#run_name_batch.append('OCR_IH_Status_Quo_ST0_Sink3')

nt_run = 100
number_barrier3d_models = 70
buffer_length = 15
All_EP_Change = []
b3d_list = []
mean_change = []
All_OV_Flux = []
All_OV_Flux_m3 = []
All_Dune_Rebuilding_TS = []
All_Sandbag_Building_TS = []
All_Road_Relocation_TS = []
All_OW_Year_TS = []
All_Total_EP_Change = []
All_OW_Years_Dict = {}
All_OW_Unique_Years_TS = {}
Total_Island_Width_Change = []
Rate_Island_Width_Change = []
Percent_Island_Width_Change = []
Road_Relocation_Years_Dict = {}
Sandbag_Presence_Years_Dict = {}
Dune_Rebuilding_TS_Dict = {}
Interior_Dune_Construction_TS_Dict = {}
Combined_Dune_Construction_TS_Dict = {}
Road_Drowning_Years_Dict = {}
Island_Drowning = {}

for k in range(0,len(run_name_batch)):
    # --------- plot ---------
    output = np.load(run_name_batch[k] + ".npz", allow_pickle=True)
    cascade = output["cascade"]
    cascade = cascade[0]
    b3d = cascade.barrier3d
    ny = np.size(b3d)

    b3d_list.append(copy.deepcopy(b3d))

    if cascade.b3d_break == 1:
        drowned_cells = {}
        for drown in range(len(b3d)):
            if b3d[drown]._drown_break == 1:
                drowned_cells[str(drown-4)] = len(b3d[drown]._InteriorWidth_AvgTS)
                years_modeled = len(b3d[drown]._InteriorWidth_AvgTS)
                final_year_index = years_modeled - 1
        Island_Drowning[run_name_batch[k]] = drowned_cells
    else:
        years_modeled = nt_run
        final_year_index = years_modeled-1
        Island_Drowning[run_name_batch[k]] = False


    directory = "/Users/frank/PycharmProjects/CASCADE/"
    # TMax_MGMT = Needed 0
    # TMAX_Sim = Last simulation year of the model 99
    TMax_Sim = nt_run  # Give length of simulation
    roadway_management_ny = [True] * ny

    barrier3d = cascade.barrier3d
    # Need to convert to be lists
    # set up the domain; here we just use the first grid, but that could break in future runs
    BarrierLength = barrier3d[0].BarrierLength
    OriginY = int(barrier3d[0].x_s_TS[0])
    total_shoreline_change = cascade._brie_coupler.brie.x_s_dt
    all_shoreline_change = cascade._brie_coupler.brie.x_s_save

    All_Year_1_Shoreline_Position = all_shoreline_change[:,1]
    All_Final_Shoreline_Position = all_shoreline_change[:,final_year_index]

    Year_1_Shoreline_Positions = All_Year_1_Shoreline_Position[buffer_length:-buffer_length]
    Year_1_Shoreline_Positions[0] = 1624
    Year_Final_Shoreline_Positions = All_Final_Shoreline_Position[buffer_length:-buffer_length]
    EP_Change = ((Year_Final_Shoreline_Positions - Year_1_Shoreline_Positions)*-1)/years_modeled
    Total_EP_Change = ((Year_Final_Shoreline_Positions - Year_1_Shoreline_Positions)*-1)


    All_EP_Change.append(copy.deepcopy(EP_Change))
    All_Total_EP_Change.append(copy.deepcopy(Total_EP_Change))
    mean_change.append(copy.deepcopy(np.mean(All_EP_Change, axis=0)))

    # Calculate the change in island width
    Width_TS = []
    Width_Percent_Change = []
    Width_Change_Rate_TS = []
    for ww in range(buffer_length, (number_barrier3d_models - buffer_length-1)):
        Year_1_Width = b3d[ww].InteriorWidth_AvgTS[0]
        Final_Year_Width = b3d[ww].InteriorWidth_AvgTS[final_year_index]
        Width_Change = Final_Year_Width - Year_1_Width
        Width_Change_Rate = Width_Change/years_modeled
        Percent_Change_Temp = (Width_Change/Year_1_Width)*100
        Width_TS.append(copy.deepcopy(Width_Change))
        Width_Percent_Change.append(copy.deepcopy(Percent_Change_Temp))
        Width_Change_Rate_TS.append(copy.deepcopy(Width_Change_Rate))

    # Save model runs values
    Total_Island_Width_Change.append(copy.deepcopy(Width_TS))
    Rate_Island_Width_Change.append(copy.deepcopy(Width_Change_Rate_TS))
    Percent_Island_Width_Change.append(copy.deepcopy(Width_Percent_Change))

    All_OV_Flux_s = []
    for l in range(buffer_length, (number_barrier3d_models - buffer_length-1)):
        total_sum = np.sum(b3d[l]._QowTS)
        All_OV_Flux_s.append(copy.deepcopy(total_sum))
    All_OV_Flux.append(copy.deepcopy(All_OV_Flux_s))

    # Filter the Overwash flux to only focus on the grids of interest and copy to a time series list
    All_OV_Flux_m3_s = []
    for klm in range(buffer_length, (number_barrier3d_models - buffer_length-1)):
        total_sum = (np.sum(b3d[klm]._QowTS))
        x = total_sum * int(BarrierLength)
        All_OV_Flux_m3_s.append(copy.deepcopy(x))
    All_OV_Flux_m3.append(copy.deepcopy(All_OV_Flux_m3_s))

    # Filter the roadbuilding TS to only show B3D grids of interest and copy to final index
    All_Dune_Rebuilding_TS_Temp = []
    Interior_Dune_Building_TS_Temp = []
    Combined_Dune_Building_TS_Temp = []
    Dune_Interior_Building_Years = {}
    Dune_Rebuilding_Years = {}
    Combined_Dune_building_years = {}
    for m in range(buffer_length, (number_barrier3d_models - buffer_length-1)):
        Temp_Years = []
        Temp_Years_Interior = []
        Combined_Temp_Years = []
        All_Dune_Rebuilding = cascade.roadways[m]._dunes_rebuilt_TS
        Interior_Dune_Building = cascade.roadways[m]._interior_dunes_built_TS
        for years in range(len(All_Dune_Rebuilding)):
            if All_Dune_Rebuilding[years] ==1:
                Temp_Years.append(copy.deepcopy(years))
            if Interior_Dune_Building[years] == 1:
                Temp_Years_Interior.append(copy.deepcopy(years))
        if len(Temp_Years) > 0:
            Dune_Rebuilding_Years[str(m-4)] = copy.deepcopy(Temp_Years)
        if len(Temp_Years_Interior) > 0:
            Dune_Interior_Building_Years[str(m-4)] = copy.deepcopy(Temp_Years_Interior)
        if len(Temp_Years) > 0 and len(Temp_Years_Interior) == 0:
            Combined_Dune_building_years[str(m-4)] = (copy.deepcopy(Temp_Years))
        elif len(Temp_Years) == 0 and len(Temp_Years_Interior) > 0:
            Combined_Dune_building_years[str(m-4)] = copy.deepcopy(Temp_Years_Interior)
        elif len(Temp_Years) > 0 and len(Temp_Years_Interior) > 0:
            Combined_Dune_building_years[str(m-4)] = copy.deepcopy(np.sort(np.append(Temp_Years,Temp_Years_Interior)))


        All_Dune_Rebuilding_TS_Temp.append(copy.deepcopy(All_Dune_Rebuilding))
    Combined_Dune_Construction_TS_Dict[str(run_name_batch[k])] = copy.deepcopy(Combined_Dune_building_years)
    Dune_Rebuilding_TS_Dict[str(run_name_batch[k])] = copy.deepcopy(Dune_Rebuilding_Years)
    All_Dune_Rebuilding_TS.append(copy.deepcopy(All_Dune_Rebuilding_TS_Temp))
    Interior_Dune_Construction_TS_Dict[str(run_name_batch[k])] = copy.deepcopy(Dune_Interior_Building_Years)

    # Filter the road relocation TS to only show B3D grids of interest and save to final index for plotting
    All_Road_Relocation_TS_Temp = []
    Road_Relocation_Years = {}
    for m in range(buffer_length, (number_barrier3d_models - buffer_length-1)):
        Temp_Years = []
        All_Road_Relocation = cascade.roadways[m]._road_relocated_TS
        for years in range(len(All_Road_Relocation)):
            if All_Road_Relocation[years] == 1:
                Temp_Years.append(copy.deepcopy(years))
        if len(Temp_Years) > 0:
            Road_Relocation_Years[str(m-4)] = (copy.deepcopy(Temp_Years))
        All_Road_Relocation_TS_Temp.append(copy.deepcopy(All_Road_Relocation))
    All_Road_Relocation_TS.append(copy.deepcopy(All_Road_Relocation_TS_Temp))
    Road_Relocation_Years_Dict[str(run_name_batch[k])] = copy.deepcopy(Road_Relocation_Years)

    # Find times the roadway broke and save the year that it did
    Road_Drowning_Years = {}
    for m in range(buffer_length, (number_barrier3d_models - buffer_length - 1)):
        Temp_Years = []
        Road_Data = cascade.roadways[m]

        if Road_Data.relocation_break == 1 or Road_Data.drown_break == True:
            Temp_Years.append(copy.deepcopy(Road_Data.time_index))

        if len(Temp_Years) > 0:
            Road_Drowning_Years[str(m - 4)] = (copy.deepcopy(Temp_Years))
    Road_Drowning_Years_Dict[str(run_name_batch[k])] = copy.deepcopy(Road_Drowning_Years)

    All_OW_Year_TS_Temp = []
    Years_OW_Unique = []
    OW_Years_Dict = {}
    for klm in range(buffer_length, (number_barrier3d_models - buffer_length-1)):
        Years_OW = []
        for years in range(len(b3d[klm]._QowTS)):
            OW_flux = b3d[klm]._QowTS[years]
            if OW_flux > 0:
                Years_OW.append(copy.deepcopy(years))
                Years_OW_Unique.append(copy.deepcopy(years))
        if len(Years_OW) > 0:
            OW_Years_Dict[str(klm-4)] = (copy.deepcopy(Years_OW))


        Years_OW_Sort = np.sort(np.unique(Years_OW))
        All_OW_Year_TS_Temp.append(copy.deepcopy(Years_OW_Unique))
    All_OW_Year_TS.append(copy.deepcopy(All_OW_Year_TS_Temp))
    All_OW_Years_Dict[str(run_name_batch[k])] = copy.deepcopy(OW_Years_Dict)
    All_OW_Unique_Years_TS[str(run_name_batch[k])] =copy.deepcopy(Years_OW_Sort)

    # Filter the road relocation TS to only show B3D grids of interest and save to final index for plotting
    All_Sandbag_Building_TS_Temp = []
    Sandbag_Years = {}
    for n in range(buffer_length, (number_barrier3d_models - buffer_length-1)):
        Temp_Years = []
        All_Sandbag_Building = cascade._sandbag_Need_TS[n]
        for years in range(len(All_Sandbag_Building)):
            if All_Sandbag_Building[years] == 1:
                Temp_Years.append(copy.deepcopy(years))
        if len(Temp_Years) > 0:
            Sandbag_Years[str(n-4)] = (copy.deepcopy(Temp_Years))
        All_Sandbag_Building_TS_Temp.append(copy.deepcopy(All_Sandbag_Building))
    All_Sandbag_Building_TS.append(copy.deepcopy(All_Sandbag_Building_TS_Temp))
    Sandbag_Presence_Years_Dict[str(run_name_batch[k])] = copy.deepcopy(Sandbag_Years)


domain_nums = range(11,50)
'''
# Set Font

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 16

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Plot all the different OCR sinks with different management options
# IL RSLR
plt.plot(domain_nums, LRR_97_20, label = 'Historic Change',color='grey')

plt.axhline(y = 0, color = 'k', linestyle = '--')
plt.plot(domain_nums, All_EP_Change[0], label= 'Natural')#,color='#1f77b4')
plt.plot(domain_nums, All_EP_Change[1], label= 'Status Quo')#,color='#ff7f0e')

plt.legend()
ax = plt.gca()
ax.set_ylim([-10.1,2.1])
plt.title('Shoreline Change Rate (IL): 2024 - 2124')
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('B3D Domain')
plt.tight_layout()
plt.show()

# Plot all the different OCR sinks with different management options
# I RSLR
plt.plot(domain_nums, LRR_97_20, label = 'Historic Change',color='grey')

plt.axhline(y = 0, color = 'k', linestyle = '--')
plt.plot(domain_nums, All_EP_Change[2], label= 'Natural')#,color='#1f77b4')
plt.plot(domain_nums, All_EP_Change[3], label= 'Status Quo')#,color='#ff7f0e')

plt.legend()
ax = plt.gca()
ax.set_ylim([-10.1,2.1])
plt.title('Shoreline Change Rate (I): 2024 - 2124')
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('B3D Domain')
plt.tight_layout()
plt.show()

# Plot all the different OCR sinks with different management options
# IL RSLR
plt.plot(domain_nums, LRR_97_20, label = 'Historic Change',color='grey')

plt.axhline(y = 0, color = 'k', linestyle = '--')
plt.plot(domain_nums, All_EP_Change[4], label= 'Natural')#,color='#1f77b4')
plt.plot(domain_nums, All_EP_Change[5], label= 'Status Quo')#,color='#ff7f0e')
plt.legend()
ax = plt.gca()
ax.set_ylim([-10.1,2.1])
plt.title('Shoreline Change Rate (IH): 2024 - 2124')
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('B3D Domain')
plt.tight_layout()
plt.show()

# Plot all Status Quo scenarios

plt.plot(domain_nums, LRR_97_20, label = 'Historic Change',color='grey')

plt.axhline(y = 0, color = 'k', linestyle = '--')
plt.plot(domain_nums, All_EP_Change[0], label= 'IL')#,color='#1f77b4')
plt.plot(domain_nums, All_EP_Change[2], label= 'I')#,color='#ff7f0e')
plt.plot(domain_nums, All_EP_Change[4], label= 'IH')#,color='#ff7f0e')

plt.legend()
plt.title('Shoreline Change Rate (Natural): 2024 - 2124')
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('B3D Domain')

ax = plt.gca()
ax.set_ylim([-10.1,2.1])
plt.tight_layout()
plt.show()

# Plot all Road management scenarios

plt.plot(domain_nums, LRR_97_20, label = 'Historic Change',color='grey')

plt.axhline(y = 0, color = 'k', linestyle = '--')
plt.plot(domain_nums, All_EP_Change[1], label= 'IL')#,color='#1f77b4')
plt.plot(domain_nums, All_EP_Change[3], label= 'I')#,color='#ff7f0e')
plt.plot(domain_nums, All_EP_Change[5], label= 'IH')#,color='#ff7f0e')

plt.legend()
plt.title('Shoreline Change Rate (Status Quo): 2024 - 2124')
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('B3D Domain')
ax = plt.gca()
ax.set_ylim([-10.1,2.1])

plt.tight_layout()
plt.show()

print(Island_Drowning)

# Plot total shoreface recession
plt.axhline(y = 0, color = 'k', linestyle = '--')
plt.plot(domain_nums, All_Total_EP_Change[0], label= 'No Management',color='#1f77b4')
plt.plot(domain_nums, All_Total_EP_Change[1], label= 'Status Quo',color='#ff7f0e')

plt.legend()
plt.title('Total Shoreline Change (IL): 2024 - 2124')
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('B3D Domain')
plt.tight_layout()
plt.show()


# Plot overwash
plt.plot(domain_nums, All_OV_Flux_m3[0], label = 'Erosional Inlet',color='#1f77b4')
plt.plot(domain_nums, All_OV_Flux_m3[1], label= 'Accretional Inlet',color='#ff7f0e')

plt.axvline(x = 40, color = 'k', linestyle = '--')
plt.axvline(x = 46, color = 'k', linestyle = '--')
plt.legend()

plt.title('Overwash Volume (m^3) (1974-1997)')
plt.ylabel('Cum Overwash Amount Per B3D Domain (m^3)')
plt.xlabel('B3D Domain')
plt.show()

# Plot roadway abondonment
All_Years = []
All_B3D_Grid_Values = []

for z in range(len(run_name_batch)):
    Year_values = []
    B3D_Grid_values = []
    B3D_Grid = copy.deepcopy(list(Road_Drowning_Years_Dict[run_name_batch[z]].keys()))
    Year = copy.deepcopy(list(Road_Drowning_Years_Dict[run_name_batch[z]].values()))

    for k in range(len(B3D_Grid)):
        Year_values.append(copy.deepcopy(Year[k][0]))
        B3D_Grid_values.append(copy.deepcopy(int(B3D_Grid[k])))
    All_Years.append(copy.deepcopy(Year_values))
    All_B3D_Grid_Values.append(copy.deepcopy(B3D_Grid_values))

# Plot roadway rebuilding as well

Rebuilding_Keys = copy.deepcopy(list(Road_Relocation_Years_Dict[run_name_batch[0]].keys()))

relocation_data = Road_Relocation_Years_Dict[run_name_batch[0]]

for keys in range(len(Rebuilding_Keys)):
    rebuild_years = relocation_data[Rebuilding_Keys[keys]]
    for y in range(len(rebuild_years)):
        plt.scatter(x=(int(Rebuilding_Keys[keys])),y=rebuild_years[y], s = 20, c= '#1f77b4', marker = '1')

Rebuilding_Keys = copy.deepcopy(list(Road_Relocation_Years_Dict[run_name_batch[1]].keys()))

relocation_data = Road_Relocation_Years_Dict[run_name_batch[1]]

for keys in range(len(Rebuilding_Keys)):
    rebuild_years = relocation_data[Rebuilding_Keys[keys]]
    for y in range(len(rebuild_years)):
        plt.scatter(x=(int(Rebuilding_Keys[keys])),y=rebuild_years[y], s = 20, c= '#ff7f0e', marker = '2')


for l in range(len(All_Years[0])):
    plt.scatter(x=All_B3D_Grid_Values[0][l], y=All_Years[0][l], s =30, c = '#1f77b4',
                marker='v')

for l2 in range(len(All_Years[1])):
    plt.scatter(x=All_B3D_Grid_Values[1][l2], y=All_Years[1][l2], s =30, c = '#ff7f0e',
                marker='^')

plt.scatter(x=All_B3D_Grid_Values[0][l], y=All_Years[0][l], s =30, c = '#1f77b4',
                marker='v', label = 'Erosional Inlet')
plt.scatter(x=All_B3D_Grid_Values[1][l2], y=All_Years[1][l2], s=30, c='#ff7f0e',
            marker='^', label='Accretional Inlet')

plt.rc('font',size=18)
plt.title('Road Abandonment: IL_RSLR')
plt.ylabel('Year')
plt.xlabel('B3D Domain')
ax = plt.gca()
ax.set_xlim([10,51])
ax.set_ylim([1,100])
plt.axvline(x = 40, color = 'k', linestyle = '--')
plt.axvline(x = 46, color = 'k', linestyle = '--')
plt.legend()
plt.tight_layout()
plt.show()



# Plot Ocracoke Storms
#plt.plot(domain_nums, LRR_74_97, label = 'Historic Change',color='grey')
plt.plot(domain_nums, LRR_97_20, label = 'Historic Change',color='grey')

plt.axhline(y = 0, color = 'k', linestyle = '--')
plt.plot(domain_nums, All_EP_Change[1], label= 'OCR: Modified', color='purple')
plt.plot(domain_nums, All_EP_Change[4], label= 'OCR: Baseline', color = 'dodgerblue')

plt.legend()
plt.title('Historic vs Modeled Change: 1997-2021')
#plt.title('Historic vs Modeled Change: 1974-1997')
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('B3D Domain')
plt.show()


# Plot Duck Storms
plt.plot(domain_nums, LRR_74_97, label = 'Historic Change',color='grey')
#plt.plot(domain_nums, LRR_97_20, label = 'Historic Change',color='grey')
plt.axhline(y = 0, color = 'k', linestyle = '--')
plt.plot(domain_nums, All_EP_Change[2], label= 'Duck: Modified',color='orange')
plt.plot(domain_nums, All_EP_Change[5], label= 'Duck: Baseline',color = 'saddlebrown')

plt.legend()
#plt.title('Historic vs Modeled Change: 1997-2021')
plt.title('Historic vs Modeled Change: 1974-1997')
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('B3D Domain')
plt.show()


# Plot NCB derived
plt.plot(domain_nums, LRR_74_97, label = 'Historic Change',color='grey')
#plt.plot(domain_nums, LRR_97_20, label = 'Historic Change',color='grey')
plt.axhline(y = 0, color = 'k', linestyle = '--')

plt.plot(domain_nums, All_EP_Change[0], label= 'NCB Modified',color='green')
plt.plot(domain_nums, All_EP_Change[3], label= 'NCB: Baseline', color ='lime')

plt.legend()
#plt.title('Historic vs Modeled Change: 1997-2021')
plt.title('Historic vs Modeled Change: 1974-1997')
plt.ylabel('Shoreline Change Rate (m/yr)')
plt.xlabel('B3D Domain')
plt.show()


# Plot Overwash values
# Plot Duck
plt.plot(domain_nums,All_OV_Flux_m3[2], label = 'Duck: Modified',color='orange')
plt.plot(domain_nums,All_OV_Flux_m3[5], label = 'Duck: Baseline',color='saddlebrown')

plt.axvline(x = 40, color = 'k', linestyle = '--')
plt.axvline(x = 46, color = 'k', linestyle = '--')
plt.legend()

plt.title('Overwash Volume (m^3)')
plt.ylabel('Cum Overwash Amount Per B3D Domain (m^3)')
plt.xlabel('B3D Domain')
plt.show()


# Plot OCR
plt.plot(domain_nums,All_OV_Flux_m3[1], label = 'OCR: Modified',color='purple')
plt.plot(domain_nums,All_OV_Flux_m3[4], label = 'OCR: Baseline',color='dodgerblue')

plt.axvline(x = 40, color = 'k', linestyle = '--')
plt.axvline(x = 46, color = 'k', linestyle = '--')
plt.legend()

plt.title('Overwash Volume (m^3)')
plt.ylabel('Cum Overwash Amount Per B3D Domain (m^3)')
plt.xlabel('B3D Domain')
plt.show()

# Plot NCB
plt.plot(domain_nums,All_OV_Flux_m3[0], label = 'NCB: Modified',color='green')
plt.plot(domain_nums,All_OV_Flux_m3[3], label = 'NCB: Baseline',color='lime')

plt.axvline(x = 40, color = 'k', linestyle = '--')
plt.axvline(x = 46, color = 'k', linestyle = '--')
plt.legend()

plt.title('Overwash Volume (m^3)')
plt.ylabel('Cum Overwash Amount Per B3D Domain (m^3)')
plt.xlabel('B3D Domain')
plt.show()
'''