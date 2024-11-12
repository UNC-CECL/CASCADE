import copy

import numpy as np
import pandas as pd
from scipy import stats as st
import os

Change_Rates = np.loadtxt('C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Revised_Offshore_Datum\\All_Shoreline_Change_Rates.csv',skiprows=1,delimiter=',')
Subset_Change_Rates = np.loadtxt('C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Revised_Offshore_Datum\\All_Annual_Change_Rates.csv',skiprows=1,delimiter=',')

#os.chdir('C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Model Runs\\Future Runs')
os.chdir('E:\\Model_Runs')


Save_Path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Model Runs\\Summary_Values\\'

Management_Name = ['Status_Quo','Natural']
RSLR_Rate = ['IL','I','IH']
Sink_Name = ['Erosional_Sink','Accretional_Sink']

Base_Name_List = []

for management in range(len(Management_Name)):
    for RSLR in range(len(RSLR_Rate)):
        #for sinks in range(len(Sink_Name)):
        Base_Name = 'OCR_' + str(RSLR_Rate[RSLR]) + '_' + str(Management_Name[management])
        Base_Name_List.append(copy.deepcopy(Base_Name))


Base_Name = Base_Name_List[0]

def Process_Batch(Base_Name,
                  Sink_Name,
                  Save_Path):
    name_list = []

    # IL
    for runs in range(0,100):
        name_list.append(str(Base_Name)+'_S'+str(runs)+'_'+str(Sink_Name))

    nt_run = 126
    number_barrier3d_models = 70
    buffer_length = 15
    All_EP_Change = []
    All_Roadway_Abandonment = []
    Island_Drowning = {}
    Years_Modeled_List = []
    Relocation_TS = []
    Frequency_TS = []
    number_sandbags_TS = []
    sandbag_duration_TS = []
    sandbag_areas_TS = []
    island_width_change_TS = []

    Model_Run_Years,Drowning_Domain_Locations, Cascade_List = Process_Data(run_name_batch = name_list)

    for runs in range(len(Cascade_List)):
        # Calculate shoreline change values
        shoreline_change = Calculate_Average_Shoreline_Change(cascade=Cascade_List[runs],
                                                              years_modeled=Model_Run_Years[runs],
                                                              buffer_length=buffer_length)
        All_EP_Change.append(copy.deepcopy(shoreline_change))
        # Calculate Roadway Abandonment metrics
        roadway_abandonment = Calculate_Roadway_Abandonmet(cascade=Cascade_List[runs],
                                                           years_modeled=Model_Run_Years[runs],
                                                           buffer_length=buffer_length,
                                                           number_barrier3d_models = number_barrier3d_models)
        All_Roadway_Abandonment.append(copy.deepcopy(roadway_abandonment))
        # Calculate roadway relocation
        roadway_relocation, relocation_frequency = Calculate_Roadway_Relocation(cascade=Cascade_List[runs],
                                                                                years_modeled=Model_Run_Years[runs],
                                                                                buffer_length=buffer_length,
                                                                                number_barrier3d_models = number_barrier3d_models)

        number_sandbags, sandbag_duration, sandbag_areas = Calculate_Sandbag_Years(cascade=Cascade_List[runs],
                                                                                   years_modeled=Model_Run_Years[runs],
                                                                                   buffer_length=buffer_length,
                                                                                   number_barrier3d_models = number_barrier3d_models)

        island_width_change = Calculate_Island_Interior_Width_Change(cascade=Cascade_List[runs],
                                                                     years_modeled=Model_Run_Years[runs],
                                                                     buffer_length=buffer_length,
                                                                     number_barrier3d_models = number_barrier3d_models)
        Relocation_TS.append(copy.deepcopy(roadway_relocation))
        Frequency_TS.append(copy.deepcopy(relocation_frequency))
        number_sandbags_TS.append(copy.deepcopy(number_sandbags))
        sandbag_duration_TS.append(copy.deepcopy(sandbag_duration))
        sandbag_areas_TS.append(copy.deepcopy(sandbag_areas))
        island_width_change_TS.append(copy.deepcopy(island_width_change))

    # Calculate the mean values for all runs
    Mean_Shoreline_Change_Rate = np.mean(All_EP_Change,axis=0)
    Mean_Roadway_Abandonment = np.mean(All_Roadway_Abandonment,axis=0)
    Mean_Roadway_Relocation = np.mean(Relocation_TS,axis=0)
    Mean_Roadway_Frequency = np.mean(Frequency_TS, axis=0)
    Mean_Sandbag_Duration = np.mean(sandbag_duration_TS, axis=0)
    Mean_Number_Sandbags = np.mean(number_sandbags_TS, axis=0)
    Mean_Island_Interior_Change = np.mean(island_width_change_TS, axis=0)


    Break_Section,Break_Domain_Location, Section_Nums = Find_Most_Common_Drowning_Area(Drowned_Domains=Drowning_Domain_Locations)

    Avg_Break_Year = np.mean(Model_Run_Years)

    # Collate all data from 100 model runs
    All_Values_Dict = {'All_EP_Change':All_EP_Change,
                       'All_Roadway_Abandonment':All_Roadway_Abandonment,
                       'Relocation_TS':Relocation_TS,
                       'Frequency_TS':Frequency_TS,
                       'Sandbag_Duration_TS':sandbag_duration_TS,
                       'Number_Sandbags_TS':number_sandbags_TS,
                       'Island_Width_Change_TS':island_width_change_TS,
                       'Model_Run_Years':Model_Run_Years,
                       #'Break_Domain_Locations':Section_Nums
                       }

    All_Values_Data_Frame = pd.DataFrame(All_Values_Dict)


    Export_Values_Dict = {
        'Mean_Shoreline_Change_Rate':Mean_Shoreline_Change_Rate,
        'Mean_Roadway_Abandonment':Mean_Roadway_Abandonment,
        'Roadway_Relocations':Mean_Roadway_Relocation,
        'Roadway_Relocation_Frequency':Mean_Roadway_Frequency,
        'Sandbag_Duration':Mean_Sandbag_Duration,
        'Number_of_Sandbag_Emplacements':Mean_Number_Sandbags,
        'Island_Interior_Change':Mean_Island_Interior_Change,
        'Island_Drown_Year':Avg_Break_Year,
        'Island_Drown_Domain':Break_Domain_Location,
        'Island_Drown_Section':Break_Section
    }

    Export_DF = pd.DataFrame(Export_Values_Dict)

    Full_Save_Path = Save_Path+Base_Name+'_'+Sink_Name+'.csv'

    #Export_DF.to_csv(Full_Save_Path)

    # Save yearly data as .pkl
    Full_Save_Path_PKL = Save_Path+Base_Name+'_'+Sink_Name+'.pkl'
    All_Values_Data_Frame.to_pickle(Full_Save_Path_PKL)


    return(Export_DF)

def Find_Most_Common_Drowning_Area(Drowned_Domains):
    # Set the numbers that comprise the 6 groups
    Greatest_Len = -50
    Most_Common_Break = -50

    All_Len = [0,
               0,
               0,
               0,
               0,
               0]

    if len(Drowned_Domains) > 0:
        Section_1 = range(11,20)
        Section_2 = range(20,30)
        Section_3 = range(30,34)
        Section_4 = range(34,40)
        Section_5 = range(40,47)
        Section_6 = range(47,50)

        # Create blank lists to be
        S1_List = []
        S2_List = []
        S3_List = []
        S4_List = []
        S5_List = []
        S6_List = []

        for drowned_cells in range(len(Drowned_Domains)):
            Domain = Drowned_Domains[drowned_cells]
            if Domain >= Section_1[0] and Domain <=Section_1[-1]:
                S1_List.append(copy.deepcopy(Domain))
            elif Domain >= Section_2[0] and Domain <=Section_2[-1]:
                S2_List.append(copy.deepcopy(Domain))
            elif Domain >= Section_3[0] and Domain <=Section_3[-1]:
                S3_List.append(copy.deepcopy(Domain))
            elif Domain >= Section_4[0] and Domain <=Section_4[-1]:
                S4_List.append(copy.deepcopy(Domain))
            elif Domain >= Section_5[0] and Domain <=Section_5[-1]:
                S5_List.append(copy.deepcopy(Domain))
            elif Domain >= Section_6[0] and Domain <=Section_6[-1]:
                S6_List.append(copy.deepcopy(Domain))

        # Find the most common drowning group
        S1_Len = len(S1_List)
        S2_Len = len(S2_List)
        S3_Len = len(S3_List)
        S4_Len = len(S4_List)
        S5_Len = len(S5_List)
        S6_Len = len(S6_List)

        All_Len = [S1_Len,
                   S2_Len,
                   S3_Len,
                   S4_Len,
                   S5_Len,
                   S6_Len]

        All_Breaks = [S1_List,
                      S2_List,
                      S3_List,
                      S4_List,
                      S5_List,
                      S6_List]

        Long_Len = 0
        for most in range(0,len(All_Len)):
            if All_Len[most] > Long_Len:
                Greatest_Len = copy.deepcopy(most+1)
                Long_Len = copy.deepcopy(All_Len[most])

        Most_Common_Break = st.mode(All_Breaks[Greatest_Len-1])[0][0]
    return(Greatest_Len,Most_Common_Break, All_Len)

def Process_Data(run_name_batch):
    cascade_list = []
    Island_Drowning_Location_List = []
    Years_Modeled_List = []
    for k in range(0,len(run_name_batch)):
        # --------- plot ---------
        output = np.load(run_name_batch[k] + ".npz", allow_pickle=True)["cascade"]
        cascade = output[0]
        #cascade = cascade[0]
        cascade_list.append(copy.deepcopy(cascade))
        b3d = cascade.barrier3d
        ny = np.size(b3d)
        print(str(k)+' is loaded. Break is equal to '+str(cascade.b3d_break))

        if cascade.b3d_break == 1:
            drowned_cells = {}
            for drown in range(len(b3d)):
                if b3d[drown]._drown_break == 1:
                    drowned_cells[str(drown-4)] = len(b3d[drown]._InteriorWidth_AvgTS)
                    years_modeled = len(b3d[drown]._InteriorWidth_AvgTS)
                    final_year_index = years_modeled - 1
                    Island_Drowning_Location_List.append(copy.deepcopy(drown - 4))
            #Island_Drowning[run_name_batch[k]] = drowned_cells
        else:
            years_modeled = cascade._nt
            final_year_index = years_modeled-1
            #Island_Drowning[run_name_batch[k]] = False
        Years_Modeled_List.append(copy.deepcopy(years_modeled))
    return(Years_Modeled_List,Island_Drowning_Location_List,cascade_list)

def Calculate_Average_Shoreline_Change(cascade, years_modeled, buffer_length):
    final_year_index = years_modeled -1
    barrier3d = cascade.barrier3d
    # Need to convert to be lists
    # set up the domain; here we just use the first grid, but that could break in future runs
    total_shoreline_change = cascade._brie_coupler.brie.x_s_dt
    all_shoreline_change = cascade._brie_coupler.brie.x_s_save

    All_Year_1_Shoreline_Position = all_shoreline_change[:, 1]
    All_Final_Shoreline_Position = all_shoreline_change[:, final_year_index]

    Year_1_Shoreline_Positions = All_Year_1_Shoreline_Position[buffer_length:-buffer_length]
    Year_1_Shoreline_Positions[0] = 1624
    Year_Final_Shoreline_Positions = All_Final_Shoreline_Position[buffer_length:-buffer_length]
    EP_Change = ((Year_Final_Shoreline_Positions - Year_1_Shoreline_Positions) * -1) / years_modeled
    return(EP_Change)

def Calculate_Island_Interior_Width_Change(cascade, years_modeled, buffer_length, number_barrier3d_models):
    final_year_index = years_modeled -1
    Width_TS = []
    Width_Percent_Change = []
    Width_Change_Rate_TS = []
    for ww in range(buffer_length, (number_barrier3d_models - buffer_length - 1)):
        b3d = cascade.barrier3d[ww]

        Year_1_Width = b3d.InteriorWidth_AvgTS[0]
        Final_Year_Width = b3d.InteriorWidth_AvgTS[final_year_index]
        Width_Change = Final_Year_Width - Year_1_Width
        Width_Change_Rate = Width_Change / years_modeled
        Percent_Change_Temp = (Width_Change / Year_1_Width) * 100
        Width_TS.append(copy.deepcopy(Width_Change))
        Width_Percent_Change.append(copy.deepcopy(Percent_Change_Temp))
        Width_Change_Rate_TS.append(copy.deepcopy(Width_Change_Rate))

    # Save model runs values
    #Total_Island_Width_Change.append(copy.deepcopy(Width_TS))
    #Rate_Island_Width_Change.append(copy.deepcopy(Width_Change_Rate_TS))
    #Percent_Island_Width_Change.append(copy.deepcopy(Width_Percent_Change))

    return(Width_Percent_Change)

def Calculate_Roadway_Abandonmet(cascade, years_modeled, buffer_length, number_barrier3d_models):
    # Find times the roadway broke and save the year that it did
    Road_Drowning_Years = []
    for m in range(buffer_length, (number_barrier3d_models - buffer_length - 1)):
        Road_Data = cascade.roadways[m]
        if Road_Data.relocation_break == 1 or Road_Data.drown_break == True:
            Road_Drowning_Years.append(copy.deepcopy(Road_Data.time_index))
        elif Road_Data.drown_break == int(0) and Road_Data.time_index == 1:
            Road_Drowning_Years.append(copy.deepcopy(1))
        else:
            Road_Drowning_Years.append(copy.deepcopy(years_modeled))
    return(Road_Drowning_Years)

def Calculate_Roadway_Relocation(cascade, years_modeled,buffer_length, number_barrier3d_models):
    Relocations = []
    Frequency = []
    for m in range(buffer_length, (number_barrier3d_models - buffer_length - 1)):
        Road_Data = cascade.roadways[m]
        Num_Relocations = np.sum(Road_Data._road_relocated_TS)
        Relocations.append(copy.deepcopy(Num_Relocations))
        if Num_Relocations == 0:
            Frequency.append(0)
        else:
            Frequency.append(copy.deepcopy(Road_Data.time_index/Num_Relocations))
    return(Relocations, Frequency)

def Calculate_Sandbag_Years(cascade, years_modeled,buffer_length, number_barrier3d_models):
    last_index = years_modeled -1
    Number_Sandbag_Emplacements_List = []
    Mean_Sandbag_Length_List = []
    Sandbag_Emplacement_Domains = []
    for m in range(buffer_length, (number_barrier3d_models - buffer_length - 1)):
        Sandbag_TS = cascade._sandbag_Need_TS[m]
        Years_Of_Sandbags = 0
        Years_Of_Sandbags_TS = []

        for years in range(1,len(Sandbag_TS)):
            Past_Sandbag = Sandbag_TS[years-1]
            Present_Sandbag = Sandbag_TS[years]
            if Past_Sandbag == Present_Sandbag and Present_Sandbag == 1:
                Years_Of_Sandbags += 1
            elif Past_Sandbag != Present_Sandbag and Present_Sandbag == 1:
                Years_Of_Sandbags = 1
            else:
                if Years_Of_Sandbags == 1:
                    Years_Of_Sandbags_TS.append(copy.deepcopy(Years_Of_Sandbags))
                elif Years_Of_Sandbags > 1:
                    Years_Of_Sandbags_TS.append(copy.deepcopy(Years_Of_Sandbags))

                Years_Of_Sandbags = 0
            if years == (len(Sandbag_TS)-1) and Years_Of_Sandbags >= 1:
                Years_Of_Sandbags_TS.append(copy.deepcopy(Years_Of_Sandbags))

        Number_Sandbag_Emplacements = len(Years_Of_Sandbags_TS)
        if Number_Sandbag_Emplacements == 0:
            Mean_Sandbag_Length = 0
        else:
            Mean_Sandbag_Length = np.mean(Years_Of_Sandbags_TS)
        if Number_Sandbag_Emplacements > 0:
            Sandbag_Emplacement_Domains.append(copy.deepcopy(m-4))
        Number_Sandbag_Emplacements_List.append(copy.deepcopy(Number_Sandbag_Emplacements))
        Mean_Sandbag_Length_List.append(copy.deepcopy(Mean_Sandbag_Length))
    return (Number_Sandbag_Emplacements_List,Mean_Sandbag_Length_List,Sandbag_Emplacement_Domains)

Output_DF = Process_Batch(Base_Name = Base_Name_List[3],Sink_Name =Sink_Name[1],Save_Path = Save_Path)

print('Hello')

'''
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
    '''

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