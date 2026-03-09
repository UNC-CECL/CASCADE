# Collate data from CASCADE runs and export as .pkl files
# BF 2/16/25

import copy

import numpy as np
import pandas as pd
from scipy import stats as st
import os

#os.chdir('C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Model Runs\\Future Runs')
os.chdir('E:\\Model_Runs')


Save_Path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Model Runs\\Summary_Values\\'

Management_Name = ['Status_Quo','Natural', 'Nourishment']
Management_Name = ['Natural']
RSLR_Rate = ['IL','I','IH']
#RSLR_Rate = ['IL']
Sink_Name = ['Erosional_Sink','Accretional_Sink']
#Sink_Name = ['Erosional_Sink']

Storm_Level = 'Baseline'

Base_Name_List = []

for management in range(len(Management_Name)):
    for RSLR in range(len(RSLR_Rate)):
        #for sinks in range(len(Sink_Name)):
        Base_Name = 'OCR_' + str(RSLR_Rate[RSLR]) + '_' + str(Management_Name[management])
        Base_Name_List.append(copy.deepcopy(Base_Name))


Base_Name = Base_Name_List[0]

def Process_Batch(Base_Name,
                  Sink_Name,
                  Save_Path,
                  Storm_Scenario):
    name_list = []

    # IL
    if Storm_Scenario == 'Baseline':
        for runs in range(0,100):
            name_list.append(str(Base_Name)+'_S'+str(runs)+'_T_'+str(Sink_Name))
    else:
        for runs in range(0,100):
            name_list.append(str(Base_Name)+'_S'+str(runs)+'_10_'+str(Sink_Name))

    nt_run = 100
    number_barrier3d_models = 70
    buffer_length = 15
    All_EP_Change = []
    All_Roadway_Abandonment = []
    All_Abandonment_Reason = []
    All_Initial_Island_Elevations = []
    All_Final_Island_Elevations = []
    All_Elevation_Changes = []
    Island_Drowning = {}
    Years_Modeled_List = []
    Relocation_TS = []
    Frequency_TS = []
    number_sandbags_TS = []
    sandbag_duration_TS = []
    sandbag_areas_TS = []
    island_width_change_TS = []
    Model_Run_Years = []
    Drowning_Domain_Locations = []
    All_Shoreline_Positions = []
    Total_Volume_TS = []
    All_Nourishment_TS = []
    Total_Overwash_TS = []
    Yearly_Overwash_TS = []
    All_Subaerial_Cells = []
    Bay_Shoreline_TS = []
    Island_Width_TS = []


    for runs in range(0,len(name_list)):
        #for runs in range(0,2):
        Model_Run_Year,Drowning_Domain_Location, Cascade_List = Process_Data(run_name_batch = name_list, load_index = runs)

        # Calculate shoreline change values
        shoreline_change, all_shoreline_positions = Calculate_Average_Shoreline_Change(cascade=Cascade_List,
                                                              years_modeled=Model_Run_Year,
                                                              buffer_length=buffer_length)
        All_EP_Change.append(copy.deepcopy(shoreline_change))
        All_Shoreline_Positions.append(copy.deepcopy(all_shoreline_positions))
        # Calculate Roadway Abandonment metrics
        roadway_abandonment,abandonmet_reason = Calculate_Roadway_Abandonmet(cascade=Cascade_List,
                                                           years_modeled=Model_Run_Year,
                                                           buffer_length=buffer_length,
                                                           number_barrier3d_models = number_barrier3d_models)
        All_Roadway_Abandonment.append(copy.deepcopy(roadway_abandonment))
        All_Abandonment_Reason.append(copy.deepcopy(abandonmet_reason))
        # Calculate roadway relocation
        roadway_relocation, relocation_frequency = Calculate_Roadway_Relocation(cascade=Cascade_List,
                                                                                years_modeled=Model_Run_Year,
                                                                                buffer_length=buffer_length,
                                                                                number_barrier3d_models = number_barrier3d_models)

        number_sandbags, sandbag_duration, sandbag_areas = Calculate_Sandbag_Years(cascade=Cascade_List,
                                                                                   years_modeled=Model_Run_Year,
                                                                                   buffer_length=buffer_length,
                                                                                   number_barrier3d_models = number_barrier3d_models)

        island_width_change,Island_Width_Vals = Calculate_Island_Interior_Width_Change(cascade=Cascade_List,
                                                                     years_modeled=Model_Run_Year,
                                                                     buffer_length=buffer_length,
                                                                     number_barrier3d_models = number_barrier3d_models)

        Total_Volume,Total_Nourishment_Per_Grid,All_Nourishment = Calculate_Nourishment_Volume(cascade=Cascade_List,
                                                                                                  years_modeled=Model_Run_Year,
                                                                                                  buffer_length=buffer_length,
                                                                                                  number_barrier3d_models=number_barrier3d_models
                                                                                                  )

        Total_OW_Volume,Yearly_Total_OW_Volume = Calculate_Overwash_Volume(cascade = Cascade_List,
                                                       years_modeled=Model_Run_Years,
                                                       buffer_length=buffer_length,
                                                       number_barrier3d_models=number_barrier3d_models)

        Domain_Bay_Shoreline = Calculate_Bay_Distance(cascade=Cascade_List,
                                                      buffer_length=buffer_length,
                                                      number_barrier3d_models=number_barrier3d_models)

        # Elevation_Change_Output, Initial_Elevation_Output, Final_Elevation_Output = Calculate_Island_Elevation_Metrics(
        #     cascade=Cascade_List,
        #     buffer_length=buffer_length,
        #     number_barrier3d_models=number_barrier3d_models)
        #
        # Subaerial_TS = Calculate_Subaerial_Elevation(
        #     cascade=Cascade_List,
        #     buffer_length=buffer_length,
        #     number_barrier3d_models=number_barrier3d_models)

        Total_Volume_TS.append(copy.deepcopy(Total_Volume))
        All_Nourishment_TS.append(copy.deepcopy(All_Nourishment))
        Relocation_TS.append(copy.deepcopy(roadway_relocation))
        Frequency_TS.append(copy.deepcopy(relocation_frequency))
        number_sandbags_TS.append(copy.deepcopy(number_sandbags))
        sandbag_duration_TS.append(copy.deepcopy(sandbag_duration))
        sandbag_areas_TS.append(copy.deepcopy(sandbag_areas))
        island_width_change_TS.append(copy.deepcopy(island_width_change))
        Model_Run_Years.append(copy.deepcopy(Model_Run_Year))
        Drowning_Domain_Locations.append(copy.deepcopy(Drowning_Domain_Location))
        Total_Overwash_TS.append(copy.deepcopy(Total_OW_Volume))
        Yearly_Overwash_TS.append(copy.deepcopy(Yearly_Total_OW_Volume))
        Island_Width_TS.append(copy.deepcopy(Island_Width_Vals))
        #All_Initial_Island_Elevations.append(copy.deepcopy(Initial_Elevation_Output))
        #All_Final_Island_Elevations.append(copy.deepcopy(Final_Elevation_Output))
        #All_Elevation_Changes.append(copy.deepcopy(Elevation_Change_Output))
        #All_Subaerial_Cells.append(copy.deepcopy(Subaerial_TS))
        Bay_Shoreline_TS.append(copy.deepcopy(Domain_Bay_Shoreline))

        z = 10


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
                       'Island_Width_TS':Island_Width_TS,
                       'Model_Run_Years':Model_Run_Years,
                       'Break_Domain_Locations':Drowning_Domain_Locations,
                       'All_Shoreline_Positions':All_Shoreline_Positions,
                       'Total_Nourishment_Volume':Total_Volume_TS,
                       'All_Nourishment_TS':All_Nourishment_TS,
                       'All_Overwash_volume':Total_Overwash_TS,
                       'All_Yearly_OW_volume':Yearly_Overwash_TS,
                       'Roadway_Abandonment_Reason':All_Abandonment_Reason,
                       'Bay_Shoreline_TS':Bay_Shoreline_TS
                       }

        #Elev_Values = {'Island_Initial_Elevation':All_Initial_Island_Elevations,
         #              'Island_Final_Elevation':All_Final_Island_Elevations,
          #             'Island_Elevation_Change':All_Elevation_Changes,
           #            'Subaerial_Elevation_Cells':All_Subaerial_Cells}'''

    All_Values_Data_Frame = pd.DataFrame(All_Values_Dict)
    #All_Elev_Values_DF = pd.DataFrame(Elev_Values)


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

    Full_Save_Path = Save_Path+Base_Name+'_T_'+Sink_Name+'_'+str(Storm_Scenario)+'.csv'

    Export_DF.to_csv(Full_Save_Path)

    # Save yearly data as .pkl
    Full_Save_Path_PKL = Save_Path+Base_Name+'_T_'+Sink_Name+'_'+str(Storm_Scenario)+'.pkl'
    #Full_Save_Path_Elev = Save_Path+Base_Name+'_'+Sink_Name+'_'+str(Storm_Scenario)+'_elev.pkl'

    All_Values_Data_Frame.to_pickle(Full_Save_Path_PKL)
    #All_Elev_Values_DF.to_pickle(Full_Save_Path_Elev)

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
        if Greatest_Len < 0:
            Most_Common_Break = -50
        else:
            Most_Common_Break = st.mode(All_Breaks[Greatest_Len-1])[0][0]
    return(Greatest_Len,Most_Common_Break, All_Len)

def Process_Data(run_name_batch,load_index):
    cascade_list = []
    Island_Drowning_Location_List = []
    Years_Modeled_List = []
    #for k in range(0,len(run_name_batch)):
    # --------- plot ---------
    output = np.load(run_name_batch[load_index] + ".npz", allow_pickle=True)["cascade"]
    cascade = output[0]
    #cascade = cascade[0]
    #cascade_list.append(copy.deepcopy(cascade))
    b3d = cascade.barrier3d
    ny = np.size(b3d)
    print(str(load_index)+' is loaded. Break is equal to '+str(cascade.b3d_break))

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
        Island_Drowning_Location_List = [0]
    Years_Modeled_List.append(copy.deepcopy(years_modeled))
    return(Years_Modeled_List[0],Island_Drowning_Location_List[0],cascade)

def Calculate_Nourishment_Volume(cascade,years_modeled,buffer_length,number_barrier3d_models):
    final_year_index = years_modeled -1
    Nourishment_Data = cascade.nourishments
    All_Nourishment_TS = []
    Total_Nourishment_Per_Grid = []
    Nourishment_Event_Volume = cascade.nourishment_volume
    for l in range(buffer_length, (number_barrier3d_models - buffer_length)):
        All_Nourishment_TS.append(copy.deepcopy(Nourishment_Data[l].nourishment_volume_TS))
        Total_Nourishment_Per_Grid.append(np.sum(Nourishment_Data[l].nourishment_volume_TS))
    Total_Volume = np.sum(Total_Nourishment_Per_Grid)
    return (Total_Volume,Total_Nourishment_Per_Grid,All_Nourishment_TS)

def Calculate_Average_Shoreline_Change(cascade, years_modeled, buffer_length):
    final_year_index = years_modeled -1
    barrier3d = cascade.barrier3d
    # Need to convert to be lists
    # set up the domain; here we just use the first grid, but that could break in future runs
    total_shoreline_change = cascade._brie_coupler.brie.x_s_dt
    all_shoreline_change = cascade._brie_coupler.brie.x_s_save
    all_focused_shoreline_change = all_shoreline_change[buffer_length:-buffer_length,:]

    All_Year_1_Shoreline_Position = all_shoreline_change[:, 1]
    All_Final_Shoreline_Position = all_shoreline_change[:, final_year_index]

    Year_1_Shoreline_Positions = All_Year_1_Shoreline_Position[buffer_length:-buffer_length]
    Year_1_Shoreline_Positions[0] = 1624
    Year_Final_Shoreline_Positions = All_Final_Shoreline_Position[buffer_length:-buffer_length]
    EP_Change = ((Year_Final_Shoreline_Positions - Year_1_Shoreline_Positions) * -1) / years_modeled
    return(EP_Change,all_focused_shoreline_change)

def Calculate_Island_Interior_Width_Change(cascade, years_modeled, buffer_length, number_barrier3d_models):
    final_year_index = years_modeled -1
    Width_TS = []
    Width_Percent_Change = []
    Width_Change_Rate_TS = []
    All_Focused_Domain_Width_TS = []
    for ww in range(buffer_length, (number_barrier3d_models - buffer_length - 1)):
        b3d = cascade.barrier3d[ww]
        Domain_Width = b3d.InteriorWidth_AvgTS
        Year_1_Width = b3d.InteriorWidth_AvgTS[0]
        Final_Year_Width = b3d.InteriorWidth_AvgTS[final_year_index]
        Width_Change = Final_Year_Width - Year_1_Width
        Width_Change_Rate = Width_Change / years_modeled
        Percent_Change_Temp = (Width_Change / Year_1_Width) * 100
        Width_TS.append(copy.deepcopy(Width_Change))
        Width_Percent_Change.append(copy.deepcopy(Percent_Change_Temp))
        Width_Change_Rate_TS.append(copy.deepcopy(Width_Change_Rate))
        All_Focused_Domain_Width_TS.append(copy.deepcopy(Domain_Width))

    # Save model runs values
    #Total_Island_Width_Change.append(copy.deepcopy(Width_TS))
    #Rate_Island_Width_Change.append(copy.deepcopy(Width_Change_Rate_TS))
    #Percent_Island_Width_Change.append(copy.deepcopy(Width_Percent_Change))
    c = 20
    return(Width_Percent_Change,All_Focused_Domain_Width_TS)

def Calculate_Roadway_Abandonmet(cascade, years_modeled, buffer_length, number_barrier3d_models):
    # Find times the roadway broke and save the year that it did
    Road_Drowning_Years = []
    Road_Drowning_Reason = []
    for m in range(buffer_length, (number_barrier3d_models - buffer_length - 1)):
        Road_Data = cascade.roadways[m]

        # Record reason for road drowning
        if Road_Data.relocation_break == 1:
            Road_Drowning_Reason.append(copy.deepcopy('No relocation room'))
        elif Road_Data.drown_break == True:
            Road_Drowning_Reason.append(copy.deepcopy('Too much water'))
        elif cascade.b3d_break == 1:
            Road_Drowning_Reason.append(copy.deepcopy('Island Drowning'))
        else:
            Road_Drowning_Reason.append(copy.deepcopy('NA'))

        # Record the year of roadway abandonment
        if Road_Data.relocation_break == 1 or Road_Data.drown_break == True:
            Road_Drowning_Years.append(copy.deepcopy(Road_Data.time_index))
        elif Road_Data.drown_break == int(0) and Road_Data.time_index == 1:
            Road_Drowning_Years.append(copy.deepcopy(1))
        else:
            Road_Drowning_Years.append(copy.deepcopy(years_modeled))
        x = 10
    return(Road_Drowning_Years,Road_Drowning_Reason)

def Calculate_Island_Elevation_Metrics(cascade, buffer_length, number_barrier3d_models):
    Domains_of_Interest = range(buffer_length, (number_barrier3d_models - buffer_length))
    b3d = cascade.barrier3d

    Elevation_Change_TS = []
    Initial_Elevation_TS = []
    Final_Elevation_TS = []

    for k in Domains_of_Interest:
        Temp_B3D = b3d[k]
        Shoreline_Change_TS = Temp_B3D._ShorelineChangeTS

        # Find the initial and final island interior elevations
        initial_elev = Temp_B3D.DomainTS[0]
        final_elev = Temp_B3D.DomainTS[Temp_B3D.time_index-1]

        # Find the initial dune crest elevation
        initial_dune = Temp_B3D.DuneDomain[0]
        final_dune_elev = Temp_B3D.DuneDomain[Temp_B3D.time_index-1]

        # Rotate np arrays to concat with
        initial_dune_elev_r = np.rot90(initial_dune,k=3)
        final_dune_elev_r = np.rot90(final_dune_elev,k=3)

        # Add dunes to front of initial elevation
        combined_initial_elev = np.concatenate((initial_dune_elev_r,initial_elev),axis=0)

        distance_traveled = int(abs(np.sum(Shoreline_Change_TS[0:int(Temp_B3D.time_index)])))

        blank_cells = np.full((distance_traveled,int(Temp_B3D.BarrierLength)),-0.3)

        updated_final_array = np.concatenate((np.concatenate((blank_cells,final_dune_elev_r),axis=0),final_elev),axis=0)

        if len(combined_initial_elev) > len(updated_final_array):
            len_difference = len(combined_initial_elev) - len(updated_final_array)
            extra_cells = np.full((len_difference,int(Temp_B3D.BarrierLength)),-0.3)
            updated_final_array = np.concatenate((updated_final_array,extra_cells),axis=0)
        elif len(combined_initial_elev) < len(updated_final_array):
            len_difference = len(updated_final_array) - len(combined_initial_elev)
            extra_cells = np.full((len_difference,int(Temp_B3D.BarrierLength)),-0.3)
            combined_initial_elev = np.concatenate((combined_initial_elev,extra_cells),axis=0)

        elev_difference = np.subtract(updated_final_array,combined_initial_elev)

        Elevation_Change_TS.append(copy.deepcopy(elev_difference))
        Initial_Elevation_TS.append(copy.deepcopy(combined_initial_elev))
        Final_Elevation_TS.append(copy.deepcopy(updated_final_array))

    return (Elevation_Change_TS,Initial_Elevation_TS,Final_Elevation_TS)

def Calculate_Subaerial_Elevation(cascade, buffer_length, number_barrier3d_models):
    Domains_of_Interest = range(buffer_length, (number_barrier3d_models - buffer_length))
    b3d = cascade.barrier3d

    Subaerial_TS = []

    for k in Domains_of_Interest:
        Temp_B3D = b3d[k]
        Domain_TS = Temp_B3D.DomainTS
        Domain_All_Year_TS = []
        for years in range(len(Domain_TS)):
            Focus_Domain = Domain_TS[years]
            if np.any(Focus_Domain):
                Subaerial_Elevs_Year_X = np.round(Focus_Domain[np.where(Focus_Domain > 0)],decimals=4)
                Domain_All_Year_TS.append(copy.deepcopy(Subaerial_Elevs_Year_X))
        Subaerial_TS.append(copy.deepcopy(Domain_All_Year_TS))
    z = 20
    return (Subaerial_TS)

def Calculate_Bay_Distance(cascade, buffer_length, number_barrier3d_models):
    Domains_of_Interest = range(buffer_length, (number_barrier3d_models - buffer_length))
    b3d = cascade.barrier3d

    Domain_TS = []

    for k in Domains_of_Interest:
        Temp_B3D = b3d[k]
        Shoreline_Change_TS = Temp_B3D._ShorelineChangeTS
        Mean_Distance = []
        for years in range((Temp_B3D.time_index-1)):
            Yearly_Vals = []
            raw_elev = Temp_B3D.DomainTS[years]
            flipped_array = copy.deepcopy(np.flip(raw_elev))
            for cols in range((raw_elev.shape[1])):
                Zero_Row_Data = (np.where(raw_elev[:,cols]>0)[0])
                if np.any(Zero_Row_Data):
                    Zero_Data_Index = Zero_Row_Data.max()
                else:
                    Zero_Data_Index = 0
                if years == 0:
                    Shoreline_Change_Dif = 0
                else:
                    Shoreline_Change_Dif = np.sum(Shoreline_Change_TS[:years])

                Zero_Data_Index_Final = Zero_Data_Index - Shoreline_Change_Dif
                Yearly_Vals.append(copy.deepcopy(Zero_Data_Index_Final))
            Mean_Distance.append(copy.deepcopy(np.mean(Yearly_Vals)))
        Domain_TS.append(copy.deepcopy(Mean_Distance))

    return (Domain_TS)


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

def Calculate_Overwash_Volume(cascade, years_modeled, buffer_length, number_barrier3d_models):
    All_OW_Year_TS_Temp = []
    Yearly_OW_TS = []
    for klm in range(buffer_length, (number_barrier3d_models - buffer_length -1)):
        b3d = cascade.barrier3d[klm]
        Total_OW = np.multiply(np.sum(b3d.QowTS),500)
        Overwash_TS = np.multiply((b3d.QowTS),500)
        All_OW_Year_TS_Temp.append(copy.deepcopy(Total_OW))
        Yearly_OW_TS.append(copy.deepcopy(Overwash_TS))
    return(All_OW_Year_TS_Temp,Yearly_OW_TS)
'''
for sinks in range(len(Sink_Name)):
    Output_DF = Process_Batch(Base_Name=Base_Name_List[7],
                          Sink_Name=Sink_Name[sinks],
                          Save_Path=Save_Path,
                          Storm_Scenario=Storm_Level)
    Output_DF = Process_Batch(Base_Name=Base_Name_List[8],
                          Sink_Name=Sink_Name[sinks],
                          Save_Path=Save_Path,
                          Storm_Scenario=Storm_Level)

'''
for base in range(len(Base_Name_List)):
    for sinks in range(len(Sink_Name)):
        Output_DF = Process_Batch(Base_Name = Base_Name_List[base],
                                  Sink_Name =Sink_Name[sinks],
                                  Save_Path = Save_Path,
                                  Storm_Scenario=Storm_Level)


print('Hello')