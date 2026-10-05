# Load data and display island widths
import copy
import pickle

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


base_data_path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Island Rollover Calculations\\'
save_path ='C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Island Rollover Calculations\\'

initial_island_widths = pd.read_csv('C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Island Rollover Calculations\\Initial_Island_Widths_Focus.csv')
initial_island_vals = np.concatenate(initial_island_widths.values)*-1

Load_Names = [
              'Data_Overwash_T.pkl', #####
              'All_Domains_Mean_Shoreline_Change_Rate_T.pkl', #####
              'Mean_Interior_Values_T.pkl', #############
              ]

name_0 = base_data_path+Load_Names[0]
name_1 = base_data_path+Load_Names[1]
name_2 = base_data_path+Load_Names[2]

# Load the data
with open(name_0, 'rb') as file_0:
    OW_Data = pickle.load(file_0)

with open(name_1, 'rb') as file_1:
   All_Ocean_Shoreline_Data = pickle.load(file_1)


with open(name_2, 'rb') as file_2:
   All_Interior_Width_Dict = pickle.load(file_2)

z = 10
Save_Names = [
    ['IL_SQ_E','IL_RR_E','IL_BN_E'], ['IL_SQ_A','IL_RR_A','IL_BN_A'],
    ['I_SQ_E','I_RR_E','I_BN_E'], ['I_SQ_A','I_RR_A','I_BN_A'],
    ['IH_SQ_E', 'IH_RR_E', 'IH_BN_E'],['IH_SQ_A','IH_RR_A','IH_BN_A']
]

RSLR_Rates = ['IL',
              'IL',
              'I',
              'I',
              'IH',
              'IH']

Inlet =['E','A','E','A','E','A']
c = 20

def Plot_Ocean_Shoreline_Interior_Width_Mean(Shoreline_Data,
                             Interior_Width,
                             Initial_Width,
                             Run_Info,
                             save_plot=False):
    z = 20
    Ocean_Year_25_TS = []
    Ocean_Year_50_TS = []
    Ocean_Year_75_TS = []
    Ocean_Year_100_TS = []

    Bay_Year_25_TS = []
    Bay_Year_50_TS = []
    Bay_Year_75_TS = []
    Bay_Year_100_TS = []

    Drown_Bay_TS = []
    Drown_Shoreline_TS = []
    Drown_Early = False

    Bay_Color = 'black'

    if Run_Info == 'IL_SQ_E_S' or Run_Info == 'IL_SQ_A_S' or Run_Info == 'I_SQ_E_S' or Run_Info == 'I_SQ_A_S' or Run_Info == 'IH_SQ_E_S' or Run_Info == 'IH_SQ_A_S' or Run_Info == 'IL_SQ_E_S_T' or Run_Info == 'IL_SQ_A_S_T' or Run_Info == 'I_SQ_E_S_T' or Run_Info == 'I_SQ_A_S_T' or Run_Info == 'IH_SQ_E_S_T' or Run_Info == 'IH_SQ_A_S_T':
        Bay_Color = '#1b9e77'

    if Run_Info == 'IL_BN_E_S' or Run_Info == 'IL_BN_A_S' or Run_Info == 'I_BN_E_S' or Run_Info == 'I_BN_A_S' or Run_Info == 'IH_BN_E_S' or Run_Info == 'IH_BN_A_S' or Run_Info == 'IL_BN_E_S_T' or Run_Info == 'IL_BN_A_S_T' or Run_Info == 'I_BN_E_S_T' or Run_Info == 'I_BN_A_S_T' or Run_Info == 'IH_BN_E_S_T' or Run_Info == 'IH_BN_A_S_T':
        Bay_Color = '#d95f02'

    if Run_Info == 'IL_RR_E_S' or Run_Info == 'IL_RR_A_S' or Run_Info == 'I_RR_E_S' or Run_Info == 'I_RR_A_S' or Run_Info == 'IH_RR_E_S' or Run_Info == 'IH_RR_A_S' or Run_Info == 'IL_RR_E_S_T' or Run_Info == 'IL_RR_A_S_T' or Run_Info == 'I_RR_E_S_T' or Run_Info == 'I_RR_A_S_T' or Run_Info == 'IH_RR_E_S_T' or Run_Info == 'IH_RR_A_S_T':
        Bay_Color = '#7570b3'


    Model_Run_Length_Ts = []
    Test_Length = Shoreline_Data[0]

    for years in range(len(Shoreline_Data[0])):
        Model_Run_Length_Ts.append(copy.deepcopy(len(Shoreline_Data[0][years])))
    Run_Lengths = np.sort(Model_Run_Length_Ts)
    Final_Year = Run_Lengths[59]+2
    if Final_Year <= 100:
        Drown_Early = True
        Focus_Version = 4

    for x in range(len(Shoreline_Data)):
        focus_data = Shoreline_Data[x]
        Ocean_Year_25_TS_Temp = []
        Ocean_Year_50_TS_Temp = []
        Ocean_Year_75_TS_Temp = []
        Ocean_Year_100_TS_Temp = []
        Ocean_Year_Drown_TS_Temp = []

        focus_bay_data = Interior_Width[x]
        Bay_Year_25_TS_Temp = []
        Bay_Year_50_TS_Temp = []
        Bay_Year_75_TS_Temp = []
        Bay_Year_100_TS_Temp = []
        Bay_Year_Drown_TS_Temp = []

        # Final_Ocean = Shoreline_Data[x][Focus_Version][-1]
        # Final_Bay = Bayshoreline_Data[x]

        for runs in range(len(focus_data)):
            focus_run = focus_data[runs]
            if len(focus_run) >= 25:
                Ocean_Year_25_TS_Temp.append(copy.deepcopy(focus_run[24]))
            if len(focus_run) >= 50:
                Ocean_Year_50_TS_Temp.append(copy.deepcopy(focus_run[49]))
            if len(focus_run) >= 75:
                Ocean_Year_75_TS_Temp.append(copy.deepcopy(focus_run[74]))
            if len(focus_run) >= 100:
                Ocean_Year_100_TS_Temp.append(copy.deepcopy(focus_run[99]))
            if Drown_Early and len(focus_run)-1 >= Final_Year:
                Ocean_Year_Drown_TS_Temp.append(copy.deepcopy(focus_run[Final_Year]))

        if len(focus_bay_data) >= 25:
            Bay_Year_25_TS_Temp.append(copy.deepcopy(focus_bay_data[24]))
        if len(focus_bay_data) >= 50:
            Bay_Year_50_TS_Temp.append(copy.deepcopy(focus_bay_data[49]))
        if len(focus_bay_data) >= 75:
            Bay_Year_75_TS_Temp.append(copy.deepcopy(focus_bay_data[74]))
        if len(focus_bay_data) >= 100:
            Bay_Year_100_TS_Temp.append(copy.deepcopy(focus_bay_data[99]))
        if Drown_Early and len(focus_bay_data) >= Final_Year:
            Bay_Year_Drown_TS_Temp.append(copy.deepcopy(focus_bay_data[Final_Year]))

        if len(Ocean_Year_25_TS_Temp) >= 20:
            Mean_Ocean_Year_25 = np.mean(Ocean_Year_25_TS_Temp)
            Ocean_Year_25_TS.append(copy.deepcopy(Mean_Ocean_Year_25))
            Mean_Bay_Year_25 = np.mean(Bay_Year_25_TS_Temp)
            Bay_Year_25_TS.append(copy.deepcopy(Mean_Bay_Year_25))
        if len(Ocean_Year_50_TS_Temp) >= 20:
            Mean_Ocean_Year_50 = np.mean(Ocean_Year_50_TS_Temp)
            Ocean_Year_50_TS.append(copy.deepcopy(Mean_Ocean_Year_50))
            Mean_Bay_Year_50 = np.mean(Bay_Year_50_TS_Temp)
            Bay_Year_50_TS.append(copy.deepcopy(Mean_Bay_Year_50))
        if len(Ocean_Year_75_TS_Temp) >= 20:
            Mean_Ocean_Year_75 = np.mean(Ocean_Year_75_TS_Temp)
            Ocean_Year_75_TS.append(copy.deepcopy(Mean_Ocean_Year_75))
            Mean_Bay_Year_75 = np.mean(Bay_Year_75_TS_Temp)
            Bay_Year_75_TS.append(copy.deepcopy(Mean_Bay_Year_75))
        if len(Ocean_Year_100_TS_Temp) >= 20:
            Mean_Ocean_Year_100 = np.mean(Ocean_Year_100_TS_Temp)
            Ocean_Year_100_TS.append(copy.deepcopy(Mean_Ocean_Year_100))
            Mean_Bay_Year_100 = np.mean(Bay_Year_100_TS_Temp)
            Bay_Year_100_TS.append(copy.deepcopy(Mean_Bay_Year_100))

        Drown_Bay_TS.append(copy.deepcopy(np.mean(Bay_Year_Drown_TS_Temp)))
        Drown_Shoreline_TS.append(copy.deepcopy(np.mean(Ocean_Year_Drown_TS_Temp)))

    Domains_of_Interest = range(31, 50)
    if len(Ocean_Year_25_TS) == 19:
        plt.plot(Domains_of_Interest, Ocean_Year_25_TS, linestyle='solid', color='black', alpha=0.8)
    if len(Ocean_Year_50_TS) == 19:
        plt.plot(Domains_of_Interest, Ocean_Year_50_TS, linestyle='dashed', color='black', alpha=0.8)
    if len(Ocean_Year_75_TS) == 19:
        plt.plot(Domains_of_Interest, Ocean_Year_75_TS, linestyle='dashdot', color='black', alpha=0.8)
    if len(Ocean_Year_100_TS) == 19:
        plt.plot(Domains_of_Interest, Ocean_Year_100_TS, linestyle='dotted', color='black', alpha=0.8)
    if len(Ocean_Year_25_TS) == 19:
        plot_25_bay = np.subtract(Ocean_Year_25_TS, np.multiply(Bay_Year_25_TS,10))
        plt.plot(Domains_of_Interest, plot_25_bay, linestyle='solid', color=Bay_Color, alpha=0.8)

    if len(Ocean_Year_50_TS) == 19:
        plot_50_bay = np.subtract(Ocean_Year_50_TS, (np.multiply(Bay_Year_50_TS,10)))
        # plt.plot(Domains_of_Interest,plot_50_bay,linestyle='dotted')
        plt.plot(Domains_of_Interest, plot_50_bay, linestyle='dashed', color=Bay_Color, alpha=0.8)

    if len(Ocean_Year_75_TS) == 19:
        plot_75_bay = np.subtract(Ocean_Year_75_TS, np.multiply(Bay_Year_75_TS,10))
        # plt.plot(Domains_of_Interest,plot_75_bay,linestyle='dashdot')
        plt.plot(Domains_of_Interest, plot_75_bay, linestyle='dashdot', color=Bay_Color, alpha=0.8)
    if len(Ocean_Year_100_TS) == 19:
        plot_100_bay = np.subtract(Ocean_Year_100_TS, np.multiply(Bay_Year_100_TS,10))
        # plt.plot(Domains_of_Interest,plot_100_bay,linestyle='dashed')
        plt.plot(Domains_of_Interest, plot_100_bay, linestyle='dotted', color=Bay_Color, alpha=0.8)

    if Drown_Early:
        plot_final_bay = np.subtract(Drown_Shoreline_TS, np.multiply(Drown_Bay_TS,10))
        plt.plot(Domains_of_Interest, plot_final_bay, linestyle='dotted', color='blue')
        plt.plot(Domains_of_Interest, Drown_Shoreline_TS, linestyle='dotted', color='brown')
        plt.title(Run_Info + ' Drowned at ' + str(Final_Year))
    else:
        plt.title(Run_Info)

    plt.xlim(31, 49)
    plt.ylim(100, -750)
    # plt.title(Run_Info)
    plt.xlabel('Domain #')
    plt.ylabel('Position Relative to Initial Ocean Shoreline (m)')
    if save_plot == True:
        base = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Island Rollover Calculations\\Figures\\'
        plt.savefig(base + str(Run_Info) + '_Island_Width_M1.pdf', format='pdf')
        #plt.savefig(base + str(Run_Info) + '_Island_Width_M1.png', format='png')
    plt.show()



   Plot_Ocean_Shoreline_Interior_Width_Mean(Shoreline_Data= All_Ocean_Shoreline_Data[Save_Names[comb][0]],
                                       Interior_Width= All_Interior_Width_Dict[Save_Names[comb][0]],
                                       Initial_Width= initial_island_vals,
                                       Run_Info= Save_Names[comb][0],
                                       save_plot=True)
   
   Plot_Ocean_Shoreline_Interior_Width_Mean(Shoreline_Data= All_Ocean_Shoreline_Data[Save_Names[comb][1]],
                                       Interior_Width= All_Interior_Width_Dict[Save_Names[comb][1]],
                                       Initial_Width= initial_island_vals,
                                       Run_Info= Save_Names[comb][1],
                                       save_plot=True)

   Plot_Ocean_Shoreline_Interior_Width_Mean(Shoreline_Data= All_Ocean_Shoreline_Data[Save_Names[comb][2]],
                                       Interior_Width= All_Interior_Width_Dict[Save_Names[comb][2]],
                                       Initial_Width= initial_island_vals,
                                       Run_Info= Save_Names[comb][2],
                                       save_plot=True)

   t = 20

