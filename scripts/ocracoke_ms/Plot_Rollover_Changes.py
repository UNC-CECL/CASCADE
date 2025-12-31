# Load data and display as heat map
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
    ['IL_SQ_E_S_T','IL_RR_E_S_T','IL_BN_E_S_T'], ['IL_SQ_A_S_T','IL_RR_A_S_T','IL_BN_A_S_T'],
    ['I_SQ_E_S_T','I_RR_E_S_T','I_BN_E_S_T'], ['I_SQ_A_S_T','I_RR_A_S_T','I_BN_A_S_T'],
    ['IH_SQ_E_S_T', 'IH_RR_E_S_T', 'IH_BN_E_S_T'],['IH_SQ_A_S_T','IH_RR_A_S_T','IH_BN_A_S_T']
]

RSLR_Rates = ['IL',
              'IL',
              'I',
              'I',
              'IH',
              'IH']

Inlet =['E','A','E','A','E','A']
c = 20

#
# Save_Names = [
#     ['IL_BN_E_S', 'IL_BN_E_S_T'],
#     ['IL_BN_A_S','IL_BN_A_S_T'],
#     ['I_BN_E_S', 'I_BN_E_S_T'],
#     ['I_BN_A_S', 'I_BN_A_S_T'],
#     #['IH_SQ_E_S', 'IH_SQ_E_S_T'],
#     #['IH_SQ_A_S', 'IH_SQ_A_S_T']
# ]
#
# Save_Names = [['I_PR_E_S','I_PR_E_S_T']]
# RSLR_Rates = ['I']
# Inlet = ['E']

# Plot data
def Plot_Data_Bay_Shoreline_Changes(SQ_Data,
                                RR_Data,
                                BN_Data,
                                RSLR_Rate,
                                Inlet,
                                Save_Fig=False):
    Title = 'Erosion Hotspot Bay Shoreline Position '+str(RSLR_Rate) +' '+Inlet
    BIGGER_SIZE = 16

    plt.rc('font', size=BIGGER_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=16)  # fontsize of the axes title
    plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=BIGGER_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=BIGGER_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=12)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    ax = plt.gca()

    plt.plot(SQ_Data, color = '#1b9e77', label='Status Quo',linewidth = 2.5)
    plt.plot(BN_Data, color = '#d95f02', label='Nourishment', linewidth = 2.5,linestyle='dashdot')
    plt.plot(RR_Data, color = '#7570b3', label='Road Removal',linewidth = 2.5,linestyle='dashed')

    plt.xlabel('Year')
    plt.title(Title)
    plt.ylabel('Total Bay Shoreline Change (m)')
    plt.ylim(-350,50)
    plt.tight_layout()
    plt.legend(loc =3)
    if Save_Fig == True:
        plt.savefig(fname=('C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Island Rollover Calculations\\Figures\\'+Title+'.pdf'),format='pdf')
        plt.savefig(fname=('C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Island Rollover Calculations\\Figures\\'+Title+'.png'),format='png')
    plt.show()

def Plot_Data_Shoreline_Changes(SQ_Data,
                                RR_Data,
                                BN_Data,
                                RSLR_Rate,
                                Inlet,
                                Save_Fig=False):
    Title = 'Erosion Hotspot Shoreline Position '+str(RSLR_Rate) +' '+Inlet
    BIGGER_SIZE = 16

    plt.rc('font', size=BIGGER_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=16)  # fontsize of the axes title
    plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=BIGGER_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=BIGGER_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=12)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    ax = plt.gca()

    plt.plot(SQ_Data, color = '#1b9e77', label='Status Quo',linewidth = 2.5)
    plt.plot(BN_Data, color = '#d95f02', label='Nourishment', linewidth = 2.5,linestyle='dashdot')
    plt.plot(RR_Data, color = '#7570b3', label='Road Removal',linewidth = 2.5,linestyle='dashed')

    plt.xlabel('Year')
    plt.title(Title)
    plt.ylabel('Total Shoreline Change (m)')
    plt.ylim(-400,25)
    plt.tight_layout()
    plt.legend(loc =3)
    if Save_Fig == True:
        plt.savefig(fname=('C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Island Rollover Calculations\\Figures\\'+Title+'.pdf'),format='pdf')
        plt.savefig(fname=(
                'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Island Rollover Calculations\\Figures\\' + Title + '.png'),
                format='png')

    plt.show()

def Plot_Data_OW(SQ_Data,
                                RR_Data,
                                BN_Data,
                                RSLR_Rate,
                                Inlet,
                                Save_Fig=False):

        color_list = ['#1b9e77',
                      '#d95f02',
                      '#7570b3']

        All_Data = [SQ_Data,BN_Data,RR_Data]
        All_OW_25 = []
        All_OW_50 = []
        All_OW_75 = []
        All_OW_100 = []
        for j in range(len(All_Data)):
            Base_Data = All_Data[j]
            OW_25_TS = []
            OW_50_TS = []
            OW_75_TS = []
            OW_100_TS = []

            for x in range(len(Base_Data)):
                focus_data = Base_Data[x]
                OW_25_TS_Temp = []
                OW_50_TS_Temp = []
                OW_75_TS_Temp = []
                OW_100_TS_Temp = []

                for runs in range(len(focus_data)):
                    focus_run = focus_data[runs]
                    if len(focus_run) >= 25:
                        OW_25_TS_Temp.append(copy.deepcopy(focus_run[24]))
                    if len(focus_run) >= 50:
                        OW_50_TS_Temp.append(copy.deepcopy(focus_run[49]))
                    if len(focus_run) >= 75:
                        OW_75_TS_Temp.append(copy.deepcopy(focus_run[74]))
                    if len(focus_run) >= 100:
                        OW_100_TS_Temp.append(copy.deepcopy(focus_run[99]))

                if len(OW_25_TS_Temp) >= 20:
                    Mean_OW_Year_25 = np.mean(OW_25_TS_Temp)
                    if Mean_OW_Year_25 > 0:
                        OW_25_TS.append(copy.deepcopy(Mean_OW_Year_25))
                    else:
                        OW_25_TS.append(copy.deepcopy(0))
                if len(OW_50_TS_Temp) >= 20:
                    Mean_OW_Year_50 = np.mean(OW_50_TS_Temp)
                    if Mean_OW_Year_50 > 0:
                        OW_50_TS.append(copy.deepcopy(Mean_OW_Year_50))
                    else:
                        OW_50_TS.append(copy.deepcopy(0))
                if len(OW_75_TS_Temp) >= 20:
                    Mean_OW_Year_75 = np.mean(OW_75_TS_Temp)
                    if Mean_OW_Year_75 > 0:
                        OW_75_TS.append(copy.deepcopy(Mean_OW_Year_75))
                    else:
                        OW_75_TS.append(0)
                if len(OW_100_TS_Temp) >= 20:
                    Mean_OW_Year_100 = np.mean(OW_100_TS_Temp)
                    if Mean_OW_Year_100 > 0:
                        OW_100_TS.append(copy.deepcopy(Mean_OW_Year_100))
                    else:
                        OW_100_TS.append(0)

            All_OW_25.append(copy.deepcopy(OW_25_TS))
            All_OW_50.append(copy.deepcopy(OW_50_TS))
            All_OW_75.append(copy.deepcopy(OW_75_TS))
            All_OW_100.append(copy.deepcopy(OW_100_TS))

        Domains_of_Interest = range(31, 50)
        for k in range(len(All_Data)):
            if len(All_OW_25[k]) == 19:
                plt.plot(Domains_of_Interest, All_OW_25[k], linestyle='solid', color=color_list[k], alpha=0.8)
            if len(All_OW_50[k]) == 19:
                plt.plot(Domains_of_Interest, All_OW_50[k], linestyle='dashed', color=color_list[k], alpha=0.8)
            if len(All_OW_75[k]) == 19:
                plt.plot(Domains_of_Interest, All_OW_75[k], linestyle='dashdot', color=color_list[k], alpha=0.8)
            if len(All_OW_100[k]) == 19:
                plt.plot(Domains_of_Interest, All_OW_100[k], linestyle='dotted', color=color_list[k], alpha=0.8)


        plt.title(RSLR_Rate+' '+Inlet+' Total OW Flux (m^3)')

        plt.xlim(31, 49)
        plt.ylim(0,300000)
        plt.xlabel('Domain #')
        plt.ylabel('Total OW Flux (m^3)')
        if Save_Fig == True:
            base = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Island Rollover Calculations\\Figures\\'
            plt.savefig(base + str(RSLR_Rate) +'_'+Inlet+'_Mean_OW_Flux.pdf', format='pdf')
            plt.savefig(base + str(RSLR_Rate) +'_'+Inlet+ '_Mean_OW_Flux.png', format='png')
        plt.show()


def Plot_Data_OW_Test(SQ_Data,
                 SQ_Data_T,
                 RSLR_Rate,
                 Inlet,
                 Save_Fig=False):
    color_list = ['#1b9e77',
                  'blue',
                  '#7570b3']

    All_Data = [SQ_Data, SQ_Data_T]
    All_OW_25 = []
    All_OW_50 = []
    All_OW_75 = []
    All_OW_100 = []
    for j in range(len(All_Data)):
        Base_Data = All_Data[j]
        OW_25_TS = []
        OW_50_TS = []
        OW_75_TS = []
        OW_100_TS = []

        for x in range(len(Base_Data)):
            focus_data = Base_Data[x]
            OW_25_TS_Temp = []
            OW_50_TS_Temp = []
            OW_75_TS_Temp = []
            OW_100_TS_Temp = []

            for runs in range(len(focus_data)):
                focus_run = focus_data[runs]
                if len(focus_run) >= 25:
                    OW_25_TS_Temp.append(copy.deepcopy(focus_run[24]))
                if len(focus_run) >= 50:
                    OW_50_TS_Temp.append(copy.deepcopy(focus_run[49]))
                if len(focus_run) >= 75:
                    OW_75_TS_Temp.append(copy.deepcopy(focus_run[74]))
                if len(focus_run) >= 100:
                    OW_100_TS_Temp.append(copy.deepcopy(focus_run[99]))

            if len(OW_25_TS_Temp) >= 20:
                Mean_OW_Year_25 = np.mean(OW_25_TS_Temp)
                if Mean_OW_Year_25 > 0:
                    OW_25_TS.append(copy.deepcopy(Mean_OW_Year_25))
                else:
                    OW_25_TS.append(copy.deepcopy(0))
            if len(OW_50_TS_Temp) >= 20:
                Mean_OW_Year_50 = np.mean(OW_50_TS_Temp)
                if Mean_OW_Year_50 > 0:
                    OW_50_TS.append(copy.deepcopy(Mean_OW_Year_50))
                else:
                    OW_50_TS.append(copy.deepcopy(0))
            if len(OW_75_TS_Temp) >= 20:
                Mean_OW_Year_75 = np.mean(OW_75_TS_Temp)
                if Mean_OW_Year_75 > 0:
                    OW_75_TS.append(copy.deepcopy(Mean_OW_Year_75))
                else:
                    OW_75_TS.append(0)
            if len(OW_100_TS_Temp) >= 20:
                Mean_OW_Year_100 = np.mean(OW_100_TS_Temp)
                if Mean_OW_Year_100 > 0:
                    OW_100_TS.append(copy.deepcopy(Mean_OW_Year_100))
                else:
                    OW_100_TS.append(0)

        All_OW_25.append(copy.deepcopy(OW_25_TS))
        All_OW_50.append(copy.deepcopy(OW_50_TS))
        All_OW_75.append(copy.deepcopy(OW_75_TS))
        All_OW_100.append(copy.deepcopy(OW_100_TS))

    Domains_of_Interest = range(31, 50)
    for k in range(len(All_Data)):
        if len(All_OW_25[k]) == 19:
            plt.plot(Domains_of_Interest, All_OW_25[k], linestyle='solid', color=color_list[k], alpha=0.25)
        if len(All_OW_50[k]) == 19:
            plt.plot(Domains_of_Interest, All_OW_50[k], linestyle='dashed', color=color_list[k], alpha=0.25)
        if len(All_OW_75[k]) == 19:
            plt.plot(Domains_of_Interest, All_OW_75[k], linestyle='dashdot', color=color_list[k], alpha=0.25)
        if len(All_OW_100[k]) == 19:
            plt.plot(Domains_of_Interest, All_OW_100[k], linestyle='dotted', color=color_list[k], alpha=0.25)

    plt.title(RSLR_Rate + ' ' + Inlet + ' Total OW Flux (m^3)')

    plt.xlim(31, 49)
    plt.ylim(0, 300000)
    plt.xlabel('Domain #')
    plt.ylabel('Total OW Flux (m^3)')
    if Save_Fig == True:
        base = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Island Rollover Calculations\\Figures\\'
        plt.savefig(base + str(RSLR_Rate) + '_' + Inlet + '_Mean_OW_Flux_PR.pdf', format='pdf')
        plt.savefig(base + str(RSLR_Rate) + '_' + Inlet + '_Mean_OW_Flux_PR.png', format='png')
    plt.show()




def Plot_Ocean_Bay_Shoreline(Shoreline_Data,
                             Bayshoreline_Data,
                             Initial_Width,
                             Run_Info,
                             save_plot = True):
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

    if Run_Info == 'IL_SQ_E_S' or Run_Info == 'IL_SQ_A_S' or Run_Info == 'I_SQ_E_S' or Run_Info == 'I_SQ_A_S' or Run_Info == 'IH_SQ_E_S' or Run_Info == 'IH_SQ_A_S':
        Bay_Color = '#1b9e77'
    if Run_Info == 'IL_BN_E_S' or Run_Info == 'IL_BN_A_S' or Run_Info == 'I_BN_E_S' or Run_Info == 'I_BN_A_S' or Run_Info == 'IH_BN_E_S' or Run_Info == 'IH_BN_A_S':
        Bay_Color = '#d95f02'
    if Run_Info == 'IL_RR_E_S' or Run_Info == 'IL_RR_A_S' or Run_Info == 'I_RR_E_S' or Run_Info == 'I_RR_A_S' or Run_Info == 'IH_RR_E_S' or Run_Info == 'IH_RR_A_S':
        Bay_Color = '#7570b3'

    Model_Run_Length_Ts = []
    Test_Length = Shoreline_Data[0]

    for years in range(len(Shoreline_Data[0])):
        Model_Run_Length_Ts.append(copy.deepcopy(len(Shoreline_Data[0][years])))
    Run_Lengths = np.sort(Model_Run_Length_Ts)
    Final_Year = Run_Lengths[59]
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

        focus_bay_data = Bayshoreline_Data[x]
        Bay_Year_25_TS_Temp = []
        Bay_Year_50_TS_Temp = []
        Bay_Year_75_TS_Temp = []
        Bay_Year_100_TS_Temp = []
        Bay_Year_Drown_TS_Temp = []

        #Final_Ocean = Shoreline_Data[x][Focus_Version][-1]
        #Final_Bay = Bayshoreline_Data[x]

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
            if Drown_Early and len(focus_run) >= Final_Year:
                Ocean_Year_Drown_TS_Temp.append(copy.deepcopy(focus_run[Final_Year-1]))



        if len(focus_bay_data) >= 25:
            Bay_Year_25_TS_Temp.append(copy.deepcopy(focus_bay_data[24]))
        if len(focus_bay_data) >= 50:
            Bay_Year_50_TS_Temp.append(copy.deepcopy(focus_bay_data[49]))
        if len(focus_bay_data) >= 75:
            Bay_Year_75_TS_Temp.append(copy.deepcopy(focus_bay_data[74]))
        if len(focus_bay_data) >= 100:
            Bay_Year_100_TS_Temp.append(copy.deepcopy(focus_bay_data[99]))
        if Drown_Early and len(focus_bay_data)>= Final_Year:
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

    Domains_of_Interest = range(31,50)
    if len(Ocean_Year_25_TS) == 19:
        plt.plot(Domains_of_Interest,Ocean_Year_25_TS,linestyle='solid',color='black',alpha=0.8)
    if len(Ocean_Year_50_TS) == 19:
        plt.plot(Domains_of_Interest,Ocean_Year_50_TS,linestyle='dashed',color='black',alpha=0.8)
    if len(Ocean_Year_75_TS) == 19:
        plt.plot(Domains_of_Interest,Ocean_Year_75_TS,linestyle='dashdot',color='black',alpha=0.8)
    if len(Ocean_Year_100_TS) == 19:
        plt.plot(Domains_of_Interest,Ocean_Year_100_TS,linestyle='dotted',color='black',alpha=0.8)
    if len(Ocean_Year_25_TS) == 19:
        plot_25_bay = np.subtract(Initial_Width,Bay_Year_25_TS)
        plt.plot(Domains_of_Interest,plot_25_bay,linestyle='solid',color=Bay_Color,alpha=0.8)

    if len(Ocean_Year_50_TS) == 19:
        plot_50_bay = np.subtract(Initial_Width,Bay_Year_50_TS)
        #plt.plot(Domains_of_Interest,plot_50_bay,linestyle='dotted')
        plt.plot(Domains_of_Interest,plot_50_bay,linestyle='dashed',color=Bay_Color,alpha=0.8)

    if len(Ocean_Year_75_TS) == 19:
        plot_75_bay = np.subtract(Initial_Width,Bay_Year_75_TS)
        #plt.plot(Domains_of_Interest,plot_75_bay,linestyle='dashdot')
        plt.plot(Domains_of_Interest,plot_75_bay,linestyle='dashdot',color=Bay_Color,alpha=0.8)
    if len(Ocean_Year_100_TS) == 19:
        plot_100_bay = np.subtract(Initial_Width,Bay_Year_100_TS)
        #plt.plot(Domains_of_Interest,plot_100_bay,linestyle='dashed')
        plt.plot(Domains_of_Interest,plot_100_bay,linestyle='dotted',color=Bay_Color,alpha=0.8)

    if Drown_Early:
        plot_final_bay = np.subtract(Initial_Width,Drown_Bay_TS)
        plt.plot(Domains_of_Interest,plot_final_bay, linestyle = 'dotted',color='blue')
        plt.plot(Domains_of_Interest,Drown_Shoreline_TS, linestyle = 'dotted',color='brown')
        plt.title(Run_Info+' Drowned at '+str(Final_Year))
    else:
        plt.title(Run_Info)

    plt.axvline(x = 40)
    plt.axvline(x = 41)
    plt.axvline(x = 42)
    plt.axvline(x = 43)


    plt.xlim(31,49)
    plt.ylim(100,-750)
    #plt.title(Run_Info)
    plt.xlabel('Domain #')
    plt.ylabel('Position Relative to Initial Ocean Shoreline (m)')
    if save_plot == True:
        base = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Island Rollover Calculations\\Figures\\'
        plt.savefig(base+str(Run_Info)+'_Island_Width.pdf',format='pdf')
        plt.savefig(base+str(Run_Info)+'_Island_Width.png',format='png')
    plt.show()

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

def Plot_Ocean_Shoreline_Interior_Width_Median(Shoreline_Data,
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

    if Run_Info == 'IL_SQ_E_S' or Run_Info == 'IL_SQ_A_S' or Run_Info == 'I_SQ_E_S' or Run_Info == 'I_SQ_A_S' or Run_Info == 'IH_SQ_E_S' or Run_Info == 'IH_SQ_A_S' or Run_Info == 'IL_SQ_E_S_T' or Run_Info == 'IL_SQ_A_S_T' or Run_Info == 'I_SQ_E_S_T' or Run_Info == 'I_SQ_A_S_T' or Run_Info == 'IH_SQ_E_S_T' or Run_Info == 'IH_SQ_A_S_T':
        Bay_Color = '#1b9e77'
    if Run_Info == 'IL_BN_E_S' or Run_Info == 'IL_BN_A_S' or Run_Info == 'I_BN_E_S' or Run_Info == 'I_BN_A_S' or Run_Info == 'IH_BN_E_S' or Run_Info == 'IH_BN_A_S':
        Bay_Color = '#d95f02'
    if Run_Info == 'IL_RR_E_S' or Run_Info == 'IL_RR_A_S' or Run_Info == 'I_RR_E_S' or Run_Info == 'I_RR_A_S' or Run_Info == 'IH_RR_E_S' or Run_Info == 'IH_RR_A_S':
        Bay_Color = '#7570b3'

    Model_Run_Length_Ts = []
    Test_Length = Shoreline_Data[0]

    for years in range(len(Shoreline_Data[0])):
        Model_Run_Length_Ts.append(copy.deepcopy(len(Shoreline_Data[0][years])))
    Run_Lengths = np.sort(Model_Run_Length_Ts)
    Final_Year = Run_Lengths[59] + 2
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
            if Drown_Early and len(focus_run) - 1 >= Final_Year:
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
        plot_25_bay = np.subtract(Ocean_Year_25_TS, np.multiply(Bay_Year_25_TS, 10))
        plt.plot(Domains_of_Interest, plot_25_bay, linestyle='solid', color=Bay_Color, alpha=0.8)

    if len(Ocean_Year_50_TS) == 19:
        plot_50_bay = np.subtract(Ocean_Year_50_TS, (np.multiply(Bay_Year_50_TS, 10)))
        # plt.plot(Domains_of_Interest,plot_50_bay,linestyle='dotted')
        plt.plot(Domains_of_Interest, plot_50_bay, linestyle='dashed', color=Bay_Color, alpha=0.8)

    if len(Ocean_Year_75_TS) == 19:
        plot_75_bay = np.subtract(Ocean_Year_75_TS, np.multiply(Bay_Year_75_TS, 10))
        # plt.plot(Domains_of_Interest,plot_75_bay,linestyle='dashdot')
        plt.plot(Domains_of_Interest, plot_75_bay, linestyle='dashdot', color=Bay_Color, alpha=0.8)
    if len(Ocean_Year_100_TS) == 19:
        plot_100_bay = np.subtract(Ocean_Year_100_TS, np.multiply(Bay_Year_100_TS, 10))
        # plt.plot(Domains_of_Interest,plot_100_bay,linestyle='dashed')
        plt.plot(Domains_of_Interest, plot_100_bay, linestyle='dotted', color=Bay_Color, alpha=0.8)

    if Drown_Early:
        plot_final_bay = np.subtract(Drown_Shoreline_TS, np.multiply(Drown_Bay_TS, 10))
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
        plt.savefig(base + str(Run_Info) + '_Island_Width_Mean.pdf', format='pdf')
        #plt.savefig(base + str(Run_Info) + '_Island_Width_Mean.png', format='png')
    plt.show()

def Plot_Island_Width_Outlier(Shoreline_Data,
                                             Interior_Width,
                                             Run_Info,
                                             save_plot=False):
    # 10th percentile
    Ocean_Year_25_TS_10 = []
    Ocean_Year_50_TS_10 = []
    Ocean_Year_75_TS_10 = []
    Ocean_Year_100_TS_10 = []

    Bay_Year_25_TS_10 = []
    Bay_Year_50_TS_10 = []
    Bay_Year_75_TS_10 = []
    Bay_Year_100_TS_10 = []

    Drown_Bay_TS_10 = []
    Drown_Ocean_TS_10 = []

    # 90th percentile
    Ocean_Year_25_TS_90 = []
    Ocean_Year_50_TS_90 = []
    Ocean_Year_75_TS_90 = []
    Ocean_Year_100_TS_90 = []

    Bay_Year_25_TS_90 = []
    Bay_Year_50_TS_90 = []
    Bay_Year_75_TS_90 = []
    Bay_Year_100_TS_90 = []

    Drown_Bay_TS_90 = []
    Drown_Ocean_TS_90 = []
    Drown_Early = False



    if Run_Info == 'IL_SQ_E_S' or Run_Info == 'IL_SQ_A_S' or Run_Info == 'I_SQ_E_S' or Run_Info == 'I_SQ_A_S' or Run_Info == 'IH_SQ_E_S' or Run_Info == 'IH_SQ_A_S':
        Bay_Color = '#1b9e77'
    if Run_Info == 'IL_BN_E_S' or Run_Info == 'IL_BN_A_S' or Run_Info == 'I_BN_E_S' or Run_Info == 'I_BN_A_S' or Run_Info == 'IH_BN_E_S' or Run_Info == 'IH_BN_A_S':
        Bay_Color = '#d95f02'
    if Run_Info == 'IL_RR_E_S' or Run_Info == 'IL_RR_A_S' or Run_Info == 'I_RR_E_S' or Run_Info == 'I_RR_A_S' or Run_Info == 'IH_RR_E_S' or Run_Info == 'IH_RR_A_S':
        Bay_Color = '#7570b3'

    Model_Run_Length_Ts = []

    for years in range(len(Shoreline_Data[0])):
        Model_Run_Length_Ts.append(copy.deepcopy(len(Shoreline_Data[0][years])))
    Run_Lengths = np.sort(Model_Run_Length_Ts)
    Final_Year = Run_Lengths[59] + 2
    if Final_Year <= 100:
        Drown_Early = True

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
            if Drown_Early and len(focus_run) - 1 >= Final_Year:
                Ocean_Year_Drown_TS_Temp.append(copy.deepcopy(focus_run[Final_Year]))

        x = 20

        if len(focus_bay_data) >= 25:
            Bay_Year_25_TS_Temp= copy.deepcopy(focus_bay_data[24])
        if len(focus_bay_data) >= 50:
            Bay_Year_50_TS_Temp=(copy.deepcopy(focus_bay_data[49]))
        if len(focus_bay_data) >= 75:
            Bay_Year_75_TS_Temp=(copy.deepcopy(focus_bay_data[74]))
        if len(focus_bay_data) >= 100:
            Bay_Year_100_TS_Temp=(copy.deepcopy(focus_bay_data[99]))
        if Drown_Early and len(focus_run) - 1 >= Final_Year:
            Bay_Year_Drown_TS_Temp=(copy.deepcopy(focus_bay_data[Final_Year]))

        if len(Ocean_Year_25_TS_Temp) >= 20:
            x = 20
            Distance_Key_Dict = {}
            Island_Width_Dict = {}
            for k in range(len(Ocean_Year_25_TS_Temp)):
                Distance_Key_Dict[str(k)] = Ocean_Year_25_TS_Temp[k]
                Island_Width_Dict[str(k)] = Bay_Year_25_TS_Temp[k]
            Sorted_Shoreline_Change = dict(sorted(Distance_Key_Dict.items(), key=lambda item: item[1]))

            Keys = list(Sorted_Shoreline_Change.keys())
            Key_90th = Keys[0:10]
            Key_10th = Keys[-10:]

            Ocean_Shoreline_10_Vals = []
            Ocean_Shoreline_90_Vals = []
            Bay_Shoreline_10_Vals = []
            Bay_Shoreline_90_Vals = []

            for ten in range(len(Key_10th)):
                Ocean_Shoreline_10_Vals.append(copy.deepcopy(Distance_Key_Dict[Key_10th[ten]]))
                Ocean_Shoreline_90_Vals.append(copy.deepcopy(Distance_Key_Dict[Key_90th[ten]]))
                Bay_Shoreline_10_Vals.append(copy.deepcopy(Island_Width_Dict[Key_10th[ten]]))
                Bay_Shoreline_90_Vals.append(copy.deepcopy(Island_Width_Dict[Key_90th[ten]]))


            Ocean_Shoreline_10_Mean = np.mean(Ocean_Shoreline_10_Vals)
            Ocean_Shoreline_90_Mean = np.mean(Ocean_Shoreline_90_Vals)

            Bay_Shoreline_10_Mean = np.mean(Bay_Shoreline_10_Vals)
            Bay_Shoreline_90_Mean = np.mean(Bay_Shoreline_90_Vals)

            Ocean_Year_25_TS_10.append(copy.deepcopy(Ocean_Shoreline_10_Mean))
            Bay_Year_25_TS_10.append(copy.deepcopy(Bay_Shoreline_10_Mean))
            Ocean_Year_25_TS_90.append(copy.deepcopy(Ocean_Shoreline_90_Mean))
            Bay_Year_25_TS_90.append(copy.deepcopy(Bay_Shoreline_90_Mean))
        if len(Ocean_Year_50_TS_Temp) >= 20:
            x = 20
            Distance_Key_Dict = {}
            Island_Width_Dict = {}
            for k in range(len(Ocean_Year_50_TS_Temp)):
                Distance_Key_Dict[str(k)] = Ocean_Year_50_TS_Temp[k]
                Island_Width_Dict[str(k)] = Bay_Year_50_TS_Temp[k]
            Sorted_Shoreline_Change = dict(sorted(Distance_Key_Dict.items(), key=lambda item: item[1]))

            Keys = list(Sorted_Shoreline_Change.keys())
            Key_90th = Keys[0:10]
            Key_10th = Keys[-10:]

            Ocean_Shoreline_10_Vals = []
            Ocean_Shoreline_90_Vals = []
            Bay_Shoreline_10_Vals = []
            Bay_Shoreline_90_Vals = []

            for ten in range(len(Key_10th)):
                Ocean_Shoreline_10_Vals.append(copy.deepcopy(Distance_Key_Dict[Key_10th[ten]]))
                Ocean_Shoreline_90_Vals.append(copy.deepcopy(Distance_Key_Dict[Key_90th[ten]]))
                Bay_Shoreline_10_Vals.append(copy.deepcopy(Island_Width_Dict[Key_10th[ten]]))
                Bay_Shoreline_90_Vals.append(copy.deepcopy(Island_Width_Dict[Key_90th[ten]]))

            Ocean_Shoreline_10_Mean = np.mean(Ocean_Shoreline_10_Vals)
            Ocean_Shoreline_90_Mean = np.mean(Ocean_Shoreline_90_Vals)

            Bay_Shoreline_10_Mean = np.mean(Bay_Shoreline_10_Vals)
            Bay_Shoreline_90_Mean = np.mean(Bay_Shoreline_90_Vals)

            Ocean_Year_50_TS_10.append(copy.deepcopy(Ocean_Shoreline_10_Mean))
            Bay_Year_50_TS_10.append(copy.deepcopy(Bay_Shoreline_10_Mean))
            Ocean_Year_50_TS_90.append(copy.deepcopy(Ocean_Shoreline_90_Mean))
            Bay_Year_50_TS_90.append(copy.deepcopy(Bay_Shoreline_90_Mean))
        if len(Ocean_Year_75_TS_Temp) >= 20:
            Distance_Key_Dict = {}
            Island_Width_Dict = {}
            for k in range(len(Ocean_Year_75_TS_Temp)):
                Distance_Key_Dict[str(k)] = Ocean_Year_75_TS_Temp[k]
                Island_Width_Dict[str(k)] = Bay_Year_75_TS_Temp[k]
            Sorted_Shoreline_Change = dict(sorted(Distance_Key_Dict.items(), key=lambda item: item[1]))

            Keys = list(Sorted_Shoreline_Change.keys())
            Key_90th = Keys[0:10]
            Key_10th = Keys[-10:]

            Ocean_Shoreline_10_Vals = []
            Ocean_Shoreline_90_Vals = []
            Bay_Shoreline_10_Vals = []
            Bay_Shoreline_90_Vals = []

            for ten in range(len(Key_10th)):
                Ocean_Shoreline_10_Vals.append(copy.deepcopy(Distance_Key_Dict[Key_10th[ten]]))
                Ocean_Shoreline_90_Vals.append(copy.deepcopy(Distance_Key_Dict[Key_90th[ten]]))
                Bay_Shoreline_10_Vals.append(copy.deepcopy(Island_Width_Dict[Key_10th[ten]]))
                Bay_Shoreline_90_Vals.append(copy.deepcopy(Island_Width_Dict[Key_90th[ten]]))

            Ocean_Shoreline_10_Mean = np.mean(Ocean_Shoreline_10_Vals)
            Ocean_Shoreline_90_Mean = np.mean(Ocean_Shoreline_90_Vals)

            Bay_Shoreline_10_Mean = np.mean(Bay_Shoreline_10_Vals)
            Bay_Shoreline_90_Mean = np.mean(Bay_Shoreline_90_Vals)

            Ocean_Year_75_TS_10.append(copy.deepcopy(Ocean_Shoreline_10_Mean))
            Bay_Year_75_TS_10.append(copy.deepcopy(Bay_Shoreline_10_Mean))
            Ocean_Year_75_TS_90.append(copy.deepcopy(Ocean_Shoreline_90_Mean))
            Bay_Year_75_TS_90.append(copy.deepcopy(Bay_Shoreline_90_Mean))
        if len(Ocean_Year_100_TS_Temp) >= 20:
            Distance_Key_Dict = {}
            Island_Width_Dict = {}
            for k in range(len(Ocean_Year_100_TS_Temp)):
                Distance_Key_Dict[str(k)] = Ocean_Year_100_TS_Temp[k]
                Island_Width_Dict[str(k)] = Bay_Year_100_TS_Temp[k]
            Sorted_Shoreline_Change = dict(sorted(Distance_Key_Dict.items(), key=lambda item: item[1]))

            Keys = list(Sorted_Shoreline_Change.keys())
            Key_90th = Keys[0:10]
            Key_10th = Keys[-10:]

            Ocean_Shoreline_10_Vals = []
            Ocean_Shoreline_90_Vals = []
            Bay_Shoreline_10_Vals = []
            Bay_Shoreline_90_Vals = []

            for ten in range(len(Key_10th)):
                Ocean_Shoreline_10_Vals.append(copy.deepcopy(Distance_Key_Dict[Key_10th[ten]]))
                Ocean_Shoreline_90_Vals.append(copy.deepcopy(Distance_Key_Dict[Key_90th[ten]]))
                Bay_Shoreline_10_Vals.append(copy.deepcopy(Island_Width_Dict[Key_10th[ten]]))
                Bay_Shoreline_90_Vals.append(copy.deepcopy(Island_Width_Dict[Key_90th[ten]]))

            Ocean_Shoreline_10_Mean = np.mean(Ocean_Shoreline_10_Vals)
            Ocean_Shoreline_90_Mean = np.mean(Ocean_Shoreline_90_Vals)

            Bay_Shoreline_10_Mean = np.mean(Bay_Shoreline_10_Vals)
            Bay_Shoreline_90_Mean = np.mean(Bay_Shoreline_90_Vals)

            Ocean_Year_100_TS_10.append(copy.deepcopy(Ocean_Shoreline_10_Mean))
            Bay_Year_100_TS_10.append(copy.deepcopy(Bay_Shoreline_10_Mean))
            Ocean_Year_100_TS_90.append(copy.deepcopy(Ocean_Shoreline_90_Mean))
            Bay_Year_100_TS_90.append(copy.deepcopy(Bay_Shoreline_90_Mean))

        if Drown_Early:
            Distance_Key_Dict = {}
            Island_Width_Dict = {}
            for k in range(len(Ocean_Year_Drown_TS_Temp)):
                Distance_Key_Dict[str(k)] = Ocean_Year_Drown_TS_Temp[k]
                Island_Width_Dict[str(k)] = Bay_Year_Drown_TS_Temp[k]
            Sorted_Shoreline_Change = dict(sorted(Distance_Key_Dict.items(), key=lambda item: item[1]))

            Keys = list(Sorted_Shoreline_Change.keys())
            Key_90th = Keys[0:10]
            Key_10th = Keys[-10:]

            Ocean_Shoreline_10_Vals = []
            Ocean_Shoreline_90_Vals = []
            Bay_Shoreline_10_Vals = []
            Bay_Shoreline_90_Vals = []

            for ten in range(len(Key_10th)):
                Ocean_Shoreline_10_Vals.append(copy.deepcopy(Distance_Key_Dict[Key_10th[ten]]))
                Ocean_Shoreline_90_Vals.append(copy.deepcopy(Distance_Key_Dict[Key_90th[ten]]))
                Bay_Shoreline_10_Vals.append(copy.deepcopy(Island_Width_Dict[Key_10th[ten]]))
                Bay_Shoreline_90_Vals.append(copy.deepcopy(Island_Width_Dict[Key_90th[ten]]))

            Ocean_Shoreline_10_Mean = np.mean(Ocean_Shoreline_10_Vals)
            Ocean_Shoreline_90_Mean = np.mean(Ocean_Shoreline_90_Vals)

            Bay_Shoreline_10_Mean = np.mean(Bay_Shoreline_10_Vals)
            Bay_Shoreline_90_Mean = np.mean(Bay_Shoreline_90_Vals)

            Drown_Ocean_TS_10.append(copy.deepcopy(Ocean_Shoreline_10_Mean))
            Drown_Bay_TS_10.append(copy.deepcopy(Bay_Shoreline_10_Mean))
            Drown_Ocean_TS_90.append(copy.deepcopy(Ocean_Shoreline_90_Mean))
            Drown_Bay_TS_90.append(copy.deepcopy(Bay_Shoreline_90_Mean))


    Domains_of_Interest = range(21, 40)
    if len(Ocean_Year_25_TS_10) == 19:
        plt.plot(Domains_of_Interest, Ocean_Year_25_TS_10, linestyle='solid', color='black', alpha=0.8)
    if len(Ocean_Year_50_TS_10) == 19:
        plt.plot(Domains_of_Interest, Ocean_Year_50_TS_10, linestyle='dashed', color='black', alpha=0.8)
    if len(Ocean_Year_75_TS_10) == 19:
        plt.plot(Domains_of_Interest, Ocean_Year_75_TS_10, linestyle='dashdot', color='black', alpha=0.8)
    if len(Ocean_Year_100_TS_10) == 19:
        plt.plot(Domains_of_Interest, Ocean_Year_100_TS_10, linestyle='dotted', color='black', alpha=0.8)
    if len(Ocean_Year_25_TS_10) == 19:
        plot_25_bay_10 = np.subtract(Ocean_Year_25_TS_10, np.multiply(Bay_Year_25_TS_10, 10))
        plt.plot(Domains_of_Interest, plot_25_bay_10, linestyle='solid', color=Bay_Color, alpha=0.8)

    if len(Ocean_Year_50_TS_10) == 19:
        plot_50_bay_10 = np.subtract(Ocean_Year_50_TS_10, (np.multiply(Bay_Year_50_TS_10, 10)))
        plt.plot(Domains_of_Interest, plot_50_bay_10, linestyle='dashed', color=Bay_Color, alpha=0.8)

    if len(Ocean_Year_75_TS_10) == 19:
        plot_75_bay_10 = np.subtract(Ocean_Year_75_TS_10, np.multiply(Bay_Year_75_TS_10, 10))
        plt.plot(Domains_of_Interest, plot_75_bay_10, linestyle='dashdot', color=Bay_Color, alpha=0.8)
    if len(Ocean_Year_100_TS_10) == 19:
        plot_100_bay_10 = np.subtract(Ocean_Year_100_TS_10, np.multiply(Bay_Year_100_TS_10, 10))
        plt.plot(Domains_of_Interest, plot_100_bay_10, linestyle='dotted', color=Bay_Color, alpha=0.8)

    if Drown_Early:
        plot_final_bay_10 = np.subtract(Drown_Ocean_TS_10, np.multiply(Drown_Bay_TS_10, 10))
        plt.plot(Domains_of_Interest, plot_final_bay_10, linestyle='dotted', color='blue')
        plt.plot(Domains_of_Interest, Drown_Ocean_TS_10, linestyle='dotted', color='brown')
        plt.title(Run_Info + ' Drowned at ' + str(Final_Year)+' 10%')
    else:
        plt.title(Run_Info+' 10%')

    plt.xlim(21, 39)
    plt.ylim(100, -750)
    # plt.title(Run_Info)
    plt.xlabel('Domain #')
    plt.ylabel('Position Relative to Initial Ocean Shoreline (m)')
    if save_plot == True:
        base = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Island Rollover Calculations\\Figures\\'
        plt.savefig(base + str(Run_Info) + '_Island_10_Mean.pdf', format='pdf')
        # plt.savefig(base + str(Run_Info) + '_Island_Width_Mean.png', format='png')
    plt.show()

    if len(Ocean_Year_25_TS_90) == 19:
        plt.plot(Domains_of_Interest, Ocean_Year_25_TS_90, linestyle='solid', color='black', alpha=0.8)
    if len(Ocean_Year_50_TS_90) == 19:
        plt.plot(Domains_of_Interest, Ocean_Year_50_TS_90, linestyle='dashed', color='black', alpha=0.8)
    if len(Ocean_Year_75_TS_90) == 19:
        plt.plot(Domains_of_Interest, Ocean_Year_75_TS_90, linestyle='dashdot', color='black', alpha=0.8)
    if len(Ocean_Year_100_TS_90) == 19:
        plt.plot(Domains_of_Interest, Ocean_Year_100_TS_90, linestyle='dotted', color='black', alpha=0.8)
    if len(Ocean_Year_25_TS_90) == 19:
        plot_25_bay_90 = np.subtract(Ocean_Year_25_TS_90, np.multiply(Bay_Year_25_TS_90, 10))
        plt.plot(Domains_of_Interest, plot_25_bay_90, linestyle='solid', color=Bay_Color, alpha=0.8)

    if len(Ocean_Year_50_TS_90) == 19:
        plot_50_bay_90 = np.subtract(Ocean_Year_50_TS_90, (np.multiply(Bay_Year_50_TS_90, 10)))
        plt.plot(Domains_of_Interest, plot_50_bay_90, linestyle='dashed', color=Bay_Color, alpha=0.8)

    if len(Ocean_Year_75_TS_90) == 19:
        plot_75_bay_90 = np.subtract(Ocean_Year_75_TS_90, np.multiply(Bay_Year_75_TS_90, 10))
        plt.plot(Domains_of_Interest, plot_75_bay_90, linestyle='dashdot', color=Bay_Color, alpha=0.8)
    if len(Ocean_Year_100_TS_90) == 19:
        plot_100_bay_90 = np.subtract(Ocean_Year_100_TS_90, np.multiply(Bay_Year_100_TS_90, 10))
        plt.plot(Domains_of_Interest, plot_100_bay_90, linestyle='dotted', color=Bay_Color, alpha=0.8)

    if Drown_Early:
        plot_final_bay_90 = np.subtract(Drown_Ocean_TS_90, np.multiply(Drown_Bay_TS_90, 10))
        plt.plot(Domains_of_Interest, plot_final_bay_90, linestyle='dotted', color='blue')
        plt.plot(Domains_of_Interest, Drown_Ocean_TS_90, linestyle='dotted', color='brown')
        plt.title(Run_Info + ' Drowned at ' + str(Final_Year)+' 90%')
        #plt.title(Run_Info + ' Drowned at ' + str(Final_Year))
    else:
        plt.title(Run_Info + ' 90%')

    plt.xlim(21, 39)
    plt.ylim(100, -750)
    # plt.title(Run_Info)
    plt.xlabel('Domain #')
    plt.ylabel('Position Relative to Initial Ocean Shoreline (m)')
    if save_plot == True:
        base = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Island Rollover Calculations\\Figures\\'
        plt.savefig(base + str(Run_Info) + '_Island_90_Mean.pdf', format='pdf')
        # plt.savefig(base + str(Run_Info) + '_Island_Width_Mean.png', format='png')
    plt.show()

    x = 20


x = 20
for comb in range(len(Save_Names)):
   # Plot_Data_Shoreline_Changes(SQ_Data=Hotspot_Shoreline_Data[Save_Names[comb][0]],
   #                              RR_Data=Hotspot_Shoreline_Data[Save_Names[comb][1]],
   #                              BN_Data=Hotspot_Shoreline_Data[Save_Names[comb][2]],
   #                              RSLR_Rate=RSLR_Rates[comb],
   #                              Inlet=Inlet[comb],
   #                              Save_Fig=True)
   #
   # Plot_Data_Bay_Shoreline_Changes(SQ_Data=Hotspot_Bay_Data[Save_Names[comb][0]],
   #                              RR_Data=Hotspot_Bay_Data[Save_Names[comb][1]],
   #                              BN_Data=Hotspot_Bay_Data[Save_Names[comb][2]],
   #                              RSLR_Rate=RSLR_Rates[comb],
   #                              Inlet=Inlet[comb],
   #                              Save_Fig=True)

   # Plot_Data_OW(SQ_Data=Hotspot_OW_Data[Save_Names[comb][0]],
   #                              RR_Data=Hotspot_OW_Data[Save_Names[comb][1]],
   #                              BN_Data=Hotspot_OW_Data[Save_Names[comb][2]],
   #                              RSLR_Rate=RSLR_Rates[comb],
   #                              Inlet=Inlet[comb],
   #                              Save_Fig=True)


   # Plot_Ocean_Bay_Shoreline(Shoreline_Data = All_Ocean_Shoreline_Data[Save_Names[comb][0]],
   #                          Bayshoreline_Data = All_Bay_Shoreline_Data[Save_Names[comb][0]],
   #                          Initial_Width = initial_island_vals,
   #                          Run_Info = Save_Names[comb][0]
   #                          )
   #
   # Plot_Ocean_Bay_Shoreline(Shoreline_Data = All_Ocean_Shoreline_Data[Save_Names[comb][1]],
   #                          Bayshoreline_Data = All_Bay_Shoreline_Data[Save_Names[comb][1]],
   #                          Initial_Width = initial_island_vals,
   #                          Run_Info = Save_Names[comb][1]
   #                          )

   # Plot_Ocean_Bay_Shoreline(Shoreline_Data = All_Ocean_Shoreline_Data[Save_Names[comb][2]],
   #                          Bayshoreline_Data = All_Bay_Shoreline_Data[Save_Names[comb][2]],
   #                          Initial_Width = initial_island_vals,
   #                          Run_Info = Save_Names[comb][2]
   #                          )

   c = 20

   # Plot_Island_Width_Outlier(Shoreline_Data= All_Ocean_Shoreline_Data[Save_Names[comb][0]],
   #                                     Interior_Width= All_Island_Interior_Vals_Dict[Save_Names[comb][0]],
   #                                     Run_Info= Save_Names[comb][0],
   #                                     save_plot=True)
   #
   # Plot_Island_Width_Outlier(Shoreline_Data= All_Ocean_Shoreline_Data[Save_Names[comb][1]],
   #                                     Interior_Width= All_Island_Interior_Vals_Dict[Save_Names[comb][1]],
   #                                     Run_Info= Save_Names[comb][1],
   #                                     save_plot=True)

   # Plot_Island_Width_Outlier(Shoreline_Data= All_Ocean_Shoreline_Data[Save_Names[comb][2]],
   #                                     Interior_Width= All_Island_Interior_Vals_Dict[Save_Names[comb][2]],
   #                                     Run_Info= Save_Names[comb][2],
   #                                     save_plot=True)
   #
   Plot_Ocean_Shoreline_Interior_Width_Mean(Shoreline_Data= All_Ocean_Shoreline_Data[Save_Names[comb][0]],
                                       Interior_Width= All_Interior_Width_Dict[Save_Names[comb][0]],
                                       Initial_Width= initial_island_vals,
                                       Run_Info= Save_Names[comb][0],
                                       save_plot=True)
   #
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
   #
   # Plot_Ocean_Shoreline_Interior_Width_Mean(Shoreline_Data= All_Ocean_Shoreline_Data[Save_Names[comb][2]],
   #                                     Interior_Width= All_Interior_Width_Dict[Save_Names[comb][2]],
   #                                     Initial_Width= initial_island_vals,
   #                                     Run_Info= Save_Names[comb][2],
   #                                     save_plot=True)

   # Plot_Data_OW(SQ_Data= OW_Data[Save_Names[comb][0]],
   #              RR_Data= OW_Data[Save_Names[comb][1]],
   #              BN_Data= OW_Data[Save_Names[comb][2]],
   #              RSLR_Rate = RSLR_Rates[comb],
   #              Inlet = Inlet[comb],
   #              Save_Fig=True)
   t = 20
   Plot_Data_OW_Test(SQ_Data= OW_Data[Save_Names[comb][1]],
                SQ_Data_T= OW_Data[Save_Names[comb][1]],
                RSLR_Rate = RSLR_Rates[comb],
                Inlet = Inlet[comb],
                Save_Fig=False)





'''   Plot_Ocean_Shoreline_Interior_Width(Shoreline_Data= All_Ocean_Shoreline_Data[Save_Names[comb][0]],
                                       Interior_Width= Median_Interior_Width_Dict[Save_Names[comb][0]],
                                       Initial_Width= initial_island_vals,
                                       Run_Info= Save_Names[comb][0],
                                       save_plot=True)

   Plot_Ocean_Shoreline_Interior_Width(Shoreline_Data= All_Ocean_Shoreline_Data[Save_Names[comb][2]],
                                       Interior_Width= Median_Interior_Width_Dict[Save_Names[comb][2]],
                                       Initial_Width= initial_island_vals,
                                       Run_Info= Save_Names[comb][2],
                                       save_plot=True)

   Plot_Ocean_Shoreline_Interior_Width(Shoreline_Data= All_Ocean_Shoreline_Data[Save_Names[comb][1]],
                                       Interior_Width= Median_Interior_Width_Dict[Save_Names[comb][1]],
                                       Initial_Width= initial_island_vals,
                                       Run_Info= Save_Names[comb][1],
                                       save_plot=True)'''


x = 20