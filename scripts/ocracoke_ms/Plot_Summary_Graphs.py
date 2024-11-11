# Plot graphs of summary plots from Future model runs
import copy
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# List the summary files to load
Summary_Values_Path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Model Runs\\Summary_Values\\'
file_names = os.listdir(Summary_Values_Path)

Run_Name_List = []
Island_Interior_Change_List = []
Roadway_Abandonment_List = []
Shoreline_Change_Rate_List = []
Number_Roadway_Relocation_List = []
Frequency_of_Roadway_Relocation_List = []
Number_Sandbag_Emplacements = []
Average_Sandbag_Duration = []
Island_Interior_Change_List = []
Island_Drown_Domain_List = []
Island_Drown_Section_List = []
Island_Drown_Year_List = []
Island_Drown_Dict = {}
Island_Section_Dict = {}

for files in range(len(file_names)):
    full_name = Summary_Values_Path+file_names[files]
    data = pd.read_csv(full_name)
    data = data.drop('Unnamed: 0', axis=1)

    start_index = (file_names[files]).index('_')+1
    end_index = (file_names[files]).index('.')

    plot_name = file_names[files][start_index:end_index]
    Run_Name_List.append(copy.deepcopy(plot_name))

    # Island Interior Width
    Island_Interior_Change_List.append(copy.deepcopy(list(data['Island_Interior_Change'])))

    # Roadway Abandonment
    Roadway_Abandonment_List.append(copy.deepcopy(list(data['Mean_Roadway_Abandonment'])))

    # Shoreline Change Rate
    Shoreline_Change_Rate_List.append(copy.deepcopy(list(data['Mean_Shoreline_Change_Rate'])))

    # Roadway Relocation Cells
    Number_Roadway_Relocation_List.append(copy.deepcopy(list(data['Roadway_Relocations'])))

    # Frequency of Roadway Relocation
    Frequency_of_Roadway_Relocation_List.append(copy.deepcopy(list(data['Roadway_Relocation_Frequency'])))

    # Number of Sandbag Emplacements
    Number_Sandbag_Emplacements.append(copy.deepcopy(list(data['Number_of_Sandbag_Emplacements'])))

    # Sandbag Duration
    Average_Sandbag_Duration.append(copy.deepcopy(list(data['Sandbag_Duration'])))

    # Island Drown Domains
    Island_Drown_Domain_List.append(copy.deepcopy(int(list(data['Island_Drown_Domain'])[0])))

    # Island Drown Section
    Island_Drown_Section_List.append(copy.deepcopy(int(list(data['Island_Drown_Section'])[0])))
    Island_Section_Dict[copy.deepcopy(plot_name)] = copy.deepcopy(int(list(data['Island_Drown_Section'])[0]))

    # Island Drown Year
    Island_Drown_Year_List.append(copy.deepcopy(int(list(data['Island_Drown_Year'])[0])))
    Island_Drown_Dict[copy.deepcopy(plot_name)] = copy.deepcopy(int(list(data['Island_Drown_Year'])[0]))

# Set the graph values
All_Domain_Values = [Island_Interior_Change_List,
    Roadway_Abandonment_List,
    Shoreline_Change_Rate_List,
    Number_Roadway_Relocation_List,
    Frequency_of_Roadway_Relocation_List,
    Number_Sandbag_Emplacements,
    Average_Sandbag_Duration]

All_Domain_Names = ['Island Interior Change',
    'Roadway Abandonment',
    'Shoreline Change Rate',
    'Number Roadway Relocation',
    'Frequency of Roadway Relocation',
    'Number Sandbag Emplacements',
    'Average Sandbag Duration']

All_Domain_Variable_Axes = ['Island Interior Change (% Change)',
    'Roadway Abandonment (yr)',
    'Shoreline Change Rate (m/yr)',
    'Number Roadway Relocations',
    'Frequency of Roadway Relocation',
    'Number Sandbag Emplacements',
    'Average Sandbag Duration (yr)']

All_Year_Values = [Island_Drown_Domain_List,
    Island_Drown_Section_List,
    Island_Drown_Year_List]

domain_nums = range(11,50)

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

pairs = [
    [0,2],
    [1,3],
    [4,6],
    [5,7],
    [8,10],
    [9,11]
]

for graphs in range(len(pairs)):
    Num_1 = pairs[graphs][0]
    Num_2 = pairs[graphs][1]
    ax = plt.gca()
    ax.set_ylim([-3.5, 1])
    plt.plot(domain_nums,Shoreline_Change_Rate_List[Num_1],label=Run_Name_List[Num_1])
    plt.plot(domain_nums,Shoreline_Change_Rate_List[Num_2],label=Run_Name_List[Num_2])
    plt.title(All_Domain_Names[2])
    plt.xlabel('B3D Domain')
    plt.ylabel(All_Domain_Variable_Axes[2])
    plt.tight_layout()
    plt.legend()
    plt.show()


print('Hello Purr')

'''
for graphs in range(0,len(All_Domain_Values)):
    print(All_Domain_Names[graphs])
    for run in range(len(Run_Name_List)):
        plt.plot(domain_nums,All_Domain_Values[graphs][run],label=Run_Name_List[run])
    plt.title(All_Domain_Names[graphs])
    plt.xlabel('B3D Domain')
    plt.ylabel(All_Domain_Variable_Axes[graphs])
    plt.tight_layout()
    plt.legend()
    plt.show()


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
'''

x = 10