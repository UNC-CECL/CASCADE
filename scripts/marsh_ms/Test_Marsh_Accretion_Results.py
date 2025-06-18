import copy

import matplotlib.pyplot as plt
import numpy as np

load_path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 2\\Updated_Marsh_Accretion_C_Release\\'
save_path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 2\\Cascade_CSV_Outputs\\Carbon Released\\'
plot_data = False

Storm_List = ['Baseline',
              '5',
              '10']

RSLR = ['IL',
        'I',
        'IH']

Base_Name_List = ['Geom_1',
                  'Geom_3',
                  'Geom_4',
                  'Geom_5']

for Geometries in range(len(Base_Name_List)):
    Base_Name = Base_Name_List[Geometries]
    for RSLR_Rates in range(len(RSLR)):
        for Storm_Mean_Intensities in range(len(Storm_List)):
            Save_Name = Base_Name+'_'+RSLR[RSLR_Rates]+'_'+Storm_List[Storm_Mean_Intensities]
            full_load_path = load_path+Save_Name

            Mean_Distance = np.load(full_load_path+'_Mean_Shoreline_Location.npy')
            Mean_C_Vals = np.load(full_load_path+'_Mean_C_Added.npy')

            Total_C_Release = []
            for years in range(len(Mean_Distance)):
                Yearly_C_Amount = Mean_C_Vals[years+1]
                Total_C_Eroded = np.sum(Yearly_C_Amount[int(Mean_Distance[years]):])
                Total_C_Release.append(copy.deepcopy(Total_C_Eroded))

            Total_C_Release_kg = np.divide(Total_C_Release,1000)
            np.savetxt(save_path+Save_Name+'_C_Release.csv',Total_C_Release_kg, delimiter=',')
            print('Saved '+str(Save_Name))
    if plot_data == True:
        plt.plot(Total_C_Release)
        plt.show()

        plt.plot(Mean_C_Vals[0])
        plt.plot(Mean_C_Vals[25])
        plt.plot(Mean_C_Vals[50])
        plt.plot(Mean_C_Vals[75])
        plt.plot(Mean_C_Vals[125])
        plt.axvline(x=Mean_Distance[0])
        plt.axvline(x=Mean_Distance[24])
        plt.axvline(x=Mean_Distance[49])
        plt.axvline(x=Mean_Distance[74])
        plt.axvline(x=Mean_Distance[99])
        plt.axvline(x=Mean_Distance[124])
        plt.xlim(4500,7000)

        plt.show()
    c = 20