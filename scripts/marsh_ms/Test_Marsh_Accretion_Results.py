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


def Calculate_D1(Input_Data):
    D1_List = []
    for k in range(1, len(Input_Data)):
        Current_Value = Input_Data[k]
        Past_Value = Input_Data[k - 1]
        D1 = Current_Value - Past_Value
        D1_List.append(copy.deepcopy(D1))
    return (D1_List)


def Calculate_Rolling_Mean(Input_Data, Seperation):
    Mean_TS = []
    Test_TS = []
    Index_Vals = list(range(1, len(Input_Data), Seperation))
    Index_Vals.append(len(Input_Data))
    for h in range(1, len(Index_Vals)):
        Subset = Input_Data[Index_Vals[h - 1]:Index_Vals[h]]
        Temp_Mean = np.mean(Subset)
        Mean_TS.append(copy.deepcopy(Temp_Mean))
        Test_TS.append(copy.deepcopy(Subset))
    return (Mean_TS, Test_TS, Index_Vals)

for Geometries in range(len(Base_Name_List)):
    Base_Name = Base_Name_List[Geometries]
    for RSLR_Rates in range(len(RSLR)):
        for Storm_Mean_Intensities in range(len(Storm_List)):
            Save_Name = Base_Name+'_'+RSLR[RSLR_Rates]+'_'+Storm_List[Storm_Mean_Intensities]
            full_load_path = load_path+Save_Name

            Mean_Distance = np.load(full_load_path+'_Mean_Shoreline_Location.npy')
            Mean_C_Vals = np.load(full_load_path+'_Mean_C_Added_Marsh.npy')

            Total_C_Release = []
            for years in range(len(Mean_Distance)):
                Yearly_C_Amount = Mean_C_Vals[years]
                Total_C_Eroded = np.sum(Yearly_C_Amount[int(Mean_Distance[years]):])
                Total_C_Release.append(copy.deepcopy(Total_C_Eroded))

            Total_C_Release_kg = np.divide(Total_C_Release,1000)
            #np.savetxt(save_path+Save_Name+'_C_Release_Marsh_Update.csv',Total_C_Release_kg, delimiter=',')
            print('Saved '+str(Save_Name))

            C_D1 = Calculate_D1(list(Total_C_Release_kg))
            Distance_D1 = Calculate_D1(list(Mean_Distance))


            Rolling_Mean_C1, Test_TS_LR_BS, Index_Vals_LR_BS = Calculate_Rolling_Mean(C_D1, 10)
            Rolling_Mean_Distance, Test_TS_IR_5S, Index_Vals_IR_5S = Calculate_Rolling_Mean(Distance_D1, 10)

            Distance_Averaged_RM = np.divide(Rolling_Mean_C1,np.absolute(Rolling_Mean_Distance))

            plt.plot(Distance_Averaged_RM)
            plt.title(Save_Name)
            plt.show()

            if plot_data == True:
                #plt.plot(Total_C_Release)
               # plt.show()

                plt.plot(Mean_C_Vals[0])
                plt.plot(Mean_C_Vals[25])
                plt.plot(Mean_C_Vals[50])
                plt.plot(Mean_C_Vals[75])
                plt.plot(Mean_C_Vals[124])
                plt.axvline(x=Mean_Distance[0])
                plt.axvline(x=Mean_Distance[24])
                plt.axvline(x=Mean_Distance[49])
                plt.axvline(x=Mean_Distance[74])
                plt.axvline(x=Mean_Distance[99])
                plt.axvline(x=Mean_Distance[119])
                plt.axvline(x=Mean_Distance[124])
                plt.xlim(4900,5500)

                plt.show()
                c = 20