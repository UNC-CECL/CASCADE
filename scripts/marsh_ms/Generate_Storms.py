from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from barrier3d import Barrier3d, Barrier3dBmi
from barrier3d.tools import input_files
import os
import sys
cwd = os.getcwd()
data_path = '/Users/ceclmac/PycharmProjects/CASCADE/data/marsh_init_data'

storm_list = np.loadtxt(data_path+"/StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",delimiter=',')

RMAX = storm_list[:,2]
Mean_Rmax = np.mean(RMAX)
Mean_Rmax_10 = Mean_Rmax*.1
Mean_Rmax_02 = Mean_Rmax*.02
Mean_Rmax_05 = Mean_Rmax*.05

'''
sys.path
for i in range(0,50):
    base_storm_series = input_files.shift_storm_intensity(
        datadir=data_path,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        shift=0,  # shift the TWL distribution to the right (i.e., increase intensity), m NAVD88, typically [-0.15, 0.15]
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=1000,
        bPlot=False,
        bSave=False,
        output_filename="StormList_10kyrs_VCR_Berm1pt9m_Slope0pt04-lowIntensity",
    )
    np.save(file='/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/Chapter 2/Storms/Baseline/StormList_'+str(i)+'_baseline.npy',arr=base_storm_series)

for i in range(0,50):
    storm_series_shifted_10 = input_files.shift_storm_intensity(
        datadir=data_path,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        shift=Mean_Rmax_10,  # shift the TWL distribution to the right (i.e., increase intensity), m NAVD88, typically [-0.15, 0.15]
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=1000,
        bPlot=False,
        bSave=False,
        output_filename="StormList_10kyrs_VCR_Berm1pt9m_Slope0pt04-lowIntensity",
    )
    np.save(file='/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/Chapter 2/Storms/Ten_Percent_Increase/StormList_'+str(i)+'_10_percent_increase.npy',arr=storm_series_shifted_10)

for i in range(0,50):
    storm_series_shifted_02 = input_files.shift_storm_intensity(
        datadir=data_path,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        shift=Mean_Rmax_02,  # shift the TWL distribution to the right (i.e., increase intensity), m NAVD88, typically [-0.15, 0.15]
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=1000,
        bPlot=False,
        bSave=False,
        output_filename="StormList_10kyrs_VCR_Berm1pt9m_Slope0pt04-lowIntensity",
    )
    np.save(file='/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/Chapter 2/Storms/Two_Percent_Increase/StormList_'+str(i)+'_2_percent_increase.npy',arr=storm_series_shifted_02)
'''
for i in range(0,50):
    storm_series_shifted_05 = input_files.shift_storm_intensity(
        datadir=data_path,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        shift=Mean_Rmax_05,  # shift the TWL distribution to the right (i.e., increase intensity), m NAVD88, typically [-0.15, 0.15]
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=1000,
        bPlot=False,
        bSave=False,
        output_filename="StormList_10kyrs_VCR_Berm1pt9m_Slope0pt04-lowIntensity",
    )
    np.save(file='/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/Chapter 2/Storms/Five_Percent_Increase/StormList_'+str(i)+'_5_percent_increase.npy',arr=storm_series_shifted_05)