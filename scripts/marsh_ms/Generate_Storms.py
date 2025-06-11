from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from barrier3d import Barrier3d, Barrier3dBmi
from barrier3d.tools import input_files
import os
import sys
cwd = os.getcwd()
data_path = '/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/Chapter 2/Storms/Historic_Storms/'

storm_list = np.loadtxt(data_path+"VCRSeaStormMatrix_1980_2020.csv",delimiter=',')

RMAX = storm_list[:,2]
Mean_Rmax = np.mean(RMAX)
Mean_Rmax_10 = Mean_Rmax*.1
Mean_Rmax_02 = Mean_Rmax*.02
Mean_Rmax_05 = Mean_Rmax*.05


sys.path
for i in range(50,100):
    base_storm_series = input_files.shift_storm_intensity(
        datadir=data_path,
        storm_list_name="VCRSeaStormMatrix_1980_2020.csv",  # can by .py or .csv
        mean_yearly_storms=42.25,
        SD_yearly_storms=4,
        shift=0,  # shift the TWL distribution to the right (i.e., increase intensity), m NAVD88, typically [-0.15, 0.15]
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.03,  # m NAVD88, just used for plotting
        model_years=125,
        bPlot=False,
        bSave=False,
        output_filename="StormList_10kyrs_VCR_Berm1pt9m_Slope0pt04-lowIntensity",
    )
    np.save(file='/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/Chapter 2/Storms/Baseline/StormList_'+str(i)+'_baseline_RSD.npy',arr=base_storm_series)

for i in range(50,100):
    storm_series_shifted_10 = input_files.shift_storm_intensity(
        datadir=data_path,
        storm_list_name="VCRSeaStormMatrix_1980_2020.csv",  # can by .py or .csv
        mean_yearly_storms=42.25,
        SD_yearly_storms=4,
        shift=Mean_Rmax_10,  # shift the TWL distribution to the right (i.e., increase intensity), m NAVD88, typically [-0.15, 0.15]
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.03,  # m NAVD88, just used for plotting
        model_years=125,
        bPlot=False,
        bSave=False,
        output_filename="StormList_10kyrs_VCR_Berm1pt9m_Slope0pt04-lowIntensity",
    )
    np.save(file='/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/Chapter 2/Storms/Ten_Percent_Increase/StormList_'+str(i)+'_10_percent_increase_RSD.npy',arr=storm_series_shifted_10)

for i in range(50,100):
    storm_series_shifted_05 = input_files.shift_storm_intensity(
        datadir=data_path,
        storm_list_name="VCRSeaStormMatrix_1980_2020.csv",  # can by .py or .csv
        mean_yearly_storms=42.25,
        SD_yearly_storms=4,
        shift=Mean_Rmax_05,  # shift the TWL distribution to the right (i.e., increase intensity), m NAVD88, typically [-0.15, 0.15]
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.03,  # m NAVD88, just used for plotting
        model_years=125,
        bPlot=False,
        bSave=False,
        output_filename="StormList_10kyrs_VCR_Berm1pt9m_Slope0pt04-lowIntensity",
    )
    np.save(file='/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/Chapter 2/Storms/Five_Percent_Increase/StormList_'+str(i)+'_5_percent_increase_RSD.npy',arr=storm_series_shifted_05)