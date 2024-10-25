from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from barrier3d import Barrier3d, Barrier3dBmi
from barrier3d.tools import input_files
import os
import sys
cwd = os.getcwd()
data_path = '/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/Chapter 3/Historic Storm Values/'

storm_list = np.loadtxt(data_path+"/SeaStormMatrix_OCR_12.csv",delimiter=',')

mean_storms = 2.45
sd_storms = 1.7
Ocracoke_MHW = 0.26
Berm_Elevation = 1.7

sys.path
for i in range(1,100):
    base_storm_series = input_files.shift_storm_intensity(
        datadir=data_path,
        storm_list_name="SeaStormMatrix_OCR_12.csv",  # can by .py or .csv
        mean_yearly_storms=mean_storms,
        SD_yearly_storms=sd_storms,
        shift=0,  # shift the TWL distribution to the right (i.e., increase intensity), m NAVD88, typically [-0.15, 0.15]
        MHW=Ocracoke_MHW,  # m NAVD88
        StormStart=0,
        BermEl=Berm_Elevation,  # m NAVD88, just used for plotting
        model_years=150,
        bPlot=False,
        bSave=False,
        output_filename=("OCR_Baseline_Storm_"+str(i)),
    )
    np.save(file='/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/Chapter 3/Storm Output/Npy files/Ocracoke_Future_StormList_'+str(i)+'_baseline.npy',arr=base_storm_series)

'''for i in range(0,1):
    storm_series_shifted_05 = input_files.shift_storm_intensity(
        datadir=data_path,
        storm_list_name="StormList_Duck.csv",  # can by .py or .csv
        mean_yearly_storms=mean_storms,
        SD_yearly_storms=sd_storms,
        shift=0.5,  # (used 2.5) shift the TWL distribution to the right (i.e., increase intensity), m NAVD88, typically [-0.15, 0.15]
        MHW=Ocracoke_MHW,  # m NAVD88
        StormStart=2,
        BermEl=Berm_Elevation,  # m NAVD88, just used for plotting
        model_years=126,
        bPlot=False,
        bSave=False,
        output_filename="StormList_10kyrs_VCR_Berm1pt9m_Slope0pt04-lowIntensity",
    )
    np.save(file='/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/Chapter 3/Storm Output/Npy files/Ocracoke_Bigger_Storms_5_m.npy',arr=storm_series_shifted_05)
'''