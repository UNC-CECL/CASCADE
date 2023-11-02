from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from barrier3d import Barrier3d, Barrier3dBmi
from barrier3d.tools import input_files
import os
import sys
cwd = os.getcwd()
data_path = cwd+'/data/marsh_init_data'

StormList = np.loadtxt(datadir / storm_list_name, delimiter=",")
storm_list = np.loadtxt(data_path+"/StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",delimiter=',')

RMAX = storm_list[:,2]
Mean_Rmax = np.mean(RMAX)
Mean_Rmax_10 = Mean_Rmax*.1
Mean_Rmax_02 = Mean_Rmax*.02


sys.path

storm_series_shifted = input_files.shift_storm_intensity(
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

Mean_base = np.mean(storm_series_shifted[:,1])
Mean_10 = np.mean(storm_series_shifted_10[:,1])
Mean_02 = np.mean(storm_series_shifted_02[:,1])

print((Mean_10-Mean_base)/Mean_base)
print((Mean_02-Mean_base)/Mean_base)

