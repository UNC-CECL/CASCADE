from pathlib import Path

import matplotlib.pyplot as plt
from tqdm import tqdm

from barrier3d import Barrier3d, Barrier3dBmi
from barrier3d.tools import input_files
from barrier3d.tools.plot import plot_dune_height

import os
import sys

number_storms = 289
years = 41
mean_storms = 6.7
sd_storms = 3.2
Ocracoke_MHW = 0.26
Berm_Elevation = 1.7

data_path = 'C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\Ocracoke_init_data'

os.getcwd()

sys.path

storm_series_normal = input_files.yearly_storms(
    datadir=data_path,
    storm_list_name="StormList_20k_Ocracoke.csv",  # can by .py or .csv
    mean_yearly_storms=mean_storms,
    SD_yearly_storms=sd_storms,
    MHW=Ocracoke_MHW,  # m NAVD88
    StormStart=2,
    BermEl=Berm_Elevation,  # m NAVD88, just used for plotting
    model_years=1000,
    bPlot=True,
    bSave=False,
    output_filename="StormList_1kyrs_VCR_Berm1pt9m_Slope0pt04",
)