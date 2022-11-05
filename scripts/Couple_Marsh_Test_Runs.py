# Nov 2 2022
# Test methods for running a version of Cascade with marsh dynamics from BarrierBMFT
# Script 1 Run 1 singular instance of Cascade

import numpy as np
import time

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import pandas as pd
import os
import imageio
import math
from scipy import signal

from scripts import CASCADE_plotters as CASCADEplt
from cascade.cascade import Cascade
from barrier3d.tools.input_files import (
    yearly_storms,
    gen_dune_height_start,
    gen_alongshore_variable_rmin_rmax,
)


# ###############################################################################
# 4 - CASCADE with only one B3D model and no human dynamics
# ###############################################################################
# Use the starting interior domain from the 10,000 yr runs for each dune growth rate and run for 1000 years
# or until the barrier drowns. All other modules (brie and human dymnamics modules) turned off. Can also use this
# run script for the 10,000 year runs.

def Single_Marsh_Plot(
    nt,
    rmin,
    rmax,
    name,
    storm_file,
    elevation_file,
    dune_file,
    marsh_dynamics,
):
    # ###############################################################################
    # 4 - CASCADE with only one B3D model and no human dynamics
    # ###############################################################################
    # Use the starting interior domain from the 10,000 yr runs for each dune growth rate and run for 1000 years
    # or until the barrier drowns. All other modules (brie and human dymnamics modules) turned off. Can also use this
    # run script for the 10,000 year runs.

    # --------- INITIALIZE ---------
    datadir = "B3D_Inputs/"
    cascade = Cascade(
        datadir,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file="RUN4-CASCADE-parameters.yaml",
        wave_height=1,
        wave_period=7,  # s (lowered from 10 s to reduce k_sf)
        wave_asymmetry=0.8,  # fraction approaching from left
        wave_angle_high_fraction=0.2,  # fraction of waves approaching from higher than 45 degrees
        sea_level_rise_rate=0.004,  # m/yr
        sea_level_rise_constant=True,  # linear SLR
        background_erosion=0.0,
        alongshore_section_count=1,  # only one B3D domain
        time_step_count=nt,
        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,
        num_cores=1,
        roadway_management_module=False,  # no roadway management
        alongshore_transport_module=False,  # no brie coupling
        beach_nourishment_module=False,  # no beach nourishment
        community_dynamics_module=False,  # no community dynamics
        marsh_dynamics = marsh_dynamics, #no marsh dynamics
    )

    # --------- LOOP ---------
    Time = time.time()

    for time_step in range(nt - 1):
        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        cascade.update()
        if cascade.b3d_break:
            break

    # --------- SAVE ---------
    save_directory = "Run_Output/"
    cascade.save(save_directory)

    return cascade


# ###############################################################################
# Example of Running Function 1
# ###############################################################################

# Specify variables to use in calling function
# Elevation file path name
e_file = "/Users/ceclmac/PycharmProjects/TestPhython/CP4.npy"
# Dune height path name
d_file = "/Users/ceclmac/PycharmProjects/CASCADE/B3D_Inputs/DuneStart_1000dam.npy"
# Storm file path name
s_file = "/Users/ceclmac/PycharmProjects/CASCADE/B3D_Inputs/StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy"

# Call function
Single_Marsh_Plot(
    nt=5,
    rmin=0.25,
    rmax=0.65,
    name="Marsh_Script_Tests",
    storm_file=s_file,
    elevation_file=e_file,
    dune_file=d_file,
    marsh_dynamics=True,
)