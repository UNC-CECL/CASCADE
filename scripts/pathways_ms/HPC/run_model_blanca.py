# Model runfile for CASCADE simulations on the HPC
# Written by K.Anarde

# remember if I move to a different computer to $ pip install -e . in the brie, CHOM, and B3D directories

import numpy as np
import time

from cascade.cascade import Cascade  # the new class


def RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
    nt,
    rmin,
    rmax,
    name,
    storm_file,
    elevation_file,
    dune_file,
):

    # ###############################################################################
    # 4 - CASCADE with only one B3D model and no human dynamics
    # ###############################################################################
    # Use the starting interior domain from the 10,000 yr runs for each dune growth rate and run for 1000 years
    # or until the barrier drowns. All other modules (brie and human dymnamics modules) turned off. Can also use this
    # run script for the 10,000 year runs.

    # --------- INITIALIZE ---------
    datadir = "cascade/data/pathways_data/"
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


cascade_10kyr_pt45_01 = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
    nt=10,
    rmin=0.25,  # rave = 0.45 (but not 0.5 spaced like in Reeves et al., 2021 -- arbitrary)
    rmax=0.65,
    name="4-CASCADE_noAST_Rave_pt45_SLR_pt004_10k-yrs_01",
    storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
    elevation_file="InitElevHog.npy",
    dune_file="DuneStart_1000dam.npy",
)
