# Lexi (Van Blunk) Fiegelist
# 09/24/2026

# code for running CASCADE for overwash and marsh dynamics
# berm elevation, slopes, and elevations are currently set for Masonboro Island, NC
# other variables set for calibration (2004-2014)


import time
import os
import numpy as np

from cascade.cascade import Cascade
from cascade.tools.plotters import plot_ElevAnimation_CASCADE

# input datadir where the 100 storms are located
datadir = r"C:\Users\agfig\model\calibration"

# make the list of dune and interior files
elev_files = []
dune_files = []
# folder that contains the dunes and interior elevations for all domains
dunes_datadir = r"C:\Users\agfig\model\calibration\dunes_2004_final_shrubs_v1_berm2pt0"
elev_datadir = r"C:\Users\agfig\model\calibration\domains_2004_final_shrubs_v1_berm2pt0"

# not including end domains in the model (1, 2, 26)
d_start_num = 3
d_end_num = 25
items = d_end_num - d_start_num + 1
for i in range(d_start_num, d_end_num+1):
    dune_name = os.path.join(dunes_datadir, 'domain_{0}_dunes_2004.npy'.format(i))
    elev_name = os.path.join(elev_datadir, 'domain_{0}_interior_2004.npy'.format(i))
    dune_files.append(dune_name)
    elev_files.append(elev_name)

# ---------------------------------- set model parameters that change per run ------------------------------------------
beach_slope = 0.06
model_duration = 10
berm_elev = 2.00  # m NAVD88
MHW = 0.421

# change r values per domain
# NOTE: has to be a list and not an array because numpy array uses type float64 instead of just a float, which raises
# an error in the load_inputs/set_yaml functions because rmin and rmax are specified as float types in configuration.py
min_dune_r = [0.25,
0.25,
0.25,
0.25,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.10,
0.10,
0.10,
0.10,
0.10]
max_dune_r = [0.45,
0.45,
0.45,
0.45,
0.35,
0.35,
0.35,
0.35,
0.35,
0.35,
0.35,
0.35,
0.35,
0.35,
0.35,
0.35,
0.35,
0.35,
0.40,
0.40,
0.40,
0.40,
0.40]


# save to results folder
save_dir = r"C:\Users\agfig\model\calibration\results"
# create the folder if it does not already exist
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# offsets
# Units: CSV is in meters, convert to decameters for CASCADE
dune_offsets_file = r"C:\Users\agfig\model\calibration\offsets\Island_Dune_Offsets_2004_CASCADE_Input.csv"
dune_offsets_raw = np.loadtxt(dune_offsets_file, skiprows=1, delimiter=',')
dune_offsets = dune_offsets_raw / 10  # Convert m → dam
print(f"✓ Loaded dune offsets: {len(dune_offsets)} domains")
print(f"  Range: {np.min(dune_offsets):.1f} to {np.max(dune_offsets):.1f} dam")

# --------------------------------- running overwash scenario for 1 storm --------------------------------------
overwash_storm = "masonboro-storms1.npy"
run_name = "calib_1_flip_domain_order"

# reverse the domain list and everything related to the domains
elev_files.reverse()
dune_files.reverse()
min_dune_r.reverse()
max_dune_r.reverse()
dune_offsets = np.flip(dune_offsets)

# initialize class
cascade_marsh = Cascade(
    datadir,
    name=run_name,
    elevation_file=elev_files,
    dune_file=dune_files,
    parameter_file="marsh-default-parameters.yaml",
    storm_file=overwash_storm,
    num_cores=1,  # cascade can run in parallel, can never specify more cores than that
    roadway_management_module=False,
    alongshore_transport_module=False,
    beach_nourishment_module=False,
    community_economics_module=False,
    outwash_module=False,
    marsh_module=False,
    alongshore_section_count=items,
    time_step_count=model_duration,
    wave_height=1,  # ---------- for BRIE and Barrier3D --------------- #
    wave_period=7,
    wave_asymmetry=0.8,
    wave_angle_high_fraction=0.2,
    bay_depth=3.0,
    s_background=0.001,
    berm_elevation=berm_elev,
    MHW=MHW,
    beta=beach_slope,
    sea_level_rise_rate=0.004,
    sea_level_rise_constant=True,
    background_erosion=0.0,
    min_dune_growth_rate=min_dune_r,
    max_dune_growth_rate=max_dune_r,
    # --------------- marsh dynamics -----------------------------------
    SSCb=0.05,
    organic_content_bay=0,
    numiterations=500,
    tidal_period=12.5,
    settling_velocity=0.05 * 10 ** (-3),
    max_biomass=2500,
    min_depth_marsh_growth=0,
    max_depth_marsh_growth=0.4,
    density_organic_matter=85,
    density_sediment=2000,
    max_depth_decomp=0.4,
    decomp_coeff=0.1,
    min_elev_marsh=-0.2,
    max_elev_marsh=0,
    tidal_amplitude=0.7,
    accretion_method=4,
    # ---------- offsets ----------- #
    enable_shoreline_offset=True,
    shoreline_offset=dune_offsets,

)

# run the time loop/update function
t0 = time.time()

for time_step in range(cascade_marsh._nt - 1):
    print("\r", "Time Step: ", time_step + 1, end="")
    cascade_marsh.update()
    if cascade_marsh.b3d_break:
        break

t1 = time.time()
t_total_seconds = t1 - t0
t_total_minutes = t_total_seconds / 60
t_total_hours = t_total_seconds / 3600

# save variables
cascade_marsh.save(save_dir)


# plot domains
plot_ElevAnimation_CASCADE(
    cascade=cascade_marsh,
    directory=r"C:\Users\agfig\model\calibration\results",
    TMAX_MGMT=0,
    name=run_name,
    TMAX_SIM=model_duration,
    ny=items,
    beach_management_ny=None,  # list of bool the length of ny, or None for all False
    roadway_management_ny=None,
    y_lim=(150, 350),
    z_lim=5,
    fig_size=(20, 5),
    fig_eps=False,
    km_on=True,
    )