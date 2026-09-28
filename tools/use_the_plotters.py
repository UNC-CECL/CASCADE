# Lexi (Van Blunk) Fiegelist
# Sept. 24, 2026
# use the plotter functions to create figures for analysis

import os
import numpy as np
from cascade.tools.plotters import plot_ElevAnimation_CASCADE, plot_dune_domain, plot_start_end_domains

# choose whether to plot and save the cascade animations since it is automatic
# the other plots just produce and output on screen
plot_and_save_elev_anim = False
plot_dunes = False
plot_domains = True

# model results
datadir = r"C:\Users\agfig\model\calibration"
run_name = "calib_1"
file = os.path.join(datadir, "results", "{0}.npz".format(run_name))
calib = np.load(file, allow_pickle=True)["cascade"][0]
calib_b3d = calib.barrier3d  # all domains
calib_tmax = calib.barrier3d[0].TMAX  # since the whole model stops when one domain drowns,
# all domains should have the same TMAX

# domains
d_start_num = 3
d_end_num = 25
items = d_end_num - d_start_num + 1

# plot domains as a full island with offsets and save as gif
if plot_and_save_elev_anim:
    plot_ElevAnimation_CASCADE(
        cascade=calib,
        directory=r"C:\Users\agfig\model\calibration\results",
        TMAX_MGMT=0,
        name=run_name,
        TMAX_SIM=calib_tmax,
        ny=items,
        beach_management_ny=None,  # list of bool the length of ny, or None for all False
        roadway_management_ny=None,
        y_lim=(150, 350),
        z_lim=5,
        fig_size=(20, 5),
        fig_eps=False,
        km_on=True,
        )

if plot_dunes:
    # plot the dune domain over time
    plot_dune_domain(b3d=calib_b3d, TMAX=calib_tmax)

if plot_domains:
    # plot 1st and last year domains to compare features
    plot_start_end_domains(
            cascade_b3d=calib_b3d,
            save_dir=r"C:\Users\agfig\model\calibration\results\{0}\domain_comparison".format(run_name),
            min_z=-3,
            max_z=5,
    )