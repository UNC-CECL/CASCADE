# Lexi (Van Blunk) Fiegelist
# Sept. 24, 2026
# use the plotter functions to create figures for analysis

import os
import numpy as np
from cascade.tools.plotters import (
    plot_ElevAnimation_CASCADE,
    plot_dune_domain,
    plot_start_end_domains,
    plot_overwash_flux,
    compare_modeled_observed_domains,
    plot_linear_dune_dif
    )

# choose whether to plot and save the cascade animations since it is automatic
# the other plots just produce and output on screen
plot_and_save_elev_anim = True
plot_dunes = True
plot_overwash = True
plot_domain_comparison = True
plot_linear_dune = True

# model results
datadir = r"C:\Users\agfig\model\calibration"
run_name = "calib_8"
file = os.path.join(datadir, "results", "{0}.npz".format(run_name))
calib = np.load(file, allow_pickle=True)["cascade"][0]
calib_b3d = calib.barrier3d  # all domains
calib_tmax = calib.barrier3d[0].TMAX  # since the whole model stops when one domain drowns,
# all domains should have the same TMAX

# 2014 datadir for comparison to observed results
dunes_2014_datadir = r"C:\Users\agfig\model\calibration\dunes_2014_final_berm2pt0"
elev_2014_datadir = r"C:\Users\agfig\model\calibration\domains_2014_final_berm2pt0"

# domains
d_start_num = 3
d_end_num = 25
items = d_end_num - d_start_num + 1
# make the list of dune and interior files
elev_files = []
dune_files = []
# load them in reverse order since that is how we modeled them and ib3d[0] corresponds to domain 25
for i in range(d_end_num, d_start_num-1, -1):
    dune_name = os.path.join(dunes_2014_datadir, 'domain_{0}_dunes_2014.npy'.format(i))  # decameters
    elev_name = os.path.join(elev_2014_datadir, 'domain_{0}_interior_2014.npy'.format(i))  # decameters above berm
    dune_files.append(dune_name)
    elev_files.append(elev_name)

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
        # y_lim=(0, 300),
        y_lim=(150, 350),
        z_lim=5,
        fig_size=(20, 5),
        fig_eps=False,
        km_on=True,
        flip_domain=True,
        invert_y=False
        )

if plot_dunes:
    # plot the dune domain over time
    dune_fig = plot_dune_domain(b3d=calib_b3d, TMAX=calib_tmax)
    # create the save folder if it does not already exist
    savedir = r"C:\Users\agfig\model\calibration\results\{0}".format(run_name)
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    dune_fig.savefig(os.path.join(savedir, "dunes.png"))


if plot_overwash:
    overwash = plot_overwash_flux(
        cascade_b3d=calib_b3d,
        save_dir=r"C:\Users\agfig\model\calibration\results\{0}".format(run_name),
        figsize=[15, 8]
    )

if plot_domain_comparison:
    compare_modeled_observed_domains(
        cascade_b3d=calib_b3d,
        save_dir=r"C:\Users\agfig\model\calibration\results\{0}\domain_comparison".format(run_name),
        observed_domains=elev_files,
        observed_dunes_list=dune_files,
        min_z=-3,
        max_z=5,
        figsize=[15, 8]
    )

if plot_linear_dune:
    plot_linear_dune_dif(
        cascade_b3d=calib_b3d,
        save_dir=r"C:\Users\agfig\model\calibration\results\{0}\dune_differences".format(run_name),
        observed_dunes_list=dune_files,
        min_z=-5,
        max_z=5,
        figsize=[15, 8]
)