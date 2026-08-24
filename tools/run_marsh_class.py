# Lexi Van Blunk
# 11/22/2025

# code for running CASCADE with marsh dynamics
import copy
import time
from cascade.marsh_dynamics import Marsh
import numpy as np
from matplotlib import pyplot as plt

# load a test domain and transect which is just an arbitrary segment of Masonboro Island
# test_domain = np.load(r"C:\Users\agfig\model\basic_cascade_run\marsh\marsh_domain_test3.npy")  # dam MHW
test_domain = np.load(r"C:\Users\agfig\model\marsh_domain_test.npy") - 0.07
n_cols = np.shape(test_domain)[1]
width = np.shape(test_domain)[0]
initial_domain = copy.deepcopy(test_domain)

# set model duration
model_duration = 1  # yrs
shoreline_change = np.zeros(model_duration)

# initialize the class
marsh_class = Marsh(
    RSLR=0.004,  # this will get replaced with cascade: self._sea_level_rise_rate which is in m/yr
    C_e=0.05,
    OCb=0,
    numiterations=500,
    P=12.5 * 3600 * 1,
    ws=0.05 * 10 ** (-3),
    n_tidal_cycles=365 * (24 / 12.5),
    Bmax=2500,
    Dmin=0,
    Dmax=1,
    rhoo=85,
    rhos=2000,
    mui=0.4,
    mki=0.1,
    m_min=-0.3,
    m_max=0,
    time_step_count=model_duration,
    alongshore_length=n_cols,
    tidal_amplitude=0.7,
    initial_width=width
    )


# run the time loop/update function
t0 = time.time()

for time_step in range(marsh_class._nt):
    print("\r", "Time Step: ", time_step + 1, end="")
    test_domain = marsh_class.update(
        interior_domain=test_domain,  # this will be the b3d interior domain at the time step
        model_year=time_step,
        shoreline_changeTS=shoreline_change
        )
    # test_domain = marsh_class._marsh_elevation[time_step]

t1 = time.time()
t_total_seconds = t1 - t0
t_total_minutes = t_total_seconds / 60
t_total_hours = t_total_seconds / 3600

# save domains
np.save(r"C:\Users\agfig\model\basic_cascade_run\marsh\results\work_comp_marsh2.npy",marsh_class._marsh_elevation)

# # save variables
# save_dir = r"C:/Users/Lexi\Documents\UNC\BarrierBMFT\tests\post_integration"
# marsh_class.save(save_dir)

plot_on = False
save_fig = False
# plot the results
if plot_on:
    # initial domain
    xlabel = "alongshore distance (m)"
    ylabel = "cross-shore distance (m)"
    minz = -3
    maxz = 5

    plot_years = np.linspace(0, len(marsh_class._marsh_elevation), 3, dtype=int)
    plot_years[-1] = plot_years[-1] - 1
    n = len(plot_years)
    plot_n = 1
    plot_width = (n-1)*10

    fig1 = plt.figure(figsize=[plot_width, 10])
    for p in plot_years:
        ax1 = fig1.add_subplot(1,n,plot_n)
        mat1 = ax1.matshow(
            marsh_class._marsh_elevation[p] * 10,
            cmap="terrain",
            vmin=minz,
            vmax=maxz,
            origin="lower"
        )
        cbar = fig1.colorbar(mat1)
        cbar.set_label('elevation (m MHW)', rotation=270, labelpad=15)
        ax1.set_xlabel(xlabel)
        ax1.set_title("t = {0}".format(p))
        fig1.gca().xaxis.tick_bottom()
        plt.gca().invert_yaxis()
        plot_n += 1
    fig1.tight_layout()

    if save_fig:
        fig1.savefig(r"C:\Users\Lexi\Documents\UNC\model\output_figures\marsh_t{0}_msl_conversion".format(marsh_class._nt))