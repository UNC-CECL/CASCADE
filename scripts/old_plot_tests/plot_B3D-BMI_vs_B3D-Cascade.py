"""
    Plot check that barrier3d call in cascade matches bmi

    NOTE to users:
        - if using Barrier3D for the first time, remember to $ pip install -e .
"""

import time

import matplotlib.pyplot as plt
from barrier3d import Barrier3dBmi
from barrier3d.tools import plot as B3Dfunc

from cascade.cascade import Cascade

# specify data directories with initial conditions
datadir = "cascade/data/pathways_data/"

nt = 100
cascade = Cascade(
    datadir,
    name="test",
    alongshore_section_count=1,
    time_step_count=nt,
    num_cores=1,
    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
    elevation_file="barrier3d-default-elevation.npy",
    dune_file="barrier3d-default-dunes.npy",
    parameter_file="barrier3d-default-parameters.yaml",
    wave_height=1,
    wave_period=7,
    wave_asymmetry=0.8,
    wave_angle_high_fraction=0.2,
    sea_level_rise_rate=0.004,
    min_dune_growth_rate=0.35,
    max_dune_growth_rate=0.85,
    roadway_management_module=False,
    alongshore_transport_module=False,
    beach_nourishment_module=False,
)

for time_step in range(nt - 1):
    # Print time step to screen (NOTE: time_index in each model is time_step+1)
    print("\r", "Time Step: ", time_step, end="")
    cascade.update()

# Barrier3D Version 2.0 ------------------------------
# create an instance of the new BMI class, which is the model
barrier3d = Barrier3dBmi()
input_file = "barrier3d-default-parameters.yaml"
barrier3d.initialize(datadir + input_file)

# increase time step
Time = time.time()
for time_step in range(1, barrier3d._model._TMAX):
    barrier3d.update()
    print("\r", "Time Step: ", time_step, end="")
SimDuration = time.time() - Time
print()
print("Elapsed Time: ", SimDuration, "sec")  # Print elapsed time of simulation

# ------------------------------

# dune Height Over Time (input in decameter)
B3Dfunc.plot_dune_height(cascade.barrier3d[0]._DuneDomain, cascade.barrier3d[0]._Dmax)
B3Dfunc.plot_dune_height(barrier3d._model._DuneDomain, barrier3d._model._Dmax)

# check that shoreline change is the same between the two and that shoreface slope starts in equilibrium
plt.figure()
plt.plot(cascade.barrier3d[0].x_s_TS, "b")
plt.plot(barrier3d._model.x_s_TS, "g")
