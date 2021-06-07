from CASCADE import Cascade
import matplotlib.pyplot as plt
from scripts import CASCADE_plotters as CASCADE_Plt

name = "bulldozer"
ny = 3
datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"  # laptop
num_cores = ny
nt = 100

# start a model without a road, and then add a road
cascade = Cascade(
    datadir,
    name=name,
    alongshore_section_count=ny,
    time_step_count=nt,
    num_cores=ny,
    roadway_management_module=False,
    alongshore_transport_module=True,
    road_ele=1.7,
    road_width=30,
    road_setback=30,
    artificial_max_dune_ele=3.7,
    artificial_min_dune_ele=2.2,
)

# advance the first 50 time steps without a road
nt = 50
for time_step in range(nt - 1):
    # Print time step to screen (NOTE: time_index in each model is time_step+1)
    print("\r", "Time Step: ", time_step, end="")
    cascade.update()

# plot domain
plt.matshow(
    cascade.barrier3d[0].InteriorDomain * 10,
    origin="lower",
    cmap="terrain",
    vmin=-1.1,
    vmax=4.0,
)
plt.colorbar()

# add a road and run for another 30 time steps
cascade.roadway_management_module = True  # turn roadway management model on
nt = 30
for time_step in range(nt - 1):
    # Print time step to screen (NOTE: time_index in each model is time_step+1)
    print("\r", "Time Step: ", time_step, end="")
    cascade.update()

directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
CASCADE_Plt.plot_ElevAnimation(
    cascade.barrier3d, ny, directory, TMAX=80, name="bulldozer_test"
)
