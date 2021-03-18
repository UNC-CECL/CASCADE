import CASCADE as CASCADE
import matplotlib.pyplot as plt
from scripts import CASCADE_plotters as CASCADE_Plt

name = "bulldozer"
ny = 3
datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"  # laptop
num_cores = ny
brie, barrier3d = CASCADE.initialize(
    datadir,
    name,
    ny=ny,
)

# advance time steps without a road
brie, barrier3d, _ = CASCADE.time_loop(
    brie,
    barrier3d,
    num_cores,
    nt=50,
    road_ele=None,
    road_width=None,
    road_setback=None,
)

# plot domain
plt.matshow(
    barrier3d[0].InteriorDomain * 10,
    origin="lower",
    cmap="terrain",
    vmin=-1.1,
    vmax=4.0,
)

# add a road
brie, barrier3d, road_setback = CASCADE.time_loop(
    brie,
    barrier3d,
    num_cores,
    nt=30,
    road_ele=1.7,  # 1.7 m NAVD88 for pt45, 1.5 m NAVD88 for pt75 (average of NC-12 is 1.3 m NAVD88), berm ele is 1.4 m NAV
    road_width=20,  # m
    road_setback=30,  # m
    artificial_max_dune_ele=3.7,  # m NAVD88, 4.7 = a 3 m dune above the roadway
    artificial_min_dune_ele=2.7,  # m NAVD88, 3.7 =  a 2 m dune above the roadway
)

directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
CASCADE_Plt.plot_ElevAnimation(barrier3d, ny, directory, TMAX=80, name="bulldozer_test")
