import CASCADE as CASCADE
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

# advance 5 time steps without a road
brie, barrier3d = CASCADE.time_loop(
    brie,
    barrier3d,
    num_cores,
    nt=50,
    road_ele=None,
    road_width=None,
    road_setback=None,
)

# plot domain
# plt.matshow(
#     barrier3d[0].InteriorDomain,
#     origin="lower",
#     cmap="terrain",
# )

# add a road
brie, barrier3d = CASCADE.time_loop(
    brie,
    barrier3d,
    num_cores,
    nt=30,
    road_ele=2.0,
    road_width=20,
    road_setback=30,
)

directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
CASCADE_Plt.plot_ElevAnimation(barrier3d, ny, directory, TMAX=80, name="bulldozer_test")
