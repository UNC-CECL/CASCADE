import numpy as np
import pandas as pd
import os
from matplotlib import pyplot as plt

os.chdir('C:/Users/Lexi/Documents/Research/ArcPro/')
# points_df = pd.read_csv("raster_points_table.csv")  # same as RawMatrix
# points_df = pd.read_csv("big_raster_points.csv")
# points_df = pd.read_csv("big_poly2_points.csv")
points_df = pd.read_csv("biggest_poly_points.csv")
elevations = points_df['grid_code']  # same as RawMatrix[,3]

num_cols_domain = 150
num_rows_domain = int(len(elevations)/num_cols_domain)  # this will be the length of the island
MHW = 0.46

elevations = elevations - MHW
# set all values less than 0 to -3 m
elevations[elevations < 0] = -3

island_width = 30
island_length = num_rows_domain

processed_elev_array = np.zeros([island_width, island_length])
dune_elev_array = np.zeros(island_length)

start = 0
end = start + num_cols_domain
for l in range(num_rows_domain):
    # elev_array[row, :] = elevations[start:end]
    z = elevations[start:end]  # z in Bentons code (one cross-shore row)
    revz = np.flip(z)
    dune_max_location = np.argmax(revz)
    dune_elev = np.max(revz)
    # Start Island Immediatly after duneline
    start_island = (dune_max_location + 1)
    end_island = (dune_max_location + island_width+1)

    one_elev_col = revz[start_island:end_island]
    # Add values to larger matrix
    processed_elev_array[:, island_length-l-1] = one_elev_col
    dune_elev_array[island_length-l-1] = dune_elev

    start += num_cols_domain
    end += num_cols_domain

smaller_length_array = processed_elev_array[:, 0:49]

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
mat = ax1.matshow(
    smaller_length_array,
    cmap="terrain",
    vmin=-3.0
    # vmax=3.0,
)
cbar = fig1.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=15)
ax1.set_title("Initial Elevation")
ax1.set_ylabel("barrier width (dam)")
ax1.set_xlabel("barrier length (dam)")
plt.gca().xaxis.tick_bottom()