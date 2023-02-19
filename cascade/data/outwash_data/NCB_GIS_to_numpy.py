import numpy as np
import pandas as pd
import os
from matplotlib import pyplot as plt

os.chdir('C:/Users/Lexi/Documents/Research/ArcPro/')
df = pd.read_csv("back_barrier_points.csv")
# df = pd.read_csv("big_raster_points.csv")
# df = pd.read_csv("big_poly2_points.csv")
# df = pd.read_csv("biggest_poly_points.csv")
elevations = df['grid_code']
elevs_numpy = np.zeros(len(elevations))

num_cols_domain = 300
num_rows_domain = int(len(elevations)/num_cols_domain)

for index, val in enumerate(elevations):
    elevs_numpy[index] = val
    if val < -0.1:
        elevs_numpy[index] = -3

elev_array = np.zeros([num_rows_domain, num_cols_domain])
start = 0
end = start + num_cols_domain
for row in range(num_rows_domain):
    elev_array[row, :] = elevs_numpy[start:end] - 0.46
    start += num_cols_domain
    end += num_cols_domain

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
mat = ax1.matshow(
    elev_array,
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


width = 40
length = 50
vertical_array = np.zeros([width, length])
start_row = 59
start_col = 0
# for col in range(length):
#     if col != 0 and col % 2 == 0:
#         start_row -= 1
#     elif col != 0 and col % 2 != 0:
#         start_col += 1
#     for i in range(width):
#         vertical_array[i, col] = elev_array[start_row + i, start_col + i]
for col in range(length):
    if col != 0:
        start_row -= 1
        start_col += 1
    for i in range(width):
        vertical_array[i, col] = elev_array[start_row + i, start_col + i]


width = 40
length = 55
vertical_array = np.zeros([width, length])
for w in range(width):
    start_col = 67
    start_row = 45 - w
    for l in range(length):
        vertical_array[width-1-w, length-1-l] = elev_array[start_row, start_col]
        start_row = start_row + 1
        start_col = start_col - 1

# width = 40
# length = 50
# vertical_array = np.zeros([width, length])
# for l in range(length):
#     start_col = 10 + l
#     start_row = 40
#     for w in range(width):
#         vertical_array[w, l] = elev_array[start_row, start_col]
#         start_row = start_row + 1
#         start_col = start_col + 1

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
mat2 = ax2.matshow(
    vertical_array,
    cmap="terrain",
    # vmin=-3.0, vmax=3.0,
)
ax2.set_xlabel('barrier length (dam)')
ax2.set_ylabel('barrier width (dam)')
ax2.set_title("Initial Elevation")
plt.gca().xaxis.tick_bottom()
cbar = fig2.colorbar(mat2)
cbar.set_label('m MHW', rotation=270, labelpad=15)