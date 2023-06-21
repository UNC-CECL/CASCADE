# load the text file from Chris
# his elevations are in m MSL

import numpy as np
from matplotlib import pyplot as plt

datadir = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/"
file_name = "sound_data.txt"  # the ADCIRC generated hydrograph from USGS
hydro_mMSL = np.loadtxt(datadir+file_name, delimiter=",")
# use Beaufort, Duke Marine Lab NC datums to convert to dam MHW: MHW is 0.47 m MSL, dam = m / 10
hydro_damMHW = np.asarray([(float(s) - 0.470)/10 for s in hydro_mMSL])

plt.rcParams.update({"font.size": 15})
# plotting the original hydrograph in m MSL
# points are every hour
fig, ax1 = plt.subplots()
x_vals = np.asarray(range(len(hydro_mMSL)))/24
ax1.set_xlabel('days')
ax1.set_title("MSL and MHW ADCIRC hydrographs")
ax1.set_ylabel('m MSL', color='g')
# ax1.plot(x_vals, hydro_mMSL, color='red')
ax1.scatter(x_vals, hydro_mMSL, color='g')
ax1.tick_params(axis='y', labelcolor='g')
ax1.set_ylim(-1, 3)
# Adding Twin Axes for m MHW
ax2 = ax1.twinx()
ax2.set_ylabel('m MHW', color='m')
# ax2.plot(x_vals, hydro_mMHW*10, color='blue', linestyle="dashed")
ax2.scatter(x_vals, hydro_damMHW*10, color='m')
ax2.tick_params(axis='y', labelcolor='m')
ax2.set_ylim(-1, 3)

# use only the highest 12 hours of data
hydro_damMHW = hydro_damMHW[21:34]
hydro_damMHW[0] = 0

plt.figure(4)
x = range(len(hydro_damMHW))
plt.plot(hydro_damMHW*10, color="darkgrey")
plt.scatter(x, hydro_damMHW*10, color="r")
plt.xlabel("hours")
plt.ylabel("m MHW")
plt.title("1 Hour Hydrograph")

freq_increase = 10  # number of points within an hour, must be a factor of 60. examples below
# 1 = 60 minute intervals (original hydrograph)
# 2 = 30 minute intervals
# 3 = 20 minute intervals
# 4 = 15 minute intervals
# 5 = 12 minute intervals
# 6 = 10 minute intervals
# 10 = 6 minute intervals
# etc.

# storms will be the sound data, so storm_series[1]
hydro_damMHW_freq_increase = []
num_new_vals = freq_increase - 1  # number of values we are adding between hourly data points, have to subtract 1
# because we are incorporating the hourly values
for s in range(len(hydro_damMHW) - 1):
    hydro_damMHW_freq_increase.append(hydro_damMHW[s])  # add the starting (hourly) value back in
    inc = (hydro_damMHW[s + 1] - hydro_damMHW[s]) / freq_increase  # increment for the substeps
    for a in range(num_new_vals):
        hydro_damMHW_freq_increase.append(hydro_damMHW[s] + (a + 1) * inc)
        # example, we want to do a substep of 3, so we are adding 2 new values
        # if our original storm series is [0, 0.15], then the increment will be 0.05
        # we will add 0 into the new storm series, then 0+1*0.05 = 0.05
        # then 0+2*0.05 = 0.1
        # the next iteration will start with 0.15, and that will be appended as well
hydro_damMHW_freq_increase.append(hydro_damMHW[-1])  # make sure to include our last value
hydro_damMHW_freq_increase = np.asarray(hydro_damMHW_freq_increase)

save_hydro = False
if save_hydro:
    name_hydro = "outwash_storms_6min.npy"
    np.save(datadir + name_hydro, hydro_damMHW_freq_increase)

# all for plotting the incrememnts on a 12-hr x axis (I think)
# new_ss_list = []
# for s in range(len(new_ss)):
#     if s == 0:
#         new_ss_list.append(new_ss[s])
#     elif s % 8 == 0:
#         new_ss_list.append(new_ss[s])
#
# new_ss = np.asarray(new_ss_list)

# plt.figure(3)
# x_3 = np.asarray(range(len(new_ss)))*8
# x_1 = []
#
# for t in x_3:
#     if t == 0:
#         x_1.append(t)
#     elif t % 60 == 0:
#         x_1.append(t)
#
# plt.plot(x_3, new_ss*10, color='darkgrey')
# plt.scatter(x_3, new_ss*10, color="b")
# plt.scatter(x_1, sound_levels*10, color="r")
# plt.xlabel("hours")
# plt.ylabel("m MHW")
# label_locs = np.arange(0, max(x_3)+1, step=120, dtype=int)
# labels = label_locs // 60
# plt.xticks(label_locs, labels)
#
# plt.title("8 Minute Hydrograph")

# plt.figure(3)
# x_15 = np.asarray(range(len(new_ss)))*15
# x_1 = []
#
# for t in x_15:
#     if t == 0:
#         x_1.append(t)
#     elif t % 60 == 0:
#         x_1.append(t)
#
# plt.plot(x_15, new_ss*10, color='darkgrey')
# plt.scatter(x_15, new_ss*10, color="b")
# plt.scatter(x_1, sound_levels*10, color="r")
# plt.xlabel("hours")
# plt.ylabel("m MHW")
# label_locs = np.arange(0, max(x_15)+1, step=120, dtype=int)
# labels = label_locs // 60
# plt.xticks(label_locs, labels)
#
# plt.title("15 Minute Hydrograph")

# creating the full storm series (years in column 1 and hydrograph in column 2)
num_simulation_years = 100
num_storms = 5
outwash_storms = np.empty([num_storms, 2], dtype=object)  # has to be object bc we are adding an array to the second col
interval = int(num_simulation_years/num_storms)  # the years the storms occur

for i in range(num_storms):
    outwash_storms[i, 0] = (i+1)*interval
    outwash_storms[i, 1] = hydro_damMHW_freq_increase

# add an additional outwash storm immediately (not included in n_storms)
storm_at_year_1 = True
if storm_at_year_1:
    year1 = np.zeros([1,2], dtype=object)
    year1[0, 0] = 1
    year1[0, 1] = hydro_damMHW_freq_increase
    outwash_storms = np.vstack((year1, outwash_storms))

save_storm_series = False
if save_storm_series:
    name = "outwash_storms_20yrs_early"
    np.save(datadir + name, outwash_storms)

