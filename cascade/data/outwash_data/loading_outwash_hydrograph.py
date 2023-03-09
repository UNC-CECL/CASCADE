# load the text file from Chris
# his elevations are in NAVD88

import numpy as np
from matplotlib import pyplot as plt

sound_levels = []

datadir = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/"
file_name = "sound_data.txt"

hydrograph = np.loadtxt(datadir+file_name, delimiter=",")  # this hydrograph is in m MSL
# using Beaufort, Duke Marine Lab NC datums to convert to dam MHW
sound_levels = np.asarray([(float(s) - 0.470)/10 for s in hydrograph])

fig, ax1 = plt.subplots()

ax1.set_xlabel('X-axis')
ax1.set_ylabel('m MSL', color='red')
ax1.plot(hydrograph, color='red')
ax1.tick_params(axis='y', labelcolor='red')
ax1.set_ylim(-1, 3)
# Adding Twin Axes
ax2 = ax1.twinx()
ax2.set_ylabel('m MHW', color='blue')
ax2.plot(sound_levels*10, color='blue', linestyle="dashed")
ax2.tick_params(axis='y', labelcolor='blue')
ax2.set_ylim(-1, 3)

sound_levels = sound_levels[21:45]
sound_levels[0] = 0


freq_increase = 10
# storms will be the sound data, so storm_series[1]
new_ss = []
num_new_vals = freq_increase - 1  # number of values we are adding between existing values
for s in range(len(sound_levels) - 1):
    new_ss.append(sound_levels[s])  # add the starting value back in
    inc = (sound_levels[s + 1] - sound_levels[s]) / freq_increase  # increment for the substeps
    for a in range(num_new_vals):
        new_ss.append(sound_levels[s] + (a + 1) * inc)
        # example, we want to do a substep of 3, so we are adding 2 new values
        # if our original storm series is [0, 0.15], then the increment will be 0.05
        # we will add 0 into the new storm series, then 0+1*0.05 = 0.05
        # then 0+2*0.05 = 0.1
        # the next iteration will start with 0.15, and that will be appended as well
new_ss.append(sound_levels[-1])  # make sure to include our last value
new_ss = np.asarray(new_ss)



plt.figure(3)
x = range(len(new_ss))
plt.scatter(x, new_ss*10, color="r")
plt.plot(new_ss*10)
plt.xlabel("hours")
plt.ylabel("m MHW")
plt.title("Bay Hydrograph")

num_years = 100
num_storms = 10
outwash_storms = np.empty([num_storms, 2], dtype=object)

interval = int(num_years/num_storms)

year1 = np.zeros([1,2], dtype=object)
year1[0, 0] = 1
year1[0, 1] = new_ss

for i in range(num_storms):
    outwash_storms[i, 0] = (i+1)*interval
    outwash_storms[i, 1] = new_ss

np.save(datadir+"outwash_storms_6min.npy", year1)