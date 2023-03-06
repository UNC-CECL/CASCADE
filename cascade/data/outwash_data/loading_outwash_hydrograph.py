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
plt.figure(3)
plt.plot(sound_levels)


num_years = 100
num_storms = 10
outwash_storms = np.empty([num_storms, 2], dtype=object)

interval = int(num_years/num_storms)

for i in range(num_storms):
    outwash_storms[i, 0] = (i+1)*interval
    outwash_storms[i, 1] = sound_levels

# np.save(datadir+"outwash_storms.npy", outwash_storms)