# Model plot file for CASCADE simulations
# Written by K.Anarde

import numpy as np
import matplotlib.pyplot as plt
import os

#===================================================
# CASCADE vs brie (LTA)



# plot_BRIE_frames
# CASCADEplt.plot_BRIE_frames(brieLTA._y,
#                             brieLTA._x_t_save,
#                             brieLTA._x_s_save,
#                             brieLTA._x_b_save,
#                             brieLTA._nt,
#                             brieLTA._ny,
#                             brieLTA._dy,
#                             brieLTA._dtsave,
#                             brieLTA._inlet_idx,
#                             'True',
#                             'test')

y = brieLTA._y
x_t = brieLTA._x_t_save
x_s = brieLTA._x_s_save
x_b = brieLTA._x_b_save
nt = brieLTA._nt
ny = brieLTA._ny
dy = brieLTA._dy
dtsave = brieLTA._dtsave
inlet_idx = brieLTA._inlet_idx
make_gif = 'True'
file_name = 'test'
directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"

# create time array
t = np.arange(0, nt, dtsave)

if make_gif:
    os.chdir(directory)
    os.mkdir(directory+'/GIF')

fig, axes = plt.subplots()
frames = []
ymax = int(math.ceil(np.max(x_b[:, -1]) / 100.0)) * 100
ymin = int(math.floor(np.min(x_t[:, 0]) / 100.0)) * 100

for i in np.arange(0, np.size(t)):

    # plot initial conditions
    axes.plot(y / 1000, x_t[:, i], color="cornflowerblue")
    axes.plot(y / 1000, x_s[:, i], color="gold")
    axes.plot(y / 1000, x_b[:, i], color="teal")
    axes.fill_between(y / 1000, ymin, x_t[:, i], color="royalblue")
    axes.fill_between(y / 1000, x_b[:, i], ymax, color="teal", alpha=0.6)
    axes.fill_between(y / 1000, x_t[:, i], x_s[:, i], color="cornflowerblue", alpha=0.6)
    axes.fill_between(y / 1000, x_s[:, i], x_b[:, i], color="gold", alpha=0.6)
    axes.legend(['x_t', 'x_s', 'x_b'])
    plt.xlabel('Alongshore (km)')
    plt.ylabel('Cross-shore (m)')
    plt.title('Time = ' + str(int(t[i])) + ' yr')
    axes.set_xlim(0, (ny-1) * dy / 1000)
    # axes.set_ylim(-2000, 6000)  # KA, will need to update later - placeholder
    # ratio = 0.4   # aspect ratio
    # axes.set_aspect(0.05)
    # axes.margins(0.5)
    # axes.set_aspect(1/axes.get_data_ratio())
    axes.set_ylim(ymin, ymax)

    # Here I chose to only mark the inlets (all indices), and not where the barrier volume is less than 0...
    # need to go back and fix
    # if np.size(inlet_idx) != 0 :
    #    axes.plot(y[np.hstack(inlet_idx)]/1000, x_s[np.hstack(inlet_idx)], 'kD')

    if make_gif:
        filename = 'GIF/year_' + str(int(t[i])) + '.png'
        fig.savefig(filename, dpi=200)

        # use imageio to create gif
        frames.append(imageio.imread(filename))
        os.remove(filename)

    axes.clear()  # alternatively, plt.pause(0.0001)

# make gif from saved png files
if make_gif:
    imageio.mimsave(name + '.gif', frames, 'GIF-FI')
    os.rmdir(directory + '/GIF')


# FIGURES TO DO
# - mean barrier width (brie and CASCADE) [future, with and without inlets]
# - barrier evolution (block diagram?, like Jaap Figure 1a) for brie and CASCADE (maybe 3 snapshots?)
# - [eventually with inlets, Qow vs Qinlet]


# 3) what is the effect of dunes on relative importance of tidal inlets (need to check brie discretization)