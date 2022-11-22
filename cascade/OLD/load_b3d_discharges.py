import numpy as np
import imageio
from matplotlib import pyplot as plt
import os
from cascade.OLD.lexi_b3d import Barrier3d

# ---------------------------------------------run lexi b3d-------------------------------------------------------------
datadir = "tests/test_params/"
b3d = Barrier3d.from_yaml(datadir)
# increase time step
TMAX = 5
for time_step in range(1, TMAX):
    # Print time step to screen
    print("\r", "Time Step: ", time_step, end="")
    b3d.update()  # update the model's main time loop
    b3d.update_dune_domain()  # now update the dune domain and increase time by one year

# ---------------------------------------------load discharges and slopes-----------------------------------------------
discharges = []
for storms in range(3):
    discharges.append(np.load("C:/Users/Lexi/Documents/Research/Outwasher/b3d_discharges2/discharge_{0}.npy".format(storms)))

slopes = []
for s in range(3):
    slopes.append(np.load("C:/Users/Lexi/Documents/Research/Outwasher/b3d_slopes/slope2_{0}.npy".format(s)))

elevations = []
for st in range(3):
    elevations.append(np.load("C:/Users/Lexi/Documents/Research/Outwasher/b3d_elev_changes/elev_change_{0}.npy".format(st)))
# ---------------------------------------------make gif of discharges---------------------------------------------------
def plot_DischargeAnimation(dis, directory, TMAX, name):
    os.chdir(directory)
    newpath = name + "/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(TMAX):
        AnimateDomain = dis[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain, origin="upper", cmap="jet_r",
            vmin=0, vmax=120,
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Discharge (dam^3/hr)")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " hrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "dis_" + str(t)
        elevFig1.savefig(name)  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(TMAX):
        filename = "dis_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    # imageio.mimsave("dis.gif", frames, fps=2)
    imageio.mimsave("dis.gif", frames, "GIF-FI")
    print()
    print("[ * dishcarge GIF {0}/{1} successfully generated * ]".format(n+1, len(discharges)))

def plot_SlopesAnimation(slopes, directory, TMAX, name):
    os.chdir(directory)
    newpath = name + "/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(TMAX):
        AnimateDomain = slopes[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain, origin="upper", cmap="jet_r",
            # vmin=0, vmax=120,
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("S2 Slopes")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " hrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "slope_" + str(t)
        elevFig1.savefig(name)  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(TMAX):
        filename = "slope_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    # imageio.mimsave("dis.gif", frames, fps=2)
    imageio.mimsave("slope2.gif", frames, "GIF-FI")
    print()
    print("[ * slope GIF {0}/{1} successfully generated * ]".format(n+1, len(discharges)))

def plot_ElevChangeAnimation(elevs, directory, TMAX, name):
    os.chdir(directory)
    newpath = name + "/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    for t in range(TMAX):
        AnimateDomain = elevs[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain, origin="upper", cmap="jet_r",
            vmin=-0.008, vmax=0.008,
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Elevation Changes")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " hrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "elev_" + str(t)
        elevFig1.savefig(name)  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(TMAX):
        filename = "elev_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    # imageio.mimsave("dis.gif", frames, fps=2)
    imageio.mimsave("elevs.gif", frames, "GIF-FI")
    print()
    print("[ * elevation GIF {0}/{1} successfully generated * ]".format(n+1, len(discharges)))
#
for n in range(3):
    dis = discharges[n]
    tmax = len(dis)
    dir1 = "C:/Users/Lexi/Documents/Research/Outwasher/b3d_discharges2"
    dir2 = "C:/Users/Lexi/Documents/Research/Outwasher/b3d_slopes"
    dir3 = "C:/Users/Lexi/Documents/Research/Outwasher/b3d_elev_changes"
    name1 = "storm_{0}_discharge".format(n)
    name2 = "storm_{0}_slopes".format(n)
    name3 = "storm_{0}_elevs".format(n)
    slope = slopes[n]
    elev = elevations[n]
    # plot_DischargeAnimation(dis, dir1, TMAX=tmax, name=name1)
    # plot_SlopesAnimation(slope, dir2, TMAX=tmax, name=name2)
    plot_ElevChangeAnimation(elev, dir3, TMAX=tmax, name=name3)

