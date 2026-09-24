# Lexi (Van Blunk) Fiegelist
# last updated: September 23, 2026

"""
plot the two equations that control dune growth rates based on user inputs
dunes grow fastest when they are at a low height and slowest when they are near their maximum height

Eqn 1: change in dune growth over time
dH/dt = r*H[1-(H/Hmax)]
----- variables ------
dH/dt - change in dune height over time
r - characteristic dune growth rate (mean of rmin and rmax)
H - current dune height
Hmax - maximum dune height

Integrate equation 1 to get Equation 2

Eqn 2: logistic curve equation
H = (Hmax*Ho) / [(Hmax-Ho)e^(-n) + Ho]
----- variables ------
H - dune height at some time
Hmax - maximum dune height
Ho - initial dune height
n - years


"""

import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams["font.size"] = 14
markers = list(Line2D.markers.keys())

# set your parameters
save_figs = True
savedir = r"C:\Users\agfig\OneDrive - University of North Carolina at Chapel Hill\UNC\data"
hmax = [3.0, 6.0]
min_r_to_plot = 0.05  # characteristic dune growth rate (NOT b3d rmin)
max_r_to_plot = 0.85  # characteristic dune growth rate (NOT b3d rmax)
interval = 0.1  # interval for r values
y_text = [0.95, 0.90]  # add more if you want to plot more hmax values
colors = ["tab:blue", "tab:pink"]  # add more if you want to plot more hmax values

# for second equation, test with different initial dune heights instead of hmax
y_text2 = [0.95, 0.90, 0.85, 0.80]  # add more if you want to plot more hmax values
colors2 = ["tab:blue", "tab:pink", "tab:green", "tab:orange"]  # add more if you want to plot more ho values
n_years = 10

# auto generated from here down
r_values = np.arange(min_r_to_plot, max_r_to_plot, interval)
r_values = np.round(r_values, 2)

# first plot: see how different characteristic growth rates impact the dune height over time
fig1 = plt.figure(figsize=[20,8])
ax1 = fig1.add_subplot(111)
c=0
for hm in hmax:
    color=colors[c]
    y=y_text[c]
    x_values = np.arange(0, hm+0.5, 0.5)
    c+=1
    m=0
    # add label for colors
    ax1.text(0.02, y=y, s="Hmax = {0}".format(hm), horizontalalignment='left',
             verticalalignment='center', transform=ax1.transAxes, c=color, size=18)
    for r in r_values:
        y_values = r * x_values * (1-(x_values/hm))
        ax1.plot(x_values, y_values,ls="solid", label=r, marker=markers[m], c=color)
        m+=1
ax1.legend(title="r values")
ax1.set_ylabel("dH/dt")
ax1.set_xlabel("dune height from 0 to Hmax")
ax1.set_xlim([0,np.max(hmax)])
# ax1.set_title(r"$\frac{dH}{dt} = rH[1-(\frac{H}{H_{max}})]$")
ax1.text(0.02, y=0.8, s=r"$\frac{dH}{dt} = rH[1-(\frac{H}{H_{max}})]$",
         transform=ax1.transAxes,
         horizontalalignment='left', verticalalignment='center', size=24)

# second plot: see how different characteristic growth rates impact the logistic curve
fig2 = plt.figure(figsize=[20,8])
ax1 = fig2.add_subplot(111)
c=0
# H = (Hmax*Ho) / [(Hmax-Ho)e^(-n) + Ho]
x_values = np.arange(0,n_years,1)
for hm in hmax:
    ho_values = np.arange(0.5, hm+0.5, 0.5)
    color=colors2[c]
    y=y_text2[c]
    c+=1
    m=0
    # add label for colors
    # ax1.text(0.02, y=y, s="Hmax = {0}".format(hm), horizontalalignment='left',
    #          verticalalignment='center', transform=ax1.transAxes, c=color, size=18)
    ax1.text(8, y=hm+0.2, s="Hmax = {0}".format(hm), horizontalalignment='left',
             verticalalignment='center', c=color, size=18)
    for ho in ho_values:
        y_values = (hm * ho) / ((hm - ho) * np.exp(-x_values) + ho)  # dune height
        ax1.plot(x_values, y_values,ls="solid", label=ho, marker=markers[m], c=color)
        m+=1
ax1.legend(title="ho values")
ax1.set_ylabel("H")
ax1.set_xlabel("years")
ax1.set_ylim([0,np.max(hmax)+0.5])
ax1.set_xlim(0,n_years)
# ax1.set_title(r"$H = \frac{H_{max}*H_o}{[(H_{max}-H_o)e^{-n} + H_o]}$", size=24)
ax1.text(6, y=1, s=r"$H = \frac{H_{max}*H_o}{[(H_{max}-H_o)e^{-n} + H_o]}$",
         horizontalalignment='left', verticalalignment='center', size=30)

if save_figs:
    fig1.savefig(os.path.join(savedir, "dune_height_change.png"))
    fig2.savefig(os.path.join(savedir, "dune_height_logistic_curve.png"))