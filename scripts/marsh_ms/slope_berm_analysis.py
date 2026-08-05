import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import lines


def calculate_slope(x1, x2, elev):
    y2 = elev[x2]
    y1 = elev[x1]
    slope = (y2-y1)/(x2-x1)
    return round(slope,3), y1, y2

# set directory to the text files folder
directory = r'C:\Users\agfig\OneDrive - University of North Carolina at Chapel Hill\UNC\data\DEM\berm_slope_analysis'
file_2004 = os.path.join(directory, "2004_transects.xls")
file_2014 = os.path.join(directory, "2014_transects.xls")
file_2020 = os.path.join(directory, "2020_transects.xls")
df_2004 = pd.read_excel(file_2004)
df_2014 = pd.read_excel(file_2014)
df_2020 = pd.read_excel(file_2020)

domains = np.arange(26,2,-1)
transects = np.unique(df_2004["transectID"].values)

# plot and save transect figures
plt.rcParams["font.size"] = 14
save_on = False
save_slopes_on = True
d = 0
beach_range = 15  # points which are spaced every 1 meter
beach_slopes = []

savedir = os.path.join(directory, "results/odd_only")
if not os.path.exists(savedir):
    os.makedirs(savedir)

for t in transects:
    if t % 2 != 0:
        # get transect elevations
        elev_2004 = df_2004[df_2004["transectID"] == t].sort_values(by='pointID').elevation.values
        elev_2014 = df_2014[df_2014["transectID"] == t].sort_values(by='pointID').elevation.values
        elev_2020 = df_2020[df_2020["transectID"] == t].sort_values(by='pointID').elevation.values
        # tide line is about halfway at each transect (except at the south end which is higher transect ID)
        tide_2004 = elev_2004[250]
        tide_2014 = elev_2014[250]
        tide_2020 = elev_2020[250]
        # assume beach elevation for streamlined beach slope calcs
        beach_elev = 1
        # find the first position of the value closest to 1 in the elev arrays
        dif_2014 = elev_2014 - 1
        # first value closest to 1
        first_beach_2014 = dif_2014[dif_2014>0][0]
        # index position of that value
        beach_2014_pos = np.where(dif_2014==first_beach_2014)[0][0]
        # use the surrounding 50 cells to calc beach slope
        pos1 = beach_2014_pos + beach_range
        pos2 = beach_2014_pos - beach_range
        # calculate slope of specific year
        beach_slope, y1, y2 = calculate_slope(x1=pos1, x2=pos2, elev=elev_2014)
        beach_slopes.append(beach_slope)
        # create the plot
        fig = plt.subplots(1, 1, figsize=(12, 6))
        # plot elevations
        plt.plot(elev_2004, label="2004")
        plt.plot(elev_2014, ls="dotted", label="2014")
        plt.plot(elev_2020, ls="dotted", label="2020")
        # plot approximate beach elevation and values used to calculate the slope
        plt.hlines(y=beach_elev, xmin=0, xmax=500, color='darkgrey', linestyle='solid', linewidth=0.25)
        plt.plot([pos1, pos2], [y1, y2], color="magenta", linewidth=0.5)
        plt.scatter([pos1, pos2], [y1, y2], color="magenta")
        # plot labels
        plt.text(10, 4.5, "beach slope = {0}".format(beach_slope))
        plt.title("transect: {0}, domain: {1}".format(t, domains[d]))
        plt.ylabel("elevation [m NAVD88]")
        plt.xlabel("alongshore position [m]")
        plt.ylim(-5, 5)
        plt.xlim(0, 500)
        # plt.hlines(y=1.95, xmin=0, xmax=500, color='magenta', linestyle='solid', linewidth=0.5)
        plt.legend(loc="upper right")
        # save and close figure
        if save_on:
            plt.savefig(os.path.join(savedir, "transect_{0}.png".format(t)))
        plt.close()
        # move to next domain
        d += 1

# plot the beach slopes per domain
fig2 = plt.figure(figsize=(15, 6))
plt.scatter(domains, beach_slopes)
avg_slope = np.mean(beach_slopes)
plt.hlines(y=avg_slope, xmin=min(domains), xmax=max(domains), linestyle="dashed", linewidth=0.75)
plt.ylabel("beach slopes")
plt.xlabel("domain")
plt.ylim(0, 0.1)
plt.text(min(domains), 0.095, "average slope = {0}".format(avg_slope))
plt.xticks(ticks=domains, labels=domains)
plt.grid(visible=True, axis="x", color="lightgrey", linewidth=0.5, linestyle="dotted")
if save_slopes_on:
    plt.savefig(os.path.join(savedir, "domain_beach_slopes.png"))
plt.show()