# Lexi (Van Blunk) Fiegelist
# Sept. 24, 2026
# initial analysis used to compare results for calibration (2004-2014) and testing (2014-2024)

# calibration: compare the following between model results and observations
# - modeled dune heights to selected 2014 dune heights
# - shoreline movement
# - marsh elevations (maybe the average)
# - overwash extents/features
# - full domains side-by-side


import os
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.size"] = 14

#--------------------- calibration -------------------------
datadir = r"C:\Users\agfig\model\calibration"
file = os.path.join(datadir, "results", "calib_1.npz")
calib = np.load(file, allow_pickle=True)["cascade"][0]

# dune heights in m
d_start_num = 3
d_end_num = 25
items = d_end_num - d_start_num + 1

# model results
all_model_2014_dunes = []
for i in range(items):
    b3d = calib.barrier3d[i]  # barrier3d class for one domain
    dunes_2014 = b3d.DuneDomain[-1]  # we only care about the last result for now
    # compare the crests only?
    dunes_2014_crest = np.max(dunes_2014, axis=1)  # max dune heights in dam
    all_model_2014_dunes.append(dunes_2014_crest)
modeled_2014_dunes = np.concatenate(all_model_2014_dunes) * 10  # dune heights in m for all domains

# load the observed 2014 dunes and combine them into a single array
all_dunes = []  # reset array each year
folder = os.path.join(datadir, 'dunes_2014_final_v2_berm2pt0')  # dune heights in dam (only one row)
for d in range(d_start_num, d_end_num+1):
    dune_name = os.path.join(folder, 'domain_{0}_dunes_2014.npy'.format(d))
    dunes = np.load(dune_name)
    all_dunes.append(dunes)  # this is now a list of arrays
# combine into a sinle array
observed_2014_dunes = np.concatenate(all_dunes) * 10  # dune heights in m for all domains

# plot both
fig1 = plt.figure(figsize=[15,5])
ylim_low = -0.1
ylim_high = 3.5
plt.plot(modeled_2014_dunes, label = "model", ls="dashed")  # dune heights in m
plt.plot(observed_2014_dunes, label = "observed")
# plot the domain edges
domain_len_dam = 50
domain_edges = np.arange(0,len(modeled_2014_dunes)+domain_len_dam,domain_len_dam)
plt.vlines(domain_edges, ymin=ylim_low, ymax=ylim_high, colors="black", ls="dotted")
# other plot features
plt.legend(loc="upper right")
plt.ylim([ylim_low, ylim_high])  # m
plt.ylabel("height above berm elev (dam)")
plt.xlabel("alongshore length south to north (dam)")
# plt.xlim([0,len(all_dunes_concat)])
# plot domain labels as text at the top of the figure, halfway between the domain limits
start_text = domain_len_dam / 2  # first label at 25
end_text = (items+1)*domain_len_dam - start_text
text_pos = np.arange(start_text, end_text, domain_len_dam)
d_text = 1
for t in range(len(text_pos)):
    # plt.text(x=text_pos[t], y=0.58, s="{0}".format(d_text), va="center", ha="center")
    plt.text(x=text_pos[t], y=3.4, s="{0}".format(d_text), va="center", ha="center")
    d_text += 1
