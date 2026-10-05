import copy

import numpy as np
import time

import matplotlib.pyplot as plt

import os
import imageio


os.chdir('C:\\Users\\frank\\PycharmProjects\\CASCADE\\Run_output')


#run_name_batch='OCR_1997_2020_Hindcast_Final'
run_name_batch='OCR_1997_2020_Hindcast_Final'

nt_run = 23
number_barrier3d_models = 70
buffer_length = 15


# --------- plot ---------
output = np.load(run_name_batch + ".npz", allow_pickle=True)
cascade = output["cascade"]
cascade = cascade[0]
b3d = cascade.barrier3d[48]
ny = np.size(b3d)
roads = cascade.roadways[48]
dunes = b3d.DuneDomain


transects_ts = []
for k in range(cascade._nt):
    dune_temp = np.mean(dunes[k],axis=0)
    dune_temp = dune_temp+ (b3d.BermEl )
    interior_mean = np.mean(b3d.DomainTS[k], axis=1)
    if k == 0:
        dummy_cells = b3d.ShorelineChangeTS[k]
    elif k == 22:
        dummy_cells = abs(np.sum(b3d.ShorelineChangeTS))
    else:
        dummy_cells = abs(np.sum(b3d.ShorelineChangeTS[:k+1]))
    if dummy_cells > 0:
        add_cells = np.zeros(int(dummy_cells))
        dune_temp = np.append(add_cells,dune_temp)
    final_transect = np.append(dune_temp,interior_mean)
    transects_ts.append(copy.deepcopy(final_transect*10))

#transect_10 = b3d.DomainTS[10][:,10]

#plt.plot(transects_ts[0][0:40],label='0')
'''ax = plt.gca()
plt.plot(transects_ts[12][0:35],label='12')
plt.plot(transects_ts[13][0:35],label='13')
ax.set_ylim([-2,4])
plt.title('Dune Rebuilding')
plt.ylabel('')
plt.xlabel('Cross-shore position')
plt.legend()
#plt.savefig(fname=('C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Figures\\Hindcast_Rebuilding_47.eps'),format='eps')
plt.show()'''


ax = plt.gca()
plt.axhline(y = 2.2, color = 'k', linestyle = '--')
plt.plot(transects_ts[22][0:35],label='22')
#plt.plot(transects_ts[15][0:40],label='15')
ax.set_ylim([-2,3])
plt.title('Dune Rebuilding')
plt.ylabel('')
plt.xlabel('Cross-shore position')
plt.legend()
plt.savefig(fname=('C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Figures\\Hindcast_Rebuilding_48_End.eps'),format='eps')
plt.show()

x = 10


transect = b3d.DomainTS[0][:,0]