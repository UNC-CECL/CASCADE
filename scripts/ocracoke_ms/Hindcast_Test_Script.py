import numpy as np
import time

import matplotlib.pyplot as plt

import os
import imageio



Change_Rates = np.loadtxt('C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Revised_Offshore_Datum\\All_Shoreline_Change_Rates.csv',skiprows=1,delimiter=',')
Linear_1974 = Change_Rates[:,0]
Linear_1988 = Change_Rates[:,1]
Endpoint_1974 = Change_Rates[:,2]
Endpoint_1988 = Change_Rates[:,3]

os.chdir('C:\\Users\\frank\\PycharmProjects\\CASCADE\\Run_output')
# run_name='Wreck_ACC_RSLR3_S3' # 5 Length
run_name = "Hindcast_Test_88"  # 4 length

name_prefix = run_name
nt_run = 32
number_barrier3d_models = 39


# --------- plot ---------
output = np.load(run_name + ".npz", allow_pickle=True)
cascade = output["cascade"]
cascade = cascade[0]
b3d = cascade.barrier3d
ny = np.size(b3d)

directory = "C:\\Users\\frank\\PycharmProjects\\CASCADE\\"
# TMax_MGMT = Needed 0
# TMAX_Sim = Last simulation year of the model 99
TMax_Sim = nt_run  # Give length of simulation
roadway_management_ny = [True] * ny




"""
NOTE THAT THE BEACH REPRESENTATION IS BASED ON A MODEL SPECIFIED BEACH WIDTH. We set the beach width for the
remaining time steps after the community has been abandoned to the last managed beach width in order to not have a
huge jump in the back-barrier position in Barrier3D. OTHERWISE, it is meaningless to the dynamics Barrier3D.
"""
barrier3d = cascade.barrier3d

# set up the domain; here we just use the first grid, but that could break in future runs
BarrierLength = barrier3d[0].BarrierLength
OriginY = int(barrier3d[0].x_s_TS[0])
total_shoreline_change = cascade._brie_coupler.brie.x_s_dt
all_shoreline_change = cascade._brie_coupler.brie.x_s_save

Year_1_Shoreline_Positions = all_shoreline_change[:,1]
Year_44_Shoreline_Positions = all_shoreline_change[:,-1]
EP_Change = ((Year_44_Shoreline_Positions - Year_1_Shoreline_Positions)*-1)/31

Dif = Endpoint_1988 - EP_Change



# Positive means observed is changing faster than modeled
# Negative means model is overpredicting change on Ocracoke

plt.plot(Linear_1988, label = 'Linear')
plt.plot(Endpoint_1988, label = 'Endpoint')
plt.axhline(y = 0, color = 'r', linestyle = '-')
plt.plot(EP_Change,label= 'Modeled')
plt.legend()
plt.show()