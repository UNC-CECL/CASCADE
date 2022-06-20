import numpy as np
from matplotlib import pyplot as plt
from lexi_b3d import Barrier3d

datadir = "tests/test_params/"
b3d = Barrier3d.from_yaml(datadir)
# increase time step
for time_step in range(1, b3d._TMAX):
    b3d.update()  # update the model's main time loop
    b3d.update_dune_domain()  # now update the dune domain and increase time by one year
    # Print time step to screen
    print("\r", "Time Step: ", time_step, end="")

#
# # discharges = []
# # for s in range(18):
# #     discharges.append(np.load("C:/Users/Lexi/Documents/Research/Outwasher/b3d_discharges_{0}.npy".format(s)))
# discharges = np.load("C:/Users/Lexi/Documents/Research/Outwasher/b3d_discharges_10.npy")
# # plot the initial full domain before sediment movement
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111)
# mat = ax1.matshow(
#     discharges[0],
#     origin="upper",
#     cmap="Greens",
#     vmin=0, vmax=0.25,
# )
# fig1.colorbar(mat)
# ax1.set_title("Initial Discharge $(dam)$")
# ax1.set_ylabel("barrier width (dam)")
# ax1.set_xlabel("barrier length (dam)")
# plt.show()