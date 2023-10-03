import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
from tests import test_marsh_dynamics
from numpy.testing import assert_array_almost_equal

from cascade import Cascade

os.chdir("/Users/ceclmac/PycharmProjects/CASCADE")
datadir = "B3D_Inputs/"

########
# Turn on marsh dynamics and test that Cascade doesn't crash
test_marsh_dynamics.run_cascade_marsh_dynamics(datadir)

########
# Make sure PyBMFT can replace topography from B3D
Test_b3d_topography = test_marsh_dynamics.check_b3d_topography_replacement_from_PYBMFT(
    datadir
)

Z0 = Test_b3d_topography[0]
Z1 = Test_b3d_topography[1]
Z2 = Test_b3d_topography[2]
Z3 = Test_b3d_topography[3]
Z4 = Test_b3d_topography[4]

plt.plot(Z1, label="TS 1")
plt.plot(Z2, label="TS 2")
plt.plot(Z3, label="TS 3")
plt.plot(Z4, label="TS 4")
plt.title("Changes to B3D Elevation from coupling to PyBMFT")
plt.xlabel("B3D Length (dm)")
plt.ylabel("B3D Elevation (m)")
plt.legend(loc="upper right")
plt.show()

###########
# Test to make sure PyBMFT correctly initalizes topography from B3D
(
    Test_Output_PyBMFT,
    Test_b3d_topography,
) = test_marsh_dynamics.check_PYBMFT_topography_replacement_from_B3D(datadir)

Z1 = Test_Output_PyBMFT[1]
B3D_Elev_1 = Test_b3d_topography[1]
Big_B3D_Elev = []


for i in range(len(B3D_Elev_1)):
    for j in range(10):
        Big_B3D_Elev.append(B3D_Elev_1[i])

Z1 -= Z1[0]
Z_Compare = Z1[158:]

plt.plot(Z_Compare, label="PyBMFT")
plt.plot(Big_B3D_Elev, label="B3D")

plt.title("Changes to PyBMFT Elevation from Altering B3D")
plt.xlabel("Length (m)")
plt.ylabel("Elevation (m)")
plt.legend(loc="upper right")
plt.show()

########
# Test that PyBMFT reflects changes in B3D elevation
(
    Test_Output_PyBMFT,
    Test_b3d_topography,
) = test_marsh_dynamics.changing_B3D_elevs_on_PyBMFT(datadir)

Z1 = Test_Output_PyBMFT[1]
Z2 = Test_Output_PyBMFT[2]
Z3 = Test_Output_PyBMFT[3]
Z4 = Test_Output_PyBMFT[4]
Z5 = Test_Output_PyBMFT[5]
Z6 = Test_Output_PyBMFT[6]


plt.plot(Z4, label="TS 4")
plt.plot(Z5, label="TS 5")
plt.plot(Z6, label="TS 6")

plt.title("Changes to PyBMFT Elevation from Altering B3D")
plt.xlabel("PyBMFT Length (dm)")
plt.ylabel("PyBMFT Elevation (m)")
plt.legend(loc="upper right")
plt.show()

########
# Test that B3D reflects changes in PyBMFT elevation
(
    Test_Output_PyBMFT,
    Test_b3d_topography,
) = test_marsh_dynamics.changing_PyBMFT_elevs_on_B3d(datadir)

Z1 = Test_b3d_topography[1]
Z2 = Test_b3d_topography[2]
Z3 = Test_b3d_topography[3]
Z4 = Test_b3d_topography[4]
Z5 = Test_b3d_topography[5]

plt.plot(Z1, label="TS 1")
plt.plot(Z2, label="TS 2")
plt.plot(Z3, label="TS 3")
plt.plot(Z4, label="TS 4")
plt.plot(Z5, label="TS 5")
plt.title("Changes to B3D Elevation from Altering PyBMFT")
plt.xlabel("B3D Length (dm)")
plt.ylabel("B3D Elevation (m)")
plt.legend(loc="upper right")
plt.show()
