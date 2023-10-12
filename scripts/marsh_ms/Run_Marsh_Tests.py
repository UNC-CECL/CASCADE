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
'''
test_marsh_dynamics.run_cascade_marsh_dynamics(datadir)
print('Test 1 complete')
'''
########

# Make sure PyBMFT can replace topography from B3D
'''
Test_b3d_topography = test_marsh_dynamics.check_b3d_topography_replacement_from_PYBMFT(
    datadir
)
'''
#Z0 = Test_b3d_topography[0]
#Z1 = Test_b3d_topography[1]
#Z2 = Test_b3d_topography[2]
#Z3 = Test_b3d_topography[3]
#Z4 = Test_b3d_topography[4]

#plt.plot(Z1, label="TS 1")
#plt.plot(Z2, label="TS 2")
#plt.plot(Z3, label="TS 3")
#plt.plot(Z4, label="TS 4")
#plt.title("Changes to B3D Elevation from coupling to PyBMFT")
#plt.xlabel("B3D Length (dm)")
#plt.ylabel("B3D Elevation (m)")
#plt.legend(loc="upper right")
#plt.show()
print('Test 2 complete')

###########
# Test to make sure PyBMFT correctly initalizes topography from B3D
b3d_grids = 1
(
    Test_Output_PyBMFT,
    Test_b3d_topography,
) = test_marsh_dynamics.check_PYBMFT_topography_replacement_from_B3D(datadir,b3d_grids)

b3d_1 = Test_b3d_topography[1]

pymft_1 = Test_Output_PyBMFT[1]
pymft_2 = Test_Output_PyBMFT[2]
pymft_3 = Test_Output_PyBMFT[3]
pymft_4 = Test_Output_PyBMFT[4]
pymft_5 = Test_Output_PyBMFT[5]
pymft_6 = Test_Output_PyBMFT[6]
pymft_7 = Test_Output_PyBMFT[7]
pymft_8 = Test_Output_PyBMFT[8]
pymft_20 = Test_Output_PyBMFT[20]



plt.plot(pymft_1)
plt.plot(pymft_2)
plt.plot(pymft_3)
plt.plot(pymft_4)
plt.plot(pymft_5)
plt.plot(pymft_6)
plt.plot(pymft_7)
plt.plot(pymft_8)
plt.plot(pymft_20)
plt.show()

plt.plot(b3d_1)
plt.show()
plt.plot(pymft_1)
plt.show()

b3d_6 = Test_b3d_topography[6]
pymft_5 = Test_Output_PyBMFT[5]

plt.plot(b3d_6)
plt.show()
plt.plot(pymft_5)
plt.show()

b3d_12 = Test_b3d_topography[12]
pymft_11 = Test_Output_PyBMFT[11]

plt.plot(b3d_12)
plt.show()
plt.plot(pymft_11)
plt.show()

plt.plot(pymft_1, label='1')
plt.plot(pymft_8, label = '8')
plt.plot(pymft_20, label = '20')
plt.legend(loc="upper right")
plt.show()

b3d_1 = Test_b3d_topography[1]
b3d_8 = Test_b3d_topography[8]
b3d_20 = Test_b3d_topography[20]

plt.plot(b3d_1, label='1')
plt.plot(b3d_8, label = '8')
plt.plot(b3d_20, label = '20')
plt.legend(loc="upper right")
plt.show()

if b3d_grids < 2:
    Z1 = Test_Output_PyBMFT[1]
    B3D_Elev_1 = Test_b3d_topography[1]
    Big_B3D_Elev = []
    for i in range(len(B3D_Elev_1)):
        for j in range(10):
            Big_B3D_Elev.append(B3D_Elev_1[i])


#if b3d_grids < 2:
#    Z1 -= Z1[0]
#    Z_Compare = Z1[158:]
    #plt.plot(Z_Compare, label="PyBMFT")
    #plt.plot(Big_B3D_Elev, label="B3D")
    #plt.title("PyBMFT Elevation vs B3D Elevation")
    #plt.xlabel("Length (m)")
    #plt.ylabel("Elevation (m)")
    #plt.legend(loc="upper right")
    #plt.show()
'''
# print('Test 3 complete')
'''
########
# Test that PyBMFT reflects changes in B3D elevation
'''
(
    Test_Output_PyBMFT,
    Test_b3d_topography,
) = test_marsh_dynamics.changing_B3D_elevs_on_PyBMFT(datadir)

#Z1 = Test_Output_PyBMFT[1]
#Z2 = Test_Output_PyBMFT[2]
#Z3 = Test_Output_PyBMFT[3]
#Z4 = Test_Output_PyBMFT[4]
#Z5 = Test_Output_PyBMFT[5]
#Z6 = Test_Output_PyBMFT[6]


#plt.plot(Z4, label="TS 4")
#plt.plot(Z5, label="TS 5")
#plt.plot(Z6, label="TS 6")

#plt.title("Changes to PyBMFT Elevation from Altering B3D")
#plt.xlabel("PyBMFT Length (dm)")
#plt.ylabel("PyBMFT Elevation (m)")
#plt.legend(loc="upper right")
#plt.show()

print('Test 4 complete')


########
# Test that B3D reflects changes in PyBMFT elevation

(
    Test_Output_PyBMFT,
    Test_b3d_topography,
) = test_marsh_dynamics.changing_PyBMFT_elevs_on_B3d(datadir)

#Z1 = Test_b3d_topography[1]
Z2 = Test_b3d_topography[2]
Z3 = Test_b3d_topography[3]
Z4 = Test_b3d_topography[4]
Z5 = Test_b3d_topography[5]

#plt.plot(Z1, label="TS 1")
plt.plot(Z2, label="TS 2")
plt.plot(Z3, label="TS 3")
plt.plot(Z4, label="TS 4")
plt.plot(Z5, label="TS 5")
plt.title("Changes to B3D Elevation from Altering PyBMFT")
plt.xlabel("B3D Length (dm)")
plt.ylabel("B3D Elevation (m)")
plt.legend(loc="upper right")
plt.show()

print('Test 5 complete')
'''