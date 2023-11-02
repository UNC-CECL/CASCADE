from cascade import Cascade
import numpy as np
import matplotlib.pyplot as plt
import os
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_equal


os.chdir("/Users/ceclmac/PycharmProjects/CASCADE")  # Note from KA: eventually this will go away
NT = 199 # KA: this too

def Compare_B3D_BMFT_Topo_No_AST_Constant_RSLR(datadir,b3d_grids=1):
    """
    check that the PyBMFT elevation is being updated based on changes to B3D elevation
    """
    cascade = Cascade(
        datadir,
        name="test_PyBMFT_topo_replacement",
        storm_file="cascade-default-storms.npy",
        elevation_file="Hog_Topo_2.npy",
        dune_file="barrier3d-default-dunes.npy",
        parameter_file="Alongshore_Test-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.006,
        sea_level_rise_constant=True,
        background_erosion=0.0,
        alongshore_section_count=b3d_grids,
        time_step_count=NT,
        min_dune_growth_rate=0.55,
        max_dune_growth_rate=0.95,  # rave = 0.75
        num_cores=1,
        roadway_management_module=False,
        alongshore_transport_module=False,
        beach_nourishment_module=False,
        community_economics_module=False,  # no community dynamics
        marsh_dynamics=True,
    )

    for i in range(150):
        cascade.update()
        print(i)

        if i > 0:
            Failed_Tests = ['Acceptable % error','Acceptable gross error']
            Test_Values = ['Acceptable % error','Acceptable gross error']
            BB_transect = np.flip(
                cascade._bmft_coupler._bmftc[0].elevation[
                cascade._bmft_coupler._bmftc[0].startyear + i - 1,
                int(
                    cascade._bmft_coupler._bmftc[0].Marsh_edge[
                        cascade._bmft_coupler._bmftc[0].startyear + i -1
                        ]
                ):,
                ]
            )

            ocean_edge = np.where(BB_transect != BB_transect[0])
            ocean_edge = ocean_edge[0][0]
            BMFT_RSLR = (cascade._bmft_coupler._bmftc[0].msl[cascade._bmft_coupler._bmftc[0].startyear - 1 +i] + cascade._bmft_coupler._bmftc[0].amp + (
                    cascade._bmft_coupler._bmftc[0].RSLRi / 1000))  # - RSLR #- 0.008
            BB_transect -= BMFT_RSLR

            B3D_Elev = (np.mean(cascade._bmft_coupler._b3d_elev_after_PyBMFT_TS[i], axis=1) * 10)
            x = np.linspace(
                1, len(B3D_Elev) * 10, num=len(B3D_Elev) * 10
            )
            xp = (
                    np.linspace(1, len(B3D_Elev), num=len(B3D_Elev))
                    * 10
            )
            xp = xp - 5
            B3D_Elev_M = np.interp(x, xp, B3D_Elev)

            if len(B3D_Elev_M) > len(BB_transect[ocean_edge:]):
                Compare_BMFT = BB_transect[ocean_edge:]
                Compare_B3D = B3D_Elev_M[:len(Compare_BMFT)]
            elif len(B3D_Elev_M) < len(BB_transect[ocean_edge:]):
                Compare_B3D = B3D_Elev_M
                Compare_BMFT = BB_transect[ocean_edge:]
                Compare_BMFT = Compare_BMFT[:len(Compare_B3D)]

            dif = abs(Compare_BMFT - Compare_B3D)
            max(dif)
            Difference_Values = np.zeros(len(dif))

            plt.plot(Compare_BMFT, label='BMFT' + str(i))
            plt.plot(Compare_B3D, label='B3D' + str(i))
            plt.legend()
            plt.show()

            for i in range(len(dif)):
                if dif[i] > 0.02:
                    Difference_Values[i] = 1

            if sum(Difference_Values) > len(dif) / 10:
                Test_Values[0] = 'Too much % error'
                print('% of incorrect transects '+str(sum(Difference_Values)/(len(dif))))

            elif sum(Difference_Values) < len(dif) / 10:
                Test_Values[0] = "Acceptable % error"

            if max(dif) > max(Compare_B3D) / 5 or max(dif) > max(Compare_BMFT) / 5:
                Test_Values[1] = 'Too much gross error'
            else:
                Test_Values[1] = 'Acceptable gross error'
            assert_equal(Test_Values,Failed_Tests)

    return

def Test_B3D_BMFT_Topo_No_AST_Accelerated_RSLR(datadir,b3d_grids=1):
    """
    check that the PyBMFT elevation is being updated based on changes to B3D elevation
    """
    cascade = Cascade(
        datadir,
        name="test_PyBMFT_topo_replacement",
        storm_file="cascade-default-storms.npy",
        elevation_file="Hog_Topo_2.npy",
        dune_file="barrier3d-default-dunes.npy",
        parameter_file="Alongshore_Test-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.25,
        sea_level_rise_rate=0.006,
        sea_level_rise_constant=False,
        background_erosion=0.0,
        alongshore_section_count=b3d_grids,
        time_step_count=NT,
        min_dune_growth_rate=0.55,
        max_dune_growth_rate=0.95,  # rave = 0.75
        num_cores=1,
        roadway_management_module=False,
        alongshore_transport_module=False,
        beach_nourishment_module=False,
        community_economics_module=False,  # no community dynamics
        marsh_dynamics=True,
    )

    for i in range(50):
        cascade.update()
        print(i)

        if i > 0:
            Failed_Tests = ['Acceptable % error','Acceptable gross error']
            Test_Values = ['Acceptable % error','Acceptable gross error']
            BB_transect = np.flip(
                cascade._bmft_coupler._bmftc[0].elevation[
                cascade._bmft_coupler._bmftc[0].startyear + i - 1,
                int(
                    cascade._bmft_coupler._bmftc[0].Marsh_edge[
                        cascade._bmft_coupler._bmftc[0].startyear + i -1
                        ]
                ):,
                ]
            )

            ocean_edge = np.where(BB_transect != BB_transect[0])
            ocean_edge = ocean_edge[0][0]
            BMFT_RSLR = (cascade._bmft_coupler._bmftc[0].msl[cascade._bmft_coupler._bmftc[0].startyear - 1 +i] + cascade._bmft_coupler._bmftc[0].amp + (
                    cascade._bmft_coupler._bmftc[0].RSLRi / 1000))  # - RSLR #- 0.008
            BB_transect -= BMFT_RSLR

            B3D_Elev = (np.mean(cascade._bmft_coupler._b3d_elev_after_PyBMFT_TS[i], axis=1) * 10)
            x = np.linspace(
                1, len(B3D_Elev) * 10, num=len(B3D_Elev) * 10
            )
            xp = (
                    np.linspace(1, len(B3D_Elev), num=len(B3D_Elev))
                    * 10
            )
            xp = xp - 5
            B3D_Elev_M = np.interp(x, xp, B3D_Elev)

            if len(B3D_Elev_M) > len(BB_transect[ocean_edge:]):
                Compare_BMFT = BB_transect[ocean_edge:]
                Compare_B3D = B3D_Elev_M[:len(Compare_BMFT)]
            elif len(B3D_Elev_M) < len(BB_transect[ocean_edge:]):
                Compare_B3D = B3D_Elev_M
                Compare_BMFT = BB_transect[ocean_edge:]
                Compare_BMFT = Compare_BMFT[:len(Compare_B3D)]

            dif = abs(Compare_BMFT - Compare_B3D)
            max(dif)
            Difference_Values = np.zeros(len(dif))

            plt.plot(Compare_BMFT, label='BMFT' + str(i))
            plt.plot(Compare_B3D, label='B3D' + str(i))
            plt.legend()
            plt.show()

            for i in range(len(dif)):
                if dif[i] > 0.02:
                    Difference_Values[i] = 1

            if sum(Difference_Values) > len(dif) / 10:
                Test_Values[0] = 'Too much % error'
                print('% of incorrect transects '+str(sum(Difference_Values)/(len(dif))))

            elif sum(Difference_Values) < len(dif) / 10:
                Test_Values[0] = "Acceptable % error"

            if max(dif) > max(Compare_B3D) / 5 or max(dif) > max(Compare_BMFT) / 5:
                Test_Values[1] = 'Too much gross error'
            else:
                Test_Values[1] = 'Acceptable gross error'
            assert_equal(Test_Values,Failed_Tests)

    return

Test_B3D_BMFT_Topo_No_AST_Accelerated_RSLR(datadir = "data/",b3d_grids=1)
