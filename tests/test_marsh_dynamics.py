import shutil
import os
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_equal

from cascade import Cascade

# os.chdir("/Users/ceclmac/PycharmProjects/CASCADE")  # Note from KA: eventually this will go away
NT = 150 # KA: this too

# test initializations -- coupled and not coupled with human dynamics and AST
def initialize_cascade_no_human_dynamics_ast_marsh(datadir):
    for data_file in datadir.iterdir():
        shutil.copy(data_file, ".")

    cascade = Cascade(
        datadir,
        name="test_marsh_dynamics",
        storm_file="cascade-default-storms.npy",  # KA: move all of these into the test directory
        elevation_file="Hog_Topo_2.npy",
        dune_file="barrier3d-dunes.npy",
        parameter_file="Alongshore_Test-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.004,
        sea_level_rise_constant=True,
        background_erosion=0.0,
        alongshore_section_count=6,
        time_step_count=150,
        min_dune_growth_rate=0.55,
        max_dune_growth_rate=0.95,  # rave = 0.75
        num_cores=1,
        roadway_management_module=False,
        alongshore_transport_module=False,
        beach_nourishment_module=False,
        community_economics_module=False,  # no community dynamics
        marsh_dynamics=True,
    )

    return cascade

def initialize_cascade_human_dynamics_marsh_no_ast(datadir):
    for data_file in datadir.iterdir():
        shutil.copy(data_file, ".")

    cascade = Cascade(
        datadir,
        name="test_marsh_dynamics",
        storm_file="cascade-default-storms.npy",  # KA: move all of these into the test directory
        elevation_file="Hog_Topo.npy",
        dune_file="barrier3d-dunes.npy",
        parameter_file="Alongshore_Test-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.004,
        sea_level_rise_constant=True,
        background_erosion=0.0,
        alongshore_section_count=1,
        time_step_count=150,
        min_dune_growth_rate=0.55,
        max_dune_growth_rate=0.95,  # rave = 0.75
        num_cores=1,
        roadway_management_module=False,
        alongshore_transport_module=False,
        beach_nourishment_module=False,
        community_economics_module=False,  # no community dynamics
        marsh_dynamics=True,
    )

    return cascade

def initialize_cascade_no_human_dynamics_marsh_no_ast(datadir):
    for data_file in datadir.iterdir():
        shutil.copy(data_file, ".")

    cascade = Cascade(
        datadir,
        name="test_marsh_dynamics",
        storm_file="cascade-default-storms.npy",  # KA: move all of these into the test directory
        elevation_file="InitElevHog.npy",
        dune_file="barrier3d-dunes.npy",
        parameter_file="Alongshore_Test-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.004,
        sea_level_rise_constant=True,
        background_erosion=0.0,
        alongshore_section_count=1,
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

    return cascade

def initialize_cascade_no_human_dynamics_marsh_no_ast_ACC_RSLR(datadir):
    for data_file in datadir.iterdir():
        shutil.copy(data_file, ".")

    cascade = Cascade(
        datadir,
        name="test_marsh_dynamics",
        storm_file="cascade-default-storms.npy",  # KA: move all of these into the test directory
        elevation_file="Hog_Topo_2.npy",
        dune_file="barrier3d-dunes.npy",
        parameter_file="Alongshore_Test-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.004,
        sea_level_rise_constant=False,
        background_erosion=0.0,
        alongshore_section_count=1,
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

    return cascade

def initialize_cascade_no_human_dynamics_no_ast_marsh_user_defined_RSLR(datadir):
    for data_file in datadir.iterdir():
        shutil.copy(data_file, ".")

    cascade = Cascade(
        datadir,
        name="test_marsh_dynamics",
        storm_file="cascade-default-storms.npy",  # KA: move all of these into the test directory
        elevation_file="Hog_Topo_2.npy",
        dune_file="barrier3d-dunes.npy",
        parameter_file="Alongshore_Test-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.004,
        sea_level_rise_constant=True,
        background_erosion=0.0,
        alongshore_section_count=1,
        time_step_count=150,
        min_dune_growth_rate=0.55,
        max_dune_growth_rate=0.95,  # rave = 0.75
        num_cores=1,
        roadway_management_module=False,
        alongshore_transport_module=False,
        beach_nourishment_module=False,
        community_economics_module=False,  # no community dynamics
        marsh_dynamics=True,
        user_inputed_RSLR=True,
        user_inputed_RSLR_rate=np.load('Low_SLR.npy'),
    )

    return cascade


# actual tests that call the functions above
def test_initialize_ast(tmp_path, datadir, monkeypatch):
    """
    test that the brie-cascade coupling hasn't been altered by the coupling with PyBMFT

    Note from KA: add a bftc.elevation (?) comparison to cascade.barrier3d.InitElev (as an example) and other variables
    that should be initialized the same

    """
    monkeypatch.chdir(tmp_path)
    CASCADE_AST_MODEL = initialize_cascade_no_human_dynamics_ast_marsh(datadir)

    # check that brie coupling is still intact
    x_t_TS, x_s_TS, x_b_TS, h_b_TS = ([] for _ in range(4))
    for iB3D in range(CASCADE_AST_MODEL.brie.ny):
        x_t_TS.append(CASCADE_AST_MODEL.barrier3d[iB3D].x_t_TS[0] * 10)
        x_s_TS.append(CASCADE_AST_MODEL.barrier3d[iB3D].x_s_TS[0] * 10)
        x_b_TS.append(CASCADE_AST_MODEL.barrier3d[iB3D].x_b_TS[0] * 10)
        h_b_TS.append(CASCADE_AST_MODEL.barrier3d[iB3D].h_b_TS[0] * 10)

    dt = CASCADE_AST_MODEL.brie.x_t - x_t_TS
    ds = CASCADE_AST_MODEL.brie.x_s - x_s_TS  # this isn't always zero due to rounding
    db = CASCADE_AST_MODEL.brie.x_b - x_b_TS
    dh = CASCADE_AST_MODEL.brie.h_b - h_b_TS

    assert_array_almost_equal([dt, ds, db, dh], np.zeros([4, 6]))

def test_B3D_BMFT_Topo_No_AST(tmp_path, datadir, monkeypatch):
    """
    Test that the PyBMFT elevation is being updated based on changes to B3D elevation
    """
    monkeypatch.chdir(tmp_path)
    cascade = initialize_cascade_human_dynamics_marsh_no_ast(datadir)

    for i in range(10):
        cascade.update()
        print(i)

        if i > 0:
            Failed_Tests = ['Acceptable % error', 'Acceptable gross error']
            Test_Values = ['Acceptable % error', 'Acceptable gross error']
            BB_transect = np.flip(
                cascade._bmft_coupler._bmftc[0].elevation[
                cascade._bmft_coupler._bmftc[0].startyear + i - 1,
                int(
                    cascade._bmft_coupler._bmftc[0].Marsh_edge[
                        cascade._bmft_coupler._bmftc[0].startyear + i - 1
                        ]
                ):,
                ]
            )

            ocean_edge = np.where(BB_transect != BB_transect[0])
            ocean_edge = ocean_edge[0][0]
            BMFT_RSLR = (cascade._bmft_coupler._bmftc[0].msl[cascade._bmft_coupler._bmftc[0].startyear - 1 + i] +
                         cascade._bmft_coupler._bmftc[0].amp + (
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

            for i in range(len(dif)):
                if dif[i] > 0.025:
                    Difference_Values[i] = 1

            if sum(Difference_Values) > len(dif) / 10:
                Test_Values[0] = 'Too much % error'
                print('% of incorrect transects ' + str(sum(Difference_Values) / (len(dif))))

            elif sum(Difference_Values) < len(dif) / 10:
                Test_Values[0] = "Acceptable % error"

            if max(dif) > max(Compare_B3D) / 5 or max(dif) > max(Compare_BMFT) / 5:
                Test_Values[1] = 'Too much gross error'
            else:
                Test_Values[1] = 'Acceptable gross error'
            assert_equal(Test_Values, Failed_Tests)

def test_B3D_BMFT_Topo_No_AST_Accelerated_RSLR(tmp_path, datadir, monkeypatch):
    """
    check that the PyBMFT elevation is being updated based on changes to B3D elevation with Accelerated RSLR rates
    """
    monkeypatch.chdir(tmp_path)
    cascade = initialize_cascade_no_human_dynamics_marsh_no_ast_ACC_RSLR(datadir)

    for i in range(10):
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


            for i in range(len(dif)):
                if dif[i] > 0.025:
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

def test_similar_RSLR_marsh_no_ast(tmp_path, datadir, monkeypatch):
    """
    check that the PyBMFT elevation is being updated based on changes to B3D elevation with Accelerated RSLR rates
    """
    monkeypatch.chdir(tmp_path)
    cascade = initialize_cascade_no_human_dynamics_no_ast_marsh_user_defined_RSLR(datadir)

    for i in range(10):
        cascade.update()
        b3d_RSLR = cascade.barrier3d[0].RSLR[i] * 10
        brie_RSLR = cascade._brie_coupler._brie.slr[i]
        bmft_RSLR = cascade._bmft_coupler._bmftc[0].RSLRi
        assert_array_almost_equal(b3d_RSLR, brie_RSLR, bmft_RSLR)
