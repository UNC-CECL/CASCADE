import shutil
import os
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from cascade import Cascade

# os.chdir("/Users/ceclmac/PycharmProjects/CASCADE")  # Note from KA: eventually this will go away
# NT = 150 # KA: this too

# test initializations -- coupled and not coupled with human dynamics and AST
def initialize_cascade_no_human_dynamics_ast_marsh(datadir):
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
        alongshore_section_count=6,
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

def initialize_cascade_human_dynamics_marsh_no_ast(datadir):
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


# def test_marsh_update_no_ast(tmp_path, datadir, monkeypatch):
#     """
#     test to make sure B3D has topography replaced by PyBMFT correctly
#     """
#     monkeypatch.chdir(tmp_path)
#     CASCADE_MARSH_MODEL = initialize_cascade_no_human_dynamics_marsh_no_ast(datadir)
#
#     PyBMFT_elevation = [[] for i in range(CASCADE_MARSH_MODEL._ny)]
#     b3d_topography = [[] for i in range(CASCADE_MARSH_MODEL._ny)]
#
#     for i in range(4):
#         CASCADE_MARSH_MODEL.update()
#         for iB3D in range(CASCADE_MARSH_MODEL._ny):
#             if CASCADE_MARSH_MODEL._ny < 2:
#                 BB_transect = np.flip(
#                     CASCADE_MARSH_MODEL._bmft_coupler._bmftc[0].elevation[
#                         CASCADE_MARSH_MODEL._bmft_coupler._bmftc[0].startyear + i - 1,
#                         int(
#                             CASCADE_MARSH_MODEL._bmft_coupler._bmftc[0].Marsh_edge[
#                                 CASCADE_MARSH_MODEL._bmft_coupler._bmftc[0].startyear + i
#                             ]
#                         ) :,
#                     ]
#                 )
#                 PyBMFT_elevation.append(BB_transect)
#                 b3d_topography.append(
#                     np.mean(CASCADE_MARSH_MODEL.barrier3d[iB3D].DomainTS[i], axis=1) * 10
#                 )
#             else:
#                 BB_transect = np.flip(
#                     CASCADE_MARSH_MODEL._bmft_coupler._bmftc[0].elevation[
#                         CASCADE_MARSH_MODEL._bmft_coupler._bmftc[iB3D].startyear + i - 1,
#                         int(
#                             CASCADE_MARSH_MODEL._bmft_coupler._bmftc[0].Marsh_edge[
#                                 CASCADE_MARSH_MODEL._bmft_coupler._bmftc[iB3D].startyear + i
#                             ]
#                         ) :,
#                     ]
#                 )
#                 PyBMFT_elevation[iB3D].append(BB_transect)
#                 b3d_topography[iB3D].append(CASCADE_MARSH_MODEL.barrier3d[iB3D].DomainTS[i])
#
#     return b3d_topography
#
#
# def check_PYBMFT_topography_replacement_from_B3D(datadir,b3d_grids=1):
#     """
#     check that the PyBMFT elevation is being updated based on changes to B3D elevation
#     """
#     cascade = Cascade(
#         datadir,
#         name="test_PyBMFT_topo_replacement",
#         storm_file="cascade-default-storms.npy",
#         elevation_file="Marsh_Test_Inputs/InitElevHog.npy",
#         dune_file="Marsh_Test_Inputs/barrier3d-dunes.npy",
#         parameter_file="Alongshore_Test-parameters.yaml",
#         wave_height=1,
#         wave_period=7,
#         wave_asymmetry=0.8,
#         wave_angle_high_fraction=0.2,
#         sea_level_rise_rate=0.004,
#         sea_level_rise_constant=True,
#         background_erosion=0.0,
#         alongshore_section_count=b3d_grids,
#         time_step_count=NT,
#         min_dune_growth_rate=0.55,
#         max_dune_growth_rate=0.95,  # rave = 0.75
#         num_cores=1,
#         roadway_management_module=False,
#         alongshore_transport_module=False,
#         beach_nourishment_module=False,
#         community_economics_module=False,  # no community dynamics
#         marsh_dynamics=True,
#     )
#
#     PyBMFT_elevation = [[] for i in range(cascade._ny)]
#     b3d_topography = [[] for i in range(cascade._ny)]
#
#     for i in range(20):
#         cascade.update()
#         for iB3D in range(cascade._ny):
#             if cascade._ny < 2:
#                 BB_transect = np.flip(
#                     cascade._bmft_coupler._bmftc[0].elevation[
#                         cascade._bmft_coupler._bmftc[0].startyear + i - 1,
#                         int(
#                             cascade._bmft_coupler._bmftc[0].Marsh_edge[
#                                 cascade._bmft_coupler._bmftc[0].startyear + i
#                             ]
#                         ) :,
#                     ]
#                 )
#                 PyBMFT_elevation.append(BB_transect)
#                 b3d_topography.append(
#                     np.mean(cascade.barrier3d[iB3D].DomainTS[i], axis=1) * 10
#                 )
#             else:
#                 BB_transect = np.flip(
#                     cascade._bmft_coupler._bmftc[0].elevation[
#                         cascade._bmft_coupler._bmftc[iB3D].startyear + i - 1,
#                         int(
#                             cascade._bmft_coupler._bmftc[0].Marsh_edge[
#                                 cascade._bmft_coupler._bmftc[iB3D].startyear + i
#                             ]
#                         ) :,
#                     ]
#                 )
#                 PyBMFT_elevation[iB3D].append(BB_transect)
#                 b3d_topography[iB3D].append(cascade.barrier3d[iB3D].DomainTS[i])
#
#     return (PyBMFT_elevation, b3d_topography)
#
#
# def changing_B3D_elevs_on_PyBMFT(datadir):
#     """
#     check that the PyBMFT elevation is being updated based on large changes to B3D elevation
#     """
#     cascade = Cascade(
#         datadir,
#         name="test_PyBMFT_topo_replacement",
#         storm_file="cascade-default-storms.npy",
#         elevation_file="Marsh_Test_Inputs/InitElevHog.npy",
#         dune_file="Marsh_Test_Inputs/barrier3d-dunes.npy",
#         parameter_file="Alongshore_Test-parameters.yaml",
#         wave_height=1,
#         wave_period=7,
#         wave_asymmetry=0.8,
#         wave_angle_high_fraction=0.2,
#         sea_level_rise_rate=0.004,
#         sea_level_rise_constant=True,
#         background_erosion=0.0,
#         alongshore_section_count=1,
#         time_step_count=NT,
#         min_dune_growth_rate=0.55,
#         max_dune_growth_rate=0.95,  # rave = 0.75
#         num_cores=1,
#         roadway_management_module=False,
#         alongshore_transport_module=False,
#         beach_nourishment_module=False,
#         community_economics_module=False,  # no community dynamics
#         marsh_dynamics=True,
#     )
#
#     PyBMFT_elevation = [[] for i in range(cascade._ny)]
#     b3d_topography = [[] for i in range(cascade._ny)]
#
#     for i in range(7):
#         cascade.update()
#         for iB3D in range(cascade._ny):
#             cascade.barrier3d[iB3D].DomainTS[i] = (
#                 cascade.barrier3d[iB3D].DomainTS[i] - 0.05
#             )
#             if cascade._ny < 2:
#                 BB_transect = np.flip(
#                     cascade._bmft_coupler._bmftc[0].elevation[
#                         cascade._bmft_coupler._bmftc[0].startyear + i - 1,
#                         int(
#                             cascade._bmft_coupler._bmftc[0].Marsh_edge[
#                                 cascade._bmft_coupler._bmftc[0].startyear + i
#                             ]
#                         ) :,
#                     ]
#                 )
#                 PyBMFT_elevation[iB3D].append(BB_transect)
#                 b3d_topography[iB3D].append(
#                     np.mean(cascade.barrier3d[iB3D].DomainTS[i], axis=1) * 10
#                 )
#             else:
#                 BB_transect = np.flip(
#                     cascade._bmft_coupler._bmftc[0].elevation[
#                         cascade._bmft_coupler._bmftc[iB3D].startyear + i - 1,
#                         int(
#                             cascade._bmft_coupler._bmftc[0].Marsh_edge[
#                                 cascade._bmft_coupler._bmftc[iB3D].startyear + i
#                             ]
#                         ) :,
#                     ]
#                 )
#                 PyBMFT_elevation[iB3D].append(BB_transect)
#                 b3d_topography[iB3D].append(cascade.barrier3d[iB3D].DomainTS[i])
#
#     return (PyBMFT_elevation, b3d_topography)
#
#
# def changing_PyBMFT_elevs_on_B3d(datadir):
#     """
#     check that the B3D elevation is being updated based on large changes to PyBMFT elevation
#     """
#     cascade = Cascade(
#         datadir,
#         name="test_PyBMFT_topo_replacement",
#         storm_file="cascade-default-storms.npy",
#         elevation_file="Marsh_Test_Inputs/InitElevHog.npy",
#         dune_file="Marsh_Test_Inputs/barrier3d-dunes.npy",
#         parameter_file="Alongshore_Test-parameters.yaml",
#         wave_height=1,
#         wave_period=7,
#         wave_asymmetry=0.8,
#         wave_angle_high_fraction=0.2,
#         sea_level_rise_rate=0.004,
#         sea_level_rise_constant=True,
#         background_erosion=0.0,
#         alongshore_section_count=2,
#         time_step_count=150,
#         min_dune_growth_rate=0.55,
#         max_dune_growth_rate=0.95,  # rave = 0.75
#         num_cores=1,
#         roadway_management_module=False,
#         alongshore_transport_module=False,
#         beach_nourishment_module=False,
#         community_economics_module=False,  # no community dynamics
#         marsh_dynamics=True,
#     )
#
#     PyBMFT_elevation = [[] for i in range(cascade._ny)]
#     b3d_topography = [[] for i in range(cascade._ny)]
#
#     for i in range(6):
#         cascade.update()
#         for iB3D in range(cascade._ny):
#             cascade._bmft_coupler._bmftc[0].elevation[
#                 cascade._bmft_coupler._bmftc[0].startyear + i - 1,
#                 int(
#                     cascade._bmft_coupler._bmftc[0].Marsh_edge[
#                         cascade._bmft_coupler._bmftc[0].startyear + i
#                     ]
#                 ) :,
#             ] += 1.75
#
#             if cascade._ny < 2:
#                 BB_transect = np.flip(
#                     cascade._bmft_coupler._bmftc[0].elevation[
#                         cascade._bmft_coupler._bmftc[0].startyear + i - 1,
#                         int(
#                             cascade._bmft_coupler._bmftc[0].Marsh_edge[
#                                 cascade._bmft_coupler._bmftc[0].startyear + i
#                             ]
#                         ) :,
#                     ]
#                 )
#                 PyBMFT_elevation[iB3D].append(BB_transect)
#                 b3d_topography[iB3D].append(
#                     np.mean(cascade.barrier3d[iB3D].DomainTS[i], axis=1) * 10
#                 )
#             else:
#                 BB_transect = np.flip(
#                     cascade._bmft_coupler._bmftc[0].elevation[
#                         cascade._bmft_coupler._bmftc[iB3D].startyear + i - 1,
#                         int(
#                             cascade._bmft_coupler._bmftc[0].Marsh_edge[
#                                 cascade._bmft_coupler._bmftc[iB3D].startyear + i
#                             ]
#                         ) :,
#                     ]
#                 )
#                 PyBMFT_elevation[iB3D].append(BB_transect)
#                 b3d_topography[iB3D].append(cascade.barrier3d[iB3D].DomainTS[i])
#
#     return (PyBMFT_elevation, b3d_topography)
