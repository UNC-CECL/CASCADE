import shutil
import os
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from cascade import Cascade

os.chdir("/Users/ceclmac/PycharmProjects/CASCADE")


NT = 180

# test to make sure Marsh_Dynamics can be run
def run_cascade_marsh_dynamics(datadir):
    # for data_file in datadir.iterdir():
    #    shutil.copy(data_file, ".")
    cascade = Cascade(
        datadir,
        name="test_marsh_dynamics",
        storm_file="cascade-default-storms.npy",
        elevation_file="Marsh_Test_Inputs/InitElevHog.npy",
        dune_file="Marsh_Test_Inputs/barrier3d-dunes.npy",
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

    for _ in range(NT - 1):
        cascade.update()
        if cascade.b3d_break:
            break

    return cascade


# test to make sure B3D has topography replaced by PyBMFT


def check_b3d_topography_replacement_from_PYBMFT(datadir):
    """
    check that the barrier3d elevation is being updated by PyBMFT
    """
    cascade = Cascade(
        datadir,
        name="test_b3d_topo_replacement",
        storm_file="cascade-default-storms.npy",
        elevation_file="Marsh_Test_Inputs/InitElevHog.npy",
        dune_file="Marsh_Test_Inputs/barrier3d-dunes.npy",
        parameter_file="Alongshore_Test-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.004,
        sea_level_rise_constant=True,
        background_erosion=0.0,
        alongshore_section_count=1,
        time_step_count=6,
        min_dune_growth_rate=0.55,
        max_dune_growth_rate=0.95,  # rave = 0.75
        num_cores=1,
        roadway_management_module=False,
        alongshore_transport_module=False,
        beach_nourishment_module=False,
        community_economics_module=False,  # no community dynamics
        marsh_dynamics=True,
    )

    PyBMFT_elevation = [[] for i in range(cascade._ny)]
    b3d_topography = [[] for i in range(cascade._ny)]

    for i in range(4):
        cascade.update()
        for iB3D in range(cascade._ny):
            if cascade._ny < 2:
                BB_transect = np.flip(
                    cascade._bmft_coupler._bmftc[0].elevation[
                        cascade._bmft_coupler._bmftc[0].startyear + i - 1,
                        int(
                            cascade._bmft_coupler._bmftc[0].Marsh_edge[
                                cascade._bmft_coupler._bmftc[0].startyear + i
                            ]
                        ) :,
                    ]
                )
                PyBMFT_elevation.append(BB_transect)
                b3d_topography.append(
                    np.mean(cascade.barrier3d[iB3D].DomainTS[i], axis=1) * 10
                )
            else:
                BB_transect = np.flip(
                    cascade._bmft_coupler._bmftc[0].elevation[
                        cascade._bmft_coupler._bmftc[iB3D].startyear + i - 1,
                        int(
                            cascade._bmft_coupler._bmftc[0].Marsh_edge[
                                cascade._bmft_coupler._bmftc[iB3D].startyear + i
                            ]
                        ) :,
                    ]
                )
                PyBMFT_elevation[iB3D].append(BB_transect)
                b3d_topography[iB3D].append(cascade.barrier3d[iB3D].DomainTS[i])

    return b3d_topography


def check_PYBMFT_topography_replacement_from_B3D(datadir):
    """
    check that the PyBMFT elevation is being updated based on changes to B3D elevation
    """
    cascade = Cascade(
        datadir,
        name="test_PyBMFT_topo_replacement",
        storm_file="cascade-default-storms.npy",
        elevation_file="Marsh_Test_Inputs/InitElevHog.npy",
        dune_file="Marsh_Test_Inputs/barrier3d-dunes.npy",
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

    PyBMFT_elevation = [[] for i in range(cascade._ny)]
    b3d_topography = [[] for i in range(cascade._ny)]

    for i in range(7):
        cascade.update()
        for iB3D in range(cascade._ny):
            if cascade._ny < 2:
                BB_transect = np.flip(
                    cascade._bmft_coupler._bmftc[0].elevation[
                        cascade._bmft_coupler._bmftc[0].startyear + i - 1,
                        int(
                            cascade._bmft_coupler._bmftc[0].Marsh_edge[
                                cascade._bmft_coupler._bmftc[0].startyear + i
                            ]
                        ) :,
                    ]
                )
                PyBMFT_elevation.append(BB_transect)
                b3d_topography.append(
                    np.mean(cascade.barrier3d[iB3D].DomainTS[i], axis=1) * 10
                )
            else:
                BB_transect = np.flip(
                    cascade._bmft_coupler._bmftc[0].elevation[
                        cascade._bmft_coupler._bmftc[iB3D].startyear + i - 1,
                        int(
                            cascade._bmft_coupler._bmftc[0].Marsh_edge[
                                cascade._bmft_coupler._bmftc[iB3D].startyear + i
                            ]
                        ) :,
                    ]
                )
                PyBMFT_elevation[iB3D].append(BB_transect)
                b3d_topography[iB3D].append(cascade.barrier3d[iB3D].DomainTS[i])

    return (PyBMFT_elevation, b3d_topography)


def changing_B3D_elevs_on_PyBMFT(datadir):
    """
    check that the PyBMFT elevation is being updated based on large changes to B3D elevation
    """
    cascade = Cascade(
        datadir,
        name="test_PyBMFT_topo_replacement",
        storm_file="cascade-default-storms.npy",
        elevation_file="Marsh_Test_Inputs/InitElevHog.npy",
        dune_file="Marsh_Test_Inputs/barrier3d-dunes.npy",
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

    PyBMFT_elevation = [[] for i in range(cascade._ny)]
    b3d_topography = [[] for i in range(cascade._ny)]

    for i in range(7):
        cascade.update()
        for iB3D in range(cascade._ny):
            cascade.barrier3d[iB3D].DomainTS[i] = (
                cascade.barrier3d[iB3D].DomainTS[i] - 0.05
            )
            if cascade._ny < 2:
                BB_transect = np.flip(
                    cascade._bmft_coupler._bmftc[0].elevation[
                        cascade._bmft_coupler._bmftc[0].startyear + i - 1,
                        int(
                            cascade._bmft_coupler._bmftc[0].Marsh_edge[
                                cascade._bmft_coupler._bmftc[0].startyear + i
                            ]
                        ) :,
                    ]
                )
                PyBMFT_elevation.append(BB_transect)
                b3d_topography.append(
                    np.mean(cascade.barrier3d[iB3D].DomainTS[i], axis=1) * 10
                )
            else:
                BB_transect = np.flip(
                    cascade._bmft_coupler._bmftc[0].elevation[
                        cascade._bmft_coupler._bmftc[iB3D].startyear + i - 1,
                        int(
                            cascade._bmft_coupler._bmftc[0].Marsh_edge[
                                cascade._bmft_coupler._bmftc[iB3D].startyear + i
                            ]
                        ) :,
                    ]
                )
                PyBMFT_elevation[iB3D].append(BB_transect)
                b3d_topography[iB3D].append(cascade.barrier3d[iB3D].DomainTS[i])

    return (PyBMFT_elevation, b3d_topography)


def changing_PyBMFT_elevs_on_B3d(datadir):
    """
    check that the B3D elevation is being updated based on large changes to PyBMFT elevation
    """
    cascade = Cascade(
        datadir,
        name="test_PyBMFT_topo_replacement",
        storm_file="cascade-default-storms.npy",
        elevation_file="Marsh_Test_Inputs/InitElevHog.npy",
        dune_file="Marsh_Test_Inputs/barrier3d-dunes.npy",
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

    PyBMFT_elevation = [[] for i in range(cascade._ny)]
    b3d_topography = [[] for i in range(cascade._ny)]

    for i in range(6):
        cascade.update()
        for iB3D in range(cascade._ny):
            cascade._bmft_coupler._bmftc[0].elevation[
                cascade._bmft_coupler._bmftc[0].startyear + i - 1,
                int(
                    cascade._bmft_coupler._bmftc[0].Marsh_edge[
                        cascade._bmft_coupler._bmftc[0].startyear + i
                    ]
                ) :,
            ] += 1.75

            if cascade._ny < 2:
                BB_transect = np.flip(
                    cascade._bmft_coupler._bmftc[0].elevation[
                        cascade._bmft_coupler._bmftc[0].startyear + i - 1,
                        int(
                            cascade._bmft_coupler._bmftc[0].Marsh_edge[
                                cascade._bmft_coupler._bmftc[0].startyear + i
                            ]
                        ) :,
                    ]
                )
                PyBMFT_elevation.append(BB_transect)
                b3d_topography.append(
                    np.mean(cascade.barrier3d[iB3D].DomainTS[i], axis=1) * 10
                )
            else:
                BB_transect = np.flip(
                    cascade._bmft_coupler._bmftc[0].elevation[
                        cascade._bmft_coupler._bmftc[iB3D].startyear + i - 1,
                        int(
                            cascade._bmft_coupler._bmftc[0].Marsh_edge[
                                cascade._bmft_coupler._bmftc[iB3D].startyear + i
                            ]
                        ) :,
                    ]
                )
                PyBMFT_elevation[iB3D].append(BB_transect)
                b3d_topography[iB3D].append(cascade.barrier3d[iB3D].DomainTS[i])

    return (PyBMFT_elevation, b3d_topography)
