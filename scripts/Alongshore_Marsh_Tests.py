# Nov 10 2022
# Script for running and plotting an instance of Cascade with back-barrier marsh dynamics


# Required Python packages
import numpy as np
import time

import matplotlib.pyplot as plt

import os
import imageio


from cascade.cascade import Cascade

os.chdir('/Users/ceclmac/PycharmProjects/CASCADE/')


def Test_Alongshore_Runs(
    nt,
    name,
    storm_file,
    alongshore_section_count,
    num_cores,
    rmin,  # the remaining variables are arrays
    rmax,
    elevation_file,
    dune_file,
    background_erosion,
    sea_level_rise_rate=0.004,  # not an array
    sea_level_constant=True,  # not an array
    enable_shoreline_offset=False,
    shoreline_offset=[0],
    marsh_dynamics = False,
):

    # ###############################################################################
    # 9 - connect cascade domains (human management) with AST
    # ###############################################################################
    print('')
    print('Nt is ')
    print(nt)
    # --------- INITIALIZE ---------
    datadir = "B3D_Inputs/"
    cascade = Cascade(
        datadir,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file="Alongshore_Test-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=sea_level_rise_rate,
        sea_level_rise_constant=sea_level_constant,
        background_erosion=background_erosion,
        alongshore_section_count=alongshore_section_count,
        time_step_count=nt,
        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,
        num_cores=num_cores,
        roadway_management_module=False,  # no roadway management
        alongshore_transport_module=True,  # Is there brie coupling?
        beach_nourishment_module=False,  # no beach nourishment
        community_dynamics_module=False,  # no community dynamics
        enable_shoreline_offset=enable_shoreline_offset,  # Bool
        shoreline_offset=shoreline_offset,
        marsh_dynamics=marsh_dynamics,
    )
    # --------- LOOP ---------
    for time_step in range(nt - 1):

        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        cascade.update(Time_step = time_step)
        if cascade.b3d_break:
            break

    # --------- SAVE ---------
    save_directory = "Run_Output/"
    cascade.save(save_directory)

    return cascade

# ###############################################################################
# New test values
# ###############################################################################

# Specify variables to use in calling function
# Elevation file path name
e_file = "/B3D_Inputs/Marsh_Test_Inputs/InitElevHog.npy"
# Dune height path name
d_file = "/B3D_Inputs/Marsh_Test_Inputs/barrier3d-dunes.npy"
# Storm file path name
s_file = "/B3D_Inputs/Default_StormTimeSeries_1000yr.npy"
c_wd = os.getcwd()
nt_run = 10 # Number of years model will run
run_name = 'Smith_Test_Elev'


number_barrier3d_models = 5
rmin = [0.55] * number_barrier3d_models
rmax = [0.95] * number_barrier3d_models
elevation_file = [
    c_wd + "/B3D_Inputs/Smith/S_M_6.npy",
    c_wd + "/B3D_Inputs/Smith/S_M_4.npy",
    c_wd + "/B3D_Inputs/Smith/S_M_3.npy",
    c_wd + "/B3D_Inputs/Smith/S_M_2.npy",
    c_wd + "/B3D_Inputs/Smith/S_M_1.npy",
]
#dune_file = [c_wd+d_file] * number_barrier3d_models

dune_file = [
    c_wd + "/B3D_Inputs/Smith/S_M_D_6.npy",
    c_wd + "/B3D_Inputs/Smith/S_M_D_4.npy",
    c_wd + "/B3D_Inputs/Smith/S_M_D_3.npy",
    c_wd + "/B3D_Inputs/Smith/S_M_D_2.npy",
    c_wd + "/B3D_Inputs/Smith/S_M_D_1.npy",
             ]

storm_file = c_wd+s_file

Test_Alongshore_Runs(
    nt=nt_run,
    name=run_name,
    storm_file=c_wd+s_file,
    alongshore_section_count=number_barrier3d_models,
    num_cores=3,
    rmin=rmin,
    rmax=rmax,
    elevation_file=elevation_file,
    dune_file=dune_file,
    background_erosion=-1.00,
    sea_level_rise_rate=0.004,  # not an array
    sea_level_constant=True,  # not an array
    enable_shoreline_offset=False,
    marsh_dynamics=True,
)

def plot_ElevAnimation_CASCADE(
    cascade,
    directory,
    TMAX_MGMT,
    name,
    TMAX_SIM,
    ny=1,
    beach_management_ny=None,  # list of booleans the length of ny
    roadway_management_ny=None,
    y_lim=None,
    z_lim=3.5,
    fig_size=None,
    fig_eps=False,
):
    """
    NOTE THAT THE BEACH REPRESENTATION IS BASED ON A MODEL SPECIFIED BEACH WIDTH. We set the beach width for the
    remaining time steps after the community has been abandoned to the last managed beach width in order to not have a
    huge jump in the back-barrier position in Barrier3D. OTHERWISE, it is meaningless to the dynamics Barrier3D.
    """
    barrier3d = cascade.barrier3d

    # set up the domain; here we just use the first grid, but that could break in future runs
    BarrierLength = barrier3d[0].BarrierLength
    if np.any(beach_management_ny):
        indices = [i for i in range(ny) if beach_management_ny[i] == 1]
        iB3D = indices[0]
        print(1)
        MaxBeachWidth = (
            np.max(cascade.nourishments[iB3D].beach_width[0 : TMAX_MGMT[iB3D]]) / 10
        )  # dam
    else:
        MaxBeachWidth = cascade._initial_beach_width[0]
        print(2)
    OriginY = int(barrier3d[0].x_s_TS[0])
    AniDomainWidth = int(
        np.amax(barrier3d[0].InteriorWidth_AvgTS)
        + MaxBeachWidth
        + np.abs(barrier3d[0]._ShorelineChange)
        + OriginY
        + 35
    )

    os.chdir(directory)
    newpath = "Output/" + name + "/SimFrames/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    maxMGMT = np.max(TMAX_MGMT)

    if np.any(beach_management_ny) or np.any(roadway_management_ny):
        for t in range(maxMGMT + 1):

            # start with plotting t=0, then plot post-storm dune, interior, shoreface, etc. before management, treating this
            # as t=0.5 (i.e., this is really the final configuration from storms at t=1,2,3,4,5,... but before human
            # modifications); for the default dune inputs, we should be able to see natural dune growth at
            # t=0.5 before rebuild and storms (i.e., t=0.5 and t=1 should not be the same)
            if 0 < t <= TMAX_SIM:

                # post-storm variables in the BeachDuneManager are: interior, dune height, x_s, s_sf, beach width
                AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1

                for iB3D in range(ny):

                    # the nourishment scenario
                    if beach_management_ny[iB3D] and (t < TMAX_MGMT[iB3D] + 1):
                        actual_shoreline_post_storm = np.hstack(
                            [
                                0,
                                cascade.nourishments[iB3D]._post_storm_x_s[
                                    1 : TMAX_MGMT[iB3D] + 1
                                ],
                            ]
                        )
                        beach_width = (
                            cascade.nourishments[iB3D]._post_storm_beach_width[t] / 10
                        )
                        Domain = cascade.nourishments[iB3D]._post_storm_interior[t] * 10
                        Dunes = (
                            cascade.nourishments[iB3D]._post_storm_dunes[t]
                            + barrier3d[iB3D].BermEl
                        ) * 10

                    elif roadway_management_ny[iB3D] and (t < TMAX_MGMT[iB3D] + 1):
                        # the roadways scenario
                        actual_shoreline_post_storm = barrier3d[iB3D].x_s_TS[
                            0 : TMAX_MGMT[iB3D] + 1
                        ]
                        beach_width = cascade._initial_beach_width[iB3D] / 10
                        Domain = cascade.roadways[iB3D]._post_storm_interior[t] * 10
                        Dunes = (
                            cascade.roadways[iB3D]._post_storm_dunes[t]
                            + barrier3d[iB3D].BermEl
                        ) * 10

                    else:
                        # the natural scenario
                        actual_shoreline_post_storm = barrier3d[iB3D].x_s_TS[0:TMAX_SIM]
                        beach_width = cascade._initial_beach_width[iB3D] / 10
                        Domain = barrier3d[iB3D].DomainTS[t] * 10
                        Dunes = (
                            barrier3d[iB3D].DuneDomain[t, :, :] + barrier3d[iB3D].BermEl
                        ) * 10

                    # Build beach elevation domain, we only show beach width decreasing in increments of 10 m and we don't
                    # illustrate a berm, just a sloping beach up to the elevation of the berm
                    cellular_dune_toe_post_storm = np.floor(
                        actual_shoreline_post_storm[t] + beach_width
                    )
                    cellular_shoreline_post_storm = np.floor(
                        actual_shoreline_post_storm[t]
                    )
                    cellular_beach_width = int(
                        cellular_dune_toe_post_storm - cellular_shoreline_post_storm
                    )  # not actual bw

                    BeachDomain = np.zeros(
                        [
                            cellular_beach_width,
                            BarrierLength,
                        ]
                    )
                    if cellular_beach_width == 0:
                        pass
                    else:
                        add = (barrier3d[iB3D].BermEl - barrier3d[iB3D].SL) / (
                            cellular_beach_width + 1
                        )
                        for i in range(0, cellular_beach_width):
                            BeachDomain[i, :] = (barrier3d[iB3D].SL + add) * (i + 1)

                    # Make animation frame domain
                    Dunes = np.rot90(Dunes)
                    Dunes = np.flipud(Dunes)
                    Beach = BeachDomain * 10
                    Domain = np.vstack([Beach, Dunes, Domain])
                    Domain[Domain < 0] = -1
                    widthTS = len(Domain)
                    OriginTstart = int(cellular_shoreline_post_storm)
                    OriginTstop = OriginTstart + widthTS
                    xOrigin = iB3D * BarrierLength
                    AnimateDomain[
                        OriginTstart:OriginTstop, xOrigin : xOrigin + BarrierLength
                    ] = Domain

                # Plot and save
                if fig_size is not None:
                    elevFig1 = plt.figure(figsize=fig_size)
                else:
                    elevFig1 = plt.figure(figsize=(7, 7))
                ax = elevFig1.add_subplot(111)
                cax = ax.pcolormesh(
                    AnimateDomain,
                    cmap="terrain",
                    vmin=-1.1,
                    vmax=z_lim,
                    # edgecolors="w",  # for debugging
                    # linewidth=0.01,
                )
                cbar = elevFig1.colorbar(cax)
                cbar.set_label("elevation (m MHW)", rotation=270)
                plt.xlabel("alongshore distance (dam)")
                plt.ylabel("cross-shore distance (dam)")
                timestr = (
                    "Time = " + str(t - 0.5) + " yrs"
                )  # we are letting the post-storm output represent 0.5 years
                if y_lim is not None:
                    plt.ylim(y_lim)
                    plt.text(3, y_lim[0] + 3, timestr, color="w")
                else:
                    plt.ylim(bottom=OriginY - 35)
                    plt.text(1, OriginY - 33, timestr)
                plt.tight_layout()
                plt.rcParams.update({"font.size": 11})
                # elevFig1.tight_layout()
                if fig_eps:
                    name = "elev_" + str(t - 1) + "pt5.eps"
                    elevFig1.savefig(name, format="eps")
                else:
                    name = "elev_" + str(t - 1) + "pt5"
                    elevFig1.savefig(name)  # dpi=200
                plt.close(elevFig1)

    for t in range(TMAX_SIM):

        # ok, now the annual time step, which incorporates human modifications to the shoreface, beach, dune, & interior
        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1

        for iB3D in range(ny):

            actual_shoreline_post_humans = barrier3d[iB3D].x_s_TS[0 : TMAX_SIM + 1]

            if beach_management_ny[iB3D]:
                beach_width = cascade.nourishments[iB3D].beach_width[t] / 10
            else:
                # both roadways and natural scenario
                beach_width = cascade._initial_beach_width[iB3D] / 10

            if beach_management_ny[iB3D] and np.isnan(beach_width):
                beach_width = (
                    cascade.nourishments[iB3D].beach_width[t - 1] / 10
                )  # this can happen if the barrier height drowns

            cellular_dune_toe_post_humans = np.floor(
                actual_shoreline_post_humans[t] + beach_width
            )

            cellular_shoreline_post_humans = np.floor(actual_shoreline_post_humans[t])
            cellular_beach_width = int(
                cellular_dune_toe_post_humans - cellular_shoreline_post_humans
            )  # not actual bw

            if t == maxMGMT + 1:
                cellular_beach_width_final = (
                    cellular_beach_width  # save this for plotting moving forward
                )
            if t > maxMGMT + 1:
                cellular_beach_width = cellular_beach_width_final

            BeachDomain = np.zeros(
                [
                    cellular_beach_width,
                    BarrierLength,
                ]
            )

            if cellular_beach_width == 0:
                pass
            else:
                add = (barrier3d[iB3D].BermEl - barrier3d[iB3D].SL) / (
                    cellular_beach_width + 1
                )
                for i in range(0, cellular_beach_width):
                    BeachDomain[i, :] = (barrier3d[iB3D].SL + add) * (i + 1)

            # Make animation frame domain
            Domain = barrier3d[iB3D].DomainTS[t] * 10
            Dunes = (barrier3d[iB3D].DuneDomain[t, :, :] + barrier3d[iB3D].BermEl) * 10
            Dunes = np.rot90(Dunes)
            Dunes = np.flipud(Dunes)
            Beach = BeachDomain * 10
            Domain = np.vstack([Beach, Dunes, Domain])
            Domain[Domain < -3] = -3
            widthTS = len(Domain)
            OriginTstart = int(cellular_shoreline_post_humans)
            OriginTstop = OriginTstart + widthTS
            xOrigin = iB3D * BarrierLength
            AnimateDomain[
                OriginTstart:OriginTstop, xOrigin : xOrigin + BarrierLength
            ] = Domain

        # Plot and save
        if fig_size is not None:
            elevFig2 = plt.figure(figsize=fig_size)
        else:
            elevFig2 = plt.figure(figsize=(7, 7))
        ax = elevFig2.add_subplot(111)
        cax = ax.pcolormesh(
            AnimateDomain,
            cmap="terrain",
            vmin=-1.1,
            vmax=z_lim,
            # edgecolors="w",  # for debugging
            # linewidth=0.01,
        )
        cbar = elevFig2.colorbar(cax)
        cbar.set_label("elevation (m MHW)", rotation=270)
        plt.xlabel("alongshore distance (dam)")
        plt.ylabel("cross-shore distance (dam)")
        timestr = "Time = " + str(t) + " yrs"
        if y_lim is not None:
            plt.ylim(y_lim)
            plt.text(3, y_lim[0] + 3, timestr, color="w")
        else:
            plt.ylim(bottom=OriginY - 35)
            plt.text(1, OriginY - 33, timestr)
        plt.tight_layout()
        plt.rcParams.update({"font.size": 11})
        # elevFig2.tight_layout()
        if fig_eps:
            name = "elev_" + str(t) + ".eps"
            elevFig2.savefig(name, format="eps")
        else:
            name = "elev_" + str(t)
            elevFig2.savefig(name)  # dpi=200
        plt.close(elevFig2)

    frames = []

    for filenum in range(TMAX_SIM):
        if fig_eps:
            filename = "elev_" + str(filenum) + ".eps"
        else:
            filename = "elev_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))

        # this might break on AST runs: need to test it (don't know if I need the first part
        for iB3D in range(ny):
            if (beach_management_ny[iB3D] or roadway_management_ny[iB3D]) and (
                filenum
                < TMAX_MGMT[iB3D]
                # and cascade.nourishments[iB3D]._post_storm_interior[TMAX_MGMT] is not None
            ):
                if fig_eps:
                    filename = "elev_" + str(filenum) + "pt5" ".eps"
                else:
                    filename = "elev_" + str(filenum) + "pt5" ".png"
                frames.append(imageio.imread(filename))

    imageio.mimsave("elev.gif", frames, "GIF-FI")
    print()
    print("[ * GIF successfully generated * ]")

os.chdir("/Users/ceclmac/PycharmProjects/CASCADE/Run_output")
name_prefix = run_name

# --------- plot ---------
output = np.load(run_name + ".npz", allow_pickle=True)
cascade = output["cascade"]
cascade = cascade[0]
b3d = cascade.barrier3d
ny = np.size(b3d)

directory = "/Users/ceclmac/PycharmProjects/CASCADE/"
# TMax_MGMT = Needed 0
# TMAX_Sim = Last simulation year of the model 99
TMax_Sim = nt_run  # Give length of simulation
TMax_MGMT = [0] * ny
beach_management_ny = [False] * ny
roadway_management_ny = [False] * ny

plot_ElevAnimation_CASCADE(
    cascade,
    ny=number_barrier3d_models,
    directory=directory,
    TMAX_MGMT=TMax_MGMT,  # an array
    name=run_name,
    TMAX_SIM=TMax_Sim,  # not an array
    beach_management_ny=beach_management_ny,
    roadway_management_ny=roadway_management_ny,
    y_lim=None,
    z_lim=None,
    fig_size=None,
    fig_eps=False,
)

