import numpy as np
import time
import matplotlib.pyplot as plt
import os
import imageio
import copy


def plot_ElevAnimation_CASCADE(
    cascade,
    directory,
    TMAX_MGMT,
    name,
    TMAX_SIM,
    Model_Grids_Of_Interest,
    ny=1,
    beach_management_ny=None,  # list of booleans the length of ny
    roadway_management_ny=None,
    y_lim=None,                # if not None, overrides auto-centering
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

    Focused_B3D = []
    for focus in Model_Grids_Of_Interest:
        Focused_B3D.append(copy.deepcopy(barrier3d[focus]))
    barrier3d = copy.deepcopy(Focused_B3D)

    # set up the domain; here we just use the first grid
    BarrierLength = barrier3d[0].BarrierLength

    if np.any(beach_management_ny):
        indices = [i for i in range(ny) if beach_management_ny[i] == 1]
        iB3D = indices[0]
        MaxBeachWidth = (
            np.max(cascade.nourishments[iB3D].beach_width[0: TMAX_MGMT[iB3D]]) / 10
        )  # dam
    else:
        MaxBeachWidth = cascade._initial_beach_width[0]

    # shoreline index (cell index of shoreline position)
    OriginY = int(barrier3d[0].x_s_TS[0])

    # maximum interior width across all grids (cells, cross-shore)
    max_interior = int(
        np.max([np.max(b.InteriorWidth_AvgTS) for b in barrier3d])
    )

    # Size of full plotting canvas (same as before)
    AniDomainWidth = int(
        np.amax(barrier3d[0].InteriorWidth_AvgTS)
        + MaxBeachWidth
        + np.abs(barrier3d[0]._ShorelineChange)
        + OriginY
        + 200  # padding
    )

    # --- NEW: default y-limits to center island better ---
    #  - 35 cells seaward of shoreline
    #  - full island interior (max_interior) plus 50-cell buffer landward
    default_ymin = OriginY - 35
    default_ymax = OriginY + max_interior + 50

    os.chdir(directory)
    newpath = "gif_output/" + name + "/SimFrames/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    maxMGMT = np.max(TMAX_MGMT)

    # ============================
    # POST-STORM (t - 0.5 yrs) PLOTS
    # ============================
    if np.any(beach_management_ny) or np.any(roadway_management_ny):
        for t in range(maxMGMT + 1):

            if 0 < t <= TMAX_SIM:

                AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1

                for iB3D in range(ny):

                    # nourishment scenario
                    if beach_management_ny[iB3D] and (t < TMAX_MGMT[iB3D] + 1):
                        actual_shoreline_post_storm = np.hstack(
                            [
                                0,
                                cascade.nourishments[iB3D]._post_storm_x_s[
                                    1: TMAX_MGMT[iB3D] + 1
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

                    # roadway scenario
                    elif roadway_management_ny[iB3D] and (t < TMAX_MGMT[iB3D] + 1):
                        actual_shoreline_post_storm = barrier3d[iB3D].x_s_TS[
                            0: TMAX_MGMT[iB3D] + 1
                        ]
                        beach_width = cascade._initial_beach_width[iB3D] / 10
                        Domain = cascade.roadways[iB3D]._post_storm_interior[t] * 10
                        Dunes = (
                            cascade.roadways[iB3D]._post_storm_dunes[t]
                            + barrier3d[iB3D].BermEl
                        ) * 10

                    # natural scenario
                    else:
                        actual_shoreline_post_storm = barrier3d[iB3D].x_s_TS[0:TMAX_SIM]
                        beach_width = cascade._initial_beach_width[iB3D] / 10
                        Domain = barrier3d[iB3D].DomainTS[t] * 10
                        Dunes = (
                            barrier3d[iB3D].DuneDomain[t, :, :]
                            + barrier3d[iB3D].BermEl
                        ) * 10

                    # Build beach elevation domain
                    cellular_dune_toe_post_storm = np.floor(
                        actual_shoreline_post_storm[t] + beach_width
                    )
                    cellular_shoreline_post_storm = np.floor(
                        actual_shoreline_post_storm[t]
                    )
                    cellular_beach_width = int(
                        cellular_dune_toe_post_storm - cellular_shoreline_post_storm
                    )

                    BeachDomain = np.zeros(
                        [
                            cellular_beach_width,
                            BarrierLength,
                        ]
                    )
                    if cellular_beach_width > 0:
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
                        OriginTstart:OriginTstop, xOrigin: xOrigin + BarrierLength
                    ] = Domain

                # Plot and save
                if fig_size is not None:
                    elevFig1 = plt.figure(figsize=fig_size)
                else:
                    elevFig1 = plt.figure(figsize=(21, 7))
                ax = elevFig1.add_subplot(111)
                cax = ax.pcolormesh(
                    AnimateDomain,
                    cmap="terrain",
                    vmin=-1.1,
                    vmax=z_lim,
                )
                cbar = elevFig1.colorbar(cax)
                cbar.set_label("elevation (m MHW)", rotation=270)
                plt.xlabel("alongshore distance (dam)")
                plt.ylabel("cross-shore distance (dam)")
                timestr = "Time = " + str(t - 0.5) + " yrs"

                # ------- NEW: centered y-limits -------
                if y_lim is not None:
                    plt.ylim(y_lim)
                    text_y = y_lim[0] + 3
                else:
                    plt.ylim(default_ymin, default_ymax)
                    text_y = default_ymin + 3
                plt.text(3, text_y, timestr, color="w")
                # ---------------------------------------

                plt.tight_layout()
                plt.rcParams.update({"font.size": 11})
                if fig_eps:
                    fname = "elev_" + str(t - 1) + "pt5.eps"
                    elevFig1.savefig(fname, format="eps")
                else:
                    fname = "elev_" + str(t - 1) + "pt5"
                    elevFig1.savefig(fname)
                plt.close(elevFig1)

    # ============================
    # ANNUAL (t yrs) PLOTS
    # ============================
    for t in range(TMAX_SIM):

        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1

        for iB3D in range(ny):

            actual_shoreline_post_humans = barrier3d[iB3D].x_s_TS[0: TMAX_SIM + 1]

            if beach_management_ny[iB3D]:
                beach_width = cascade.nourishments[iB3D].beach_width[t] / 10
            else:
                beach_width = cascade._initial_beach_width[iB3D] / 10

            if beach_management_ny[iB3D] and np.isnan(beach_width):
                beach_width = cascade.nourishments[iB3D].beach_width[t - 1] / 10

            cellular_dune_toe_post_humans = np.floor(
                actual_shoreline_post_humans[t] + beach_width
            )
            cellular_shoreline_post_humans = np.floor(
                actual_shoreline_post_humans[t]
            )
            cellular_beach_width = int(
                cellular_dune_toe_post_humans - cellular_shoreline_post_humans
            )

            if t == maxMGMT + 1:
                cellular_beach_width_final = cellular_beach_width
            if t > maxMGMT + 1:
                cellular_beach_width = cellular_beach_width_final

            BeachDomain = np.zeros(
                [
                    cellular_beach_width,
                    BarrierLength,
                ]
            )

            if cellular_beach_width > 0:
                add = (barrier3d[iB3D].BermEl - barrier3d[iB3D].SL) / (
                    cellular_beach_width + 1
                )
                for i in range(0, cellular_beach_width):
                    BeachDomain[i, :] = (barrier3d[iB3D].SL + add) * (i + 1)

            Domain = barrier3d[iB3D].DomainTS[t] * 10
            Dunes = (
                barrier3d[iB3D].DuneDomain[t, :, :] + barrier3d[iB3D].BermEl
            ) * 10
            Dunes = np.rot90(Dunes)
            Dunes = np.flipud(Dunes)
            Beach = BeachDomain * 10
            Domain = np.vstack([Beach, Dunes, Domain])
            Domain = np.fliplr(Domain)
            Domain[Domain < -3] = -3
            Domain[Domain > 6] = 6
            widthTS = len(Domain)
            OriginTstart = int(cellular_shoreline_post_humans)
            OriginTstop = OriginTstart + widthTS
            xOrigin = iB3D * BarrierLength

            AnimateDomain[
                OriginTstart:OriginTstop, xOrigin: xOrigin + BarrierLength
            ] = Domain

        # Plot and save
        if fig_size is not None:
            elevFig2 = plt.figure(figsize=fig_size)
        else:
            elevFig2 = plt.figure(figsize=(21, 7))
        ax = elevFig2.add_subplot(111)
        cax = ax.pcolormesh(
            AnimateDomain,
            cmap="terrain",
            vmin=-1.1,
            vmax=z_lim,
        )
        cbar = elevFig2.colorbar(cax, pad=0.01)
        cbar.set_label("elevation (m MHW)", rotation=270)

        plt.xlabel("alongshore distance (km)")
        plt.ylabel("cross-shore distance (km)")
        timestr = "Time = " + str(t) + " yrs"

        # ------- NEW: centered y-limits again -------
        if y_lim is not None:
            plt.ylim(y_lim)
            text_y = y_lim[0] + 3
        else:
            plt.ylim(default_ymin, default_ymax)
            text_y = default_ymin + 3
        plt.text(3, text_y, timestr, color="w")
        # --------------------------------------------

        # Convert decameter labels to kilometers (unchanged)
        xticks = plt.xticks()
        yticks = plt.yticks()

        xtick_loc = copy.deepcopy(xticks[0])
        ytick_loc = copy.deepcopy(yticks[0])

        xmin = 0
        xmax = len(AnimateDomain[0])

        focused_x_tick_location = []
        for k in range(len(xtick_loc)):
            if xmin <= xtick_loc[k] <= xmax:
                focused_x_tick_location.append(copy.deepcopy(xtick_loc[k]))

        New_X_Labels = np.divide(copy.deepcopy(focused_x_tick_location), 100)
        plt.xticks(focused_x_tick_location, labels=New_X_Labels)

        Y_vals = plt.ylim()
        Y_min = min(Y_vals)
        Y_max = max(Y_vals)

        focused_y_tick_location = []
        for l in range(len(ytick_loc)):
            if Y_min <= ytick_loc[l] <= Y_max:
                focused_y_tick_location.append(copy.deepcopy(ytick_loc[l]))

        New_Y_Labels = np.divide(copy.deepcopy(focused_y_tick_location), 100)
        plt.yticks(focused_y_tick_location, labels=New_Y_Labels)

        plt.tight_layout()
        plt.subplots_adjust(right=1.12)
        plt.rcParams.update({"font.size": 11})

        if fig_eps:
            fname = "elev_" + str(t) + ".eps"
            elevFig2.savefig(fname, format="eps")
        else:
            fname = "elev_" + str(t)
            elevFig2.savefig(fname)
        plt.close(elevFig2)

    # Build GIF (unchanged)
    frames = []
    for filenum in range(TMAX_SIM):
        if fig_eps:
            filename = "elev_" + str(filenum) + ".eps"
        else:
            filename = "elev_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))

        for iB3D in range(ny):
            if (beach_management_ny[iB3D] or roadway_management_ny[iB3D]) and (
                filenum < TMAX_MGMT[iB3D]
            ):
                if fig_eps:
                    filename = "elev_" + str(filenum) + "pt5.eps"
                else:
                    filename = "elev_" + str(filenum) + "pt5.png"
                frames.append(imageio.imread(filename))

    imageio.mimsave("elev.gif", frames, "GIF", fps=2)
    print()
    print("[ * GIF successfully generated * ]")


# ----------------- DRIVER CODE (same as yours) -----------------

os.chdir(r'/output/raw_runs/HAT_1978_1997_Natural_State')
run_name = "HAT_1978_1997_Natural_State"
name_prefix = run_name
nt_run = 20  # Number of years to plot
number_barrier3d_models = 105  # Number of B3D domains to plot
Model_Grids_Of_Interest = range(15, 120)  # Specific B3D domains to plot

output = np.load(run_name + ".npz", allow_pickle=True)
cascade = output["cascade"][0]
b3d = cascade.barrier3d
ny = np.size(b3d)

directory = r"C:\Users\hanna\PycharmProjects\CASCADE\output"
TMax_Sim = nt_run
TMax_MGMT = [0] * ny
beach_management_ny = [False] * ny
roadway_management_ny = [False] * ny

plot_ElevAnimation_CASCADE(
    cascade,
    ny=number_barrier3d_models,
    directory=directory,
    TMAX_MGMT=TMax_MGMT,
    name=run_name,
    TMAX_SIM=TMax_Sim,
    beach_management_ny=beach_management_ny,
    roadway_management_ny=roadway_management_ny,
    y_lim=None,          # or pass a tuple like (200, 600) if you want explicit limits
    z_lim=3.5,
    fig_size=None,
    fig_eps=False,
    Model_Grids_Of_Interest=Model_Grids_Of_Interest,
)
