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
        beach_management_ny=None,
        roadway_management_ny=None,
        y_lim=None,
        z_lim=3.5,
        fig_size=None,
        fig_eps=False,
        show_domain_labels=True,   # Control domain labeling
        domain_label_interval=10,  # Label every N domains
        first_domain_index=0,      # Index of first plotted domain in full cascade
):
    """
    Plot CASCADE elevation animations with proper scaling and domain labels.

    Key changes:
    - Ensures full island is visible (no cutoff)
    - Adds domain number labels at specified intervals
    - Allows plotting a subset of domains (e.g., hide buffers)
    """
    barrier3d_full = cascade.barrier3d

    # Subset barrier3d to the domains of interest
    Focused_B3D = []
    for focus in Model_Grids_Of_Interest:
        Focused_B3D.append(copy.deepcopy(barrier3d_full[focus]))
    barrier3d = copy.deepcopy(Focused_B3D)

    # Make ny consistent with what we're actually plotting
    ny = len(barrier3d)

    # Set up the domain
    BarrierLength = barrier3d[0].BarrierLength

    if np.any(beach_management_ny):
        indices = [i for i in range(ny) if beach_management_ny[i] == 1]
        iB3D = indices[0]
        MaxBeachWidth = (
                np.max(cascade.nourishments[iB3D].beach_width[0: TMAX_MGMT[iB3D]]) / 10
        )
    else:
        MaxBeachWidth = cascade._initial_beach_width[0]

    OriginY = int(barrier3d[0].x_s_TS[0])

    # Calculate domain width with extra padding to prevent cutoff
    AniDomainWidth = int(
        np.amax(barrier3d[0].InteriorWidth_AvgTS)
        + MaxBeachWidth
        + np.abs(barrier3d[0]._ShorelineChange)
        + OriginY
        + 250  # Increased padding to prevent cutoff
    )

    os.chdir(directory)
    newpath = "gif_output/" + name + "/SimFrames/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    maxMGMT = np.max(TMAX_MGMT)

    # -------------------------------------------------------------------------
    # Post-storm frames (if any management)
    # -------------------------------------------------------------------------
    if np.any(beach_management_ny) or np.any(roadway_management_ny):
        for t in range(maxMGMT + 1):
            if 0 < t <= TMAX_SIM:
                AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1

                for iB3D in range(ny):
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

                    else:
                        actual_shoreline_post_storm = barrier3d[iB3D].x_s_TS[0:TMAX_SIM]
                        beach_width = cascade._initial_beach_width[iB3D] / 10
                        Domain = barrier3d[iB3D].DomainTS[t] * 10
                        Dunes = (
                                        barrier3d[iB3D].DuneDomain[t, :, :]
                                        + barrier3d[iB3D].BermEl
                                ) * 10

                    cellular_dune_toe_post_storm = np.floor(
                        actual_shoreline_post_storm[t] + beach_width
                    )
                    cellular_shoreline_post_storm = np.floor(
                        actual_shoreline_post_storm[t]
                    )
                    cellular_beach_width = int(
                        cellular_dune_toe_post_storm - cellular_shoreline_post_storm
                    )

                    BeachDomain = np.zeros([cellular_beach_width, BarrierLength])
                    if cellular_beach_width == 0:
                        pass
                    else:
                        add = (barrier3d[iB3D].BermEl - barrier3d[iB3D].SL) / (
                                cellular_beach_width + 1
                        )
                        for i in range(0, cellular_beach_width):
                            BeachDomain[i, :] = (barrier3d[iB3D].SL + add) * (i + 1)

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

                # Plot post-storm frame
                _plot_frame(
                    AnimateDomain, t - 0.5, OriginY, z_lim, y_lim,
                    fig_size, fig_eps, name, BarrierLength, ny,
                    show_domain_labels, domain_label_interval,
                    is_post_storm=True,
                    first_domain_index=first_domain_index,
                )

    # -------------------------------------------------------------------------
    # Annual time step frames
    # -------------------------------------------------------------------------
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

            BeachDomain = np.zeros([cellular_beach_width, BarrierLength])

            if cellular_beach_width == 0:
                pass
            else:
                add = (barrier3d[iB3D].BermEl - barrier3d[iB3D].SL) / (
                        cellular_beach_width + 1
                )
                for i in range(0, cellular_beach_width):
                    BeachDomain[i, :] = (barrier3d[iB3D].SL + add) * (i + 1)

            Domain = barrier3d[iB3D].DomainTS[t] * 10
            Dunes = (
                            barrier3d[iB3D].DuneDomain[t, :, :]
                            + barrier3d[iB3D].BermEl
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

        # Plot annual frame
        _plot_frame(
            AnimateDomain, t, OriginY, z_lim, y_lim,
            fig_size, fig_eps, name, BarrierLength, ny,
            show_domain_labels, domain_label_interval,
            is_post_storm=False,
            first_domain_index=first_domain_index,
        )

    # -------------------------------------------------------------------------
    # Create GIF
    # -------------------------------------------------------------------------
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
                    filename = "elev_" + str(filenum) + "pt5" ".eps"
                else:
                    filename = "elev_" + str(filenum) + "pt5" ".png"
                frames.append(imageio.imread(filename))

    imageio.mimsave("elev.gif", frames, "GIF", fps=2)
    print()
    print("[ * GIF successfully generated * ]")


def _plot_frame(AnimateDomain, time_val, OriginY, z_lim, y_lim,
                fig_size, fig_eps, run_name, BarrierLength, ny,
                show_domain_labels, domain_label_interval,
                is_post_storm=False,
                first_domain_index=0):
    """
    Helper function to plot and save individual frames with proper formatting.
    """
    if fig_size is not None:
        elevFig = plt.figure(figsize=fig_size)
    else:
        elevFig = plt.figure(figsize=(21, 7))

    ax = elevFig.add_subplot(111)
    cax = ax.pcolormesh(
        AnimateDomain,
        cmap="terrain",
        vmin=-1.1,
        vmax=z_lim,
    )
    cbar = elevFig.colorbar(cax, pad=0.01)
    cbar.set_label("elevation (m MHW)", rotation=270, labelpad=20)

    plt.xlabel("alongshore distance (km)")
    plt.ylabel("cross-shore distance (km)")
    timestr = "Time = " + str(time_val) + " yrs"

    # Set y-limits to show full island (prevent cutoff)
    if y_lim is not None:
        plt.ylim(y_lim)
        text_y_pos = y_lim[0] + 3
    else:
        # Calculate appropriate y-limits to show all data
        data_extent = np.where(AnimateDomain > -1)
        if len(data_extent[0]) > 0:
            y_min = max(0, np.min(data_extent[0]) - 50)  # Add padding
            y_max = np.max(data_extent[0]) + 50
            plt.ylim(y_min, y_max)
            text_y_pos = y_min + 3
        else:
            plt.ylim(bottom=OriginY - 50)
            text_y_pos = OriginY - 48

    plt.text(3, text_y_pos, timestr, color='w', fontsize=12, weight='bold')

    # Convert decameter labels to kilometers
    xticks = plt.xticks()
    yticks = plt.yticks()

    xtick_loc = copy.deepcopy(xticks[0])
    ytick_loc = copy.deepcopy(yticks[0])

    # X-axis (alongshore) conversion
    xmin = 0
    xmax = len(AnimateDomain[0])

    focused_x_tick_location = []
    for k in range(len(xtick_loc)):
        if xmin <= xtick_loc[k] <= xmax:
            focused_x_tick_location.append(copy.deepcopy(xtick_loc[k]))

    New_X_Labels = np.divide(copy.deepcopy(focused_x_tick_location), 100)
    plt.xticks(focused_x_tick_location, labels=New_X_Labels)

    # Y-axis (cross-shore) conversion
    Y_vals = plt.ylim()
    Y_min = min(Y_vals)
    Y_max = max(Y_vals)

    focused_y_tick_location = []
    for l in range(len(ytick_loc)):
        if Y_min <= ytick_loc[l] <= Y_max:
            focused_y_tick_location.append(copy.deepcopy(ytick_loc[l]))

    New_Y_Labels = np.divide(copy.deepcopy(focused_y_tick_location), 100)
    plt.yticks(focused_y_tick_location, labels=New_Y_Labels)

    # Add domain number labels
    if show_domain_labels:
        # Add secondary x-axis for domain numbers
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())

        # Calculate domain label positions
        # Assuming domains start at x=0 and each domain is BarrierLength wide
        domain_positions = []
        domain_numbers = []

        for i in range(0, ny, domain_label_interval):
            # Position at center of domain
            domain_center = i * BarrierLength + BarrierLength / 2
            if domain_center < xmax:
                domain_positions.append(domain_center)
                # Account for which domain index in the full cascade we're starting from
                actual_domain_num = first_domain_index + i + 1
                domain_numbers.append(f'D{actual_domain_num}')

        ax2.set_xticks(domain_positions)
        ax2.set_xticklabels(domain_numbers, fontsize=9)
        ax2.set_xlabel("Domain Number", fontsize=11)

    plt.tight_layout()
    plt.subplots_adjust(right=1.12)
    plt.rcParams.update({"font.size": 11})

    # Save figure
    if is_post_storm:
        if fig_eps:
            filename = "elev_" + str(int(time_val)) + "pt5.eps"
            elevFig.savefig(filename, format="eps")
        else:
            filename = "elev_" + str(int(time_val)) + "pt5"
            elevFig.savefig(filename)
    else:
        if fig_eps:
            filename = "elev_" + str(int(time_val)) + ".eps"
            elevFig.savefig(filename, format="eps")
        else:
            filename = "elev_" + str(int(time_val))
            elevFig.savefig(filename)

    plt.close(elevFig)


# =============================================================================
# MAIN EXECUTION
# =============================================================================

os.chdir(r'C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_1978_1997_mod_storms_be1')
run_name = "HAT_1978_1997_mod_storms_be1"
name_prefix = run_name
nt_run = 19  # Number of years to plot

# Load CASCADE output
output = np.load(run_name + ".npz", allow_pickle=True)
cascade = output["cascade"][0]
b3d = cascade.barrier3d
total_domains = len(b3d)
print("Total domains in cascade:", total_domains)

# === CONTROL THIS FLAG TO SHOW/HIDE BUFFERS ===
include_buffers = False          # set False to hide buffers
n_buffer_each_side = 15         # adjust if different

if include_buffers:
    # Plot ALL domains (real + buffers)
    Model_Grids_Of_Interest = range(26, 33) #range(total_domains)
    first_domain_index = 0
else:
    # Plot ONLY real domains (no buffers)
    start_real = n_buffer_each_side
    end_real = total_domains - n_buffer_each_side   # exclusive
    Model_Grids_Of_Interest = range(start_real, end_real)
    first_domain_index = start_real                 # for domain labels

ny = len(Model_Grids_Of_Interest)  # number of domains actually plotted

directory = r"C:\Users\hanna\PycharmProjects\CASCADE\output"
TMax_Sim = nt_run
TMax_MGMT = [0] * ny
beach_management_ny = [False] * ny
roadway_management_ny = [False] * ny

plot_ElevAnimation_CASCADE(
    cascade,
    directory=directory,
    TMAX_MGMT=TMax_MGMT,
    name=run_name,
    TMAX_SIM=TMax_Sim,
    Model_Grids_Of_Interest=Model_Grids_Of_Interest,
    ny=ny,
    beach_management_ny=beach_management_ny,
    roadway_management_ny=roadway_management_ny,
    y_lim=None,
    z_lim=3.5,
    fig_size=None,
    fig_eps=False,
    show_domain_labels=True,
    domain_label_interval=10,
    first_domain_index=first_domain_index,
)
