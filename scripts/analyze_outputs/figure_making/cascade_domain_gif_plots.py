import numpy as np
import time
import matplotlib.pyplot as plt
import os
import imageio.v2 as imageio  # quiets deprecation warnings
import copy
import matplotlib.patheffects as pe


def plot_ElevAnimation_CASCADE(
    cascade,
    directory,
    TMAX_MGMT,
    name,
    TMAX_SIM,
    Model_Grids_Of_Interest,
    ny=None,
    beach_management_ny=None,  # list of booleans length == ny
    roadway_management_ny=None,
    y_lim=None,
    z_lim=3.5,
    fig_size=None,
    fig_eps=False,
    *,
    # Alongshore orientation helpers
    reverse_alongshore=False,  # reverse the order of domains left↔right
    flip_each_domain=False,    # mirror each domain's columns left↔right
    draw_domain_seams=True,    # show faint vertical lines between domains
    draw_domain_labels=True,   # <-- NEW: overlay domain numbers on the figure
):
    """
    Render elevation frames (and optional post-storm frames) for a subset of Barrier3D domains,
    then stitch them into a GIF.

    - Distances are plotted in decameters (dam) internally; ticks relabeled to km.
    - Elevations are meters MHW.
    - For natural-only runs, keep both management lists all False and TMAX_MGMT all zeros.
    """
    # --- Normalize inputs ---
    if ny is None:
        ny = len(Model_Grids_Of_Interest)
    if beach_management_ny is None:
        beach_management_ny = [False] * ny
    if roadway_management_ny is None:
        roadway_management_ny = [False] * ny
    assert len(beach_management_ny) == ny, "beach_management_ny length must equal ny"
    assert len(roadway_management_ny) == ny, "roadway_management_ny length must equal ny"
    assert len(TMAX_MGMT) == ny, "TMAX_MGMT length must equal ny"

    barrier3d_full = cascade.barrier3d

    # Focus only on requested alongshore grids (deep copies to avoid side-effects)
    Focused_B3D = [copy.deepcopy(barrier3d_full[idx]) for idx in Model_Grids_Of_Interest]
    # Labels should reflect original ids in displayed order
    display_order_labels = list(Model_Grids_Of_Interest)
    if reverse_alongshore:
        Focused_B3D = Focused_B3D[::-1]
        display_order_labels = display_order_labels[::-1]
    barrier3d = copy.deepcopy(Focused_B3D)

    # Domain setup (use first focused grid for canvas sizing)
    BarrierLength = barrier3d[0].BarrierLength

    # Max beach width for canvas sizing (mgmt vs natural)
    if np.any(beach_management_ny):
        indices = [i for i in range(ny) if beach_management_ny[i]]
        iB3D = indices[0]
        MaxBeachWidth = (np.max(cascade.nourishments[iB3D].beach_width[0:TMAX_MGMT[iB3D]]) / 10.0)  # dam
    else:
        MaxBeachWidth = cascade._initial_beach_width[0]  # dam
    OriginY = int(barrier3d[0].x_s_TS[0])
    AniDomainWidth = int(
        np.amax(barrier3d[0].InteriorWidth_AvgTS)
        + MaxBeachWidth
        + np.abs(barrier3d[0]._ShorelineChange)
        + OriginY
        + 50  # padding
    )

    # Helper to draw seams and labels (call after pcolormesh + y-limits are set)
    def _decorate_axes(ax):
        if draw_domain_seams:
            for k in range(1, ny):
                ax.axvline(x=k * BarrierLength, color=(1, 1, 1, 0.25), linewidth=0.8)
        if draw_domain_labels:
            # place labels near the top of the current y-range
            y_top = ax.get_ylim()[1]
            y_text = y_top - 8  # dam from top; tweak if needed
            for k in range(ny):
                x_center = k * BarrierLength + BarrierLength / 2.0
                label = str(display_order_labels[k])
                ax.text(
                    x_center, y_text, label,
                    ha="center", va="top", color="white", fontsize=9,
                    path_effects=[pe.withStroke(linewidth=2, foreground="black", alpha=0.6)]
                )

    # Output folder
    cwd_before = os.getcwd()
    try:
        os.chdir(directory)
        newpath = os.path.join("Output", name, "SimFrames")
        os.makedirs(newpath, exist_ok=True)
        os.chdir(newpath)

        maxMGMT = int(np.max(TMAX_MGMT)) if len(TMAX_MGMT) > 0 else 0

        # ---------- (t - 0.5) post-storm frames (only if any management flags are True) ----------
        if np.any(beach_management_ny) or np.any(roadway_management_ny):
            for t in range(maxMGMT + 1):
                if 0 < t <= TMAX_SIM:
                    AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1
                    for iB3D in range(ny):
                        # nourishment
                        if beach_management_ny[iB3D] and (t < TMAX_MGMT[iB3D] + 1):
                            actual_shoreline_post_storm = np.hstack(
                                [0, cascade.nourishments[iB3D]._post_storm_x_s[1 : TMAX_MGMT[iB3D] + 1]]
                            )
                            beach_width = cascade.nourishments[iB3D]._post_storm_beach_width[t] / 10.0
                            Domain = cascade.nourishments[iB3D]._post_storm_interior[t] * 10.0
                            Dunes = (cascade.nourishments[iB3D]._post_storm_dunes[t] + barrier3d[iB3D].BermEl) * 10.0
                        # roadway
                        elif roadway_management_ny[iB3D] and (t < TMAX_MGMT[iB3D] + 1):
                            actual_shoreline_post_storm = barrier3d[iB3D].x_s_TS[0 : TMAX_MGMT[iB3D] + 1]
                            beach_width = cascade._initial_beach_width[iB3D] / 10.0
                            Domain = cascade.roadways[iB3D]._post_storm_interior[t] * 10.0
                            Dunes = (cascade.roadways[iB3D]._post_storm_dunes[t] + barrier3d[iB3D].BermEl) * 10.0
                        # natural
                        else:
                            actual_shoreline_post_storm = barrier3d[iB3D].x_s_TS[0:TMAX_SIM]
                            beach_width = cascade._initial_beach_width[iB3D] / 10.0
                            Domain = barrier3d[iB3D].DomainTS[t] * 10.0
                            Dunes = (barrier3d[iB3D].DuneDomain[t, :, :] + barrier3d[iB3D].BermEl) * 10.0

                        # beach wedge
                        cellular_dune_toe_post_storm = np.floor(actual_shoreline_post_storm[t] + beach_width)
                        cellular_shoreline_post_storm = np.floor(actual_shoreline_post_storm[t])
                        cellular_beach_width = int(cellular_dune_toe_post_storm - cellular_shoreline_post_storm)

                        BeachDomain = np.zeros([max(0, cellular_beach_width), BarrierLength])
                        if cellular_beach_width > 0:
                            add = (barrier3d[iB3D].BermEl - barrier3d[iB3D].SL) / (cellular_beach_width + 1)
                            for i in range(cellular_beach_width):
                                BeachDomain[i, :] = (barrier3d[iB3D].SL + add) * (i + 1)

                        # assemble
                        Dunes_plot = np.flipud(np.rot90(Dunes))
                        Beach = BeachDomain * 10.0
                        Domain_stack = np.vstack([Beach, Dunes_plot, Domain])
                        Domain_stack[Domain_stack < 0] = -1

                        if flip_each_domain:
                            Domain_stack = Domain_stack[:, ::-1]

                        widthTS = len(Domain_stack)
                        OriginTstart = int(cellular_shoreline_post_storm)
                        OriginTstop = OriginTstart + widthTS
                        xOrigin = iB3D * BarrierLength
                        AnimateDomain[OriginTstart:OriginTstop, xOrigin : xOrigin + BarrierLength] = Domain_stack

                    # plot & save
                    elevFig1 = plt.figure(figsize=(7, 7) if fig_size is None else fig_size)
                    ax = elevFig1.add_subplot(111)
                    cax = ax.pcolormesh(AnimateDomain, cmap="terrain", vmin=-1.1, vmax=z_lim)
                    cbar = elevFig1.colorbar(cax)
                    cbar.set_label("elevation (m MHW)", rotation=270)

                    # Seams + labels
                    if draw_domain_seams:
                        for k in range(1, ny):
                            ax.axvline(x=k * BarrierLength, color=(1, 1, 1, 0.25), linewidth=0.8)

                    plt.xlabel("alongshore distance (dam)")
                    plt.ylabel("cross-shore distance (dam)")
                    timestr = f"Time = {t - 0.5} yrs"
                    if y_lim is not None:
                        plt.ylim(y_lim)
                        plt.text(3, y_lim[0] + 3, timestr, color="w")
                    else:
                        plt.ylim(bottom=OriginY - 35)
                        plt.text(1, OriginY - 33, timestr, color="w")

                    # Labels after ylim is set
                    _decorate_axes(ax)

                    plt.tight_layout()
                    plt.rcParams.update({"font.size": 11})
                    fname = f"elev_{t - 1}pt5"
                    elevFig1.savefig(fname + (".eps" if fig_eps else ".png"), format=("eps" if fig_eps else None))
                    plt.close(elevFig1)

        # ---------- annual (t) frames ----------
        for t in range(TMAX_SIM):
            AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1
            for iB3D in range(ny):
                actual_shoreline_post_humans = barrier3d[iB3D].x_s_TS[0 : TMAX_SIM + 1]

                if beach_management_ny[iB3D]:
                    beach_width = cascade.nourishments[iB3D].beach_width[t] / 10.0
                    if np.isnan(beach_width):
                        beach_width = cascade.nourishments[iB3D].beach_width[t - 1] / 10.0
                else:
                    beach_width = cascade._initial_beach_width[iB3D] / 10.0

                cellular_dune_toe_post_humans = np.floor(actual_shoreline_post_humans[t] + beach_width)
                cellular_shoreline_post_humans = np.floor(actual_shoreline_post_humans[t])
                cellular_beach_width = int(cellular_dune_toe_post_humans - cellular_shoreline_post_humans)

                if t == maxMGMT + 1:
                    cellular_beach_width_final = cellular_beach_width
                if t > maxMGMT + 1:
                    cellular_beach_width = cellular_beach_width_final

                BeachDomain = np.zeros([max(0, cellular_beach_width), BarrierLength])
                if cellular_beach_width > 0:
                    add = (barrier3d[iB3D].BermEl - barrier3d[iB3D].SL) / (cellular_beach_width + 1)
                    for i in range(cellular_beach_width):
                        BeachDomain[i, :] = (barrier3d[iB3D].SL + add) * (i + 1)

                Domain = barrier3d[iB3D].DomainTS[t] * 10.0
                Dunes = (barrier3d[iB3D].DuneDomain[t, :, :] + barrier3d[iB3D].BermEl) * 10.0
                Dunes = np.flipud(np.rot90(Dunes))
                Beach = BeachDomain * 10.0
                Domain_stack = np.vstack([Beach, Dunes, Domain])
                Domain_stack[Domain_stack < -3] = -3
                Domain_stack[Domain_stack > 6] = 6

                if flip_each_domain:
                    Domain_stack = Domain_stack[:, ::-1]

                widthTS = len(Domain_stack)
                OriginTstart = int(cellular_shoreline_post_humans)
                OriginTstop = OriginTstart + widthTS
                xOrigin = iB3D * BarrierLength
                AnimateDomain[OriginTstart:OriginTstop, xOrigin : xOrigin + BarrierLength] = Domain_stack

            # plot & save
            elevFig2 = plt.figure(figsize=(7, 7) if fig_size is None else fig_size)
            ax = elevFig2.add_subplot(111)
            cax = ax.pcolormesh(AnimateDomain, cmap="terrain", vmin=-1.1, vmax=z_lim)
            cbar = elevFig2.colorbar(cax)
            cbar.set_label("elevation (m MHW)", rotation=270)

            # Seams first
            if draw_domain_seams:
                for k in range(1, ny):
                    ax.axvline(x=k * BarrierLength, color=(1, 1, 1, 0.25), linewidth=0.8)

            plt.xlabel("alongshore distance (km)")
            plt.ylabel("cross-shore distance (km)")
            timestr = f"Time = {t} yrs"
            if y_lim is not None:
                plt.ylim(y_lim)
                plt.text(3, y_lim[0] + 3, timestr, color="w")
            else:
                plt.ylim(bottom=OriginY - 35)
                plt.text(1, OriginY - 33, timestr, color="w")

            # Labels after ylim is set
            _decorate_axes(ax)

            # relabel ticks (dam -> km)
            xticks = plt.xticks()
            yticks = plt.yticks()
            xtick_loc = list(xticks[0])
            ytick_loc = list(yticks[0])

            xmin, xmax = 0, len(AnimateDomain[0])
            focused_x_tick_location = [x for x in xtick_loc if xmin <= x <= xmax]
            New_X_Labels = np.divide(focused_x_tick_location, 100.0)
            plt.xticks(focused_x_tick_location, labels=New_X_Labels)

            Y_vals = plt.ylim()
            Y_min, Y_max = min(Y_vals), max(Y_vals)
            focused_y_tick_location = [y for y in ytick_loc if Y_min <= y <= Y_max]
            New_Y_Labels = np.divide(focused_y_tick_location, 100.0)
            plt.yticks(focused_y_tick_location, labels=New_Y_Labels)

            plt.tight_layout()
            plt.rcParams.update({"font.size": 11})
            fname = f"elev_{t}"
            elevFig2.savefig(fname + (".eps" if fig_eps else ".png"), format=("eps" if fig_eps else None))
            plt.close(elevFig2)

        # ---------- assemble GIF ----------
        frames = []
        for filenum in range(TMAX_SIM):
            filename = f"elev_{filenum}" + (".eps" if fig_eps else ".png")
            frames.append(imageio.imread(filename))
            for iB3D in range(ny):
                if (beach_management_ny[iB3D] or roadway_management_ny[iB3D]) and (filenum < TMAX_MGMT[iB3D]):
                    filename = f"elev_{filenum}pt5" + (".eps" if fig_eps else ".png")
                    frames.append(imageio.imread(filename))

        imageio.mimsave("elev.gif", frames, "GIF", fps=2)
        print("\n[ * GIF successfully generated * ]")
    finally:
        os.chdir(cwd_before)


# ------------------------ HATTERAS DRIVER (natural-only) ------------------------

NPZ_PATH = r"/output/Hindcast_1978_1997_buffers/Hindcast_1978_1997_buffers.npz"
DIRECTORY = r"C:\Users\hanna\PycharmProjects\CASCADE"

output = np.load(NPZ_PATH, allow_pickle=True)
cascade = output["cascade"][0]
b3d = cascade.barrier3d

# Hatteras: buffers 0..14; real 15..106; buffers 107..121
Model_Grids_Of_Interest = list(range(55, 65))  # 10 real domains; adjust as desired
ny = len(Model_Grids_Of_Interest)

TMAX_SIM = 23
TMAX_MGMT = [0] * ny
beach_management_ny = [False] * ny
roadway_management_ny = [False] * ny

run_name = "Hindcast_1978_1997_buffers"

plot_ElevAnimation_CASCADE(
    cascade=cascade,
    directory=DIRECTORY,
    TMAX_MGMT=TMAX_MGMT,
    name=run_name,
    TMAX_SIM=TMAX_SIM,
    Model_Grids_Of_Interest=Model_Grids_Of_Interest,
    ny=ny,
    beach_management_ny=beach_management_ny,
    roadway_management_ny=roadway_management_ny,
    y_lim=None,
    z_lim=3.5,
    fig_size=None,
    fig_eps=False,
    reverse_alongshore=True,   # your fix
    flip_each_domain=False,
    draw_domain_seams=True,
    draw_domain_labels=True,   # <-- turn labels on
)
