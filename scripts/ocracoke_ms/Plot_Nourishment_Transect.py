import os
import numpy as np
import matplotlib.pyplot as plt
import math

def fig3_initialCNH_topo(
    cascade_model_list,  # must be from a nourishment simulation (i.e., have a beach width)
    km_on=True,
    year_wanted = 0,
):
    fig, axs = plt.subplots(1, 4, figsize=(10, 3), sharey=True, sharex=True)

    # make the mat image of the beach to the back-barrier; stole this code from the animation plots above
    for iCascade in range(len(cascade_model_list)):
        cascade = cascade_model_list[iCascade]
        barrier3d = cascade.barrier3d[35:52]
        domain_num = 35
        BarrierLength = barrier3d[0]._BarrierLength
        OriginY = 16
        AniDomainWidth = 120  # end Y
        z = int(year_wanted)
        ny = 1
        iB3D = 0
        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1
        scts = 0  # make shoreline position relative to zero starting

        # Build beach elevation domain, we only show beach width decreasing in increments of 10 m and we don't
        # illustrate a berm, just a sloping beach up to the elevation of the berm
        BeachWidth = math.floor(
            cascade.nourishments[domain_num].beach_width[0] / 10
        )  # switched from ceil
        BeachDomain = np.zeros(
            [
                BeachWidth,
                BarrierLength,
            ]
        )
        # beach width must be greater than zero
        add = (barrier3d[iB3D].BermEl - barrier3d[iB3D]._SL) / (BeachWidth + 1)
        for i in range(0, BeachWidth):
            BeachDomain[i, :] = (barrier3d[iB3D]._SL + add) * (i + 1)

        # Make animation frame domain
        Domain = barrier3d[iB3D]._DomainTS[z] * 10  # m MHW
        Dunes = (
            barrier3d[iB3D]._DuneDomain[z, :, :] + barrier3d[iB3D]._BermEl
        ) * 10  # m MHW
        Dunes = np.rot90(Dunes)
        Dunes = np.flipud(Dunes)
        Beach = BeachDomain * 10
        Domain = np.vstack([Beach, Dunes, Domain])
        Domain[Domain < 0] = -1  # anything underwater
        widthTS = len(Domain)
        OriginTstart = OriginY + math.floor(scts)  # ceil
        OriginTstop = OriginTstart + widthTS
        xOrigin = iB3D * BarrierLength
        AnimateDomain[
            OriginTstart:OriginTstop, xOrigin : xOrigin + BarrierLength
        ] = Domain

        # plot
        print(np.max(AnimateDomain))
        cax = axs[iCascade].matshow(
            AnimateDomain, origin="lower", cmap="terrain", vmin=-1, vmax=3.0
        )  # , interpolation='gaussian') # analysis:ignore
        axs[iCascade].xaxis.set_ticks_position("bottom")
        axs[iCascade].set(xlabel="alongshore distance (dam)")
        if km_on:
            axs[iCascade].set(xlabel="alongshore distance (km)")
        axs[iCascade].set_ylim([40, 110])
        axs[iCascade].set_xlim([-1, 50])

    axs[0].set(ylabel="cross-shore distance (dam)")
    # cbar = fig.colorbar(cax)
    # cbar.set_label("elevation (m MHW)", rotation=270)
    # plt.tight_layout()
    if km_on:
        locs, _ = plt.yticks()
        plt.yticks(locs, locs / 100)
        locs, _ = plt.xticks()
        plt.xticks(locs[1:], locs[1:] / 100)
        axs[0].set(ylabel="cross-shore distance (km)")

    # now make the cross-section; stole this from the cross-section code above and modified
    v = 10  # just use the 10th transect
    fig, axs = plt.subplots(1, 2, figsize=(10, 3), sharey=True, sharex=True)

    for iCascade in range(len(cascade_model_list)):
        cascade = cascade_model_list[iCascade]

        sea_level = barrier3d[iB3D]._SL

        # Create data points
        shoreface_toe_x = (
            barrier3d[iB3D].x_t_TS[z] - barrier3d[iB3D].x_t_TS[z]
        )
        shoreface_toe_y = (sea_level - barrier3d[iB3D].DShoreface) * 10  # m
        shoreline_x = (
            barrier3d[iB3D].x_s_TS[z] - barrier3d[iB3D].x_t_TS[z]
        )
        shoreline_y = sea_level * 10  # m
        bay_y = (sea_level - barrier3d[iB3D]._BayDepth) * 10  # m
        end_of_bay_y = bay_y

        berm_x = shoreline_x + (
            cascade.nourishments[domain_num].beach_width[z] / 10
        )  # beach width (in dam)
        berm_y = (
            barrier3d[iB3D]._BermEl * 10
        ) + shoreline_y  # convert to meters
        dune_toe_x = berm_x
        dune_toe_y = berm_y

        interior_y = barrier3d[iB3D]._DomainTS[0]
        interior_y = interior_y[:, v]
        print(
            np.max(barrier3d[iB3D]._DuneDomain[z, v, :] * 10)
        )  # max dune height
        dunes_y = (
            barrier3d[iB3D]._DuneDomain[z, v, :]
            + barrier3d[iB3D]._BermEl
        )
        cross_barrier_y = np.insert(interior_y, 0, dunes_y)
        cross_barrier_y = (cross_barrier_y * 10) + shoreline_y  # Convert to meters
        cross_barrier_x = np.arange(0, len(cross_barrier_y), 1) + dune_toe_x

        end_of_bay_x = (
            cross_barrier_x[-1] + 20
        )  # just add a buffer to the end of the plt

        x = np.hstack(
            [
                shoreface_toe_x,
                shoreline_x,
                berm_x,
                dune_toe_x,
                cross_barrier_x,
                end_of_bay_x,
            ]
        )
        y = np.hstack(
            [
                shoreface_toe_y,
                shoreline_y,
                berm_y,
                dune_toe_y,
                cross_barrier_y,
                end_of_bay_y,
            ]
        )

        # Plot
        if iCascade < 2:  # just the 0.45 cases
            axs[0].plot(x, y)
            axs[0].hlines(
                sea_level * 10, shoreface_toe_x, end_of_bay_x, colors="dodgerblue"
            )
            # axs[0].set_xlim([0, 110])
            axs[0].set(xlabel="cross-shore distance (dam)")
            axs[0].set(ylabel="elevation (m MHW)")
        else:  # 0.75 cases
            axs[1].plot(x, y)
            axs[1].hlines(
                sea_level * 10, shoreface_toe_x, end_of_bay_x, colors="dodgerblue"
            )
            # axs[1].set_xlim([0, 110])
            axs[1].set(xlabel="cross-shore distance (dam)")
    plt.tight_layout()
    plt.savefig(fname=('C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Figures\\Nourishment_11.eps'),format='eps')
    plt.show()

    if km_on:
        locs, _ = plt.xticks()
        plt.xticks(locs[1:], locs[1:] / 100)
        axs[0].set(xlabel="cross-shore distance (km)")
        axs[1].set(xlabel="cross-shore distance (km)")
        axs[0].set_xlim([-1, 141])
        axs[1].set_xlim([-1, 141])
        axs[0].legend(["profile A", "profile B"])
        axs[1].legend(["profile C", "profile D"])



os.chdir('C:\\Users\\frank\\PycharmProjects\\CASCADE\\Run_output')


run_name_batch='OCR_IH_Nourishment_S49_Erosional_Sink'

number_barrier3d_models = 70
buffer_length = 15


# --------- plot ---------
output = np.load(run_name_batch + ".npz", allow_pickle=True)
cascade = output["cascade"]
cascade = cascade[0]

fig3_initialCNH_topo(cascade_model_list=[cascade],
                     year_wanted = 11)