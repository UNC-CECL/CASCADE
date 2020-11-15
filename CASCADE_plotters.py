# Plotting functions for

# ~******* CASCADE ********~

"""----------------------------------------------------
Copyright (C) 2020 Katherine Anarde
Full copyright notice located in main CASCADE.py file
----------------------------------------------------"""

# Simulation Functions Included:
#   Brie_gif            Creates a gif of barrier island evolution from BRIE outputs

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import numpy as np
import os
import imageio
import math
from scipy import signal

#===================================================
# 1: Animation Frames of Barrier and Dune Elevation (#4 in Barrier3D_Functions, modified here for CASCADE)
#
#       inputs:         - barrier3d (a list containing BMI objects)
#                       - ny (the number of barrier3d subgrids that you want to plot)
#                       - directory (for saving)
#                       - TMAX (the last time index for plotting)
#                       - name (for saving new directory)
#       outputs:        - gif

def plot_ElevAnimation(barrier3d, ny, directory, TMAX, name):

    BarrierLength = barrier3d[0]._BarrierLength

    BeachWidth = 6
    OriginY = 10
    AniDomainWidth = int(np.amax(barrier3d[0]._InteriorWidth_AvgTS) + BeachWidth + np.abs(barrier3d[0]._ShorelineChange) + OriginY + 35)

    os.chdir(directory)
    newpath = 'Output/' + name + '/SimFrames/'
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)


    for t in range(TMAX-1):

        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1

        for iB3D in range(ny):
            # Build beach elevation domain
            BeachDomain = np.zeros([BeachWidth, BarrierLength])
            berm = math.ceil(BeachWidth * 0.65)
            BeachDomain[berm:BeachWidth + 1, :] = barrier3d[iB3D]._BermEl
            add = (barrier3d[iB3D]._BermEl - barrier3d[iB3D]._SL) / berm
            for i in range(berm):
                BeachDomain[i, :] = barrier3d[iB3D]._SL + add * i

            # Make animation frame domain
            Domain = barrier3d[iB3D]._DomainTS[t] * 10
            Dunes = (barrier3d[iB3D]._DuneDomain[t, :, :] + barrier3d[iB3D]._BermEl) * 10
            Dunes = np.rot90(Dunes)
            Dunes = np.flipud(Dunes)
            Beach = BeachDomain * 10
            Domain = np.vstack([Beach, Dunes, Domain])
            Domain[Domain < 0] = -1
            widthTS = len(Domain)
            scts = [(x - barrier3d[iB3D]._x_s_TS[0]) for x in barrier3d[iB3D]._x_s_TS]
            if scts[t] >= 0:
                OriginTstart = OriginY + math.floor(scts[t])
            else:
                OriginTstart = OriginY + math.ceil(scts[t])
            OriginTstop = OriginTstart + widthTS
            #AnimateDomain[OriginTstart:OriginTstop, 0: barrier3d[iB3D]._BarrierLength] = Domain
            xOrigin = iB3D * BarrierLength
            AnimateDomain[OriginTstart:OriginTstop, xOrigin : xOrigin + BarrierLength] = Domain


        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(AnimateDomain, origin='lower', cmap='terrain', vmin=-1.1,
                         vmax=4.0)  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position('bottom')
        # cbar = elevFig1.colorbar(cax)
        plt.xlabel('Alongshore Distance (dam)')
        plt.ylabel('Cross-Shore Diatance (dam)')
        plt.title('Interior Elevation')
        plt.tight_layout()
        timestr = 'Time = ' + str(t) + ' yrs'
        plt.text(1, 1, timestr)
        name = 'elev_' + str(t)
        elevFig1.savefig(name)  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(TMAX-1):
        filename = 'elev_' + str(filenum) + '.png'
        frames.append(imageio.imread(filename))
    imageio.mimsave('elev.gif', frames, 'GIF-FI')
    print()
    print('[ * GIF successfully generated * ]')


# ===================================================
# 2: Shoreline positions over time (#6 in Barrier3D_Functions)
#
#       inputs:         - x_s_TS, x_b_TS (shoreline and back-barrier time series for one B3D subgrid)
#       outputs:        - fig

def plot_ShorelinePositions(x_s_TS, x_b_TS):
    scts = [(x - x_s_TS[0]) * -10 for x in x_s_TS]
    bscts = [(x - x_s_TS[0]) * -10 for x in x_b_TS]
    shorelinefig = plt.figure()
    plt.plot(scts, 'b')
    plt.plot(bscts, 'g')
    fig = plt.gcf()
    fig.set_size_inches(14, 8)
    plt.ylabel('Shoreline Position (m)')
    plt.xlabel('Year')
    plt.show()
    name = 'Output/Shorelines'
    #shorelinefig.savefig(name)


# ===================================================
# #3 B3D cross-Shore Transect for one subgrid every 100 m for last time step (#5 in Barrier3D_Functions)
#
#       inputs:         - barrier3d (a list containing BMI objects)
#                       - TMAX (the last time index for plotting)
#       outputs:        - fig

def plot_XShoreTransects(barrier3d, TMAX):
    # Build beach elevation
    BW = 6  # beach width (dam) for illustration purposes
    BeachX = np.zeros(BW)
    berm = math.ceil(BW * 0.5)
    BeachX[berm:BW + 1] = barrier3d._BermEl
    add = (barrier3d._BermEl - barrier3d._SL) / berm
    for i in range(berm):
        BeachX[i] = barrier3d._SL + add * i
        # Plot full tranects
    plt.figure()

    for v in range(0, barrier3d._BarrierLength, 20):
        CrossElev = barrier3d._InteriorDomain[:, v]
        Dunes = barrier3d._DuneDomain[TMAX, v, :] + barrier3d._BermEl
        CrossElev1 = np.insert(CrossElev, 0, Dunes)
        CrossElev2 = np.insert(CrossElev1, 0, BeachX)
        CrossElev = CrossElev2 * 10  # Convert to meters
        plt.plot(CrossElev)
    fig = plt.gcf()
    fig.set_size_inches(14, 6)
    plt.hlines(barrier3d._SL, -1, len(CrossElev + 1), colors='dodgerblue')
    plt.xlabel('Cross-Shore Distance (dam)')
    plt.ylabel('Elevation (m)')
    plt.title('Cross-shore Topo Transects')
    plt.show()
    name = 'Output/Profiles'
    # fig.savefig(name)

# ===================================================
# #4: Cross-shore transects for both brie and B3d
#
#       inputs:         - barrier3d (a list of B3D objects)
#                       - brieLTA (a brie object with LTA model on)
#                       - time_step (a list of time steps to plot)
#                       - iB3D (an integer corresponding to the B3D subdomain / brie alongshore grid)
#       outputs:        - fig

def plot_ModelTransects(b3d, brieLTA, time_step, iB3D):

    plt.figure(figsize=(10, 5))
    colors = mpl.cm.viridis(np.linspace(0, 1, b3d[0]._TMAX))

    plt.subplot(2, 1, 1)

    for t in time_step:

        # Sea level
        SL = (b3d[iB3D]._SL + (t * b3d[iB3D]._RSLR[t]))

        # Create data points
        Tx = b3d[iB3D]._x_t_TS[t]
        Ty = (SL - b3d[iB3D]._DShoreface) * 10
        Sx = b3d[iB3D]._x_s_TS[t]
        Sy = SL * 10
        Bx = b3d[iB3D]._x_b_TS[t]
        By = (SL - b3d[iB3D]._BayDepth) * 10
        My = By

        BW = 6  # beach width (dam) for illustration purposes
        BeachX = np.zeros(BW)
        berm = math.ceil(BW * 0.5)
        BeachX[berm:BW + 1] = b3d[iB3D]._BermEl
        add = (b3d[iB3D]._BermEl - b3d[iB3D]._SL) / berm
        for i in range(berm):
            BeachX[i] = b3d[iB3D]._SL + add * i

        v = 0  # just use the first transect
        CrossElev = b3d[iB3D]._DomainTS[t]
        CrossElev = CrossElev[:, v]
        Dunes = b3d[iB3D]._DuneDomain[t, v, :] + b3d[iB3D]._BermEl
        CrossElev1 = np.insert(CrossElev, 0, Dunes)
        CrossElev2 = np.insert(CrossElev1, 0, BeachX)
        CrossElev = (CrossElev2 * 10) + Sy  # Convert to meters
        xCrossElev = np.arange(0, len(CrossElev), 1) + Sx

        Mx = xCrossElev[-1] + 20  # just add a buffer to the end of the plt

        x = np.hstack([Tx, Sx, xCrossElev, Mx])
        y = np.hstack([Ty, Sy, CrossElev, My])

        # Plot
        plt.plot(x, y, color=colors[t])
        plt.hlines(SL * 10, Tx, Mx, colors='dodgerblue')
        # plt.hlines((b3d[iB3D]._SL + (t * b3d[iB3D]._RSLR[t])) * 10, Tx, Sx, colors='dodgerblue')
        # plt.hlines((b3d[iB3D]._SL + (t * b3d[iB3D]._RSLR[t])) * 10, Bx + 2, xCrossElev[-1]+20, colors='dodgerblue')  # KA: scrappy fix

        plt.ylabel('Elevation (m)')
        plt.title('Profile Evolution')
        plt.show()

    plt.subplot(2, 1, 2)

    for t in time_step:
        # BRIE

        # Sea level
        SL = (b3d[iB3D]._SL + (t * b3d[iB3D]._RSLR[t]))

        # Create data points
        Tx = brieLTA._x_t_save[iB3D, t] / 10  # convert to dam
        Ty = (SL * 10) - brieLTA._d_sf  # units of m
        Sx = brieLTA._x_s_save[iB3D, t] / 10
        Sy = SL * 10
        Bx = brieLTA._x_b_save[iB3D, t] / 10
        By = (SL * 10) - brieLTA._bb_depth
        Hx1 = Sx
        Hy1 = brieLTA._h_b_save[iB3D, t] + Sy
        Hx2 = Bx
        Hy2 = brieLTA._h_b_save[iB3D, t] + Sy
        Mx = Bx + 20
        My = By

        x = [Tx, Sx, Hx1, Hx2, Bx, Mx]
        y = [Ty, Sy, Hy1, Hy2, By, My]

        # Plot
        plt.plot(x, y, color=colors[t])
        plt.hlines(SL * 10, Tx, Mx, colors='dodgerblue')
        # plt.hlines((b3d[iB3D]._SL + (t * b3d[iB3D]._RSLR[t])) * 10, Tx, Sx, colors='dodgerblue')
        # plt.hlines((b3d[iB3D]._SL + (t * b3d[iB3D]._RSLR[t])) * 10, Bx, Mx, colors='dodgerblue')

        plt.xlabel('Alongshore Distance (dam)')
        plt.ylabel('Elevation (m)')
        plt.show()

# ===================================================
# #5: Statistics from B3D
#
#       inputs:         - b3d (a list of B3D objects)
#                       - iB3D (an integer corresponding to the B3D subdomain / brie alongshore grid)
#                       - TMAX (the last time index for plotting)
#                       - iB3D
#       outputs:        - fig

def plot_statistics(b3d, iB3D, TMAX):

    colors = mpl.cm.viridis(np.linspace(0, 1, 10))
    plt.figure(figsize=(10, 10))

    # A Shoreface Slope
    plt.subplot(3, 2, 1)
    ssfTS = b3d[iB3D]._s_sf_TS
    plt.plot(ssfTS, color=colors[1])
    plt.hlines(b3d[iB3D]._s_sf_eq, 0, TMAX, colors='black', linestyles='dashed')
    #plt.ylabel(r'$\alpha$')
    plt.ylabel('Shoreface Slope')
    plt.legend(['B3D sub-grid #'+str(iB3D)])
    plt.rcParams["legend.loc"] = 'lower right'

    # A Interior Width
    plt.subplot(3, 2, 2)
    aiw = [a * 10 for a in b3d[iB3D]._InteriorWidth_AvgTS]
    plt.plot(aiw, color=colors[2])
    plt.ylabel('Avg. Width (m)')  # Average Interior Width
    plt.legend(['B3D sub-grid #'+str(iB3D)])
    plt.rcParams["legend.loc"] = 'lower right'

    # A Shoreline Change
    scts = [(x - b3d[iB3D]._x_s_TS[0]) * 10 for x in b3d[iB3D]._x_s_TS]
    plt.subplot(3, 2, 3)
    plt.plot(scts, color=colors[3])
    plt.ylabel('Shoreline Position (m)')
    plt.legend(['B3D sub-grid #'+str(iB3D)])
    plt.rcParams["legend.loc"] = 'lower right'

    # A Dune Height
    aHd = [a * 10 for a in b3d[iB3D]._Hd_AverageTS]
    plt.subplot(3, 2, 4)
    plt.plot(aHd, color=colors[4])
    plt.ylabel('Avg. Dune Height (m)')  # Average Dune Height
    plt.legend(['B3D sub-grid #'+str(iB3D)])
    plt.rcParams["legend.loc"] = 'lower right'

    # Qoverwash for entire B3D grid
    plt.subplot(3, 2, 5)
    Qoverwash = np.zeros(np.size(b3d[0]._QowTS[0:TMAX]))
    for iGrid in range(len(b3d)):
        Qoverwash = Qoverwash + (np.array(b3d[iGrid]._QowTS[0:TMAX]) * (b3d[iGrid]._BarrierLength * 10)) # m^3/yr
    Qoverwash = Qoverwash / (len(b3d) * b3d[0]._BarrierLength * 10) # m^3/m/yr
    plt.plot(Qoverwash, color=colors[5])
    #    movingavg = np.convolve(QowTS, np.ones((50,))/50, mode='valid')
    #    movingavg = [i * 10 for i in movingavg]
    #    plt.plot(movingavg, 'r--')
    plt.ylabel(r'Qow ($m^3/m/yr$)')
    plt.xlabel('Year')

    # Shoreface flux for entire B3D grid
    plt.subplot(3, 2, 6)
    Qsf = np.zeros(np.size(b3d[0]._QsfTS[0:TMAX]))
    for iGrid in range(len(b3d)):
        Qsf = Qsf + (np.array(b3d[iGrid]._QsfTS[0:TMAX]) * (b3d[iGrid]._BarrierLength * 10)) # m^3/yr
    Qsf = Qsf / (len(b3d) * b3d[0]._BarrierLength * 10) # m^3/m/yr
    plt.plot(Qsf, color=colors[6])
    plt.ylabel(r'Qsf ($m^3/m/yr$)')
    plt.xlabel('Year')

    plt.show()

# ===================================================
# #6: Statistics from BRIE with LTA
#
#       inputs:         - brieLTA (brie model with LTA on)
#                       - iB3D (an integer corresponding to the B3D subdomain / brie alongshore grid)
#                       - TMAX (the last time index for plotting and calculations)
#       outputs:        - fig

def plot_statistics_BRIE(brieLTA, iB3D, TMAX):

    colors = mpl.cm.viridis(np.linspace(0, 1, 10))
    plt.figure(figsize=(10, 10))

    # A Shoreface Slope
    plt.subplot(3, 2, 1)
    ssfTS = brieLTA._s_sf_save[iB3D,:]
    plt.plot(ssfTS, color=colors[1])
    plt.hlines(brieLTA._s_sf_eq, 0, len(brieLTA._s_sf_save[iB3D,:]), colors='black', linestyles='dashed')
    #plt.ylabel(r'$\alpha$')
    plt.ylabel('Shoreface Slope')
    plt.legend(['BRIE sub-grid #'+str(iB3D)])
    plt.rcParams["legend.loc"] = 'lower right'

    # A Interior Width
    plt.subplot(3, 2, 2)
    aiw = brieLTA._x_b_save[iB3D,:] - brieLTA._x_s_save[iB3D,:]
    plt.plot(aiw, color=colors[2])
    plt.ylabel('Barrier Width (m)')
    plt.legend(['BRIE sub-grid #'+str(iB3D)])
    plt.rcParams["legend.loc"] = 'lower right'

    # A Shoreline Change
    scts = [(x - brieLTA._x_s_save[iB3D,0]) for x in brieLTA._x_s_save[iB3D,:]]
    plt.subplot(3, 2, 3)
    plt.plot(scts, color=colors[3])
    plt.ylabel('Shoreline Position (m)')
    plt.legend(['BRIE sub-grid #'+str(iB3D)])
    plt.rcParams["legend.loc"] = 'lower right'

    # Barrier Height
    aHd = brieLTA._h_b_save[iB3D,:]
    plt.subplot(3, 2, 4)
    plt.plot(aHd, color=colors[4])
    plt.ylabel('Barrier Height (m)')  # Average Dune Height
    plt.legend(['BRIE sub-grid #'+str(iB3D)])
    plt.rcParams["legend.loc"] = 'upper right'
    plt.xlabel('Year')

    # Qoverwash for entire BRIE grid
    plt.subplot(3, 2, 5)
    QoverwashLTA = brieLTA._Qoverwash[0:TMAX] / (brieLTA._ny * brieLTA._dy)  # from brie in m^3/yr --> m^3/m/yr
    plt.plot( brieLTA._t[0:TMAX], QoverwashLTA, color=colors[5])
    #    movingavg = np.convolve(QowTS, np.ones((50,))/50, mode='valid')
    #    movingavg = [i * 10 for i in movingavg]
    #    plt.plot(movingavg, 'r--')
    plt.ylabel(r'Qow ($m^3/m/yr$)')
    plt.xlabel('Year')

    # # Shoreface flux for entire B3D grid
    # plt.subplot(3, 2, 6)
    # Qsf = np.zeros(np.size(b3d[0]._QsfTS[0:TMAX]))
    # for iGrid in range(len(b3d)):
    #     Qsf = Qsf + (np.array(b3d[iGrid]._QsfTS[0:TMAX]) * (b3d[iGrid]._BarrierLength * 10)) # m^3/yr
    # Qsf = Qsf / (len(b3d) * b3d[0]._BarrierLength * 10) # m^3/m/yr
    # plt.plot(Qsf, color=colors[6])
    # plt.ylabel(r'Qsf ($m^3/m/yr$)')
    # plt.xlabel('Year')

    plt.show()

# ===================================================
# #7: Calculate shoreline change periodicity
#
#       inputs:         - x_s_TS (an single shoreline timeseries for one B3D subdomain)
#       outputs:        - Periodicity (period of punctuated retreat)
#                       - AvgFastDur (average duration of fast periods)
#                       - AvgSlowDur (average duration of slow periods)
#                       - Punc (boolean for punctuated retreat)

def calc_ShorelinePeriodicity(x_s_TS):

    # Shoreline Change & Change Rate Over Time
    scts = [(x - x_s_TS[0]) * 10 for x in x_s_TS]

    # Filter
    win = 31
    poly = 3
    der1 = (signal.savgol_filter(scts, win, poly, deriv=1))

    HitDown = []  # Slow-downs
    HitUp = []  # Speed-ups

    window1 = 3  # Maximum allowed length for gaps in slow periods
    window2 = 30  # Minimum length required for slow periods, including gaps
    buffer = 3
    thresh1 = 0.5  # Max slope for slow periods (if shoreline change rate is below this, immobile, otherwise transgressive)
    thresh2 = 1

    # Find slow periods
    der_under = np.where(der1 < thresh1)[0]

    if len(der_under) > 0:

        gaps = np.diff(der_under) > window1
        peak_start = np.insert(der_under[1:][gaps], 0, der_under[0])
        peak_stop = np.append(der_under[:-1][gaps], der_under[-1])

        for n in range(len(peak_stop)):
            if peak_stop[n] - peak_start[n] > window2:
                if len(HitDown) == 0:
                    if peak_start[n] > buffer:
                        HitDown.append(peak_start[n])
                    if peak_stop[n] < len(scts) - buffer:
                        HitUp.append(peak_stop[n])
                else:
                    gap_length = peak_start[n] - HitUp[-1]
                    gap_slope = (scts[peak_start[n]] - scts[HitUp[-1]]) / gap_length
                    if gap_length < window2 and gap_slope < thresh2:
                        if peak_stop[n] < len(scts) - buffer:
                            HitUp[-1] = peak_stop[n]
                        else:
                            del HitUp[-1]
                    else:
                        if peak_start[n] > buffer:
                            HitDown.append(peak_start[n])
                        if peak_stop[n] < len(scts) - buffer:
                            HitUp.append(peak_stop[n])

    ### CALCULATE STATS

    Jumps = len(HitDown)
    Slows = len(HitUp)
    SlowDur = []
    FastDur = []

    if Jumps > 0 and Slows > 0:
        DownFirst = HitDown[0] < HitUp[0]
    elif Jumps == 0 and Slows > 1:
        DownFirst = True
    else:
        DownFirst = False

    if Jumps >= 2 or Slows >= 2:
        if Jumps >= 2 and Slows >= 2:
            Periodicity = (np.mean(np.diff(HitDown)) + np.mean(np.diff(HitUp))) / 2
        elif Jumps >= Slows:
            Periodicity = np.mean(np.diff(HitDown))
        else:
            Periodicity = np.mean(np.diff(HitUp))
        if DownFirst:
            for n in range(Slows):
                SlowDur.append(HitUp[n] - HitDown[n])
            for n in range(Jumps - 1):
                FastDur.append(HitDown[n + 1] - HitUp[n])
        else:
            for n in range(Slows - 1):
                SlowDur.append(HitUp[n + 1] - HitDown[n])
            for n in range(Jumps):
                FastDur.append(HitDown[n] - HitUp[n])
    else:
        Periodicity = 0
        if Jumps == 1 and Slows == 1:
            if DownFirst:
                SlowDur.append(HitUp[0] - HitDown[0])
            else:
                FastDur.append(HitDown[0] - HitUp[0])

    AvgFastDur = np.mean(FastDur)
    if np.isnan(AvgFastDur):
        AvgFastDur = 0
    AvgSlowDur = np.mean(SlowDur)
    if np.isnan(AvgSlowDur):
        AvgSlowDur = 0

    if len(SlowDur) >= 2 and len(FastDur) >= 2:
        Punc = 1
    else:
        Punc = 0

    return Periodicity, AvgFastDur, AvgSlowDur, Punc

# ===================================================
# #8: Dune Height Over Time for CASCADE
#
#       inputs:         - barrier3d (a list containing BMI objects)
#                       - TMAX (the last time index for plotting)
#       outputs:        - fig (dune domain for all sub-grids)

def plot_dune_domain(b3d, TMAX):

    DuneCrest = []
    Dmax = []

    for iB3D in range(len(b3d)):
        sub_domain = b3d[iB3D]._DuneDomain[0:TMAX, :, :]
        DuneCrest.append(sub_domain.max(axis=2))
        Dmax.append(b3d[iB3D]._Dmax)

    DuneCrest = np.hstack(DuneCrest).astype(float)
    Dmax = np.max(Dmax)

    duneFig = plt.figure(figsize=(14, 8))
    plt.rcParams.update({'font.size': 13})
    ax = duneFig.add_subplot(111)
    ax.matshow((DuneCrest) * 10, origin='lower', cmap='bwr', aspect='auto', vmin=0, vmax=Dmax * 10)
    cax = ax.xaxis.set_ticks_position('bottom')  # analysis:ignore
    # cbar = duneFig.colorbar(cax)
    # cbar.set_label('Dune Height Above Berm Elevation (m)', rotation=270)
    plt.xlabel('Alongshore Distance (dam)')
    plt.ylabel('Year')
    plt.title('Dune Height (m)')