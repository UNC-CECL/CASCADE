# Lexi (Van Blunk) Fiegelist
# last updated: June 04, 2026

"""Simulate marsh dynamics in CASCADE

This module provides functions for modifying a barrier segment from Barrier3D --
consisting of an interior grid with marsh cells to replicate marsh dynamics
including:

1. marsh platform accretion based on the sum of mineral and organic contributions of mud-sized particles
2. marsh platform decomposition based on organic sediment within the marsh soil profile


References
----------

.. [1]

.. [2]

.. [3]


Notes
-----
Here, we assume all contributions from suspended sediment in the bay are mineral

"""

import math
import numpy as np
from matplotlib import pyplot as plt
import copy

def evolvemarsh(
        marshelevation,     # elevation domain, here it is ONLY marsh cells
        msl,                # msl in m
        C_e,                # SSC at marsh edge (kg/m3)
        OCb,                # organic content of uppermost layer of bay sediment (%), currently set to 0.05
        tr,                 # tidal range, m
        numiterations,      # tidal iterations, currently set to 500
        P,                  # tidal period, 12.5 hours for semi-diurnal tide, converted to seconds
        ws,                 # settling velocity (m/s), currently set to 5 x 10-5
        timestep,           # 365 * (24 / 12.5)  # number of tidal cycles per year
        Bmax,               # g/m2 (I think this is maximum biomass)
        Dmin,               # minimum depth below high water that marsh vegetation can grow (m)
        Dmax,               # maximum depth below high water that marsh vegetation can grow (m)
        rhoo,               # organic matter bulk density (kg/m3), set to 85
        rhos,               # sediment bulk density (kg/m3), set to 2000
        plot,
):
    """Calculates biomass and mineral and organic deposition for each cell in the marsh as a function of flooding
    frequency; calculates the total flux of sediment onto the marsh from the bay."""
    plt.rcParams.update({"font.size": 16})

    # updates from original:
    # 1. convert elevation into meters
    # 2. we now have multiple transects that represent the marsh, so loop through all columns of the interior domain
    # 3. we need to identify the marsh cells

    dt = P/numiterations  # inundation time, seconds

    # np.save(r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\marsh_elevation.npy", marshelevation)
    L = len(marshelevation)
    sedimentcycle = np.zeros([numiterations, L])
    depth = np.zeros([numiterations, L])
    C = np.zeros([L])

    min_dep_method = 1  # see options below
    # 1 for original
    # 2 for fixed marsh edge
    # 2.5 for fixed marsh edge plus smoothing
    # 3 for no pond
    plot_on = plot

    # Loop through a tidal cycle to determine flooding duration for each point in the marsh
    # time_submerged is used to calculate floodfraction but that is never used itself,
    # BUT the depth variable is used to calculate sediment deposition, so we definitely need that still
    for i in range(1, numiterations):
        # tidal cycle depth + depth of water above marsh elevation for this msl
        # some cells are below msl so they are negative depths
        depth_iteration = 0.5 * tr * math.sin(2 * math.pi * ((i + 1) * dt / P)) + (msl - marshelevation[0: L + 1])  # [m] Depth at each position in the marsh
        depth[i, :] = depth_iteration  # Store depths in array

    # -------------------------
    # Belowground Productivity
    # Creates a biomass curve (Mariotti & Carr, 2014) where peak biomass occurs at a depth halfway between the maximum
    # depth for vegetation and the minimum (here, mean high water level)
    # we might want to just set this value instead of using tidal range
    # so instead of using HWL as our baseline, use SL
    # check the paper they reference
    dm = msl + tr / 2 - marshelevation  # [m] Depth of the marsh surface below HWL at any given point

    if plot_on:
        # plot marsh elevation and marsh depth on same figure
        fig, ax1 = plt.subplots(1, 1, figsize=(10, 5), sharex="all")
        x_values = np.arange(0, L, 1)
        # Plot marsh elevation on the left y-axis
        ax1.plot(x_values, marshelevation, 'r', label='marsh elevation')
        ax1.set_xlabel('distance from marsh edge (m)')
        ax1.set_ylabel('marsh elevation (m MSL)', color='r')
        ax1.tick_params(axis='y', labelcolor='r')
        # Create a second y-axis
        ax2 = ax1.twinx()
        ax2.plot(x_values, dm, 'b', label='dm')
        xmin = 0
        xmax = L
        plt.hlines(Dmax, xmin, xmax, ls="dashed", colors="grey")
        plt.hlines(Dmin, xmin, xmax, ls="dashed", colors="grey")
        ax2.text(xmax-1.5, Dmax, "Dmax", size="x-small", va="bottom", c="grey")
        ax2.text(xmax-1.5, Dmin, "Dmin", size="x-small", va="bottom", c="grey")
        ax2.set_ylabel('marsh depth (m)', color='b')
        ax2.tick_params(axis='y', labelcolor='b')

    bgb = np.zeros([L])  # Belowground biomass
    organic_autoch = np.zeros([L])  # Autochthonous organic material

    for ii in range(L):
        if dm[ii] > Dmax:  # If depth is above vegetation maximum, there is no production
            bgb[ii] = 0  # [g] Mudflat
            organic_autoch[ii] = 0  # [g] No autochthonous material stored in the soil; we do not multiply by lingin content (as in earlier version of the model) because we subtract mass due to decomposition in the 'decompose' function
        elif dm[ii] <= Dmin:  # If depth is below vegetation minimum, there is very little belowground productivity
            bgb[ii] = 0.1  # [g] "Forest" - really high marsh
            organic_autoch[ii] = bgb[ii]  # [g] Forest organic matter stored in soil
        else:
            bgb[ii] = Bmax * (dm[ii] - Dmax) * (dm[ii] - Dmin) / (0.25 * (-Dmin - Dmax) * (Dmax - 3 * Dmin))  # [g] Marsh
            organic_autoch[ii] = bgb[ii]  # [g] As mentioned above, no longer multiply by lingin content

    if plot_on:
        # FIRST SUBPLOT
        # plot marsh elevation and bgb on same figure
        fig, (ax1, ax3, ax5) = plt.subplots(3, 1, figsize=(10, 15), sharex="all")
        x_values = np.arange(0, L, 1)
        # Plot marsh elevation on the left y-axis
        ax1.plot(x_values, marshelevation, 'r', label='marsh elevation')
        ax1.set_ylabel('marsh elevation (m MSL)', color='r')
        ax1.tick_params(axis='y', labelcolor='r')
        # Create a second y-axis
        ax2 = ax1.twinx()
        ax2.plot(x_values, bgb, 'b', label='bgb')
        ax2.set_ylabel('biomass (g)', color='b')
        ax2.tick_params(axis='y', labelcolor='b')
        ax2.set_ylim([0, 2600])

    if min_dep_method == 1:
        # -------------------------
        # Mineral Deposition - Original
        coeff = -0.002  # Coefficient of -0.0031 is a fitted parameter for realistic marsh topography
        distance = 0  # [m] Initialize, distance from marsh edge
        pond_index = []
        pond_y = []
        pond = False
        for xx in range(L):
            if bgb[xx] > 0:
                pond = False
                distance += 1  # [m]
                C[xx] = C_e * math.exp(coeff * distance)  # [kg/m3] Concentration at each marsh cell. Coefficient of -0.0031 is a fitted parameter for realistic marsh topography
            else:
                if not pond:  # Allow sediment suply to beyond first ponded cell: this is an ALTERATION/NEW ADDITION not included in original Matlab CoLT version
                    C_e = C_e * 0.9  # Decrease concentration at the new "marsh edge" by 10% with each subsequent pond formation
                    pond = True
                distance = 1  # [m]
                C[xx] = C_e
                pond_index.append(xx)
                pond_y.append(C_e)

    elif min_dep_method == 2:
        # Mineral Deposition Option 2 - fixed marsh edge
        coeff = -0.002  # Coefficient of -0.0031 is a fitted parameter for realistic marsh topography
        distance = 0  # [m] Initialize, distance from marsh edge
        pond_index = []
        pond_y = []
        pond = False
        for xx in range(L):
            distance += 1  # [m]
            if bgb[xx] > 0:
                pond = False
                # distance += 1  # [m]
                C[xx] = C_e * math.exp(coeff * distance)  # [kg/m3] Concentration at each marsh cell. Coefficient of -0.0031 is a fitted parameter for realistic marsh topography
            else:
                if not pond:  # Allow sediment suply to beyond first ponded cell: this is an ALTERATION/NEW ADDITION not included in original Matlab CoLT version
                    C_e = C_e * 0.9  # Decrease concentration at the new "marsh edge" by 10% with each subsequent pond formation
                    pond = True
                # distance = 1  # [m]
                C[xx] = C_e
                pond_index.append(xx)
                pond_y.append(C_e)

    elif min_dep_method == 2.5:
        # Mineral Deposition Option 2.5 - fixed marsh edge
        coeff = -0.002  # Coefficient of -0.0031 is a fitted parameter for realistic marsh topography
        distance = 0  # [m] Initialize, distance from marsh edge
        pond_index = []
        pond_y = []
        pond = False
        for xx in range(L):
            distance += 1  # [m]
            if bgb[xx] > 0:
                pond = False
                # distance += 1  # [m]
                C[xx] = C_e * math.exp(coeff * distance)  # [kg/m3] Concentration at each marsh cell. Coefficient of -0.0031 is a fitted parameter for realistic marsh topography
            else:
                if not pond:  # Allow sediment suply to beyond first ponded cell: this is an ALTERATION/NEW ADDITION not included in original Matlab CoLT version
                    if xx == 0:
                        C_e = C_e * 0.9  # Decrease concentration at the new "marsh edge" by 10% with each subsequent pond formation
                    else:
                        C_e = C[xx-1] * 0.9  # Decrease concentration at the new "marsh edge" by 10% with each subsequent pond formation
                    pond = True
                # distance = 1  # [m]
                C[xx] = C_e
                pond_index.append(xx)
                pond_y.append(C_e)

    elif min_dep_method == 3:
        # Mineral Deposition Option 3 - no pond
        coeff = -0.002  # Coefficient of -0.0031 is a fitted parameter for realistic marsh topography
        distance = 0  # [m] Initialize, distance from marsh edge
        for xx in range(L):
            if bgb[xx] > 0:
                distance += 1  # [m]
                C[xx] = C_e * math.exp(
                    coeff * distance)  # [kg/m3] Concentration at each marsh cell. Coefficient of -0.0031 is a fitted parameter for realistic marsh topography

    else:
        print("not a valid mineral deposition method")

    if plot_on:
        #-------------------------------------------------------
        # SECOND SUBPLOT
        # plot marsh elevtion and Ce on same figure
        # Plot marsh elevation on the left y-axis
        ax3.plot(x_values, marshelevation, 'r', label='marsh')
        ax3.set_ylabel('marsh elevation (m MSL)', color='r')
        ax3.tick_params(axis='y', labelcolor='r')
        # Create a second y-axis
        ax4 = ax3.twinx()
        ax4.plot(x_values, C, 'b', label='SSC (kg/m3)')
        ax4.set_ylabel('SSC (kg/m3)', color='b')
        ax4.tick_params(axis='y', labelcolor='b')
        ax4.set_ylim([0, 0.051])
        # plot the ponds as stars
        if min_dep_method != 3:
            ax4.scatter(pond_index, pond_y, marker="*", c="b")

    for i in range(1, numiterations):
        tempdepth = depth[i, :]
        tempy = np.zeros([L])
        tempy[tempdepth > 0] = C[tempdepth > 0] * ws * dt  # [kg] mass of mineral sediment deposited, where depth > 0
        # not entirely sure how tempy is in kg bc C is kg/m3, ws is m/s, and dt is s, so it seems like it is kg/m2
        # except that it is a transect, so maybe assuming 1m x 1m?
        sedimentcycle[i, :] = tempy

    susp_dep = np.sum(sedimentcycle, axis=0) * timestep * 1000  # [g/yr] suspended sediment deposition in an entire year, in each cell
    mineral = susp_dep * (1 - OCb)  # [g] Mineral deposition of suspended sediment in a given year is equal determined by the organic content of the bay sediment
    organic_alloch = susp_dep * OCb  # [g] Organic deposition of suspended sediment is equal determined by the organic content of bay sed

    # -------------------------
    # Calculations & Conversions

    # # Total flux of sediment onto the marsh from the bay
    # Fm_min = np.sum(mineral) / 1000  # [kg/yr] Flux of mineral sediment from the bay
    # Fm_org = np.sum(organic_alloch) / 1000  # [kg/yr] Flux of organic sediment from the bay

    # Calculate thickness of new sediment (mineral+organic) based off LOI and its effect on density
    loi = (organic_autoch + organic_alloch) / (mineral + organic_autoch + organic_alloch)
    density = 1 / ((loi / rhoo) + ((1 - loi) / rhos)) * 1000  # [g/m3] Bulk density is calculated according to Morris et al. (2016)
    density[np.isnan(density)] = 1  # If there is zero deposition, loi calculation will divide by zero and make density nan. Set density equal to 1 in this case, so that accretion is zero, instead of nan.

    accretion = (mineral + organic_autoch + organic_alloch) / density  # [m] accretion in a given year

    if plot_on:
        # THIRD SUBPLOT
        # plot marsh elevation and bgb on same figure
        # Plot marsh elevation on the left y-axis
        ax5.plot(x_values, marshelevation, 'r', label='marsh elevation')
        ax5.set_xlabel('distance from marsh edge (m)')
        ax5.set_ylabel('marsh elevation (m MSL)', color='r')
        ax5.tick_params(axis='y', labelcolor='r')
        # Create a second y-axis
        ax6 = ax5.twinx()
        ax6.plot(x_values, accretion, 'b', label='accretion')
        ax6.set_ylabel('accretion (m)', color='b')
        ax6.tick_params(axis='y', labelcolor='b')
        ax6.set_ylim([0, 0.061])
        # Add title and show plot
        fig.tight_layout()
        plt.show()

        # save the figure
        if min_dep_method == 2.5:
            fig.savefig(r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\output_figures\option_2pt5_dmax_{0}".format(Dmax))
        else:
            fig.savefig(
                r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\output_figures\option_{0}_dmax_{1}".format(min_dep_method, Dmax))

    # Update elevation
    marshelevation += accretion

    return marshelevation, accretion, organic_autoch


def decompose(
        x_m,  # marsh position
        x_f,  # forest position
        yr,   # current year
        organic_dep_autoch,
        elevation,  # elevation domain
        B,  # transect width (marsh, bay, and forest)
        mui,  # depth in marsh where decomp goes to 0 (m)
        mki,  # decomp coefficient
        rhoo,  # organic matter bulk density (kg/m3)
):
    """Decomposes all of the organic sediment within the marsh soil profile at a rate determined by depth."""

    compaction = np.zeros([B])
    Fd = 0

    # Decompose the marsh sediment
    for x in range(x_m, x_f):  # Loop through each marsh cell in the domain
        decomp = np.zeros([yr + 1])
        for tempyr in range(yr, 0, -1):  # Loop through each layer of sediment in each cell, starting at the most recently deposited layer of sediment at the surface
            # but the most recent layer is actually on the "bottom" of the array because row number increases as you go down in the array
            depth = elevation[yr, x] - elevation[tempyr, x]  # Depth of sediment layer below the surface
            if depth > mui:  # Maximum depth at which decomposition occurs
                decomp[tempyr] = 0
                break
            else:
                decomp[tempyr] = organic_dep_autoch[tempyr, x] * (mki * math.exp(-depth / mui))  # [g] Mass of organic material decomposed from a given "pocket" of sediment
                # decomp[tempyr] = organic_dep_autoch[tempyr][x] * (mki * math.exp(-depth / mui))  # [g] Mass of organic material decomposed from a given "pocket" of sediment
                # organic_dep_autoch[tempyr][x] -= decomp[tempyr]  # [g] Autochthanous organic material in a given "pocket" of sediment updated for deomposition
                organic_dep_autoch[tempyr, x] -= decomp[tempyr]  # [g] Autochthanous organic material in a given "pocket" of sediment updated for deomposition
        compaction[x] = np.sum(decomp) / 1000 / rhoo  # [m] Total compaction in a given cell is a result of the sum of all decomposition in that cell
        Fd += np.sum(decomp)  # [kg] Flux of organic matter out of the marsh due to decomposition

    # Update the elevation of only the most recent layer
    elevation[yr, x_m:x_f] -= compaction[x_m:x_f]

    return elevation[yr, x_m:x_f], compaction, Fd, organic_dep_autoch


class Marsh:
    """apply marsh dynamics

    Examples
    --------
    # >>> from cascade.marsh_dynamics import Marsh
    # >>> marsh = Marsh()
    # >>> marsh.update(INPUTS)
    """

    def __init__(
        self,
        RSLR,               # rate of SLR in m/yr
        C_e,                # SSC at marsh edge (kg/m3)
        OCb,                # organic content of uppermost layer of bay sediment (%), currently set to 0.05
        numiterations,      # tidal iterations, currently set to 500
        P,                  # tidal period, 12.5 hours for semi-diurnal tide, converted to seconds
        ws,                 # settling velocity (m/s), currently set to 5 x 10-5
        n_tidal_cycles,     # 365 * (24 / 12.5)  # number of tidal cycles per year
        Bmax,               # g/m2 (I think this is maximum biomass)
        Dmin,               # minimum depth below high water that marsh vegetation can grow (m)
        Dmax,               # maximum depth below high water that marsh vegetation can grow (m)
        rhoo,               # organic matter bulk density (kg/m3), set to 85
        rhos,               # sediment bulk density (kg/m3), set to 2000
        mui,                # depth in marsh where decomp goes to 0 (m)
        mki,                # decomp coefficient
        m_min,              # minimum depth that is considered marsh (- 0.3 dam MHW)
        m_max,              # maximum depth that is considered marsh (0 dam MHW)
        time_step_count,    # total model time steps
        alongshore_length,  # alongshore length of the barrier (500 m or 50 cells)
        tidal_amplitude,    # tidal amplitude (m)
    ):
        """

        """

        # initial variables
        self._RSLR = RSLR  # m/yr
        self._SSCb = C_e  # kg/m3
        self._OCb = OCb  # %
        self._numiterations = numiterations
        self._P = P  # seconds
        self._ws = ws  # m/s
        self._n_tides = n_tidal_cycles
        self._Bmax = Bmax  # g/m2
        self._Dmin = Dmin  # m
        self._Dmax = Dmax  # m
        self._rhoo = rhoo  # kg/m3
        self._rhos = rhos  # kg/m3
        self._mui = mui  # m
        self._mki = mki
        self._m_min = m_min
        self._m_max = m_max
        self._nt = time_step_count
        self._amp = tidal_amplitude  # m
        self._tr = tidal_amplitude * 2  # m

        # initialize other variables that do not have inputs
        self._inundation_time = P / numiterations  # seconds

        # initialize empty arrays to keep track of variables at each model time step
        self._marsh_elevation = [np.nan] * self._nt
        self._yearly_decomp_elev_TS = [np.nan] * alongshore_length
        self._yearly_organic_autoch_TS = [np.nan] * alongshore_length
        self._accretion_TS = [np.nan] * self._nt
        self._compaction_TS = [np.nan] * self._nt
        self._msl = np.linspace(
            1, self._nt, num=self._nt) * self._RSLR  # [m] Mean sea level over time relative to start

        # initialize arrays for decomp so the rows are the total model duration and columns are barrier width
        # NOTE: not sure what we do about changing barrier width throughout time, so I made these extra long and set
        # the values to nans or maybe I will change to bay elevation?
        initial_width = 500
        for c in range(alongshore_length):
            self._yearly_decomp_elev_TS[c] = np.ones([self._nt, initial_width]) * np.nan  # time is years and width is dam
            self._yearly_organic_autoch_TS[c] = np.zeros([self._nt, initial_width])  # zeros are better for this array bc nothing
            # happens when this value is zero

        # initialize the time index variable (changed to b3d_time_step in the update function)
        self._time_index = 0

    def update(self, interior_domain, model_year):

        self._time_index = model_year

        # bmft assumes the domain starts at the marsh edge and ends at the higher elevation "forest", so we need
        # to flip our domain so that 0 is the marsh edge (top is marsh edge, bottom is barrier edge)
        interior_domain = np.flip(interior_domain, axis=0)

        # we will need to convert the interior domain from dam MHW to m MSL
        interior_domain = interior_domain * 10  # convert to m
        interior_domain = interior_domain + self._msl[model_year] + self._amp + self._RSLR  # convert to MSL
        m_max_msl = (self._m_max * 10) + self._msl[model_year] + self._amp + self._RSLR
        m_min_msl = (self._m_min * 10) + self._msl[model_year] + self._amp + self._RSLR

        n_cols = np.shape(interior_domain)[1]
        barrier_width = np.shape(interior_domain)[0]

        # initialize accretion and compaction arrays
        self._accretion_TS[model_year] = np.zeros(np.shape(interior_domain))
        self._compaction_TS[model_year] = np.zeros(np.shape(interior_domain))

        for c in range(n_cols):
            # initial transect
            transect = interior_domain[:, c]
            if model_year == 1:
                self._yearly_decomp_elev_TS[c][0, 0:barrier_width] = copy.deepcopy(transect)  # set the first row to the initial transect elevation
            self._yearly_decomp_elev_TS[c][model_year, 0:barrier_width] = copy.deepcopy(transect)  # m MSL, this model year's transect
            # identify the marsh cells
            marsh_cells = np.where((transect <= m_max_msl) & (transect > m_min_msl))[0]  # if none, all cells are too high or too low to be marsh
            if len(marsh_cells) == 0:
                marsh_transect = []
            else:
                start_marsh_cell = np.min(marsh_cells)
                end_marsh_cell = np.max(marsh_cells)
                marsh_transect = transect[start_marsh_cell:end_marsh_cell+1]

            # if marsh_transect is empty, there are no marsh cells and we skip the calcs
            if len(marsh_transect) != 0:
                # marsh accretion
                marsh_transect, accretion, organic_autoch = evolvemarsh(
                    marshelevation=marsh_transect,
                    msl=self._msl[model_year],
                    C_e=self._SSCb,
                    OCb=self._OCb,
                    tr=self._tr,
                    numiterations=self._numiterations,
                    P=self._P,
                    ws=self._ws,
                    timestep=self._n_tides,
                    Bmax=self._Bmax,
                    Dmin=self._Dmin,
                    Dmax=self._Dmax,
                    rhoo=self._rhoo,
                    rhos=self._rhos,
                    plot=False,
                    )

                # store accretion values
                # currently oriented with marsh on top, dunes/ocean on bottom
                self._accretion_TS[model_year][start_marsh_cell:end_marsh_cell+1, c] = accretion

                # oriented marsh at column 0, barrier dunes at last column, newest elevation layer on bottom, olders on top
                # add the most recent marsh transect to the marsh layers
                self._yearly_decomp_elev_TS[c][model_year, start_marsh_cell:end_marsh_cell+1] = copy.deepcopy(marsh_transect)  # m MSL
                # previous time periods as well as marsh_transect updates
                self._yearly_organic_autoch_TS[c][model_year, start_marsh_cell:end_marsh_cell+1] = organic_autoch

                marsh_transect, compaction, Fd, organic_dep_autoch = decompose(
                        x_m=start_marsh_cell,
                        x_f=end_marsh_cell+1,
                        yr=model_year,
                        organic_dep_autoch=self._yearly_organic_autoch_TS[c],
                        elevation=self._yearly_decomp_elev_TS[c],
                        B=len(transect),
                        mui=self._mui,
                        mki=self._mki,
                        rhoo=self._rhoo,
                    )

                # store compaction values
                # currently oriented with marsh on top, dunes/ocean on bottom
                self._compaction_TS[model_year][:, c] = compaction  # compaction is the full transect, not just marsh

        # re-flip the domain so it is oriented correctly for b3d (marsh/bay on bottomw, barrier/ocean on top)
        interior_domain = np.flip(interior_domain, axis=0)
        # convert back to dam MHW
        interior_domain = interior_domain - (self._msl[model_year] + self._amp + self._RSLR)  # convert to MHW
        # in barrierbmft, Ian does not subtract RSLR and I think it is bc barrier3d lowers the elevation of the barrier
        # by RLSR, so it is already going to be accounted for?

        # keep in same orientation as b3d domain for comparison
        self._accretion_TS[model_year] = np.flip(self._accretion_TS[model_year], axis=0)
        self._compaction_TS[model_year] = np.flip(self._compaction_TS[model_year], axis=0)

        interior_domain = interior_domain / 10  # convert to dam
        self._marsh_elevation[model_year] = copy.deepcopy(interior_domain)  # dam MHW

        return interior_domain
