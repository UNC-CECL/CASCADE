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

def evolvemarsh(
        marshelevation,     # elevation domain, here it is ONLY marsh cells
        msl,                # msl in m
        C_e,                # SSC at marsh edge (kg/m3)
        OCb,                # organic content of uppermost layer of bay sediment (%), currently set to 0.05
        tr,                 # tidal range, m
        numiterations,      # tidal iterations, currently set to 500
        P,                  # tidal period, 12.5 hours for semi-diurnal tide, converted to seconds
        dt,                 # tidal dt ( a.k.a. inundation time) = (P/numiterations), seconds?
        ws,                 # settling velocity (m/s), currently set to 5 x 10-5
        timestep,           # 365 * (24 / 12.5)  # number of tidal cycles per year
        BMax,               # g/m2 (I think this is maximum biomass)
        Dmin,               # minimum depth below high water that marsh vegetation can grow (m)
        Dmax,               # maximum depth below high water that marsh vegetation can grow (m)
        rhoo,               # organic matter bulk density (kg/m3), set to 85
        rhos,               # sediment bulk density (kg/m3), set to 2000
):
    """Calculates biomass and mineral and organic deposition for each cell in the marsh as a function of flooding
    frequency; calculates the total flux of sediment onto the marsh from the bay."""

    L = len(marshelevation)
    time_submerged = np.zeros([numiterations, L])
    sedimentcycle = np.zeros([numiterations, L])
    depth = np.zeros([numiterations, L])
    C = np.zeros([L])

    # Loop through a tidal cycle to determine flooding duration for each point in the marsh
    # time_submerged is used to calculate floodfraction but that is never used itself,
    # BUT the depth variable is used to calculate sediment deposition, so we definitely need that still
    for i in range(1, numiterations):
        depth_iteration = 0.5 * tr * math.sin(2 * math.pi * ((i + 1) * dt / P)) + (msl - marshelevation[0: L + 1])  # [m] Depth at each position in the marsh
        depth[i, :] = depth_iteration  # Store depths in array
        temp = np.zeros([L])
        temp[depth_iteration > 0] = dt
        time_submerged[i, :] = temp  # Inundation for a single flood cycle recorded for each cell

    # -------------------------
    # Belowground Productivity

    # Creates a biomass curve (Mariotti & Carr, 2014) where peak biomass occurs at a depth halfway between the maximum depth for vegetation and the minimum (here, mean high water level)

    dm = msl + tr / 2 - marshelevation  # [m] Depth of the marsh surface below HWL at any given point

    bgb = np.zeros([L])  # Belowground biomass
    agb = np.zeros([L])  # Aboveground biomass
    organic_autoch = np.zeros([L])  # Autochthonous organic material

    for ii in range(L):
        if dm[ii] > Dmax:  # If depth is above vegetation maximum, there is no production
            agb[ii] = 0  # [g] Mudflat
            bgb[ii] = 0  # [g] Mudflat
            organic_autoch[ii] = 0  # [g] No autochthonous material stored in the soil; we do not multiply by lingin content (as in earlier version of the model) because we subtract mass due to decomposition in the 'decompose' function
        elif dm[ii] <= Dmin:  # If depth is below vegetation minimum, there is very little belowground productivity
            if ii > 6000:
                agb[ii] = 100  # [g] Forest aboveground - constant for now, should depend on elevation
                bgb[ii] = 0.00001  # [g] Forest
                organic_autoch[ii] = bgb[ii]  # [g] Forest organic matter stored in soil
            else:
                agb[ii] = 1 * (BMax * (dm[ii] - Dmax) * (dm[ii] - Dmin) / (0.25 * (-Dmin - Dmax) * (Dmax - 3 * Dmin)))
                bgb[ii] = 0.1  # [g] "Forest" - really high marsh
                organic_autoch[ii] = bgb[ii]  # [g] Forest organic matter stored in soil
        else:
            agb[ii] = 1 * (BMax * (dm[ii] - Dmax) * (dm[ii] - Dmin) / (0.25 * (-Dmin - Dmax) * (Dmax - 3 * Dmin)))  # [g] Marsh
            bgb[ii] = BMax * (dm[ii] - Dmax) * (dm[ii] - Dmin) / (0.25 * (-Dmin - Dmax) * (Dmax - 3 * Dmin))  # [g] Marsh
            organic_autoch[ii] = bgb[ii]  # [g] As mentioned above, no longer multiply by lingin content

    # -------------------------
    # Mineral Deposition
    coeff = -0.002  # Coefficient of -0.0031 is a fitted parameter for realistic marsh topography
    distance = 0  # [m] Initialize, distance from marsh edge
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

    floodfraction = np.sum(time_submerged, axis=0) / P  # Portion of the tidal cycle that each point is submerged

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

    # Total flux of sediment onto the marsh from the bay
    Fm_min = np.sum(mineral) / 1000  # [kg/yr] Flux of mineral sediment from the bay
    Fm_org = np.sum(organic_alloch) / 1000  # [kg/yr] Flux of organic sediment from the bay

    # Calculate thickness of new sediment (mineral+organic) based off LOI and its effect on density
    loi = (organic_autoch + organic_alloch) / (mineral + organic_autoch + organic_alloch)
    density = 1 / ((loi / rhoo) + ((1 - loi) / rhos)) * 1000  # [g/m3] Bulk density is calculated according to Morris et al. (2016)
    density[np.isnan(density)] = 1  # If there is zero deposition, loi calculation will divide by zero and make density nan. Set density equal to 1 in this case, so that accretion is zero, instead of nan.

    accretion = (mineral + organic_autoch + organic_alloch) / density  # [m] accretion in a given year

    # Update elevation
    marshelevation += accretion

    return marshelevation, organic_autoch, organic_alloch, mineral, Fm_min, Fm_org, bgb, accretion, agb


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
    for x in range(x_m, x_f):  # Loop through each marsh and upland cell in the domain (really just the marsh since it stops at the forest)
        decomp = np.zeros([yr + 1])
        for tempyr in range(yr, 0, -1):  # Loop through each pocket of sediment in each cell, starting at the most recently deposited packet of sediment at the surface
            depth = elevation[yr, x] - elevation[tempyr, x]  # Depth of sediment pocket below the surface
            if depth > mui:  # Maximum depth at which decomposition occurs
                decomp[tempyr] = 0
                break
            else:
                decomp[tempyr] = organic_dep_autoch[tempyr, x] * (mki * math.exp(-depth / mui))  # [g] Mass of organic material decomposed from a given "pocket" of sediment
                organic_dep_autoch[tempyr, x] -= decomp[tempyr]  # [g] Autochthanous organic material in a given "pocket" of sediment updated for deomposition
        compaction[x] = np.sum(decomp) / 1000 / rhoo  # [m] Total compaction in a given cell is a result of the sum of all decomposition in that cell
        Fd += np.sum(decomp)  # [kg] Flux of organic matter out of the marsh due to decomposition

    return compaction, Fd, organic_dep_autoch


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
        interior_elev,      # interior elevation domain
        msl,                # msl in m
        C_e,                # SSC at marsh edge (kg/m3)
        OCb,                # organic content of uppermost layer of bay sediment (%), currently set to 0.05
        tr,                 # tidal range, m
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
    ):
        """

        """

        # initial variables
        self._elevation = interior_elev,  # dam
        self._msl = msl,  # m
        self._SSC = C_e,  # kg/m3
        self._OCb = OCb,  # %
        self._tr = tr,  # m
        self._numiterations = numiterations,
        self._P = P,  # seconds
        self._ws = ws,  # m/s
        self._n_tides = n_tidal_cycles,
        self._BMax = Bmax,  # g/m2
        self._Dmin = Dmin,  # m
        self._Dmax = Dmax,  # m
        self._rhoo = rhoo,  # kg/m3
        self._rhos = rhos,  # kg/m3
        self._mui = mui,  # m
        self._mki = mki,

        # initialize other variables that do not have inputs
        self._inundation_time = P / numiterations  # seconds
        self._organic_dep_autoch = np.zeros(
            np.shape(self._elevation), dtype=object)  # grams, this is an output from growing the marsh,
            # initializing at 0 here


        # initialize the time index variable (changed to b3d_time_step in the update function)
        self._time_index = 0


    def update(self, b3d):

        self._time_index = b3d.time_index

        # marsh deposition


        # marsh decomposition


        # update Barrier3D variables

        # interior domain: remove all rows of bay without any deposition from
        # the domain we check the first row rather than the last row (which is
        # used in B3D) because we have not flipped the domain yet
        check = 1
        while check == 1:
            if all(
                    x <= -self._bay_depth for x in post_outwash_interior_domain[0, :]
                    ):  # bay_depth = 0.3
                post_outwash_interior_domain = np.delete(
                    post_outwash_interior_domain, 0, axis=0
                    )
            else:
                check = 0

        new_ave_interior_height = np.average(
            post_outwash_interior_domain[
                post_outwash_interior_domain >= b3d.SL
                ]  # all in dam MHW
            )
        b3d.h_b_TS[-1] = new_ave_interior_height
        # flip the Interior and Dune Domains to put back into Barrier3D
        b3d.InteriorDomain = np.flip(post_outwash_interior_domain)
        b3d.DomainTS[self._time_index - 1] = np.flip(post_outwash_interior_domain)


        # final update to outwash and domain variables ---
        # recalculate and save DomainWidth and InteriorWidth
        _, _, new_ave_interior_width = b3d.FindWidths(b3d.InteriorDomain, b3d.SL)
        b3d.InteriorWidth_AvgTS[-1] = new_ave_interior_width

        # replace value of x_b_TS with new InteriorWidth_Avg
        b3d.x_b_TS[-1] = b3d.x_s_TS[-1] + new_ave_interior_width
