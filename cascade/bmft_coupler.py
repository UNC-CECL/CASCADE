from bmftc import Bmftc
import numpy as np
import math

"""BMFT Coupler

This module couples barrier3d with BMFT to allow for basic marsh dynamics

References
----------

Notes
---------

"""

###############################################################################
# initialize marsh dynamic modules
###############################################################################


def B3d_PyBMFT_equal(barrier3d, ny, bmftc):
    for iB3D in range(ny):
        barrier3d[iB3D]._TMAX = bmftc[iB3D].dur + 1  # [yrs] Duration of simulation
        barrier3d[iB3D]._RSLR = (
            np.ones([len(barrier3d[iB3D].RSLR) + 1]) * (bmftc[iB3D].RSLRi / 1000) / 10
        )  # [m/yr] Relative sea-level rise rate, converted units
        barrier3d[iB3D]._BayDepth = (
            bmftc[iB3D].Bay_depth[bmftc[iB3D].startyear - 1] / 10
        )  # [yrs] Initial depth of bay
    return barrier3d


class BMFTCoupler:
    def __init__(self, ny, nt, barrier3d, name):

        """The AlongshoreCoupler module.

        Parameters
        ----------
        name: string, mandatory
            Name of simulation
        ny: int, mandatory
            The number of alongshore Barrier3D domains for simulation in BRIE
        nt: int, mandatory
            Number of time steps.
        barrier3d: list, mandatory
            Barrier3d subdomains created in earlier functions
        """
        self._ny = ny
        self._nt = nt + 1

        # PyBMFT Variables
        self._bmftc = []
        self._BMFTC_Break = False  # Initialize BMFTC break variable
        # Initialize blank PyBMFT variables as lists
        self._name = []
        self._x_b_TS = []
        self._LandscapeTypeWidth_TS = []
        self._bay_overwash_carryover = (
            []
        )  # [m^3] Volume of overwash deposition into back-barrier bay from previous year that did not fill new cell up to sea level; is added to overwash bay dep in following year
        self._x_s_offset = (
            []
        )  # Initial location of B in PyBMFT-C relative to x_s_initial in Barrier3D
        self._cumul_len_change = []
        self._delta_fetch_TS = []
        self._OWspread = (
            0  # [%] Percentage of overwash past marsh edge that is spread across bay
        )
        self._b3d_transect_difs = []

        # initialize PyBMFT models (number set by brie ny above)
        for iB3D in range(self._ny):
            self._bmftc.append(
                Bmftc(
                    name="back-barrier",
                    time_step_count=self._nt,
                    # relative_sea_level_rise=barrier3d[iB3D]._RSLR[1] * 1000,
                    relative_sea_level_rise=12,
                    reference_concentration=60,
                    slope_upland=0.005,
                    bay_fetch_initial=5000,
                    forest_width_initial_fixed=False,
                    forest_width_initial=5000,  # 5000 accomodates 250 yrs at R=15 and S=0.001
                    wind_speed=6,
                    forest_on=False,
                    filename_equilbaydepth="/Users/ceclmac/PycharmProjects/PyBMFT-C/Input/PyBMFT-C/Equilibrium Bay Depth.mat",
                    filename_marshspinup="/Users/ceclmac/PycharmProjects/PyBMFT-C/Input/PyBMFT-C/MarshStrat_all_RSLR1_CO50_width500.mat",
                    marsh_width_initial=500,
                )
            )

        for iB3D in range(self._ny):
            # ===========================================
            # Add initial barrier topography from Barrier3D to initial "forest" (i.e., subaerial) portion of PyBMFT-C transect
            b3d_transect = (
                np.mean(barrier3d[iB3D].InteriorDomain, axis=1) * 10
            )  # Take average across alongshore dimension, convert to m (vertical dimension)
            x = np.linspace(1, len(b3d_transect) * 10, num=len(b3d_transect) * 10)
            xp = np.linspace(1, len(b3d_transect), num=len(b3d_transect)) * 10
            xp = xp - 5
            b3d_transect = np.interp(
                x, xp, b3d_transect
            )  # Interpolate from dam to m (horizontal dimension)
            x_f = (
                np.where(b3d_transect < (barrier3d[iB3D].SL * 10))[0][0] - 1
            )  # [m] Distance of first interior (subaerial) cell from B3D ocean shoreline (excluding dunes/beach)
            b3d_transect = b3d_transect[:x_f]
            b3d_transect = np.flip(b3d_transect)
            b3d_transect = (
                b3d_transect
                + self._bmftc[iB3D].msl[self._bmftc[iB3D].startyear - 1]
                + self._bmftc[iB3D].amp
                + (self._bmftc[iB3D].RSLRi / 1000)
            )  # Convert vertical datums

            # Adjust size of Barrier3D topo to fit PyBMFT-C "forest" section
            BB_forest_len = len(
                self._bmftc[iB3D].elevation[
                    self._bmftc[iB3D].startyear, self._bmftc[iB3D].x_f :
                ]
            )
            if len(b3d_transect) > BB_forest_len:
                subtract = len(b3d_transect) - BB_forest_len
                b3d_transect = b3d_transect[:-subtract]
            elif len(b3d_transect) < BB_forest_len:
                add = np.ones([BB_forest_len - len(b3d_transect)]) * (
                    self._bmftc[iB3D].msl[self._bmftc[iB3D].startyear]
                    + self._bmftc[iB3D].amp
                )
                b3d_transect = np.append(b3d_transect, add)

            # Replace initial subaerial elevation in PyBMFT-C with Barrier3D initial barrier elevation
            self._bmftc[iB3D].elevation[
                self._bmftc[iB3D].startyear - 1, self._bmftc[iB3D].x_f :
            ] = b3d_transect  # Replace!

            # ===========================================
            # Populate blank PyBMFT list variables
            self._name.append(name)
            self._x_b_TS.append(np.zeros([self._bmftc[iB3D].dur]))
            self._LandscapeTypeWidth_TS.append(np.zeros([self._bmftc[iB3D].dur, 4]))
            self._bay_overwash_carryover.append(
                0
            )  # [m^3] Volume of overwash deposition into back-barrier bay from previous year that did not fill new cell up to sea level; is added to overwash bay dep in following year
            initial_subaerial_width = self._bmftc[iB3D].B - self._bmftc[iB3D].x_f
            self._x_s_offset.append(
                initial_subaerial_width - (barrier3d[iB3D].InteriorWidth_AvgTS[-1] * 10)
            )  # Initial location of B in PyBMFT-C relative to x_s_initial in Barrier3D
            self._cumul_len_change.append([0])
            # self._OWspread.append(0)  # [%] Percentage of overwash past marsh edge that is spread across bay
            self._delta_fetch_TS.append([])
            # self._OWspread.append([0])

            B3d_PyBMFT_equal(ny=self._ny, barrier3d=barrier3d, bmftc=self._bmftc)

    ###############################################################################
    # Backbarrier marsh module
    ###############################################################################
    def updateMarsh(self, ny, time_step, barrier3d):
        # Get new Barrier3D geometry
        time_step = time_step - 2
        # Just barrier3d look at roadway manager class

        # ~~~~~~~~~~~~~~ PyBMFT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Run PyBMFT module to represent marsh growth and erosion from the back-bay
        for iB3D in range(ny):
            """Update BarrierBMFT by one time step"""
            # ===================================================================================================================================================================================================================================
            # ===================================================================================================================================================================================================================================
            # Advance PyBMFT-C back-barrier marshes
            # Update RSLR rate from Barrier3D
            self._bmftc[iB3D]._RSLRi = barrier3d[iB3D].RSLR[time_step] * 10
            self._bmftc[iB3D]._RSLR = (
                self._bmftc[iB3D].RSLRi * 10 ** (-3) / (3600 * 24 * 365)
            )  # Convert from mm/yr to m/s
            self._bmftc[iB3D].update()

            # Check if marsh has completely drowned or basin is completely full
            if self._bmftc[iB3D].drown_break == 1:
                self._bmftc[iB3D]._dur = time_step
                self._bmftc[iB3D]._endyear = self._bmftc[iB3D].startyear + time_step
                self._BMFTC_Break = True
                print(
                    "PyBMFT-C Simulation Break: marsh has completely drowned or basin is completely full"
                )
                return  # If so, end simulation
            # ===================================================================================================================================================================================================================================
            # ===================================================================================================================================================================================================================================
            # Update fetch and marsh point locations from PyBMFT-C bay erosion/deposition processes

            # Calculate change in fetch from erosion of both marshes
            delta_fetch_BB = (
                self._bmftc[iB3D].bfo
                - self._bmftc[iB3D].fetch[self._bmftc[iB3D].startyear + time_step - 1]
            )  # [m] Back-barrier marsh

            self._delta_fetch_TS[iB3D].append(delta_fetch_BB)

            # Determine change in x_b location
            self._x_b_TS[iB3D][time_step] = self._bmftc[iB3D].x_b  # Save to array

            # ===================================================================================================================================================================================================================================
            # ===================================================================================================================================================================================================================================
            # Adjust bay depth in Barrier3D according to depth calculated in PyBMFT-C
            barrier3d[iB3D]._BayDepth = np.mean([self._bmftc[iB3D].db]) / 10

            # ===================================================================================================================================================================================================================================
            # ===================================================================================================================================================================================================================================
            # Add marsh from PyBMFT-C to Barrier3D
            previous_elevation = barrier3d[iB3D].DomainTS[time_step]

            # Extract and convert marsh elevation from PyBMFT-C
            marsh_transect = self._bmftc[iB3D].elevation[
                self._bmftc[iB3D].startyear + time_step,
                self._bmftc[iB3D].x_m : self._bmftc[iB3D].x_f + 1,
            ]  # Marsh elevation from PyBMFT-C
            if len(marsh_transect) >= 1:
                len_marsh_transect = 10 * (
                    (len(marsh_transect) + 5) // 10
                )  # Cross-shore length of marsh rounded to nearest dam
                self._cumul_len_change[iB3D].append(
                    self._cumul_len_change[iB3D][-1]
                    + (len_marsh_transect - len(marsh_transect))
                )
                x = np.linspace(
                    1, len(marsh_transect) / 10, num=int((len_marsh_transect / 10))
                )
                xp = np.linspace(
                    1, len(marsh_transect) / 10, num=int(len(marsh_transect))
                )
                marsh_transect = np.interp(
                    x, xp, marsh_transect
                )  # Interpolate marsh elevation from m to dam in the horizontal dimension
                marsh_transect = marsh_transect - (
                    self._bmftc[iB3D].msl[self._bmftc[iB3D].startyear + time_step]
                    + self._bmftc[iB3D].amp
                )  # Make marsh elevation relative to MHW datum
                marsh_transect = (
                    marsh_transect / 10
                )  # Convert from m to dam in the vertial dimension
                marsh_transect = np.flip(marsh_transect)

            StartDomainWidth = np.shape(barrier3d[iB3D].InteriorDomain)[
                0
            ]  # Width of interior domain from last time step

            # Find barrier interior widths for each dam alongshore
            InteriorWidth = [0] * barrier3d[iB3D].BarrierLength
            for bl in range(barrier3d[iB3D].BarrierLength):
                width = next(
                    (
                        index
                        for index, value in enumerate(
                            barrier3d[iB3D].InteriorDomain[:, bl]
                        )
                        if value <= barrier3d[iB3D].SL
                    ),
                    StartDomainWidth,
                )
                width = width - 1
                if width < 0:
                    width = 0
                InteriorWidth[bl] = width

            # Update Barrier3D Domain Sizes
            Target_width_barriermarsh = (
                self._bmftc[iB3D].B - self._bmftc[iB3D].x_m - self._x_s_offset[iB3D]
            )  # [m] Target width of barrier-marsh
            Target_width_barriermarsh = math.ceil(
                Target_width_barriermarsh / 10
            )  # [dam]
            addRows = (
                Target_width_barriermarsh - StartDomainWidth + 1
            )  # Number of rows to add (if positive) or subtract (if negative) from Barrier3D domain

            if addRows > 0:
                # Update interior domain size
                Marsh_Addition = (
                    np.ones([addRows, barrier3d[iB3D].BarrierLength])
                    * -barrier3d[iB3D]._BayDepth
                )
                Zero_Addition = np.zeros([addRows, barrier3d[iB3D].BarrierLength])
                NewDomain = np.vstack([barrier3d[iB3D].InteriorDomain, Marsh_Addition])
                previous_elevation = np.append(
                    np.mean(barrier3d[iB3D].DomainTS[time_step], axis=1),
                    np.mean(Marsh_Addition, axis=1),
                )
                # Update size of shrub domains, too
                barrier3d[iB3D]._ShrubDomainFemale = np.vstack(
                    [barrier3d[iB3D]._ShrubDomainFemale, Zero_Addition]
                )
                barrier3d[iB3D]._ShrubDomainMale = np.vstack(
                    [barrier3d[iB3D]._ShrubDomainMale, Zero_Addition]
                )
                barrier3d[iB3D]._ShrubDomainDead = np.vstack(
                    [barrier3d[iB3D]._ShrubDomainDead, Zero_Addition]
                )
                barrier3d[iB3D]._ShrubPercentCover = np.vstack(
                    [barrier3d[iB3D]._ShrubPercentCover, Zero_Addition]
                )
                barrier3d[iB3D]._DeadPercentCover = np.vstack(
                    [barrier3d[iB3D]._DeadPercentCover, Zero_Addition]
                )
                barrier3d[iB3D]._BurialDomain = np.vstack(
                    [barrier3d[iB3D]._BurialDomain, Zero_Addition]
                )
                barrier3d[iB3D]._ShrubDomainAll = (
                    barrier3d[iB3D]._ShrubDomainFemale
                    + barrier3d[iB3D]._ShrubDomainMale
                )
            elif addRows < 0:
                # Update interior domain size
                NewDomain = barrier3d[iB3D].InteriorDomain[:addRows, :]
                previous_elevation = previous_elevation[:addRows]
                # Update size of shrub domains, too
                barrier3d[iB3D]._ShrubDomainFemale = barrier3d[iB3D]._ShrubDomainFemale[
                    :addRows, :
                ]
                barrier3d[iB3D]._ShrubDomainMale = barrier3d[iB3D]._ShrubDomainMale[
                    :addRows, :
                ]
                barrier3d[iB3D]._ShrubDomainDead = barrier3d[iB3D]._ShrubDomainDead[
                    :addRows, :
                ]
                barrier3d[iB3D]._ShrubPercentCover = barrier3d[iB3D]._ShrubPercentCover[
                    :addRows, :
                ]
                barrier3d[iB3D]._DeadPercentCover = barrier3d[iB3D]._DeadPercentCover[
                    :addRows, :
                ]
                barrier3d[iB3D]._BurialDomain = barrier3d[iB3D]._BurialDomain[
                    :addRows, :
                ]
                barrier3d[iB3D]._ShrubDomainAll = (
                    barrier3d[iB3D]._ShrubDomainFemale
                    + barrier3d[iB3D]._ShrubDomainMale
                )
            else:
                NewDomain = barrier3d[iB3D].InteriorDomain  # Domains stay same size

            if len(marsh_transect) >= 1:
                # Update Marsh In Barrier3D
                x_marsh = (
                    Target_width_barriermarsh + 1
                )  # [dam] Cross-shore location of marsh edge relative to interior domain
                for w in range(barrier3d[iB3D].BarrierLength):
                    width_diff = x_marsh - (InteriorWidth[w] + len(marsh_transect))
                    if width_diff < 0:
                        MarshTransect = marsh_transect[: -int(abs(width_diff))]  # [dam]
                    elif width_diff > 0:
                        add = (
                            np.ones([int(abs(width_diff))]) * marsh_transect[-1]
                        )  # Set additional marsh cells to elevation of last marsh
                        MarshTransect = np.append(marsh_transect, add)  # [dam]
                    else:
                        MarshTransect = marsh_transect  # [dam]

                    InteriorTransect = NewDomain[: InteriorWidth[w], w]  # [dam]
                    BarrierMarshTransect = np.append(
                        InteriorTransect, MarshTransect
                    )  # Combine interior and marsh

                    NewDomain[: len(BarrierMarshTransect), w] = BarrierMarshTransect
                    NewDomain[len(BarrierMarshTransect) :, w] = (
                        barrier3d[iB3D].SL - np.mean([self._bmftc[iB3D].db])
                    ) / 10

            barrier3d[iB3D].InteriorDomain = NewDomain
            # ===================================================================================================================================================================================================================================
            # Update PyBMFT-C transect elevation based on Barrier3D elevation change
            old_elev = previous_elevation * 10

            shoreline_change = barrier3d[iB3D].x_s_TS[-1] - barrier3d[iB3D].x_s_TS[-2]
            self._x_s_offset[iB3D] = self._x_s_offset[iB3D] + (shoreline_change * 10)
            if time_step < 3:
                start_b3d = (
                    np.mean(NewDomain, axis=1) * 10
                )  # Barrier3D domain before update, averaged across alongshore dimension, converted to m (vertical dimension)
                start_b3d = start_b3d + (barrier3d[iB3D].RSLR[time_step] * 10)
            else:
                start_b3d = (
                    np.mean(barrier3d[iB3D].DomainTS[time_step - 1], axis=1) * 10
                )

            if time_step < 3:
                end_b3d = (
                    np.mean(barrier3d[iB3D].InteriorDomain, axis=1) * 10
                )  # Barrier3D domain after update, averaged across alongshore dimension, converted to m (vertical dimension)
            else:
                end_b3d = np.mean(barrier3d[iB3D].DomainTS[time_step], axis=1) * 10

            # Update start domain size to match end domain
            sc_b3d = int(
                barrier3d[iB3D].ShorelineChangeTS[time_step] # NEW
                #barrier3d[iB3D].ShorelineChangeTS[-1] # OLD
            )  # Shoreline change [dam] from Barrier3D model update (this timestep)
            if sc_b3d < 0:  # Shoreline erosion
                start_b3d = start_b3d[abs(sc_b3d) :]  # Trim off front
            elif sc_b3d > 0:  # Shoreline progradation
                add = np.zeros([sc_b3d])
                start_b3d = np.append(add, start_b3d)  # Add zeros to front

            if len(start_b3d) < len(end_b3d):
                #add = np.ones([len(end_b3d) - len(start_b3d)]) * np.mean([self._bmftc[iB3D].db]) * -1  # Add bay cells OLD
                #add = np.ones([len(end_b3d) - len(start_b3d)]) * end_b3d[-1]  # Add marsh cells OLD
                add = end_b3d[-(len(end_b3d) - len(start_b3d)) :] + (
                    barrier3d[iB3D].RSLR[time_step] * 10
                ) # NEW
                start_b3d = np.append(start_b3d, add)
            elif len(start_b3d) > len(end_b3d):
                subtract = len(end_b3d) - len(start_b3d)
                start_b3d = start_b3d[:subtract]

            # Calculate change in elevation from Barrier3D update
            '''
            Find total amount of movement
            '''
            # Calculate the amount of total shoreline change
            if time_step>2:
                Current_B3D_Shoreline_Change = int(sum(barrier3d[iB3D].ShorelineChangeTS[:time_step]))
                Previous_B3D_Shoreline_Change = int(sum(barrier3d[iB3D].ShorelineChangeTS[:time_step]))

                # Add amount of shoreline erosion to the front of islands
                Current_Offset_B3D_Elevations = np.zeros(int(abs(Current_B3D_Shoreline_Change)))
                Previous_Offset_B3D_Elevations = np.zeros(int(abs(Previous_B3D_Shoreline_Change)))

                B3D_Elev_Current = np.mean(barrier3d[iB3D].DomainTS[time_step], axis=1)*10
                B3D_Elev_Previous = np.mean(barrier3d[iB3D].DomainTS[time_step-1], axis=1)*10

                # Convert to list
                Current_Offset_B3D_Elevations = Current_Offset_B3D_Elevations.tolist()
                Previous_Offset_B3D_Elevations = Previous_Offset_B3D_Elevations.tolist()
                B3D_Elev_Current = B3D_Elev_Current.tolist()
                B3D_Elev_Previous = B3D_Elev_Previous.tolist()

                # Combine
                Current_Offset_B3D_Elevations.extend(B3D_Elev_Current)
                Previous_Offset_B3D_Elevations.extend(B3D_Elev_Previous)

                if len(Current_Offset_B3D_Elevations) > len(Previous_Offset_B3D_Elevations):
                    Dif_Len = len(Current_Offset_B3D_Elevations) - len(Previous_Offset_B3D_Elevations)
                    Zeros_to_back = np.zeros(Dif_Len)
                    Zeros_to_back = Zeros_to_back.tolist()
                    Previous_Offset_B3D_Elevations.extend(Zeros_to_back)
                elif len(Current_Offset_B3D_Elevations) < len(Previous_Offset_B3D_Elevations):
                    Dif_Len = len(Previous_Offset_B3D_Elevations) - len(Current_Offset_B3D_Elevations)
                    Zeros_to_back = np.zeros(Dif_Len)
                    Zeros_to_back = Zeros_to_back.tolist()
                    Current_Offset_B3D_Elevations.extend(Zeros_to_back)

                # Return to NP.array for math
                Current_Offset_B3D_Elevations = np.array(Current_Offset_B3D_Elevations)
                Previous_Offset_B3D_Elevations = np.array(Previous_Offset_B3D_Elevations)




            start_b3d = start_b3d #- (
                #barrier3d[iB3D].RSLR[time_step] * 10
            #)  # Offset sea-level rise from Barrier3D so that it isn't counted twice (i.e. RSLR already taken into account in PyBMFT-C)
            end_b3d = end_b3d # + (barrier3d[iB3D].RSLR[time_step] * 10)  # Offset sea-level rise from Barrier3D so that it isn't counted twice (i.e. RSLR already taken into account in PyBMFT-C)
            elevation_change_b3d = (
                end_b3d - start_b3d
            )
            #print('')
            #print('Method 1 '+str(elevation_change_b3d))

            # Checking methods are the same
            if time_step > 2:
                dif_elev_change = Current_Offset_B3D_Elevations - Previous_Offset_B3D_Elevations
                #Any_Actual_Diffence = elevation_change_b3d - dif_elev_change
                #print('')
                #print('Difference method 2 '+str(dif_elev_change))

            # Change in elevation across transect after Barrier3d update; [dam] horizontal dimension, [m] vertical dimentsion
            #print('Elevation Change is '+str(elevation_change_b3d))
            self._b3d_transect_difs.append(elevation_change_b3d)
            # Interpolate from dam to m (horizontal dimension)
            x = np.linspace(
                1, len(elevation_change_b3d) * 10, num=len(elevation_change_b3d) * 10
            )
            xp = (
                np.linspace(1, len(elevation_change_b3d), num=len(elevation_change_b3d))
                * 10
            )
            xp = xp - 5
            elevation_change_b3d = np.interp(x, xp, elevation_change_b3d)
            off = int(
                abs(math.floor(self._x_s_offset[iB3D]))
            )  # [m] Offset of barrier shoreline and B

            # Incorporate elevation change from Barrier3D into back-barrier instance of PyBMFT-C
            if int(math.floor(self._x_s_offset[iB3D])) < 0:
                elevation_change_b3d = np.flip(
                    elevation_change_b3d[off:]
                )  # Flip orientation
                marsh_barrier_width = self._bmftc[iB3D].B - self._bmftc[iB3D].x_m
                x_m_change = abs(
                    math.floor(len(elevation_change_b3d) - marsh_barrier_width)
                )  # Location of marsh edge within elevation_change_b3d

                # Add subaerial elevation change
                self._bmftc[iB3D].elevation[
                    self._bmftc[iB3D].startyear + time_step,
                    -len(elevation_change_b3d[x_m_change:]) :,
                ] += elevation_change_b3d[
                    x_m_change:
                ]  # Store mass of overwash mineral sediment deposited across transect
                self._bmftc[iB3D].mineral_dep[
                    self._bmftc[iB3D].startyear + time_step,
                    -len(elevation_change_b3d[x_m_change:]) :,
                ] += (
                    elevation_change_b3d[x_m_change:] * self._bmftc[iB3D].rhos * 1000
                )  # [g] Mass of pure mineral sediment deposited by overwash

                # Determine volume of sed deposited past initial marsh edge and into bay
                sum_marsh_dep = np.sum(elevation_change_b3d[:x_m_change]) * (
                    1 - self._OWspread
                )  # [m^3] Volume of overwash deposition landward of marsh edge deposited as marsh
                sum_bay_dep = (
                    np.sum(elevation_change_b3d[:x_m_change]) * self._OWspread
                )  # [m^3] Volume of overwash deposition landward of marsh edge deposited across bay bottom
                self._bmftc[iB3D]._Fow_min = max(
                    0, sum_bay_dep * self._bmftc[iB3D].rhos
                )  # [kg/yr] Overwash deposition into bay, volume converted to mass

                # Add volume of carryover from last time step
                sum_marsh_dep += (
                    self._bay_overwash_carryover
                )  # [m^3] Bay deposition from previous time step that wasn't enough to fully fill bay cell up to sea level

                # Calculate height of deposition needed to bring bay bottom up to avg marsh elevation
                new_marsh_height = self._bmftc[iB3D].db
                # Determine distance of marsh progradation from overwash deposition
                progradation_actual = (
                    sum_marsh_dep / new_marsh_height
                )  # [m] Amount of marsh progradation, in which all overwash dep in bay fills first bay cell, then second, and so on until no more sediment. Assumes overwash is not spread out over bay.

                progradation = int(
                    max(math.floor(progradation_actual[iB3D]), 0)
                )  # Round to nearest FULL meter
                self._bay_overwash_carryover = (
                    progradation_actual - progradation
                ) * new_marsh_height  # Save leftover volume of sediment to be added to sum_bay_dep in following time step

                if progradation > 0:
                    # Add subaqueous elevation change
                    self._bmftc[iB3D].elevation[
                        self._bmftc[iB3D].startyear + time_step,
                        self._bmftc[iB3D].x_m - progradation : self._bmftc[iB3D].x_m,
                    ] += new_marsh_height
                    # Store mass of overwash mineral sediment deposited across transect
                    self._bmftc[iB3D].mineral_dep[
                        self._bmftc[iB3D].startyear + time_step,
                        self._bmftc[iB3D].x_m - progradation : self._bmftc[iB3D].x_m,
                    ] += (
                        new_marsh_height * self._bmftc[iB3D].rhos * 1000
                    )  # [g] Mass of pure mineral sediment deposited by overwash

                # Spread overwash bay flux evenly across bay bottom
                bay_accrete = sum_bay_dep / (
                    self._bmftc[iB3D].bfo - progradation
                )  # [m] Vertical accretion of bay bottom from overwash deposition in bay
                self._bmftc[iB3D]._db = (
                    self._bmftc[iB3D].db + bay_accrete
                )  # Update bay depth

            elif int(math.floor(self._x_s_offset[iB3D])) > 0:
                elevation_change_b3d = np.flip(elevation_change_b3d)
                marsh_barrier_width = self._bmftc[iB3D].B - self._bmftc[iB3D].x_m
                x_m_change = abs(
                    math.floor(len(elevation_change_b3d) - (marsh_barrier_width - off))
                )  # Location of marsh edge within elevation_change_b3d

                # Add subaerial elevation change
                self._bmftc[iB3D].elevation[
                    self._bmftc[iB3D].startyear + time_step,
                    self._bmftc[iB3D].B
                    - off
                    - len(elevation_change_b3d[x_m_change:]) : self._bmftc[iB3D].B
                    - off,
                ] += elevation_change_b3d[x_m_change:]
                # Store mass of overwash mineral sediment deposited across transect
                self._bmftc[iB3D].mineral_dep[
                    self._bmftc[iB3D].startyear + time_step,
                    self._bmftc[iB3D].B
                    - off
                    - len(elevation_change_b3d[x_m_change:]) : self._bmftc[iB3D].B
                    - off,
                ] += (
                    elevation_change_b3d[x_m_change:] * self._bmftc[iB3D].rhos * 1000
                )  # [g] Mass of pure mineral sediment deposited by overwash

                # Determine volume of sed deposited past initial marsh edge and into bay
                sum_marsh_dep = np.sum(elevation_change_b3d[:x_m_change]) * (
                    1 - self._OWspread
                )  # [m^3] Volume of overwash deposition landward of marsh edge deposited as marsh
                sum_bay_dep = (
                    np.sum(elevation_change_b3d[:x_m_change]) * self._OWspread
                )  # [m^3] Volume of overwash deposition landward of marsh edge deposited across bay bottom
                self._bmftc[iB3D]._Fow_min = max(
                    0, sum_bay_dep * self._bmftc[iB3D].rhos
                )  # [kg/yr] Overwash deposition into bay, volume converted to mass

                # Add volume of carryover from last time step
                sum_marsh_dep += (
                    self._bay_overwash_carryover
                )  # [m^3] Bay deposition from previous time step that wasn't enough to fully fill bay cell up to sea level

                # Calculate height of deposition needed to bring bay bottom up to avg marsh elevation
                new_marsh_height = self._bmftc[iB3D].db

                # Determine distance of marsh progradation from overwash deposition
                progradation_actual = (
                    sum_marsh_dep / new_marsh_height
                )  # [m] Amount of marsh progradation, in which all overwash dep in bay fills first bay cell, then second, and so on until no more sediment. Assumes overwash is not spread out over bay.
                progradation = int(
                    max(math.floor(progradation_actual[iB3D]), 0)
                )  # Round to nearest FULL meter
                self._bay_overwash_carryover = (
                    progradation_actual - progradation
                ) * new_marsh_height  # Save leftover volume of sediment to be added to sum_bay_dep in following time step

                if progradation > 0:
                    # Add subaqueous elevation change
                    self._bmftc[iB3D].elevation[
                        self._bmftc[iB3D].startyear + time_step,
                        self._bmftc[iB3D].x_m - progradation : self._bmftc[iB3D].x_m,
                    ] += new_marsh_height
                    # Store mass of overwash mineral sediment deposited across transect
                    self._bmftc[iB3D].mineral_dep[
                        self._bmftc[iB3D].startyear + time_step,
                        self._bmftc[iB3D].x_m - progradation : self._bmftc[iB3D].x_m,
                    ] += (
                        new_marsh_height * self._bmftc[iB3D].rhos * 1000
                    )  # [g] Mass of pure mineral sediment deposited by overwash

                # Spread 50% of overwash bay flux evenly across bay bottom
                bay_accrete = sum_bay_dep / (
                    self._bmftc[iB3D].bfo - progradation
                )  # [m] Vertical accretion of bay bottom from overwash deposition in bay
                self._bmftc[iB3D]._db = (
                    self._bmftc[iB3D].db + bay_accrete
                )  # Update bay depth

                # Remove barrier at front and set to msl to account for shoreline change
                self._bmftc[iB3D].elevation[
                    self._bmftc[iB3D].startyear + time_step, -off:
                ] = (
                    self._bmftc[iB3D].msl[self._bmftc[iB3D].startyear + time_step]
                    + self._bmftc[iB3D].amp
                )

            else:
                elevation_change_b3d = np.flip(elevation_change_b3d)
                marsh_barrier_width = self._bmftc[iB3D].B - self._bmftc[iB3D].x_m
                x_m_change = abs(
                    math.floor(len(elevation_change_b3d) - marsh_barrier_width)
                )  # Location of marsh edge within elevation_change_b3d

                # Add subaerial elevation change
                self._bmftc[iB3D].elevation[
                    self._bmftc[iB3D].startyear + time_step,
                    -len(elevation_change_b3d[x_m_change:]) :,
                ] += elevation_change_b3d[x_m_change:]
                # Store mass of overwash mineral sediment deposited across transect
                self._bmftc[iB3D].mineral_dep[
                    self._bmftc[iB3D].startyear + time_step,
                    -len(elevation_change_b3d[x_m_change:]) :,
                ] += (
                    elevation_change_b3d[x_m_change:] * self._bmftc[iB3D].rhos * 1000
                )  # [g] Mass of pure mineral sediment deposited by overwash

                # Determine volume of sed deposited past initial marsh edge and into bay
                sum_marsh_dep = np.sum(elevation_change_b3d[:x_m_change]) * (
                    1 - self._OWspread
                )  # [m^3] Volume of overwash deposition landward of marsh edge deposited as marsh
                sum_bay_dep = (
                    np.sum(elevation_change_b3d[:x_m_change]) * self._OWspread
                )  # [m^3] Volume of overwash deposition landward of marsh edge deposited across bay bottom
                self._bmftc[iB3D]._Fow_min = max(
                    0, sum_bay_dep * self._bmftc[iB3D].rhos
                )  # [kg/yr] Overwash deposition into bay, volume converted to mass

                # Add volume of carryover from last time step
                sum_marsh_dep += (
                    self._bay_overwash_carryover
                )  # [m^3] Bay deposition from previous time step that wasn't enough to fully fill bay cell up to sea level

                # Calculate height of deposition needed to bring bay bottom up to avg marsh elevation
                new_marsh_height = self._bmftc[iB3D].db

                # Determine distance of marsh progradation from overwash deposition
                progradation_actual = (
                    sum_marsh_dep / new_marsh_height
                )  # [m] Amount of marsh progradation, in which all overwash dep in bay fills first bay cell, then second, and so on until no more sediment. Assumes overwash is not spread out over bay.
                progradation = int(
                    max(math.floor(progradation_actual[iB3D]), 0)
                )  # Round to nearest FULL meter
                self._bay_overwash_carryover = (
                    progradation_actual - progradation
                ) * new_marsh_height  # Save leftover volume of sediment to be added to sum_bay_dep in following time step

                if progradation > 0:
                    # Add subaqueous elevation change
                    self._bmftc[iB3D].elevation[
                        self._bmftc[iB3D].startyear + time_step,
                        self._bmftc[iB3D].x_m - progradation : self._bmftc[iB3D].x_m,
                    ] += new_marsh_height
                    # Store mass of overwash mineral sediment deposited across transect
                    self._bmftc[iB3D].mineral_dep[
                        self._bmftc[iB3D].startyear + time_step,
                        self._bmftc[iB3D].x_m - progradation : self._bmftc[iB3D].x_m,
                    ] += (
                        new_marsh_height * self._bmftc[iB3D].rhos * 1000
                    )  # [g] Mass of pure mineral sediment deposited by overwash

                # Spread 50% of overwash bay flux evenly across bay bottom
                bay_accrete = sum_bay_dep / (
                    self._bmftc[iB3D].bfo - progradation
                )  # [m] Vertical accretion of bay bottom from overwash deposition in bay
                self._bmftc[iB3D]._db = (
                    self._bmftc[iB3D].db + bay_accrete
                )  # Update bay depth

            # Calculate new marsh and "forest" edge positions after overwash
            self._bmftc[iB3D]._x_m = self._bmftc[iB3D].x_m - progradation

            try:
                self._bmftc[iB3D]._x_f = max(
                    self._bmftc[iB3D].x_m + 1,
                    np.where(
                        self._bmftc[iB3D].elevation[
                            self._bmftc[iB3D].startyear + time_step, :
                        ]
                        > self._bmftc[iB3D].msl[self._bmftc[iB3D].startyear + time_step]
                        + self._bmftc[iB3D].amp
                        - self._bmftc[iB3D].Dmin
                        + 0.03
                    )[0][0],
                )
            except IndexError:
                self._bmftc[iB3D]._x_f = self._bmftc[iB3D].B
                # If x_f can't be found, barrier has drowned
                self._bmftc[iB3D]._dur = time_step
                self._bmftc[iB3D]._endyear = self._bmftc[iB3D].startyear + time_step
                self._BMFTC_Break = True
                print(
                    "PyBMFT-C Simulation Break: marsh has completely drowned or basin is completely full"
                )
                return  # End simulation

            # Store new positions
            self._bmftc[iB3D].Marsh_edge[
                self._bmftc[iB3D].startyear + time_step
            ] = self._bmftc[
                iB3D
            ].x_m  # Save to array
            self._bmftc[iB3D].Forest_edge[
                self._bmftc[iB3D].startyear + time_step
            ] = self._bmftc[
                iB3D
            ].x_f  # Save to array

            # Determine change in x_b location
            self._x_b_TS[iB3D][time_step] = self._bmftc[iB3D].x_b  # Save to array

            # Determine new fetch based on change in opposite marsh - both fetches should be exactly the same!
            self._bmftc[iB3D]._bfo = self._bmftc[iB3D].bfo - progradation
            self._bmftc[iB3D].fetch[
                self._bmftc[iB3D].startyear + time_step
            ] = self._bmftc[
                iB3D
            ].bfo  # Save to array

            # Update marsh scarp height parameter
            self._bmftc[iB3D]._dmo = (
                self._bmftc[iB3D].msl[self._bmftc[iB3D].startyear + time_step]
                + self._bmftc[iB3D].amp
                - self._bmftc[iB3D].elevation[
                    self._bmftc[iB3D].startyear + time_step, self._bmftc[iB3D].x_m
                ]
            )

            # Store landscape type widths for this time step
            if int(math.floor(self._x_s_offset[iB3D])) < 0:
                barrier_width = (
                    len(
                        self._bmftc[iB3D].elevation[
                            self._bmftc[iB3D].startyear + time_step,
                            self._bmftc[iB3D].x_f :,
                        ]
                    )
                    + off
                )
            elif int(math.floor(self._x_s_offset[iB3D])) > 0:
                barrier_width = (
                    len(
                        self._bmftc[iB3D].elevation[
                            self._bmftc[iB3D].startyear + time_step,
                            self._bmftc[iB3D].x_f :,
                        ]
                    )
                    - off
                )
            else:
                barrier_width = len(
                    self._bmftc[iB3D].elevation[
                        self._bmftc[iB3D].startyear + time_step, self._bmftc[iB3D].x_f :
                    ]
                )
            BB_marsh_width = (
                self._bmftc[iB3D].elevation[
                    self._bmftc[iB3D].startyear + time_step,
                    self._bmftc[iB3D].x_m : self._bmftc[iB3D].x_f,
                ]
                > self._bmftc[iB3D].msl[self._bmftc[iB3D].startyear + time_step]
                + self._bmftc[iB3D].amp
                - self._bmftc[iB3D].Dmax
            ).sum()
            BB_marsh_pond_width = (
                self._bmftc[iB3D].elevation[
                    self._bmftc[iB3D].startyear + time_step,
                    self._bmftc[iB3D].x_m : self._bmftc[iB3D].x_f,
                ]
                < self._bmftc[iB3D].msl[self._bmftc[iB3D].startyear + time_step]
                + self._bmftc[iB3D].amp
                - self._bmftc[iB3D].Dmax
            ).sum()
            self._LandscapeTypeWidth_TS[iB3D][time_step, :] = [
                barrier_width,
                BB_marsh_width,
                self._bmftc[iB3D].bfo,
                BB_marsh_pond_width,
            ]
