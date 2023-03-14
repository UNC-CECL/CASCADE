"""Couple CHOM with Barrier3D (in development)

This module couples Barrier3D with the Coastal Home Ownership Model (CHOM), an agent-based model for coastal real estate
markets.

References
----------

.. [1]


Notes
---------
The user specifies the number of communities to be described by the Barrier3D models, and the initialization of
the coupler determines the (sequential) Barrier3d indices that are associated with each CHOM model instance. For
example, if 6 Barrier3D models are used to describe the alongshore domain (each 500 m long), and 2 communities are
specified, then indices [0,1,2] and [3,4,5] are associated with each community. Barrier3D variables for coupling with
CHOM are then averaged for each set of indices before passing to CHOM.


"""
import numpy as np
from chom import Chom

from .roadway_manager import rebuild_dunes

dm3_to_m3 = 1000  # convert from cubic decameters to cubic meters


def create_communities(ny, number_of_communities):
    community_index_initial = np.arange(ny)
    alongshore_community_length = int(np.ceil(ny / number_of_communities))
    community_index = [
        community_index_initial[x : x + alongshore_community_length]
        for x in range(0, len(community_index_initial), alongshore_community_length)
    ]

    if len(community_index) < number_of_communities:
        number_of_communities = len(community_index)
        print(
            "Number of communities simulated in CHOM updated based on number of Barrier3D domains"
        )

    return number_of_communities, community_index


def community_initial_statistics(
    community_indices, barrier3d, dune_design_elevation, alongshore_length_b3d
):
    avg_beach_width = []
    avg_barrier_height = []
    avg_shoreface_depth = []
    avg_dune_design_height = []
    avg_interior_width = []
    avg_dune_height = []
    avg_dune_width = []

    for iB3D in community_indices:
        # Barrier3D in decameters --> convert to meters for CHOM
        bh_array = np.array(barrier3d[iB3D].DomainTS[0]) * 10
        avg_barrier_height.append(bh_array[bh_array > 0].mean())

        interior_width = np.array(barrier3d[iB3D].InteriorWidth_AvgTS[0]) * 10
        avg_interior_width.append(interior_width)

        avg_shoreface_depth.append(barrier3d[iB3D].DShoreface * 10)
        avg_dune_design_height.append(
            dune_design_elevation[iB3D] - (barrier3d[iB3D].BermEl * 10)
        )
        initial_beach_width = (
            int(barrier3d[iB3D].BermEl / barrier3d[iB3D]._beta) * 10
        )  # m
        avg_beach_width.append(initial_beach_width)

        # maximum height of each (alongshore) row in DuneDomain (the largest foredune height), then average
        dune_height = barrier3d[iB3D].DuneDomain[0, :, :].max(axis=1)
        avg_dune_height.append(np.mean(dune_height) * 10)

        avg_dune_width.append(barrier3d[iB3D]._DuneWidth * 10)

        # community length is the Barrier3D discretization * # of Barrier3D domains describing community
        alongshore_community_length = len(community_indices) * alongshore_length_b3d

    avg_beach_width = np.mean(avg_beach_width)
    avg_barrier_height = np.mean(avg_barrier_height)
    avg_shoreface_depth = np.mean(avg_shoreface_depth)
    avg_dune_design_height = np.mean(avg_dune_design_height)
    avg_interior_width = np.mean(avg_interior_width)
    avg_dune_height = np.mean(avg_dune_height)
    avg_dune_width = np.mean(avg_dune_width)

    return (
        avg_beach_width,
        avg_barrier_height,
        avg_shoreface_depth,
        avg_dune_design_height,
        avg_interior_width,
        avg_dune_height,
        avg_dune_width,
        alongshore_community_length,
    )


def community_update_statistics(
    community_indices, barrier3d, time_index_b3d, nourishments, dune_design_elevation
):
    # once averaged, these are saved as a time series in CHOM
    (
        avg_barrier_height_msl,
        avg_change_shoreline_position,
        avg_beach_width,
        avg_dune_height,
        total_dune_sand_volume_rebuild,
    ) = ([] for _ in range(5))

    for iB3D in community_indices:
        # Barrier3D in dam MHW --> convert to m MSL for CHOM; b/c time_index in B3D is updated at the end
        # of the time loop, time_index-1 is the current time step for passing variable to CHOM
        bh_array = np.array(barrier3d[iB3D].DomainTS[time_index_b3d - 1]) * 10  # m MHW
        avg_barrier_height = bh_array[bh_array > 0].mean()
        avg_barrier_height_msl.append(
            avg_barrier_height + barrier3d[iB3D]._MHW
        )  # m NAVD88 (~MSL)

        change_in_shoreline_position = (
            barrier3d[iB3D].x_s_TS[-1] - barrier3d[iB3D].x_s_TS[-2]
        ) * 10  # meters
        avg_change_shoreline_position.append(change_in_shoreline_position)

        # nourishments.beach_width not updated until we call nourishments (after chom_coupler), so just calculate here
        beach_width = (
            nourishments[iB3D].beach_width[time_index_b3d - 2]
            - change_in_shoreline_position
        )
        avg_beach_width.append(beach_width)

        # maximum height of each row in DuneDomain, then average
        dune_height = barrier3d[iB3D].DuneDomain[time_index_b3d - 1, :, :].max(axis=1)
        avg_dune_height.append(np.mean(dune_height) * 10)

        # volume needed to rebuild dunes up to design height
        artificial_max_dune_height = dune_design_elevation[iB3D] - (
            barrier3d[iB3D].BermEl * 10
        )
        _, rebuild_dune_volume = rebuild_dunes(
            barrier3d[iB3D].DuneDomain[time_index_b3d - 1, :, :],  # dam
            max_dune_height=artificial_max_dune_height,  # in m
            min_dune_height=artificial_max_dune_height,  # in m
            dz=10,  # specifies dune domain is in dam
            rng=True,  # adds stochasticity to dune height (seeded)
        )
        total_dune_sand_volume_rebuild.append(rebuild_dune_volume * dm3_to_m3)

    avg_beach_width = np.mean(avg_beach_width)
    avg_barrier_height_msl = np.mean(avg_barrier_height_msl)
    avg_dune_height = np.mean(avg_dune_height)
    avg_change_shoreline_position = np.mean(avg_change_shoreline_position)
    total_dune_sand_volume_rebuild = np.sum(total_dune_sand_volume_rebuild)

    return (
        avg_beach_width,
        avg_barrier_height_msl,
        avg_dune_height,
        avg_change_shoreline_position,
        total_dune_sand_volume_rebuild,
    )


class ChomCoupler:
    """Couple Barrier3d with C-HOM to simulate community dynamics

    Examples
    --------
    # >>> from cascade.chom_coupler import ChomCoupler
    # >>> chom_coupler = ChomCoupler(barrier3d=, dune_design_elevation=)
    # >>> chom_coupler.update()
    """

    def __init__(
        self,
        barrier3d,
        dune_design_elevation,
        name="default",
        total_time=100,
        alongshore_length_b3d=500,
        number_of_communities=1,
        sand_cost=10,
        taxratio_oceanfront=1,
        external_housing_market_value_oceanfront=6e5,
        external_housing_market_value_nonoceanfront=4e5,
        fixed_cost_beach_nourishment=2e6,
        fixed_cost_dune_nourishment=2e5,
        nourishment_cost_subsidy=10e6,
        house_footprint_x=15,
        house_footprint_y=20,
        beach_full_cross_shore=70,
    ):
        """The ChomCoupler module.

        Parameters
        ----------
        barrier3d: class
            Class of barrier3d models
        dune_design_elevation: list of floats
            Elevation for rebuilding dunes [m NAVD88]; one value for each Barrier3D domain
        name: string, optional
            Name of simulation
        total_time: int, optional,
            Number of time steps [yr]
        alongshore_length_b3d: int, optional
            Static length of island segment in Barrier3D (comprised of 10x10 cells)
        number_of_communities: int, optional
            Number of communities (CHOM model instances) described by the alongshore section count (Barrier3D grids)
        sand_cost: int, optional
            Unit cost of sand $/m^3
        taxratio_oceanfront: float, optional
            The proportion of tax burden placed on oceanfront row for beach management
        external_housing_market_value_oceanfront:  , optional
            Value of comparable housing option outside the coastal system
        external_housing_market_value_nonoceanfront:  , optional
            Value of comparable housing option outside the coastal system
        fixed_cost_beach_nourishment: int, optional
            Fixed cost of 1 nourishment project
        fixed_cost_dune_nourishment: int, optional
            Fixed cost of building dunes once
        nourishment_cost_subsidy: int, optional
            Subsidy on cost of entire nourishment plan
        house_footprint_x: int, optional
            Length of house footprint in the cross-shore (meters)
        house_footprint_y: int, optional
            Length of house footprint in the alongshore (meters)
        beach_full_cross_shore: float, optional
            The cross-shore extent (meters) of fully nourished beach (i.e., the community desired beach width) [m]

        """
        self._dune_design_elevation = dune_design_elevation
        self._chom = []
        ny = len(barrier3d)

        # create groups of barrier3d indices for each community
        [self._number_of_communities, self._community_index] = create_communities(
            ny, number_of_communities
        )

        # find the average statistics for each community and then initialize CHOM models
        for iCommunity in range(self._number_of_communities):
            # calculate the average environmental statistics for each community
            [
                avg_beach_width,
                avg_barrier_height_msl,
                avg_shoreface_depth,
                avg_dune_design_height,
                avg_interior_width,
                avg_dune_height,
                avg_dune_width,
                self._alongshore_community_length,
            ] = community_initial_statistics(
                community_indices=self._community_index[iCommunity],
                barrier3d=barrier3d,
                dune_design_elevation=dune_design_elevation,
                alongshore_length_b3d=alongshore_length_b3d,
            )

            # initialize CHOM model instance -- ZACHK, are these still correct or are there more?
            self._chom.append(
                Chom(
                    name=name,
                    total_time=total_time,
                    average_interior_width=avg_interior_width,
                    barrier_island_height=avg_barrier_height_msl,
                    beach_width=avg_beach_width,
                    dune_height=avg_dune_height,
                    shoreface_depth=avg_shoreface_depth,
                    dune_width=avg_dune_width,
                    dune_height_build=avg_dune_design_height,
                    alongshore_domain_extent=self._alongshore_community_length,
                    shoreline_retreat_rate=0,  # no knowledge of shoreline retreat rate at t=0
                    sand_cost=sand_cost,
                    taxratio_oceanfront=taxratio_oceanfront,
                    external_housing_market_value_oceanfront=external_housing_market_value_oceanfront,
                    external_housing_market_value_nonoceanfront=external_housing_market_value_nonoceanfront,
                    fixed_cost_beach_nourishment=fixed_cost_beach_nourishment,
                    fixed_cost_dune_nourishment=fixed_cost_dune_nourishment,
                    nourishment_cost_subsidy=nourishment_cost_subsidy,
                    house_footprint_x=house_footprint_x,
                    house_footprint_y=house_footprint_y,
                    beach_full_cross_shore=beach_full_cross_shore,
                )
            )

    def update(self, barrier3d, nourishments, community_break):
        # time is updated at the end of the CHOM and Barrier3d model loop
        time_index_b3d = barrier3d[0].time_index
        time_index_chom = self._chom[0].time_index
        nourish_now, rebuild_dune_now, nourishment_volume = (
            np.zeros(len(barrier3d)) for _ in range(3)
        )

        # calculate physical variables needed to update CHOM for each Barrier3D model, then group by community
        for iCommunity in range(self._number_of_communities):
            community_indices = self._community_index[iCommunity]

            [
                avg_beach_width,
                avg_barrier_height_msl,
                avg_dune_height,
                avg_change_shoreline_position,
                total_dune_sand_volume_rebuild,
            ] = community_update_statistics(
                community_indices=community_indices,
                barrier3d=barrier3d,
                time_index_b3d=time_index_b3d,
                nourishments=nourishments,
                dune_design_elevation=self.dune_design_elevation,
            )

            # update CHOM model variables before advancing one time step (but only if the community hasn't
            # been abandoned)
            if any(community_break[community_indices[0] : community_indices[-1] + 1]):
                pass
            else:
                self._chom[
                    iCommunity
                ].height_above_msl = avg_barrier_height_msl  # m MSL
                self._chom[iCommunity].bw_erosion_rate[
                    time_index_chom
                ] = avg_change_shoreline_position
                self._chom[iCommunity].beach_width[time_index_chom] = avg_beach_width
                self._chom[iCommunity].dune_height[time_index_chom] = avg_dune_height
                self._chom[iCommunity].dune_sand_volume[
                    time_index_chom
                ] = total_dune_sand_volume_rebuild
                # NOTE: dune design ele cannot yet be updated at each time step in CHOM -- potential future update

                self._chom[iCommunity].update()

                # specify for each community, and iB3D cell, if it is a nourishment and/or dune rebuilding year
                nourish_now[community_indices] = self._chom[iCommunity].nourish_now[
                    time_index_chom
                ]
                rebuild_dune_now[community_indices] = self._chom[
                    iCommunity
                ].rebuild_dune_now[time_index_chom]
                nourishment_volume[community_indices] = (
                    self._chom[iCommunity].nourishment_volume[time_index_chom]
                    / self._chom[iCommunity].lLength
                )  # this is in m^3 in CHOM, divide by the length of the community to get m^3/m

        return (nourish_now, rebuild_dune_now, nourishment_volume)

    @property
    def chom(self):
        return self._chom

    @property
    def dune_design_elevation(self):
        return self._dune_design_elevation

    @dune_design_elevation.setter
    def dune_design_elevation(self, value):
        self._dune_design_elevation = value
