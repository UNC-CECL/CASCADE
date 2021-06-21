"""Chome Coupler

This module couples Barrier3d with CHOME to simulate community dynamics. Longer explanation with references go here.
User specifies the number of communities to be described by the Barrier3D model instances, and the initialization of the
coupler determines the (sequential) Barrier3d indices that are associated with each CHOME model instance. For example,
if 6 Barrier3D models are used to describe the alongshore domain (each 500 m long), and 2 communities are specified,
then indices [0,1,2] and [3,4,5] are associated with each community. Barrier3D variables for coupling with CHOME are
then averaged for each set of indices before passing to CHOME.

References
----------

.. [1] Zachary Williams citation


Notes
---------

"""
import numpy as np

from chome import Chome
from .roadway_manager import rebuild_dunes

dm3_to_m3 = 1000  # convert from cubic decameters to cubic meters


class ChomeCoupler:
    """Couple Barrier3d with CHOME to simulate community dynamics

    Examples
    --------
    >>> from cascade.chome_coupler import ChomeCoupler
    >>> chome_coupler = ChomeCoupler(barrier3d=, dune_design_elevation=)
    >>> chome_coupler.update()
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
        house_footprint=15,
    ):
        """The ChomeBuyer module.

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
            Number of communities (CHOME model instances) described by the alongshore section count (Barrier3D grids)
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

        """
        self._number_of_communities = number_of_communities
        self._dune_design_elevation = dune_design_elevation
        self._chome = []
        ny = len(barrier3d)

        # create groups of barrier3d indices for each community
        self._community_index = np.arange(ny)
        alongshore_community_length = int(np.ceil(ny / self._number_of_communities))
        self._community_index = [
            self._community_index[x: x + alongshore_community_length]
            for x in range(0, len(self._community_index), alongshore_community_length)
        ]

        if len(self._community_index) < self._number_of_communities:
            self._number_of_communities = len(self._community_index)
            print(
                "Number of communities simulated in CHOME updated based on number of Barrier3D domains"
            )

        for iCommunity in range(self._number_of_communities):
            # community length is the Barrier3D discretization * # of Barrier3D domains describing community
            self._alongshore_community_length = (
                len(self._community_index[iCommunity]) * alongshore_length_b3d
            )
            avg_beach_width = []
            avg_barrier_height = []
            avg_shoreface_depth = []
            avg_dune_design_height = []
            avg_interior_width = []
            avg_dune_height = []
            avg_dune_width = []

            for iB3D in self._community_index[iCommunity]:
                # Barrier3D in decameters --> convert to meters for CHOME
                bh_array = np.array(barrier3d[iB3D].DomainTS[0]) * 10
                avg_barrier_height.append(bh_array[bh_array > 0].mean())

                interior_width = np.array(barrier3d[iB3D].InteriorWidth_AvgTS[0]) * 10
                avg_interior_width.append(interior_width)

                avg_shoreface_depth.append(barrier3d[iB3D].DShoreface * 10)
                avg_dune_design_height.append(
                    self._dune_design_elevation[iB3D] - (barrier3d[iB3D].BermEl * 10)
                )
                initial_beach_width = (
                    int(barrier3d[iB3D].BermEl / barrier3d[iB3D]._beta) * 10
                )  # m
                avg_beach_width.append(initial_beach_width)

                # maximum height of each row in DuneDomain, then average
                dune_height = barrier3d[iB3D].DuneDomain[0, :, :].max(axis=1)
                avg_dune_height.append(np.mean(dune_height) * 10)

                avg_dune_width.append(barrier3d[iB3D]._DuneWidth * 10)

            # initialize chome model instance
            self._chome.append(
                Chome(
                    name=name,
                    total_time=total_time,
                    average_interior_width=np.mean(avg_interior_width),
                    barrier_island_height=np.mean(avg_barrier_height),
                    beach_width=np.mean(avg_beach_width),
                    dune_height=np.mean(avg_dune_height),
                    shoreface_depth=np.mean(avg_shoreface_depth),
                    dune_width=np.mean(avg_dune_width),
                    dune_height_build=np.mean(avg_dune_design_height),
                    alongshore_domain_extent=self._alongshore_community_length,
                    shoreline_retreat_rate=0,  # no knowledge of shoreline retreat rate at t=0
                    sand_cost=sand_cost,
                    taxratio_oceanfront=taxratio_oceanfront,
                    external_housing_market_value_oceanfront=external_housing_market_value_oceanfront,
                    external_housing_market_value_nonoceanfront=external_housing_market_value_nonoceanfront,
                    fixed_cost_beach_nourishment=fixed_cost_beach_nourishment,
                    fixed_cost_dune_nourishment=fixed_cost_dune_nourishment,
                    nourishment_cost_subsidy=nourishment_cost_subsidy,
                    house_footprint=house_footprint,
                )
            )

    def update(self, barrier3d, nourishments):
        # NOTE: dune design ele can not yet be updated at each time step in CHOME -- potential future functionality

        # time is updated at the end of the chome and barrier3d model loop
        time_index_b3d = barrier3d[0].time_index
        time_index_chome = self._chome[0].time_index
        nourish_now, rebuild_dune_now, nourishment_volume = [[] for _ in range(3)]

        # calculate physical variables needed to update CHOME for each Barrier3D cell, then group by community
        for iCommunity in range(self._number_of_communities):
            # once averaged, these are saved as a time series in CHOME
            (
                avg_barrier_height,
                avg_change_shoreline_position,
                avg_beach_width,
                avg_dune_height,
                total_dune_sand_volume_rebuild,
            ) = [[] for _ in range(5)]

            for iB3D in self._community_index[iCommunity]:
                # Barrier3D in decameters --> convert to meters for CHOME; b/c time_index in B3D is updated at the end
                # of the time loop, time_index-1 is the current time step for passing variable to CHOME
                bh_array = np.array(barrier3d[iB3D].DomainTS[time_index_b3d - 1]) * 10  # meters
                avg_barrier_height.append(bh_array[bh_array > 0].mean())

                change_in_shoreline_position = (
                    barrier3d[iB3D].x_s_TS[-1] - barrier3d[iB3D].x_s_TS[-2]
                ) * 10
                avg_change_shoreline_position.append(change_in_shoreline_position)

                # nourishments.beach_width not updated until we call nourishments, so just calculate here
                beach_width = (
                    nourishments[iB3D].beach_width[time_index_b3d - 2]
                    - change_in_shoreline_position
                )
                avg_beach_width.append(beach_width)

                # maximum height of each row in DuneDomain, then average
                dune_height = (
                    barrier3d[iB3D].DuneDomain[time_index_b3d - 1, :, :].max(axis=1)
                )
                avg_dune_height.append(np.mean(dune_height) * 10)

                # volume needed to rebuild dunes up to design height
                artificial_max_dune_height = self._dune_design_elevation[iB3D] - (
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

            # update chome model variables before advancing one time step
            self._chome[iCommunity].barr_elev[time_index_chome] = np.mean(
                avg_barrier_height
            )
            self._chome[iCommunity].bw_erosion_rate[time_index_chome] = np.mean(
                avg_change_shoreline_position
            )
            self._chome[iCommunity].beach_width[time_index_chome] = np.mean(
                avg_beach_width
            )
            self._chome[iCommunity].dune_height[time_index_chome] = np.mean(
                avg_dune_height
            )
            self._chome[iCommunity].dune_sand_volume[time_index_chome] = np.sum(
                total_dune_sand_volume_rebuild
            )

            self._chome[iCommunity].update()

            # specify for each community, and iB3D cell, if it is a nourishment and/or dune rebuilding year
            nourish_now.append(
                [
                    self._chome[iCommunity].nourish_now[time_index_chome]
                    for _ in self._community_index[iCommunity]
                ]
            )
            rebuild_dune_now.append(
                [
                    self._chome[iCommunity].rebuild_dune_now[time_index_chome]
                    for _ in self._community_index[iCommunity]
                ]
            )
            nourishment_volume.append(
                [
                    # this is in m^3 in chome, divide by the length of the community to get m^3/m
                    (
                        self._chome[iCommunity].nourishment_volume[
                            time_index_chome
                        ]
                        / self._chome[iCommunity].lLength
                    )
                    for _ in self._community_index[iCommunity]
                ]
            )

        return (
            np.hstack(nourish_now),
            np.hstack(rebuild_dune_now),
            np.hstack(nourishment_volume),
        )

    @property
    def chome(self):
        return self._chome

    @property
    def dune_design_elevation(self):
        return self._dune_design_elevation

    @dune_design_elevation.setter
    def dune_design_elevation(self, value):
        self._dune_design_elevation = value