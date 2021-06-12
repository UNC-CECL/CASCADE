"""Chome Coupler

This module couples Barrier3d with CHOME to simulate community dynamics. Longer explanation with references go here.

References
----------

.. [1] Zachary Williams citation


Notes
---------

"""
import numpy as np
from chome import Chome


class ChomeCoupler:
    """Couple Barrier3d with CHOME to simulate community dynamics

    Examples
    --------
    >>> from chome_coupler import ChomeCoupler
    >>> chome_coupler = ChomeCoupler(barrier3d, nourishments)
    >>> chome_coupler.update()
    """

    def __init__(
            self,
            barrier3d,
            dune_design_elevation,
            dy=500,
            number_of_communities=1,
            name="default",
    ):
        """The ChomeBuyer module.

        Parameters
        ----------
        name: string, optional
            Name of simulation
        dy: int, optional
            Static length of island segment in Barrier3D (comprised of 10x10 cells)
        barrier3d: class
            Class of barrier3d models
        number_of_communities: int, optional
            Number of communities contained within the Barrier3D domain
        dune_design_elevation: list of floats, optional
            Elevation for rebuilding dunes [m NAVD88]; one value for each Barrier3D domain
        """
        self._number_of_communities = number_of_communities
        self._communities = []
        ny = len(barrier3d)

        # create groups of barrier3d indices for each community
        community_index = np.arange(ny)
        alongshore_community_length = int(np.ceil(ny / self._number_of_communities))
        self._community_index = [
            community_index[x: x + alongshore_community_length]
            for x in range(0, len(community_index), alongshore_community_length)
        ]

        if len(self._community_index) < self._number_of_communities:
            self._number_of_communities = len(self._community_index)
            print("Number of communities updated based on number of Barrier3D domains")

        for iCommunity in range(self._number_of_communities):
            average_barrier_height = []
            average_shoreface_depth = []
            average_dune_design_height = []

            for iB3D in self._community_index[iCommunity]:
                bh_array = np.array(barrier3d[iB3D].DomainTS[0]) * 10  # in meters

                average_barrier_height.append(bh_array[bh_array > 0].mean())
                average_shoreface_depth.append(
                    barrier3d[iB3D].DShoreface * 10
                )  # should all be the same, meters
                average_dune_design_height.append(
                    dune_design_elevation[iB3D] - (barrier3d[iB3D].BermEl * 10)
                )  # again, meters

            # initialize chome model instance
            self._communities.append(
                Chome(
                    name=name,
                    barrier_island_height=np.mean(average_barrier_height),
                    shoreface_depth=np.mean(average_shoreface_depth),
                    alongshore_domain_extent=dy,  # 500 m for Barrier3D
                    dune_height_build=np.mean(average_dune_design_height),
                )
            )

    def update(self, barrier3d, nourishments):

        time_index_b3d = barrier3d.time_index
        time_index_chome = self._communities[
            0
        ].time_index  # NOTE: time is updated at the end of the chome model loop
        nourish_now = []
        rebuild_dune_now = []
        nourishment_volume = []

        # calculate physical variables needed to update CHOME for each Barrier3D cell, then group by community
        for iCommunity in range(self._number_of_communities):
            # x_t_TS, x_s_TS, x_b_TS, h_b_TS = [[] for _ in range(4)]
            avg_barrier_height = (
                []
            )  # once averaged, these are saved as a time series in CHOME
            avg_change_shoreline_position = []
            avg_beach_width = []
            avg_interior_width = []
            avg_dune_height = []
            avg_dune_sand_volume_for_rebuild = []

            for iB3D in self._community_index[iCommunity]:
                bh_array = np.array(barrier3d[iB3D].DomainTS[-1]) * 10  # meters
                change_in_shoreline_position = (
                                                       barrier3d[iB3D].x_s_TS[-1] - barrier3d[iB3D].x_s_TS[-2]
                                               ) * 10  # meters
                beach_width = (
                        nourishments.beach_width[time_index_b3d - 2]
                        - change_in_shoreline_position
                )  # meters
                interior_width = barrier3d[iB3D].InteriorWidth_AvgTS * 10  # meters

                avg_barrier_height.append(bh_array[bh_array > 0].mean())
                avg_change_shoreline_position.append(change_in_shoreline_position)
                avg_beach_width.append(beach_width)
                avg_interior_width.append(interior_width)

            # update chome model variables before advancing one time step
            self._communities[iCommunity].barr_elev[time_index_chome] = np.mean(
                avg_barrier_height
            )
            self._communities[iCommunity].bw_erosion_rate = np.mean(
                avg_change_shoreline_position
            )
            self._communities[iCommunity].beach_width[time_index_chome] = np.mean(
                avg_beach_width
            )
            # self._communities[iCommunity].interior_width[time_index_chome] = np.mean(avg_interior_width)
            self._communities[iCommunity].dune_height[time_index_chome] = np.mean(
                avg_dune_height
            )
            # self._communities[iCommunity].dune_sand_volume[time_index_chome] = np.mean(avg_dune_sand_volume_for_rebuild)

            self._communities[iCommunity].update()

            # specify for each community, and iB3D cell, if it is a nourishment and/or dune rebuilding year
            nourish_now.append(
                [
                    self._communities[iCommunity].nourish_now[time_index_chome]
                    for _ in self._community_index[iCommunity]
                ]
            )
            rebuild_dune_now.append(
                [
                    self._communities[iCommunity].rebuild_dune_now[time_index_chome]
                    for _ in self._community_index[iCommunity]
                ]
            )
            nourishment_volume.append(
                [
                    self._communities[iCommunity].nourishment_volume[time_index_chome]
                    for _ in self._community_index[iCommunity]
                ]
            )

        return (
            np.hstack(nourish_now),
            np.hstack(rebuild_dune_now),
            np.hstack(nourishment_volume),
        )

    @property
    def communities(self):
        return self._communities
