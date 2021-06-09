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
            dy=500,
            dune_design_elevation=3.7,
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
        dune_design_elevation: float, optional
            Elevation for rebuilding dunes [m NAVD88]
        """
        self._number_of_communities = number_of_communities
        self._communities = []
        ny = len(barrier3d)

        # create groups of barrier3d indices for each community
        iB3D_communities = np.arange(ny)
        alongshore_community_length = int(np.ceil(ny / self._number_of_communities))
        self._iB3D_communities = [iB3D_communities[x:x + alongshore_community_length]
                                  for x in range(0, len(iB3D_communities), alongshore_community_length)]

        if len(self._iB3D_communities) < self._number_of_communities:
            self._number_of_communities = len(self._iB3D_communities)
            print("Number of communities updated based on number of Barrier3D domains")

        for iCommunity in range(self._number_of_communities):
            average_barrier_height = []
            average_shoreface_depth = []
            average_dune_design_height = []

            for iB3D in self._iB3D_communities[iCommunity]:
                bh_array = np.array(barrier3d[iB3D].DomainTS[0]) * 10  # in meters

                average_barrier_height.append(bh_array[bh_array > 0].mean())
                average_shoreface_depth.append(barrier3d[iB3D].DShoreface * 10)  # should all be the same, meters
                average_dune_design_height.append(
                    dune_design_elevation - (barrier3d[iB3D].BermEl * 10))  # again, meters

            # initialize chome model instance
            self._communities.append(
                Chome(
                    name=name,
                    barrier_island_height=np.mean(average_barrier_height),
                    # shoreface_depth=np.mean(average_shoreface_depth),
                    beach_nourishment_fill_width=dy,  # 500 m for Barrier3D
                    dune_height_build=np.mean(average_dune_design_height),
                )
            )

    def update(self, barrier3d, nourishments):

        time_index = barrier3d.time_index

        # calculate physical variables needed to update CHOME for each Barrier3D cell, then group by community
        for iCommunity in range(self._number_of_communities):
            average_barrier_height = []  # I NEED TO SAVE THESE AS TIME SERIES
            average_change_shoreline_position = []
            average_beach_width = []
            average_interior_width = []
            average_dune_height = []
            average_dune_sand_volume_for_rebuild = []

            for iB3D in self._iB3D_communities[iCommunity]:
                bh_array = np.array(barrier3d[iB3D].DomainTS[-1]) * 10  # in meters
                change_in_shoreline_position = (barrier3d[iB3D].x_s_TS[-1] - barrier3d[iB3D].x_s_TS[-2]) * 10  # meters
                beach_width = nourishments.beach_width[time_index - 2] - change_in_shoreline_position  # meters

                average_barrier_height.append(bh_array[bh_array > 0].mean())
                average_change_shoreline_position.append(change_in_shoreline_position)  # meters
                average_beach_width.append(beach_width)

            # initialize chome model instance
            self._communities[iCommunity].barr_elev = np.mean(average_barrier_height)

            self._communities[iCommunity].update()
