"""Bulldozer

This module provides functions modifying a 3D grid (X,Y,Z) for human coastal management decision,
including 1) bulldozing sand from roadways and adding to the dune line ....

References
----------

.. [1] NCDOT (North Carolina Department of Transportation), 2008. NC 12 Transportation Improvements. http://www.ncdot.org/projects/
NC12/ (accessed October 24, 2007).
.. [2]
.. [3] P.D. Komar, 1998, Beach processes and sedimentation: Upper Saddle River, New Jersey, Prentice Hall , 544 p.
.. [4] Jaap H. Nienhuis, Jorge Lorenzo Trueba; Simulating barrier island response to sea level rise with the barrier
    island and inlet environment (BRIE) model v1.0 ; Geosci. Model Dev., 12, 4013–4030, 2019;
    https://doi.org/10.5194/gmd-12-4013-2019


Notes
---------

"""
import numpy as np
from scipy.interpolate import interp2d

dm3_to_m3 = 1000  # convert from cubic decameters to cubic meters


def bulldoze(
    xyz_interior_grid,
    yxz_dune_grid,
    road_ele=2.0,
    road_width=20,
    road_setback=30,
    dx=10,
    dy=10,
    dz=10,
    drown_threshold=0,
):
    r"""Remove overwash from roadway and put it back on the dune. Spreads sand evenly across dune cells.

    Parameters
    ----------
    xyz_interior_grid: array
        Interior barrier island topography [z units specified by dz; if dz=10, decameters NAVD88]
    yxz_dune_grid: array
        Dune topography [z units specified by dz; if dz=10, decameters]
    road_ele: float
        Road elevation [m; needs to be in same reference frame as xyz and yxz, typically NAVD88]
    road_width: int
        Width of roadway [m]
    road_setback: int
        Setback distance of roadway from edge of interior domain [m]
    dx: int
        Cross-shore discretization of x [m]
    dy: int
        Alongshore discretization of y [m]
    dz: int
        Vertical discretization of z [m]
    drown_threshold: float
        threshold for roadway drowning [m; needs to be in same reference frame as xyz and yxz, typically NAVD88]

    Returns
    -------
    array of float
        new_road_domain: in units of dx, dy, dz
        new_dune_domain: in units of dx, dy, dz
        overwash removal: in units of dz^3
    int
        roadway_drown: flag for if roadway borders water

    """

    # road parameters
    road_start = int(
        road_setback / dy
    )  # grid index for start of roadway in interior domain (for B3D, convert to dm)
    road_width = int(road_width / dx)
    road_end = road_start + road_width
    road_ele = (
        road_ele / dz
    )  # convert to units of grid (NOTE: in B3D default simulation, berm is 1.44 m)

    # remove sand from roadway (only account for positive values)
    old_road_domain = xyz_interior_grid[road_start:road_end, :]
    new_road_domain = np.zeros((road_width, np.size(old_road_domain, 1))) + road_ele
    road_overwash_removal = sum(old_road_domain - new_road_domain)
    road_overwash_removal[road_overwash_removal < 0] = 0

    # spread overwash removed from roadway equally over all dune cells
    number_dune_cells = np.size(yxz_dune_grid, 1)
    overwash_to_dune = np.transpose(
        [road_overwash_removal / number_dune_cells] * number_dune_cells
    )
    new_dune_domain = yxz_dune_grid + overwash_to_dune

    xyz_interior_grid[
        road_start:road_end, :
    ] = new_road_domain  # update interior domain

    # check if any water cells border the road on either side
    if (
        road_start > 0
        and (np.min(xyz_interior_grid[road_start - 1, :] * dz) <= drown_threshold)
    ) or (np.min(xyz_interior_grid[road_end + 1, :] * dz) <= drown_threshold):
        roadway_drown = True
    else:
        roadway_drown = False

    return (
        new_dune_domain,
        xyz_interior_grid,
        np.sum(road_overwash_removal),
        roadway_drown,
    )


def rebuild_dunes(
    yxz_dune_grid, max_dune_height=3.0, min_dune_height=2.4, dz=10, rng=True
):
    r"""Construct artificial dunes after each overwash or inundation event. Along NC-12, artificial dune heights range
    from 2.4 to 4.6 m (NCDOT, 2008; Overton, Judge, and Fisher, 2000; Magliocca et al., 2011). From a more recent pub
    (Valasquez et al., 2020), the average elevation of the road along NC-12 is 1.3 m (NAVD88); they find that in order
    for the road to not be vulnerable to overwash, the dune crest must be higher than 4.3 m (NAVD88), so here the
    default max_dune_height is 3 m (here, dune height is measured as the height above the dune toe).

    Note that while artificial dune geometry of a given height is constrained by the angle of repose of beach, here we
    just assume the dunes are built to a width capable of surviving impacts from a "moderate" storm (25-m wide along
    NC-12, corresponding to a 50-year event; NCDOT, 2008, Overton, Judge, and Fisher, 2000).

    For the given dune grid, we apply a linear gradient from the first to last dune row given the max and min dune
    height, with small random perturbations.

    Parameters
    ----------
    yxz_dune_grid: float or array of float
        Dune topography [z units specified by dz; if dz=10, decameters]
    max_dune_height: float
        Maximum dune height [m]
    min_dune_height: float
        Minimum dune height [m]
    dz: int
        Vertical discretization of z [m]
    rng: boolean
        If True, add random perturbations alongshore to dune height

    Returns
    -------
    new_dune_domain: float or array of float
        New yxz dune domain in z units specified by dz

    """

    # convert from m to grid z discretization
    max_height = max_dune_height / dz
    min_height = min_dune_height / dz
    ny = np.size(yxz_dune_grid, 0)
    nx = np.size(yxz_dune_grid, 1)

    if rng:
        # add some random perturbations to dune heights
        RNG = np.random.default_rng(seed=1973)

        dune_start_max = np.ones([ny]) * (
            # max_height + (-0.01 + (0.01 - (-0.01)) * np.random.rand(ny))
            max_height
            + (-0.01 + (0.01 - (-0.01)) * RNG.random(ny))
        )
        dune_start_min = np.ones([ny]) * (
            # min_height + (-0.01 + (0.01 - (-0.01)) * np.random.rand(ny))
            min_height
            + (-0.01 + (0.01 - (-0.01)) * RNG.random(ny))
        )
    else:
        dune_start_max = np.ones([ny]) * max_height
        dune_start_min = np.ones([ny]) * min_height

    # linearly interpolate (2D) from front row (max height) to back row (min height)
    x = [0, nx - 1]
    y = [np.arange(0, ny, 1)]
    z = np.transpose([dune_start_max, dune_start_min])
    f = interp2d(x, y, z)
    new_dune_domain = f(np.arange(0, nx, 1), np.arange(0, ny, 1))

    return new_dune_domain


def set_growth_parameters(
    yxz_dune_grid,
    Dmax,
    growthparam,
    original_growth_param=None,
    rmin=0.35,
    rmax=0.85,
):
    r"""Set dune growth rate to zero for next time step if the dune elevation (front row) is larger than the natural eq.
    dune height (DMax) -- we know from modeling work (Duran and Moore 2013) that the dune will not grow above
    the natural equilibrium height because of interactions between the wind field and the dune: too high and no sand
    flux

    Parameters
    ----------
    yxz_dune_grid: float or array of float
        Dune topography [decameters]
    Dmax: float
        Maximum natural equilibrium dune height [decameters]
    growthparam: float
        growth parameters from last time step, used in Houser formulation [unitless]
    original_growth_param: float, optional
        growth parameters from first time step, before humans interfered [unitless]
    rmin: float, optional
        Minimum growth rate - used if now original_growth_parm provided [unitless]
    rmax: float, optional
        Maximum growth rate - used if now original_growth_parm provided [unitless]

    Returns
    -------
    new_growth_param: array of float
        New growth parameter array that accounts for human modifications to dune height above/below the equilibrium

    """

    ny = np.size(yxz_dune_grid, 0)
    new_growth_param = np.copy(growthparam)
    rng = np.random.default_rng(seed=1973)

    for idune in range(ny):

        if yxz_dune_grid[idune, 0] > Dmax:  # if dune height above dmax, don't grow
            new_growth_param[0, idune] = 0

        else:
            # if dune is now below the Dmax (was formerly above), make sure it has a growth rate either the same as
            # before (if original growth rate provided) or a random number between rmin and rmax
            if growthparam[0, idune] == 0:
                if original_growth_param is not None:
                    new_growth_param[0, idune] = original_growth_param[0, idune]
                else:
                    new_growth_param[0, idune] = rmin + (rmax - rmin) * rng.random()

    return new_growth_param


class RoadwayManager:
    #
    # """Manage the road.
    #
    # Examples
    # --------
    # >>> from cascade.bulldozer import RoadwayManager
    # >>> roadways = RoadwayManager()
    # >>> roadways.update()
    # """

    def __init__(
        self,
        road_elevation,
        road_width,
        road_setback,
        road_relocation_setback,
        artificial_max_dune_ele,
        artificial_min_dune_ele,
        original_growth_param,
    ):
        """The RoadwayManager module.

        Parameters
        ----------

        """

        self._dunes_rebuilt = None
        self._road_overwash_removal = None
        self._road_width = road_width
        self._road_setback = road_setback
        self._road_ele = road_elevation
        self._road_relocation_setback = road_relocation_setback
        self._artificial_max_dune_ele = artificial_max_dune_ele
        self._artificial_min_dune_ele = artificial_min_dune_ele
        self._original_growth_param = original_growth_param
        self._drown_break = 0
        self._narrow_break = 0

    def update(self, barrier3d):

        # if dune line moved by one grid cell, substract that amount from the setback to keep road in the same place
        if (
            barrier3d.dune_migration
        ):  # (NOTE: I may have been able to just use the B3D variable shoreline_change here!)
            self._road_setback = self._road_setback - 10

            # with this shoreline change and dune migration, check if the roadway needs to be relocated
            if self._road_setback < 0:
                # relocate the road only if the width of the island allows it
                if (
                    self._road_relocation_setback
                    + (2 * self._road_width)  # bay shoreline buffer of 2x the roadway
                    > barrier3d._InteriorWidth_AvgTS[-1] * 10
                ):
                    self._narrow_break = 1
                    print(
                        "Island is to narrow for roadway to be relocated. Roadway drowned at {time} years".format(
                            time=barrier3d.time_index
                        )
                    )
                    return barrier3d  # exit program
                else:
                    self._road_setback = self._road_relocation_setback

        # now check the other shoreline: if the roadway gets eaten up by the back-barrier, then it is lost
        if (barrier3d._InteriorWidth_AvgTS[-1] * 10) <= self._road_setback:
            self._drown_break = 1
            print(
                "Roadway drowned at {time} years from the back-bay".format(
                    time=barrier3d.time_index
                )
            )
            return barrier3d  # exit program

        # manage that roadway!
        else:
            # bulldoze that road and put bulldozed sand back on the dunes
            (
                new_dune_domain,
                new_xyz_interior_domain,
                self._road_overwash_removal,
                roadway_drown,
            ) = bulldoze(
                road_ele=self._road_ele,  # m NAVD88
                road_width=self._road_width,
                road_setback=self._road_setback,
                xyz_interior_grid=barrier3d.InteriorDomain,  # interior domain from this last time step
                yxz_dune_grid=barrier3d.DuneDomain[
                    barrier3d.time_index - 1, :, :
                ],  # dune domain from this last time step
                dx=10,
                dy=10,
                dz=10,
            )

            # another check on the location of the back-barrier shoreline: if the roadway touches a water cell, drown!
            if roadway_drown:
                self._drown_break = 1
                print(
                    "Roadway drowned at {time} years from the back-bay".format(
                        time=barrier3d.time_index
                    )
                )
                return barrier3d

            self._road_overwash_removal = (
                self._road_overwash_removal * dm3_to_m3
            )  # convert from dm^3 to m^3

            # dune management: rebuild artifical dunes!
            if (
                self._artificial_max_dune_ele is None
                or self._artificial_min_dune_ele is None
            ):
                pass
            else:
                # in B3D, dune height is the height above the berm crest (keep this in m)
                self._artificial_max_dune_height = self._artificial_max_dune_ele - (
                    barrier3d.BermEl * 10
                )
                self._artificial_min_dune_height = self._artificial_min_dune_ele - (
                    barrier3d.BermEl * 10
                )

                # if any dune cell in the front row of dunes is less than a minimum threshold -- as measured above the berm
                # crest -- then rebuild artificial dunes (first row up to the max dune height, second row to the minimum?
                # -- nah, lets make it the max as well)
                if np.min(new_dune_domain[:, 0]) < (
                    self._artificial_min_dune_height / 10
                ):  # dune domain is in decameters; previously had np.max() and got a lot of small dune breaches that didn't get
                    # repaired -- not realistic
                    new_dune_domain = rebuild_dunes(  # decameters
                        new_dune_domain,
                        max_dune_height=self._artificial_max_dune_height,  # in m
                        min_dune_height=self._artificial_max_dune_height,  # in m
                        dz=10,  # specifies dune domain is in decameters
                        rng=True,  # adds stochasticity to dune height (seeded)
                    )
                    self._dunes_rebuilt = 1  # track when dunes are rebuilt
                else:
                    self._dunes_rebuilt = 0

            # set dune growth rate to zero for next time step if the dune elevation (front row) is larger than the natural eq.
            # dune height (Dmax)
            new_growth_parameters = set_growth_parameters(
                new_dune_domain,
                barrier3d.Dmax,
                barrier3d.growthparam,
                original_growth_param=self._original_growth_param,  # use original growth rates for resetting values
            )

            # update class variables
            barrier3d.DuneDomain[barrier3d.time_index - 1, :, :] = new_dune_domain
            barrier3d.InteriorDomain = new_xyz_interior_domain  # I think these class variables are altered above, but again here for safety
            barrier3d.DomainTS[barrier3d.time_index - 1] = new_xyz_interior_domain
            barrier3d.growthparam = new_growth_parameters

        return barrier3d
