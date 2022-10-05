"""BRIE Coupler

This module couples barrier3d with brie via alongshore sediment transport

References
----------

.. [1] Jaap H. Nienhuis, Jorge Lorenzo Trueba; Simulating barrier island response to sea level rise with the barrier
    island and inlet environment (BRIE) model v1.0 ; Geosci. Model Dev., 12, 4013–4030, 2019;
    https://doi.org/10.5194/gmd-12-4013-2019


Notes
---------

"""
import numpy as np
import math

from yaml import full_load, dump
from brie import Brie
from barrier3d import Barrier3d


def set_yaml(var_name, new_vals, file_name):
    with open(file_name) as f:
        doc = full_load(f)
    doc[var_name] = new_vals
    with open(file_name, "w") as f:
        dump(doc, f)


def batchB3D(subB3D):

    """parallelize update function for each B3D sub-grid (spread overwash routing algorithms to different cores)"""

    subB3D.update()

    # calculate the diff in shoreface toe, shoreline, and height of barrier (dam)
    sub_x_t_dt = (subB3D.x_t_TS[-1] - subB3D.x_t_TS[-2]) * 10
    sub_x_s_dt = (subB3D.x_s_TS[-1] - subB3D.x_s_TS[-2]) * 10
    sub_h_b_dt = (subB3D.h_b_TS[-1] - subB3D.h_b_TS[-2]) * 10

    return sub_x_t_dt, sub_x_s_dt, sub_h_b_dt, subB3D


def initialize_equal(
    datadir,
    brie,
    slr_constant,
    rmin,
    rmax,
    background_erosion,
    parameter_file,
    storm_file,
    dune_file,
    elevation_file,
    MHW=0.46,
    beta=0.04
):
    """
    for each B3D subgrid, set the initial shoreface geometry equal to what is set in brie (some random
    # perturbations); all other B3D variables are set equal

    :param datadir: the directory containing the Barrier3D parameter file
    :param brie: the brie class
    :param slr_constant: specifies whether SLR will be constant
    :param rmin: the minimum dune growth rate
    :param rmax: the maximum dune growth rate
    :param background_erosion: Rate of shoreline retreat attributed to gradients in alongshore transport
    :param parameter_file: name of the Barrier3D parameter file
    :param storm_file: name of the Barrier3D storms file
    :param dune_file: name of the Barrier3D dunes file
    :param elevation_file: name of the Barrier3D elevation file
    :param MHW: elevation of mean high water [m NAVD88]
    :param beta: beach slope for runup calculations

    :return: barrier3d

    """
    barrier3d = []

    for iB3D in range(brie.ny):

        fid = datadir + parameter_file
        parameter_file_prefix = parameter_file.replace("-parameters.yaml", "")

        # update yaml file (these are the only variables that I'm like to change from default)
        set_yaml("Shrub_ON", 0, fid)  # make sure that shrubs are turned off
        set_yaml(
            "TMAX", brie.nt, fid
        )  # [yrs] Duration of simulation (if brie._dt = 1 yr, set to ._nt)
        set_yaml(
            "BarrierLength", brie._dy, fid
        )  # [m] Static length of island segment (comprised of 10x10 cells)
        set_yaml(
            "DShoreface", brie.d_sf, fid
        )  # [m] Depth of shoreface (set to brie depth, function of wave height)
        set_yaml(
            "LShoreface", float(brie.x_s[iB3D] - brie.x_t[iB3D]), fid
        )  # [m] Length of shoreface (calculate from brie variables, shoreline - shoreface toe)
        set_yaml(
            "ShorefaceToe", float(brie.x_t[iB3D]), fid
        )  # [m] Start location of shoreface toe
        set_yaml(
            "BayDepth", brie._bb_depth, fid
        )  # [m] Depth of bay behind island segment (set to brie bay depth)
        set_yaml(
            "DuneParamStart", True, fid
        )  # Dune height will come from external file
        # Rate of shoreline retreat attributed to gradients in alongshore transport; (-) = erosion, (+) = acc [m / y]
        if np.size(background_erosion) > 1:
            set_yaml("Rat", background_erosion[iB3D], fid)
        else:
            set_yaml("Rat", background_erosion, fid)
        set_yaml(
            "RSLR_Constant", slr_constant, fid
        )  # Relative sea-level rise rate will be constant, otherwise logistic growth function used
        set_yaml(
            "RSLR_const", brie._slr[0], fid
        )  # [m / y] Relative sea-level rise rate; initialized in brie, but saved as time series, so we use index 0
        set_yaml(
            "k_sf", float(brie.k_sf), fid
        )  # [m^3 / m / y] Shoreface flux rate constant (function of wave parameters from brie)
        set_yaml(
            "s_sf_eq", float(brie.s_sf_eq), fid
        )  # Equilibrium shoreface slope (function of wave and sediment parameters from brie)
        set_yaml(
            "GrowthParamStart", False, fid
        )  # Dune growth parameter WILL NOT come from external file
        if np.size(rmin) > 1:
            set_yaml(
                "rmin", rmin[iB3D], fid
            )  # Minimum growth rate for logistic dune growth
        else:
            set_yaml("rmin", rmin, fid)  # Minimum growth rate for logistic dune growth
        if np.size(rmax) > 1:
            set_yaml(
                "rmax", rmax[iB3D], fid
            )  # Maximum growth rate for logistic dune growth
        else:
            set_yaml("rmax", rmax, fid)  # Maximum growth rate for logistic dune growth

        # external file names used for initialization
        set_yaml("storm_file", storm_file, fid)
        if np.size(dune_file) > 1:
            set_yaml("dune_file", dune_file[iB3D], fid)
        else:
            set_yaml("dune_file", dune_file, fid)
        if np.size(elevation_file) > 1:
            set_yaml("elevation_file", elevation_file[iB3D], fid)
        else:
            set_yaml("elevation_file", elevation_file, fid)

        # the following parameters CANNOT be changed or else the MSSM storm list & storm time series needs to be remade
        set_yaml("MHW", MHW, fid)  # [m] Elevation of Mean High Water
        set_yaml("beta", beta, fid)  # Beach slope for runup calculations
        set_yaml(
            "BermEl", float(brie._h_b_crit), fid
        )  # [m] Static elevation of berm

        barrier3d.append(Barrier3d.from_yaml(datadir, prefix=parameter_file_prefix))

        # now update brie back barrier position, height of barrier, and SLRR time series with that from B3D so all
        # the initial conditions are the same! The rate of slr can only be constant in brie, whereas it can accelerate
        # in barrier3d, so by replacing the slr time series in brie we enable new functionality!
        # NOTE: interestingly here we don't need to have a "setter" in the property class for x_b, h_b, etc. because
        # we are only replacing certain indices but added for completeness
        brie.x_b[iB3D] = (
            barrier3d[iB3D].x_b_TS[0] * 10
        )  # the shoreline position + average interior width
        brie.h_b[iB3D] = (
            barrier3d[iB3D].h_b_TS[0] * 10
        )  # average height of the interior domain
        brie.x_b_save[iB3D, 0] = brie.x_b[iB3D]
        brie.h_b_save[iB3D, 0] = brie.h_b[iB3D]

    brie.slr = (
        np.array(barrier3d[0].RSLR) * 10
    )  # same for all b3d domains, just use first

    return barrier3d


class BrieCoupler:
    """Couple brie and barrier3d

    Examples
    --------
    # >>> from cascade.brie_coupler import BrieCoupler
    # >>> brie_coupler = BrieCoupler()
    # >>> brie_coupler.update_inlets(brie, barrier3d)
    # >>> brie_coupler.update_ast(brie, barrier3d)
    """

    def __init__(
        self,
        name="default",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.004,
        back_barrier_depth=3.0,
        s_background=0.001,
        h_b_crit=1.9,
        ny=1,
        nt=200,
    ):
        """The AlongshoreCoupler module.

        Parameters
        ----------
        name: string, optional
            Name of simulation
        wave_height: float, optional
            Mean offshore significant wave height [m].
        wave_period: float, optional
            Mean wave period [s].
        wave_asymmetry: float, optional
            Fraction of waves approaching from left (looking onshore).
        wave_angle_high_fraction: float, optional
            Fraction of waves approaching from angles higher than 45 degrees.
        sea_level_rise_rate: float, optional
            Rate of sea_level rise [m/yr].
        back_barrier_depth: float, optional
            Depth of back-barrier bay [m]
        ny: int, optional
            The number of alongshore Barrier3D domains for simulation in BRIE
        nt: int, optional
            Number of time steps.
        s_background: float, optional
            background slope (for shoreface toe position, back-barrier & inlet calculations)
        h_b_crit: float, optional
            critical barrier height for overwash [m], used also to calc shoreline diffusivity; we set = to B3D berm ele

        """
        ###############################################################################
        # initial conditions for BRIE
        ###############################################################################

        # parameters that we need to initialize in BRIE for coupling (not necessarily default values), but I won't be
        # modifying often (or ever) for CASCADE
        brie_ast_model = True  # shoreface formulations on
        brie_barrier_model = False  # LTA14 overwash model off
        brie_inlet_model = False  # inlet model off -- once coupled, will make this more functional
        b3d_barrier_model = True  # B3d overwash model on
        z = 10.0  # initial sea level (for tracking SL, Eulerian reference frame)

        # inlet parameters (use default; these are here to remind me later that they are important and I can change)
        Jmin = 10000  # minimum inlet spacing [m]
        a0 = 0.5  # amplitude of tide [m]
        marsh_cover = 0.5  # % of backbarrier covered by marsh and therefore does not contribute to tidal prism

        # grid and time step params
        dy = 500  # m, length of alongshore section (same as B3D) -- do not change
        dt = 1  # yr, timestep (same as B3D) -- do not change
        dtsave = 1  # save spacing (every year) -- do not change

        # start by initializing BRIE b/c it has parameters related to wave climate that we use to initialize B3D
        self._brie = Brie(
            name=name,
            ast_model=brie_ast_model,
            barrier_model=brie_barrier_model,
            inlet_model=brie_inlet_model,
            b3d=b3d_barrier_model,
            wave_height=wave_height,
            wave_period=wave_period,
            wave_asymmetry=wave_asymmetry,
            wave_angle_high_fraction=wave_angle_high_fraction,
            sea_level_rise_rate=sea_level_rise_rate,
            sea_level_initial=z,
            barrier_height_critical=h_b_crit,
            tide_amplitude=a0,
            back_barrier_marsh_fraction=marsh_cover,
            back_barrier_depth=back_barrier_depth,
            xshore_slope=s_background,
            inlet_min_spacing=Jmin,
            alongshore_section_length=dy,
            alongshore_section_count=ny,
            time_step=dt,
            time_step_count=nt,
            save_spacing=dtsave,
        )  # initialize class

    def update_ast(self, barrier3d, x_t_dt, x_s_dt, h_b_dt):
        """Pass shoreline and shoreface values from B3D subdomains to brie for use in second time step

        Parameters
        ----------
        barrier3d: list
            Barrier3D classes
        x_t_dt: list of floats
            Change in shoreface toe [m]
        x_s_dt: list of floats
            Change in shoreline position [m]
        h_b_dt: list of floats
            Change in barrier height [m]

        """
        self._brie.x_t_dt = x_t_dt  # this accounts for RSLR
        self._brie.x_s_dt = x_s_dt
        self._brie.x_b_dt = 0  # we set x_b below
        self._brie.h_b_dt = h_b_dt

        # update brie one time step (this is time_index = 2 at start of loop)
        self._brie.update()

        for iB3D in range(self._brie.ny):
            # pass shoreline position back to B3D from Brie (convert from m to dam)
            barrier3d[iB3D].x_s = self._brie.x_s[iB3D] / 10
            barrier3d[iB3D].x_s_TS[-1] = self._brie.x_s[iB3D] / 10

            # update dune domain in B3D (erode/prograde) based on shoreline change from Brie
            barrier3d[iB3D].update_dune_domain()

            # update back-barrier shoreline location in BRIE based on new shoreline + average interior width in B3D
            # NOTE: all of these positions also get updated at the end of the human management time loop
            self._brie.x_b[iB3D] = barrier3d[iB3D].x_b_TS[-1] * 10

            # this (I think) prevents a Barrier3D drowning error
            if math.isnan(self._brie.x_b[iB3D]):
                pass
            else:
                self._brie.x_b_save[iB3D, self._brie.time_index - 1] = self._brie.x_b[
                    iB3D
                ]

    def update_brie_for_human_modifications(self, x_t, x_s, x_b, h_b, s_sf):
        """Reset all of the save variables and current variables for the barrier geometry to account for human
        modifications in the cascade modules

        Parameters
        ----------
        x_t: list of floats
            Shoreface toe [dam]
        x_s: list of floats
            Shoreline [dam]
        x_b: list of floats
            Back-barrier shoreline [dam]
        h_b: list of floats
            Height of barrier [dam]
        s_sf: list of floats
            Slope of shoreface

        """
        for iB3D in range(self._brie.ny):
            self._brie.x_s[iB3D] = x_s[iB3D] * 10  # convert from dam to meters
            self._brie.x_b[iB3D] = x_b[iB3D] * 10
            self._brie.x_t[iB3D] = x_t[iB3D] * 10
            self._brie.h_b[iB3D] = h_b[iB3D] * 10

            self._brie.x_s_save[iB3D, self._brie.time_index - 1] = self._brie.x_s[iB3D]
            self._brie.x_b_save[iB3D, self._brie.time_index - 1] = self._brie.x_b[iB3D]
            self._brie.x_t_save[iB3D, self._brie.time_index - 1] = self._brie.x_t[iB3D]
            self._brie.h_b_save[iB3D, self._brie.time_index - 1] = self._brie.h_b[iB3D]
            self._brie.s_sf_save[iB3D, self._brie.time_index - 1] = s_sf[
                iB3D
            ]  # slope, so no unit change

    @property
    def brie(self):
        return self._brie

    # def update_tidal_inlets(self):
    #     # just a reminder that when we couple inlets, we're going to have to reconcile the sloping back-barrier
    #     # (vs not sloping in barrier3d) for basin_width -- maybe replace basin_width in the coupled version?
