# Lexi (Van Blunk) Fiegelist
# last updated: June 10, 2026

"""test the marsh dynamics module

Step 1. import the functions and test them with a full B3D domain
Step 2. write tests to ensure modules are functioning as we expect
Step 3. integrate into Cascade and test running cascade

"""
import copy
import numpy as np
from numpy.testing import assert_array_almost_equal

# marsh versions
from cascade.marsh_dynamics import Marsh
from cascade.marsh_dynamics import evolvemarsh as evolve_cascade
from cascade.marsh_dynamics import decompose as decomp_cascade

# BMFTC versions
from evolvemarsh import evolvemarsh as evolve_bmft
from decompose import decompose as decomp_bmft

# cascade
from cascade.cascade import Cascade

# BarrierBMFT
from barrierbmft.barrierbmft import BarrierBMFT

# -------------------------------------------------------------------------------------------------------------
# --------------------------------------- test the individual functions ---------------------------------------
# -------------------------------------------------------------------------------------------------------------
def test_evolvemarsh():
    """
    test to make sure the marsh growth function behaves the same as the barrierbmft growth function
    """
    # create deep copy of marsh transect because the function alters marsh_transect internally so we would initialize
    # the bmft with the results of the casc (I think)
    marsh_transect = np.load(r"C:\Users\agfig\model\marsh_transect_test.npy")
    marsh_transect2 = copy.deepcopy(marsh_transect)
    msl = 0
    C_e = 0.05
    OCb = 0
    numiterations = 500
    P = 12.5 * 3600 * 1
    ws = 0.05 * 10 ** (-3)
    n_tidal_cycles = 365 * (24 / 12.5)
    Bmax = 2500
    Dmin = 0
    Dmax = 1
    rhoo = 85
    rhos = 2000
    tidal_amplitude = 0.7
    tidal_range = tidal_amplitude * 2
    plot=False

    # marsh accretion cascade version
    marshelevation_casc, accretion_casc, organic_autoch_casc = evolve_cascade(
        marshelevation=marsh_transect,
        msl=msl,
        C_e=C_e,
        OCb=OCb,
        tr=tidal_range,
        numiterations=numiterations,
        P=P,
        ws=ws,
        timestep=n_tidal_cycles,
        Bmax=Bmax,
        Dmin=Dmin,
        Dmax=Dmax,
        rhoo=rhoo,
        rhos=rhos,
        plot=plot,
        )

    # marsh accretion bmft version
    marshelevation_bmft, organic_autoch_bmft, _, _, _, _, _, accretion_bmft, _ = evolve_bmft(
        marshelevation=marsh_transect2,
        msl=msl,
        C_e=C_e,
        OCb=OCb,
        tr=tidal_range,
        numiterations=numiterations,
        P=P,
        dt=P/numiterations,
        ws=ws,
        timestep=n_tidal_cycles,
        BMax=Bmax,
        Dmin=Dmin,
        Dmax=Dmax,
        rhoo=rhoo,
        rhos=rhos,
        plot=plot,
        )

    assert_array_almost_equal(marshelevation_casc, marshelevation_bmft)
    assert_array_almost_equal(accretion_casc, accretion_bmft)
    assert_array_almost_equal(organic_autoch_casc, organic_autoch_bmft)


def test_decompose():
    """
    test to make sure the marsh decomposition function behaves the same as the barrierbmft decompose function
    """
    marsh_transect = np.load(r"C:\Users\agfig\model\marsh_transect_test.npy")
    model_year = 4
    if model_year == 1:
        marsh_elev = np.zeros([2, len(marsh_transect)])
        marsh_elev[1, :] = copy.deepcopy(marsh_transect)
        autoch_values = np.random.rand(2, len(marsh_transect)) * 1000
        autoch_values_initial = copy.deepcopy(autoch_values)
    elif model_year == 4:
        marsh_elev = np.zeros([5, len(marsh_transect)])
        marsh_elev[1, :] = copy.deepcopy(marsh_transect)
        for row in range(2,5):
            marsh_elev[row, :] = marsh_elev[row-1,:] + 0.1  # just add 0.1 meters for each year
        autoch_values = np.random.rand(5, len(marsh_transect)) * 1000
        autoch_values_initial = copy.deepcopy(autoch_values)
    x_m = 0
    mui = 0.4
    mki = 0.1
    rhoo = 85

    # marsh decomposition cascade
    marsh_transect, compaction_casc, _, organic_dep_autoch_casc = decomp_cascade(
        x_m=x_m,
        x_f=len(marsh_transect),
        yr=model_year,
        organic_dep_autoch=autoch_values,
        elevation=marsh_elev,
        B=len(marsh_transect),
        mui=mui,
        mki=mki,
        rhoo=rhoo,
        )
    elevation_casc = copy.deepcopy(marsh_elev)

    # DO NOT RE-GENERATE AUTOCH VALUES, THEY MUST BE THE SAME AS ABOVE, SO USE THE autoch_values_initial VARIABLE
    marsh_transect = np.load(r"C:\Users\agfig\model\marsh_transect_test.npy")
    if model_year == 1:
        marsh_elev = np.zeros([2, len(marsh_transect)])
        marsh_elev[1, :] = copy.deepcopy(marsh_transect)
    elif model_year == 4:
        marsh_elev = np.zeros([5, len(marsh_transect)])
        marsh_elev[1, :] = copy.deepcopy(marsh_transect)
        for row in range(2,5):
            marsh_elev[row, :] = marsh_elev[row-1,:] + 0.1  # just add 0.1 meters for each year

    # marsh decomposition bmft
    compaction_bmft, _, organic_dep_autoch_bmft = decomp_bmft(
        x_m=x_m,
        x_f=len(marsh_transect),
        yr=model_year,
        organic_dep_autoch=autoch_values_initial,
        elevation=marsh_elev,
        B=len(marsh_transect),
        mui=mui,
        mki=mki,
        rhoo=rhoo,
        )
    marsh_elev[model_year] = marsh_elev[model_year] - compaction_bmft
    elevation_bmft = marsh_elev

    assert_array_almost_equal(elevation_casc, elevation_bmft)
    assert_array_almost_equal(compaction_casc, compaction_bmft)


# ----------------------------------------------------------------------------------------------------
# --------------------------------------- test the marsh class ---------------------------------------
# ----------------------------------------------------------------------------------------------------
# load a test domain and transect which is just an arbitrary segment of Masonboro Island
test_domain = np.load(r"C:\Users\agfig\model\marsh_domain_test.npy")
n_cols = np.shape(test_domain)[1]
initial_width=np.shape(test_domain)[0]
def run_marsh_dynamics(model_domain=test_domain, cols=n_cols, time=5, width=initial_width):

    # initialize the class
    marsh_class = Marsh(
        RSLR=0.012,  # this will get replaced with cascade: self._sea_level_rise_rate which is in m/yr
        C_e=0.05,
        OCb=0,
        numiterations=500,
        P=12.5 * 3600 * 1,
        ws=0.05 * 10 ** (-3),
        n_tidal_cycles=365 * (24 / 12.5),
        Bmax=2500,
        Dmin=0,
        Dmax=1,
        rhoo=85,
        rhos=2000,
        mui=0.4,
        mki=0.1,
        m_min=-0.3,
        m_max=0,
        time_step_count=time,
        alongshore_length=cols,
        tidal_amplitude=0.7,
        initial_width=width
        )

    # run the time loop/update function
    for time_step in range(marsh_class._nt):
        print("\r", "Time Step: ", time_step + 1, end="")
        marsh_class.update(
            interior_domain=model_domain,  # this will be the b3d interior domain at the time step
            model_year=time_step,
            shoreline_changeTS=np.zeros(marsh_class._nt)  # there is no shoreline change since there is no cascade/b3d
            )
        model_domain = marsh_class._marsh_elevation[time_step]

    return marsh_class, model_domain


def test_domain_orientation():
    """
    we have to convert the domain from dam MHW to m MSL which requires many flips and edits
    this test checks to make sure that the input domain orientation matches the output domain orientation
    """
    input_domain_values = np.linspace(0.5, 1, 5)  # dam MHW (high enough values that no accretion/erosion occurs)
    input_domain = np.zeros([2, 5])
    input_domain[0, :] = input_domain_values
    input_domain[1, :] = input_domain_values + 0.5
    n_cols = np.shape(input_domain)[1]
    init_width = np.shape(input_domain)[0]
    output_domain = copy.deepcopy(input_domain)  # what we expect it to be

    model_duration = 1

    # function includes class initialization and time loop
    marsh_class, model_domain = run_marsh_dynamics(
        model_domain=input_domain,
        cols=n_cols,
        time=model_duration,
        width=init_width,
        )

    assert_array_almost_equal(output_domain, model_domain)

def test_marsh_accretion(domain=test_domain, cols=n_cols, width=initial_width):
    """
    this test ensures that all cells above the water line never accrete AND therefore the elevation of those cells
    should not change over time either
    """

    model_duration = 100
    elev_limit = 0.2
    high_cells = domain[domain > elev_limit]

    # function includes class initialization and time loop
    marsh_class, model_domain = run_marsh_dynamics(
        model_domain=domain,
        cols=cols,
        time=model_duration,
        width=width
        )

    # are the high cells unchanged?
    high_cells_post_update = model_domain[model_domain > elev_limit]
    assert_array_almost_equal(high_cells, high_cells_post_update, decimal=5)

    # loop through all accretion arrays and make sure high cells are approximately 0
    for t in range(len(marsh_class._accretion_TS)):
        accretion_t = marsh_class._accretion_TS[t]
        compaction_t = marsh_class._compaction_TS[t]
        accretion_test = accretion_t[model_domain > elev_limit]
        compaction_test = compaction_t[model_domain > elev_limit]
        zero_array = np.zeros(np.shape(accretion_test))
        assert_array_almost_equal(accretion_test, zero_array, decimal=5)
        assert_array_almost_equal(compaction_test, zero_array, decimal=5)  # mostly curious if this is 0


def test_strat_compaction():
    """
    I ran two barrierbmft models, one with compaction caluclated normally (marsh_elev_full_strat) and one which only
    uses the top layer to calculate compaction (marsh_elev_no_strat)
    I just wanted to see how/if they differ
    """
    # full_strat = np.load(r"C:\Users\Lexi\Documents\UNC\model\marsh_elevation_full_strat.npy")
    # edited_strat = np.load(r"C:\Users\Lexi\Documents\UNC\model\marsh_elevation_no_strat_5yrs.npy")
    # assert_array_almost_equal(full_strat, edited_strat, decimal=5)
    b3d1 = np.load(r"C:\Users\agfig\model\barrierbmft\b3d_barrierbmft1.npz", allow_pickle=True)["cascade"][0]
    b3d2 = np.load(r"C:\Users\agfig\model\barrierbmft\b3d_barrierbmft2.npz", allow_pickle=True)["cascade"][0]
    elev1 = b3d1.InteriorDomain
    elev2 = b3d2.InteriorDomain
    if np.array_equal(elev1, elev2):
        print(
            "decomposition with a stratigraphic record produces very similar results than only decomposing the top layer")
    else:
        print(
            "decomposition with a stratigraphic record produces different results than only decomposing the top layer")

# ----------------------------------------------------------------------------------------------------------
# --------------------------------------- test CASCADE incorporation ---------------------------------------
# ----------------------------------------------------------------------------------------------------------
def test_initialization():
    """
    test to see if overlapping variables are initialized the same way
    variables we care about:
    - msl
    - RSLR
    - duration
    - tmax
    - SL
    """

    model_duration = 10
    slr_m_yr = 0.04
    slr_mm_yr = slr_m_yr * 1000

    # cascade init
    casc = Cascade(
        datadir=r"C:\Users\agfig\model\basic_cascade_run",
        time_step_count=model_duration,
        sea_level_rise_rate=slr_m_yr,
        parameter_file="barrier3d-default-parameters.yaml",
        marsh_module=True,
        alongshore_section_count=1
        )
    casc_RSLR = casc._sea_level_rise_rate   # m/yr
    casc_b3d_RSLR = casc.barrier3d[0].RSLR  # dam/yr for barrier3d variables
    casc_marsh_RSLR = casc.marsh[0]._RSLR   # m/yr but only one value, not an array. used to create marsh msl array
    casc_marsh_msl = casc.marsh[0]._msl     # m
    casc_b3d_tmax = casc.barrier3d[0].TMAX
    casc_marsh_tmax = casc.marsh[0]._nt
    casc_dur = casc._nt
    casc_b3d_SL = casc.barrier3d[0].SL

    # barrierbmft init
    barrierbmft = BarrierBMFT(
        time_step_count=model_duration,
        relative_sea_level_rise=slr_mm_yr,  # mm/yr
        parameter_file="barrier3d-parameters.yaml",
        )
    # barrierbmft creates historic data, so we only want to look at the sizes and values of the model duration
    start = barrierbmft.bmftc.startyear
    bmft_b3d_RSLR = barrierbmft.barrier3d.model.RSLR[0:-1]  # dam/yr, RSLR is initialized to dur + 1 for some reason...
    bmft_msl = barrierbmft.bmftc.msl[start:]                # m
    bmft_b3d_tmax = barrierbmft.barrier3d.model.TMAX
    bmft_dur = barrierbmft.bmftc.dur
    bmft_b3d_SL = barrierbmft.barrier3d.model.SL

    # tests
    assert_array_almost_equal(casc_RSLR, casc_marsh_RSLR)      # check that marsh is same as cascade
    assert_array_almost_equal(casc_b3d_RSLR, bmft_b3d_RSLR)    # check that b3ds are equal
    assert_array_almost_equal(casc_marsh_msl, bmft_msl)        # check that marsh msl is initialized the same as bmft msl
    assert_array_almost_equal(casc_marsh_tmax, bmft_b3d_tmax)  # check to make sure durations are the same
    assert_array_almost_equal(casc_b3d_tmax, bmft_b3d_tmax)    # check to make sure durations are the same
    assert_array_almost_equal(casc_dur, bmft_dur)              # check to make sure durations are the same
    assert_array_almost_equal(casc_b3d_SL, bmft_b3d_SL)        # check to make sure sea levels are the same (always 0)

def test_casc_run():
    """
    test to make sure the marsh is doing something!
    run cascade with and without marsh and compare final elevations
    """

    datadir = r"C:\Users\agfig\model\basic_cascade_run\marsh"
    min_dune_r = 0.05
    max_dune_r = 0.45
    beach_slope = 0.06
    model_duration = 3
    slr_m_yr = 0.04

    overwash_storm = "cascade-default-storms.npy"

    # marsh class on
    cascade_marsh = Cascade(
        datadir,
        name="marsh_run",
        elevation_file=f"marsh_domain_test3.npy",
        dune_file=f"marsh_dunes_test3.npy",
        parameter_file="marsh-default-parameters.yaml",
        storm_file=overwash_storm,
        num_cores=1,  # cascade can run in parallel, can never specify more cores than that
        roadway_management_module=False,
        alongshore_transport_module=False,
        beach_nourishment_module=False,
        community_economics_module=False,
        outwash_module=False,
        marsh_module=True,
        alongshore_section_count=1,
        time_step_count=model_duration,
        # ---------- for BRIE and Barrier3D --------------- #
        beta=beach_slope,
        sea_level_rise_rate=slr_m_yr,
        min_dune_growth_rate=min_dune_r,
        max_dune_growth_rate=max_dune_r,
        # --------------- marsh dynamics -----------------------------------
        SSCb=0.05,
        organic_content_bay=0,
        numiterations=500,
        tidal_period=12.5,
        settling_velocity=0.05 * 10 ** (-3),
        max_biomass=2500,
        min_depth_marsh_growth=0,
        max_depth_marsh_growth=0.4,
        density_organic_matter=85,
        density_sediment=2000,
        max_depth_decomp=0.4,
        decomp_coeff=0.1,
        min_elev_marsh=-0.3,
        max_elev_marsh=0,
        tidal_amplitude=0.7,

        )

    # run the time loop/update function
    for time_step in range(cascade_marsh._nt - 1):
        print("\r", "Time Step: ", time_step + 1, end="")
        cascade_marsh.update()
        if cascade_marsh.b3d_break:
            break


    # run normal - marsh class off
    cascade = Cascade(
        datadir,
        name="overwash_only",
        elevation_file=f"marsh_domain_test3.npy",
        dune_file=f"marsh_dunes_test3.npy",
        parameter_file="marsh-default-parameters.yaml",
        storm_file=overwash_storm,
        num_cores=1,  # cascade can run in parallel, can never specify more cores than that
        roadway_management_module=False,
        alongshore_transport_module=False,
        beach_nourishment_module=False,
        community_economics_module=False,
        outwash_module=False,
        marsh_module=False,
        alongshore_section_count=1,
        time_step_count=model_duration,
        # ---------- for BRIE and Barrier3D --------------- #
        beta=beach_slope,
        sea_level_rise_rate=slr_m_yr,
        min_dune_growth_rate=min_dune_r,
        max_dune_growth_rate=max_dune_r,
        )

    # run the time loop/update function
    for time_step in range(cascade._nt - 1):
        print("\r", "Time Step: ", time_step + 1, end="")
        cascade.update()
        if cascade.b3d_break:
            break

    # final domains should be different!
    casc_marsh_domain = cascade_marsh.barrier3d[0].InteriorDomain
    casc_domain = cascade.barrier3d[0].InteriorDomain
    dif = casc_marsh_domain - casc_domain
    # they should be different!!
    assert not np.array_equal(casc_marsh_domain, casc_domain)
    assert np.any(dif)  # dif should have non-zero values

def test_casc_accretion_erosion():
    """
    test to make sure we are still changing the correct cells
    only cells identified as marsh should be changing elevation after applying marsh dynamics
    """

    datadir = r"C:\Users\agfig\model\basic_cascade_run\marsh"
    min_dune_r = 0.05
    max_dune_r = 0.45
    beach_slope = 0.06
    model_duration = 15
    slr_m_yr = 0.004
    berm_elev = 1.95  # m NAVD88
    MHW = 0.421

    overwash_storm = "cascade-default-storms.npy"

    # marsh class on
    cascade_marsh = Cascade(
        datadir,
        name="marsh_run",
        elevation_file=f"marsh_domain_test2_short.npy",
        dune_file=f"marsh_dunes_test3.npy",
        parameter_file="marsh-default-parameters.yaml",
        storm_file=overwash_storm,
        num_cores=1,  # cascade can run in parallel, can never specify more cores than that
        roadway_management_module=False,
        alongshore_transport_module=False,
        beach_nourishment_module=False,
        community_economics_module=False,
        outwash_module=False,
        marsh_module=True,
        alongshore_section_count=1,
        time_step_count=model_duration,
        # ---------- for BRIE and Barrier3D --------------- #
        beta=beach_slope,
        sea_level_rise_rate=slr_m_yr,
        min_dune_growth_rate=min_dune_r,
        max_dune_growth_rate=max_dune_r,
        berm_elevation=berm_elev,
        MHW=MHW,
        # --------------- marsh dynamics -----------------------------------
        SSCb=0.05,
        organic_content_bay=0,
        numiterations=500,
        tidal_period=12.5,
        settling_velocity=0.05 * 10 ** (-3),
        max_biomass=2500,
        min_depth_marsh_growth=0,
        max_depth_marsh_growth=0.4,
        density_organic_matter=85,
        density_sediment=2000,
        max_depth_decomp=0.4,
        decomp_coeff=0.1,
        min_elev_marsh=-0.3,
        max_elev_marsh=0,
        tidal_amplitude=0.7,
        )

    # run the time loop/update function
    for time_step in range(cascade_marsh._nt - 1):
        print("\r", "Time Step: ", time_step + 1, end="")
        cascade_marsh.update()
        if cascade_marsh.b3d_break:
            break

    # min and max marsh elevations in dam MHW
    m_max_msl = cascade_marsh.marsh[0]._m_max
    m_min_msl = cascade_marsh.marsh[0]._m_min
    # marsh class stores pre and post marsh elevations in dam MHW
    pre_elev = cascade_marsh.marsh[0]._pre_marsh_elev
    post_elev = cascade_marsh.marsh[0]._marsh_elevation
    for t in range(1, model_duration):
        non_marsh_cells = np.where((pre_elev[t] > m_max_msl) | (pre_elev[t] < m_min_msl))  # | is or operator
        assert_array_almost_equal(post_elev[t][non_marsh_cells], pre_elev[t][non_marsh_cells])  # actual, desired


