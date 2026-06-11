# Lexi (Van Blunk) Fiegelist
# last updated: June 10, 2026

"""test the marsh dynamics module

Step 1. import the functions and test them with a full B3D domain
Step 2. write tests to ensure modules are functioning as we expect
Step 3. integrate into Cascade and test running cascade

"""
import copy
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

# marsh versions
from cascade.marsh_dynamics import Marsh
from cascade.marsh_dynamics import evolvemarsh as evolve_cascade
from cascade.marsh_dynamics import decompose2 as decomp_cascade

# BMFTC versions
from evolvemarsh import evolvemarsh as evolve_bmft
from decompose import decompose as decomp_bmft

# load a test domain and transect which is just an arbitrary segment of Masonboro Island
test_domain = np.load(r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\marsh_domain_test.npy")
n_cols = np.shape(test_domain)[1]

# variables used in tests
RSLR = 0.012  # m/yr
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
mui = 0.4
mki = 0.1
m_min = -0.3  # dam MHW
m_max = 0  # dam MHW
time_step_count = 20
alongshore_length = n_cols
tidal_amplitude = 0.7
tidal_range = tidal_amplitude * 2

def run_marsh_dynamics(model_domain=test_domain, cols=n_cols):
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
        time_step_count=20,
        alongshore_length=cols,
        tidal_amplitude=0.7,
        )

    # run the time loop/update function
    for time_step in range(marsh_class._nt - 1):
        print("\r", "Time Step: ", time_step + 1, end="")
        marsh_class.update(
            b3d=None,
            interior_domain=model_domain,  # this will be the b3d interior domain at the time step
            model_year=time_step
            )
        model_domain = marsh_class._marsh_elevation[time_step]

    return marsh_class


def test_evolvemarsh():
    # create deep copy of marsh transect because the function alters marsh_transect internally so we would initialize
    # the bmft with the results of the casc (I think)
    marsh_transect = np.load(r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\marsh_transect_test.npy")
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
    marsh_transect = np.load(r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\marsh_transect_test.npy")
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
    elevation_casc, compaction_casc, _, organic_dep_autoch_casc = decomp_cascade(
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

    # DO NOT RE-GENERATE AUTOCH VALUES, THEY MUST BE THE SAME AS ABOVE, SO USE THE autoch_values_initial VARIABLE
    marsh_transect = np.load(r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\marsh_transect_test.npy")
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

