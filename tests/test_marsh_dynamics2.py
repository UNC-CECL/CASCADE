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


# i think we will want to test this for multiple years bc it looks at all the layers of years
marsh_transect = np.load(r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\marsh_transect_test.npy") * 10 # conver to m MHW
marsh_elev1 = np.zeros([2, len(marsh_transect)])
marsh_elev1[1,:] = copy.deepcopy(marsh_transect)
marsh_elev4 = np.zeros([5, len(marsh_transect)])
marsh_elev4[1,:] = copy.deepcopy(marsh_transect)
for row in range(2,5):
    marsh_elev4[row, :] = marsh_elev4[row-1,:] + 0.1  # just add 0.1 meters for each year
autoch_values_1yr = np.random.rand(2, len(marsh_transect)) * 1000
autoch_values_4yr = np.random.rand(5, len(marsh_transect)) * 1000
@pytest.mark.parametrize("yr", [1,4])
@pytest.mark.parametrize("marsh_elevations", [marsh_elev1, marsh_elev4])
@pytest.mark.parametrize("organic_dep_autoch", [autoch_values_1yr, autoch_values_4yr])  # this is going to be list one and list 2 eventually
def test_decompose(yr, marsh_elevations, organic_dep_autoch):
    marsh_transect1 = np.load(r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\marsh_transect_test.npy")
    marsh_transect2 = np.load(r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\marsh_transect_test.npy")
    x_m = 0
    model_year = yr
    organic_dep_autoch1 = copy.deepcopy(organic_dep_autoch)
    organic_dep_autoch2 = copy.deepcopy(organic_dep_autoch)
    elevation1 = copy.deepcopy(marsh_elevations)
    elevation2 = copy.deepcopy(marsh_elevations)
    mui = 0.4
    mki = 0.1
    rhoo = 85

    # marsh decomposition cascade
    elevation_casc, compaction_casc, _, organic_dep_autoch_casc = decomp_cascade(
        x_m=x_m,
        x_f=len(marsh_transect1),
        yr=model_year,
        organic_dep_autoch=organic_dep_autoch1,
        elevation=elevation1,
        B=len(marsh_transect1),
        mui=mui,
        mki=mki,
        rhoo=rhoo,
        )

    # marsh decomposition bmft
    compaction_bmft, _, organic_dep_autoch_bmft, elevation_bmft = decomp_bmft(
        x_m=x_m,
        x_f=len(marsh_transect2),
        yr=model_year,
        organic_dep_autoch=organic_dep_autoch2,
        elevation=elevation2,
        B=len(marsh_transect2),
        mui=mui,
        mki=mki,
        rhoo=rhoo,
        )
    # elevation_bmftc = marsh_transect2 - compaction_bmft

    assert_array_almost_equal(elevation_casc, elevation_bmft)
    assert_array_almost_equal(compaction_casc, compaction_bmft)

