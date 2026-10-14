"""Tests for cascade.groin: the corrected groin's schedule, its zero-strength identity and its direction.

The callback reads BRIE's state through ``cascade._brie_coupler._brie``, so a stand-in
carrying exactly the attributes it reads is enough to pin its arithmetic. Whether it
tracks a real CASCADE run is answered by full runs (strength 0 reproduces the no-groin
run to 0 m; the runner reproduces the experiment runs it was fitted with to 0 m).
"""

from types import SimpleNamespace

import numpy as np
import pytest
from scipy.stats import uniform

from cascade.groin import CorrectedGroinCallback

NY = 40
UP, DOWN = 21, 20          # updrift is the higher index, as at Buxton
DY, DT = 500.0, 1.0


def fake_cascade(x_s, diffusivity=8000.0, drift=1.0e5):
    """A stand-in with the BRIE attributes the callback reads; the drift table is uniform."""
    wave_climl = 180
    angle_array = np.linspace(-90, 90, wave_climl)
    ny = len(x_s)
    di = np.r_[ny, np.arange(2, ny + 1), np.arange(1, ny + 1), np.arange(1, ny), 1] - 1
    dj = np.r_[1, np.arange(1, ny), np.arange(1, ny + 1), np.arange(2, ny + 1), ny] - 1
    brie = SimpleNamespace(
        x_s=np.asarray(x_s, dtype=float),
        _coast_diff=np.full(wave_climl + 1, diffusivity),
        _coast_qs=np.full(wave_climl, drift),
        _angles=uniform(loc=-90, scale=180),
        _angle_array=angle_array,
        _wave_climl=wave_climl,
        _dy=DY, _dt=DT, _h_b_crit=2.0, _d_sf=10.0,
        _di=di, _dj=dj,
    )
    return SimpleNamespace(_brie_coupler=SimpleNamespace(_brie=brie))


def groin(history=((1969, 0.25), (2017, 0.0)), start=2000):
    return CorrectedGroinCallback(UP, DOWN, list(history), start_year=start, install_year=1969, n_domains=NY)


def test_strength_follows_the_history():
    g = groin()
    assert [g._effective_trapping_rate(y) for y in (1968, 1969, 2016, 2017, 2030)] == [0.0, 0.25, 0.25, 0.0, 0.0]


def test_history_must_be_shares():
    with pytest.raises(ValueError):
        groin(history=((1969, 1.5),))
    with pytest.raises(ValueError):
        groin(history=())


def test_pads_must_be_adjacent():
    with pytest.raises(ValueError):
        CorrectedGroinCallback(UP + 1, DOWN, [(1969, 0.25)], start_year=2000, install_year=1969, n_domains=NY)


def test_zero_strength_changes_nothing():
    x = np.linspace(0.0, -500.0, NY)
    cascade = fake_cascade(x)
    g = groin(history=((1969, 0.0),))
    for _ in range(3):
        out = g(cascade, np.zeros(NY))
        assert np.array_equal(out, np.zeros(NY))
    assert g.footprint_max_m_TS == [0.0, 0.0, 0.0]


def test_the_groin_moves_sand_from_downdrift_to_updrift_and_conserves_it():
    x = np.linspace(0.0, -500.0, NY)
    g = groin()
    out = g(fake_cascade(x), np.zeros(NY))
    # x_s is landward +: the side the drift comes from advances (negative), the other retreats, by the same amount
    assert out[DOWN] == -out[UP] != 0.0
    assert abs(out.sum()) < 1e-9


def test_off_after_the_end_year():
    x = np.linspace(0.0, -500.0, NY)
    g = groin(start=2017)
    out = g(fake_cascade(x), np.zeros(NY))
    assert g.blocking_applied_TS == [0.0]
    assert np.array_equal(out, np.zeros(NY))
