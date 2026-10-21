"""Tests for cascade.groin: the corrected groin's schedule, its zero-strength identity, its direction and its spread.

The callback reads BRIE's state through ``cascade._brie_coupler._brie``, so a stand-in
carrying exactly the attributes it reads is enough to pin its arithmetic. Whether it
tracks a real CASCADE run is answered by full runs (strength 0 reproduces the no-groin
run to 0 m; the runner reproduces the experiment runs it was fitted with to 0 m).
"""

from types import SimpleNamespace

import numpy as np
import pytest
from scipy.stats import uniform

from cascade.groin import GroinCallback, GroinOutsideTestedRange

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
    return GroinCallback(UP, DOWN, list(history), start_year=start, n_domains=NY)


def test_strength_follows_the_history():
    g = groin()
    assert [g.strength(y) for y in (1968, 1969, 2016, 2017, 2030)] == [0.0, 0.25, 0.25, 0.0, 0.0]


def test_history_must_be_shares():
    with pytest.raises(ValueError):
        groin(history=((1969, 1.5),))
    with pytest.raises(ValueError):
        groin(history=())


def test_pads_must_be_adjacent():
    with pytest.raises(ValueError):
        GroinCallback(UP + 1, DOWN, [(1969, 0.25)], start_year=2000, n_domains=NY)


def test_zero_strength_changes_nothing():
    x = np.linspace(0.0, -500.0, NY)
    cascade = fake_cascade(x)
    g = groin(history=((1969, 0.0),))
    for _ in range(3):
        out = g(cascade, np.zeros(NY))
        assert np.array_equal(out, np.zeros(NY))
    assert [r["footprint_max_m"] for r in g.record] == [0.0, 0.0, 0.0]


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
    assert [r["strength_applied"] for r in g.record] == [0.0]
    assert np.array_equal(out, np.zeros(NY))


def test_a_face_past_89_5_degrees_does_not_index_past_the_table():
    x = np.zeros(NY)
    x[5] = 1.0e6                                  # a ~90 degree face far from the groin
    g = groin()
    cascade = fake_cascade(x)
    brie = cascade._brie_coupler._brie
    brie._coast_diff = np.full(brie._wave_climl, 8000.0)   # as long as BRIE's: _wave_climl entries
    for _ in range(2):                            # the second call re-steps last year's footprint
        g(cascade, np.zeros(NY))


def test_the_correction_zone_must_fit_inside_the_domains():
    with pytest.raises(GroinOutsideTestedRange):
        GroinCallback(3, 2, [(1969, 0.25)], start_year=2000, n_domains=NY)
    GroinCallback(10, 9, [(1969, 0.25)], start_year=2000, n_domains=NY)


def test_a_changed_wave_climate_stops_the_run():
    x = np.linspace(0.0, -500.0, NY)
    g = groin()
    g(fake_cascade(x), np.zeros(NY))
    with pytest.raises(GroinOutsideTestedRange, match="wave climate"):
        g(fake_cascade(x, diffusivity=9000.0), np.zeros(NY))


def test_a_face_steeper_than_the_tested_range_stops_the_run():
    x = np.linspace(0.0, -500.0, NY)
    x[UP:] -= 1000.0                              # the groin face at about 63 degrees
    with pytest.raises(GroinOutsideTestedRange, match="groin face"):
        groin()(fake_cascade(x), np.zeros(NY))


def test_a_groin_that_makes_or_loses_sand_stops_the_run():
    x = np.linspace(0.0, -500.0, NY)
    g = groin()
    cascade = fake_cascade(x)
    g(cascade, np.zeros(NY))
    g._footprint_m = np.zeros(NY)
    g._footprint_m[DOWN] = 2 * g.SAND_TOLERANCE_M   # a footprint that is all loss: past the tolerance
    g._last_year = None                             # keep this footprint rather than re-stepping it
    with pytest.raises(GroinOutsideTestedRange, match="lost"):
        g(cascade, np.zeros(NY))


SET_E = (0.0, 2 / 3, 1 / 3)       # Buxton's fit: none of the loss on the updrift pad, 2/3 and 1/3 beyond it


def stepped_coast(updrift, step_m=200.0):
    """A level coast with the updrift side standing step_m seaward (x_s landward +) of the downdrift side."""
    x = np.zeros(NY)
    if updrift == UP:
        x[UP:] -= step_m
    else:
        x[:DOWN + 1] -= step_m
    return x


@pytest.mark.parametrize("updrift, downdrift, away", [(UP, DOWN, 1), (DOWN, UP, -1)])
def test_the_spread_moves_the_updrift_loss_to_the_pads_beyond_and_conserves_it(updrift, downdrift, away):
    x = stepped_coast(updrift)
    plain = GroinCallback(updrift, downdrift, [(1969, 0.25)], start_year=2000, n_domains=NY)
    spread = GroinCallback(updrift, downdrift, [(1969, 0.25)], start_year=2000, n_domains=NY,
                           spread_weights=SET_E)
    difference = spread(fake_cascade(x), np.zeros(NY)) - plain(fake_cascade(x), np.zeros(NY))
    loss_m = spread.record[-1]["spread_loss_m"]
    assert loss_m > 0.0 and plain.record[-1]["spread_loss_m"] == 0.0
    # x_s is landward +: the updrift pad is given its loss back, the two pads beyond take 2/3 and 1/3 of it
    expected = np.zeros(NY)
    expected[updrift] = -loss_m
    expected[updrift + away] = 2 / 3 * loss_m
    expected[updrift + 2 * away] = 1 / 3 * loss_m
    assert np.allclose(difference, expected, atol=1e-12)
    assert abs(difference.sum()) < 1e-9


def test_the_spread_is_off_while_the_groin_is():
    x = stepped_coast(UP)
    g = GroinCallback(UP, DOWN, [(2010, 0.25)], start_year=2000, n_domains=NY, spread_weights=SET_E)
    assert np.array_equal(g(fake_cascade(x), np.zeros(NY)), np.zeros(NY))
    assert g.record[-1]["spread_loss_m"] == 0.0


@pytest.mark.parametrize("weights", [(0.5, 0.5, 0.5), (-0.1, 0.6, 0.5), (1.0, 0.0)])
def test_spread_weights_must_be_three_shares_summing_to_one(weights):
    with pytest.raises(ValueError, match="spread_weights"):
        GroinCallback(UP, DOWN, [(1969, 0.25)], start_year=2000, n_domains=NY, spread_weights=weights)
