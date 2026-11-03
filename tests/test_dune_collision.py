"""Tests for cascade/dune_collision.py (DuneCollisionManager)

The functions are checked against hand-worked values; the manager on Barrier3D's default domain
(tests/test_coupling), against plain Barrier3D run from the same state and seeds.
"""

import copy
import pickle
from pathlib import Path

import numpy as np
import pytest
from barrier3d import Barrier3d

from cascade.dune_collision import (
    DuneCollisionManager,
    DuneCollisionOutsideTestedRange,
    beach_factor,
    collision_cells,
    larson_volume,
    lower_dunes,
)

DATA_DIR = Path(__file__).parent / "test_coupling"
YEARS = 12


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def new_barrier3d():
    np.random.seed(0)
    return Barrier3d.from_path(DATA_DIR, prefix="barrier3d-default")


def run(barrier3d, years=YEARS, before_update=None):
    """Advance Barrier3D, seeding each year so two runs from one state draw the same numbers."""
    for _ in range(years):
        if before_update is not None:
            before_update(barrier3d)
        np.random.seed(1000 + barrier3d.time_index)
        barrier3d.update()
        barrier3d.update_dune_domain()
    return barrier3d


def low_toe(barrier3d):
    """A toe just above the berm, so every storm above it reaches the dunes."""
    return np.full(barrier3d.BarrierLength, barrier3d.BermEl + 0.01)


def assert_same_run(a, b):
    np.testing.assert_array_equal(a.DuneDomain, b.DuneDomain)
    np.testing.assert_array_equal(a.InteriorDomain, b.InteriorDomain)
    np.testing.assert_array_equal(a.x_s_TS, b.x_s_TS)


@pytest.fixture(scope="module")
def plain_run():
    return run(new_barrier3d())


# --------------------------------------------------------------------------- #
# beach_factor
# --------------------------------------------------------------------------- #
def test_beach_factor_falls_linearly_to_zero():
    w = np.array([0.0, 2.25, 4.5, 9.0])
    np.testing.assert_allclose(beach_factor(w, 4.5), [1.0, 0.5, 0.0, 0.0])


# --------------------------------------------------------------------------- #
# collision_cells
# --------------------------------------------------------------------------- #
def test_collision_cells_needs_toe_below_water_below_crest():
    r_high = 0.5
    crest_pre = np.array([0.6, 0.4, 0.6, 0.6, 0.6, 0.6])
    crest = np.array([0.6, 0.6, 0.5, 0.6, 0.6, 0.6])
    toe = np.array([0.3, 0.3, 0.3, 0.5, 0.3, 0.3])
    beach = np.array([1.0, 1.0, 1.0, 1.0, 0.0, 0.2])
    hit = collision_cells(crest_pre, crest, toe, r_high, beach)
    # 0 in collision; 1 overwashed by Barrier3D; 2 crest now at the water; 3 water at the toe;
    # 4 behind a beach too wide; 5 in collision, part beach
    np.testing.assert_array_equal(hit, [True, False, False, False, False, True])


def test_collision_cells_crest_prestorm_at_water_is_not_overwashed():
    # Barrier3D overwashes where crest < Rhigh, so crest == Rhigh is left to the module
    hit = collision_cells(np.array([0.5]), np.array([0.6]), np.array([0.3]), 0.5, np.ones(1))
    assert hit[0]


def test_collision_cells_r_dune_sets_only_the_toe_test():
    crest = np.array([0.6, 0.6])
    toe = np.array([0.3, 0.3])
    r_dune = np.array([0.4, 0.25])  # the second falls short of its toe
    hit = collision_cells(crest, crest, toe, 0.5, np.ones(2), r_dune)
    np.testing.assert_array_equal(hit, [True, False])


# --------------------------------------------------------------------------- #
# larson_volume
# --------------------------------------------------------------------------- #
def test_larson_volume_by_hand():
    # 1 m above the toe, Cs 1e-4, 10 h, T 10 s: 4 * 1e-4 * 1^2 * 36000 / 10 = 1.44 m3/m
    v = larson_volume(0.3, np.array([0.2, 0.2]), 10.0, 10.0, 1e-4, np.array([1.0, 0.5]))
    np.testing.assert_allclose(v, [1.44, 0.72])


def test_larson_volume_scales_with_square_of_excess():
    v = larson_volume(0.4, np.array([0.3, 0.2]), 10.0, 10.0, 1e-4, np.ones(2))
    np.testing.assert_allclose(v[1] / v[0], 4.0)


# --------------------------------------------------------------------------- #
# lower_dunes
# --------------------------------------------------------------------------- #
BERM, WIDTH, RESTART = 0.1, 2.0, 0.0075


def test_lower_dunes_by_hand():
    # 20 m3/m over 20 m of dune lowers the crest 1 m: top row 0.4 -> 0.3 dam, the other row
    # in proportion, 0.2 -> 0.15
    new = lower_dunes(np.array([[0.4, 0.2]]), np.array([20.0]), np.array([0.2]), BERM, WIDTH, RESTART)
    np.testing.assert_allclose(new, [[0.3, 0.15]])


def test_lower_dunes_stops_at_the_toe():
    rows = np.array([[0.4, 0.2]])
    new = lower_dunes(rows, np.array([1e4]), np.array([0.3]), BERM, WIDTH, RESTART)
    assert new.max() + BERM == pytest.approx(0.3)


def test_lower_dunes_stops_at_dune_restart():
    # a toe below the berm would take the dune below zero height
    new = lower_dunes(np.array([[0.4, 0.2]]), np.array([1e4]), np.array([0.0]), BERM, WIDTH, RESTART)
    np.testing.assert_allclose(new, RESTART)


def test_lower_dunes_never_raises_a_row():
    rng = np.random.default_rng(1)
    rows = rng.uniform(0, 0.6, (200, 2))
    rows[:20, 1] = 0.001  # rows below DuneRestart stay as they are
    new = lower_dunes(rows, rng.uniform(0, 80, 200), rng.uniform(0.1, 0.4, 200), BERM, WIDTH, RESTART)
    assert np.all(new <= rows)


def test_lower_dunes_zero_volume_is_exact():
    rng = np.random.default_rng(2)
    rows = rng.uniform(0.01, 0.6, (500, 2))
    new = lower_dunes(rows, np.zeros(500), np.full(500, 0.12), 0.144, WIDTH, RESTART)
    np.testing.assert_array_equal(new, rows)


# --------------------------------------------------------------------------- #
# DuneCollisionManager: construction
# --------------------------------------------------------------------------- #
def test_manager_installs_itself():
    b3d = new_barrier3d()
    m = DuneCollisionManager(b3d, low_toe(b3d), cs=1e-4)
    assert b3d.DuneGrowth.__self__ is m
    assert b3d.DiffuseDunes.__self__ is m


@pytest.mark.parametrize("name", ["toe_elevation", "beach_width"])
def test_manager_rejects_wrong_lengths(name):
    b3d = new_barrier3d()
    kwds = {"toe_elevation": low_toe(b3d), "beach_width": np.ones(b3d.BarrierLength)}
    kwds[name] = np.ones(b3d.BarrierLength + 1)
    with pytest.raises(ValueError, match=name):
        DuneCollisionManager(b3d, cs=1e-4, **kwds)


@pytest.mark.parametrize(
    "option", [{"beach_follows_shoreline": True}, {"beach_runup_loss": 0.05}]
)
def test_manager_rejects_beach_options_without_a_beach(option):
    b3d = new_barrier3d()
    with pytest.raises(ValueError, match="needs a beach_width"):
        DuneCollisionManager(b3d, low_toe(b3d), cs=1e-4, **option)


# --------------------------------------------------------------------------- #
# the tested range
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("cs", [-1e-5, 2.6e-4, 1.4e-3])
def test_manager_rejects_cs_outside_the_tested_range(cs):
    b3d = new_barrier3d()
    with pytest.raises(DuneCollisionOutsideTestedRange, match="cs"):
        DuneCollisionManager(b3d, low_toe(b3d), cs=cs)


@pytest.mark.parametrize("cs", [0.0, 2.0e-5, 1.24e-4, 2.5e-4])
def test_manager_accepts_cs_inside_the_tested_range(cs):
    b3d = new_barrier3d()
    DuneCollisionManager(b3d, low_toe(b3d), cs=cs)


def test_manager_rejects_a_toe_below_mhw():
    b3d = new_barrier3d()
    toe = low_toe(b3d)
    toe[3] = -0.01
    with pytest.raises(DuneCollisionOutsideTestedRange, match="below MHW"):
        DuneCollisionManager(b3d, toe, cs=1e-4)


@pytest.mark.parametrize(
    "storm, match",
    [
        (dict(period=10.0, duration=250.0), "250 h"),
        (dict(period=25.0, duration=24.0), "25 s"),
        (dict(period=3.0, duration=24.0), "3 s"),
    ],
)
def test_storm_outside_the_tested_range_is_rejected(storm, match):
    b3d = new_barrier3d()
    m = DuneCollisionManager(b3d, low_toe(b3d), cs=1e-4)
    with pytest.raises(DuneCollisionOutsideTestedRange, match=match):
        m._check_storm(np.array([0.30]), np.array([0.20]), **storm)


def test_water_too_far_above_the_toe_is_rejected():
    b3d = new_barrier3d()
    m = DuneCollisionManager(b3d, low_toe(b3d), cs=1e-4)
    m._check_storm(np.array([0.70]), np.array([0.21]), period=10.0, duration=24.0)  # 4.9 m: inside
    with pytest.raises(DuneCollisionOutsideTestedRange, match="5.10 m above a dune toe"):
        m._check_storm(np.array([0.70]), np.array([0.19]), period=10.0, duration=24.0)


def test_a_long_storm_that_lowers_nothing_is_not_checked():
    # the storm checks apply only to storms that lower dunes: toes above every storm, storms lengthened
    b3d = new_barrier3d()
    b3d._StormSeries[:, 4] = 500.0
    toe = np.full(b3d.BarrierLength, b3d._StormSeries[:, 1].max() + 0.01)
    DuneCollisionManager(b3d, toe, cs=1e-4)
    run(b3d, years=3)


def test_a_long_storm_that_lowers_dunes_stops_the_run():
    b3d = new_barrier3d()
    b3d._StormSeries[:, 4] = 500.0
    DuneCollisionManager(b3d, low_toe(b3d), cs=1e-4)
    with pytest.raises(DuneCollisionOutsideTestedRange, match="500 h"):
        run(b3d)


# --------------------------------------------------------------------------- #
# DuneCollisionManager in Barrier3D
# --------------------------------------------------------------------------- #
def test_zero_cs_reproduces_barrier3d_exactly(plain_run):
    b3d = new_barrier3d()
    m = DuneCollisionManager(b3d, low_toe(b3d), cs=0.0)
    run(b3d)
    assert_same_run(plain_run, b3d)
    assert not any(m.volume_TS) and not any(m.crest_lowering_TS)


def test_toe_above_every_storm_reproduces_barrier3d_exactly(plain_run):
    b3d = new_barrier3d()
    m = DuneCollisionManager(b3d, np.full(b3d.BarrierLength, 10.0), cs=1e-4)
    run(b3d)
    assert_same_run(plain_run, b3d)
    assert not any(m.volume_TS)


def test_collision_lowers_dunes(plain_run):
    b3d = new_barrier3d()
    m = DuneCollisionManager(b3d, low_toe(b3d), cs=1.24e-4)
    run(b3d)
    assert sum(m.crest_lowering_TS) > 0
    assert len(m.crest_lowering_TS) == len(m.volume_TS) == YEARS
    # crest lowering is the larger row's change only, so no more than the volume (both rows)
    assert sum(m.crest_lowering_TS) <= sum(m.volume_TS) + 1e-12


def test_storms_handled_match_the_storm_series():
    b3d = new_barrier3d()
    m = DuneCollisionManager(b3d, low_toe(b3d), cs=1.24e-4)
    handled = {}
    update = m.update

    def counting(barrier3d, DuneDomain, t):
        update(barrier3d, DuneDomain, t)
        handled[t] = m._storm

    m.update = counting
    run(b3d)
    in_series = {
        t: int(np.sum(b3d._StormSeries[:, 0] == t))
        for t in range(b3d._StormStart, b3d.time_index)
    }
    assert handled == {t: n for t, n in in_series.items() if n > 0}
    assert handled  # the default series has storms in these years


def test_overwash_lowering_is_barrier3d_own(monkeypatch):
    """One storm that overwashes some cells and leaves others in collision: at the dune
    diffusion, the overwashed cells match plain Barrier3D, and only collision cells differ."""

    start = run(new_barrier3d(), 3)
    t = start.time_index

    # the pre-storm crest this year, from a dry run with no storms
    dry = copy.deepcopy(start)
    dry._StormSeries = np.empty((0, 5))
    np.random.seed(7)
    dry.update()
    crest_pre = np.maximum(dry.DuneDomain[t].max(axis=1), dry._DuneRestart) + dry._BermEl
    r_high = float(np.median(crest_pre))
    assert np.any(crest_pre < r_high) and np.any(crest_pre > r_high)

    series = np.array([[t, r_high, start._BermEl, 10.0, 24.0]])
    at_diffusion = []
    diffuse = Barrier3d.DiffuseDunes

    def capture(self, DuneDomain, time_index):
        at_diffusion.append(DuneDomain[time_index].copy())
        return diffuse(self, DuneDomain, time_index)

    monkeypatch.setattr(Barrier3d, "DiffuseDunes", capture)

    plain = copy.deepcopy(start)
    plain._StormSeries = series
    np.random.seed(7)
    plain.update()

    managed = copy.deepcopy(start)
    managed._StormSeries = series
    m = DuneCollisionManager(managed, low_toe(managed), cs=1.24e-4)
    np.random.seed(7)
    managed.update()

    before, after = at_diffusion
    overwashed = crest_pre < r_high
    np.testing.assert_array_equal(before[overwashed], after[overwashed])
    changed = np.any(before != after, axis=1)
    assert changed.any()
    hit = collision_cells(crest_pre, before.max(axis=1) + start._BermEl, m.toe_elevation, r_high, np.ones(len(crest_pre)))
    np.testing.assert_array_equal(changed, hit)


def test_manager_survives_pickling():
    # joblib (n_jobs > 1) pickles each domain every year; the manager must come back on the copy
    b3d = new_barrier3d()
    DuneCollisionManager(b3d, low_toe(b3d), cs=1.24e-4)
    copied = pickle.loads(pickle.dumps(b3d))
    m = copied.DuneGrowth.__self__
    assert copied.DiffuseDunes.__self__ is m
    assert m._barrier3d is copied
    run(copied, 6)
    assert len(m.volume_TS) == 6 and sum(m.volume_TS) > 0


def test_beach_follows_shoreline():
    b3d = new_barrier3d()
    width = np.linspace(0.0, 6.0, b3d.BarrierLength)
    m = DuneCollisionManager(
        b3d, low_toe(b3d), cs=1.24e-4, beach_width=width, beach_follows_shoreline=True
    )
    x0 = b3d.x_s_TS[-1]
    shorelines = []
    run(b3d, before_update=lambda b: shorelines.append(b.x_s_TS[-1]))
    expected = [np.maximum(width - (x - x0), 0).mean() for x in shorelines]
    np.testing.assert_allclose(m.beach_width_TS, expected)


def test_fixed_beach_keeps_its_width():
    b3d = new_barrier3d()
    m = DuneCollisionManager(b3d, low_toe(b3d), cs=1.24e-4, beach_width=np.ones(b3d.BarrierLength))
    run(b3d, 4)
    assert m.beach_width_TS == []
    np.testing.assert_array_equal(m.beach_width, 1.0)


def test_runup_loss_zero_is_no_beach():
    runs = []
    for runup in (False, True):
        b3d = new_barrier3d()
        extra = (
            {"beach_width": np.full(b3d.BarrierLength, 3.0), "beach_runup_loss": 0.0}
            if runup
            else {}
        )
        DuneCollisionManager(b3d, low_toe(b3d), cs=1.24e-4, **extra)
        runs.append(run(b3d))
    assert_same_run(*runs)


def test_runup_loss_beyond_every_toe_is_plain_barrier3d(plain_run):
    b3d = new_barrier3d()
    m = DuneCollisionManager(
        b3d,
        low_toe(b3d),
        cs=1.24e-4,
        beach_width=np.full(b3d.BarrierLength, 10.0),
        beach_runup_loss=1.0,  # the beach takes 10 dam of water level
    )
    run(b3d)
    assert_same_run(plain_run, b3d)
    assert not any(m.volume_TS)


def test_runup_loss_replaces_the_beach_factor():
    # with beach_runup_loss, a beach wider than beach_width_limit no longer stops collision
    b3d = new_barrier3d()
    m = DuneCollisionManager(
        b3d,
        low_toe(b3d),
        cs=1.24e-4,
        beach_width=np.full(b3d.BarrierLength, 10.0),
        beach_width_limit=4.5,
        beach_runup_loss=1e-4,
    )
    run(b3d)
    assert sum(m.volume_TS) > 0
