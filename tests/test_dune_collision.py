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
    collision_cells,
    erode_laterally,
    larson_volume,
    lower_dunes,
    refill_deficit,
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
# collision_cells
# --------------------------------------------------------------------------- #
def test_collision_cells_needs_toe_below_water_below_crest():
    r_high = 0.5
    crest_pre = np.array([0.6, 0.4, 0.6, 0.6])
    crest = np.array([0.6, 0.6, 0.5, 0.6])
    toe = np.array([0.3, 0.3, 0.3, 0.5])
    hit = collision_cells(crest_pre, crest, toe, r_high)
    # 0 in collision; 1 overwashed by Barrier3D; 2 crest now at the water; 3 water at the toe
    np.testing.assert_array_equal(hit, [True, False, False, False])


def test_collision_cells_crest_prestorm_at_water_is_not_overwashed():
    # Barrier3D overwashes where crest < Rhigh, so crest == Rhigh is left to the module
    hit = collision_cells(np.array([0.5]), np.array([0.6]), np.array([0.3]), 0.5)
    assert hit[0]


# --------------------------------------------------------------------------- #
# larson_volume
# --------------------------------------------------------------------------- #
def test_larson_volume_by_hand():
    # 1 m above the toe, Cs 1e-4, 10 h, T 10 s: 4 * 1e-4 * 1^2 * 36000 / 10 = 1.44 m3/m
    v = larson_volume(0.3, np.array([0.2]), 10.0, 10.0, 1e-4)
    np.testing.assert_allclose(v, [1.44])


def test_larson_volume_scales_with_square_of_excess():
    v = larson_volume(0.4, np.array([0.3, 0.2]), 10.0, 10.0, 1e-4)
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


def test_lower_dunes_front_by_hand():
    # 20 m3/m over the front row's 10 m lowers it 2 m: 0.4 -> 0.2 dam; the back row is unchanged
    new = lower_dunes(np.array([[0.4, 0.2]]), np.array([20.0]), np.array([0.2]), BERM, WIDTH, RESTART,
                      rows="front")
    np.testing.assert_allclose(new, [[0.2, 0.2]])


def test_lower_dunes_equal_by_hand():
    # 20 m3/m over 20 m of dune lowers every row 1 m: 0.4 -> 0.3 and 0.2 -> 0.1 dam
    new = lower_dunes(np.array([[0.4, 0.2]]), np.array([20.0]), np.array([0.1]), BERM, WIDTH, RESTART,
                      rows="equal")
    np.testing.assert_allclose(new, [[0.3, 0.1]])


@pytest.mark.parametrize("rows", ["front", "equal"])
def test_lower_dunes_front_and_equal_conserve_volume(rows):
    # away from the toe and DuneRestart, the sand taken from the rows is the eroded volume
    dune_rows = np.array([[0.6, 0.5], [0.5, 0.6]])
    volume = np.array([10.0, 30.0])  # m3/m
    new = lower_dunes(dune_rows, volume, np.full(2, 0.11), BERM, WIDTH, RESTART, rows=rows)
    taken = (dune_rows - new).sum(axis=1) * 10 * 10 * (WIDTH / 2) * 10 / 10  # dam of height x 10 m rows -> m3/m
    np.testing.assert_allclose(taken, volume)


@pytest.mark.parametrize("rows", ["front", "equal"])
def test_lower_dunes_other_modes_stop_at_the_toe_and_never_raise(rows):
    rng = np.random.default_rng(3)
    dune_rows = rng.uniform(0, 0.6, (200, 2))
    toe = rng.uniform(0.1, 0.4, 200)
    new = lower_dunes(dune_rows, rng.uniform(0, 80, 200), toe, BERM, WIDTH, RESTART, rows=rows)
    assert np.all(new <= dune_rows)
    floor = np.maximum(toe - BERM, RESTART)
    lowered = new < dune_rows
    assert np.all(new[lowered] >= np.broadcast_to(floor[:, None], new.shape)[lowered] - 1e-12)


@pytest.mark.parametrize("rows", ["front", "equal"])
def test_lower_dunes_other_modes_zero_volume_is_exact(rows):
    rng = np.random.default_rng(4)
    dune_rows = rng.uniform(0.01, 0.6, (100, 2))
    new = lower_dunes(dune_rows, np.zeros(100), np.full(100, 0.12), 0.144, WIDTH, RESTART, rows=rows)
    np.testing.assert_array_equal(new, dune_rows)


# lateral: the toe at the berm, so a row holds its height x 10 m: front 0.4 dam -> 40 m3/m, back 0.2 -> 20
LAT_ROWS, LAT_TOE = np.array([[0.4, 0.2]]), np.array([BERM])


def test_erode_laterally_keeps_the_height_until_the_front_row_is_gone():
    rows, lost = erode_laterally(LAT_ROWS, np.zeros((1, 2)), np.array([30.0]), LAT_TOE, BERM, RESTART, 1.0)
    np.testing.assert_allclose(rows, LAT_ROWS)
    np.testing.assert_allclose(lost, [[30.0, 0.0]])


def test_erode_laterally_removes_the_front_row_and_carries_the_rest_back():
    rows, lost = erode_laterally(LAT_ROWS, np.zeros((1, 2)), np.array([50.0]), LAT_TOE, BERM, RESTART, 1.0)
    np.testing.assert_allclose(rows, [[RESTART, 0.2]])
    np.testing.assert_allclose(lost, [[0.0, 10.0]])


def test_erode_laterally_carries_the_deficit_between_storms():
    rows, lost = erode_laterally(LAT_ROWS, np.zeros((1, 2)), np.array([30.0]), LAT_TOE, BERM, RESTART, 1.0)
    rows, lost = erode_laterally(rows, lost, np.array([20.0]), LAT_TOE, BERM, RESTART, 1.0)
    np.testing.assert_allclose(rows, [[RESTART, 0.2]])
    np.testing.assert_allclose(lost, [[0.0, 10.0]])


def test_erode_laterally_both_rows_gone_stops_at_the_toe():
    toe = np.array([0.15])  # 0.05 dam above the berm: the rows drop to it, the sand below it stays
    rows, lost = erode_laterally(LAT_ROWS, np.zeros((1, 2)), np.array([1e4]), toe, BERM, RESTART, 1.0)
    np.testing.assert_allclose(rows, [[0.05, 0.05]])
    np.testing.assert_allclose(lost, 0.0)


def test_erode_laterally_zero_volume_is_exact():
    rng = np.random.default_rng(5)
    dune_rows, lost = rng.uniform(0.01, 0.6, (100, 2)), rng.uniform(0, 5, (100, 2))
    rows, new_lost = erode_laterally(dune_rows, lost, np.zeros(100), np.full(100, BERM), BERM, RESTART, 1.0)
    np.testing.assert_array_equal(rows, dune_rows)
    np.testing.assert_array_equal(new_lost, lost)


def test_refill_deficit_pays_back_growth_first():
    # 0.05 dam of growth over a 10 m row is 5 m3/m: it repays 5 of the 30 lost and the row stays at 0.4
    rows, lost = refill_deficit(np.array([[0.45, 0.25]]), LAT_ROWS, np.array([[30.0, 0.0]]), 1.0)
    np.testing.assert_allclose(rows, [[0.4, 0.25]])
    np.testing.assert_allclose(lost, [[25.0, 0.0]])


def test_refill_deficit_beyond_the_deficit_raises_the_row():
    rows, lost = refill_deficit(np.array([[0.45, 0.2]]), LAT_ROWS, np.array([[2.0, 0.0]]), 1.0)
    np.testing.assert_allclose(rows, [[0.43, 0.2]])
    np.testing.assert_allclose(lost, 0.0)


def test_lateral_deficit_moves_with_the_rows_when_the_shoreline_erodes():
    b3d = new_barrier3d()
    m = DuneCollisionManager(b3d, low_toe(b3d), cs=1e-4, rows="lateral")
    m._deficit[:, 0], m._deficit[:, 1] = 7.0, 3.0
    b3d.migrate_dunes(-1.0, 1, 0)  # erosion: the front row goes, the back row moves to the front
    np.testing.assert_allclose(m.deficit[:, 0], 3.0)
    np.testing.assert_allclose(m.deficit[:, 1], 0.0)


def test_lateral_deficit_moves_with_the_rows_when_the_shoreline_progrades():
    b3d = new_barrier3d()
    m = DuneCollisionManager(b3d, low_toe(b3d), cs=1e-4, rows="lateral")
    m._deficit[:, 0], m._deficit[:, 1] = 7.0, 3.0
    b3d.migrate_dunes(1.0, 1, 0)  # progradation: a new front row, the old front moves back
    np.testing.assert_allclose(m.deficit[:, 0], 0.0)
    np.testing.assert_allclose(m.deficit[:, 1], 7.0)


def test_lateral_runs_and_keeps_a_valid_deficit():
    b3d = new_barrier3d()
    m = DuneCollisionManager(b3d, low_toe(b3d), cs=1.24e-4, rows="lateral")
    run(b3d)
    assert np.isfinite(m.deficit).all() and (m.deficit >= 0).all()
    assert sum(m.volume_TS) > 0
    pickle.loads(pickle.dumps(m))


@pytest.mark.parametrize("rows", ["proportional", "lateral"])
def test_per_cell_series_add_up(rows):
    b3d = new_barrier3d()
    m = DuneCollisionManager(b3d, low_toe(b3d), cs=1.24e-4, rows=rows)
    run(b3d)
    n = b3d.BarrierLength
    assert m.collision_volume_TS.shape == (b3d._TMAX, n)
    assert m.dune_row_loss_TS.shape == (b3d._TMAX, n, b3d._DuneWidth)
    assert (m.collision_volume_TS >= 0).all() and (m.dune_row_loss_TS >= 0).all()
    assert m.collision_storms_TS.sum() > 0
    # no volume in a cell-year without a collision storm
    assert not (m.collision_volume_TS[m.collision_storms_TS == 0] > 0).any()
    if rows == "lateral":
        # the island mean of the per-cell volume is the existing yearly series
        yearly = m.collision_volume_TS.sum(axis=1) / n
        np.testing.assert_allclose(yearly[1:1 + len(m.volume_TS)], m.volume_TS, rtol=1e-9, atol=1e-12)
        assert m.collision_deficit_TS.shape == (b3d._TMAX, n, b3d._DuneWidth)


def test_a_flattened_row_is_not_counted_again():
    # a row already at its toe, topped up by a sliver, is not a new collapse when a storm takes the sliver
    b3d = new_barrier3d()
    m = DuneCollisionManager(b3d, low_toe(b3d), cs=1e-4, rows="lateral")
    hit = np.zeros(b3d.BarrierLength, bool)
    hit[:2] = True
    floor = max(low_toe(b3d)[0] - b3d.BermEl, b3d._DuneRestart)
    old = np.array([[floor + 0.01, 0.2], [floor + 0.2, 0.2]])   # a sliver above the toe, and a real dune
    new = np.array([[floor, 0.2], [floor, 0.2]])
    m._track(1, hit, np.array([5.0, 50.0]), old, new)
    np.testing.assert_array_equal(m.dune_row_loss_TS[1, :2, 0], [0, 1])


def test_lateral_records_a_collapse_when_a_front_row_is_eroded_through():
    b3d = new_barrier3d()
    m = DuneCollisionManager(b3d, low_toe(b3d), cs=2.5e-4, rows="lateral")
    run(b3d)
    assert m.dune_row_loss_TS[..., 0].sum() > 0


def test_lower_dunes_rejects_an_unknown_row_mode():
    with pytest.raises(ValueError, match="rows must be one of"):
        lower_dunes(np.array([[0.4, 0.2]]), np.array([20.0]), np.array([0.2]), BERM, WIDTH, RESTART,
                    rows="back")


# --------------------------------------------------------------------------- #
# DuneCollisionManager: construction
# --------------------------------------------------------------------------- #
def test_manager_installs_itself():
    b3d = new_barrier3d()
    m = DuneCollisionManager(b3d, low_toe(b3d), cs=1e-4)
    assert b3d.DuneGrowth.__self__ is m
    assert b3d.DiffuseDunes.__self__ is m


def test_manager_rejects_an_unknown_row_mode():
    with pytest.raises(ValueError, match="rows must be one of"):
        DuneCollisionManager(new_barrier3d(), low_toe(new_barrier3d()), cs=1e-4, rows="back")


def test_front_mode_leaves_the_back_row_to_barrier3d(plain_run):
    # collision takes sand only from the front row, so the back row matches plain Barrier3D's until
    # the alongshore diffusion or overwash moves sand between cells; at least it is never lower
    # than with both rows lowered
    both, front = new_barrier3d(), new_barrier3d()
    DuneCollisionManager(both, low_toe(both), cs=1.24e-4)
    DuneCollisionManager(front, low_toe(front), cs=1.24e-4, rows="front")
    run(both), run(front)
    assert front.DuneDomain[:, :, 1].sum() > both.DuneDomain[:, :, 1].sum()


def test_manager_rejects_a_wrong_toe_length():
    b3d = new_barrier3d()
    with pytest.raises(ValueError, match="toe_elevation"):
        DuneCollisionManager(b3d, np.ones(b3d.BarrierLength + 1), cs=1e-4)


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

    def counting(DuneDomain, t):
        update(DuneDomain, t)
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
    hit = collision_cells(crest_pre, before.max(axis=1) + start._BermEl, m.toe_elevation, r_high)
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

