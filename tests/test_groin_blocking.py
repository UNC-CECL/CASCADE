"""Tests for cascade.groin: the blocking groin, and the shared schedule.

The blocking groin reads BRIE's state through ``cascade._brie_coupler._brie``,
so a minimal stand-in carrying exactly the attributes it reads is enough to
pin its arithmetic. Whether the callback form tracks a real CASCADE run is an
integration question answered by a full run (b = 0 must reproduce the
no-groin run), not here.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from cascade.groin import BlockingGroinCallback, GroinCallback

NY = 40
UP, DOWN = 21, 20          # updrift is the higher index, as at Buxton
DY, DT = 500.0, 1.0


def fake_cascade(x_s, diffusivity=8000.0):
    """A stand-in exposing the BRIE attributes the callback reads."""
    brie = SimpleNamespace(
        x_s=np.asarray(x_s, dtype=float),
        _coast_diff=np.full(181, diffusivity),
        _wave_climl=180,
        _dy=DY,
        _dt=DT,
    )
    return SimpleNamespace(_brie_coupler=SimpleNamespace(_brie=brie))


def stepped_shoreline(step=-250.0):
    """Straight coast with a step at the DOWN|UP face (UP seaward if step < 0)."""
    x = np.zeros(NY)
    x[UP:] += step
    return x


def blocking(b=0.6, **kw):
    args = dict(updrift_pad=UP, downdrift_pad=DOWN, blocking_fraction=b,
                start_year=1996, install_year=1969, n_domains=NY)
    args.update(kw)
    return BlockingGroinCallback(**args)


def test_zero_blocking_leaves_x_s_dt_untouched():
    cb = blocking(b=0.0)
    x_s_dt = [0.1 * i for i in range(NY)]
    before = list(x_s_dt)
    assert cb(fake_cascade(stepped_shoreline()), x_s_dt) == before


def test_inert_before_install_year():
    cb = blocking(b=0.8, start_year=1960, install_year=1969)
    x_s_dt = [0.0] * NY
    assert cb(fake_cascade(stepped_shoreline()), x_s_dt) == [0.0] * NY
    assert cb.active_TS == [False]


def test_cancels_b_of_the_explicit_face_transfer():
    b, D = 0.6, 8000.0
    x = stepped_shoreline(-250.0)
    cb = blocking(b=b)
    x_s_dt = cb(fake_cascade(x, D), [0.0] * NY)

    r = D * DT / 2 / DY ** 2
    theta_down = np.degrees(np.arctan2(x[UP] - x[DOWN], DY))
    assert theta_down < 0     # the step leans the face, as at Cape Hatteras
    offset = x[UP] - x[DOWN]
    assert x_s_dt[DOWN] == pytest.approx(-b * 2 * r * offset)
    assert x_s_dt[UP] == pytest.approx(b * 2 * r * offset)
    # Only the two flanking cells change.
    assert all(v == 0.0 for i, v in enumerate(x_s_dt) if i not in (UP, DOWN))


def test_holds_the_step_rather_than_closing_it():
    # Diffusion across the face would move DOWN seaward and UP landward,
    # closing the step. The groin must push the other way.
    x = stepped_shoreline(-250.0)            # UP sits 250 m seaward of DOWN
    x_s_dt = blocking(b=0.6)(fake_cascade(x), [0.0] * NY)
    assert x_s_dt[DOWN] > 0                  # DOWN held landward
    assert x_s_dt[UP] < 0                    # UP held seaward


def test_diffusivity_is_brie_clipped_at_zero():
    cb = blocking(b=1.0)
    x_s_dt = cb(fake_cascade(stepped_shoreline(), diffusivity=-119.0), [0.0] * NY)
    assert x_s_dt == [0.0] * NY
    assert cb.r_ipl_updrift_TS == [0.0] and cb.r_ipl_downdrift_TS == [0.0]


def test_instant_failure_schedule():
    cb = blocking(b=0.6, deterioration_mode="instant",
                  deterioration_delay_years=2004 - 1969,
                  deterioration_fraction=0.3)
    assert cb._effective_trapping_rate(2003) == 0.6
    assert cb._effective_trapping_rate(2004) == pytest.approx(0.18)
    assert cb._effective_trapping_rate(2020) == pytest.approx(0.18)


def test_diagnostics_record_equivalent_trapping_rate():
    cb = blocking(b=0.5)
    cas = fake_cascade(stepped_shoreline())
    for _ in range(3):
        cb(cas, [0.0] * NY)
    rows = cb.diagnostics_frame()
    assert [r["model_year"] for r in rows] == [1996, 1997, 1998]
    assert rows[0]["trapping_rate_applied_m_yr"] == pytest.approx(
        abs(rows[0]["applied_dx_updrift_m"]))
    assert cb.summary()["kind"] == "blocking"
    assert cb.mean_trapping_rate_m_yr > 0


@pytest.mark.parametrize("bad", [
    dict(downdrift_pad=UP - 2),                    # not adjacent
    dict(blocking_fraction=1.2),
    dict(deterioration_mode="instant", deterioration_ramp_years=7.0),
])
def test_rejects_bad_configuration(bad):
    with pytest.raises(ValueError):
        blocking(**bad)


def test_dipole_schedule_arithmetic_unchanged():
    # The schedule moved into a shared helper; the dipole must still compute
    # M - taper * (M - M*f) exactly as before.
    kw = dict(updrift_pad=UP, downdrift_pad=DOWN, trapping_rate_m_yr=60.0,
              start_year=1984, install_year=1969, n_domains=NY,
              deterioration_delay_years=27, deterioration_mode="linear_ramp",
              deterioration_ramp_years=7.0, deterioration_fraction=0.6)
    cb = GroinCallback(**kw)
    for year in range(1984, 2025):
        if year < 1996:
            expected = 60.0
        else:
            taper = min(1.0, (year - 1996) / 7.0)
            expected = 60.0 - taper * (60.0 - 60.0 * 0.6)
        assert cb._effective_trapping_rate(year) == expected
    assert cb.kind == "dipole"
