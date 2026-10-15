#!/usr/bin/env python3
"""
Try to break the runner's groin: CorrectedGroinCallback under forcing it was never fitted to.

    python stress_test.py   ->  tables/stress_runs.csv, tables/stress_by_year.csv.gz

The runner's own callback (cascade.groin.CorrectedGroinCallback, imported unchanged) on the
straight-coast rig of ../2-straight-coast/1-run/: 61 domains of 500 m, one groin at the middle
face, BRIE's implicit alongshore step alone on an infinitely long tilted coast. The callback reads a
stand-in for BRIE that carries a real Brie instance's tables and the rig's shoreline, as
CASCADE hands it x_s_dt before BRIE's solve. Every run is paired with the same run with the groin
off, and scored on the difference. Families, each set before the runs (README.md):

    0 rig check       f 0.6, 5-cell zone: must match version 6 of 1-run/ on a straight coast
    1 waves           one wave parameter at a time, far outside option A
    2 strength        0 to 1; 0 must be the no-groin run exactly
    3 notch           a sharp step already at the groin when it starts (Buxton's is 256 m)
    4 switch-off      on for 50 years, then off for 100
    5 cycling         on and off every 5 years
    6 variable waves  a new wave climate each year, drift reversing in some; the callback as written
                      (drift table read once) and with its drift table refreshed each year
    7 background      cross-shore change in x_s_dt, as Barrier3D passes it: steady retreat,
                      storm-like noise per domain, a storm pulse every 10 years
    8 placement       the groin near the end of the array
    9 long run        500 years

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-15
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path[:0] = [str(HERE.parent / "2-straight-coast" / "1-run"), str(REPO)]
import straight_coast_test as sc  # noqa: E402
from brie import Brie  # noqa: E402
from cascade.groin import CorrectedGroinCallback, GroinOutsideTestedRange  # noqa: E402

DIFFUSION_NUMBER = CorrectedGroinCallback._diffusion_number   # BRIE's r at each face, as the callback reads it

# --- CONFIG ------------------------------------------------------------------
NY, DY_M, DT_YR = sc.NY, sc.DY_M, sc.DT_YR      # 61 domains of 500 m, 1-year step
LO = sc.LO                                      # the groin is the face LO|LO+1
OPTION_A = dict(Hs=2.0, Tp=7.5, asym=0.6, ahf=0.5)
HALF_WIDTH = 10                                 # the runner's correction zone
STRENGTHS = (0.25, 1.0)                         # Buxton's whole-life strength, and all of the drift
ORIENTATIONS = (-30, -20, -10, 0, 10, 20, 30)
YEARS = 100
WAVE_SWEEP = dict(ahf=(0.1, 0.3, 0.5, 0.6, 0.7), asym=(0.5, 0.55, 0.6, 0.7, 0.8, 0.9),
                  Hs=(1.0, 1.5, 2.0, 2.5, 3.0, 3.5))
SWEEP_STRENGTHS = (0.0, 0.1, 0.25, 0.5, 0.75, 1.0)
NOTCHES_M = (-300, -150, -50, 50, 150, 300)     # starting gap, updrift seaward of downdrift +
NOTCH_DECAY_CELLS = 3.0
VARIABLE_ASYM = np.round(np.arange(0.40, 0.801, 0.05), 2)   # below 0.5 the drift reverses
VARIABLE_HS = (1.5, 2.0, 2.5)
VARIABLE_SEEDS = (1, 2, 3)
RETREAT_M_YR = 2.0                              # steady landward retreat everywhere
NOISE_SD_M = 10.0                               # storm-like change per domain per year
PULSE_M, PULSE_EVERY = 30.0, 10                 # a storm cutting every domain back
PLACEMENTS = (0, 1, 2, 5, 9, 30, 58, 59)        # LO of the face LO|LO+1
LONG_YEARS = 500
# The pass rules, set before the runs
CONSERVE_ABS_M, CONSERVE_FRAC = 5.0, 0.05      # |net| <= max(5 m, 5% of gross), both against the no-groin run
MAX_ANGLE_DEG = 60.0                            # steepest face angle, below where BRIE's table runs out (90)
STEADY_FRAC, STEADY_ABS_M = 0.10, 5.0           # gap change over the last 20% <= 10% of |gap| + 5 m
TABLES = HERE / "tables"
# -----------------------------------------------------------------------------

_BRIE = {}


# The runner's callback, its tested-range guard recording when it would stop the run instead of stopping it
class Recording(CorrectedGroinCallback):
    def __init__(self, *args, **kw):
        self.guard_year, self.guard_reason, self.trace = None, None, []
        super().__init__(*args, **kw)

    def _check_zone(self, *args):
        try:
            super()._check_zone(*args)
        except GroinOutsideTestedRange as e:
            self.guard_year, self.guard_reason = 0, str(e).split(";")[0]

    def _check_tested_range(self, brie, x_s, near):
        self.trace.append((float(np.sum(self._a)), float(np.sum(np.abs(self._a))),
                           float(abs(np.degrees(np.arctan2(x_s[self._hi] - x_s[self._lo], brie._dy))))))
        try:
            super()._check_tested_range(brie, x_s, near)
        except GroinOutsideTestedRange as e:
            if self.guard_year is None:
                self.guard_year, self.guard_reason = self.start_year + self._call_count - 1, str(e).split(";")[0]


# A stand-in for BRIE with a real instance's tables, one per wave climate
def brie_like(climate):
    key = tuple(sorted(climate.items()))
    if key not in _BRIE:
        b = Brie(ast_model=True, barrier_model=False, inlet_model=False, b3d=True,
                 wave_height=climate["Hs"], wave_period=climate["Tp"],
                 wave_asymmetry=climate["asym"], wave_angle_high_fraction=climate["ahf"],
                 alongshore_section_length=DY_M, alongshore_section_count=NY,
                 time_step=DT_YR, time_step_count=10, save_spacing=1)
        _BRIE[key] = SimpleNamespace(
            _dy=b._dy, _dt=b._dt, _wave_climl=b._wave_climl, _coast_diff=b._coast_diff.copy(),
            _di=b._di.copy(), _dj=b._dj.copy(), _h_b_crit=b._h_b_crit, _d_sf=b._d_sf,
            _angles=b._angles, _angle_array=b._angle_array, _coast_qs=b._coast_qs.copy(), x_s=None)
    return _BRIE[key]


# BRIE's implicit alongshore step (brie.py:1293-1326) on the infinite tilted coast: eta is the departure from the tilt
def brie_step(b, slope, eta, dx):
    step = slope * DY_M + (np.roll(eta, -1) - eta)
    theta = 180.0 * np.arctan2(step, DY_M) / np.pi
    k = np.maximum(1, np.minimum(b._wave_climl, np.round(90 - theta).astype(int)))
    r = np.maximum(0.0, b._coast_diff[k] * DT_YR / 2 / DY_M ** 2)
    dv = np.r_[-r[-1], -r[1:], 1 + 2 * r, -r[0:-1], -r[0]]
    rhs = eta + r * (np.roll(eta, -1) - 2 * eta + np.roll(eta, 1)) + dx
    return spsolve(csr_matrix((dv, (b._di, b._dj))), rhs), theta


# Net drift toward higher index on the undisturbed coast, m3/yr, as the callback reads it
def drift0(b, theta0):
    return float(CorrectedGroinCallback._drift_table(b)[int(np.clip(round(90 - theta0), 1, b._wave_climl - 1))])


# One run: the callback in front of BRIE's step every year; returns one row per year
def run(theta0, years=YEARS, strength=0.25, history=None, climate=OPTION_A, climate_of_year=None,
        refresh_drift=False, half_width=HALF_WIDTH, lo=LO, notch_m=0.0, background=None, seed=0,
        groin=True):
    rng = np.random.default_rng(seed)
    slope = np.tan(np.deg2rad(theta0))
    origin = (lo + 1 + NY // 2) % NY             # where the periodic tilt jumps: opposite the groin
    tilt = slope * DY_M * ((np.arange(NY) - origin) % NY)
    b0 = brie_like(climate)
    up_hi = drift0(b0, theta0) > 0               # drift toward higher index: the updrift cell is lo
    up, down = (lo, lo + 1) if up_hi else (lo + 1, lo)
    side = np.where((np.arange(NY) <= lo) == up_hi, 1.0, -1.0)   # +1 on the updrift side
    dist = np.minimum(np.abs(np.arange(NY) - (lo + 0.5)), NY - np.abs(np.arange(NY) - (lo + 0.5))) - 0.5
    eta = -0.5 * notch_m * side * np.exp(-dist / NOTCH_DECAY_CELLS)   # landward +: updrift seaward by notch/2
    eta0 = eta.copy()
    history = history if history is not None else [(1, strength)]
    cb = Recording(updrift_pad=up, downdrift_pad=down, strength_history=history,
                                start_year=1, install_year=1, n_domains=NY, half_width=half_width)
    fake = SimpleNamespace(_brie_coupler=SimpleNamespace(_brie=None))
    bg_cum = np.zeros(NY)
    rows, crashed = [], None
    for year in range(1, years + 1):
        b = brie_like(climate_of_year(year, rng) if climate_of_year else climate)
        b.x_s = tilt + eta
        bg = background(year, rng) if background else np.zeros(NY)
        x_s_dt = list(bg)                         # CASCADE passes Barrier3D's change as a list
        r_now = DIFFUSION_NUMBER(b, b.x_s)
        zone = np.abs(np.arange(NY) - (lo + 0.5)) <= half_width
        if groin:
            if refresh_drift:
                cb._q_net = None
            fake._brie_coupler._brie = b
            try:
                x_s_dt = cb(fake, x_s_dt)
            except IndexError:                    # an angle past BRIE's table: the callback stops the run
                crashed = year
                break
        dx = np.asarray(x_s_dt, dtype=float)
        eta, theta = brie_step(b, slope, eta, dx)
        bg_cum += bg
        dev = eta - eta0 - bg_cum                  # change the background does not explain
        rows.append(dict(year=year, gap_m=float(eta[down] - eta[up]),    # updrift seaward of downdrift, m
                         updrift_seaward_m=float(-(eta - eta0)[up]),
                         downdrift_seaward_m=float(-(eta - eta0)[down]),
                         net_m=float(-dev.sum()), gross_m=float(np.abs(dev).sum()),
                         max_angle_deg=float(np.abs(theta).max()),
                         finite=bool(np.isfinite(eta).all()),
                         strength=cb.blocking_applied_TS[-1] if groin else 0.0,
                         held_m=-cb.applied_dx_updrift_TS[-1] if groin else 0.0,
                         correction_m=cb.leak_correction_m_TS[-1] if groin else 0.0,
                         drift_m3_yr=cb.net_drift_m3_yr_TS[-1] if groin else np.nan,
                         r_face=float(r_now[lo]), r_zone_max=float(r_now[zone].max())))
        if groin and cb.trace:
            rows[-1].update(footprint_net_m=cb.trace[-1][0], footprint_gross_m=cb.trace[-1][1],
                            groin_face_deg=cb.trace[-1][2])
        if not np.isfinite(eta).all():
            break
    frame = pd.DataFrame(rows)
    frame.attrs["eta"] = eta - eta0 - bg_cum
    frame.attrs["crashed"] = crashed
    frame.attrs["guard"] = (cb.guard_year, cb.guard_reason) if groin else (None, None)
    return frame


# The groin run against the same run with the groin off, scored on the rules in CONFIG
def score(g, n, steady=True):
    gross = float(np.abs(g.attrs["eta"] - n.attrs["eta"]).sum())   # the groin's own footprint, summed
    n = n.set_index("year").reindex(g.year).reset_index()
    net = g.net_m.iloc[-1] - n.net_m.iloc[-1]
    gap = g.gap_m.iloc[-1] - n.gap_m.iloc[-1]
    tail = g.gap_m - n.gap_m
    k = max(1, int(round(0.8 * len(tail))))
    drift_tail = float(tail.iloc[-1] - tail.iloc[k - 1])
    out = dict(years_run=len(g), crashed_year=g.attrs["crashed"],
               guard_year=g.attrs["guard"][0], guard_reason=g.attrs["guard"][1],
               finite=bool(g.attrs["crashed"] is None and g.finite.all() and np.isfinite(g.gap_m).all()),
               gap_vs_nogroin_m=gap, net_vs_nogroin_m=net, gross_m=gross,
               max_angle_deg=float(g.max_angle_deg.max()), gap_change_last20pct_m=drift_tail,
               mean_held_m_yr=float(g.held_m.mean()),
               max_gap_diff_m=float((g.gap_m - n.gap_m).abs().max()),
               r_face_min=float(g.r_face.min()), r_face_max=float(g.r_face.max()),
               r_zone_max=float(g.r_zone_max.max()))
    out["conserves"] = bool(abs(net) <= max(CONSERVE_ABS_M, CONSERVE_FRAC * gross))
    out["angle_ok"] = bool(out["max_angle_deg"] < MAX_ANGLE_DEG)
    out["steady"] = bool(abs(drift_tail) <= STEADY_FRAC * abs(gap) + STEADY_ABS_M) if steady else None
    out["passes"] = bool(out["finite"] and out["conserves"] and out["angle_ok"] and (out["steady"] is not False))
    return out


# BRIE's diffusion number r = D dt / (2 dy^2) at the start angle, for the case's (first) climate
def diffusion_number(kw):
    b = brie_like(kw.get("climate", OPTION_A))
    k = int(np.clip(round(90 - kw["theta0"]), 1, len(b._coast_diff) - 1))
    return float(max(0.0, b._coast_diff[k]) * DT_YR / 2 / DY_M ** 2)


# One case of a family: the groin run, its no-groin pair, the score and both year tables
def case(family, label, steady=True, **kw):
    g = run(**kw)
    n = run(groin=False, **{k: v for k, v in kw.items() if k not in ("strength", "history", "refresh_drift")})
    params = {k: (v if np.isscalar(v) else str(v)) for k, v in kw.items()
              if k not in ("climate_of_year", "background")}
    row = dict(family=family, case=label, **params, r_groin=diffusion_number(kw), **score(g, n, steady))
    g.attrs, n.attrs = {}, {}
    by_year = pd.concat([g.assign(groin=True), n.assign(groin=False)]).assign(family=family, case=label,
                                                                                theta0=kw["theta0"])
    return row, by_year


# Family 0: the callback on a straight coast at f 0.6 and a 5-cell zone is version 6 of 1-run/
def rig_check():
    t = sc.brie_tables()
    rows = []
    for th in sc.ORIENTATIONS:
        ref, _ = sc.run(t, th, "corrected", 50)
        g = run(th, years=50, strength=sc.BLOCKING_FRACTION, half_width=sc.CORRECTION_HALF_WIDTH)
        rows.append(dict(theta0=th, updrift_1run_m=ref.updrift_seaward_m.iloc[-1],
                         updrift_callback_m=g.updrift_seaward_m.iloc[-1],
                         reach_1run_m=ref.reach_change_m.iloc[-1], reach_callback_m=g.net_m.iloc[-1]))
    frame = pd.DataFrame(rows)
    frame["max_diff_m"] = np.maximum((frame.updrift_1run_m - frame.updrift_callback_m).abs(),
                                     (frame.reach_1run_m - frame.reach_callback_m).abs())
    return frame


# A wave climate drawn each year from VARIABLE_ASYM x VARIABLE_HS
def variable_climate(year, rng):
    return dict(OPTION_A, asym=float(rng.choice(VARIABLE_ASYM)), Hs=float(rng.choice(VARIABLE_HS)))


BACKGROUNDS = {
    "retreat": lambda year, rng: np.full(NY, RETREAT_M_YR),
    "noise": lambda year, rng: rng.normal(0.0, NOISE_SD_M, NY),
    "pulse": lambda year, rng: np.full(NY, PULSE_M if year % PULSE_EVERY == 0 else 0.0),
}


# Run: every family, every case; write the run table and the year table
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--only", help="run one family by name (waves, strength, ...)")
    a = ap.parse_args()
    TABLES.mkdir(parents=True, exist_ok=True)
    runs, years = [], []

    def add(item):
        runs.append(item[0])
        years.append(item[1])

    want = lambda f: a.only in (None, f)  # noqa: E731
    families(want, add)
    if want("rig"):
        chk = rig_check()
        chk.to_csv(TABLES / "stress_rig_check.csv", index=False, float_format="%.4f")
        print(chk.round(3).to_string(index=False))
    _write(runs, years, a)


# Every family
def families(want, add):
    at = case

    for th in ORIENTATIONS:
        for s in STRENGTHS:
            if want("waves"):
                for p, values in WAVE_SWEEP.items():
                    for v in values:
                        add(at("waves", f"{p}={v}", theta0=th, strength=s, climate=dict(OPTION_A, **{p: v})))
            if want("notch"):
                for nm in NOTCHES_M:
                    add(at("notch", f"notch={nm}", theta0=th, strength=s, notch_m=float(nm)))
            if want("switch-off"):
                add(at("switch-off", "on 50, off 100", steady=False, theta0=th, strength=s, years=150,
                         history=[(1, s), (51, 0.0)]))
            if want("cycling"):
                add(at("cycling", "5 on, 5 off", theta0=th, strength=s,
                         history=[(y, s if (y // 5) % 2 == 0 else 0.0) for y in range(1, YEARS + 1, 5)]))
            if want("background"):
                for name, fn in BACKGROUNDS.items():
                    add(at("background", name, steady=False, theta0=th, strength=s, background=fn, seed=7))
            if want("placement"):
                for lo in PLACEMENTS:
                    add(at("placement", f"lo={lo}", theta0=th, strength=s, lo=lo))
            if want("long"):
                add(at("long", f"{LONG_YEARS} yr", theta0=th, strength=s, years=LONG_YEARS))
        if want("strength"):
            for s in SWEEP_STRENGTHS:
                add(at("strength", f"s={s}", theta0=th, strength=s))
    if want("variable"):
        for th in (-20, 0, 20):
            for s in STRENGTHS:
                for seed in VARIABLE_SEEDS:
                    for refresh in (False, True):
                        add(at("variable", f"seed {seed}, {'refreshed' if refresh else 'as written'}",
                                 steady=False, theta0=th, strength=s, climate_of_year=variable_climate,
                                 refresh_drift=refresh, seed=seed))


# Write the run table and the year table, and print the summary
def _write(runs, years, a):
    if not runs:
        return
    table = pd.DataFrame(runs)
    suffix = f"_{a.only}" if a.only else ""
    table.to_csv(TABLES / f"stress_runs{suffix}.csv", index=False, float_format="%.4f")
    pd.concat(years).to_csv(TABLES / f"stress_by_year{suffix}.csv.gz", index=False, float_format="%.3f")
    pd.set_option("display.width", 220)
    print(table.groupby(["family", "strength"], dropna=False)[["passes", "finite", "conserves", "angle_ok"]]
          .mean().round(2).to_string())
    bad = table[~table.passes]
    print(f"\n{len(bad)} of {len(table)} runs fail")
    cols = ["family", "case", "theta0", "strength", "gap_vs_nogroin_m", "net_vs_nogroin_m", "max_angle_deg",
            "gap_change_last20pct_m", "conserves", "angle_ok", "steady"]
    print(bad[[c for c in cols if c in bad]].round(1).to_string(index=False))


if __name__ == "__main__":
    main()
