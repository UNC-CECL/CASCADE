"""
Which collision-regime dune erosion rule fits the 1996-2009 lidar best, and is the simplest that does.

    python scripts/real-world-cases/Hatteras/dune-collision/3-module-tests/1-rule-comparison/HAT_collision_rule_comparison.py

Offline, no model runs. The cells are HAT_dune_collision.py calibrate's: per survey interval, the 10 m cells the
total water level never overtopped (in collision or swash only), with their earlier crest, toe (per-cell median
over the surveys, 1998 left out) and crest loss to the next survey. Every rule is fit the same way, on
domain medians per interval over the fit intervals, with the same net growth term g (m/yr) the model's dune
growth stands in for: observed loss = median(rule's erosion) - g * years. The rules, each from the model's storm
series (peak Rhigh over the duration, as Barrier3D sees it) unless hourly:

    null              no collision erosion, growth only
    long10_storm      Long, de Bakker & Plant (2014): each storm that reaches the toe takes 10% of D_high
    long_storm        the same, the fraction fit
    long10_year       10% of D_high once in each calendar year a storm reaches the toe
    long_year         the same, the fraction fit
    long_ramp_storm   per storm, a fraction of D_high rising from 0 at the toe to a at the crest,
                      a (Rhigh - toe) / (D_high - toe): Long et al.'s ~10% in collision, up to ~50% near overtopping
    larson_storm      Larson et al. (2004) per storm, 4 Cs (Rhigh - toe)^2 dur / T / 20 m (the module now)
    larson_hourly     Larson hour by hour through the hourly TWL (Palmsten & Holman 2012; Splinter & Palmsten 2012)
    larson_general_storm  Larson et al.'s general solution (Eqs. 19-20) per storm: the dune foot rises up the
                      foreshore (tan b, MHW to the toe, per-cell median) as the face retreats, so the erosion
                      slows and stops when the foot reaches Rhigh; larson_storm is its small-retreat limit (Eq. 25)
    larson_cap_storm  Eq. 25 per storm, but no more than the general solution's limit at the foot = Rhigh
    *_beach           the rule times max(0, 1 - W / W0), W the cell's median beach width over the surveys
                      (1998 left out) and W0 fit; larson_storm_w45 holds W0 at 45 m, the module's setting

Losses compound per storm for the fraction rules, D (1 - prod(1 - f_i)), and no rule takes a crest below its toe.
D_high is the crest in m NAVD88, as in HAT_dune_collision.py empirical. The hourly window runs to the end of the
last storm the interval holds: the 1999 survey's month starts on 1 Sep, two days into Dennis (calibrate's
hourly exposure stops there). Scores: RMSE and r on domain medians, the same for Dennis (1998-1999), the
leave-one-interval-out RMSE (each interval predicted from a fit on the others; the test of a rule, the others
reward parameters), AIC, and the erosion predicted on wide beaches (domain median W > 60 m), where the lidar
shows none. Writes tables/rule_comparison.csv, tables/domain_predictions.csv, figures/rule_comparison.png and
figures/observed_vs_predicted.png under the workspace's dune-collision/3-module-tests/1-rule-comparison/, beside
its NOTE.md (offline, no runs: under output/raw_runs/dune-storm-response/2026-10-29-collision-rule-comparison/
until 2026-11-03)

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-29
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

_HERE = Path(__file__).resolve()
# HAT_dune_collision.py (the toes, storms and calibration cells) sits with the 1996-2009 calibration
_DSR = next(_p for _p in _HERE.parents if (_p / "4-hindcast").is_dir())
sys.path.insert(0, str(_DSR / "4-hindcast" / "1-calibration-1996-2009"))
import HAT_dune_collision as H  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "2026-10-29-collision-rule-comparison"
OUT = H.PROJECT_ROOT / "dune-collision" / "3-module-tests" / "1-rule-comparison"
DENNIS = "1998-1999"
WIDE_BEACH_M = 60.0
LONG_FRACTION = 0.10               # Long et al. (2014): collision crest change within ~10%
W0_MODULE_M = 45.0                 # the module's beach_width_limit on the beach to the crest, until 2026-11-03
W0_DRY_M = 25.0                    # the module's beach_width_limit on the dry beach (to the toe), since 2026-11-03
MODULE_RULE = f"larson_storm_dry{W0_DRY_M:g}"
FORESHORE_MIN_M = 10.0             # shortest MHW-to-toe run for the foreshore slope (one cell)
TANB_RANGE = (0.01, 0.3)           # foreshore slope bounds; Larson et al.'s tanks span 0.06-0.26
# Per-cell runup (2026-11-03): each cell's own peak water level in place of the storm series' island-wide one,
# from the hourly water level, Hs and Tp, and the series' own Stockdon et al. (2006) R2% (H.twl) with
#   stockdon   the cell's beach slope, beach start (0.5 m MHW) to the toe, rise and run from the same point
#   parkcox    Park & Cox (2016, Coastal Eng. 115, 67-78) Eq. 23: R2% x their berm reduction factor C(W_B / L0),
#              mean (Eq. 20) or conservative (Eq. 21) curve, the slope moving from the beach slope to the dune slope
#              as the surge rises over the berm (Eqs. 17-18); W_B the beach start to the toe
# Their Eqs. 6 and 23 print Stockdon's 0.004 as 0.0004; the series' 0.06^2 = 0.0036 is kept, so that the island
# slope reproduces the storm series exactly (checked when the features are built)
BEACH_START_M_NAVD88 = H.MHW_NAVD88_M + 0.5
BERM_M_NAVD88 = 1.7                # Barrier3D's berm
DUNE_SLOPE_RANGE = (0.05, 1.0)
G = 9.81
# -----------------------------------------------------------------------------


# Stockdon et al. (2006) R2% as the storm series computes it (H.twl), for slope b and deep-water Hs, Tp
def stockdon_r2(b, hs, tp):
    hl = hs * G * tp ** 2 / (2 * np.pi)
    return 1.1 * (0.35 * b * np.sqrt(hl) + 0.5 * np.sqrt(hl * (0.5625 * b ** 2 + 0.0036)))


# Park & Cox (2016) berm reduction factor: Eq. 20 (mean) or Eq. 21 (conservative), W_B* = W_B / L0
def parkcox_c(w_b, tp, curve):
    wstar = w_b / (G * tp ** 2 / (2 * np.pi))
    return 0.8 - 0.4 * np.tanh(2 * wstar) if curve == "mean" else 1.0 - 0.5 * np.tanh(2 * wstar)


# Per cell and interval, everything the rules need
def features():
    c = pd.read_csv(H.EXP_DIR / "collision_calibration_cells.csv")
    c = c[c.interval.isin(H.FIT_INTERVALS)].copy()
    st = H.survey_table()
    # Two beach widths, per-cell medians over the surveys (1998 left out): w, beach start to crest, the dune face
    # included (the module's until 2026-11-03), and w_dry, the dry beach, beach start to toe (the module's since):
    # 0 where the crest stands within one cell of the beach start and no toe is picked, as HAT_dune_toes.py
    s0 = st[~st.survey.isin(H.TOE_SURVEYS_EXCLUDED)].copy()
    s0["w_dry"] = np.where(s0.toe_to_crest_m.notna(), s0.beach_width_m - s0.toe_to_crest_m,
                           np.where(s0.beach_width_m < 20, 0.0, np.nan))
    w = s0.groupby(["gis", "cell"]).beach_width_m.median()
    c = c.join(w.rename("w"), on=["gis", "cell"]).dropna(subset=["w"]).reset_index(drop=True)
    c = c.join(s0.groupby(["gis", "cell"]).w_dry.median(), on=["gis", "cell"])
    c["w_dry"] = c.w_dry.fillna(c.groupby("gis").w_dry.transform("median"))
    # Foreshore slope, MHW to the toe, per-cell median over the same surveys
    s = st[~st.survey.isin(H.TOE_SURVEYS_EXCLUDED)]
    run = (s.beach_width_m - s.toe_to_crest_m).clip(lower=FORESHORE_MIN_M)
    tanb = ((s.toe_m_navd88 - H.MHW_NAVD88_M) / run).clip(*TANB_RANGE)
    c = c.join(tanb.groupby([s.gis, s.cell]).median().rename("tanb"), on=["gis", "cell"])
    c["tanb"] = c.tanb.fillna(c.tanb.median())
    # For the per-cell runup rules: the beach slope with rise and run from the beach start, the dry beach width
    # to the toe (Park & Cox's W_B) and the dune face slope (their tan b_D), per-cell medians
    face = ((s.toe_m_navd88 - BEACH_START_M_NAVD88) / run).clip(*TANB_RANGE)
    dune = ((s.crest_m_navd88 - s.toe_m_navd88) / s.toe_to_crest_m.clip(lower=FORESHORE_MIN_M)).clip(*DUNE_SLOPE_RANGE)
    w_toe = (s.beach_width_m - s.toe_to_crest_m).clip(lower=0)
    for name, v in (("tanb_face", face), ("tanb_dune", dune), ("w_toe", w_toe)):
        c = c.join(v.groupby([s.gis, s.cell]).median().rename(name), on=["gis", "cell"])
        c[name] = c[name].fillna(c[name].median())
    flown = pd.read_csv(H.SURVEY_CRESTS).groupby("survey").flown.first()
    storms = H.storm_table()
    hours = pd.read_csv(H.EXP_DIR / "twl_hourly_1996_2009.csv", index_col="time", parse_dates=True)
    c["larson_storm_x"] = c["larson_hourly_x"] = 0.0
    runup_cols = ("x_island", "x_stockdon", "x_parkcox_mean", "x_parkcox_cons")
    for col in runup_cols:
        c[col] = 0.0
    worst_check = 0.0
    c["n_storm"] = c["n_year"] = 0
    ramp = {}                       # interval -> (row index, cells x storms excess / (D - toe))
    general = {}                    # interval -> (row index, excess, crest - Rhigh, t / T per storm)
    for iv, x in c.groupby("interval"):
        a, b = (int(y) for y in iv.split("-"))
        t0, t1 = H.flown_span(flown[a])[1], H.flown_span(flown[b])[0]
        ss = storms[(storms.StartTime >= t0) & (storms.StartTime <= t1)]
        toe, d = x.toe_cell_m_navd88.to_numpy(), x.crest_m_navd88.to_numpy()
        ex = np.clip(ss.rhigh_m_navd88.to_numpy()[None, :] - toe[:, None], 0, None)
        c.loc[x.index, "larson_storm_x"] = 4 * (ex ** 2 * (ss.duration * 3600.0 / ss.period).to_numpy()).sum(1) / H.DUNE_WIDTH_M
        hit = ex > 0
        c.loc[x.index, "n_storm"] = hit.sum(1)
        yrs = ss.StartTime.dt.year.to_numpy()
        c.loc[x.index, "n_year"] = np.column_stack([hit[:, yrs == y].any(1) for y in np.unique(yrs)]).sum(1)
        ramp[iv] = (x.index.to_numpy(), np.clip(ex / (d - toe)[:, None], 0, 1))
        # For Larson's general solution: excess over the toe, crest above Rhigh, number of impacts
        general[iv] = (x.index.to_numpy(), ex, np.clip(d[:, None] - ss.rhigh_m_navd88.to_numpy()[None, :], 0, None),
                       (ss.duration * 3600.0 / ss.period).to_numpy())
        # Each cell's own peak water level per storm, held over the storm as the series' is (Larson per storm)
        bf, bd, wb = (x[k].to_numpy()[:, None] for k in ("tanb_face", "tanb_dune", "w_toe"))
        peaks = {col: [] for col in runup_cols}
        for _, st in ss.iterrows():
            hs_ = hours.loc[st.StartTime:st.EndTime]
            wl, hs, tp = (hs_[k].to_numpy()[None, :] for k in ("water_level", "Hs", "Tp"))
            peaks["x_island"].append(np.max(wl + stockdon_r2(H.BEACH_SLOPE, hs, tp), axis=1).repeat(len(x)))
            peaks["x_stockdon"].append(np.max(wl + stockdon_r2(bf, hs, tp), axis=1))
            alpha = np.clip((wl - BERM_M_NAVD88 + 0.7 * hs) / (1.4 * hs), 0, 1)       # Park & Cox Eq. 18
            b_r = (1 - alpha) * bf + alpha * bd                                        # Eq. 17
            for curve in ("mean", "cons"):
                peaks[f"x_parkcox_{curve}"].append(
                    np.max(wl + parkcox_c(wb, tp, "mean" if curve == "mean" else "cons") * stockdon_r2(b_r, hs, tp),
                           axis=1))
        if len(ss):
            imp = (ss.duration * 3600.0 / ss.period).to_numpy()
            for col in runup_cols:
                r = np.column_stack(peaks[col])                                        # cells x storms
                c.loc[x.index, col] = 4 * (np.clip(r - toe[:, None], 0, None) ** 2 * imp).sum(1) / H.DUNE_WIDTH_M
            worst_check = max(worst_check, float(np.abs(np.column_stack(peaks["x_island"])[0]
                                                        - ss.rhigh_m_navd88.to_numpy()).max()))
        end = max(t1, ss.EndTime.max()) if len(ss) else t1
        hh = hours.loc[t0:end]
        hh = hh[hh.TWL > toe.min()]
        ehr = np.clip(hh.TWL.to_numpy()[None, :] - toe[:, None], 0, None)
        c.loc[x.index, "larson_hourly_x"] = 4 * (ehr ** 2 * 3600.0 / hh.Tp.to_numpy()).sum(1) / H.DUNE_WIDTH_M
        print(f"{iv}: {len(x)} cells, {len(ss)} storms (peak {ss.rhigh_m_navd88.max():.2f} m), hourly to {end:%Y-%m-%d}",
              flush=True)
    # Check: the island slope through the same code reproduces the storm series, peaks and exposure
    print(f"per-cell runup check: island-slope peaks against the series, max |difference| {worst_check:.4f} m; "
          f"exposure max |difference| {np.abs(c.x_island - c.larson_storm_x).max():.4g}", flush=True)
    c["room"] = (c.crest_m_navd88 - c.toe_cell_m_navd88).clip(lower=0)
    ramp["general"] = general
    return c, ramp


# Larson et al. (2004) general solution, one storm at constant Rhigh: the foot rises up the foreshore as the
# face retreats, dz/dt = 4 Cs tanb (R - z)^2 / (T (D - z)), so erosion slows and stops at z = R. With a = R - z
# and h = D - R this integrates to h (1/a - 1/a0) - ln(a / a0) = 4 Cs tanb t / T (Eq. 19, datum free), and the
# volume is ((h + a0)^2 - (h + a)^2) / (2 tanb) (Eq. 20). cap: Eq. 25, no more than that volume at a = 0.
def larson_general_volume(cs, a0, h, tanb, n_impacts, cap_only=False):
    a0 = np.asarray(a0, float)
    v_max = ((h + a0) ** 2 - h ** 2) / (2 * tanb)
    if cap_only:
        return np.minimum(4 * cs * a0 ** 2 * n_impacts, v_max)
    kt = 4 * cs * tanb * n_impacts
    on = a0 > 0
    la0 = np.log(np.where(on, a0, 1.0))
    lo, hi = la0 - 40.0, la0.copy()
    for _ in range(60):             # bisection on ln a; the left side falls monotonically in a
        mid = 0.5 * (lo + hi)
        f = h * (np.exp(-mid) - np.exp(-la0)) - (mid - la0) - kt
        lo, hi = np.where(f > 0, mid, lo), np.where(f > 0, hi, mid)
    a = np.exp(0.5 * (lo + hi))
    return np.where(on, ((h + a0) ** 2 - (h + a) ** 2) / (2 * tanb), 0.0)


# Erosion per cell (m), before the growth term
def erosion(rule, p, c, ramp):
    width = c[rule.get("wcol", "w")].to_numpy()          # w (to the crest) unless the rule names w_dry
    beach = np.clip(1 - width / p["w0"], 0, None) if "w0" in p else 1.0
    d = c.crest_m_navd88.to_numpy()
    kind = rule["kind"]
    if kind == "null":
        e = np.zeros(len(c))
    elif kind == "larson_storm":
        e = p["cs"] * c.larson_storm_x.to_numpy() * beach
    elif kind == "larson_runup":
        e = p["cs"] * c[rule["col"]].to_numpy() * beach
    elif kind == "larson_hourly":
        e = p["cs"] * c.larson_hourly_x.to_numpy() * beach
    elif kind in ("larson_general", "larson_cap"):
        e = np.zeros(len(c))
        b = beach if np.ndim(beach) else np.full(len(c), beach)
        tanb = c.tanb.to_numpy()
        for iv in c.interval.unique():
            idx, a0, h, n_imp = ramp["general"][iv]
            keep = c.index.get_indexer(idx)
            sel = keep >= 0
            pos = keep[sel]
            v = larson_general_volume(p["cs"], a0[sel], h[sel], tanb[pos, None], n_imp[None, :],
                                      cap_only=kind == "larson_cap")
            e[pos] = v.sum(1) * b[pos] / H.DUNE_WIDTH_M
    elif kind in ("long_storm", "long_year"):
        n = c.n_storm.to_numpy() if kind == "long_storm" else c.n_year.to_numpy()
        e = d * (1 - (1 - p["f"] * beach) ** n)
    elif kind == "long_ramp_storm":
        e = np.zeros(len(c))
        b = beach if np.ndim(beach) else np.full(len(c), beach)
        for iv in c.interval.unique():
            idx, frac = ramp[iv]
            pos = c.index.get_indexer(idx)
            e[pos] = d[pos] * (1 - np.prod(1 - p["a"] * frac * b[pos, None], axis=1))
    return np.minimum(e, c.room.to_numpy())


RULES = [
    dict(name="null", kind="null", fit=[], fixed={}),
    dict(name="long10_storm", kind="long_storm", fit=[], fixed={"f": LONG_FRACTION}),
    dict(name="long_storm", kind="long_storm", fit=["f"], fixed={}),
    dict(name="long10_year", kind="long_year", fit=[], fixed={"f": LONG_FRACTION}),
    dict(name="long_year", kind="long_year", fit=["f"], fixed={}),
    dict(name="long_ramp_storm", kind="long_ramp_storm", fit=["a"], fixed={}),
    dict(name="larson_storm", kind="larson_storm", fit=["cs"], fixed={}),
    dict(name="larson_hourly", kind="larson_hourly", fit=["cs"], fixed={}),
    dict(name="long10_storm_beach", kind="long_storm", fit=["w0"], fixed={"f": LONG_FRACTION}),
    dict(name="long_storm_beach", kind="long_storm", fit=["f", "w0"], fixed={}),
    dict(name="long_year_beach", kind="long_year", fit=["f", "w0"], fixed={}),
    dict(name="long_ramp_storm_beach", kind="long_ramp_storm", fit=["a", "w0"], fixed={}),
    dict(name="larson_storm_w45", kind="larson_storm", fit=["cs"], fixed={"w0": W0_MODULE_M}),
    dict(name="larson_storm_beach", kind="larson_storm", fit=["cs", "w0"], fixed={}),
    dict(name="larson_hourly_beach", kind="larson_hourly", fit=["cs", "w0"], fixed={}),
    dict(name="larson_general_storm", kind="larson_general", fit=["cs"], fixed={}),
    dict(name="larson_general_storm_w45", kind="larson_general", fit=["cs"], fixed={"w0": W0_MODULE_M}),
    dict(name="larson_general_storm_beach", kind="larson_general", fit=["cs", "w0"], fixed={}),
    dict(name="larson_cap_storm_w45", kind="larson_cap", fit=["cs"], fixed={"w0": W0_MODULE_M}),
    # Per-cell runup in place of the beach factor (2026-11-03); *_w45 asks whether the ramp still adds anything
    dict(name="larson_stockdon", kind="larson_runup", col="x_stockdon", fit=["cs"], fixed={}),
    dict(name="larson_parkcox_mean", kind="larson_runup", col="x_parkcox_mean", fit=["cs"], fixed={}),
    dict(name="larson_parkcox_cons", kind="larson_runup", col="x_parkcox_cons", fit=["cs"], fixed={}),
    dict(name="larson_stockdon_w45", kind="larson_runup", col="x_stockdon", fit=["cs"], fixed={"w0": W0_MODULE_M}),
    dict(name="larson_parkcox_mean_w45", kind="larson_runup", col="x_parkcox_mean", fit=["cs"],
         fixed={"w0": W0_MODULE_M}),
    # The beach factor on the dry beach, beach start to toe (2026-11-03; the rules above use the beach to the crest)
    dict(name="larson_storm_dry", kind="larson_storm", wcol="w_dry", fit=["cs", "w0"], fixed={}),
    dict(name=MODULE_RULE, kind="larson_storm", wcol="w_dry", fit=["cs"], fixed={"w0": W0_DRY_M}),
    dict(name="long10_storm_dry", kind="long_storm", wcol="w_dry", fit=["w0"], fixed={"f": LONG_FRACTION}),
    dict(name="larson_general_storm_dry", kind="larson_general", wcol="w_dry", fit=["cs", "w0"], fixed={}),
    dict(name="larson_hourly_dry", kind="larson_hourly", wcol="w_dry", fit=["cs", "w0"], fixed={}),
]
START = {"cs": 5e-5, "f": 0.05, "a": 0.2, "w0": 45.0}
BOUNDS = {"cs": (0, 1e-2), "f": (0, 1), "a": (0, 1), "w0": (5, 500)}
SCALE = {"cs": 1e-4, "f": 0.1, "a": 0.1, "w0": 10.0}


# Domain medians per interval of a per-cell array
def dom_median(c, v):
    return pd.Series(v, index=c.index).groupby([c.interval, c.gis]).median()


def fit(rule, c, ramp, keep):
    cc = c[keep]
    obs = dom_median(cc, cc.crest_loss_m.to_numpy())
    years = cc.groupby(["interval", "gis"]).years.first()
    names = rule["fit"] + ["g"]

    def params(v):
        return {**rule["fixed"], **dict(zip(names, v))}

    def resid(v):
        p = params(v)
        return (dom_median(cc, erosion(rule, p, cc, ramp)) - p["g"] * years - obs).to_numpy()

    x0 = [START[n] for n in rule["fit"]] + [0.0]
    lo = [BOUNDS[n][0] for n in rule["fit"]] + [-1.0]
    hi = [BOUNDS[n][1] for n in rule["fit"]] + [1.0]
    sc = [SCALE[n] for n in rule["fit"]] + [0.01]
    best = None
    # A second start from a stronger rule, so a fraction or Cs stuck at 0 is a fit, not a start
    for k in (1.0, 4.0):
        r = least_squares(resid, [v * k if i < len(rule["fit"]) and names[i] != "w0" else v for i, v in enumerate(x0)],
                          bounds=(lo, hi), x_scale=sc)
        if best is None or r.cost < best.cost:
            best = r
    return params(best.x)


def predict(rule, p, c, ramp, keep):
    cc = c[keep]
    years = cc.groupby(["interval", "gis"]).years.first()
    e = dom_median(cc, erosion(rule, p, cc, ramp))
    out = pd.DataFrame({"obs": dom_median(cc, cc.crest_loss_m.to_numpy()), "erosion": e, "pred": e - p["g"] * years,
                        "w": cc.groupby(["interval", "gis"]).w.median(),
                        "w_dry": cc.groupby(["interval", "gis"]).w_dry.median()})
    return out


def score(rule, c, ramp):
    p = fit(rule, c, ramp, np.ones(len(c), bool))
    d = predict(rule, p, c, ramp, np.ones(len(c), bool))
    cv = []
    for iv in H.FIT_INTERVALS:
        hold = (c.interval == iv).to_numpy()
        q = fit(rule, c, ramp, ~hold)
        cv.append(predict(rule, q, c, ramp, hold))
    cv = pd.concat(cv)
    res = d.pred - d.obs
    n, k = len(d), len(rule["fit"]) + 1
    den = d.xs(DENNIS, level="interval")
    wide = d[d.w > WIDE_BEACH_M]
    row = dict(rule=rule["name"], params=k,
               cs=p.get("cs"), fraction=p.get("f", p.get("a")), w0_m=p.get("w0"), growth_m_yr=p["g"],
               rmse=np.sqrt((res ** 2).mean()), r=np.corrcoef(d.pred, d.obs)[0, 1] if d.pred.std() > 0 else np.nan,
               cv_rmse=np.sqrt(((cv.pred - cv.obs) ** 2).mean()),
               dennis_rmse=np.sqrt(((den.pred - den.obs) ** 2).mean()),
               dennis_r=np.corrcoef(den.pred, den.obs)[0, 1] if den.pred.std() > 0 else np.nan,
               dennis_cv_rmse=np.sqrt(((lambda x: x.pred - x.obs)(cv.xs(DENNIS, level="interval")) ** 2).mean()),
               aic=n * np.log((res ** 2).mean()) + 2 * k,
               wide_beach_erosion_m=wide.erosion.mean(), wide_beach_obs_m=wide.obs.mean(), n=n)
    print(f"{rule['name']:24s} cv {row['cv_rmse']:.3f}  r {row['r']:.2f}  Dennis r {row['dennis_r']:.2f}  "
          f"{ {k: float(f'{v:.3g}') for k, v in p.items()} }", flush=True)
    return row, d.assign(rule=rule["name"]), cv.assign(rule=rule["name"])


def plots(tab, pred):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import apply_style, record_caption
    apply_style()
    fam = {"null": ("#7A7A7A", "No collision"), "long": ("#D95F02", "Long et al. (2014) fraction"),
           "larson": ("#1B7837", "Larson et al. (2004)")}

    def family(name):
        return "null" if name == "null" else name.split("_")[0].replace("long10", "long")

    t = tab.sort_values("cv_rmse", ascending=False).reset_index(drop=True)
    fig, axs = plt.subplots(1, 2, figsize=(10, 5.2), sharey=True)
    y = np.arange(len(t))
    for ax, col, label in ((axs[0], "cv_rmse", "Leave-one-interval-out RMSE (m)"),
                           (axs[1], "dennis_r", "Dennis (1998–99) r, domain medians")):
        for i, row in t.iterrows():
            colr = fam[family(row.rule)][0]
            beach = "beach" in row.rule or "w45" in row.rule
            ax.plot(row[col], i, "o", ms=8, color=colr, mfc=colr if beach else "white", mew=2)
        ax.set_xlabel(label)
        ax.grid(axis="x", color="0.9")
    axs[0].axvline(t.loc[t.rule == "null", "cv_rmse"].iloc[0], color="0.6", lw=1, ls="--")
    axs[0].set_yticks(y, [f"{r}  ({k})" for r, k in zip(t.rule, t.params)], fontsize=8)
    axs[0].set_title("Lower is better", loc="left", fontsize=9)
    axs[1].set_title("Higher is better", loc="left", fontsize=9)
    hs = [plt.Line2D([], [], marker="o", ls="", color=c, mfc=c, mew=2, label=l) for c, l in fam.values()]
    hs.append(plt.Line2D([], [], marker="o", ls="", color="0.3", mfc="white", mew=2, label="no beach factor"))
    hs.append(plt.Line2D([], [], marker="o", ls="", color="0.3", mfc="0.3", mew=2, label="× beach factor"))
    axs[1].legend(handles=hs, fontsize=7, loc="lower right")
    fig.tight_layout()
    out = OUT / "figures" / "rule_comparison.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, facecolor="white")
    record_caption(out, (
        "Collision-regime dune erosion rules against the 1996–2009 lidar crest change, domain medians per survey "
        "interval of the cells the water never overtopped. Left: RMSE when each interval is predicted from a fit "
        "on the others (dashed: no collision erosion, growth only); right: correlation on the Dennis interval. "
        "In brackets, the fitted parameters, the growth term included. Filled: the rule times "
        f"max(0, 1 − W / W0) for beach width W. Experiment {TAG}."))
    plt.close(fig)

    show = ["null", "long10_storm", "long_year_beach", "larson_storm", "larson_hourly",
            "larson_storm_w45", "larson_hourly_beach", "long_ramp_storm_beach"]
    show = [s for s in show if s in set(pred.rule)]
    fig, axs = plt.subplots(2, 4, figsize=(12, 6.4), sharex=True, sharey=True)
    lim = (-1.5, 3.0)
    for ax, name in zip(axs.flat, show):
        d = pred[pred.rule == name]
        den = d.index.get_level_values("interval") == DENNIS
        ax.plot(d.pred[~den], d.obs[~den], "o", ms=3, color="0.65", label="other intervals")
        ax.plot(d.pred[den], d.obs[den], "o", ms=4, color="#B2182B", label="Dennis, 1998–99")
        ax.plot(lim, lim, color="k", lw=0.8)
        row = tab.set_index("rule").loc[name]
        ax.set_title(f"{name}\nCV RMSE {row.cv_rmse:.2f} m, r {row.r:.2f}", loc="left", fontsize=8.5)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
    for ax in axs.flat[len(show):]:
        ax.axis("off")
    axs[0, 0].legend(fontsize=7)
    for ax in axs[-1]:
        ax.set_xlabel("Predicted crest lowering (m)")
    for ax in axs[:, 0]:
        ax.set_ylabel("Lidar crest lowering (m)")
    fig.tight_layout()
    out = OUT / "figures" / "observed_vs_predicted.png"
    fig.savefig(out, dpi=200, facecolor="white")
    record_caption(out, (
        "Lidar against predicted crest lowering, domain medians per survey interval, for eight of the rules, each "
        f"fit on all the intervals (Dennis red); line 1:1. Experiment {TAG}."))
    plt.close(fig)
    print(f"wrote {OUT / 'figures'}")


# Why Larson alone is not enough (2026-11-03): Larson per storm with and without the beach factor, predicted against
# the lidar on domain medians, and for Dennis by beach width. Reads the tables main() writes
def larson_alone_figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import INK, INK_MUTED, apply_style, record_caption
    apply_style()
    tab = pd.read_csv(OUT / "tables" / "rule_comparison.csv").set_index("rule")
    pred = pd.read_csv(OUT / "tables" / "domain_predictions.csv")
    alone, ramp_ = "larson_storm", MODULE_RULE
    obs_c, alone_c, ramp_c, dennis_c = INK, "#2166AC", "#D94801", "#B2182B"
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.9), gridspec_kw={"width_ratios": [1, 1, 1.25]})
    lim = (-1.5, 3.0)
    for ax, name, title in ((axs[0], alone, "Larson alone"), (axs[1], ramp_, "Larson × beach factor")):
        d = pred[pred.rule == name]
        den = d.interval == DENNIS
        ax.plot(d.pred[~den], d.obs[~den], "o", ms=3, color="0.72", label="other intervals")
        ax.plot(d.pred[den], d.obs[den], "o", ms=4.5, color=dennis_c, label="Dennis, 1998–99")
        ax.plot(lim, lim, color=INK, lw=0.8)
        row = tab.loc[name]
        ax.set_title(f"{title}\nCV RMSE {row.cv_rmse:.2f} m · r {row.r:.2f} · Dennis r {row.dennis_r:.2f} · "
                     f"Cs {row.cs:.1e}", loc="left", fontsize=8.5)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xlabel("Predicted crest lowering (m)")
    axs[0].set_ylabel("Lidar crest lowering (m)")
    axs[0].legend(fontsize=7, loc="upper left")
    ax = axs[2]
    bins = np.arange(0, 101, 10)
    mids = (bins[:-1] + bins[1:]) / 2
    for name, col, label, mk in ((None, obs_c, "Lidar", "o"), (alone, alone_c, "Larson alone", "s"),
                                 (ramp_, ramp_c, "Larson × beach factor", "^")):
        d = pred[(pred.rule == (name or alone)) & (pred.interval == DENNIS)]
        g = d.groupby(pd.cut(d.w_dry, bins), observed=False)
        n = g.size().to_numpy()
        v = (g.obs if name is None else g.pred).median().to_numpy()
        ok = n >= 3
        ax.plot(mids[ok], v[ok], mk + "-", color=col, lw=2.0 if name is None else 1.6, ms=5, label=label)
    ax.axhline(0, color=INK_MUTED, lw=0.7)
    ax.axvline(W0_DRY_M, color=INK_MUTED, lw=0.7, ls=":")
    ax.set_xlim(0, 100)
    ax.set_xlabel("Dry beach width, beach start to toe, domain median (m)")
    ax.set_ylabel("Crest lowering, Dennis (m)")
    ax.set_title("Dennis, by beach width (domains with ≥ 3)", loc="left", fontsize=8.5)
    ax.legend(fontsize=7, loc="upper right")
    for i, ax in enumerate(axs):
        ax.text(0.0, 1.13, f"({'abc'[i]})", transform=ax.transAxes, fontweight="bold", fontsize=10, va="bottom")
    fig.tight_layout()
    out = OUT / "figures" / "larson_alone.png"
    fig.savefig(out, dpi=250, facecolor="white")
    plt.close(fig)
    record_caption(out, (
        f"Larson et al. (2004) per storm, with and without the beach factor max(0, 1 − W/{W0_DRY_M:g} m) on the dry beach (beach start to toe), against the "
        "1996–2009 lidar crest change on the cells the water never overtopped; domain medians per survey interval, "
        "each rule fit on all the intervals (Cs and a growth term). (a, b) Lidar against predicted crest lowering, "
        "Dennis red; line 1:1. (c) Dennis only, by each domain's median beach width (10 m bins with at least three "
        f"domains): the lidar and the two predictions; dotted, {W0_DRY_M:g} m. Experiment {TAG}."))
    print(f"wrote {out}")


def main():
    c, ramp = features()
    rows, preds = [], []
    for rule in RULES:
        row, d, _ = score(rule, c, ramp)
        rows.append(row)
        preds.append(d)
    tab = pd.DataFrame(rows)
    pred = pd.concat(preds)
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    tab.round(6).to_csv(OUT / "tables" / "rule_comparison.csv", index=False)
    pred.round(4).to_csv(OUT / "tables" / "domain_predictions.csv")
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(tab.sort_values("cv_rmse").round(4).to_string(index=False))
    plots(tab, pred)
    larson_alone_figure()


if __name__ == "__main__":
    main()
