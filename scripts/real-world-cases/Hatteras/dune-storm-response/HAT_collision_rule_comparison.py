"""
Which collision-regime dune erosion rule fits the 1996-2009 lidar best, and is the simplest that does.

    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_collision_rule_comparison.py

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
    *_beach           the rule times max(0, 1 - W / W0), W the cell's median beach width over the surveys
                      (1998 left out) and W0 fit; larson_storm_w45 holds W0 at 45 m, the module's setting

Losses compound per storm for the fraction rules, D (1 - prod(1 - f_i)), and no rule takes a crest below its toe.
D_high is the crest in m NAVD88, as in HAT_dune_collision.py empirical. The hourly window runs to the end of the
last storm the interval holds: the 1999 survey's month starts on 1 Sep, two days into Dennis (calibrate's
hourly exposure stops there). Scores: RMSE and r on domain medians, the same for Dennis (1998-1999), the
leave-one-interval-out RMSE (each interval predicted from a fit on the others; the test of a rule, the others
reward parameters), AIC, and the erosion predicted on wide beaches (domain median W > 60 m), where the lidar
shows none. Writes tables/rule_comparison.csv, tables/domain_predictions.csv, figures/rule_comparison.png and
figures/observed_vs_predicted.png under dune-storm-response/2026-10-29-collision-rule-comparison/

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
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_collision as H  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "2026-10-29-collision-rule-comparison"
OUT = H.PROJECT_ROOT / "output" / "raw_runs" / "dune-storm-response" / TAG
DENNIS = "1998-1999"
WIDE_BEACH_M = 60.0
LONG_FRACTION = 0.10               # Long et al. (2014): collision crest change within ~10%
W0_MODULE_M = 45.0                 # the module's beach_width_limit
# -----------------------------------------------------------------------------


# Per cell and interval, everything the rules need
def features():
    c = pd.read_csv(H.EXP_DIR / "collision_calibration_cells.csv")
    c = c[c.interval.isin(H.FIT_INTERVALS)].copy()
    st = H.survey_table()
    w = st[~st.survey.isin(H.TOE_SURVEYS_EXCLUDED)].groupby(["gis", "cell"]).beach_width_m.median()
    c = c.join(w.rename("w"), on=["gis", "cell"]).dropna(subset=["w"]).reset_index(drop=True)
    flown = pd.read_csv(H.SURVEY_CRESTS).groupby("survey").flown.first()
    storms = H.storm_table()
    hours = pd.read_csv(H.EXP_DIR / "twl_hourly_1996_2009.csv", index_col="time", parse_dates=True)
    c["larson_storm_x"] = c["larson_hourly_x"] = 0.0
    c["n_storm"] = c["n_year"] = 0
    ramp = {}                       # interval -> (row index, cells x storms excess / (D - toe))
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
        end = max(t1, ss.EndTime.max()) if len(ss) else t1
        hh = hours.loc[t0:end]
        hh = hh[hh.TWL > toe.min()]
        ehr = np.clip(hh.TWL.to_numpy()[None, :] - toe[:, None], 0, None)
        c.loc[x.index, "larson_hourly_x"] = 4 * (ehr ** 2 * 3600.0 / hh.Tp.to_numpy()).sum(1) / H.DUNE_WIDTH_M
        print(f"{iv}: {len(x)} cells, {len(ss)} storms (peak {ss.rhigh_m_navd88.max():.2f} m), hourly to {end:%Y-%m-%d}",
              flush=True)
    c["room"] = (c.crest_m_navd88 - c.toe_cell_m_navd88).clip(lower=0)
    return c, ramp


# Erosion per cell (m), before the growth term
def erosion(rule, p, c, ramp):
    beach = np.clip(1 - c.w.to_numpy() / p["w0"], 0, None) if "w0" in p else 1.0
    d = c.crest_m_navd88.to_numpy()
    kind = rule["kind"]
    if kind == "null":
        e = np.zeros(len(c))
    elif kind == "larson_storm":
        e = p["cs"] * c.larson_storm_x.to_numpy() * beach
    elif kind == "larson_hourly":
        e = p["cs"] * c.larson_hourly_x.to_numpy() * beach
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
                        "w": cc.groupby(["interval", "gis"]).w.median()})
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


if __name__ == "__main__":
    main()
