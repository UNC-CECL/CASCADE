"""
Which storm duration rule should the hindcast use?

    python scripts/hatteras_ms/experiments/HAT_storm_length_selection.py build
    python scripts/hatteras_ms/experiments/HAT_storm_length_selection.py run --workers 4
    python scripts/hatteras_ms/experiments/HAT_storm_length_selection.py validate
    python scripts/hatteras_ms/experiments/HAT_storm_length_selection.py score
    python scripts/hatteras_ms/experiments/HAT_storm_length_selection.py ends --variants drop72 trim24
    python scripts/hatteras_ms/experiments/HAT_storm_length_selection.py figures

Every event kept and trimmed to L hours (or dropped past 72 h, the control), run in
both windows and scored on overwash first. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-29
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(_HERE.parent))
import HAT_storm_max_duration as MD  # noqa: E402  (builder functions, launcher)

# --- CONFIG ------------------------------------------------------------------
TAG = "storms-and-overwash/2026-09-28-storm-length-selection"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
MATRIX = PROJECT_ROOT / "output" / "raw_runs" / "matrix"
CONTROL_MATRIX = PROJECT_ROOT / "output" / "raw_runs" / "archive" / "2026-09-28-loess10-ends" / "matrix"   # where the drop72 controls scored here now live
WINDOWS = MD.WINDOWS
TRIMS = (24, 36, 48, 72, 96, 120, 168)
VARIANTS = ["drop72"] + [f"trim{L}" for L in TRIMS] + ["full"]
RUN_VARIANTS = [v for v in VARIANTS if v != "drop72"]
SCENARIOS = ("full_management", "natural")
CONTROLS = MD.CONTROLS
THRESHOLD = 0.0          # m3/m of shared overwash in the window
THRESHOLDS = (0.0, 1.0, 5.0)
# -----------------------------------------------------------------------------


# A window's folder name, start_end
def wtag(w):
    return f"{w[0]}_{w[1]}"


# A window's storm folder
def storm_dir(w):
    return EXP_DIR / "storms" / wtag(w)


# A variant's storm file
def storm_file(w, v):
    return storm_dir(w) / f"{wtag(w)}_storms_{v}.npy"


# A variant's storm summary
def summary_csv(w, v):
    return storm_dir(w) / f"{wtag(w)}_storms_{v}_summary.csv"


# Build

# Build the trimmed series from the builder's own functions
def build():
    fn = MD.builder_functions()
    from site_layer import hat_env_forcings as env
    for w in WINDOWS:
        storm_dir(w).mkdir(parents=True, exist_ok=True)
        df = MD.merged_record(fn, w)
        for v, mx in (("drop72", 72), ("full", 10 ** 6)):
            fn["create_storms"](df_merged=df, berm_elevation=MD.BERM, weather_grouping=MD.GROUPING, MHW=MD.MHW,
                                min_storm_dur=MD.MIN_DUR, max_storm_dur=mx, save_dfs=True,
                                save_dir=str(storm_dir(w)), save_name=f"{wtag(w)}_storms_{v}",
                                window_start_year=w[0])
        full = pd.read_csv(summary_csv(w, "full"), parse_dates=["StartTime", "EndTime"])
        for L in TRIMS:
            s = MD.trim_long_events(df, full, w, limit=L)
            s.to_csv(summary_csv(w, f"trim{L}"), index=False)
            np.save(storm_file(w, f"trim{L}"), s[["time", "Rhigh", "Rlow", "period", "duration"]].to_numpy())
        same = np.array_equal(np.load(env.storm_series_file(*w)), np.load(storm_file(w, "drop72")))
        print(f"{wtag(w)}: drop72 rebuild == committed series: {same}")
        if not same:
            raise SystemExit("rebuilt control differs from the committed series; stopping")
        for v in VARIANTS:
            a = np.load(storm_file(w, v))
            print(f"  {v:8s} {len(a):4d} storms  max {a[:, 4].max():4.0f} h  storm-hours {a[:, 4].sum():6.0f}  "
                  f"max Rhigh {a[:, 1].max() * 10:.2f} m")


# Runs

# Where a variant's run lands: matrix ends at stage 1, ends_<variant>/ at stage 2
def run_path(variant, scenario, w, ends=None):
    group = f"{variant}_{scenario}" if ends is None else f"ends_{variant}/{ends}_{scenario}"
    base = EXP_DIR / "runs" / group / wtag(w) / "edgeBE"
    hits = sorted(base.glob("HAT_*")) if base.exists() else []
    return hits[-1] if hits else None


# Run one job in its own process, logged, via this file's _launch entry
def launch(w, variant, scenario, group, override=None, logname=None):
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    log = logs / f"{logname or group.replace('/', '__')}_{wtag(w)}.log"
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(w[0]),
                "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": scenario,
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{group}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8"})
    if override:
        env["HAT_BE_OVERRIDE"] = override
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(MD._HERE), "_launch", str(w[0]), str(storm_file(w, variant))],
                           env=env, cwd=MD.HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    rec = dict(window=wtag(w), variant=variant, scenario=scenario, group=group, override=override,
               returncode=p.returncode, minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{wtag(w)} {group:40s} {override or '':22s} exit {p.returncode} {rec['minutes']} min", flush=True)
    return rec


# Every variant x scenario x window, in parallel
def run(workers):
    jobs = [(w, v, s) for w in WINDOWS for v in RUN_VARIANTS for s in SCENARIOS]
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(lambda j: launch(j[0], j[1], j[2], f"{j[1]}_{j[2]}"), jobs))


# Which storm made the overwash

# Barrier3D's gap discharge per cell, dam^3/hr (barrier3d.py, gap loop)
def _qdune(rexcess_dam):
    return math.sqrt(2 * 9.8 * (rexcess_dam * 10)) / 10 * rexcess_dam * 3600


# Each domain-year's overwash shared among that year's storms that reached a gap
def storm_shares(c, summ):
    out = {}
    for p, b in enumerate(c.barrier3d):
        q = np.asarray(b.QowTS)
        for t in range(1, len(q)):
            if q[t] <= 0:
                continue
            dd = np.array(b.DuneDomain, dtype=float, copy=True)
            dd, _, _ = type(b).DuneGrowth(b, dd, t)
            crest = dd[t].max(axis=1)
            crest[crest < b._DuneRestart] = b._DuneRestart
            rows = summ.index[summ.time == t]
            weights = {}
            for i in rows:
                rh = summ.at[i, "Rhigh"]
                dow = [k for k, v in enumerate(crest + b._BermEl) if v < rh]
                gaps = b.DuneGaps(crest, dow, b._BermEl, rh)
                wgt = sum((g[1] - g[0] + 1) * _qdune(g[2]) for g in gaps if g[2] > 0) * summ.at[i, "duration"]
                if wgt > 0:
                    weights[i] = wgt
            tot = sum(weights.values())
            for i, wgt in weights.items():
                out[(p, i)] = q[t] * wgt / tot
    return out


# A run's saved CASCADE object
def load_state(d):
    return np.load(d / f"{d.name}.npz", allow_pickle=True)["cascade"][0]


# overwash_vs_model.py, loaded by path
def overwash_module():
    path = PROJECT_ROOT / "scripts" / "input_prep" / "8-overwash-analysis" / "4-vs-model" / "overwash_vs_model.py"
    spec = importlib.util.spec_from_file_location("overwash_vs_model", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# Observed and modelled overwash, one row per assessed image and domain
def overwash_cells(c, w, summ, obs, ovm):
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    shares = storm_shares(c, summ)
    dated = pd.to_datetime(summ.EndTime) - ovm.GRACE
    rows = []
    for _, im in ovm.images_in(w, obs).iterrows():
        o = obs[obs.Obs_ID == im.Obs_ID].set_index("domain").overwash
        in_win = set(summ.index[(dated > im["from"]) & (dated <= im.date)])
        for gis in range(1, 91):
            p = DOM.gis_to_pad(gis)
            vol = sum(v for (pp, i), v in shares.items() if pp == p and i in in_win)
            rows.append(dict(obs_id=im.Obs_ID, image=im.date, gis=gis,
                             observed=o.get(gis, np.nan), model_m3_per_m=vol))
    return pd.DataFrame(rows)


# Hit, miss and false-alarm scores for one overwash threshold
def overwash_scores(cells, thr=THRESHOLD):
    a = cells.dropna(subset=["observed"]).copy()
    o, m = a.observed.astype(int), (a.model_m3_per_m > thr).astype(int)
    hits, miss = int(((o == 1) & (m == 1)).sum()), int(((o == 1) & (m == 0)).sum())
    fa, cn = int(((o == 0) & (m == 1)).sum()), int(((o == 0) & (m == 0)).sum())
    pod, pofd = hits / max(hits + miss, 1), fa / max(fa + cn, 1)
    a["m"] = m
    per_img = a.groupby("obs_id").agg(o=("observed", "sum"), m=("m", "sum"))
    per_dom = a.groupby("gis").agg(o=("observed", "mean"), m=("m", "mean"))
    r = lambda x, y: float(np.corrcoef(x, y)[0, 1]) if x.std() > 0 and y.std() > 0 else np.nan  # noqa: E731
    return dict(cells=len(a), hits=hits, misses=miss, false_alarms=fa, correct_negatives=cn,
                POD=pod, POFD=pofd, PSS=pod - pofd, CSI=hits / max(hits + miss + fa, 1),
                timing_r=r(per_img.o, per_img.m), space_r=r(per_dom.o, per_dom.m),
                model_cells=int(m.sum()), observed_cells=int(o.sum()))


_TARGETS = {}                        # (start year, LOWESS domains) -> CoastSat target, built once


# The CoastSat LRR target for a start year at a LOWESS window, cached
def _target(start, lowess_domains):
    import HAT_metres_1_offset_units as O
    if (start, lowess_domains) not in _TARGETS:
        _TARGETS[(start, lowess_domains)] = O.coastsat_target(start, lowess_domains=lowess_domains)
    return _TARGETS[(start, lowess_domains)]


# Interior bias and RMSE (LRR and endpoint) at today's window, after the run's stored ones are reproduced at its own
def rescored_skill(d):
    import re
    import HAT_metres_1_offset_units as O
    meta = json.loads((d / f"{d.name}_run_metadata.json").read_text(encoding="utf-8"))
    sk = meta["skill"]
    start = int(re.match(r"HAT_(\d{4})_", d.name).group(1))
    window = int(re.search(r"(\d+)-domain", sk["target"]).group(1))
    for k, v in O.rate_skill(d, _target(start, window)).items():
        if not np.isclose(v, float(sk[k]), rtol=1e-3, atol=5e-5):      # metadata stores 4 decimals
            raise ValueError(f"{d}: {k} against the rebuilt {window}-domain target does not "
                             f"match the runner's -- not the same target")
    return O.rate_skill(d, _target(start, O.SMOOTH_DOMAINS))


# A run's interior skill against CoastSat, at today's LOWESS window
def shoreline_scores(d):
    sk = rescored_skill(d)
    return dict(rmse_interior_m_yr=sk["rmse_interior_m_yr"],
                bias_interior_m_yr=sk["mean_bias_interior_m_yr"])


# Score one stage's runs: overwash against the imagery, then shoreline
def score(stage="1"):
    from site_layer import hat_overwash as ow
    ovm = overwash_module()
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    rows, cells_all = [], []
    ends = load_ends() if stage == "2" else {}
    for w in WINDOWS:
        for v in VARIANTS:
            if stage == "2" and v not in ends.get(wtag(w), {}):
                continue
            summ = pd.read_csv(summary_csv(w, v))
            for s in SCENARIOS:
                if stage == "1":
                    d = (CONTROL_MATRIX / wtag(w) / "edgeBE" / CONTROLS[(w[0], s)]) if v == "drop72" else run_path(v, s, w)
                else:
                    d = run_path(v, s, w, ends="solved")
                if d is None:
                    print(f"  missing run: {wtag(w)} {v} {s}")
                    continue
                cells = overwash_cells(load_state(d), w, summ, obs, ovm)
                cells_all.append(cells.assign(window=wtag(w), variant=v, scenario=s))
                for thr in THRESHOLDS:
                    rows.append(dict(window=wtag(w), variant=v, scenario=s, threshold=thr,
                                     **overwash_scores(cells, thr), **shoreline_scores(d)))
                print(f"  scored {wtag(w)} {v} {s}", flush=True)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t = pd.DataFrame(rows)
    t.to_csv(EXP_DIR / "tables" / f"stage{stage}_scores.csv", index=False)
    pd.concat(cells_all).to_csv(EXP_DIR / "tables" / f"stage{stage}_cells.csv", index=False)
    h = t[t.threshold == THRESHOLD]
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(h[["window", "scenario", "variant", "POD", "POFD", "PSS", "CSI", "timing_r", "space_r",
                 "model_cells", "observed_cells", "rmse_interior_m_yr", "bias_interior_m_yr"]]
              .round(3).to_string(index=False))


# Validation of the storm sharing

# Check the storm sharing against a Barrier3D replay of the same domain-years
def validate(n_per_window=12, seed=0):
    sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "figure_making" / "model"))
    from storm_replay import replay
    rng = np.random.default_rng(seed)
    rows = []
    for w in WINDOWS:
        for v in ("drop72", "full"):
            d = (MATRIX / wtag(w) / "edgeBE" / CONTROLS[(w[0], "natural")]) if v == "drop72" else run_path(v, "natural", w)
            c = load_state(d)
            summ = pd.read_csv(summary_csv(w, v))
            shares = storm_shares(c, summ)
            keys = {}
            for (p, i), vol in shares.items():
                keys.setdefault((p, int(summ.at[i, "time"])), []).append(i)
            multi = [k for k, ii in keys.items() if len(ii) > 1]
            pick = [multi[j] for j in rng.choice(len(multi), size=min(n_per_window, len(multi)), replace=False)]
            for p, t in pick:
                b = c.barrier3d[p]
                year = summ[summ.time == t]
                storms = [(r.Rhigh * 10, r.Rlow * 10, r.period, int(r.duration)) for r in year.itertuples()]
                got, _ = replay(b, t, storms)
                for i, g in zip(year.index, got):
                    rows.append(dict(window=wtag(w), variant=v, pad=p, t=t, storm=i,
                                     replay_m3_per_m=g["owloss"], shared_m3_per_m=shares.get((p, i), 0.0)))
            print(f"  validated {wtag(w)} {v}: {len(pick)} domain-years", flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(EXP_DIR / "tables" / "sharing_validation.csv", index=False)
    tot = t.groupby(["window", "variant", "pad", "t"]).agg(r=("replay_m3_per_m", "sum"))
    t = t.join(tot, on=["window", "variant", "pad", "t"])
    t["replay_frac"] = t.replay_m3_per_m / t.r.where(t.r > 0)
    t["shared_frac"] = t.shared_m3_per_m / t.groupby(["window", "variant", "pad", "t"]).shared_m3_per_m.transform("sum")
    same_storm = t.loc[t.groupby(["window", "variant", "pad", "t"]).replay_m3_per_m.idxmax()]
    ok = (same_storm.shared_frac > 0).mean()
    err = (t.assign(d=(t.replay_frac - t.shared_frac).abs())
           .groupby(["window", "variant", "pad", "t"]).d.sum() / 2)
    print(f"storm carrying most replayed overwash also gets a share: {ok:.0%}; "
          f"volume assigned to the wrong storms: median {err.median():.0%}, 90th pct {err.quantile(0.9):.0%}")


# Stage 2: the ends

# be_edge_domain_solve.py, loaded by path
def _solver():
    path = PROJECT_ROOT / "scripts" / "input_prep" / "7-source-sink" / "2-calibrate" / "be_edge_domain_solve.py"
    spec = importlib.util.spec_from_file_location("be_edge_domain_solve", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# Where the solved end rates are kept
def ends_file():
    return EXP_DIR / "tables" / "ends.json"


# The solved end rates so far
def load_ends():
    return json.loads(ends_file().read_text(encoding="utf-8")) if ends_file().exists() else {}


# Secant solve of the GIS 1 and 90 end rates on LRR, for one variant
def solve_ends(w, v, tol=0.05, max_probes=6):
    from site_layer.hatteras_site_config import HATTERAS_BE_EDGE_ONLY
    solver = _solver()
    target = solver.load_target(*w)
    x = {g: float(HATTERAS_BE_EDGE_ONLY[w[0]][g]) for g in (1, 90)}
    hist = {1: [], 90: []}
    best = None
    for k in range(max_probes):
        override = f"1={x[1]:.4f},90={x[90]:.4f}"
        group = f"ends_{v}/probe{k}_full_management"
        launch(w, v, "full_management", group, override)
        d = run_path(v, "full_management", w, ends=f"probe{k}")
        model = solver.read_model(d / "tables" / "shoreline_change_rate.csv")
        res = {g: model[g] - target[g] for g in (1, 90)}
        size = max(abs(res[1]), abs(res[90]))
        if best is None or size < best[0]:
            best = (size, dict(x), dict(res), k)
        print(f"  {wtag(w)} {v} probe {k}: ends {x[1]:+.3f} / {x[90]:+.3f}  residual "
              f"{res[1]:+.3f} / {res[90]:+.3f}", flush=True)
        if size <= tol:
            break
        for g in (1, 90):
            hist[g].append((x[g], model[g]))
            if len(hist[g]) >= 2 and hist[g][-1][0] != hist[g][-2][0]:
                (x0, y0), (x1, y1) = hist[g][-2], hist[g][-1]
                gain = (y1 - y0) / (x1 - x0)
                if abs(gain) < 1e-3:
                    gain = solver.NOMINAL_GAIN
            else:
                gain = solver.NOMINAL_GAIN
            x[g] = x[g] - res[g] / gain
    return dict(ends=best[1], residual=best[2], probe=best[3], converged=best[0] <= tol)


# Solve the end rates for each variant and window, in parallel
def ends(variants, workers):
    solved = load_ends()
    jobs = [(w, v) for w in WINDOWS for v in variants]

    def one(job):
        w, v = job
        r = solve_ends(w, v)
        return w, v, r

    with ThreadPoolExecutor(workers) as ex:
        for w, v, r in ex.map(one, jobs):
            solved.setdefault(wtag(w), {})[v] = r
    ends_file().parent.mkdir(parents=True, exist_ok=True)
    ends_file().write_text(json.dumps(solved, indent=2), encoding="utf-8")
    # both scenarios at the solved ends
    jobs = [(w, v, s) for w in WINDOWS for v in variants for s in SCENARIOS]

    def final(job):
        w, v, s = job
        e = solved[wtag(w)][v]["ends"]
        return launch(w, v, s, f"ends_{v}/solved_{s}", f"1={e['1'] if '1' in e else e[1]:.4f},"
                      f"90={e['90'] if '90' in e else e[90]:.4f}")

    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(final, jobs))


# Figures

LENGTH = {"drop72": None, **{f"trim{L}": L for L in TRIMS}, "full": 240}


# One stage's score and per-image figures
def figures(stage="1"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import (apply_style, C, C_1984, C_1997, INK_MUTED, figsize, save,
                                             record_caption, _title, open_frame)
    apply_style()
    t = pd.read_csv(EXP_DIR / "tables" / f"stage{stage}_scores.csv")
    t = t[t.threshold == THRESHOLD]
    cells = pd.read_csv(EXP_DIR / "tables" / f"stage{stage}_cells.csv", parse_dates=["image"])
    lines = {("1996_2010", "full_management"): (C_1984, "-"), ("1996_2010", "natural"): (C_1984, (0, (3, 2))),
             ("2010_2024", "full_management"): (C_1997, "-"), ("2010_2024", "natural"): (C_1997, (0, (3, 2)))}
    xs = [L for L in TRIMS] + [240]
    fig, axes = plt.subplots(2, 3, figsize=figsize("double", height=5.2), constrained_layout=True)
    panels = [("PSS", "Peirce skill (POD − POFD)"), ("POD", "hit rate (POD)"), ("POFD", "false-alarm rate (POFD)"),
              ("timing_r", "timing r (per image)"), ("space_r", "space r (per domain)"),
              ("rmse_interior_m_yr", "interior RMSE (m/yr)")]
    for k, (col, label) in enumerate(panels):
        ax = axes.flat[k]
        for (win, scen), (colr, ls) in lines.items():
            d = t[(t.window == win) & (t.scenario == scen)].set_index("variant")
            y = [d.at[f"trim{L}", col] if f"trim{L}" in d.index else np.nan for L in TRIMS] +                 [d.at["full", col] if "full" in d.index else np.nan]
            ax.plot(xs, y, color=colr, ls=ls, lw=1.2, marker="o", ms=3,
                    label=f"{win.replace('_', '–')} {'managed' if scen == 'full_management' else 'natural'}")
            if "drop72" in d.index:
                ax.axhline(d.at["drop72", col], color=colr, ls=ls, lw=0.6, alpha=0.6)
        ax.set_xscale("log")
        ax.set_xticks(xs)
        ax.set_xticklabels([str(L) for L in TRIMS] + ["full"], fontsize=7)
        ax.minorticks_off()
        if k >= 3:
            ax.set_xlabel("storm events trimmed to (h)")
        open_frame(ax)
        _title(ax, k, label)
    axes.flat[0].legend(frameon=False, fontsize=6.5, loc="best")
    out = save(fig, EXP_DIR / "figures" / f"stage{stage}_scores.png")
    plt.close(fig)
    record_caption(out[0],
        "Storm-length candidates against the observed overwash record and CoastSat. Every candidate keeps "
        "every storm event and trims those longer than L hours to the L hours around their peak; 'full' "
        "trims nothing (a 240 h limit). Thin horizontal lines: the committed series, which drops events "
        "over 72 h. Red 1996-2010, blue 2010-2024; solid managed (full_management), dashed natural. "
        "Overwash scores use every image x domain cell the imagery assessed; a model year's overwash is "
        "shared among the storms that reached a dune gap and dated by their end (7-day grace). "
        f"{'Edge rates as the matrix (solved on the committed series), so RMSE is not yet a fair comparison.' if stage == '1' else 'Edge rates re-solved on each candidate.'}")

    # per-image counts, managed, for the committed series and the candidates
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=3.2), constrained_layout=True)
    for ax, win in zip(axes, ("1996_2010", "2010_2024")):
        d = cells[(cells.window == win) & (cells.scenario == "full_management")].dropna(subset=["observed"])
        imgs = sorted(d.image.unique())
        x = np.arange(len(imgs))
        obs_n = [d[(d.image == i) & (d.variant == "drop72")].observed.sum() for i in imgs]
        ax.bar(x, obs_n, color="0.8", width=0.8, label="observed")
        for v, colr in (("drop72", INK_MUTED), ("trim48", C["ADDED"]), ("trim72", C["ACCENT"]), ("full", C_1984)):
            m = [(d[(d.image == i) & (d.variant == v)].model_m3_per_m > THRESHOLD).sum() for i in imgs]
            if any(m):
                ax.plot(x, m, marker="o", ms=3, lw=1.1, color=colr, label=v)
        ax.set_xticks(x)
        ax.set_xticklabels([pd.Timestamp(i).strftime("%Y-%m") for i in imgs], rotation=60, fontsize=7)
        ax.set_ylabel("domains with overwash")
        open_frame(ax)
        _title(ax, list(axes).index(ax), f"{win.replace('_', '–')}, managed")
    axes[0].legend(frameon=False, fontsize=7)
    out = save(fig, EXP_DIR / "figures" / f"stage{stage}_per_image.png")
    plt.close(fig)
    record_caption(out[0],
        "Domains with overwash in each image (grey bars, observed; only domains the image assessed) against "
        "the model's count for the same image window, managed runs, for the committed series and three "
        "candidates. Timing agreement is what the per-image correlation in stage1_scores measures.")


# Run: the subprocess entry, or the action asked for
def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        MD._launch(sys.argv[2], sys.argv[3])
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["build", "run", "validate", "score", "ends", "score2", "figures", "figures2"])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--variants", nargs="*", default=[])
    a = ap.parse_args()
    if a.action == "build":
        build()
    elif a.action == "run":
        run(a.workers)
    elif a.action == "validate":
        validate()
    elif a.action == "score":
        score("1")
    elif a.action == "ends":
        ends(a.variants, a.workers)
    elif a.action == "score2":
        score("2")
    elif a.action == "figures":
        figures("1")
    elif a.action == "figures2":
        figures("2")


if __name__ == "__main__":
    main()
