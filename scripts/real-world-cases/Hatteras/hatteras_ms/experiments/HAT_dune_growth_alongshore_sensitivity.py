"""
Does an alongshore dune growth pattern, measured from the surveys, change the model's skill?

    python scripts/hatteras_ms/experiments/HAT_dune_growth_alongshore_sensitivity.py pattern
    python scripts/hatteras_ms/experiments/HAT_dune_growth_alongshore_sensitivity.py run
    python scripts/hatteras_ms/experiments/HAT_dune_growth_alongshore_sensitivity.py score

Why (2026-10-21, Hannah): the survey fit (HAT_dune_growth_alongshore_survey_fit.py
--widened --intervals) found an alongshore pattern in the dune growth rate that
repeats between the windows once smoothed (7-domain running medians r 0.74),
but nothing physical explains it (beach width runs the wrong way;
HAT_dune_growth_beach_width_check.py). The island-wide rate is kept; this is
the sensitivity test that says how much the choice matters.

Pattern, per window: the survey-fit domain medians, with the domains a fill
reached inside the window and Rodanthe / the S-curves (GIS 76-88, rebuilt by
the road managers) blanked, the gaps filled linearly from their neighbours,
then the 7-domain running median. Divided by its island mean, it is a factor f
per domain (mean 1), so only the alongshore distribution changes, not the
island-wide growth: rmin, rmax = f * the period's adopted bounds.

Variants (each both windows, full management, edgeBE, groin on, relocations
off, the runner's settings otherwise):
    island          no override: the adopted rate everywhere (rerun here so all three share one runner state)
    own_pattern     each window's own pattern
    shared_pattern  the mean of the two windows' factors, the same in both windows
Bounds reach the runner through HAT_GROWTH_OVERRIDE (added 2026-10-21; inert
when unset). Scores: interior shoreline RMSE and bias, dune crest against each
survey (island median and 10th percentile, and how well the domain medians
follow the surveyed ones alongshore), overwash skill.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-21
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(_p for _p in _HERE.parents  if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
HINDCAST = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"
sys.path[:0] = [str(_HAT_DIR), str(_HERE.parent)]

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-21-alongshore-growth-sensitivity"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
FIT_TABLES = (PROJECT_ROOT / "output" / "raw_runs" / "experiments" / "dunes"
              / "2026-10-21-alongshore-dune-growth-survey-fit" / "tables")
FIT_SUFFIX = "_widened_at_edge_intervals"
STARTS = {1996: "1996_2009", 2009: "2009_2025"}
ADOPTED = {1996: (0.01, 0.03), 2009: (0.02, 0.04)}      # the runner's DUNE_GROWTH_RATES
REBUILT = range(76, 89)
VARIANTS = ("island", "own_pattern", "shared_pattern")
BERM_NAVD88 = 1.7
WORKERS = 4
# -----------------------------------------------------------------------------


# The survey-fit factor per GIS domain 1-90, mean 1
def factor(start):
    t = pd.read_csv(FIT_TABLES / f"domain_rates_{start}_start{FIT_SUFFIX}.csv").set_index("gis")
    r = t.r_median.where(~t.fill_in_window & ~t.index.isin(REBUILT))
    r = r.interpolate(limit_direction="both")
    s = r.rolling(7, center=True, min_periods=4).median()
    return s / s.mean()


def factors():
    f = {s: factor(s) for s in STARTS}
    shared = (f[1996] + f[2009]) / 2
    return {"own_pattern": f, "shared_pattern": {s: shared / shared.mean() for s in STARTS}}


def override(variant, start):
    if variant == "island":
        return ""
    f = factors()[variant][start]
    lo, hi = ADOPTED[start]
    return ",".join(f"{g}={lo * f[g]:.5f}:{hi * f[g]:.5f}" for g in f.index)


def pattern():
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    fs = factors()
    t = pd.DataFrame({f"{v}_{s}": fs[v][s] for v in fs for s in STARTS}).round(3)
    for s in STARTS:
        lo, hi = ADOPTED[s]
        t[f"own_pattern_{s}_rmin_rmax"] = [f"{lo * x:.4f}-{hi * x:.4f}" for x in fs["own_pattern"][s]]
    t.index.name = "gis"
    t.to_csv(EXP_DIR / "tables" / "growth_factor_per_domain.csv")
    with pd.option_context("display.width", 200, "display.max_rows", 100):
        print(t.to_string())
    return t


def run_dir(variant, start):
    base = EXP_DIR / variant / STARTS[start] / "edgeBE"
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


def launch(job):
    variant, start = job
    if run_dir(variant, start) is not None:
        return job, "done"
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(start), "HAT_SOURCE_SINK_PRESET": "edgeBE",
                "HAT_SCENARIO": "full_management", "HAT_RELOCATIONS": "False", "HAT_GROIN_ENABLED": "True",
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/{variant}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8", "MPLBACKEND": "Agg"})
    ov = override(variant, start)
    if ov:
        env["HAT_GROWTH_OVERRIDE"] = ov
    t0 = time.time()
    with open(logs / f"{variant}_{start}.log", "w", encoding="utf-8") as fh:
        code = subprocess.run([sys.executable, str(HINDCAST)], env=env, cwd=HINDCAST.parent,
                              stdout=fh, stderr=subprocess.STDOUT).returncode
    msg = f"{variant} {start}: exit {code}, {(time.time() - t0) / 60:.1f} min"
    print(msg, flush=True)
    return job, msg


# overwash_vs_model's comparison on this run's overwash (as HAT_low_growth_test_period.cells does), then its skill
def overwash_skill(c, obs, start, end):
    import overwash_vs_model as ovm
    import storms_vs_overwash as SV
    q = np.array([np.asarray(b.QowTS) for b in c.barrier3d]).T
    orig = ovm.load_qow
    ovm.load_qow = lambda w, arm: (q, "run")
    try:
        d = ovm.compare((start, end), "managed", obs, ovm.HEADLINE_THR).dropna(subset=["observed"])
    finally:
        ovm.load_qow = orig
    o, m = d.observed.astype(int), d.model.astype(int)
    d["outcome"] = np.select([(o == 1) & (m == 1), (o == 1) & (m == 0), (o == 0) & (m == 1)],
                             ["hit", "miss", "false_alarm"], "neither")
    return SV.skill(d)


def score():
    import HAT_dune_ceiling_uniform as U
    import HAT_dune_crest_vs_surveys as V
    from site_layer import hat_overwash as ow
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    rows = []
    for start in STARTS:
        first = V.PERIODS[start][1]
        sv = V.surveys(start)
        sv = sv[sv.survey != start]
        for variant in VARIANTS:
            rd = run_dir(variant, start)
            if rd is None:
                print(f"  missing {variant} {start}")
                continue
            meta = json.loads(next(rd.glob("*_run_metadata.json")).read_text(encoding="utf-8"))
            if (variant == "island") != (meta["dunes and roadway"].get("growth_override", "none") == "none"):
                raise SystemExit(f"{rd}: growth override does not match the variant")
            c = U.load(rd)
            crest = np.stack(U.cell_crests(c, DOM), axis=1) + BERM_NAVD88          # years x 90 x 50
            med_err, p10_err, along_r = [], [], []
            for year, g in sv.groupby("survey"):
                k = int(g.model_state.iloc[0][:4]) - first
                if not 0 < k < crest.shape[0]:
                    continue
                o = g.set_index(["gis", "cell"]).crest.reindex(
                    pd.MultiIndex.from_product([range(1, 91), range(50)])).to_numpy().reshape(90, 50)
                m = crest[k]
                ok = np.isfinite(o)
                med_err.append(np.median(m[ok]) - np.median(o[ok]))
                p10_err.append(np.percentile(m[ok], 10) - np.percentile(o[ok], 10))
                full = ok.sum(1) >= 25
                along_r.append(np.corrcoef(np.nanmedian(o[full], 1), np.median(m[full], 1))[0, 1])
            sk = overwash_skill(c, obs, start, U.PERIODS[start][0])
            rows.append(dict(window=STARTS[start], variant=variant,
                             shoreline_rmse=round(float(meta["skill"]["rmse_interior_m_yr"]), 3),
                             shoreline_bias=round(float(meta["skill"]["mean_bias_interior_m_yr"]), 3),
                             crest_median_rmse=round(float(np.sqrt(np.mean(np.square(med_err)))), 2),
                             crest_p10_rmse=round(float(np.sqrt(np.mean(np.square(p10_err)))), 2),
                             crest_alongshore_r=round(float(np.mean(along_r)), 2),
                             overwash_skill=round(sk["skill"], 3), overwash_hits=sk["hits"],
                             overwash_false_alarms=sk["false_alarms"], run=str(rd.relative_to(EXP_DIR))))
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "alongshore_growth_sensitivity_scores.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 20):
        print(t.drop(columns="run").to_string(index=False))


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "score"
    if cmd == "pattern":
        pattern()
    elif cmd == "run":
        pattern()
        jobs = [(v, s) for s in STARTS for v in VARIANTS]
        with ThreadPoolExecutor(WORKERS) as ex:
            list(ex.map(launch, jobs))
        score()
    else:
        score()


if __name__ == "__main__":
    main()
