"""
Which dune ceiling and growth rate together reproduce the 1996-2009 dunes and overwash best?

    python scripts/hatteras_ms/experiments/HAT_dune_ceiling_growth_grid.py run      # all stages, in order
    python scripts/hatteras_ms/experiments/HAT_dune_ceiling_growth_grid.py score    # rank from runs on disk

Stage 1, grid: Dmaxel 5.0 / 5.5 / 6.0 / 6.5 m NAVD88 x rmin-rmax 0.15-0.30 /
0.25-0.45 / 0.35-0.55, plus the control (5.5 m, 0.55-0.95), 1996-2009 only.
Stage 2, seeds: the top three by skill and the control, each re-run with the
growth rates drawn from three other seeds, to measure how much the random draw
alone moves the scores. Stage 3, selection on 1996-2009 only (Hannah, 2026-10-19):
overwash skill first; candidates within the seed noise of the best are a tie,
broken by the modelled 2009 crest against the lidar (median offset, then the
alongshore r); the shoreline RMSE must not exceed the control's by more than
GUARD_RMSE. Stage 4: the winner and the control run once on 2009-2025, and both
periods with the groin on. Managed edgeBE, in-process only.
Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-19
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_ceiling_uniform as U  # noqa: E402
import HAT_dune_ceiling_dem_growth as DG  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-19-dune-ceiling-growth-grid"
EXP_DIR = U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
CEILINGS = (5.0, 5.5, 6.0, 6.5)                       # m NAVD88
RATES = ((0.15, 0.30), (0.25, 0.45), (0.35, 0.55))    # rmin, rmax
CONTROL = (5.5, (0.55, 0.95))                         # the parameter file as adopted
TOP_N = 3                                             # grid candidates carried to the seed stage
SEEDS = (1, 2, 3)                                     # 0 is Barrier3D's own draw (seed 1973)
TIE_SD = 2.0                                          # within this many seed SDs of the best skill = a tie
GUARD_RMSE = 0.005                                    # m/yr the shoreline RMSE may exceed the control's
WORKERS = 4
# -----------------------------------------------------------------------------


def name(ceiling, rates, seed=0, groin=False):
    base = "control" if (ceiling, rates) == CONTROL else \
        f"c{round(ceiling * 100)}_r{round(rates[0] * 100):03d}_{round(rates[1] * 100):03d}"
    return f"{base}_s{seed}" + ("_groin" if groin else "")


def parse(variant):
    m = re.fullmatch(r"(?:control|c(\d+)_r(\d+)_(\d+))_s(\d+)(_groin)?", variant)
    if m.group(1) is None:
        ceiling, rates = CONTROL
    else:
        ceiling, rates = int(m.group(1)) / 100, (int(m.group(2)) / 100, int(m.group(3)) / 100)
    return ceiling, rates, int(m.group(4)), bool(m.group(5))


# Barrier3D's draw (load_input.py), one RNG per domain with the same seed, so only the seed changes
def redraw_growth(cascade, rates, seed):
    lo, hi = rates
    for b in cascade.barrier3d:
        g = b._growthparam
        g[...] = lo + (hi - lo) * np.random.default_rng(seed=seed).random(g.shape)   # in place: the managers share it


# Child process: set rates, seed and ceiling on the built Cascade, then run the unchanged hindcast
def _launch(variant):
    import runpy
    import cascade_pipeline.hindcast as H
    ceiling, rates, seed, _ = parse(variant)
    orig = H.build_cascade

    def build_cascade(*a, **k):
        k["rmin"] = [rates[0]] * len(k["rmin"])
        k["rmax"] = [rates[1]] * len(k["rmax"])
        c = orig(*a, **k)
        if seed:
            redraw_growth(c, rates, seed)
        if ceiling != CONTROL[0]:
            U.set_uniform(c, ceiling)
        print(f"EXPERIMENT dune ceiling {ceiling:.2f} m NAVD88, rmin/rmax {rates[0]}/{rates[1]}, "
              f"growth-rate seed {seed or 'Barrier3D default (1973)'}", flush=True)
        return c

    H.build_cascade = build_cascade
    sys.path.insert(0, str(U.HINDCAST.parent))
    sys.argv = [str(U.HINDCAST)]
    runpy.run_path(str(U.HINDCAST), run_name="__main__")


def run_dir(start, variant):
    base = EXP_DIR / "runs" / variant / f"{start}_{U.PERIODS[start][0]}" / U.PRESET
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


def launch(job):
    start, variant = job
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(start),
                "HAT_SOURCE_SINK_PRESET": U.PRESET, "HAT_SCENARIO": U.SCENARIO,
                "HAT_GROIN_ENABLED": str(parse(variant)[3]),
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{variant}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8", "MPLBACKEND": "Agg"})
    t0 = time.time()
    with open(logs / f"{variant}_{start}.log", "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", variant],
                           env=env, cwd=U.HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    rec = dict(start=start, variant=variant, returncode=p.returncode, minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{start} {variant:22s} exit {p.returncode}  {rec['minutes']} min", flush=True)


def run_all(jobs):
    jobs = [j for j in jobs if run_dir(*j) is None]
    with ThreadPoolExecutor(WORKERS) as ex:
        list(ex.map(launch, jobs))


# Scores

_LIDAR = {}


# 2009 starting crest, per-domain median, m above the berm: the 1996 run's end is compared with it
def lidar():
    if not _LIDAR:
        sys.path.insert(0, str(U.PROJECT_ROOT / "scripts" / "figure_making" / "pipeline" / "1-barrier3d-domains"))
        import dune_crest_figures as DC
        d = DC.datum()
        c, *_ = DC.crests(2009, d["berm_mhw"])
        _LIDAR["med"] = np.median(c - d["berm_mhw"], axis=1)
    return _LIDAR["med"]


def metrics(start, variant):
    from site_layer import hat_overwash as ow
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    rd = run_dir(start, variant)
    if rd is None:
        return None
    c = U.load(rd)
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    sk = DG.overwash_skill(c, obs, DOM, U.PERIODS[start][1])
    meta = json.loads((rd / f"{rd.name}_run_metadata.json").read_text(encoding="utf-8"))
    ceiling, rates, seed, groin = parse(variant)
    cc = U.cell_crests(c, DOM)
    row = dict(period=f"{start}_{U.PERIODS[start][0]}", variant=variant, ceiling_m_navd88=ceiling,
               rmin=rates[0], rmax=rates[1], seed=seed, groin=groin,
               skill=round(sk["skill"], 3), hits=sk["hits"], misses=sk["misses"], false_alarms=sk["false_alarms"],
               rmse_interior_m_yr=float(meta["skill"]["rmse_interior_m_yr"]),
               dunes_rebuilt=int(pd.read_csv(rd / "tables" / "road_management.csv").dunes_rebuilt.sum()))
    if start == 1996:
        med_end = np.array([np.median(x[-1]) for x in cc])
        row.update(crest_minus_lidar_m=round(float(np.median(med_end - lidar())), 3),
                   crest_vs_lidar_r=round(float(np.corrcoef(med_end, lidar())[0, 1]), 3))
    return row


def score_runs(jobs):
    rows = [m for m in (metrics(*j) for j in jobs) if m]
    return pd.DataFrame(rows)


def grid_jobs():
    return [(1996, name(cl, r)) for cl in CEILINGS for r in RATES] + [(1996, name(*CONTROL))]


def seed_jobs(top):
    return [(1996, name(cl, r, s)) for cl, r in top for s in SEEDS]


# Stage 3: the agreed ranking, on 1996-2009 only
def select(grid, seeds):
    allr = pd.concat([grid, seeds], ignore_index=True)
    allr["combo"] = allr.variant.str.replace(r"_s\d+$", "", regex=True)
    g = allr.groupby("combo").agg(ceiling_m_navd88=("ceiling_m_navd88", "first"), rmin=("rmin", "first"),
                                  rmax=("rmax", "first"), n_seeds=("seed", "size"),
                                  skill_mean=("skill", "mean"), skill_sd=("skill", "std"),
                                  hits_mean=("hits", "mean"), false_alarms_mean=("false_alarms", "mean"),
                                  crest_minus_lidar_m=("crest_minus_lidar_m", "mean"),
                                  crest_vs_lidar_r=("crest_vs_lidar_r", "mean"),
                                  rmse_mean=("rmse_interior_m_yr", "mean"))
    seeded = g[g.n_seeds > 1]
    noise = float(np.sqrt((seeded.skill_sd ** 2).mean())) if len(seeded) else np.nan   # pooled seed SD
    control_rmse = g.loc[name(*CONTROL)[:-3], "rmse_mean"]
    g["passes_guard"] = g.rmse_mean <= control_rmse + GUARD_RMSE
    cand = g[(g.n_seeds > 1) & g.passes_guard].copy()
    best = cand.skill_mean.max()
    cand["tied_with_best"] = cand.skill_mean >= best - TIE_SD * noise
    tied = cand[cand.tied_with_best].copy()
    tied["abs_offset"] = tied.crest_minus_lidar_m.abs()
    tied = tied.sort_values(["abs_offset", "crest_vs_lidar_r"], ascending=[True, False])
    winner = tied.index[0]
    g["seed_noise_sd"] = noise
    g["tied_with_best"] = cand.tied_with_best.reindex(g.index).fillna(False)
    g["winner"] = g.index == winner
    return g.sort_values("skill_mean", ascending=False).round(4), winner, noise


def write(df, fname):
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    df.to_csv(EXP_DIR / "tables" / fname, index=bool(df.index.name))


def stages(run):
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        if run:
            run_all(grid_jobs())
        grid = score_runs(grid_jobs())
        write(grid, "stage1_grid_1996.csv")
        print(grid.to_string(index=False), flush=True)
        ranked = grid[grid.variant != name(*CONTROL)].sort_values("skill", ascending=False)
        top = [parse(v)[:2] for v in ranked.variant.head(TOP_N)] + [CONTROL]
        if run:
            run_all(seed_jobs(top))
        seeds = score_runs(seed_jobs(top))
        write(seeds, "stage2_seeds_1996.csv")
        sel, winner, noise = select(grid, seeds)
        write(sel, "stage3_selection_1996.csv")
        print(sel.to_string(), flush=True)
        print(f"seed noise (pooled SD of skill) {noise:.3f}; winner {winner}", flush=True)
        wc, wr, _, _ = parse(winner + "_s0")
        final = [(2009, name(wc, wr)), (2009, name(*CONTROL)),
                 (1996, name(wc, wr, groin=True)), (2009, name(wc, wr, groin=True)),
                 (1996, name(*CONTROL, groin=True)), (2009, name(*CONTROL, groin=True))]
        if run:
            run_all(final)
        fin = score_runs(final + [(1996, name(wc, wr)), (1996, name(*CONTROL))])
        write(fin, "stage4_winner_vs_control.csv")
        print(fin.to_string(index=False), flush=True)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(sys.argv[2])
        return
    stages(run=len(sys.argv) > 1 and sys.argv[1] == "run")


if __name__ == "__main__":
    main()
