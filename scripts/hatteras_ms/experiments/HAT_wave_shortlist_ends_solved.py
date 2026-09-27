"""Wave shortlist with the end domains solved per setting (2026-09-26).

Hannah, 2026-09-26: "what about when you solve for the ends, are there
different wave parameters that perform the best?" The 2026-09-25 wave grid
ran zeroBE (nothing imposed at GIS 1 or 90), and the stored edgeBE values
(HATTERAS_BE_EDGE_ONLY) were solved at Hs 2.5 on the old /10 offset, so they
do not carry over. What the ends must carry depends on the waves, so each
wave setting gets its own solve. Chosen with Hannah:

    shortlist  the top 10 zeroBE settings (smoothed score) per window x
               scenario from wave-climate/2026-09-25-wave-grid-smoothed-score
    scenarios  natural and full management, both windows (40 chains)
    ends       solved against each window's CoastSat LRR, as the matrix end
               values were: GIS 1 against the raw domain mean, GIS 90 against
               the LOESS-10 value (the target table's own splice)
    solve      step 0 is the setting's zeroBE grid run (ends 0, 0). Each end
               is stepped on its own (they are 89 domains apart): step 1 from
               the 2026-09-11 response (about 0.09 m/yr of residual per m/yr
               imposed at GIS 1, 0.13 at GIS 90), then the secant through the
               last two probes. Lockstep: every chain runs its next probe
               before any is solved again. Converged at |residual| <= 0.02
               m/yr at both ends, or stop after MAX_STEPS probes
    score      the converged (or last) run, as the grid: share of the
               alongshore variation explained by the model smoothed like the
               CoastSat target, interior GIS 2-89; bias, r, the raw score, the
               imposed ends and their residuals beside it
    fixed      metres offset (dune line v1), edgeBE preset with HAT_BE_OVERRIDE,
               no groin, no relocations, Barrier3D route_overwash fix

WHERE: output/raw_runs/experiments/wave-climate/2026-09-26-wave-shortlist-ends-solved/
    README.md, tables/{shortlist,solve_log,all_runs}.csv, figures/,
    logs/<scenario>/step<k>/<period>_<settings>.log
    runs/<scenario>_step<k>/<period>/edgeBE/<run_name>/     (on disk only)

    python scripts/hatteras_ms/experiments/HAT_wave_shortlist_ends_solved.py run [--jobs 8]
    python scripts/hatteras_ms/experiments/HAT_wave_shortlist_ends_solved.py score
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_wave_grid_smoothed_score as grid  # noqa: E402

common, step2 = grid.common, grid.step2
TAG = "wave-climate/2026-09-26-wave-shortlist-ends-solved"
STUDY_DIR = grid.RAW_RUNS / "experiments" / TAG
TABLES_DIR, LOGS_DIR = STUDY_DIR / "tables", STUDY_DIR / "logs"
KEYS = list(grid.KEYS)
TOP_N = 10
MAX_STEPS = 5
TOL = 0.02                                 # m/yr, both ends
FIRST_GAIN = {1: 0.09, 90: 0.13}           # d(residual)/d(imposed), the 09-11 solve
ENDS = (1, 90)


def shortlist():
    import pandas as pd
    t = pd.read_csv(grid.TABLES_DIR / "all_runs.csv")
    t = t[t.status == "scored"].drop_duplicates(["scenario", "period_start", *KEYS])
    out = []
    for (sc, p), x in t.groupby(["scenario", "period_start"]):
        top = x.nlargest(TOP_N, "smoothed_variance_explained").copy()
        top["zerobe_rank"] = range(1, len(top) + 1)
        out.append(top)
    s = pd.concat(out)[["scenario", "period_start", *KEYS, "zerobe_rank",
                        "smoothed_variance_explained", "bias_m_yr", "run_dir"]]
    return s.rename(columns={"smoothed_variance_explained": "zerobe_smoothed_variance_explained",
                             "bias_m_yr": "zerobe_bias_m_yr", "run_dir": "zerobe_run_dir"})


def label(s):
    return grid.label({k: float(s[k]) for k in KEYS})


def residuals(run_dir, period, targets):
    rates = common.run_rates(run_dir)
    return {g: float(rates[g] - targets[period][g]) for g in ENDS}


def log_path(sc, step, period, s):
    return LOGS_DIR / sc / f"step{step}" / f"{period}_{label(s)}.log"


def run_env(sc, step, period, s, ends):
    e = grid.run_env("x", sc, period, {k: float(s[k]) for k in KEYS})
    e["HAT_SOURCE_SINK_PRESET"] = "edgeBE"
    e["HAT_BE_OVERRIDE"] = ",".join(f"{g}={ends[g]:.4f}" for g in ENDS)
    e["HAT_RUN_TAG"] = f"{TAG}/runs/{sc}_step{step}"
    return e


def run_dir_of(sc, step, period):
    return STUDY_DIR / "runs" / f"{sc}_step{step}" / grid.window(period) / "edgeBE"


def launch(job):
    sc, step, period, s, ends = job
    log = log_path(sc, step, period, s)
    if grid.finished(log):
        return
    log.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    p = subprocess.run([sys.executable, str(grid.HINDCAST)], env=run_env(sc, step, period, s, ends),
                       cwd=str(grid.PROJECT_ROOT), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=grid.RUN_TIMEOUT_S)
    log.write_text(f'HAT_BE_OVERRIDE="{run_env(sc, step, period, s, ends)["HAT_BE_OVERRIDE"]}"\n'
                   + (p.stdout or "") + "\n--- STDERR ---\n" + (p.stderr or ""), encoding="utf-8")
    what = "done" if p.returncode == 0 else f"FAILED ({step2.stop_reason(log)})"
    print(f"{what} {sc} step{step} {grid.window(period)} {label(s)} "
          f"ends {ends[1]:+.2f}/{ends[90]:+.2f} in {(time.perf_counter() - t0) / 60:.1f} min",
          flush=True)


def find_run(sc, step, period, s):
    """The run a probe made: its folder under the step's tag, matched on the
    wave settings in its metadata (the run name leaves defaults out)."""
    for md in run_dir_of(sc, step, period).glob("*/*_run_metadata.json"):
        w = json.loads(md.read_text(encoding="utf-8"))["wave climate"]
        got = (float(w["wave_height_m"]), float(w["wave_period_s"]),
               float(w["wave_asymmetry"]), float(w["wave_angle_high_frac"]))
        if np.allclose(got, [float(s[k]) for k in KEYS]):
            return md.parent
    return None


def cmd_run(a):
    import pandas as pd
    grid.check_barrier3d()
    common.keep_awake()
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    sl = shortlist()
    sl.to_csv(TABLES_DIR / "shortlist.csv", index=False)
    targets = {p: common.coastsat_target(p) for p in grid.PERIODS}
    chains = []
    for _, r in sl.iterrows():
        d0 = grid.STUDY_DIR / r.zerobe_run_dir
        chains.append({"s": r, "sc": r.scenario, "p": int(r.period_start),
                       "hist": [({1: 0.0, 90: 0.0}, residuals(d0, int(r.period_start), targets), d0)],
                       "state": "live"})
    print(f"{len(chains)} chains, lockstep, {a.jobs} at a time", flush=True)
    log_rows = []
    for step in range(1, MAX_STEPS + 1):
        jobs = []
        for c in chains:
            if c["state"] != "live":
                continue
            e, r, _ = c["hist"][-1]
            if all(abs(r[g]) <= TOL for g in ENDS):
                c["state"] = "converged"
                continue
            nxt = {}
            for g in ENDS:
                if len(c["hist"]) == 1:
                    nxt[g] = e[g] - r[g] / FIRST_GAIN[g]
                else:
                    (e0, r0, _), (e1, r1, _) = c["hist"][-2], c["hist"][-1]
                    slope = (r1[g] - r0[g]) / (e1[g] - e0[g]) if e1[g] != e0[g] else FIRST_GAIN[g]
                    if not np.isfinite(slope) or abs(slope) < 1e-3:
                        slope = FIRST_GAIN[g]
                    nxt[g] = e1[g] - r1[g] / slope
            c["next"] = nxt
            jobs.append((c["sc"], step, c["p"], c["s"], nxt))
        if not jobs:
            break
        print(f"=== step {step}: {len(jobs)} probes", flush=True)
        with ThreadPoolExecutor(max_workers=a.jobs) as pool:
            list(pool.map(launch, jobs))
        for c in chains:
            if c["state"] != "live" or "next" not in c:
                continue
            d = find_run(c["sc"], step, c["p"], c["s"])
            if d is None:
                c["state"] = f"failed at step {step}: {step2.stop_reason(log_path(c['sc'], step, c['p'], c['s']))}"
                continue
            r = residuals(d, c["p"], targets)
            c["hist"].append((c.pop("next"), r, d))
            log_rows.append({"scenario": c["sc"], "period_start": c["p"], **{k: c["s"][k] for k in KEYS},
                             "step": step, "gis1_imposed": c["hist"][-1][0][1],
                             "gis90_imposed": c["hist"][-1][0][90],
                             "gis1_residual": r[1], "gis90_residual": r[90]})
        pd.DataFrame(log_rows).to_csv(TABLES_DIR / "solve_log.csv", index=False)
    for c in chains:
        if c["state"] == "live":
            e, r, _ = c["hist"][-1]
            c["state"] = "converged" if all(abs(r[g]) <= TOL for g in ENDS) else "not converged"
    json.dump([{"scenario": c["sc"], "period_start": c["p"], **{k: float(c["s"][k]) for k in KEYS},
                "state": c["state"], "steps": len(c["hist"]) - 1,
                "run_dir": str(c["hist"][-1][2]), "ends": c["hist"][-1][0],
                "residuals": c["hist"][-1][1]} for c in chains],
              (TABLES_DIR / "chains.json").open("w", encoding="utf-8"), indent=1)
    return cmd_score(a)


def cmd_score(_=None):
    import pandas as pd
    chains = json.loads((TABLES_DIR / "chains.json").read_text(encoding="utf-8"))
    sl = shortlist()
    targets = {p: common.coastsat_target(p) for p in grid.PERIODS}
    rows = []
    for c in chains:
        z = sl[(sl.scenario == c["scenario"]) & (sl.period_start == c["period_start"])
               & np.logical_and.reduce([np.isclose(sl[k], c[k]) for k in KEYS])].iloc[0]
        rec = {**{k: c[k] for k in ("scenario", "period_start", *KEYS, "state", "steps")},
               "gis1_imposed": c["ends"]["1"] if "1" in c["ends"] else c["ends"][1],
               "gis90_imposed": c["ends"]["90"] if "90" in c["ends"] else c["ends"][90],
               "gis1_residual": list(c["residuals"].values())[0],
               "gis90_residual": list(c["residuals"].values())[1],
               "zerobe_rank": int(z.zerobe_rank),
               "zerobe_smoothed_variance_explained": z.zerobe_smoothed_variance_explained,
               "zerobe_bias_m_yr": z.zerobe_bias_m_yr}
        d = Path(c["run_dir"])
        if c["steps"] > 0 and not c["state"].startswith("failed"):
            rec.update(grid.score_run(d, targets[c["period_start"]]))
            rec["run_dir"] = str(d.relative_to(STUDY_DIR)).replace("\\", "/")
        rows.append(rec)
    t = pd.DataFrame(rows)
    t["ends_solved_rank"] = t.groupby(["scenario", "period_start"]).smoothed_variance_explained.rank(
        ascending=False).astype("Int64")
    t = t.sort_values(["scenario", "period_start", "ends_solved_rank"])
    t.to_csv(TABLES_DIR / "all_runs.csv", index=False)
    pd.set_option("display.width", 250)
    cols = ["scenario", "period_start", *KEYS, "state", "steps", "gis1_imposed", "gis90_imposed",
            "zerobe_rank", "zerobe_smoothed_variance_explained", "ends_solved_rank",
            "smoothed_variance_explained", "bias_m_yr"]
    print(t[[c for c in cols if c in t]].round(3).to_string(index=False))
    return 0


# =============================================================================
# RESUME WITH A SAFEGUARDED STEP (added 2026-09-26)
# =============================================================================
# The first pass (plain secant, 5 probes) converged 6 of 40 chains. GIS 90
# was fine; GIS 1 in 2010-2024 is not smooth -- imposed 0-25 m/yr gives
# residuals near zero, anything above ~25 gives +5 to +15 whatever the value
# -- so the secant took steps to -466 and +335 m/yr, and four probes drowned
# the barrier. The resume keeps every probe already run and steps each end on
# its own:
#   bracketed  (a probe on each side of the target): interpolate between the
#              closest pair, held at least 10% inside it so it always shrinks
#   otherwise  secant through the two latest probes, capped at +-STEP_CAP
#   converged  an end within TOL keeps its value
# A drowned probe carries no residual and is left out of the history.
STEP_CAP = 30.0
EXTRA_STEPS = 4


def history(chains_row, log, targets):
    """[(ends, residuals)] for one chain: its zeroBE run, then every probe
    that produced a run (from solve_log.csv)."""
    s, sc, p = chains_row, chains_row.scenario, int(chains_row.period_start)
    h = [({1: 0.0, 90: 0.0}, residuals(grid.STUDY_DIR / s.zerobe_run_dir, p, targets))]
    x = log[(log.scenario == sc) & (log.period_start == p)
            & np.logical_and.reduce([np.isclose(log[k], s[k]) for k in KEYS])].sort_values("step")
    for _, r in x.iterrows():
        h.append(({1: r.gis1_imposed, 90: r.gis90_imposed}, {1: r.gis1_residual, 90: r.gis90_residual}))
    return h


def safeguarded_next(h, g):
    pts = [(e[g], r[g]) for e, r in h]
    e_last, r_last = pts[-1]
    if abs(r_last) <= TOL:
        return e_last
    neg = [pt for pt in pts if pt[1] < 0]
    pos = [pt for pt in pts if pt[1] > 0]
    if neg and pos:
        lo = max(neg, key=lambda pt: pt[1])          # closest to zero from below
        hi = min(pos, key=lambda pt: pt[1])          # closest to zero from above
        x = lo[0] - lo[1] * (hi[0] - lo[0]) / (hi[1] - lo[1])
        a, b = sorted((lo[0], hi[0]))
        return float(np.clip(x, a + 0.1 * (b - a), b - 0.1 * (b - a)))
    (e0, r0), (e1, r1) = pts[-2], pts[-1]
    slope = (r1 - r0) / (e1 - e0) if e1 != e0 else FIRST_GAIN[g]
    if not np.isfinite(slope) or slope <= 1e-3:      # the response is positive
        slope = FIRST_GAIN[g]
    return float(e1 - np.clip(r1 / slope, -STEP_CAP, STEP_CAP))


def cmd_resume(a):
    import pandas as pd
    grid.check_barrier3d()
    common.keep_awake()
    sl = pd.read_csv(TABLES_DIR / "shortlist.csv")
    log = pd.read_csv(TABLES_DIR / "solve_log.csv")
    targets = {p: common.coastsat_target(p) for p in grid.PERIODS}
    chains = []
    for _, s in sl.iterrows():
        h = history(s, log, targets)
        chains.append({"s": s, "sc": s.scenario, "p": int(s.period_start), "hist": h,
                       "state": "converged" if all(abs(h[-1][1][g]) <= TOL for g in ENDS) else "live"})
    start = int(log.step.max()) + 1
    print(f"resume: {sum(c['state'] == 'live' for c in chains)} of {len(chains)} chains live, "
          f"steps {start}-{start + a.extra - 1}", flush=True)
    rows = log.to_dict("records")
    for step in range(start, start + a.extra):
        jobs = []
        for c in chains:
            if c["state"] != "live":
                continue
            if all(abs(c["hist"][-1][1][g]) <= TOL for g in ENDS):
                c["state"] = "converged"
                continue
            c["next"] = {g: safeguarded_next(c["hist"], g) for g in ENDS}
            jobs.append((c["sc"], step, c["p"], c["s"], c["next"]))
        if not jobs:
            break
        print(f"=== step {step}: {len(jobs)} probes", flush=True)
        with ThreadPoolExecutor(max_workers=a.jobs) as pool:
            list(pool.map(launch, jobs))
        for c in chains:
            if "next" not in c:
                continue
            nxt = c.pop("next")
            d = find_run(c["sc"], step, c["p"], c["s"])
            if d is None:                            # drowned: no residual, keep going
                print(f"  no run: {c['sc']} {c['p']} {label(c['s'])} at {nxt}", flush=True)
                continue
            r = residuals(d, c["p"], targets)
            c["hist"].append((nxt, r))
            rows.append({"scenario": c["sc"], "period_start": c["p"], **{k: c["s"][k] for k in KEYS},
                         "step": step, "gis1_imposed": nxt[1], "gis90_imposed": nxt[90],
                         "gis1_residual": r[1], "gis90_residual": r[90]})
        pd.DataFrame(rows).to_csv(TABLES_DIR / "solve_log.csv", index=False)
    out = []
    for c in chains:
        e, r = c["hist"][-1]
        state = "converged" if all(abs(r[g]) <= TOL for g in ENDS) else "not converged"
        step = int(max([0] + [x["step"] for x in rows if x["scenario"] == c["sc"]
                              and x["period_start"] == c["p"]
                              and all(np.isclose(x[k], c["s"][k]) for k in KEYS)
                              and np.isclose(x["gis1_imposed"], e[1]) and np.isclose(x["gis90_imposed"], e[90])]))
        d = find_run(c["sc"], step, c["p"], c["s"]) if step else grid.STUDY_DIR / c["s"].zerobe_run_dir
        out.append({"scenario": c["sc"], "period_start": c["p"], **{k: float(c["s"][k]) for k in KEYS},
                    "state": state, "steps": step, "run_dir": str(d), "ends": {1: e[1], 90: e[90]},
                    "residuals": {1: r[1], 90: r[90]}})
    json.dump(out, (TABLES_DIR / "chains.json").open("w", encoding="utf-8"), indent=1)
    return cmd_score(a)


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--jobs", type=int, default=8)
    sub.add_parser("score")
    rs = sub.add_parser("resume", help="continue unconverged chains, safeguarded")
    rs.add_argument("--jobs", type=int, default=8)
    rs.add_argument("--extra", type=int, default=EXTRA_STEPS)
    a = ap.parse_args()
    return {"run": cmd_run, "score": cmd_score, "resume": cmd_resume}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
