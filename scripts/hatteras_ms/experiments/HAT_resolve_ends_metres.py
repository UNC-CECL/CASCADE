"""Re-solve the two end domains under the metres offset (2026-09-27).

Hannah, 2026-09-26/27: sweep the waves with the end source/sink terms FIXED,
"but we might need to resolve the ends now that we fixed the offset". The
stored values (HATTERAS_BE_EDGE_ONLY: 1996 +32.2 / +10.0, 2010 +72.6 / +31.3
m/yr) were solved at Hs 2.5 on the /10 offset. Chosen with Hannah:

    reference  the step-2 baseline wave climate, Hs 1.0 m, Tp 8 s,
               asymmetry 0.8, high-angle 0.45: a neutral setting, not the
               winner of either search, so the fixed ends do not pre-favour
               the sweep that uses them
    scenario   full management (road, beach and dune management, fills; no
               relocations, no groin), as the matrix end values always were;
               the pair is then used for natural and managed runs alike
    target     each window's CoastSat LRR: GIS 1 against the raw domain mean,
               GIS 90 against the LOWESS value (10 domains until 2026-09-28,
               7 since: common.SMOOTH_DOMAINS)
    solve      step 0 a fresh zeroBE run; then HAT_wave_shortlist_ends_solved's
               safeguarded step (secant capped at +-30 m/yr until a probe lies
               on each side of the target, then interpolation held inside the
               bracket), first step from the metres response measured on
               2026-09-26 (about 0.24 m/yr of residual per m/yr imposed at
               GIS 1, 0.20 at GIS 90); converged at |residual| <= 0.02 m/yr,
               at most MAX_STEPS probes
    output     tables/ends.json -- {period: {"1": rate, "90": rate}} -- read by
               HAT_wave_grid_fixed_ends.py. The config (HATTERAS_BE_EDGE_ONLY)
               is NOT changed.

WHERE: output/raw_runs/experiments/end-domain-boundaries/2026-09-27-ends-resolved-metres-offset/

    python scripts/hatteras_ms/experiments/HAT_resolve_ends_metres.py

    --tag   file the solve under another study (2026-09-28: the LOWESS-7 re-solve)
    --seed  "1996=4.8394,17.545;2010=18.8,24.535": step 1 probes these ends
            instead of the first-gain guess from zeroBE, so a re-solve near a
            known answer starts there

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_wave_shortlist_ends_solved as E  # noqa: E402

grid, common, step2 = E.grid, E.common, E.step2
TAG = "end-domain-boundaries/2026-09-27-ends-resolved-metres-offset"
STUDY_DIR = grid.RAW_RUNS / "experiments" / TAG
REFERENCE = {"hs": 1.0, "wave_period_s": 8.0, "wave_asymmetry": 0.8,
             "wave_angle_high_fraction": 0.45}
SCENARIO = "full_management"
MAX_STEPS = 8
E.FIRST_GAIN = {1: 0.24, 90: 0.20}          # measured under metres, 2026-09-26

# point the shared launch/find helpers at this study
E.TAG, E.STUDY_DIR = TAG, STUDY_DIR
E.TABLES_DIR, E.LOGS_DIR = STUDY_DIR / "tables", STUDY_DIR / "logs"


def zerobe_run(period):
    """Step 0: the reference setting on zeroBE, run here on the fixed Barrier3D."""
    import subprocess
    s = pd.Series(REFERENCE)
    log = E.LOGS_DIR / SCENARIO / "step0" / f"{period}_{E.label(s)}.log"
    if not grid.finished(log):
        env = grid.run_env("x", SCENARIO, period, REFERENCE)
        env["HAT_RUN_TAG"] = f"{TAG}/runs/{SCENARIO}_step0"
        log.parent.mkdir(parents=True, exist_ok=True)
        p = subprocess.run([sys.executable, str(grid.HINDCAST)], env=env, cwd=str(grid.PROJECT_ROOT),
                           capture_output=True, text=True, encoding="utf-8", errors="replace",
                           timeout=grid.RUN_TIMEOUT_S)
        log.write_text((p.stdout or "") + "\n--- STDERR ---\n" + (p.stderr or ""), encoding="utf-8")
        print(f"step0 {period}: exit {p.returncode}", flush=True)
    # Matched on the wave settings (fixed 2026-09-27): the folder holds the
    # step-0 run of every reference solved so far, and taking the first one
    # found gave the Hs-2 and Hs-2.5 solves the Hs-1 run's residuals as their
    # step 0. The values those solves adopted were each confirmed by a probe
    # run at the right waves, so only the solver's path was affected.
    d = STUDY_DIR / "runs" / f"{SCENARIO}_step0" / grid.window(period) / "zeroBE"
    want = [REFERENCE[k] for k in E.KEYS]
    for md in d.glob("*/*_run_metadata.json"):
        w = json.loads(md.read_text(encoding="utf-8"))["wave climate"]
        got = [float(w["wave_height_m"]), float(w["wave_period_s"]), float(w["wave_asymmetry"]),
               float(w["wave_angle_high_frac"])]
        if all(abs(g - v) < 1e-9 for g, v in zip(got, want)):
            return md.parent
    raise SystemExit(f"step 0 for {period} left no run: {step2.stop_reason(log)}")


def main():
    # Re-solve one window at other waves (2026-09-27, Hannah: "re-solve the
    # 2010 ends at Hs 2"): --periods and the wave settings; the result is
    # merged into ends.json, the replaced value kept under "history". The
    # closest probe (smallest worst-end residual) is the answer, and --accept
    # is the residual the caller will act on (the 0.02 target is not reachable
    # at GIS 1 in 2010-2024, whose response is not monotonic below ~0.1).
    import argparse
    global REFERENCE
    ap = argparse.ArgumentParser()
    ap.add_argument("--periods", nargs="+", type=int, default=list(grid.PERIODS))
    ap.add_argument("--hs", type=float, default=REFERENCE["hs"])
    ap.add_argument("--tp", type=float, default=REFERENCE["wave_period_s"])
    ap.add_argument("--asym", type=float, default=REFERENCE["wave_asymmetry"])
    ap.add_argument("--ahf", type=float, default=REFERENCE["wave_angle_high_fraction"])
    ap.add_argument("--accept", type=float, default=E.TOL)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--seed", default=None)
    a = ap.parse_args()
    global TAG, STUDY_DIR
    if a.tag:
        TAG, STUDY_DIR = a.tag, grid.RAW_RUNS / "experiments" / a.tag
        E.TAG, E.STUDY_DIR = TAG, STUDY_DIR
        E.TABLES_DIR, E.LOGS_DIR = STUDY_DIR / "tables", STUDY_DIR / "logs"
    seed = {}
    for part in (a.seed or "").split(";"):
        if part.strip():
            yr, vals = part.split("=")
            g1, g90 = (float(v) for v in vals.split(","))
            seed[int(yr)] = {1: g1, 90: g90}
    REFERENCE = {"hs": a.hs, "wave_period_s": a.tp, "wave_asymmetry": a.asym,
                 "wave_angle_high_fraction": a.ahf}
    periods = tuple(a.periods)
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    grid.check_barrier3d()
    common.keep_awake()
    E.TABLES_DIR.mkdir(parents=True, exist_ok=True)
    targets = {p: common.coastsat_target(p) for p in periods}
    s = pd.Series(REFERENCE)
    print(f"solving {periods} at {E.label(s)}", flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:
        d0 = dict(zip(periods, pool.map(zerobe_run, periods)))
    chains = {p: [({1: 0.0, 90: 0.0}, E.residuals(d0[p], p, targets), d0[p])] for p in periods}
    rows = []
    for p, h in chains.items():
        print(f"{p}: zeroBE residuals GIS1 {h[0][1][1]:+.3f}, GIS90 {h[0][1][90]:+.3f}", flush=True)
    for step in range(1, MAX_STEPS + 1):
        jobs = []
        for p, h in chains.items():
            if all(abs(h[-1][1][g]) <= E.TOL for g in E.ENDS):
                continue
            hist = [(e, r) for e, r, _ in h]
            if len(h) == 1 and p in seed:
                nxt = dict(seed[p])
            elif len(h) == 1:
                nxt = {g: h[-1][0][g] - h[-1][1][g] / E.FIRST_GAIN[g] for g in E.ENDS}
            else:
                nxt = {g: E.safeguarded_next(hist, g) for g in E.ENDS}
            jobs.append((SCENARIO, step, p, s, nxt))
        if not jobs:
            break
        with ThreadPoolExecutor(max_workers=2) as pool:
            list(pool.map(E.launch, jobs))
        for job in jobs:
            p, nxt = job[2], job[4]
            d = E.find_run(SCENARIO, step, p, s)
            if d is None:
                print(f"step {step} {p}: no run ({step2.stop_reason(E.log_path(SCENARIO, step, p, s))})",
                      flush=True)
                continue
            r = E.residuals(d, p, targets)
            chains[p].append((nxt, r, d))
            rows.append({"period_start": p, "step": step, "gis1_imposed": nxt[1], "gis90_imposed": nxt[90],
                         "gis1_residual": r[1], "gis90_residual": r[90]})
            print(f"step {step} {p}: ends {nxt[1]:+.3f} / {nxt[90]:+.3f} -> residuals "
                  f"{r[1]:+.3f} / {r[90]:+.3f}", flush=True)
        pd.DataFrame(rows).to_csv(E.TABLES_DIR / f"solve_log_{'_'.join(map(str, periods))}_{E.label(s)}.csv",
                                  index=False)
    ends, report = {}, []
    for p, h in chains.items():
        e, r, d = min(h[1:] or h, key=lambda x: max(abs(x[1][g]) for g in E.ENDS))
        ok = all(abs(r[g]) <= E.TOL for g in E.ENDS)
        ends[str(p)] = {"1": round(e[1], 4), "90": round(e[90], 4)}
        report.append({"period": grid.window(p), "gis1": e[1], "gis90": e[90], "gis1_residual": r[1],
                       "gis90_residual": r[90], "converged": ok, "steps": len(h) - 1,
                       "waves": E.label(s),
                       "run_dir": str(Path(d).relative_to(STUDY_DIR)).replace("\\", "/")})
    f = E.TABLES_DIR / "ends.json"
    old = json.loads(f.read_text(encoding="utf-8")) if f.is_file() else {}
    new = dict(old) if old else {"scenario": SCENARIO, "target": "CoastSat LRR per window",
                                 "ends_m_yr": {}}
    new.setdefault("history", [])
    new.setdefault("reference_waves_by_period", {})
    for p, v in ends.items():
        if p in new["ends_m_yr"]:
            new["history"].append({"period": p, "replaced": new["ends_m_yr"][p],
                                   "waves": new["reference_waves_by_period"].get(p, old.get("reference_waves")),
                                   "note": old.get("accepted", {}).get(p, "")})
        new["ends_m_yr"][p] = v
        new["reference_waves_by_period"][p] = dict(REFERENCE)
        new.setdefault("accepted", {}).pop(p, None)
        rep = next(x for x in report if x["period"].startswith(p))
        if not rep["converged"]:
            new["accepted"][p] = (f"closest probe, residuals GIS1 {rep['gis1_residual']:+.3f} / GIS90 "
                                  f"{rep['gis90_residual']:+.3f} m/yr (target 0.02; accepted at <= {a.accept})")
    f.write_text(json.dumps(new, indent=1), encoding="utf-8")
    pd.DataFrame(report).to_csv(E.TABLES_DIR / f"ends_{'_'.join(map(str, periods))}_{E.label(s)}.csv", index=False)
    print(pd.DataFrame(report).round(3).to_string(index=False))
    if not all(max(abs(x["gis1_residual"]), abs(x["gis90_residual"])) <= a.accept for x in report):
        raise SystemExit(f"not within {a.accept} m/yr: the fixed-ends sweep must not start on these values")
    return 0


if __name__ == "__main__":
    sys.exit(main())
