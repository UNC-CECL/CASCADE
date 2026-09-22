"""
be_dune_edgesolve_loop.py
==============================================================================
Drive the dune-line end-domain solve to convergence: for each (window,
reading) chain, ask be_edge_domain_solve.py for the next probe, run it,
and repeat until both ends sit within --tol of the dune-line target. Written
2026-09-18 for the re-solve after the 1997, 2009 and 2023 dune lines were
re-digitized; the 09-16 solve was stepped by hand.

WHAT IT DOES NOT CHANGE
    The arithmetic is be_edge_domain_solve.py's (target from the stored
    5-scr/3-rates/duneline/endpoint product, the local secant through the last
    two runs, --estimator endpoint). This only runs the probes it prints.

LOCKSTEP
    Every live chain runs its next probe AT THE SAME TIME (one runner process
    each), then all are waited for, then every solver is read. The solver
    reads each run's imposed rates from run_index.csv, so reading only after a
    whole step has finished keeps it from reading an index a still-finishing
    run is rewriting.

BRACKETS (step 0, not re-run)
    1996, 2010   the current matrix zeroBE and edgeBE full-management runs
    2004         the 09-16 brackets (experiments/2026-09-16-dune-edgesolve/
                 brackets): the 2004-start inputs did not change on 09-18

OUTPUT   output/raw_runs/experiments/<exp>/<reading>/step<k>/<window>/edgeBE/<run>/
         output/raw_runs/experiments/<exp>/logs/<reading>_step<k>_<start>.log
         output/raw_runs/experiments/<exp>/loop_log.csv   one row per step
                                                          per chain

COASTSAT TARGET (2026-09-19)
    --target coastsat runs the same loop against the CoastSat target, as the
    matrix end values were solved (model lrr_m_yr against target_lrr_m_yr,
    GIS 1 raw, GIS 90 LOESS-10). One chain per window, filed under the
    reading name "coastsat"; --smooth is ignored. E.g. --exp
    2026-09-19-edgesolve-2010 --windows 2010 --target coastsat.

USAGE
    python be_dune_edgesolve_loop.py --exp 2026-09-18-dune-edgesolve \\
        --windows 1996 2004 2010 --smooth raw mean3
==============================================================================
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
SOLVER = Path(__file__).with_name("be_edge_domain_solve.py")
RUNNER = PROJECT_ROOT / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
RUN_ROOT = PROJECT_ROOT / "output" / "raw_runs"
BRACKET_EXP = "2026-09-16-dune-edgesolve"

END = {1996: 2010, 2004: 2024, 2010: 2024}
SUFFIX = {1996: "road_bdm", 2004: "road_bdm_nourish", 2010: "road_bdm_nourish"}


def brackets(start):
    """[(run_name, kind, tag), ...] for zeroBE then edgeBE, full management."""
    name = lambda p: f"HAT_{start}_{END[start]}_{p}_{SUFFIX[start]}_nogroin"  # noqa: E731
    if start == 2004:
        tag = f"{BRACKET_EXP}/brackets"
        return [(name("zeroBE"), "experiment", tag), (name("edgeBE"), "experiment", tag)]
    return [(name("zeroBE"), "matrix", ""), (name("edgeBE"), "matrix", "")]


COASTSAT_WINDOW = None   # set by --coastsat-window


def solve(start, smooth, runs):
    """Run the solver over `runs`; return (residual GIS 1, residual GIS 90,
    override string, full text). smooth == "coastsat" is the CoastSat target."""
    if smooth == "coastsat":
        cmd = [sys.executable, str(SOLVER), "--period", str(start),
               "--target", "coastsat", "--estimator", "lrr"]
        if COASTSAT_WINDOW:
            cmd += ["--coastsat-window", COASTSAT_WINDOW]
    else:
        cmd = [sys.executable, str(SOLVER), "--period", str(start),
               "--target", "duneline", "--dune-smooth", smooth, "--estimator", "endpoint"]
    for name, kind, tag in runs:
        cmd += ["--run", name, "--kind", kind, "--tag", tag]
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    out = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                         env=env, cwd=PROJECT_ROOT)
    text = out.stdout + out.stderr
    if out.returncode != 0:
        raise RuntimeError(f"solver failed for {start} {smooth}:\n{text}")
    resid = []
    for block in re.split(r"\n(?=GIS \d+)", text)[1:3]:
        rows = [l for l in block.splitlines() if l.strip().startswith("HAT_")
                and not l.strip().startswith("HAT_BE_OVERRIDE")]
        resid.append(float(rows[-1].split()[-1]))
    m = re.search(r'HAT_BE_OVERRIDE="([^"]+)"', text)
    return resid[0], resid[1], (m.group(1) if m else None), text


def imposed(run_dir):
    """{gis: rate} a finished run imposed at the two ends, from its log line
    'GIS <n>  <preset> -> <imposed> m/yr' or, when not overridden, the preset."""
    import json
    meta = next(Path(run_dir).glob("*_run_metadata.json"))
    text = meta.read_text(encoding="utf-8")
    out = {}
    for gis, key in ((1, "be_rate_gis1_m_yr"), (90, "be_rate_gis90_m_yr")):
        m = re.search(rf'"{key}": *(-?[0-9.eE+-]+)', text)
        out[gis] = float(m.group(1)) if m else None
    return out


def merge_override(override, last):
    """The solver prints only the ends it wants MOVED; an end it leaves out
    keeps the value the chain's last run imposed. Passing the solver's string
    through as-is (the first version of this driver, 2026-09-18) reset such an
    end to the edgeBE PRESET -- +32.2 at GIS 1 for 1996 -- for one probe."""
    pairs = dict(p.split("=") for p in (override or "").split(",") if p)
    for gis, val in last.items():
        if str(gis) not in pairs and val is not None:
            pairs[str(gis)] = f"{val:g}"
    return ",".join(f"{g}={pairs[g]}" for g in sorted(pairs, key=int))


def existing_steps(exp, start, smooth):
    """[(run_name, 'experiment', tag, run_dir), ...] for every FINISHED step
    already on disk, oldest first, stopping at the first gap."""
    out, k = [], 1
    while True:
        tag = f"{exp}/{smooth}/step{k}"
        base = RUN_ROOT / "experiments" / exp / smooth / f"step{k}" / f"{start}_{END[start]}" / "edgeBE"
        runs = [d for d in base.glob("HAT_*") if list(d.glob("*_run_metadata.json"))] if base.is_dir() else []
        if not runs:
            return out
        out.append((runs[0].name, "experiment", tag, runs[0]))
        k += 1


def launch(start, smooth, step, override, exp):
    tag = f"{exp}/{smooth}/step{step}"
    env = dict(os.environ,
               PYTHONIOENCODING="utf-8", HAT_IGNORE_SETTINGS="1",
               HAT_START_YEAR=str(start), HAT_SOURCE_SINK_PRESET="edgeBE",
               HAT_SCENARIO="full_management", HAT_RELOCATIONS="False",
               HAT_GROIN_ENABLED="False", HAT_RUN_KIND="experiment",
               HAT_RUN_TAG=tag, HAT_SAVE_MODEL_STATE="False",
               HAT_OVERWRITE="False", HAT_BE_OVERRIDE=override)
    logs = RUN_ROOT / "experiments" / exp / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    log = open(logs / f"{smooth}_step{step}_{start}.log", "w", encoding="utf-8")
    proc = subprocess.Popen([sys.executable, str(RUNNER)], stdout=log,
                            stderr=subprocess.STDOUT, env=env, cwd=PROJECT_ROOT)
    return proc, log, tag


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", required=True)
    ap.add_argument("--windows", type=int, nargs="+", default=[1996, 2004, 2010])
    ap.add_argument("--smooth", nargs="+", default=["raw", "mean3"])
    ap.add_argument("--tol", type=float, default=0.01, help="m/yr, both ends")
    ap.add_argument("--max-steps", type=int, default=6)
    ap.add_argument("--target", choices=("duneline", "coastsat"), default="duneline")
    ap.add_argument("--coastsat-window", default=None, metavar="START_END",
                    help="with --target coastsat: solve against this LRR window "
                         "instead of the run's own (e.g. 1996_2024)")
    ap.add_argument("--resume", action="store_true",
                    help="pick every chain up from the finished steps on disk")
    a = ap.parse_args(argv)
    if a.target == "coastsat":
        a.smooth = ["coastsat"]
        global COASTSAT_WINDOW
        COASTSAT_WINDOW = a.coastsat_window

    # chain: [(run_name, kind, tag), ...]; last: {gis: rate} its last run imposed
    chains, last = {}, {}
    for w in a.windows:
        for s in a.smooth:
            chain = brackets(w)
            last[(w, s)] = {1: None, 90: None}
            if a.resume:
                for name, kind, tag, run_dir in existing_steps(a.exp, w, s):
                    chain.append((name, kind, tag))
                    last[(w, s)] = imposed(run_dir)
            chains[(w, s)] = chain
    done = {}
    loop_log = RUN_ROOT / "experiments" / a.exp / "loop_log.csv"
    loop_log.parent.mkdir(parents=True, exist_ok=True)
    if not (a.resume and loop_log.is_file()):
        with open(loop_log, "w", newline="", encoding="utf-8") as fh:
            csv.writer(fh).writerow(["window", "smooth", "step", "run_tag",
                                     "resid_gis1", "resid_gis90", "next_override"])

    def nsteps(key):
        return len(chains[key]) - 2          # the two brackets are step 0

    for _round in range(a.max_steps + 1):
        live = [k for k in chains if k not in done and nsteps(k) <= a.max_steps]
        if not live:
            break
        procs = []
        for w, s in live:
            key = (w, s)
            r1, r90, override, _ = solve(w, s, chains[key])
            k = nsteps(key)
            with open(loop_log, "a", newline="", encoding="utf-8") as fh:
                csv.writer(fh).writerow([w, s, k, chains[key][-1][2], r1, r90, override])
            print(f"step {k}  {w} {s:<5}  residual {r1:+.4f} / {r90:+.4f}"
                  f"   next {override}", flush=True)
            if k >= 1 and abs(r1) < a.tol and abs(r90) < a.tol:
                done[key] = k
                print(f"          {w} {s}: CONVERGED at step {k}", flush=True)
                continue
            if k >= a.max_steps:
                print(f"          {w} {s}: NOT CONVERGED after {k} steps", flush=True)
                continue
            full = merge_override(override, last[key])
            proc, log, tag = launch(w, s, k + 1, full, a.exp)
            procs.append((key, proc, log, tag, full))
        for key, proc, log, tag, full in procs:
            rc = proc.wait()
            log.close()
            if rc != 0:
                raise RuntimeError(f"run failed: {key} {tag}, see its log")
            name = chains[key][-1][0].replace("zeroBE", "edgeBE")
            chains[key].append((name, "experiment", tag))
            last[key] = {int(g): float(v) for g, v in
                         (p.split("=") for p in full.split(","))}
        time.sleep(2)

    print(chr(10) + "solved: " + " ".join(f"{w}:{s}:{k}" for (w, s), k in sorted(done.items())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
