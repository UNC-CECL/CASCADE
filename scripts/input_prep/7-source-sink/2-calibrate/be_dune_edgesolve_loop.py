"""
Drive the dune-line end-domain solve to convergence: ask for the next probe, run it, repeat.

    python scripts/input_prep/7-source-sink/2-calibrate/be_dune_edgesolve_loop.py --exp <experiment>

For each (window, reading) chain until both ends converge; runs land
under the experiment's folder. Details: scripts/input_prep/7-source-sink/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
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
# --- CONFIG ------------------------------------------------------------------
SOLVER = Path(__file__).with_name("be_edge_domain_solve.py")
RUNNER = PROJECT_ROOT / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
RUN_ROOT = PROJECT_ROOT / "output" / "raw_runs"
BRACKET_EXP = "end-domain-boundaries/2026-09-16-end-domains-solved-on-duneline"

# Window labels from the config (1996_2015, 2010_2026 since 2026-10-02); were hardcoded 2010/2024 ends
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from site_layer.hatteras_site_config import HATTERAS_PERIODS  # noqa: E402
END = {p: v["end_year"] for p, v in HATTERAS_PERIODS.items()}
# Every window now carries a fill, so the full-management run is the _nourish one
SUFFIX = {1996: "road_bdm_nourish", 2004: "road_bdm_nourish", 2010: "road_bdm_nourish"}


# The matrix run names carry the offset token since the metres offset became the default (2026-09-24
OFFSET_TOKEN = "offsetmetres"
# -----------------------------------------------------------------------------


# [(run_name, kind, tag), ...] for zeroBE then edgeBE, full management
def brackets(start):
    name = lambda p: f"HAT_{start}_{END[start]}_{p}_{OFFSET_TOKEN}_{SUFFIX[start]}_nogroin"  # noqa: E731
    if start == 2004:
        tag = f"{BRACKET_EXP}/brackets"
        return [(name("zeroBE"), "experiment", tag), (name("edgeBE"), "experiment", tag)]
    return [(name("zeroBE"), "matrix", ""), (name("edgeBE"), "matrix", "")]


COASTSAT_WINDOW = None   # set by --coastsat-window


# Run the solver over `runs`
def solve(start, smooth, runs):
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


# {gis
def imposed(run_dir):
    import json
    meta = next(Path(run_dir).glob("*_run_metadata.json"))
    text = meta.read_text(encoding="utf-8")
    out = {}
    for gis, key in ((1, "be_rate_gis1_m_yr"), (90, "be_rate_gis90_m_yr")):
        m = re.search(rf'"{key}": *(-?[0-9.eE+-]+)', text)
        out[gis] = float(m.group(1)) if m else None
    return out


# The solver prints only the ends it wants MOVED
def merge_override(override, last):
    pairs = dict(p.split("=") for p in (override or "").split(",") if p)
    for gis, val in last.items():
        if str(gis) not in pairs and val is not None:
            pairs[str(gis)] = f"{val:g}"
    return ",".join(f"{g}={pairs[g]}" for g in sorted(pairs, key=int))


# Every finished step already on disk as (run_name, 'experiment', tag, run_dir), oldest first
def existing_steps(exp, start, smooth):
    out, k = [], 1
    while True:
        tag = f"{exp}/{smooth}/step{k}"
        base = RUN_ROOT / "experiments" / exp / smooth / f"step{k}" / f"{start}_{END[start]}" / "edgeBE"
        runs = [d for d in base.glob("HAT_*") if list(d.glob("*_run_metadata.json"))] if base.is_dir() else []
        if not runs:
            return out
        out.append((runs[0].name, "experiment", tag, runs[0]))
        k += 1


# Run one probe of the hindcast with the solved end rates
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


# Run: every chain until both ends converge
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
    ap.add_argument("--hold", type=int, nargs="*", default=[],
                    help="end domains kept at their config value and left out of convergence (e.g. --hold 1)")
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
            if k >= 1 and (1 in a.hold or abs(r1) < a.tol) and (90 in a.hold or abs(r90) < a.tol):
                done[key] = k
                print(f"          {w} {s}: CONVERGED at step {k}", flush=True)
                continue
            if k >= a.max_steps:
                print(f"          {w} {s}: NOT CONVERGED after {k} steps", flush=True)
                continue
            full = merge_override(override, last[key])
            if a.hold:
                # A held end keeps the config value exactly, whatever the solver printed
                from site_layer.hatteras_site_config import HATTERAS_BE_RATES_EDGE
                pairs = dict(p.split("=") for p in full.split(","))
                for g in a.hold:
                    pairs[str(g)] = f"{HATTERAS_BE_RATES_EDGE[w][g]:g}"
                full = ",".join(f"{g}={pairs[g]}" for g in sorted(pairs, key=int))
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
