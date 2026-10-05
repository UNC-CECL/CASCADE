"""
Solve the GIS 1 and GIS 90 end rates of one period against its DEM-to-DEM net-change target.

    python scripts/hatteras_ms/HAT_end_solve_net_change.py solve --period 1996
    python scripts/hatteras_ms/HAT_end_solve_net_change.py final --period 1996 --ends 14.6 39.9

The secant starts from the period's zeroBE matrix run (both ends at 0) and steps
through be_edge_domain_solve --target net_change until both ends are within
0.02 m/yr times the run length. Probes are experiment runs under
end-domain-boundaries/<date>-ends-solved-on-net-change-<window>/; only the final
edgeBE run goes to the matrix. Full management, no groin, relocations off.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-05
"""
from __future__ import annotations

import argparse
import io
import os
import subprocess
import sys
from contextlib import redirect_stdout
from pathlib import Path

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from site_layer.hatteras_site_config import HATTERAS_PERIODS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
HINDCAST = _HERE.parent / "HAT_hindcast_1984_2024.py"
SOLVER_DIR = PROJECT_ROOT / "scripts" / "input_prep" / "7-source-sink" / "2-calibrate"
SCENARIO = "full_management"
SOLVE_DATE = "2026-10-05"
MAX_STEPS = 8
LOG_ROOT = PROJECT_ROOT / "output" / "logs" / "driver" / "dem_to_dem"
EXP_ROOT = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / "end-domain-boundaries"
# -----------------------------------------------------------------------------


def window(period):
    return f"{period}_{HATTERAS_PERIODS[period]['end_year']}"


def solve_tag(period):
    return f"end-domain-boundaries/{SOLVE_DATE}-ends-solved-on-net-change-{window(period)}"


# The zeroBE base run every solve starts from, as the matrix files it
def base_run(period):
    fill = "_nourish" if HATTERAS_PERIODS[period]["enable_nourishment"] else ""
    return f"HAT_{window(period)}_zeroBE_offsetmetres_road_bdm{fill}_nogroin"


# One runner process; matrix unless a tag is given
def run(period, preset, label, override=None, tag=None):
    env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUNBUFFERED="1", MPLBACKEND="Agg",
               HAT_IGNORE_SETTINGS="1", HAT_START_YEAR=str(period),
               HAT_SOURCE_SINK_PRESET=preset, HAT_SCENARIO=SCENARIO, HAT_RELOCATIONS="False",
               HAT_GROIN_ENABLED="False", HAT_OVERWRITE="False", HAT_MAKE_GIFS="False",
               HAT_SHOW_FIGURES="False")
    if tag:
        env.update(HAT_RUN_KIND="experiment", HAT_RUN_TAG=tag, HAT_SAVE_MODEL_STATE="False")
    else:
        env.update(HAT_RUN_KIND="matrix")
    if override:
        env["HAT_BE_OVERRIDE"] = override
    log = LOG_ROOT / f"end_solve_{window(period)}_{label}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    print(f"start {label} {override or ''} -> {log.relative_to(PROJECT_ROOT)}", flush=True)
    with open(log, "w", encoding="utf-8") as f:
        code = subprocess.run([sys.executable, str(HINDCAST)], stdout=f,
                              stderr=subprocess.STDOUT, env=env, cwd=PROJECT_ROOT).returncode
    if code:
        raise SystemExit(f"run {label} failed ({code}); see {log}")


# The run name a tag holds, read from the index
def tagged_run(tag):
    from cascade_pipeline.run_registry import load_run_index
    idx = load_run_index(PROJECT_ROOT / "output" / "raw_runs" / "run_index.csv")
    rows = idx[(idx["kind"] == "experiment") & (idx["tag"] == tag)]
    if rows.empty:
        raise SystemExit(f"no run with tag {tag} in run_index.csv")
    return str(rows.iloc[-1]["run_name"])


# The solver's next step; None once both ends close
def solve_step(period, steps):
    sys.path.insert(0, str(SOLVER_DIR))
    import be_edge_domain_solve as solver
    names = [base_run(period)] + [tagged_run(t) for t in steps]
    kinds = ["matrix"] + ["experiment"] * len(steps)
    buf = io.StringIO()
    with redirect_stdout(buf):
        nxt = solver.report(period, names, None, kinds, [""] + steps, target_source="net_change")
    text = buf.getvalue()
    print(text, flush=True)
    out = EXP_ROOT / solve_tag(period).split("/", 1)[1]
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "solve.txt", "a", encoding="utf-8") as f:
        f.write(text)
    return None if text.count("CONVERGED") >= 2 else nxt


def cmd_solve(a):
    base = solve_tag(a.period)
    steps = []
    while (EXP_ROOT / base.split("/", 1)[1] / "runs" / f"step{len(steps) + 1}").is_dir():
        steps.append(f"{base}/runs/step{len(steps) + 1}")
    while True:
        nxt = solve_step(a.period, steps)
        if not nxt or len(steps) >= MAX_STEPS:
            print("solved" if not nxt else "stopped at MAX_STEPS", flush=True)
            break
        override = ",".join(f"{g}={r:.4f}" for g, r in sorted(nxt.items()))
        tag = f"{base}/runs/step{len(steps) + 1}"
        run(a.period, "edgeBE", f"step{len(steps) + 1}", override, tag=tag)
        steps.append(tag)
    print(f"last step: {steps[-1] if steps else 'base run'}")


def cmd_final(a):
    run(a.period, "edgeBE", "final", f"1={a.ends[0]},90={a.ends[1]}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 2)[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("solve", "final"):
        p = sub.add_parser(name)
        p.add_argument("--period", type=int, required=True, choices=sorted(HATTERAS_PERIODS))
        if name == "final":
            p.add_argument("--ends", type=float, nargs=2, metavar=("GIS1", "GIS90"), required=True)
    a = ap.parse_args(argv)
    {"solve": cmd_solve, "final": cmd_final}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
