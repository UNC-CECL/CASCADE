"""
One run over the whole record: 1996 start, through December 2025 (30 model years).

    python scripts/hatteras_ms/HAT_full_window_1996_2025.py check
    python scripts/hatteras_ms/HAT_full_window_1996_2025.py zero
    python scripts/hatteras_ms/HAT_full_window_1996_2025.py solve
    python scripts/hatteras_ms/HAT_full_window_1996_2025.py final

Runs the unchanged runner with the 1996 period patched in-process (end 2025,
its own storms, sea-level rate and CoastSat LRR); HATTERAS_PERIODS is not
edited. zero and final file in matrix/1996_2025/; the end-solve probes file as
an experiment. Full management, no groin, no relocation.

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

import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT))

# --- CONFIG ------------------------------------------------------------------
HINDCAST = PROJECT_ROOT / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
SOLVER_DIR = PROJECT_ROOT / "scripts" / "input_prep" / "7-source-sink" / "2-calibrate"
START, END = 1996, 2025
# The window label is inclusive: the run steps 1996..2025 and ends on 1 Jan 2026
LAST_MODEL_YEAR = 2025
SCENARIO = "full_management"
# Solve probes are not production: they file as an experiment, only zero and final go to matrix/
SOLVE_TAG = "end-domain-boundaries/2026-10-05-ends-solved-on-1996_2025"
SOLVE_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / SOLVE_TAG
LOG_DIR = PROJECT_ROOT / "output" / "logs" / "driver" / "full_window_1996_2025"
MAX_STEPS = 6
# -----------------------------------------------------------------------------


# The sea-level rate the config would carry for this window: the 0.001-rounded Duck fit
def rslr_rate():
    from site_layer import hat_env_forcings as _env
    rates = pd.read_csv(_env.RSLR_FITS_DIR / "duck_rslr_rates.csv")
    row = rates[rates["window"] == f"{START}_{END}"]
    if row.empty:
        raise SystemExit(f"no {START}_{END} row in duck_rslr_rates.csv; run duck_rslr_analysis.py")
    return float(row["config_m_yr"].iloc[0])


# Patch the 1996 period to the full window; the runner then reads storms, RSLR and target from it
def patch_period():
    from site_layer import hatteras_site_config as sc
    from site_layer import hat_env_forcings as _env
    period = sc.HATTERAS_PERIODS[START]
    period["end_year"] = END
    period["last_model_year"] = LAST_MODEL_YEAR
    period["storm_file"] = _env.init_relpath(_env.storm_series_file(START, END))
    period["sea_level_rise_rate"] = rslr_rate()
    return period


# Every input the run reads for the window, present and the length it should be
def cmd_check(_args):
    import numpy as np
    from site_layer import hatteras_site_config as sc
    from site_layer import hat_env_forcings as _env
    from site_layer.hat_observed_rates import lrr_csv
    period = patch_period()
    storms = np.load(_env.storm_series_file(START, END))
    print(f"window        {START}-{END}, last model year {sc.last_model_year(START)}, "
          f"{sc.run_years(START)} transitions")
    print(f"storms        {period['storm_file']}: {len(storms)} events, "
          f"model years {int(storms[:, 0].min())}-{int(storms[:, 0].max())}")
    print(f"RSLR          {period['sea_level_rise_rate']} m/yr")
    print(f"target        {lrr_csv(START, END)}")
    print(f"offset        {period['island_offset_file']}")
    print(f"setback       {period['road_setback_file']}")
    print(f"topography    {period['topo_product']}")
    if int(storms[:, 0].max()) != sc.run_years(START):
        raise SystemExit("storm file and run length disagree")


# Subprocess entry: patch the period, then run the unchanged runner
def _launch():
    import runpy
    period = patch_period()
    print(f"full window: {START}-{END}; storms {period['storm_file']}; "
          f"RSLR {period['sea_level_rise_rate']}")
    sys.argv = [str(HINDCAST)]
    runpy.run_path(str(HINDCAST), run_name="__main__")


# One run through _launch in its own process; matrix unless a tag is given
def run(preset, label, override=None, tag=None):
    env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUNBUFFERED="1",
               HAT_IGNORE_SETTINGS="1",
               HAT_START_YEAR=str(START), HAT_SOURCE_SINK_PRESET=preset,
               HAT_SCENARIO=SCENARIO, HAT_RELOCATIONS="False",
               HAT_GROIN_ENABLED="False", HAT_OVERWRITE="False")
    if tag:
        env.update(HAT_RUN_KIND="experiment", HAT_RUN_TAG=tag, HAT_SAVE_MODEL_STATE="False")
    else:
        env.update(HAT_RUN_KIND="matrix")
    if override:
        env["HAT_BE_OVERRIDE"] = override
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"{preset}_{label}.log"
    print(f"start {preset} {label} -> {log_path.relative_to(PROJECT_ROOT)}", flush=True)
    with open(log_path, "w", encoding="utf-8") as log:
        code = subprocess.run([sys.executable, str(_HERE), "_launch"],
                              stdout=log, stderr=subprocess.STDOUT, env=env,
                              cwd=PROJECT_ROOT).returncode
    if code:
        raise SystemExit(f"run {preset} {label} failed ({code}); see {log_path}")
    print(f"done  {preset} {label}", flush=True)


# The run name a tag holds, read from the index
def run_name(tag):
    from cascade_pipeline.run_registry import load_run_index
    idx = load_run_index(PROJECT_ROOT / "output" / "raw_runs" / "run_index.csv")
    rows = idx[(idx["kind"] == "experiment") & (idx["tag"] == tag)]
    if rows.empty:
        raise SystemExit(f"no run with tag {tag} in run_index.csv")
    return str(rows.iloc[-1]["run_name"])


def cmd_zero(_args):
    run("zeroBE", "zero")


# The solver's next step with the period patched; None once both ends close
def solve_step(steps):
    patch_period()
    sys.path.insert(0, str(SOLVER_DIR))
    import be_edge_domain_solve as solver
    names = [run_name(t) for t in steps]
    buf = io.StringIO()
    with redirect_stdout(buf):
        nxt = solver.report(START, names, "edgeBE", ["experiment"] * len(steps),
                            steps, target_source="coastsat", coastsat_window=(START, END))
    text = buf.getvalue()
    print(text)
    SOLVE_DIR.mkdir(parents=True, exist_ok=True)
    with open(SOLVE_DIR / "solve.txt", "a", encoding="utf-8") as f:
        f.write(text)
    # The solver flags each end under 0.02 m/yr; both flagged is solved
    return None if text.count("CONVERGED") >= 2 else nxt


# Secant from the 1996-2015 ends until both ends close
def cmd_solve(args):
    from site_layer.hatteras_site_config import HATTERAS_BE_EDGE_ONLY
    steps = []
    k = 0
    while (SOLVE_DIR / "runs" / f"step{k}").is_dir():
        steps.append(f"{SOLVE_TAG}/runs/step{k}")
        k += 1
    if not steps:
        g1, g90 = args.seed if args.seed else HATTERAS_BE_EDGE_ONLY[START]
        tag = f"{SOLVE_TAG}/runs/step0"
        run("edgeBE", "step0", f"1={g1},90={g90}", tag=tag)
        steps.append(tag)
    while True:
        nxt = solve_step(steps)
        if not nxt or len(steps) > MAX_STEPS:
            print("solved" if not nxt else "stopped at MAX_STEPS")
            break
        override = ",".join(f"{g}={r:.4f}" for g, r in sorted(nxt.items()))
        tag = f"{SOLVE_TAG}/runs/step{len(steps)}"
        run("edgeBE", f"step{len(steps)}", override, tag=tag)
        steps.append(tag)
    print(f"final step: {steps[-1]}")


# The solved ends as a matrix run
def cmd_final(args):
    run("edgeBE", "final", f"1={args.ends[0]},90={args.ends[1]}")


def main(argv=None):
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch()
        return 0
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 2)[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check")
    sub.add_parser("zero")
    s = sub.add_parser("solve")
    s.add_argument("--seed", type=float, nargs=2, metavar=("GIS1", "GIS90"),
                   help="first-step end rates, default the 1996-2015 ones")
    f = sub.add_parser("final")
    f.add_argument("--ends", type=float, nargs=2, metavar=("GIS1", "GIS90"), required=True)
    a = ap.parse_args(argv)
    {"check": cmd_check, "zero": cmd_zero, "solve": cmd_solve, "final": cmd_final}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
