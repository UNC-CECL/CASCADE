"""
Solve the GIS 1 and GIS 90 end rates of one period against its DEM-to-DEM net-change target.

    python scripts/real-world-cases/Hatteras/input_prep/7-source-sink/2-calibrate/be_end_solve_net_change.py base --period 1996 [--overwrite]
    python scripts/real-world-cases/Hatteras/input_prep/7-source-sink/2-calibrate/be_end_solve_net_change.py solve --period 1996
    python scripts/real-world-cases/Hatteras/input_prep/7-source-sink/2-calibrate/be_end_solve_net_change.py final --period 1996 --ends 14.6 39.9
    python scripts/real-world-cases/Hatteras/input_prep/7-source-sink/2-calibrate/be_end_solve_net_change.py base --period 1996 --hs 2.5

The secant starts from the period's zeroBE matrix run (both ends at 0) and steps
through be_edge_domain_solve --target net_change until both ends are within
0.02 m/yr times the run length. Probes are experiment runs under
end-domain-boundaries-<window>/<date>-ends-solved-on-net-change-<window>/; only the final
edgeBE run goes to the matrix. Full management, no groin, relocations off.
With --hs, every run uses that wave height and the whole solve, base and final
included, is an experiment under <date>-<name>-hs<Hs>-<window>/: the matrix
refuses off-default waves.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-12
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
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in _HERE.parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
sys.path.insert(0, str(_HAT_DIR))

from site_layer.hatteras_site_config import HATTERAS_PERIODS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
HINDCAST = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"
SOLVER_DIR = _HAT_DIR / "input_prep" / "7-source-sink" / "2-calibrate"
SCENARIO = "full_management"
# The solve folder and log date; 2026-10-05 solved the first targets, 2026-10-09 the calendar-2008 end and 17-yr test,
# 2026-10-12 the +/-6-month offset (v3) and start windows and the 2025 +/-6-month test end;
# SOLVE_NAME keeps a second solve on one day apart (2026-10-12 evening: the calendar-year setup);
# 2026-10-13 the dune settings (design height 3.94 m MHW, one island-wide Dmaxel 5.5 m NAVD88);
# 2026-10-20 the v3_240 storms, beach slope 0.06, Dmaxel 12.2 m with growth 0.01-0.03, and the
# 1996 start with the ALACE offset corrected (dune-topo v3); 2026-10-21 the 2009-2025 ends for that window's
# own growth rate, 0.02-0.04 (runner DUNE_GROWTH_RATES; stopped after step 1 for the buffer fix), then both
# windows for the quintic 30-buffer ring closure (hindcast.close_offset_ring) and Cape Point growth
# 0.10-0.20 at GIS 1-5 (runner CAPE_GROWTH_RATE); the same day re-solved with the cape's own rate dropped
# (one dune growth rate island-wide); 2026-10-26 the 1996-2009 ends on CASCADE-obx (the obx-cascade
# branch), at Hs 2.0 and, with --hs 2.5, at Hs 2.5 to test whether the sweep's Hs gain survives its own ends
SOLVE_DATE = "2026-10-26"
SOLVE_NAME = "ends-solved-on-obx"
MAX_STEPS = 8
LOG_ROOT = PROJECT_ROOT / "output" / "logs" / "driver" / "dem_to_dem"
EXP_ROOT = PROJECT_ROOT / "output" / "raw_runs" / "experiments"
# The wave height every run uses; None is the calibration default (set from --hs)
HS = None
# -----------------------------------------------------------------------------


def window(period):
    return f"{period}_{HATTERAS_PERIODS[period]['end_year']}"


def solve_name():
    return SOLVE_NAME + (f"-hs{HS:g}".replace(".", "p") if HS is not None else "")


def solve_tag(period):
    # one topic per window since 2026-10-20 (experiments/RENAMES.md)
    topic = "end-domain-boundaries-" + window(period).replace("_", "-")
    return f"{topic}/{SOLVE_DATE}-{solve_name()}-{window(period)}"


# The zeroBE base run every solve starts from, as the matrix files it
def base_run(period):
    fill = "_nourish" if HATTERAS_PERIODS[period]["enable_nourishment"] else ""
    return f"HAT_{window(period)}_zeroBE_offsetmetres_road_bdm{fill}_nogroin"


# One runner process; matrix unless a tag is given (or --hs, which keeps base and final out of the matrix)
def run(period, preset, label, override=None, tag=None, overwrite=False):
    env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUNBUFFERED="1", MPLBACKEND="Agg",
               HAT_IGNORE_SETTINGS="1", HAT_START_YEAR=str(period),
               HAT_SOURCE_SINK_PRESET=preset, HAT_SCENARIO=SCENARIO, HAT_RELOCATIONS="False",
               HAT_GROIN_ENABLED="False", HAT_OVERWRITE=str(overwrite), HAT_MAKE_GIFS="False",
               HAT_SHOW_FIGURES="False")
    if HS is not None:
        env["HAT_HS"] = repr(HS)
        tag = tag or f"{solve_tag(period)}/runs/{label}"
    if tag:
        env.update(HAT_RUN_KIND="experiment", HAT_RUN_TAG=tag, HAT_SAVE_MODEL_STATE="False")
    else:
        env.update(HAT_RUN_KIND="matrix")
    if override:
        env["HAT_BE_OVERRIDE"] = override
    log = LOG_ROOT / f"end_solve_{SOLVE_DATE}_{solve_name()}_{window(period)}_{label}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    print(f"start {label} {override or ''} -> {log.relative_to(PROJECT_ROOT)}", flush=True)
    with open(log, "w", encoding="utf-8") as f:
        code = subprocess.run([sys.executable, str(HINDCAST)], stdout=f,
                              stderr=subprocess.STDOUT, env=env, cwd=PROJECT_ROOT).returncode
    if code:
        raise SystemExit(f"run {label} failed ({code}); see {log}")


# The run name a tag holds, read from the runs' metadata
def tagged_run(tag):
    from cascade_pipeline.run_registry import scan_runs
    idx = scan_runs(PROJECT_ROOT / "output" / "raw_runs")
    rows = idx[(idx["kind"] == "experiment") & (idx["tag"] == tag)]
    if rows.empty:
        raise SystemExit(f"no run with tag {tag} under output/raw_runs")
    return str(rows.iloc[-1]["run_name"])


# The solver's next step; None once both ends close
def solve_step(period, steps):
    sys.path.insert(0, str(SOLVER_DIR))
    import be_edge_domain_solve as solver
    if HS is None:
        base, base_kind, base_tag = base_run(period), "matrix", ""
    else:
        base_tag = f"{solve_tag(period)}/runs/base"
        base, base_kind = tagged_run(base_tag), "experiment"
    names = [base] + [tagged_run(t) for t in steps]
    kinds = [base_kind] + ["experiment"] * len(steps)
    buf = io.StringIO()
    with redirect_stdout(buf):
        nxt = solver.report(period, names, None, kinds, [base_tag] + steps, target_source="net_change")
    text = buf.getvalue()
    print(text, flush=True)
    out = EXP_ROOT / solve_tag(period)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "solve.txt", "a", encoding="utf-8") as f:
        f.write(text)
    return None if text.count("CONVERGED") >= 2 else nxt


def cmd_solve(a):
    base = solve_tag(a.period)
    steps = []
    while (EXP_ROOT / base / "runs" / f"step{len(steps) + 1}").is_dir():
        steps.append(f"{base}/runs/step{len(steps) + 1}")
    rates = {}
    while True:
        nxt = solve_step(a.period, steps)
        if not nxt or len(steps) >= MAX_STEPS:
            print("solved" if not nxt else "stopped at MAX_STEPS", flush=True)
            break
        # An end the solver has no new step for keeps its last rate; left out of the
        # override, it would fall back to the preset's value (2026-10-20, test step 8)
        rates.update(nxt)
        override = ",".join(f"{g}={r:.4f}" for g, r in sorted(rates.items()))
        tag = f"{base}/runs/step{len(steps) + 1}"
        run(a.period, "edgeBE", f"step{len(steps) + 1}", override, tag=tag)
        steps.append(tag)
    print(f"last step: {steps[-1] if steps else 'base run'}")


# The zeroBE matrix run the solve starts from (both ends at 0)
# --overwrite replaces a base made under earlier settings (the matrix run is the current model's)
def cmd_base(a):
    run(a.period, "zeroBE", "base", overwrite=a.overwrite)


def cmd_final(a):
    run(a.period, "edgeBE", "final", f"1={a.ends[0]},90={a.ends[1]}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 2)[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("base", "solve", "final"):
        p = sub.add_parser(name)
        p.add_argument("--period", type=int, required=True, choices=sorted(HATTERAS_PERIODS))
        if name == "base":
            p.add_argument("--overwrite", action="store_true")
        if name == "final":
            p.add_argument("--ends", type=float, nargs=2, metavar=("GIS1", "GIS90"), required=True)
        p.add_argument("--hs", type=float, default=None,
                       help="wave height for every run; files the whole solve as an experiment")
    a = ap.parse_args(argv)
    global HS
    HS = a.hs
    {"base": cmd_base, "solve": cmd_solve, "final": cmd_final}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
