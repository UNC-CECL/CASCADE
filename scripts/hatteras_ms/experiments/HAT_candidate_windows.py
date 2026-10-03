"""
What does the model give on the candidate windows 1996-2015 and 2010-2026?

    python scripts/hatteras_ms/experiments/HAT_candidate_windows.py build
    python scripts/hatteras_ms/experiments/HAT_candidate_windows.py forcing
    python scripts/hatteras_ms/experiments/HAT_candidate_windows.py zero
    python scripts/hatteras_ms/experiments/HAT_candidate_windows.py solve --period 1996
    python scripts/hatteras_ms/experiments/HAT_candidate_windows.py solve --period 2010

Runs the unchanged runner with the window swapped in, in-process: 1996 ends at
2015 on its own storms and sea-level rate, and each period is graded against
its candidate window's CoastSat table. 2010 ends at 2026 on forcing records
extended through 2025 (forcing). Details: output/raw_runs/experiments/candidate-windows/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-02
"""
from __future__ import annotations

import argparse
import io
import json
import os
import subprocess
import sys
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT))

# --- CONFIG ------------------------------------------------------------------
HINDCAST = PROJECT_ROOT / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
BUILDER = (PROJECT_ROOT / "scripts" / "input_prep" / "3-env-forcings" / "3-storms"
           / "historical_storm_creation_v3_HAT.py")
SOLVER_DIR = PROJECT_ROOT / "scripts" / "input_prep" / "7-source-sink" / "2-calibrate"
TAG = "candidate-windows/2026-10-02-candidate-windows"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
INPUTS = EXP_DIR / "inputs"

# Period start -> (model end year, CoastSat window it is graded and solved against)
# 2010 ran to 2024 until Hannah moved it to 2026 (2026-10-02); those runs stay under runs/2010_*
PERIODS = {1996: (2015, (1996, 2015)), 2010: (2026, (2010, 2026))}
# The run-folder prefix of a period: the 2026 runs need their own beside the superseded 2024 ones
RUN_KEY = {1996: "1996", 2010: "2010_2026"}
# The forcing records extended through 2025, beside the 1984-2024 originals
FORCING_END = "2025-12-31 23:00:00"
WIS_ARCHIVE = "https://chldata.erdc.dren.mil/thredds/dodsC/wis/Atlantic/ST63228"
# The adopted storm rules (split12, trim24); the series lands in hindcast_storms/1996_2015/
STORM_ARGS = ["--long-events", "trim", "--max-duration", "24", "--split-gap", "12"]
# The scenario every run here uses: the hindcast baseline, no groin, no relocation
SCENARIO = "full_management"
MAX_STEPS = 6
# -----------------------------------------------------------------------------


# The 1996-2015 sea-level rate, fitted the way duck_rslr_analysis.py fits every window
def fit_rslr(start, end):
    sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "input_prep" / "3-env-forcings"
                           / "2-rslr"))
    import duck_rslr_analysis as dr
    with redirect_stdout(io.StringIO()):
        df = dr.load_noaa_meantrend(dr.DATA_FILE)
    fit = dr.fit_linear_trend(df, start, end)
    return dict(window=f"{start}_{end}", slope_m_yr=float(fit["slope_m_yr"]),
                ci95_m_yr=float(fit["ci95_m_yr"]), n_months=int(fit["n"]),
                config_m_yr=round(float(fit["slope_m_yr"]), dr.CONFIG_DECIMALS))


# Storms and sea-level rate for 1996-2015, into the experiment's inputs/
def cmd_build(_args):
    INPUTS.mkdir(parents=True, exist_ok=True)
    # The builder's own window folder: the runner records storm paths relative to data/hatteras_init
    subprocess.run([sys.executable, str(BUILDER), "--start-year", "1996",
                    "--end-year", "2015", *STORM_ARGS],
                   check=True, cwd=PROJECT_ROOT)
    rslr = fit_rslr(1996, 2015)
    (INPUTS / "rslr_1996_2015.json").write_text(json.dumps(rslr, indent=2))
    print(json.dumps(rslr, indent=2))


# The two extended records: (Duck water levels, WIS waves), both through 2025
def extended_records():
    from site_layer import hat_env_forcings as _env
    end = FORCING_END[:10].replace("-", "")
    duck = _env.WATER_LEVEL_DIR / f"8651370_DUCK_19840101_{end}_NAVD.csv"
    wis = _env.WIS_FILE.parent / f"ST63228_19840101_{end}_export_plus_archive.csv"
    return duck, wis


# Parse an OPeNDAP .ascii response into {name: array}
def _opendap(url):
    import urllib.request
    txt = urllib.request.urlopen(url, timeout=180).read().decode()
    body = txt.split("---------------------------------------------", 1)[1]
    out, name, buf = {}, None, []
    for line in body.strip().split("\n") + [""]:
        line = line.strip()
        if not line or "[" in line or (", " in line and line.split(",")[0].isalpha()):
            if name and buf:
                out[name] = np.array([float(v) for v in ",".join(buf).split(",")
                                      if v.strip()])
            name, buf = None, []
            if "[" in line:
                name = line.split("[")[0]
            elif line:
                k, v = line.split(",", 1)
                out[k.strip()] = np.array([float(v)])
            continue
        buf.append(line)
    return out


# One month of WIS station 63228 from the archive, in the export's columns
def wis_month(y, m):
    a = _opendap(f"{WIS_ARCHIVE}/{y}/WIS-ocean_waves_ST63228_{y}{m:02d}.nc.ascii?"
                 "time,waveHs,waveTp,waveMeanDirection,latitude,longitude")
    t = pd.to_datetime(np.round(a["time"] / 3600) * 3600, unit="s")
    return pd.DataFrame({"time": t.strftime("%Y-%m-%d %H:%M:%S"),
                         "lat": a["latitude"][0], "lon": a["longitude"][0],
                         "waveTp": a["waveTp"], "waveMeanDirection": a["waveMeanDirection"],
                         "waveHs": a["waveHs"]})


# Extend both records through 2025, then build the 2010-2026 storms and sea-level rate
def cmd_forcing(_args):
    import importlib.util
    from site_layer import hat_env_forcings as _env
    duck, wis = extended_records()
    INPUTS.mkdir(parents=True, exist_ok=True)

    # Water levels: the unchanged downloader with its end moved; cached months are not refetched
    spec = importlib.util.spec_from_file_location(
        "dl", PROJECT_ROOT / "scripts" / "input_prep" / "3-env-forcings" / "1-records"
        / "HAT_download_water_levels.py")
    dl = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(dl)
    dl.END = FORCING_END
    dl.main()
    if not duck.exists():
        raise SystemExit(f"the downloader did not write {duck.name}")

    # Waves: the export as it is, plus every 2025 month from the archive
    old = pd.read_csv(_env.WIS_FILE)
    check = wis_month(2024, 12)
    ref = old[old["time"].str.startswith("2024-12")].reset_index(drop=True)
    same = (len(check) == len(ref) and (check["time"] == ref["time"]).all() and all(
        np.allclose(check[c], ref[c]) for c in ("waveHs", "waveTp", "waveMeanDirection")))
    if not same:
        raise SystemExit("the archive's Dec 2024 differs from the export; stopping")
    months = [wis_month(2025, m) for m in range(1, 13)]
    ext = pd.concat([old] + months, ignore_index=True)
    ext = ext.drop_duplicates("time", keep="first")
    ext.to_csv(wis, index=False)
    print(f"WIS: {len(old)} export rows + {sum(map(len, months))} archive rows "
          f"-> {wis.name} ({ext['time'].iloc[0]} to {ext['time'].iloc[-1]})")

    # Storms: the unchanged builder reading the extended records
    import runpy
    _env.DUCK_GAUGE_FILE, _env.WIS_FILE = duck, wis
    argv = sys.argv
    sys.argv = [str(BUILDER), "--start-year", "2010", "--end-year", "2026", *STORM_ARGS]
    try:
        runpy.run_path(str(BUILDER), run_name="__main__")
    finally:
        sys.argv = argv

    rslr = fit_rslr(2010, 2026)
    (INPUTS / "rslr_2010_2026.json").write_text(json.dumps(rslr, indent=2))
    print(json.dumps(rslr, indent=2))


# Subprocess entry: patch the period and the target, then run the unchanged runner
def _launch(start):
    import runpy
    from site_layer import hatteras_site_config as sc
    from site_layer import hat_env_forcings as _env
    import cascade_pipeline.coastsat_lowess as cl
    from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT
    end, (w0, w1) = PERIODS[start]
    period = sc.HATTERAS_PERIODS[start]
    if end != period["end_year"]:
        period["end_year"] = end
        period["storm_file"] = _env.init_relpath(_env.storm_series_file(start, end))
        rslr = json.loads((INPUTS / f"rslr_{start}_{end}.json").read_text())
        period["sea_level_rise_rate"] = rslr["config_m_yr"]
    original = cl.CoastSatDataset

    # The runner lists its four targets by start year; this one's points at the candidate window
    def dataset(label, period_start, csv_path, **kw):
        if period_start == start:
            label = f"CoastSat LRR ({w0}-{w1})"
            csv_path = str(COASTSAT_LRR_ROOT / f"{w0}_{w1}" / "transect_lrr_full.csv")
        return original(label=label, period_start=period_start, csv_path=csv_path, **kw)

    cl.CoastSatDataset = dataset
    print(f"candidate window: {start}-{end}, graded on {w0}_{w1}; "
          f"storms {period['storm_file']}; RSLR {period['sea_level_rise_rate']}")
    sys.argv = [str(HINDCAST)]
    runpy.run_path(str(HINDCAST), run_name="__main__")


# One run through _launch in its own process; returns its tag
def run(start, preset, label, override=None):
    tag = f"{TAG}/runs/{RUN_KEY[start]}_{preset}_{label}"
    env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUNBUFFERED="1",
               HAT_IGNORE_SETTINGS="1",
               HAT_START_YEAR=str(start), HAT_SOURCE_SINK_PRESET=preset,
               HAT_SCENARIO=SCENARIO, HAT_RELOCATIONS="False",
               HAT_GROIN_ENABLED="False", HAT_RUN_KIND="experiment",
               HAT_RUN_TAG=tag, HAT_SAVE_MODEL_STATE="False", HAT_OVERWRITE="False")
    if override:
        env["HAT_BE_OVERRIDE"] = override
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    with open(logs / f"{RUN_KEY[start]}_{preset}_{label}.log", "w", encoding="utf-8") as log:
        code = subprocess.run([sys.executable, str(_HERE), "_launch", str(start)],
                              stdout=log, stderr=subprocess.STDOUT, env=env,
                              cwd=PROJECT_ROOT).returncode
    if code:
        raise SystemExit(f"run {tag} failed ({code}); see {log.name}")
    print(f"done  {tag}")
    return tag


# The run name a tag holds, read from the index
def run_name(tag):
    from cascade_pipeline.run_registry import load_run_index
    idx = load_run_index(PROJECT_ROOT / "output" / "raw_runs" / "run_index.csv")
    rows = idx[(idx["kind"] == "experiment") & (idx["tag"] == tag)]
    if rows.empty:
        raise SystemExit(f"no run with tag {tag} in run_index.csv")
    return str(rows.iloc[-1]["run_name"])


# zeroBE on the candidate periods (all, or --period)
def cmd_zero(args):
    for start in ([args.period] if args.period else PERIODS):
        run(start, "zeroBE", "zero")


# The solver's next step, in-process, with the period's end year patched; None once both ends close
def solve_step(start, steps):
    from site_layer import hatteras_site_config as sc
    sys.path.insert(0, str(SOLVER_DIR))
    import be_edge_domain_solve as solver
    end, window = PERIODS[start]
    sc.HATTERAS_PERIODS[start]["end_year"] = end
    names = [run_name(t) for t in steps]
    buf = io.StringIO()
    with redirect_stdout(buf):
        nxt = solver.report(start, names, "edgeBE", ["experiment"] * len(steps),
                            steps, target_source="coastsat", coastsat_window=window)
    text = buf.getvalue()
    print(text)
    with open(EXP_DIR / "logs" / f"{RUN_KEY[start]}_solve.txt", "a", encoding="utf-8") as f:
        f.write(text)
    # The solver flags each end under 0.02 m/yr; both flagged is solved
    return None if text.count("CONVERGED") >= 2 else nxt


# Secant from the adopted ends until both ends close to 0.02 m/yr
def cmd_solve(args):
    from site_layer.hatteras_site_config import HATTERAS_BE_EDGE_ONLY
    start = args.period
    steps = []
    k = 0
    while (EXP_DIR / "runs" / f"{RUN_KEY[start]}_edgeBE_step{k}").is_dir():
        steps.append(f"{TAG}/runs/{RUN_KEY[start]}_edgeBE_step{k}")
        k += 1
    if not steps:
        # --seed starts from a known solve instead of the adopted ends; the secant lands on the same answer
        g1, g90 = args.seed if args.seed else HATTERAS_BE_EDGE_ONLY[start]
        steps.append(run(start, "edgeBE", "step0", f"1={g1},90={g90}"))
    while True:
        nxt = solve_step(start, steps)
        if not nxt or len(steps) > MAX_STEPS:
            print("solved" if not nxt else "stopped at MAX_STEPS")
            break
        override = ",".join(f"{g}={r:.4f}" for g, r in sorted(nxt.items()))
        steps.append(run(start, "edgeBE", f"step{len(steps)}", override))
    print(f"final step: {steps[-1]}")


# Run: one subcommand
def main(argv=None):
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(int(sys.argv[2]))
        return 0
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 2)[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("build")
    sub.add_parser("forcing")
    z = sub.add_parser("zero")
    z.add_argument("--period", type=int, choices=sorted(PERIODS))
    s = sub.add_parser("solve")
    s.add_argument("--period", type=int, choices=sorted(PERIODS), required=True)
    s.add_argument("--seed", type=float, nargs=2, metavar=("GIS1", "GIS90"),
                   help="first-step end rates, default the adopted ones")
    a = ap.parse_args(argv)
    {"build": cmd_build, "forcing": cmd_forcing, "zero": cmd_zero,
     "solve": cmd_solve}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
