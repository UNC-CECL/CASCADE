"""
Which wave climate fits each window best, on the current dunes, storms, buffers and groin?

    python scripts/hatteras_ms/experiments/HAT_wave_grid.py run      # waits for the matrix driver, then solves + grid
    python scripts/hatteras_ms/experiments/HAT_wave_grid.py score

Why (2026-10-21): the last wave sweep (10-15) ran before the v3_240 storms, the new dunes
(Dmaxel 12.2, growth 0.01-0.03), the quintic 30-buffer ring, the spread groin and beach slope
0.06, and its model states are gone. Hannah asked for every combination, best per period.

Grid (Hannah, 10-21): Hs 1.5 / 2.0 / 2.5 / 3.0 x asymmetry 0.5 / 0.55 / 0.6 / 0.65 x high-angle
fraction 0.45 / 0.5 / 0.55, Tp fixed at 7.5 (flat in all four earlier sweeps); 48 cells per window.
The end rates are valid only at the Hs they were solved at, so each Hs gets its own ends:
Hs 2.0 uses the config's (HATTERAS_BE_EDGE_ONLY, solved 10-21); 1.5, 2.5 and 3.0 are solved
here the way HAT_end_solve_net_change.py solves them (zeroBE base, then secant per end on the
DEM-to-DEM net change, full management, no groin, relocations off, option-A asymmetry and
high-angle). Asymmetry and high-angle cells use their Hs's ends as they stand.
Grid cells: edgeBE, full management, groin on, relocations off -- the matrix baseline. Scores
come from each run's own metadata (interior LRR RMSE against CoastSat LOWESS 7, bias, endpoint
RMSE, roads drowned), so a run_index race between parallel runs cannot change them.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-21
"""
from __future__ import annotations

import itertools
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve()
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in _HERE.parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
HINDCAST = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"
sys.path.insert(0, str(_HAT_DIR))
sys.path.insert(0, str(_HAT_DIR / "input_prep" / "7-source-sink" / "2-calibrate"))
from site_layer.hatteras_site_config import HATTERAS_BE_EDGE_ONLY, HATTERAS_PERIODS  # noqa: E402
import be_edge_domain_solve as solver  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "wave-climate/2026-10-21-wave-grid"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
DRIVER_LOCK = PROJECT_ROOT / "output" / "logs" / "driver" / "driver.lock"
STARTS = (1996, 2009)
HS = (1.5, 2.0, 2.5, 3.0)
ASYM = (0.5, 0.55, 0.6, 0.65)
AHF = (0.45, 0.5, 0.55)
HS_DEFAULT = 2.0             # its ends are the config's
MAX_STEPS = 8                # as HAT_end_solve_net_change.py
WORKERS = 4
# -----------------------------------------------------------------------------


def window(start):
    return f"{start}_{HATTERAS_PERIODS[start]['end_year']}"


def token(x):
    return f"{x:g}".replace(".", "p")


def cell_name(hs, asym, ahf):
    return f"Hs{token(hs)}_asym{token(asym)}_ahf{token(ahf)}"


def run_dir(tag, start, preset):
    base = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / tag / window(start) / preset
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


# One runner process; skipped when the run is already on disk (so `run` resumes)
def launch(tag, start, preset, groin, hs, asym=0.6, ahf=0.5, override=None):
    done = run_dir(tag, start, preset)
    if done is not None:
        return done
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(start), "HAT_SOURCE_SINK_PRESET": preset,
                "HAT_SCENARIO": "full_management", "HAT_RELOCATIONS": "False",
                "HAT_GROIN_ENABLED": str(groin), "HAT_HS": repr(hs), "HAT_WAVE_ASYMMETRY": repr(asym),
                "HAT_WAVE_ANGLE_HIGH_FRACTION": repr(ahf),
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": tag,
                "HAT_SAVE_MODEL_STATE": "0", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8", "MPLBACKEND": "Agg"})
    if override:
        env["HAT_BE_OVERRIDE"] = override
    log = logs / f"{tag.rsplit('/', 2)[-2]}_{tag.rsplit('/', 1)[-1]}_{window(start)}.log"
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        code = subprocess.run([sys.executable, str(HINDCAST)], env=env, cwd=HINDCAST.parent,
                              stdout=fh, stderr=subprocess.STDOUT).returncode
    if code:
        raise RuntimeError(f"{tag} {window(start)} failed ({code}); see {log}")
    print(f"  done {tag.split('/', 2)[-1]} {window(start)}  {(time.time() - t0) / 60:.1f} min", flush=True)
    return run_dir(tag, start, preset)


# Secant per end against the DEM-to-DEM net change, as be_edge_domain_solve.report steps it
def solve_ends(job):
    start, hs = job
    target, _ = solver.load_net_change_target(start, HATTERAS_PERIODS[start]["end_year"])
    years = HATTERAS_PERIODS[start]["last_model_year"] - start + 1
    tol = solver.CONVERGED_M_YR * years
    member = f"{TAG}/ends-Hs{token(hs)}"
    model = solver.read_model_net_change(launch(f"{member}/base", start, "zeroBE", False, hs))
    hist = {g: [(0.0, model[g] - target[g])] for g in target}
    for step in range(1, MAX_STEPS + 1):
        rates = {}
        for g, pts in hist.items():
            # A closed end keeps its rate; stepping it again is noise (and a flat secant, 10-21 crash)
            if len(pts) > 1 and abs(pts[-1][1]) < tol:
                rates[g] = pts[-1][0]
                continue
            gain = solver.NOMINAL_GAIN * years
            if len(pts) > 1:
                (x0, r0), (x1, r1) = pts[-2], pts[-1]
                if x1 != x0 and r1 != r0:
                    gain = (r1 - r0) / (x1 - x0)
            rates[g] = pts[-1][0] - pts[-1][1] / gain
        override = ",".join(f"{g}={r:.4f}" for g, r in sorted(rates.items()))
        rd = launch(f"{member}/step{step}", start, "edgeBE", False, hs, override=override)
        model = solver.read_model_net_change(rd)
        # The rates the run actually imposed: a resumed step on disk may predate this code
        imposed = json.loads(next(rd.glob("*_run_metadata.json")).read_text(encoding="utf-8"))["source/sink"]
        for g in target:
            hist[g].append((float(imposed[f"rate_gis{g}_m_yr"]), model[g] - target[g]))
        line = (f"{window(start)} Hs {hs} step{step} " + "  ".join(
            f"GIS {g} {hist[g][-1][0]:+.4f} resid {hist[g][-1][1]:+.2f} m" for g in sorted(target)))
        print(line, flush=True)
        with open(EXP_DIR / "solve.txt", "a", encoding="utf-8") as fh:
            fh.write(line + "\n")
        if all(abs(hist[g][-1][1]) < tol for g in target):
            break
    # Each end at its best probe (2009 GIS 1 is noisy and may not close inside MAX_STEPS)
    best = {g: min(hist[g][1:], key=lambda p: abs(p[1])) for g in target}
    out = dict(start=start, hs=hs, tolerance_m=tol,
               **{f"gis{g}_m_yr": round(best[g][0], 4) for g in sorted(target)},
               **{f"gis{g}_resid_m": round(best[g][1], 3) for g in sorted(target)})
    with open(EXP_DIR / "ends.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(out) + "\n")
    return out


def ends_for(start, hs):
    if hs == HS_DEFAULT:
        return HATTERAS_BE_EDGE_ONLY[start]
    path = EXP_DIR / "ends.jsonl"
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()] if path.exists() else []
    rows = [r for r in rows if r["start"] == start and r["hs"] == hs]
    return (rows[-1]["gis1_m_yr"], rows[-1]["gis90_m_yr"]) if rows else None


def grid_cell(job):
    start, hs, asym, ahf = job
    if ends_for(start, hs) is None:
        raise RuntimeError(f"no solved ends for {window(start)} Hs {hs}; skipped {cell_name(hs, asym, ahf)}")
    gis1, gis90 = ends_for(start, hs)
    return launch(f"{TAG}/grid/{cell_name(hs, asym, ahf)}", start, "edgeBE", True, hs, asym, ahf,
                  override=f"1={gis1:.4f},90={gis90:.4f}")


def wait_for_matrix():
    while DRIVER_LOCK.exists():
        print(f"waiting for the matrix driver ({DRIVER_LOCK.read_text().strip()})", flush=True)
        time.sleep(300)


def run():
    EXP_DIR.mkdir(parents=True, exist_ok=True)
    wait_for_matrix()
    solves = [(s, hs) for hs in HS if hs != HS_DEFAULT for s in STARTS
              if ends_for(s, hs) is None]
    with ThreadPoolExecutor(WORKERS) as ex:
        for job, f in [(j, ex.submit(solve_ends, j)) for j in solves]:
            try:
                f.result()
            except Exception as err:                # its Hs's grid cells fail below; the rest still run
                print(f"  SOLVE FAILED {job}: {err!r}", flush=True)
    cells = [(s, *c) for s in STARTS for c in itertools.product(HS, ASYM, AHF)]
    with ThreadPoolExecutor(WORKERS) as ex:
        futures = [ex.submit(grid_cell, c) for c in cells]
        for f in futures:
            try:
                f.result()
            except Exception as err:                # one failed cell does not stop the grid
                print(f"  FAILED {err}", flush=True)


def score():
    rows = []
    for start in STARTS:
        for hs, asym, ahf in itertools.product(HS, ASYM, AHF):
            rd = run_dir(f"{TAG}/grid/{cell_name(hs, asym, ahf)}", start, "edgeBE")
            if rd is None:
                continue
            meta = json.loads(next(rd.glob("*_run_metadata.json")).read_text(encoding="utf-8"))
            sk, ix = meta["skill"], meta.get("index row", {})
            ends = ends_for(start, hs)
            rows.append(dict(window=window(start), hs=hs, asym=asym, ahf=ahf,
                             rmse_interior=float(sk["rmse_interior_m_yr"]),
                             bias_interior=float(sk["mean_bias_interior_m_yr"]),
                             endpoint_rmse_interior=float(sk["endpoint_rmse_interior_m_yr"]),
                             roads_drowned=ix.get("roads_drowned"),
                             gis1_m_yr=ends[0], gis90_m_yr=ends[1],
                             # pinned settings: a cell that differs ran on a changed runner
                             rmin_rmax=meta["dunes and roadway"]["rmin / rmax"],
                             groin=json.dumps([meta["groin"].get("strength_history"),
                                               meta["groin"].get("spread_weights")]),
                             run=str(rd.relative_to(EXP_DIR))))
    t = pd.DataFrame(rows)
    if t.empty:
        print("no grid cells on disk yet")
        return
    t["rank_in_window"] = t.groupby("window").rmse_interior.rank(method="min").astype(int)
    t = t.sort_values(["window", "rank_in_window"])
    # The other window's score for the same climate, so a per-window best can be cross-checked
    other = t.set_index(["window", "hs", "asym", "ahf"]).rmse_interior
    wins = {window(s): window(o) for s, o in zip(STARTS, STARTS[::-1])}
    t["rmse_other_window"] = [other.get((wins[w], h, a, f)) for w, h, a, f in
                              zip(t.window, t.hs, t.asym, t.ahf)]
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "wave_grid_scores.csv", index=False)
    for col in ("rmin_rmax", "groin"):
        for w, g in t.groupby("window"):
            if g[col].nunique() > 1:
                print(f"WARNING {w}: cells differ in {col}: {sorted(g[col].unique())}")
    with pd.option_context("display.width", 250, "display.max_columns", 20):
        for w, g in t.groupby("window"):
            print(f"\n{w}: {len(g)} cells, top 10 by interior RMSE (m/yr)")
            print(g.head(10)[["hs", "asym", "ahf", "rmse_interior", "bias_interior", "endpoint_rmse_interior",
                              "roads_drowned", "rmse_other_window"]].to_string(index=False))
            base = g[(g.hs == 2.0) & (g.asym == 0.6) & (g.ahf == 0.5)]
            if not base.empty:
                print(f"  option A: rank {int(base.rank_in_window.iloc[0])}, "
                      f"RMSE {base.rmse_interior.iloc[0]:.3f}")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        run()
    score()


if __name__ == "__main__":
    main()
