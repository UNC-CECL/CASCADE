"""
Do smoother or longer buffers stop the ring closure eroding the island's end domains?

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/ends/HAT_buffer_closure_test.py run
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/ends/HAT_buffer_closure_test.py score

Why (2026-10-21): BRIE's alongshore domain is periodic, so the buffers carry an
invented coast 5.4 km back across the shore from GIS 90 to GIS 1. The cubic
Hermite closure (hindcast.close_offset_ring) bends most right at the joins: the
buffers bend up to 4.1 deg a domain against the real coast's 1.1, the north
buffer loses 82 / 133 m (1996-2009 / 2009-2025, no source/sink) and drags GIS
88-90 and, through the wrap, GIS 1-3 down with it; in 2009-2025 the buffers lose
3.4x the real island's shoreline area. The end source/sink rates have been
compensating for that sink.

Variants (Hannah, 2026-10-21), all zeroBE, full management, groin on, today's
settings from the runner, both windows:
    A   the current cubic closure, 15 buffers a side (reference)
    B   a quintic Hermite closure: slope matched at each join to the mean of the
        last SLOPE_STEPS real steps, curvature zero there, so the coast leaves the
        island straight and turns inside the buffer; 15 buffers a side
    C   B with 30 buffers a side (150 padded domains): the same turn over twice
        the length
In-process only: close_offset_ring is replaced in the run's own process (B, C),
and for C the site geometry and build_island_offset (which re-closes the
120-value offset file's real domains over 30 buffers). Nothing in the main code
changes. Scores: mean shoreline change at GIS 1-3 and 88-90 and in each buffer,
the buffers' shoreline area change against the real island's, the sharpest bend
at each join, interior shoreline RMSE and bias.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-21
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

_HERE = Path(__file__).resolve()
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in _HERE.parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
HINDCAST = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"

# --- CONFIG ------------------------------------------------------------------
TAG = "source-sink/2026-10-21-buffer-closure"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
VARIANTS = {"A_cubic15": ("cubic", 15), "B_quintic15": ("quintic", 15), "C_quintic30": ("quintic", 30)}
STARTS = (1996, 2009)
SLOPE_STEPS = 3                     # real steps averaged for the join slope (B, C)
OFFSET_FILE_BUFFERS = 15            # the padded offset files hold 15 buffers a side
PRESET = "zeroBE"
WORKERS = 4
DOMAIN_M = 500.0
# -----------------------------------------------------------------------------


# The quintic closure: real[-1] -> right buffer -> (wrap) -> left buffer -> real[0], slope matched, curvature 0
def quintic_closure(real_m, buffers):
    import numpy as np
    n = 2 * buffers + 1
    p0, p1 = real_m[-1], real_m[0]
    m0 = np.mean(np.diff(real_m[-(SLOPE_STEPS + 1):])) * n
    m1 = np.mean(np.diff(real_m[:SLOPE_STEPS + 1])) * n
    t = np.arange(1, n) / n
    h0 = 1 - 10 * t ** 3 + 15 * t ** 4 - 6 * t ** 5
    h1 = t - 6 * t ** 3 + 8 * t ** 4 - 3 * t ** 5
    h3 = 10 * t ** 3 - 15 * t ** 4 + 6 * t ** 5
    h4 = -4 * t ** 3 + 7 * t ** 4 - 3 * t ** 5
    path = h0 * p0 + h1 * m0 + h3 * p1 + h4 * m1       # the second-derivative terms are zero
    return path[:buffers], path[buffers:]


# Child process: patch the closure (and for C the geometry and offset re-closure), then run the unchanged hindcast
def _launch(variant):
    import dataclasses
    import runpy
    import numpy as np
    sys.path.insert(0, str(_HAT_DIR))
    closure, buffers = VARIANTS[variant]
    import site_layer.hatteras_site_config as cfg
    if buffers != cfg.HATTERAS_DOMAINS.num_buffer_domains:
        cfg.HATTERAS_DOMAINS = dataclasses.replace(cfg.HATTERAS_DOMAINS, num_buffer_domains=buffers)
    import cascade_pipeline.hindcast as H
    if closure == "quintic":
        H.close_offset_ring = quintic_closure
    if buffers != OFFSET_FILE_BUFFERS:
        def build_island_offset(offset_path, geometry):
            offset_m = np.loadtxt(offset_path, skiprows=1, delimiter=",")
            real_m = offset_m[OFFSET_FILE_BUFFERS:offset_m.size - OFFSET_FILE_BUFFERS]
            if real_m.size != geometry.num_real_domains:
                raise ValueError(f"expected {geometry.num_real_domains} real domains")
            return H.pad_offset_ring(real_m, geometry.num_buffer_domains)
        H.build_island_offset = build_island_offset
    print(f"EXPERIMENT buffer closure {closure}, {buffers} buffers a side "
          f"({cfg.HATTERAS_DOMAINS.total_domains} padded domains)", flush=True)
    sys.path.insert(0, str(HINDCAST.parent))
    sys.argv = [str(HINDCAST)]
    runpy.run_path(str(HINDCAST), run_name="__main__")


def run_dir(variant, start):
    sys.path.insert(0, str(_HAT_DIR))
    from site_layer.hatteras_site_config import HATTERAS_PERIODS
    base = EXP_DIR / "runs" / variant / f"{start}_{HATTERAS_PERIODS[start]['end_year']}" / PRESET
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


def launch(job):
    variant, start = job
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(start),
                "HAT_SOURCE_SINK_PRESET": PRESET, "HAT_SCENARIO": "full_management", "HAT_GROIN_ENABLED": "True",
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{variant}",
                "HAT_SAVE_MODEL_STATE": "0", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8", "MPLBACKEND": "Agg"})
    t0 = time.time()
    with open(logs / f"{variant}_{start}.log", "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", variant],
                           env=env, cwd=HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    r = dict(variant=variant, start=start, returncode=p.returncode, minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(r) + "\n")
    print(f"{start} {variant} exit {p.returncode}  {r['minutes']} min", flush=True)


def score_run(variant, start):
    import numpy as np
    rd = run_dir(variant, start)
    m = np.load(next(rd.glob("*_shoreline_matrix.npy")))       # states x padded domains, x_s landward +
    nb = (m.shape[1] - 90) // 2
    a, b = nb, nb + 90
    x0 = m[0]
    ang = np.degrees(np.arctan(np.diff(np.r_[x0, x0[0]]) / DOMAIN_M))
    bend = np.diff(np.r_[ang[-1], ang])
    chg = -(m[-1] - m[0])                                         # seaward +
    meta = json.loads(next(rd.glob("*_run_metadata.json")).read_text(encoding="utf-8"))
    return dict(window=f"{start}", variant=variant, buffers=nb,
                gis1_3_m=round(float(chg[a:a + 3].mean()), 1), gis88_90_m=round(float(chg[b - 3:b].mean()), 1),
                south_buffer_m=round(float(chg[:a].mean()), 1), north_buffer_m=round(float(chg[b:].mean()), 1),
                buffer_area_km2=round(float((chg[:a].sum() + chg[b:].sum()) * DOMAIN_M / 1e6), 2),
                real_area_km2=round(float(chg[a:b].sum() * DOMAIN_M / 1e6), 2),
                max_bend_south_join=round(float(np.abs(bend[a - 3:a + 3]).max()), 1),
                max_bend_north_join=round(float(np.abs(bend[b - 3:b + 3]).max()), 1),
                max_bend_buffers=round(float(np.abs(np.r_[bend[:a], bend[b:]]).max()), 1),
                shoreline_rmse=round(float(meta["skill"]["rmse_interior_m_yr"]), 3),
                shoreline_bias=round(float(meta["skill"]["mean_bias_interior_m_yr"]), 3))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        return _launch(sys.argv[2])
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        todo = [(v, s) for v in VARIANTS for s in STARTS if run_dir(v, s) is None]
        with ThreadPoolExecutor(WORKERS) as ex:
            list(ex.map(launch, todo))
    import pandas as pd
    rows = [score_run(v, s) for s in STARTS for v in VARIANTS if run_dir(v, s) is not None]
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "buffer_closure_scores.csv", index=False)
    with pd.option_context("display.width", 300, "display.max_columns", 30):
        print(t.to_string(index=False))


if __name__ == "__main__":
    main()
