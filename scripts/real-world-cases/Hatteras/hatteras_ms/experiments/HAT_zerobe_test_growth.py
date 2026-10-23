"""
Does the test-period growth rate (0.02-0.04) improve 2009-2025 with no source/sink, full management and the groin?

    python scripts/hatteras_ms/experiments/HAT_zerobe_test_growth.py run
    python scripts/hatteras_ms/experiments/HAT_zerobe_test_growth.py score

Asked by Hannah 2026-10-21. The 2009 start, zeroBE (no source/sink anywhere, so the
ends cannot hide anything), full management, groin on (the wave-driven
GroinCallback), Dmaxel 12.2 m, today's setup; rmin-rmax 0.02-0.04 (the runner's
2009 rate since 2026-10-21) against 0.01-0.03 (the calibration rate). The ceiling
and rates are set by the grid study's build wrapper, so both runs differ in the
rates only. Scores: interior shoreline RMSE and bias, overwash hits and
model-only cases, the 2014-2019 median dune crest error. The runner draws each
run's position-change and rate figures.

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

import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_ceiling_uniform as U  # noqa: E402
import HAT_dune_ceiling_growth_grid as G  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-21-zerobe-test-growth"
EXP_DIR = U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
RATES = {"0.02-0.04 (test-period rate)": (0.02, 0.04), "0.01-0.03 (calibration rate)": (0.01, 0.03)}
CEILING = 12.2
START = 2009
PRESET = "zeroBE"
# -----------------------------------------------------------------------------

G.TAG, G.EXP_DIR = TAG, EXP_DIR


def run_dir(variant):
    base = EXP_DIR / "runs" / variant / f"{START}_{U.PERIODS[START][0]}" / PRESET
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


def launch(variant):
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(START),
                "HAT_SOURCE_SINK_PRESET": PRESET, "HAT_SCENARIO": "full_management", "HAT_GROIN_ENABLED": "True",
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{variant}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8", "MPLBACKEND": "Agg"})
    t0 = time.time()
    with open(logs / f"{variant}.log", "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(G._HERE), "_launch", variant],
                           env=env, cwd=U.HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    print(f"{variant} exit {p.returncode}  {round((time.time() - t0) / 60, 1)} min", flush=True)


def main():
    variants = {label: G.name(CEILING, r, groin=True) for label, r in RATES.items()}
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        todo = [v for v in variants.values() if run_dir(v) is None]
        with ThreadPoolExecutor(2) as ex:
            list(ex.map(launch, todo))
    import HAT_low_growth_test_period as T
    import HAT_test_period_growth_rate as P
    rows = []
    for label, v in variants.items():
        G.TAG, G.EXP_DIR = TAG, EXP_DIR
        G.run_dir = lambda start, variant: run_dir(variant)   # zeroBE layout
        rd = run_dir(v)
        if rd is None:
            print(f"missing {v}")
            continue
        meta = json.loads(next(rd.glob("*_run_metadata.json")).read_text(encoding="utf-8"))
        e = P.survey_errors(v)
        d = T.cells(START, v)
        o, mo = d.observed.astype(int), d.model.astype(int)
        rows.append(dict(rate=label, shoreline_rmse=round(float(meta["skill"]["rmse_interior_m_yr"]), 3),
                         shoreline_bias=round(float(meta["skill"]["mean_bias_interior_m_yr"]), 3),
                         hits=int(((o == 1) & (mo == 1)).sum()), observed=int(o.sum()),
                         model_only=int(((o == 0) & (mo == 1)).sum()),
                         irene=T.major(d, "2011-08-27"),
                         median_crest_err_2014_2019=" / ".join(f"{x:+.2f}" for x in e.median_err),
                         figures=str(rd / "figures")))
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "zerobe_test_growth.csv", index=False)
    with pd.option_context("display.width", 300, "display.max_columns", 20, "display.max_colwidth", 200):
        print(t.drop(columns="figures").to_string(index=False))


if __name__ == "__main__":
    main()
