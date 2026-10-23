"""
How much do the dune growth rate bounds (rmin, rmax) change the dunes and the overwash?

    python scripts/hatteras_ms/experiments/HAT_dune_growth_rate.py run
    python scripts/hatteras_ms/experiments/HAT_dune_growth_rate.py score

Sets rmin/rmax on the build_cascade call, in-process only, to 0.25-0.65,
0.35-0.85 and 0.75-1.15; managed edgeBE runs of both periods. The current
0.55-0.95 control is the per-cell run of HAT_dune_ceiling_uniform.py, scored
here alongside. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-13
"""
from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_ceiling_uniform as U  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-13-dune-growth-rate"
RANGES = {"r025_065": (0.25, 0.65), "r035_085": (0.35, 0.85), "r075_115": (0.75, 1.15)}
# -----------------------------------------------------------------------------


# Child process: override rmin/rmax on the build_cascade call, then run the unchanged hindcast
def _launch(variant):
    import runpy
    import cascade_pipeline.hindcast as H
    lo, hi = RANGES[variant]
    orig = H.build_cascade

    def build_cascade(*a, **k):
        k["rmin"] = [lo] * len(k["rmin"])
        k["rmax"] = [hi] * len(k["rmax"])
        print(f"EXPERIMENT rmin/rmax = {lo}/{hi} (runner's 0.55/0.95 overridden)", flush=True)
        return orig(*a, **k)

    H.build_cascade = build_cascade
    sys.path.insert(0, str(U.HINDCAST.parent))
    sys.argv = [str(U.HINDCAST)]
    runpy.run_path(str(U.HINDCAST), run_name="__main__")


# Re-point the uniform study at this one: its launcher, run folders and scorer, with these variants
def adopt():
    uni_dir = U.EXP_DIR
    U.TAG, U.EXP_DIR, U._HERE = TAG, U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG, _HERE
    U.VARIANTS = {"r025_065": None, "r035_085": None, "percell": None, "r075_115": None}
    own = U.run_dir

    def run_dir(start, variant):
        if variant in RANGES:
            return own(start, variant)
        base = uni_dir / "runs" / variant / f"{start}_{U.PERIODS[start][0]}" / U.PRESET
        hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
        return hits[-1] if hits else None

    U.run_dir = run_dir


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(sys.argv[2])
        return
    adopt()
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        from concurrent.futures import ThreadPoolExecutor
        jobs = [(s, v) for s in U.PERIODS for v in RANGES if U.run_dir(s, v) is None]
        with ThreadPoolExecutor(U.WORKERS) as ex:
            list(ex.map(U.launch, jobs))
    else:
        U.score()


if __name__ == "__main__":
    main()
