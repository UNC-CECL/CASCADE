"""
Does one dune ceiling per domain, from that domain's own starting crests, beat the per-cell ceilings?

    python scripts/hatteras_ms/experiments/HAT_dune_ceiling_domain_value.py run
    python scripts/hatteras_ms/experiments/HAT_dune_ceiling_domain_value.py score

Replaces each domain's 50 per-cell ceilings with one value, the median or 75th
percentile of its starting crests, in-process only; managed edgeBE runs of both
periods. The per-cell control and the island-wide values are the runs of
HAT_dune_ceiling_uniform.py, scored here alongside. Details:
scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-13
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_ceiling_uniform as U  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "storms-overwash-and-dunes/2026-10-13-dune-ceiling-domain-value"
PERCENTILES = {"dom50": 50, "dom75": 75}      # of the domain's starting cell crests
FLOOR_DAM = 0.05                              # DuneCeilingFloor, 0.5 m above the berm
# -----------------------------------------------------------------------------


# Every cell in a domain grows toward one value: a percentile of that domain's starting crests
def set_domain_value(cascade, q):
    for b in cascade.barrier3d:
        crest = np.asarray(b._DuneDomain[0]).max(axis=1)              # dam above the berm, per cell
        h = max(float(np.percentile(crest, q)), FLOOR_DAM)
        b._DuneCeiling = None
        b._Dmax = h
        b._Dmaxel = h + b._BermEl


# Child process: wrap build_cascade to set the ceilings, then run the unchanged hindcast
def _launch(variant):
    import runpy
    import cascade_pipeline.hindcast as H
    orig = H.build_cascade

    def build_cascade(*a, **k):
        c = orig(*a, **k)
        set_domain_value(c, PERCENTILES[variant])
        print(f"EXPERIMENT per-domain dune ceiling = P{PERCENTILES[variant]} of each domain's starting crests",
              flush=True)
        return c

    H.build_cascade = build_cascade
    sys.path.insert(0, str(U.HINDCAST.parent))
    sys.argv = [str(U.HINDCAST)]
    runpy.run_path(str(U.HINDCAST), run_name="__main__")


# Re-point the uniform study at this one: its launcher, run folders and scorer, with these variants
def adopt():
    own_dir = U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
    uni_dir = U.EXP_DIR
    U.TAG, U.EXP_DIR, U._HERE = TAG, own_dir, _HERE
    U.VARIANTS = {**{v: None for v in PERCENTILES}, "percell": None, "uni50": 5.0, "uni55": 5.5}
    own = U.run_dir

    def run_dir(start, variant):
        if variant in PERCENTILES:
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
        jobs = [(s, v) for s in U.PERIODS for v in PERCENTILES if U.run_dir(s, v) is None]
        with ThreadPoolExecutor(U.WORKERS) as ex:
            list(ex.map(U.launch, jobs))
    else:
        U.score()


if __name__ == "__main__":
    main()
