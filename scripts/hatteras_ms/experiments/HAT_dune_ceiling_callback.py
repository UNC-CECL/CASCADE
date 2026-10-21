"""
Can a CASCADE-side callback reproduce the per-cell dune ceilings on base Barrier3D dune growth?

    python scripts/hatteras_ms/experiments/HAT_dune_ceiling_callback.py run
    python scripts/hatteras_ms/experiments/HAT_dune_ceiling_callback.py score

Turns the fork's per-cell ceilings off and puts Barrier3D on its own scalar
ceiling (each domain's median cell ceiling), then rescales each cell's growth
rate every year so the front dune row grows as if toward its own starting
crest. In-process only; managed edgeBE runs of both periods, scored against the
fork's per-cell runs of HAT_dune_ceiling_uniform.py. Details:
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
TAG = "dunes/2026-10-13-dune-ceiling-callback"
FLOOR_DAM = 0.05             # lowest cell ceiling, 0.5 m above the berm (the fork's DuneCeilingFloor)
SCALAR = {"callback": "median", "callback2x": "twice_max"}   # the domain's scalar Barrier3D ceiling
EPS = 1e-9                   # below this the scalar logistic term is zero and cannot carry a growth rate
# -----------------------------------------------------------------------------


# Per-cell dune ceilings through the growth rate, on base Barrier3D's scalar ceiling
class PerCellCeilingCallback:
    """Barrier3D grows each dune row by r * D * (1 - D/Dmax) with one Dmax per
    domain. Setting a cell's rate to r * (1 - D/K) / (1 - D/Dmax) for one year
    makes its front row grow by r * D * (1 - D/K): logistic toward its own
    ceiling K. The rescaled rate is in place only during the Barrier3D step;
    the rest of the year growthparam holds the plain rates, so the road and
    beach managers read and reset what they always have."""

    def __init__(self, scalar="median"):
        self.scalar = scalar

    # Each cell's ceiling from the starting dunes; the domain's scalar ceiling is their median,
    # or twice their maximum so that no dune reaches it and the rescaling stays well conditioned
    def attach(self, cascade):
        for b in cascade.barrier3d:
            K = np.maximum(np.asarray(b.DuneDomain[0], float).max(axis=1), FLOOR_DAM)
            b._DuneCeiling = None                     # base Barrier3D's scalar path
            b._Dmax = float(np.median(K)) if self.scalar == "median" else 2.0 * float(K.max())
            b._Dmaxel = b._Dmax + b._BermEl
            b._hat_ceiling = K
            b._hat_growth = np.array(b.growthparam, float)
            self.rescale(b)

    # The growth rates that make next year's front-row growth logistic toward each cell's ceiling
    def rescale(self, b):
        # the heights growth will see: SeaLevel() first lowers the dunes by next year's RSLR
        rslr = b._RSLR[b.time_index] if b.time_index < len(b._RSLR) else 0.0
        D = np.asarray(b.DuneDomain[b.time_index - 1], float)[:, 0] - rslr
        den = 1 - D / b._Dmax
        f = np.where(np.abs(den) > EPS, (1 - D / b._hat_ceiling) / np.where(np.abs(den) > EPS, den, 1), 0.0)
        b._hat_rescaled = b._hat_growth * f

    # Barrier3D step with the rescaled rates in place, plain rates restored after
    def wrap_step(self, step):
        def batch(b):
            b.growthparam = b._hat_rescaled
            out = step(b)
            out[3].growthparam = out[3]._hat_growth
            return out
        return batch

    # After the managers: keep their rates as the plain rates, rescale for next year
    def after_update(self, cascade):
        for b in cascade.barrier3d:
            b._hat_growth = np.array(b.growthparam, float)
            self.rescale(b)


# Child process: base scalar ceilings plus the callback, then the unchanged hindcast
def _launch(variant):
    import runpy
    import cascade.cascade_groin as CG
    import cascade_pipeline.hindcast as H
    cb = PerCellCeilingCallback(SCALAR[variant])
    CG.batchB3D = cb.wrap_step(CG.batchB3D)
    orig_update = CG.Cascade.update

    def update(self):
        orig_update(self)
        cb.after_update(self)

    CG.Cascade.update = update
    orig_build = H.build_cascade

    def build_cascade(*a, **k):
        c = orig_build(*a, **k)
        cb.attach(c)
        print(f"EXPERIMENT per-cell ceilings via growth-rate callback, scalar {cb.scalar} (fork ceilings off)",
              flush=True)
        return c

    H.build_cascade = build_cascade
    sys.path.insert(0, str(U.HINDCAST.parent))
    sys.argv = [str(U.HINDCAST)]
    runpy.run_path(str(U.HINDCAST), run_name="__main__")


# Re-point the uniform study at this one: its launcher, run folders and scorer, with these variants
def adopt():
    uni_dir = U.EXP_DIR
    U.TAG, U.EXP_DIR, U._HERE = TAG, U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG, _HERE
    U.VARIANTS = {"percell": None, "callback": None, "callback2x": None}
    own = U.run_dir

    def run_dir(start, variant):
        if variant in SCALAR:
            return own(start, variant)
        base = uni_dir / "runs" / variant / f"{start}_{U.PERIODS[start][0]}" / U.PRESET
        hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
        return hits[-1] if hits else None

    U.run_dir = run_dir


# How far the callback's dunes and shoreline sit from the fork's, cell by cell
def compare():
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    for start in U.PERIODS:
        a = U.run_dir(start, "percell")
        for v in SCALAR:
            b = U.run_dir(start, v)
            if a is None or b is None:
                continue
            ca, cb = U.cell_crests(U.load(a), DOM), U.cell_crests(U.load(b), DOM)
            d = np.concatenate([(y - x).ravel() for x, y in zip(ca, cb)])
            sa = np.load(a / f"{a.name}_shoreline_matrix.npy")
            sb = np.load(b / f"{b.name}_shoreline_matrix.npy")
            print(f"  {start} {v}: crest difference, every cell-year: median {np.median(np.abs(d)):.3f} m, "
                  f"P95 {np.percentile(np.abs(d), 95):.3f} m, max {np.abs(d).max():.2f} m; "
                  f"shoreline max {np.nanmax(np.abs(sb - sa)):.2f} m", flush=True)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(sys.argv[2])
        return
    adopt()
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        from concurrent.futures import ThreadPoolExecutor
        jobs = [(s, v) for s in U.PERIODS for v in SCALAR if U.run_dir(s, v) is None]
        with ThreadPoolExecutor(U.WORKERS) as ex:
            list(ex.map(U.launch, jobs))
    else:
        U.score()
        compare()


if __name__ == "__main__":
    main()
