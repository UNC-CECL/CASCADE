#!/usr/bin/env python3
"""
Scores each arm's modelled 2004 road setback against the measured 2004 setback.

    python scripts/hatteras_ms/experiments/HAT_score_road_position.py --arms pea1989base

Complements the relocation-timing score, which cannot rank the setback sources. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-27
"""

from __future__ import annotations

import argparse
import csv
import glob
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
# Repo root, found by searching upward
REPO = next(_p for _p in HERE.parents if (_p / 'pyproject.toml').exists())
sys.path.insert(0, str(REPO / "scripts"))

# --- CONFIG ------------------------------------------------------------------
BUFFER = 15
OUT = REPO / "output" / "raw_runs" / "experiments" / "topography-and-domains" / "2026-09-02-pea-island-row-insert-control" / "results"
# -----------------------------------------------------------------------------


# An arm's 1984-2004 calibBE run
def load(arm):
    hits = glob.glob(str(REPO / "output" / "raw_runs" / arm / "1984_2004"
                         / "calibBE" / "*" / "*.npz"))
    if not hits:
        raise SystemExit("\nno run for arm {!r}\n".format(arm))
    return np.load(hits[0], allow_pickle=True)["cascade"][0]


# Run: score every arm, write the table
def main() -> None:
    ap = argparse.ArgumentParser()
    # pea1989base is the only experiment arm with runs left
    ap.add_argument("--arms", default="pea1989base")
    ap.add_argument("--suffix", default="")
    args = ap.parse_args()

    from site_layer.hatteras_site_config import HATTERAS_RELOCATION_CHECK_2004 as MEASURED

    rows = []
    print("MEASURED 2004 setback is the observation; modelled is the run's final year.\n")
    for arm in args.arms.split(","):
        arm = arm.strip()
        if not arm:
            continue
        c = load(arm + args.suffix)
        print("=" * 70)
        print("ARM {}".format(arm + args.suffix))
        print("=" * 70)
        print("{:>4} {:>10} {:>10} {:>9}".format("GIS", "modelled", "measured", "error"))
        errs = []
        for D in sorted(MEASURED):
            mgr = c.roadways[D + BUFFER - 1]
            if mgr is None:
                continue
            sb = np.asarray(getattr(mgr, "_road_setback_TS", []), dtype=float)
            if sb.size == 0:
                continue
            model = float(sb[-1])
            meas = float(MEASURED[D])
            err = model - meas
            errs.append(err)
            print("{:>4} {:>9.0f}m {:>9.0f}m {:>+8.0f}m".format(D, model, meas, err))
            rows.append({"arm": arm, "domain": D, "modelled_2004_m": model,
                         "measured_2004_m": meas, "error_m": err})
        e = np.array(errs)
        print("\n  n {}  |  mean error {:+.1f} m  median {:+.1f}  MAE {:.1f} m  "
              "RMSE {:.1f} m\n".format(e.size, e.mean(), np.median(e),
                                       np.abs(e).mean(), np.sqrt((e ** 2).mean())))

    OUT.mkdir(parents=True, exist_ok=True)
    # Named for the arms scored; a fixed name overwrote earlier arms.
    arms_tag = "-".join(a.strip().replace("blocks", "")
                        for a in args.arms.split(",") if a.strip())
    p = OUT / "HAT_road_position_score_{}.csv".format(arms_tag)
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote {}".format(p))


if __name__ == "__main__":
    main()
