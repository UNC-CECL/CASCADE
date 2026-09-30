#!/usr/bin/env python3
"""
Scores each arm's predicted NC-12 relocation year against the 1989 and 1999 events.

    python scripts/hatteras_ms/experiments/HAT_score_relocation_timing.py --arms pea1989base
    python scripts/hatteras_ms/experiments/HAT_score_relocation_timing.py --holdout 1989

Only meaningful with the prescribed road events off; --holdout scores one
event as an out-of-sample test. Details: scripts/hatteras_ms/experiments/README.md.

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
# --- CONFIG ------------------------------------------------------------------
BUFFER = 15
START, END = 1984, 2004
EVENTS = {**{d: 1999 for d in range(9, 15)},
          **{d: 1989 for d in range(84, 88)}}
OUT = REPO / "output" / "raw_runs" / "experiments" / "topography-and-domains" / "2026-09-02-pea-island-row-insert-control" / "results"
# -----------------------------------------------------------------------------


# An arm's 1984-2004 calibBE run
def load(arm):
    hits = glob.glob(str(REPO / "output" / "raw_runs" / arm / "1984_2004"
                         / "calibBE" / "*" / "*.npz"))
    if not hits:
        raise SystemExit("\nno run for arm {!r}\n".format(arm))
    return np.load(hits[0], allow_pickle=True)["cascade"][0]


# A domain's first relocation year and its starting setback
def first_relocation(c, gis):
    mgr = c.roadways[gis + BUFFER - 1]
    if mgr is None:
        return None, np.nan
    rel = np.asarray(getattr(mgr, "_road_relocated_TS", []), dtype=float)
    sb = np.asarray(getattr(mgr, "_road_setback_TS", []), dtype=float)
    t0 = float(sb[0]) if sb.size else np.nan
    hit = np.flatnonzero(np.nan_to_num(rel) > 0)
    return (START + int(hit[0]) if hit.size else None), t0


# Run: score every arm, write the table named for the arms
def main() -> None:
    ap = argparse.ArgumentParser()
    # pea1989base is the only experiment arm with runs left
    ap.add_argument("--arms", default="pea1989base")
    ap.add_argument("--suffix", default="noreloc")
    ap.add_argument("--holdout", type=int, choices=(1989, 1999), default=None,
                    help="score ONLY the other event, as an out-of-sample test")
    args = ap.parse_args()
    arms = [a.strip() for a in args.arms.split(",") if a.strip()]

    domains = sorted(d for d, y in EVENTS.items()
                     if args.holdout is None or y != args.holdout)
    if args.holdout:
        print("HOLDOUT: {} withheld; scoring only the {} block\n".format(
            args.holdout, 1999 if args.holdout == 1989 else 1989))

    rows = []
    for arm in arms:
        c = load(arm + args.suffix)
        print("=" * 78)
        print("ARM {}".format(arm + args.suffix))
        print("=" * 78)
        print("{:>4} {:>10} {:>11} {:>8} {:>8}".format(
            "GIS", "setback t0", "predicted", "actual", "error"))
        errs, censored = [], 0
        for D in domains:
            pred, t0 = first_relocation(c, D)
            actual = EVENTS[D]
            if pred is None:
                censored += 1
                print("{:>4} {:>9.0f}m {:>11} {:>8} {:>8}".format(
                    D, t0, ">2004", actual, "censored"))
                err = np.nan
            else:
                err = pred - actual
                errs.append(err)
                print("{:>4} {:>9.0f}m {:>11} {:>8} {:>+8}".format(
                    D, t0, pred, actual, err))
            rows.append({"arm": arm, "domain": D, "setback_t0_m": t0,
                         "predicted": pred if pred else "",
                         "actual": actual, "error_yr": err,
                         "censored": pred is None})
        e = np.array(errs, dtype=float)
        if e.size:
            print("\n  n scored {}  censored {}  |  mean error {:+.1f} yr  "
                  "median {:+.1f}  MAE {:.1f}".format(
                      e.size, censored, e.mean(), np.median(e), np.abs(e).mean()))
        else:
            print("\n  nothing relocated inside the run; all {} censored".format(
                censored))
        print()

    OUT.mkdir(parents=True, exist_ok=True)
    # Named for the arms scored, so a later run cannot overwrite these scores
    tag = "_holdout{}".format(args.holdout) if args.holdout else ""
    tag += "_" + "-".join(a.replace("blocks", "") for a in arms)
    p = OUT / "HAT_relocation_timing_score{}.csv".format(tag)
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote {}".format(p))


if __name__ == "__main__":
    main()
