#!/usr/bin/env python3
r"""
HAT_score_relocation_timing.py
==============================================================================
Scores each arm's PREDICTED NC-12 relocation year against the two documented
events: 1989 Pea Island (GIS 84-87) and 1999 inter-village (GIS 9-14).

ONLY MEANINGFUL WITH THE PRESCRIBED EVENTS OFF. With HATTERAS_ROAD_EVENTS on,
1989 and 1999 are inputs, and scoring the model against its own input is
circular. Pass arms built with --relocations 0.

WHAT IS BEING DISCRIMINATED
    The two independent measurements of the 1984 offset disagree by a factor of
    ~3 at GIS 85 -- the digitized dune line says 65.9 m, the DSAS shoreline
    record says 19.5 m. Neither can be preferred on its own terms. But they
    imply very different relocation dates, and the relocation dates are
    observed. So the timing test is the tie-breaker the two measurements cannot
    provide for each other.

THE HOLDOUT, AND WHY IT MATTERS
    There are TWO events, so an arm can be judged on one and tested on the
    other. `--holdout 1989` scores only the 1999 block, and vice versa. An N
    chosen because it reproduces 1989 has no claim on 1999, and that is the
    check worth having: fitting to both at once produces a better number and no
    way to know whether it means anything.

CENSORING
    A domain whose road never relocates inside 1984-2004 is RIGHT-CENSORED, not
    an error of +20 years. It is reported as ">2004" and excluded from the mean
    error, with the count stated -- averaging a censored value in would quietly
    reward an arm for never relocating anything.

USAGE
    python HAT_score_relocation_timing.py --arms pea1989base
    (the insert arms this compared -- blocksv4, blocksduneline, blocksdsas... --
     lost their run outputs on 2026-09-07; only unmodified topography is kept)
    python HAT_score_relocation_timing.py --holdout 1989
==============================================================================
"""

from __future__ import annotations

import argparse
import csv
import glob
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

BUFFER = 15
START, END = 1984, 2004
EVENTS = {**{d: 1999 for d in range(9, 15)},
          **{d: 1989 for d in range(84, 88)}}
OUT = REPO / "output" / "experiments" / "pea1989_crest"


def load(arm):
    hits = glob.glob(str(REPO / "output" / "raw_runs" / arm / "1984_2004"
                         / "calibBE" / "*" / "*.npz"))
    if not hits:
        raise SystemExit("\nno run for arm {!r}\n".format(arm))
    return np.load(hits[0], allow_pickle=True)["cascade"][0]


def first_relocation(c, gis):
    mgr = c.roadways[gis + BUFFER - 1]
    if mgr is None:
        return None, np.nan
    rel = np.asarray(getattr(mgr, "_road_relocated_TS", []), dtype=float)
    sb = np.asarray(getattr(mgr, "_road_setback_TS", []), dtype=float)
    t0 = float(sb[0]) if sb.size else np.nan
    hit = np.flatnonzero(np.nan_to_num(rel) > 0)
    return (START + int(hit[0]) if hit.size else None), t0


def main() -> None:
    ap = argparse.ArgumentParser()
    # blocksv4 was the default until 2026-09-07, when every insert arm's run
    # output was deleted; pea1989base(noreloc) is the one experiment arm left.
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
    # NAMED FOR THE ARMS SCORED. A fixed filename meant each run
    # silently replaced the previous arms' scores -- the
    # blocksduneline/dsas/minimum comparison was lost that way and had
    # to be re-derived from the runs.
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
