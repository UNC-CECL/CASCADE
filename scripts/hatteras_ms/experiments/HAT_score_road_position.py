#!/usr/bin/env python3
r"""
HAT_score_road_position.py
==============================================================================
Scores each arm's MODELLED 2004 road setback against the MEASURED 2004 setback.

WHY THIS EXISTS -- the timing test cannot do it
    HAT_score_relocation_timing.py scores the predicted relocation YEAR. It
    ranks the DSAS-derived N far above the dune-line-derived one (MAE 1.2 yr
    against 11.0). That ranking is not trustworthy on its own, because the
    relocation year is a function of TWO unknowns that trade off exactly:

        small initial setback + trigger at zero
        large initial setback + trigger at a maintenance buffer

    Both reproduce 1989. One observation cannot resolve two unknowns, so a good
    timing score is evidence about the PAIR, not about N.

    The road's POSITION breaks the degeneracy. Where the road physically ended
    up by 2004 is measured -- RoadOffset_2004_domains.csv, an independent
    same-year measurement on the 2004-start topography -- and the two N
    estimates predict different positions regardless of when the move happened.

REQUIRES THE PRESCRIBED RELOCATIONS ON
    The real NC-12 was moved by NCDOT in 1989 and 1999. A model run with the
    events off has not been given those moves, so its 2004 position answers a
    different question. Score arms built with --relocations 1.

WHAT A GOOD SCORE DOES AND DOES NOT MEAN
    Agreement here says the modelled road ends the period where the real one
    did. It does NOT validate the relocation year, the dune history, or the
    fabricated land -- those need their own observations. It is one number
    against one measurement, which is exactly why it is worth having alongside
    the timing test rather than instead of it.

USAGE
    python HAT_score_road_position.py --arms pea1989base
    (the insert arms this compared -- blocksv4, blocksduneline, blocksdsas... --
     lost their run outputs on 2026-09-07; only unmodified topography is kept)
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
# Anchored by SEARCHING UPWARD for the project root rather than by
# counting parent directories (2026-09-13). A counted depth is correct
# only while the file stays where it was written, and these moved into
# subfolders of hatteras_ms. Six files here already did it this way.
REPO = next(_p for _p in HERE.parents if (_p / 'pyproject.toml').exists())
sys.path.insert(0, str(REPO / "scripts"))

BUFFER = 15
OUT = REPO / "output" / "experiments" / "pea1989_crest"


def load(arm):
    hits = glob.glob(str(REPO / "output" / "raw_runs" / arm / "1984_2004"
                         / "calibBE" / "*" / "*.npz"))
    if not hits:
        raise SystemExit("\nno run for arm {!r}\n".format(arm))
    return np.load(hits[0], allow_pickle=True)["cascade"][0]


def main() -> None:
    ap = argparse.ArgumentParser()
    # blocksv4 was the default until 2026-09-07, when every insert arm's run
    # output was deleted; pea1989base is the one experiment arm left.
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
