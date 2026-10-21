#!/usr/bin/env python3
"""
The pinned and drift groins rerun under the 2026-10-12 setup, so the figures compare like with like.

    python same_setup_runs.py [--streams 4]     # 4 runs, about 3 min each

The pinned groin (the runner's BlockingGroinCallback, b 0.6, f 0.6, instant failure at the 2004 step) and the
10-08 drift groin (b 0.1, f 0.4) in both periods, filed under groin-wave-driven/2026-10-12-same-setup/ so their 10-08 runs
are untouched. Each is launched exactly as its own study launches it (schedule_refit.py, drift_groin.py); only
the filing location differs. results_figures.py reads them through pinned_file() and drift_file().

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-12
"""
from __future__ import annotations

import argparse
import sys
import threading

import corrected_groin as cgm  # noqa: E402  (sets the import paths)
import drift_groin as dg  # noqa: E402
import schedule_refit as sr  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
STUDY = "groin-wave-driven/2026-10-12-same-setup"
PIN = (0.6, 0.6)                     # the runner's pinned groin
DRIFT = (0.1, 0.4)                   # the drift groin's 10-08 best cell
LOGS = cgm.LOGS / "same_setup"       # not the two studies' own logs/, which hold their 10-08 runs' logs
# -----------------------------------------------------------------------------
sr.LOGS = dg.LOGS = LOGS


def pinned_file(period):
    sr.STUDY = STUDY
    return sr.shoreline_file("instant2004", period, *PIN)


def drift_file(period):
    return dg.shoreline_file(f"{STUDY}/drift", period, *DRIFT)


def _pinned(period):
    sr.STUDY = STUDY
    return sr.run("instant2004", period, *PIN)


def _drift(period):
    dg.STUDY = f"{STUDY}/drift"
    return dg.run(period, *DRIFT)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--streams", type=int, default=4)
    ap.parse_args()
    jobs = [(_pinned, p) for p in (cgm.CALIBRATION, cgm.TEST) if pinned_file(p) is None]
    jobs += [(_drift, p) for p in (cgm.CALIBRATION, cgm.TEST) if drift_file(p) is None]
    codes = []
    threads = [threading.Thread(target=lambda f=f, p=p: codes.append(f(p))) for f, p in jobs]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    print(f"{len(jobs)} runs, failed: {sum(1 for c in codes if c)}")
    sys.exit(1 if any(codes) else 0)


if __name__ == "__main__":
    main()
