#!/usr/bin/env python3
"""
Split-sample validation of the groin's strength: fit on one hindcast window, test unchanged on the other.

    python split_sample.py run   [--streams 4]   # 8 runs through the unchanged runner
    python split_sample.py score                 # -> split_sample_scores.csv, printed table

The window strengths (0.16 for 1996-2009, 0.10 for 2009-2025) were each fitted on their own window
(../2026-10-14-window-refit/), so neither window is an independent test of the groin. Here each window is also run
at the other window's strength, and at no groin, on today's runner with the setup of the window refit (edgeBE,
full_management, no relocation, off 2017). The whole-life 0.25 is run too, for context: it was fitted on
1967-2025, which contains both windows, so it is not a test. Running and scoring are window_refit.py's, pointed
at this study. Skill is against the no-groin run: 1 - MSE / MSE(no groin), on the annual gap across GIS 5|6.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-15
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-10-14-window-refit"))
import window_refit as wr  # noqa: E402  (its run and score, pointed at this study below)

# --- CONFIG ------------------------------------------------------------------
FIT = {1996: 0.16, 2009: 0.10}            # each window's own fit (2026-10-14-window-refit)
WHOLE_LIFE = 0.25                         # 1967-2025 fit (2026-10-13-pregroin-planform); not independent
OTHER = {1996: 2009, 2009: 1996}
wr.STUDY = "groin-wave-driven/2026-10-15-split-sample"
wr.RAW = wr.REPO / "output" / "raw_runs" / "experiments" / wr.STUDY
wr.LOGS = HERE / "logs"
wr.STRENGTHS = {p: (None, FIT[p], FIT[OTHER[p]], WHOLE_LIFE) for p in FIT}
# -----------------------------------------------------------------------------


# The role of a strength in a window's test
def role(period, s):
    if s is None:
        return "no groin"
    if s == FIT[period]:
        return "fitted on this window"
    if s == FIT[OTHER[period]]:
        return f"fitted on {wr.WINDOW[OTHER[period]]} (independent test)"
    return "whole life, 1967-2025 (not independent)"


def cmd_score(_a):
    rows = []
    for period, ss in wr.STRENGTHS.items():
        obs = wr.cg.observed_change(period)
        for s in ss:
            p = wr.shoreline_file(period, s)
            if p is None:
                print(f"missing: {wr.WINDOW[period]} {wr.member(s)}")
                continue
            rows.append(dict(window=wr.WINDOW[period], strength=0.0 if s is None else s, role=role(period, s),
                             **wr.cg.score_annual(p, period, obs)))
    t = pd.DataFrame(rows)
    base = t[t.role == "no groin"].set_index("window").rmse_annual
    t["skill_vs_no_groin"] = 1 - (t.rmse_annual / t.window.map(base)) ** 2
    fit = t[t.role == "fitted on this window"].set_index("window").rmse_annual
    # Share of the improvement the window's own fit makes over no groin that the transferred strength keeps
    t["share_of_fit_improvement"] = (t.window.map(base) - t.rmse_annual) / (t.window.map(base) - t.window.map(fit))
    t.to_csv(HERE / "split_sample_scores.csv", index=False, float_format="%.3f")
    pd.set_option("display.width", 200)
    print(t[["window", "strength", "role", "rmse_annual", "bias_annual", "skill_vs_no_groin",
             "share_of_fit_improvement"]].round(3).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--streams", type=int, default=4)
    sub.add_parser("score")
    a = ap.parse_args()
    {"run": wr.cmd_run, "score": cmd_score}[a.cmd](a)


if __name__ == "__main__":
    main()
