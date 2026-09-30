"""Dune line vs shoreline as the island offset, on the OLD ÷10 offset (2026-09-28).

Asked by Hannah on 2026-09-28: rebuild the offset-source experiment on the ÷10
island offset from before the metres fix, clearly labelled, "because I am
trying to track the changes as I was experimenting".

    offsets   the builds from before the metres fix, as the runs of the time
              read them: offset_mode "asrun" (offset / 10, the historical unit
              error) on superseded_20260924_pre-metres/v1 of each source.
              asrun divides the buffers too, so only the build a ÷10 run was
              made from reproduces it (1996/shoreline/v1/PROVENANCE.md)
    waves     the ÷10-era defaults, the calibration the model ran at before
              2026-09-24: Hs 2.5 m, Tp 8 s, asymmetry 0.7, high-angle 0.1
              (Hannah's choice). Offset AND waves differ from the 09-25 study,
              so the two cannot separate their effects
    ends      zeroBE; relocations and groins off, as on 09-25
    scope     natural and full management, 1996-2010; 4 runs
    Barrier3D the current one, with the route_overwash fix. The ÷10-era runs
              used the unfixed router, which segfaulted on some storms
    figures   the 09-25 driver's house form: metres of net change, observations
              LOWESS 7, no scores on the figures

WHERE: output/raw_runs/experiments/island-offset/2026-09-28-div10-offset-duneline-vs-shoreline-rebuild/

    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_div10.py run --jobs 4
    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_div10.py score
    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_div10.py plot

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_offset_source_comparison as base  # noqa: E402

TAG = "island-offset/2026-09-28-div10-offset-duneline-vs-shoreline-rebuild"
STUDY_DIR = base.grid.RAW_RUNS / "experiments" / TAG
PRE_METRES = "superseded_20260924_pre-metres/v1"
DIV10_WAVES = {"hs": 2.5, "wave_period_s": 8.0, "wave_asymmetry": 0.7}
HEADLINE = 0.1

# The 09-25 functions read these module globals at call time.
base.TAG, base.STUDY_DIR = TAG, STUDY_DIR
base.TABLES_DIR, base.LOGS_DIR, base.FIG = (STUDY_DIR / "tables", STUDY_DIR / "logs",
                                            STUDY_DIR / "figures")
base.BASE, base.HIGH_ANGLE, base.HEADLINE = DIV10_WAVES, (HEADLINE,), HEADLINE
base.WAVE_NOTE = "the ÷10-era defaults, before 2026-09-24"
base.OFFSET_NOTE = (f"÷10: offset_mode asrun (offset / 10) on {PRE_METRES} of each "
                    "source, as before the metres fix (2010 shoreline: its v1, the only "
                    "build, which differs from a pre-metres one only in the buffers)")
base.FORM_NOTE = ("A REBUILD (2026-09-28), not a record: the 09-25 design on the ÷10 offset; "
                  "compare ../2026-09-25-metres-offset-duneline-vs-shoreline-waves-hs1-tp8/ (metres, "
                  "09-25 waves).")
base.LEGEND_TITLE = "REBUILD on the ÷10 island offset (before the metres fix); no source/sink correction at the ends"

_env = base.env


def env(src, sc, s):
    """The 09-25 run environment, on the ÷10 offset."""
    e = _env(src, sc, s)
    e["HAT_OFFSET_MODE"] = "asrun"
    e["HAT_OFFSET_VERSION_1996"] = PRE_METRES
    e["HAT_OFFSET_VERSION_1996_SHORELINE"] = PRE_METRES
    # 2010: the dune line has a pre-metres build; the 2010 shoreline offset was
    # first built 2026-09-28, so asrun reads its v1. On the 1996 builds the two
    # differ only in the 15+15 buffer domains, not the real ones.
    e["HAT_OFFSET_VERSION_2010"] = PRE_METRES
    return e


base.env = env


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--jobs", type=int, default=4)
    sub.add_parser("score")
    sub.add_parser("plot")
    r = sub.add_parser("run-2010")
    r.add_argument("--jobs", type=int, default=2)
    a = ap.parse_args()
    return {"run": base.cmd_run, "score": base.cmd_score, "plot": base.cmd_plot,
            "run-2010": base.cmd_run_2010}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
