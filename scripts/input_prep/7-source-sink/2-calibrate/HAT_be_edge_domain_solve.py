# ==============================================================================
# HAT_be_edge_domain_solve.py
#
# What background-erosion rate should the two LOCKED END DOMAINS carry, for one
# hindcast period?
#
# WHAT THE END VALUES ARE FOR
#   GIS 1 and GIS 90 sit on the open boundaries of the modelled reach and
#   absorb the artefact there. They are boundary-artefact absorbers, not a
#   sediment budget -- see the end-domain note in hatteras_site_config.py,
#   which this script does not restate and must not contradict.
#
# WHY A SCRIPT
#   The solve is a Newton iteration: run, read the residual at the two ends,
#   step, run again. It was done by hand for 1984 and 2004, and the arithmetic
#   between runs -- which target column, which estimator, which secant -- was
#   carried in a person's head. Two periods were added on 2026-09-11 and both
#   need the same solve, so the arithmetic is written down here.
#
#   It does NOT run the model. It reads runs that have already happened and
#   prints the next probe, so every step stays a deliberate act.
#
# THE TWO TARGETS ARE DIFFERENT ESTIMATORS, DELIBERATELY
#   GIS 1  the raw per-domain transect mean. LoessConfig.skip_southern_domains
#          is 10, so D1-D10 are drawn raw rather than smoothed.
#   GIS 90 the LOESS-10 value, which is what is drawn everywhere north of D10.
#   That splice is what the rate-comparison figure draws, so fitting against
#   the same table means fit and figure cannot disagree. Both come out of
#   build_target_table, so neither is computed here.
#
# THE GAIN IS SMALL AND NOT CONSTANT
#   d(LRR)/d(BE) ran 0.092 to 0.123 across the four solved cases: only about a
#   tenth of an imposed edge rate survives in that domain's own shoreline, the
#   rest being diffused alongshore by BRIE. So each value is roughly ten times
#   the misfit it closes, and a single global slope should not be assumed --
#   which is why the second step uses the LOCAL secant through two real runs
#   rather than the nominal gain again.
#
# USAGE
#   One run -- report the residual and a first step at the nominal gain:
#       python 2-calibrate/HAT_be_edge_domain_solve.py --period 1996 --run <run_name>
#
#   Two or more -- local secant through the last two, and the next step:
#       python 2-calibrate/HAT_be_edge_domain_solve.py --period 1996 --run <first> --run <second>
#
#   Runs are named as they appear in run_index.csv. The rate each was run
#   under is read from that index, not retyped.
#
# Author: Hannah A. Henry, UNC CECL
# ==============================================================================

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents
                    if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT))

from hatteras_site_config import (                       # noqa: E402
    HATTERAS_PERIODS, HATTERAS_DOMAINS, HATTERAS_BE_EDGE_DOMAINS)
from cascade_pipeline.hindcast import build_target_table  # noqa: E402
from cascade_pipeline.coastsat_loess import (            # noqa: E402
    CoastSatDataset, LoessConfig, build_coastsat_series)
from hat_observed_rates import lrr_csv                   # noqa: E402

RUN_ROOT = PROJECT_ROOT / "output" / "raw_runs"
RUN_INDEX = RUN_ROOT / "run_index.csv"
# Resolved through hat_observed_rates, which owns the location.

# Section 8 of the runner builds the target this way. Kept identical rather
# than imported from it, because importing that file RUNS a hindcast.
LOESS_CONFIG = LoessConfig(window_domains=(10,), skip_southern_domains=10)
TARGET_WINDOW = 10

# The model side of the residual. Must be the OLS slope, matching the
# observed side; change_rate_m_yr is the endpoint difference and every preset
# fitted against it before 2026-08-22 is not reproducible from this pipeline.
RATE_COLUMN = "lrr_m_yr"

# d(LRR)/d(BE), used ONLY for the first step, when there is nothing to take a
# secant through. Mid-range of the four solved cases.
NOMINAL_GAIN = 0.105

# The index column holding the rate each run imposed, per end domain.
INDEX_RATE_COLUMN = {1: "be_rate_gis1_m_yr", 90: "be_rate_gis90_m_yr"}


def load_target(start_year, end_year):
    """The target table the runner grades against, as {gis: rate}."""
    csv_path = lrr_csv(start_year, end_year)   # raises, listing windows
    series = build_coastsat_series(
        [CoastSatDataset(label="CoastSat {0}".format(start_year),
                         period_start=start_year, csv_path=str(csv_path))],
        start_year, LOESS_CONFIG, domains=HATTERAS_DOMAINS)
    table = build_target_table(series[0], LOESS_CONFIG, HATTERAS_DOMAINS,
                               TARGET_WINDOW)
    return dict(zip([int(g) for g in table["gis_domain"]],
                    [float(r) for r in table["target_lrr_m_yr"]]))


def find_run(run_name, start_year, end_year):
    """The rate CSV of one run, located under the period's own directory."""
    period_dir = RUN_ROOT / "{0}_{1}".format(start_year, end_year)
    hits = sorted(period_dir.glob("*/{0}/tables/shoreline_change_rate.csv"
                                  .format(run_name)))
    if not hits:
        raise FileNotFoundError(
            "no rate CSV for run {0!r} under {1}. Check the name against "
            "run_index.csv.".format(run_name, period_dir))
    return hits[0]


def imposed_rates(run_name):
    """What this run imposed at each end domain, from run_index.csv."""
    index = pd.read_csv(RUN_INDEX)
    rows = index[index["run_name"] == run_name]
    if rows.empty:
        raise ValueError("{0!r} is not in run_index.csv".format(run_name))
    row = rows.iloc[-1]
    return {gis: float(row[column])
            for gis, column in INDEX_RATE_COLUMN.items()}


def read_model(csv_path):
    frame = pd.read_csv(csv_path)
    if RATE_COLUMN not in frame.columns:
        raise KeyError(
            "{0} has no {1!r} column -- it predates the LRR estimator. "
            "Re-run, or backfill with HAT_backfill_run_lrr.py.".format(
                csv_path.name, RATE_COLUMN))
    return frame.set_index("gis_domain")[RATE_COLUMN].to_dict()


def report(period, runs):
    start_year = period
    end_year = HATTERAS_PERIODS[period]["end_year"]
    target = load_target(start_year, end_year)

    states = []
    for run_name in runs:
        model = read_model(find_run(run_name, start_year, end_year))
        states.append({"run": run_name,
                       "imposed": imposed_rates(run_name),
                       "model": model})

    print("\n" + "=" * 74)
    print("END-DOMAIN SOLVE   {0}-{1}".format(start_year, end_year))
    print("=" * 74)

    suggestion = {}
    for gis in HATTERAS_BE_EDGE_DOMAINS:
        want = target[gis]
        print("\nGIS {0}   target {1:+.3f} m/yr".format(gis, want))
        print("  {0:<46} {1:>9} {2:>9} {3:>9}".format(
            "run", "imposed", "model", "residual"))
        points = []
        for state in states:
            imposed = state["imposed"][gis]
            got = state["model"][gis]
            residual = got - want
            points.append((imposed, got))
            print("  {0:<46} {1:>+9.2f} {2:>+9.3f} {3:>+9.3f}".format(
                state["run"][:46], imposed, got, residual))

        last_imposed, last_model = points[-1]
        residual = last_model - want

        if len(points) >= 2:
            (x0, y0), (x1, y1) = points[-2], points[-1]
            if x1 == x0:
                print("  the last two runs imposed the same rate; "
                      "no secant to take")
                continue
            gain = (y1 - y0) / (x1 - x0)
            source = "local secant through the last two runs"
        else:
            gain = NOMINAL_GAIN
            source = "nominal gain, no second run to take a secant through"

        step = -residual / gain
        suggestion[gis] = last_imposed + step
        print("  d(LRR)/d(BE) {0:.4f}   ({1})".format(gain, source))
        print("  next         {0:+.1f} m/yr   "
              "(step {1:+.1f} to close {2:+.3f})".format(
                  suggestion[gis], step, residual))
        if abs(residual) < 0.02:
            print("  CONVERGED at this tolerance; another step is noise")

    if suggestion:
        override = ",".join("{0}={1:.1f}".format(gis, rate)
                            for gis, rate in sorted(suggestion.items()))
        print("\nnext probe:")
        print('  HAT_BE_OVERRIDE="{0}"'.format(override))
    print()


def main():
    parser = argparse.ArgumentParser(
        description="solve the two locked end domains for one period")
    parser.add_argument("--period", type=int, required=True,
                        help="period start year")
    parser.add_argument("--run", action="append", required=True,
                        help="run name, oldest first; repeatable")
    args = parser.parse_args()

    if args.period not in HATTERAS_PERIODS:
        parser.error("no such period {0}; have {1}".format(
            args.period, sorted(HATTERAS_PERIODS)))
    report(args.period, args.run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
