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
# AN EXTENDED GEOMETRY (2026-09-16, the Pea Island extension experiment)
#   The ends are wherever HATTERAS_DOMAINS puts them -- GIS 1 and 115 under
#   HAT_GEOMETRY=n115 -- and the target for a
#   domain beyond GIS 90 comes from the window's extension rate table
#   (coastsat_lrr/<window>/ext/transect_lrr_with_base.csv, the surveyed
#   transects plus the extension's, one LOESS over the whole reach), which is
#   exactly the table the runner grades that geometry's ends against. Run
#   this script with the SAME HAT_GEOMETRY as the runs it reads. A run that
#   imposed nothing at an end (a zeroBE probe) reads as 0.0 there.
#
# THE DUNE LINE AS THE TARGET (2026-09-16, Hannah: "what if we used the dune
# line change" to set the ends)
#   --target duneline reads the observation from the two digitised dune lines
#   of the window instead of CoastSat: the end vintage's line minus the start
#   vintage's, per domain, over the survey interval (coastsat_vs_duneline
#   .KNOWN_SURVEY_DATES; a missing date is mid-year), seaward positive --
#   exactly what HAT_rate_windows.py draws under vs_duneline/endpoint_net_change. Two readings of
#   it at an end domain, --dune-smooth raw (the domain's own value) and mean3
#   (the mean of it and its two inward neighbours, GIS 1-3 / 88-90). A dune
#   line is two surveys, so --estimator endpoint reads change_rate_m_yr on
#   the model side; the default lrr keeps the CoastSat protocol. The runs of
#   that solve are under output/raw_runs/experiments/2026-09-16-dune-edgesolve/.
#
# USAGE
#   One run -- report the residual and a first step at the nominal gain:
#       python 2-calibrate/HAT_be_edge_domain_solve.py --period 1996 --run <run_name>
#
#   Two or more -- local secant through the last two, and the next step:
#       python 2-calibrate/HAT_be_edge_domain_solve.py --period 1996 --run <first> --run <second>
#
#   Runs are named as they appear in run_index.csv and located through
#   run_registry.find_run_dir: --kind and --tag say where (the matrix by
#   default; an experiment's runs need --kind experiment --tag <set>/<member>,
#   one --tag per --run when the members differ). The rate each was run
#   under is read from the index, not retyped.
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

from site_layer.hatteras_site_config import (                       # noqa: E402
    HATTERAS_PERIODS, HATTERAS_DOMAINS, HATTERAS_BE_EDGE_DOMAINS,
    HATTERAS_GEOMETRY, HATTERAS_GEOMETRY_EXTENDED)
from cascade_pipeline.hindcast import build_target_table  # noqa: E402
from cascade_pipeline.coastsat_loess import (            # noqa: E402
    CoastSatDataset, LoessConfig, build_coastsat_series)
from cascade_pipeline.run_layout import resolve as resolve_run_file  # noqa: E402
from cascade_pipeline.run_registry import (              # noqa: E402
    MATRIX_KIND, find_run_dir, load_run_index)
from site_layer.hat_observed_rates import lrr_csv, lrr_csv_ext      # noqa: E402

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

# The index column holding the rate each run imposed, per end domain. The
# runner writes be_rate_gis<N>_m_yr for each of HATTERAS_BE_EDGE_DOMAINS.
INDEX_RATE_COLUMN = {gis: "be_rate_gis{0}_m_yr".format(gis)
                     for gis in HATTERAS_BE_EDGE_DOMAINS}


def load_target(start_year, end_year, window=None):
    """The target table the runner grades against, as {gis: rate}. `window`
    (start, end) grades against another CoastSat LRR window instead, e.g.
    the full 1996-2024 rate for a 1996-2010 run (2026-09-19, Hannah)."""
    if window is not None:
        start_year, end_year = window
    # raises, listing windows, if absent; the extension table if the
    # geometry reaches beyond GIS 90, as the runner's section 8 does
    csv_path = (lrr_csv_ext(start_year, end_year) if HATTERAS_GEOMETRY_EXTENDED
                else lrr_csv(start_year, end_year))
    series = build_coastsat_series(
        [CoastSatDataset(label="CoastSat {0}".format(start_year),
                         period_start=start_year, csv_path=str(csv_path))],
        start_year, LOESS_CONFIG, domains=HATTERAS_DOMAINS)
    table = build_target_table(series[0], LOESS_CONFIG, HATTERAS_DOMAINS,
                               TARGET_WINDOW)
    return dict(zip([int(g) for g in table["gis_domain"]],
                    [float(r) for r in table["target_lrr_m_yr"]]))


def load_dune_target(start_year, end_year, smooth):
    """The dune-line endpoint rate per domain, seaward positive, read at each
    end domain raw or as a three-domain mean. Returns ({gis: rate}, note).

    READ FROM the stored product 5-scr/3-rates/duneline/endpoint/<window>/
    (2026-09-18), the same numbers HAT_rate_windows.py draws, rather than
    recomputed here from the raw offsets."""
    from site_layer.hat_observed_rates import dune_endpoint_csv
    dom = pd.read_csv(dune_endpoint_csv(start_year, end_year, "domain"))
    meta = pd.read_csv(dune_endpoint_csv(start_year, end_year, "transect")).iloc[0]
    rate = dom.set_index("domain_number")["mean_rate_m_yr"]
    v0, v1 = int(meta["start_vintage"]), int(meta["end_vintage"])
    d0, d1 = pd.Timestamp(meta["start_date"]), pd.Timestamp(meta["end_date"])
    a0, a1 = bool(meta["start_date_assumed"]), bool(meta["end_date_assumed"])
    years = float(meta["interval_yr"])
    first, last = int(rate.index.min()), int(rate.index.max())
    out = {}
    for gis in HATTERAS_BE_EDGE_DOMAINS:
        if smooth == "raw":
            out[gis] = float(rate.loc[gis])
        elif smooth == "mean3":
            # the end domain and its two INWARD neighbours
            sel = rate.loc[gis - 2:gis] if gis == last else rate.loc[gis:gis + 2]
            out[gis] = float(sel.mean())
        else:
            raise ValueError(smooth)
    note = ("dune line {0} ({1}{2}) to {3} ({4}{5}), {6:.2f} yr, {7}".format(
        v0, d0.date(), " assumed" if a0 else "", v1, d1.date(),
        " assumed" if a1 else "", years,
        "raw domain value" if smooth == "raw" else "three-domain mean"))
    return out, note


def find_run(run_name, start_year, end_year, preset, kind, tag):
    """The rate CSV of one run, located the way every other reader does."""
    run_dir = find_run_dir(RUN_ROOT, run_name, (start_year, end_year), preset,
                           kind=kind, tag=tag)
    return Path(resolve_run_file(run_dir, "rate_csv", run_name))


def index_row(run_name, kind, tag):
    """The run's row in run_index.csv, by its full identity."""
    index = load_run_index(RUN_INDEX)
    rows = index[(index["run_name"] == run_name) & (index["kind"] == kind)
                 & (index["tag"] == tag)]
    if rows.empty:
        raise ValueError("{0!r} (kind {1}, tag {2!r}) is not in run_index.csv"
                         .format(run_name, kind, tag))
    return rows.iloc[-1]


def imposed_rates(row):
    """What this run imposed at each end domain, from its index row.

    A run from a preset that names no rate at an end (zeroBE, or an edgeBE
    whose column predates this end) reads as 0.0 there, which is what it
    imposed.
    """
    out = {}
    for gis, column in INDEX_RATE_COLUMN.items():
        value = row.get(column, "")
        out[gis] = float(value) if value not in ("", None) else 0.0
    return out


ESTIMATOR_COLUMN = {"lrr": RATE_COLUMN, "endpoint": "change_rate_m_yr"}


def read_model(csv_path, column=RATE_COLUMN):
    frame = pd.read_csv(csv_path)
    if column not in frame.columns:
        raise KeyError(
            "{0} has no {1!r} column -- it predates the LRR estimator. "
            "Re-run, or backfill with HAT_backfill_run_lrr.py.".format(
                csv_path.name, column))
    return frame.set_index("gis_domain")[column].to_dict()


def report(period, runs, preset, kinds, tags, target_source="coastsat",
           dune_smooth="raw", estimator="lrr", coastsat_window=None):
    start_year = period
    end_year = HATTERAS_PERIODS[period]["end_year"]
    if target_source == "coastsat":
        target = load_target(start_year, end_year, coastsat_window)
        target_note = "CoastSat LRR{0}: raw mean at GIS 1, LOESS-10 at GIS 90".format(
            " {0}-{1}".format(*coastsat_window) if coastsat_window else "")
    else:
        target, target_note = load_dune_target(start_year, end_year, dune_smooth)
    column = ESTIMATOR_COLUMN[estimator]

    states = []
    for run_name, kind, tag in zip(runs, kinds, tags):
        row = index_row(run_name, kind, tag)
        # The preset folder the run sits under is the preset it ran, which
        # the index knows: a zeroBE stage-0 run and the edgeBE probes after
        # it belong to one solve and are read together.
        folder = preset or str(row["source_sink_preset"])
        model = read_model(find_run(run_name, start_year, end_year, folder, kind, tag),
                           column)
        states.append({"run": run_name, "tag": tag,
                       "imposed": imposed_rates(row),
                       "model": model})

    print("\n" + "=" * 74)
    print("END-DOMAIN SOLVE   {0}-{1}   geometry {2} (ends GIS {3} and {4})".format(
        start_year, end_year, HATTERAS_GEOMETRY, *HATTERAS_BE_EDGE_DOMAINS))
    print("target   {0}".format(target_note))
    print("model    {0} ({1})".format(estimator, column))
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
            label = state["run"] if len(set(tags)) == 1 else "{0} [{1}]".format(
                state["run"], state["tag"].rsplit("/", 1)[-1])
            print("  {0:<46} {1:>+9.2f} {2:>+9.3f} {3:>+9.3f}".format(
                label[:46], imposed, got, residual))

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
    return suggestion


def main():
    parser = argparse.ArgumentParser(
        description="solve the two locked end domains for one period")
    parser.add_argument("--period", type=int, required=True,
                        help="period start year")
    parser.add_argument("--run", action="append", required=True,
                        help="run name, oldest first; repeatable")
    parser.add_argument("--preset", default=None,
                        help="the preset folder the runs sit under; by default "
                             "each run's own preset, read from the index")
    parser.add_argument("--kind", action="append", default=None,
                        help="matrix (default), experiment, sensitivity, "
                             "version; one for all runs, or one per --run")
    parser.add_argument("--tag", action="append", default=None,
                        help="the run's tag; one for all runs, or one per --run")
    parser.add_argument("--target", choices=("coastsat", "duneline"),
                        default="coastsat",
                        help="the observation the ends are solved against")
    parser.add_argument("--dune-smooth", choices=("raw", "mean3"), default="raw",
                        help="with --target duneline: the end domain's own "
                             "value, or the mean of it and its two neighbours")
    parser.add_argument("--estimator", choices=("lrr", "endpoint"), default="lrr",
                        help="model column: lrr_m_yr (the OLS slope, default) "
                             "or change_rate_m_yr (endpoint over run years)")
    parser.add_argument("--coastsat-window", default=None, metavar="START_END",
                        help="with --target coastsat: grade against this CoastSat "
                             "LRR window instead of the run's own, e.g. 1996_2024")
    args = parser.parse_args()

    if args.period not in HATTERAS_PERIODS:
        parser.error("no such period {0}; have {1}".format(
            args.period, sorted(HATTERAS_PERIODS)))

    def per_run(values, default, what):
        values = values or [default]
        if len(values) == 1:
            values = values * len(args.run)
        if len(values) != len(args.run):
            parser.error("give one --{0}, or one per --run".format(what))
        return values

    tags = per_run(args.tag, "", "tag")
    kinds = per_run(args.kind, MATRIX_KIND, "kind")
    report(args.period, args.run, args.preset, kinds, tags,
           target_source=args.target, dune_smooth=args.dune_smooth,
           estimator=args.estimator,
           coastsat_window=(tuple(int(x) for x in args.coastsat_window.split("_"))
                            if args.coastsat_window else None))
    return 0


if __name__ == "__main__":
    sys.exit(main())
