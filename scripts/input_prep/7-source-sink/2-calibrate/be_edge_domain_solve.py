"""
Solve the background-erosion rate the two locked end domains (GIS 1 and 90) carry, for one period.

    python scripts/input_prep/7-source-sink/2-calibrate/be_edge_domain_solve.py --period 1996 --run <run_name>
    python scripts/input_prep/7-source-sink/2-calibrate/be_edge_domain_solve.py --period 1996 --run <first> --run <second>

One run gives the residual and a first step at the nominal gain; two or more
a secant step. Targets the CoastSat LRR, the dune line (--target duneline), or the
DEM-to-DEM net change in metres (--target net_change). Details: scripts/input_prep/7-source-sink/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-05
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
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
from cascade_pipeline.coastsat_lowess import (            # noqa: E402
    CoastSatDataset, LowessConfig, build_coastsat_series)
from cascade_pipeline.run_layout import resolve as resolve_run_file  # noqa: E402
from cascade_pipeline.run_registry import (              # noqa: E402
    MATRIX_KIND, find_run_dir, load_run_index)
from site_layer.hat_observed_rates import (                         # noqa: E402
    NET_CHANGE_WINDOWS, lrr_csv, lrr_csv_ext, net_change_domain_csv)

# --- CONFIG ------------------------------------------------------------------
RUN_ROOT = PROJECT_ROOT / "output" / "raw_runs"
RUN_INDEX = RUN_ROOT / "run_index.csv"
# Resolved through hat_observed_rates, which owns the location.

# Section 8 of the runner builds the target this way
LOWESS_CONFIG = LowessConfig(window_domains=(7,), skip_southern_domains=10)
TARGET_WINDOW = 7   # 10 until 2026-09-28, with the runner

# The model side of the residual
RATE_COLUMN = "lrr_m_yr"

# d(LRR)/d(BE) for the first step only, before a secant exists (mid-range of four solves)
NOMINAL_GAIN = 0.105

# Convergence on a rate; a net-change solve scales it by the run length
CONVERGED_M_YR = 0.02

# The index column holding the rate each run imposed, per end domain
INDEX_RATE_COLUMN = {gis: "be_rate_gis{0}_m_yr".format(gis)
                     for gis in HATTERAS_BE_EDGE_DOMAINS}
# -----------------------------------------------------------------------------


# The target table the runner grades against, as {gis
def load_target(start_year, end_year, window=None):
    if window is not None:
        start_year, end_year = window
    # Raises, listing windows, if absent
    csv_path = (lrr_csv_ext(start_year, end_year) if HATTERAS_GEOMETRY_EXTENDED
                else lrr_csv(start_year, end_year))
    series = build_coastsat_series(
        [CoastSatDataset(label="CoastSat {0}".format(start_year),
                         period_start=start_year, csv_path=str(csv_path))],
        start_year, LOWESS_CONFIG, domains=HATTERAS_DOMAINS)
    table = build_target_table(series[0], LOWESS_CONFIG, HATTERAS_DOMAINS,
                               TARGET_WINDOW)
    return dict(zip([int(g) for g in table["gis_domain"]],
                    [float(r) for r in table["target_lrr_m_yr"]]))


# The dune-line endpoint rate per domain (seaward positive), read raw or as a three-domain end mean
def load_dune_target(start_year, end_year, smooth):
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


# The observed net change at the ends, as the runner smooths a target: raw at GIS 1, LOWESS-7 at GIS 90
def load_net_change_target(start_year, end_year):
    if (start_year, end_year) not in NET_CHANGE_WINDOWS:
        raise SystemExit("no net-change target for {0}-{1}; have {2}".format(
            start_year, end_year, sorted(NET_CHANGE_WINDOWS)))
    dom = pd.read_csv(net_change_domain_csv(start_year, end_year), index_col=0)
    first, last = HATTERAS_BE_EDGE_DOMAINS
    out = {first: float(dom.loc[first, "net_change_m"]),
           last: float(dom.loc[last, "net_change_lowess7_m"])}
    (s0, s1), (e0, e1) = NET_CHANGE_WINDOWS[(start_year, end_year)]
    note = ("CoastSat net change, mean {0} to {1} minus mean {2} to {3}: raw at GIS {4}, "
            "LOWESS-{5} at GIS {6}".format(e0, e1, s0, s1, first, TARGET_WINDOW, last))
    return out, note


# A run's end-minus-start shoreline per GIS domain, seaward + (x_s grows landward)
def read_model_net_change(run_dir):
    m = np.load(next(Path(run_dir).glob("*_shoreline_matrix.npy")))
    D = HATTERAS_DOMAINS
    change = -(m[-1] - m[0])[D.start_real_index:D.end_real_index]
    return dict(zip(range(D.first_gis_id, D.last_gis_id + 1), map(float, change)))


# The rate CSV of one run, located the way every other reader does
def find_run(run_name, start_year, end_year, preset, kind, tag):
    run_dir = find_run_dir(RUN_ROOT, run_name, (start_year, end_year), preset,
                           kind=kind, tag=tag)
    return Path(resolve_run_file(run_dir, "rate_csv", run_name))


# The run's row in run_index.csv, by its full identity
def index_row(run_name, kind, tag):
    index = load_run_index(RUN_INDEX)
    rows = index[(index["run_name"] == run_name) & (index["kind"] == kind)
                 & (index["tag"] == tag)]
    if rows.empty:
        raise ValueError("{0!r} (kind {1}, tag {2!r}) is not in run_index.csv"
                         .format(run_name, kind, tag))
    return rows.iloc[-1]


# What this run imposed at each end domain, from its index row
def imposed_rates(row):
    out = {}
    for gis, column in INDEX_RATE_COLUMN.items():
        value = row.get(column, "")
        out[gis] = float(value) if value not in ("", None) else 0.0
    return out


ESTIMATOR_COLUMN = {"lrr": RATE_COLUMN, "endpoint": "change_rate_m_yr"}


# One run's modelled rate column per domain
def read_model(csv_path, column=RATE_COLUMN):
    frame = pd.read_csv(csv_path)
    if column not in frame.columns:
        raise KeyError(
            "{0} has no {1!r} column -- it predates the LRR estimator. "
            "Re-run, or backfill with backfill_run_lrr.py.".format(
                csv_path.name, column))
    return frame.set_index("gis_domain")[column].to_dict()


# The ends' residuals and the next step for a set of runs
def report(period, runs, preset, kinds, tags, target_source="coastsat",
           dune_smooth="raw", estimator="lrr", coastsat_window=None):
    start_year = period
    end_year = HATTERAS_PERIODS[period]["end_year"]
    net = target_source == "net_change"
    unit = "m" if net else "m/yr"
    years = HATTERAS_PERIODS[period]["last_model_year"] - start_year + 1
    tolerance = CONVERGED_M_YR * (years if net else 1)
    if net:
        target, target_note = load_net_change_target(start_year, end_year)
        estimator = "net_change"
    elif target_source == "coastsat":
        target = load_target(start_year, end_year, coastsat_window)
        target_note = "CoastSat LRR{0}: raw mean at GIS 1, LOWESS-{1} at GIS 90".format(
            " {0}-{1}".format(*coastsat_window) if coastsat_window else "", TARGET_WINDOW)
    else:
        target, target_note = load_dune_target(start_year, end_year, dune_smooth)
    column = ESTIMATOR_COLUMN.get(estimator, "end minus start shoreline, m")

    states = []
    for run_name, kind, tag in zip(runs, kinds, tags):
        row = index_row(run_name, kind, tag)
        # The preset folder the run sits under is the preset it ran, which the index knows
        folder = preset or str(row["source_sink_preset"])
        model = (read_model_net_change(find_run_dir(RUN_ROOT, run_name, (start_year, end_year),
                                                    folder, kind=kind, tag=tag))
                 if net else read_model(find_run(run_name, start_year, end_year, folder,
                                                 kind, tag), column))
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
        print("\nGIS {0}   target {1:+.3f} {2}".format(gis, want, unit))
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
            gain = NOMINAL_GAIN * (years if net else 1)
            source = "nominal gain, no second run to take a secant through"

        step = -residual / gain
        suggestion[gis] = last_imposed + step
        print("  d({0})/d(BE) {1:.4f}   ({2})".format("net" if net else "LRR", gain, source))
        print("  next         {0:+.1f} m/yr   "
              "(step {1:+.1f} to close {2:+.3f})".format(
                  suggestion[gis], step, residual))
        if abs(residual) < tolerance:
            print("  CONVERGED at this tolerance; another step is noise")

    if suggestion:
        override = ",".join("{0}={1:.1f}".format(gis, rate)
                            for gis, rate in sorted(suggestion.items()))
        print("\nnext probe:")
        print('  HAT_BE_OVERRIDE="{0}"'.format(override))
    print()
    return suggestion


# Run: read the runs, report the residual and the next step
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
    parser.add_argument("--target", choices=("coastsat", "duneline", "net_change"),
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
