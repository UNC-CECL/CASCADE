#!/usr/bin/env python3
"""Solves the D2 / D12 edge background-erosion correction for the 1967 rig.

WHY THIS EXISTS
    `HAT_groin_hindcast_1967_2017.py` carries

        EDGE_BE_RATES_GIS = {2: 5.0, 12: 10.0}   # "(solve for this)"

    -- placeholders that were never solved. They are not a minor detail here.
    The reduced reach is 11 real domains inside 41, so D2 sits THREE domains
    from the groin pair at D5/D6: an imposed rate at the edge diffuses into the
    pair within a few years and lands directly on the signal the sweep is
    fitting. The observed edge rates over 1967-2023 are +1.18 and +1.35 m/yr,
    so the placeholders are 4-7x too strong.

WHAT IS SOLVED, AND AGAINST WHAT
    The same method the production hindcast uses for GIS 1 / GIS 90 -- the site
    config records those as "solved on the edgeBE road_bdm base run". A trial
    edge rate is imposed, the base run is stepped, and the modelled change at
    that domain is compared against the surveyed change from the fixed 1967
    datum. The rate is then adjusted and the run repeated.

    GROIN OFF, ALWAYS. If the edge were solved with the groin attached, the
    correction could absorb groin signal and the sweep would afterwards be
    fitting M against a background that had already eaten part of its effect.
    The whole point of an edge correction is that it is STRUCTURAL -- a fix for
    the buffer's artificial orientation -- so it must be solved with the
    structure of interest switched off.

    Both edges are solved together rather than one at a time. They are nearly
    independent (opposite ends of the reach) but not exactly, since each one's
    signal diffuses inward, so a joint secant step converges without pretending
    the coupling is zero.

METHOD
    Secant iteration on each edge independently. The shoreline response to an
    imposed background rate is close to linear over this range, so two
    evaluations bracket the root and each further step refines it. Converges in
    2-3 iterations; the cap exists so a non-convergent case stops rather than
    spinning.

Usage:
    python HAT_solve_edge_be.py [--max-iter 4] [--tol 1.0]

Prints the solved rates and writes them to edge_be_solved.json for the sweep
to read. Nothing else may run while this does -- every CASCADE construction
writes the shared Hatteras-CASCADE-parameters.yaml, and concurrent writers
corrupt it.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys

import numpy as np
import pandas as pd

sys.argv = [sys.argv[0]]          # hc parses nothing, but be explicit
import HAT_groin_hindcast_1967_2017 as hc   # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
SOLVED_JSON = HERE / "edge_be_solved.json"

EDGE_GIS = (2, 12)
DATUM_YEAR = 1967
FIT_YEAR = 2023                   # last surveyed wet/dry year

WETDRY_CHANGE_TABLE = (
    hc.PROJECT_BASE_DIR / "hard-structures" / "groin" / "HAT-groin-test-output"
    / "shoreline_position_output" / "Change_from_wetdry_1967_D2_D12.csv"
    if isinstance(hc.PROJECT_BASE_DIR, pathlib.Path)
    else pathlib.Path(hc.PROJECT_BASE_DIR) / "hard-structures" / "groin"
    / "HAT-groin-test-output" / "shoreline_position_output"
    / "Change_from_wetdry_1967_D2_D12.csv")

STORM_FILE_1967_2024 = (
    pathlib.Path(hc.PROJECT_BASE_DIR) / "hard-structures" / "groin"
    / "HAT-groin-test-input" / "groin_init" / "storms" / "1967_2024"
    / "1967_2024_groin_storms.npy")


def observed_edge_change():
    """Surveyed change at each edge domain, 1967 -> FIT_YEAR, in metres.

    Returns:
        {gis: change_m}, landward-positive.
    """
    frame = pd.read_csv(WETDRY_CHANGE_TABLE).set_index("Domain_ID")
    # The SECOND year in the name is the survey year; the first is the 1967
    # datum. Matching the first one silently returns the datum for every column.
    column = None
    for candidate in frame.columns:
        match = re.match(r"change_from_wetdry_1967_wetdry_(\d{4})", candidate)
        if match and int(match.group(1)) == FIT_YEAR:
            column = candidate
            break
    if column is None:
        raise ValueError(f"no wetdry column for {FIT_YEAR} in {WETDRY_CHANGE_TABLE.name}")
    return {gis: float(frame.loc[gis, column]) for gis in EDGE_GIS}


def modelled_edge_change(edge_rates):
    """Runs the no-groin base at `edge_rates` and returns its edge changes.

    Args:
        edge_rates: {gis: m/yr} imposed at the edge domains.

    Returns:
        {gis: modelled change in metres over the run}, landward-positive.
    """
    hc.MAKE_FIGURES = False
    hc.MAKE_RUN_GIF = False
    hc.RUN_MATRIX = ["no_groin"]
    hc.APPLY_EDGE_BE_CORRECTION = True
    hc.EDGE_BE_RATES_GIS = dict(edge_rates)
    hc.STORM_FILE = str(STORM_FILE_1967_2024)
    hc.END_YEAR = 2025            # exclusive; RUN_YEARS = 58, states 1967..2025

    hc.check_inputs_exist()
    offsets = hc.load_island_offset_dam()
    elevation_files, dune_files = hc.build_file_lists()
    # The 1971/73 fills are history, applied to every run in the matrix -- the
    # runner requires them, and omitting them would push their signal into the
    # solved edge rate.
    nourish_on, nourish_vol = hc.build_nourishment_arrays_from_manual_inputs()

    # run_one returns a run NAME and writes the matrix to disk; it does not
    # hand back a Cascade. Same extraction the sweep worker uses.
    run_name = hc.run_one("no_groin", offsets, elevation_files, dune_files,
                          nourish_on, nourish_vol)
    matrix = np.load(pathlib.Path(hc.OUTPUT_BASE_DIR) / run_name
                     / f"{run_name}_shoreline_matrix.npy")

    fit_state = FIT_YEAR - hc.START_YEAR
    if not (0 <= fit_state < matrix.shape[0]):
        raise ValueError(
            f"FIT_YEAR {FIT_YEAR} is state {fit_state}, outside this run's "
            f"{matrix.shape[0]} states -- check END_YEAR and the storm file.")

    gis_axis = list(range(hc.FIRST_FILE_NUMBER, hc.LAST_FILE_NUMBER + 1))
    real0 = matrix[0][hc.START_REAL_INDEX:hc.END_REAL_INDEX]
    realN = matrix[fit_state][hc.START_REAL_INDEX:hc.END_REAL_INDEX]
    change = realN - real0
    return {gis: float(change[gis_axis.index(gis)]) for gis in EDGE_GIS}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--max-iter", type=int, default=4)
    parser.add_argument("--tol", type=float, default=1.0,
                        help="convergence tolerance in metres of change")
    args = parser.parse_args()

    target = observed_edge_change()
    print("=" * 70)
    print("EDGE BE SOLVE -- 1967 rig, groin OFF")
    print("=" * 70)
    for gis, value in target.items():
        print(f"  target D{gis}: {value:+.1f} m over {FIT_YEAR - DATUM_YEAR} yr "
              f"({value / (FIT_YEAR - DATUM_YEAR):+.2f} m/yr)")

    # Two brackets to start the secant: the placeholders, and half of them.
    trials = [{2: 5.0, 12: 10.0}, {2: 2.5, 12: 5.0}]
    history = []
    for rates in trials:
        got = modelled_edge_change(rates)
        history.append((rates, got))
        print(f"\n  rates {rates} -> modelled "
              + ", ".join(f"D{g}={v:+.1f} m" for g, v in got.items()))

    for iteration in range(args.max_iter):
        (r0, g0), (r1, g1) = history[-2], history[-1]
        nxt, resid = {}, {}
        for gis in EDGE_GIS:
            e0, e1 = g0[gis] - target[gis], g1[gis] - target[gis]
            resid[gis] = abs(e1)
            if abs(e1 - e0) < 1e-9:
                nxt[gis] = r1[gis]
            else:
                nxt[gis] = r1[gis] - e1 * (r1[gis] - r0[gis]) / (e1 - e0)
            nxt[gis] = float(np.clip(nxt[gis], -50.0, 50.0))
        if max(resid.values()) <= args.tol:
            print(f"\n  converged: max residual {max(resid.values()):.2f} m "
                  f"<= tol {args.tol} m")
            break
        got = modelled_edge_change(nxt)
        history.append((nxt, got))
        print(f"\n  iter {iteration + 1}: rates "
              + ", ".join(f"D{g}={v:+.2f}" for g, v in nxt.items())
              + " -> " + ", ".join(f"D{g}={v:+.1f} m" for g, v in got.items())
              + "  (residual " + ", ".join(f"{abs(got[g]-target[g]):.1f}" for g in EDGE_GIS) + " m)")

    best_rates, best_got = history[-1]
    print("\n" + "=" * 70)
    for gis in EDGE_GIS:
        print(f"  SOLVED D{gis}: {best_rates[gis]:+.2f} m/yr  "
              f"(modelled {best_got[gis]:+.1f} m vs observed {target[gis]:+.1f} m)")
    SOLVED_JSON.write_text(json.dumps(
        {"edge_be_rates_gis": {str(g): best_rates[g] for g in EDGE_GIS},
         "observed_change_m": {str(g): target[g] for g in EDGE_GIS},
         "modelled_change_m": {str(g): best_got[g] for g in EDGE_GIS},
         "fit_year": FIT_YEAR, "groin": "off"}, indent=2), encoding="utf-8")
    print(f"\n  written {SOLVED_JSON}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
