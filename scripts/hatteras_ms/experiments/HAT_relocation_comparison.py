#!/usr/bin/env python3
"""
Does CASCADE relocate NC-12 where and when history did?

    python scripts/hatteras_ms/experiments/HAT_relocation_comparison.py
    python scripts/hatteras_ms/experiments/HAT_relocation_comparison.py --period 1996 --preset edgeBE
    python scripts/hatteras_ms/experiments/HAT_relocation_comparison.py --arm-a DIR --arm-b DIR

Relocations off (the module decides) against relocations on (the measured events):
timing, footprint, near misses and animations, one set per preset. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import argparse
import contextlib
import dataclasses
import datetime
import glob
import io
import json
import os
import sys
from pathlib import Path

import re

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
# Repo root, found by searching upward
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents if (_p / 'pyproject.toml').exists())
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live in scripts/hatteras_ms/.")
for _path in (SCRIPTS_DIR, _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from cascade_pipeline import roadway as roadway_module          # noqa: E402
from site_layer.hatteras_site_config import (                              # noqa: E402
    HATTERAS_DOMAINS,
    HATTERAS_FIRST_ROAD_DOMAIN,
    HATTERAS_LAST_ROAD_DOMAIN,
    HATTERAS_PERIODS,
    HATTERAS_RELOCATION_CHECK_2004,
    HATTERAS_ROAD_EVENTS,
    resolve_be_preset,
)
from cascade_pipeline.roadway import RelocationEvent            # noqa: E402
from cascade_pipeline.run_info import RunInfo                   # noqa: E402
from cascade_pipeline.run_registry import preset_dir_for        # noqa: E402
from cascade_pipeline.plotting.road_relocation_gif import (     # noqa: E402
    make_road_relocation_gif,
    make_topography_gif,
)
from cascade_pipeline.plotting.shoreline_gif import (           # noqa: E402
    DEFAULT_GIF_CONFIG,
)

# --- CONFIG ------------------------------------------------------------------
# The period: module globals, set once by set_period() from --period
DEFAULT_PERIOD = 1984
START_YEAR = DEFAULT_PERIOD
END_YEAR = HATTERAS_PERIODS[DEFAULT_PERIOD]["end_year"]

# The one surveyed road position: the end of 1984-2004, the middle of 1996-2010
CHECK_YEAR = 2004
# -----------------------------------------------------------------------------


# Points the module at one hindcast window
def set_period(start_year):
    global START_YEAR, END_YEAR, OUTPUT_ROOT
    if start_year not in HATTERAS_PERIODS:
        raise SystemExit(f"--period {start_year}: not a hindcast period "
                         f"(have {sorted(HATTERAS_PERIODS)})")
    START_YEAR = int(start_year)
    END_YEAR = int(HATTERAS_PERIODS[START_YEAR]["end_year"])
    OUTPUT_ROOT = (PROJECT_BASE_DIR / "output" / "comparisons"
                   / "relocation" / f"{START_YEAR}_{END_YEAR}")

# The scenario both arms run: full management, the only one where relocation means anything
ARM_SCENARIO_TOKENS = ("road", "bdm", "nogroin")

DEFAULT_PRESET = "zeroBE"


# The two run directory names this comparison reads, for one preset
def arm_names(preset):
    canonical, _ = resolve_be_preset(preset)
    head, tail = ARM_SCENARIO_TOKENS[0], ARM_SCENARIO_TOKENS[1:]
    stem = f"HAT_{START_YEAR}_{END_YEAR}_{canonical}_{head}"
    return (f"{stem}_" + "_".join(tail),
            f"{stem}_reloc_" + "_".join(tail))

# Filename-safe form of a window name
def _slug(text):
    return re.sub(r"[^A-Za-z0-9]+", "", str(text))


# One folder per preset, so one preset's run cannot overwrite another's
OUTPUT_ROOT = PROJECT_BASE_DIR / "output" / "comparisons" / "relocation" / f"{START_YEAR}_{END_YEAR}"

# Two tolerance windows, because the answer is sensitive to it
TOLERANCE_YEARS = (2, 5)

# 1 fps, so each frame can be read: one second per model year
GIF_CONFIG = dataclasses.replace(DEFAULT_GIF_CONFIG, fps=1)

# Windows padded past each event, so the relocating stretch shows against road that is not
GIF_WINDOWS = (
    ("full island", HATTERAS_FIRST_ROAD_DOMAIN, HATTERAS_LAST_ROAD_DOMAIN),
    ("1999 event (GIS 9-14)", 9, 20),
    ("1989 event (GIS 84-87)", 80, 90),
)

# Windows for the topographic raster; one folder per place in a set
PLACE_DIR = {
    "full island": "1-island",
    "1999 event (GIS 9-14)": "2-event-1999_GIS9-14",
    "1989 event (GIS 84-87)": "3-event-1989_GIS84-87",
}
GIF_FILES = {"lines": "dune-and-road.gif", "topography": "topography.gif"}

TOPO_WINDOWS = (
    ("full island", HATTERAS_FIRST_ROAD_DOMAIN, HATTERAS_LAST_ROAD_DOMAIN),
    ("1999 event (GIS 9-14)", 6, 20),
    ("1989 event (GIS 84-87)", 82, 90),
)


# The windows worth drawing for this period
def windows_in_period(windows, targets):
    years = {str(y) for y in targets.values()}
    return tuple(w for w in windows
                 if w[0] == "full island" or w[0].split()[0] in years)

# Written on every topographic frame: the run carries a tenth of the planform
PLANFORM_NOTE = ("planform note: this run's shoreline offset is 1/10 of the "
                 "measured island curvature (shoreline_offset unit mismatch, under review)")


# Loading

# Loads the plan-view shoreline matrix a run saved beside its model
def load_shoreline_matrix(run_dir):
    hits = sorted(glob.glob(os.path.join(run_dir, "*_shoreline_matrix.npy")))
    return np.load(hits[0]) if hits else None


# Back-barrier shoreline per domain per year, metres, raw convention
def back_barrier_matrix(cascade):
    b3d = getattr(cascade, "barrier3d", None)
    if not b3d or not hasattr(b3d[0], "x_b_TS"):
        return None
    n = min(len(bb.x_b_TS) for bb in b3d)
    return np.array([[float(bb.x_b_TS[t]) * 10.0 for bb in b3d]
                     for t in range(n)])


# Loads the pickled Cascade a run wrote with `cascade.save(run_dir)`
def load_cascade(run_dir):
    matches = sorted(glob.glob(os.path.join(run_dir, "*.npz")))
    if not matches:
        raise FileNotFoundError(
            f"no .npz model state in {run_dir}. Re-run with "
            f"HAT_SAVE_MODEL_STATE=1 -- the roadway managers' time series only "
            f"exist inside the saved model.")
    # Cascade.save() changes the working directory, so every path here is absolute
    with np.load(matches[0], allow_pickle=True) as handle:
        return handle["cascade"][0]


# Pulls the per-domain roadway time series out of a finished run
def road_series(cascade, geometry, first_gis, last_gis):
    roadways = getattr(cascade, "roadways", None)
    management = getattr(cascade, "roadway_management_module", None)
    if roadways is None:
        return {}

    out = {}
    for gis in range(first_gis, last_gis + 1):
        pad = geometry.gis_to_pad(gis)
        if not 0 <= pad < len(roadways) or roadways[pad] is None:
            continue
        if management is not None and not management[pad]:
            continue
        manager = roadways[pad]
        out[gis] = {
            "setback": np.asarray(manager._road_setback_TS, dtype=float),
            "relocated": np.asarray(manager._road_relocated_TS, dtype=float),
            "elevation": np.asarray(manager._road_ele_TS, dtype=float),
            "drowned": bool(getattr(manager, "drown_break", 0)),
            "relocation_blocked": bool(getattr(manager, "relocation_break", 0)),
        }
    return out


# Maps each domain a relocation event moves to that event's year
def historical_targets(start_year, end_year):
    return {gis: event.year
            for event in HATTERAS_ROAD_EVENTS
            if isinstance(event, RelocationEvent) and event.enabled
            and start_year <= event.year <= end_year
            for gis in event.displacement_m}


# Scoring

# Calendar year of the first modelled relocation, or None
def first_relocation_year(relocated_ts, start_year):
    fired = np.flatnonzero(relocated_ts > 0)
    return int(start_year + fired[0]) if fired.size else None


# Every calendar year the domain relocated
def relocation_years(relocated_ts, start_year):
    return [int(start_year + i) for i in np.flatnonzero(relocated_ts > 0)]


# How close a road came to firing the relocation trigger, and never did
def relocation_margin(entry, start_year):
    sb = np.asarray(entry["setback"], dtype=float)
    idx = int(np.argmin(sb))
    closest = float(sb[idx])
    return dict(
        min_setback_m=closest,
        closest_year=start_year + idx,
        cells_remaining=int(round(closest / 10.0)),
        migration_needed_m=closest + 10.0,
        years_at_min=int((sb <= closest + 1e-9).sum()),
        end_setback_m=float(sb[-1]),
    )


# Every managed domain that never relocated, ranked by how close it came
def near_miss_table(series, targets, start_year):
    rows = []
    for gis, entry in sorted(series.items()):
        if np.any(np.asarray(entry["relocated"]) > 0):
            continue
        m = relocation_margin(entry, start_year)
        rows.append(dict(
            gis=gis,
            kind="historical" if gis in targets else "control",
            historical_year=targets.get(gis),
            start_setback_m=float(np.asarray(entry["setback"], dtype=float)[0]),
            **m))
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["min_setback_m", "gis"]).reset_index(drop=True)
    return df


# Builds the first-relocation-year table -- the primary result
def score_first_year(series, targets, start_year):
    rows = []
    for gis, event_year in sorted(targets.items()):
        entry = series.get(gis)
        if entry is None:
            rows.append(dict(gis=gis, historical_year=event_year,
                             modelled_first_year=None, error_years=None,
                             n_relocations=0, managed=False,
                             outcome="not managed in this run"))
            continue
        first = first_relocation_year(entry["relocated"], start_year)
        margin = relocation_margin(entry, start_year)
        rows.append(dict(
            gis=gis,
            historical_year=event_year,
            modelled_first_year=first,
            error_years=None if first is None else first - event_year,
            n_relocations=int(np.sum(entry["relocated"] > 0)),
            managed=True,
            outcome=("never relocated" if first is None else
                     "drowned" if entry["drowned"] else
                     "relocation blocked" if entry["relocation_blocked"] else
                     "relocated"),
            # Near misses only where the trigger never fired: a reset setback is not a near miss
            **({k: None for k in ("min_setback_m", "closest_year",
                                  "cells_remaining", "migration_needed_m",
                                  "years_at_min", "end_setback_m")}
               if first is not None else margin),
        ))
    return pd.DataFrame(rows)


# Hit/miss matrix at one tolerance window
def score_confusion(series, targets, start_year, tolerance):
    hits, misses = [], []
    for gis, event_year in targets.items():
        entry = series.get(gis)
        if entry is None:
            misses.append(gis)
            continue
        years = relocation_years(entry["relocated"], start_year)
        if any(abs(y - event_year) <= tolerance for y in years):
            hits.append(gis)
        else:
            misses.append(gis)

    control = sorted(set(series) - set(targets))
    false_pos = [gis for gis in control
                 if np.any(series[gis]["relocated"] > 0)]
    return {
        "tolerance_years": tolerance,
        "historical_domains": len(targets),
        "hits": len(hits),
        "misses": len(misses),
        "hit_domains": hits,
        "miss_domains": misses,
        "control_domains": len(control),
        "false_positives": len(false_pos),
        "false_positive_domains": false_pos,
        "recall": len(hits) / len(targets) if targets else float("nan"),
        "false_positive_rate": (len(false_pos) / len(control)
                                if control else float("nan")),
    }


# Index of the last year the manager actually ran for this domain
def last_managed_index(entry):
    written = np.flatnonzero(np.asarray(entry["elevation"], dtype=float) != 0)
    return int(written[-1]) if written.size else -1


# Per-domain setback trajectory comparison between the two arms
def score_trajectories(series_a, series_b, targets, start_year, check_2004):
    summary, long = [], []
    for gis in sorted(set(series_a) | set(series_b)):
        a = series_a.get(gis, {}).get("setback")
        b = series_b.get(gis, {}).get("setback")
        if a is None or b is None:
            continue
        n = min(len(a), len(b))
        a, b = a[:n], b[:n]
        last_a = last_managed_index(series_a[gis])
        last_b = last_managed_index(series_b[gis])
        common = min(last_a, last_b)
        if common < 0:
            continue

        years = start_year + np.arange(n)
        for i in range(n):
            long.append(dict(gis=gis, year=int(years[i]),
                             setback_free_m=float(a[i]),
                             setback_prescribed_m=float(b[i]),
                             managed_free=bool(i <= last_a),
                             managed_prescribed=bool(i <= last_b)))

        wa, wb = a[:common + 1], b[:common + 1]
        # The modelled position in the check year per arm, or None if management had stopped
        k = CHECK_YEAR - start_year
        at_a = float(a[k]) if 0 <= k <= last_a else None
        at_b = float(b[k]) if 0 <= k <= last_b else None
        summary.append(dict(
            gis=gis,
            historical=gis in targets,
            event_year=targets.get(gis),
            free_start_m=float(a[0]),
            free_last_m=float(a[last_a]),
            free_last_year=int(start_year + last_a),
            prescribed_start_m=float(b[0]),
            prescribed_last_m=float(b[last_b]),
            prescribed_last_year=int(start_year + last_b),
            compared_through=int(start_year + common),
            difference_at_common_m=float(wb[-1] - wa[-1]),
            rmse_m=float(np.sqrt(np.mean((wa - wb) ** 2))),
            check_year=CHECK_YEAR,
            free_at_check_m=at_a,
            prescribed_at_check_m=at_b,
            measured_2004_m=check_2004.get(gis),
        ))
    return pd.DataFrame(summary), pd.DataFrame(long)


# Confirms the two arms are identical before the first prescribed event
def check_determinism(series_a, series_b, first_event_year, start_year):
    cutoff = first_event_year - start_year        # exclusive
    offenders = []
    for gis in sorted(set(series_a) & set(series_b)):
        for key in ("setback", "relocated", "elevation"):
            a = series_a[gis][key][:cutoff]
            b = series_b[gis][key][:cutoff]
            if not np.array_equal(a, b):
                offenders.append(gis)
                break
    return not offenders, offenders


# Confirms index i really is calendar year start_year + i
def check_event_indexing(series_b, targets, start_year):
    rows = []
    for gis, event_year in sorted(targets.items()):
        entry = series_b.get(gis)
        if entry is None:
            continue
        setback = entry["setback"]
        # Search only up to the last managed year, or unwritten years add a false jump
        stop = last_managed_index(entry) + 1 or len(setback)
        jumps = np.diff(setback[:stop])
        if not jumps.size:
            continue
        idx = int(np.argmax(jumps))
        rows.append(dict(
            gis=gis, expected_year=event_year,
            largest_jump_year=int(start_year + idx + 1),
            largest_jump_m=float(jumps[idx]),
        ))
    df = pd.DataFrame(rows)
    if not df.empty:
        df["matches"] = df["largest_jump_year"] == df["expected_year"]
    return df


# Report

# Mirrors everything printed to the console into a buffer
class _Tee:

    def __init__(self):
        self._real = sys.stdout
        self._buf = io.StringIO()

    def write(self, text):
        self._real.write(text)
        self._buf.write(text)
        return len(text)

    def flush(self):
        self._real.flush()

    def getvalue(self):
        return self._buf.getvalue()


# One line describing which run a comparison arm actually read
def _arm_provenance(run_dir):
    hits = sorted(glob.glob(os.path.join(str(run_dir), "*_run_metadata.json")))
    if not hits:
        return f"{os.path.basename(str(run_dir))}   (no run metadata found)"
    with open(hits[0], "r", encoding="utf-8") as fh:
        ident = json.load(fh).get("identity", {})
    dirty = "  DIRTY TREE" if ident.get("git_dirty") else ""
    return (f"{ident.get('run_name', '?')}\n"
            f"      run {ident.get('timestamp', '?')}   "
            f"topo {ident.get('topo_product', '?')}/"
            f"{ident.get('topo_dune_version', '?')}   "
            f"commit {str(ident.get('git_commit', '?'))[:12]}{dirty}")


# The dune-topo version a run was made on, from its metadata ('v?' if none)
def _topo_version(run_dir):
    hits = sorted(glob.glob(os.path.join(str(run_dir), "*_run_metadata.json")))
    if not hits:
        return "v?"
    with open(hits[0], "r", encoding="utf-8") as fh:
        return str(json.load(fh).get("identity", {}).get("topo_dune_version", "v?"))


# The output folder, named for the dune-topo version and preset
def default_out_dir(arm_a, arm_b, preset):
    va, vb = _topo_version(arm_a), _topo_version(arm_b)
    if va != vb:
        raise SystemExit(f"the two arms are on different dune-topo versions ({va} vs {vb}); "
                         f"a comparison across versions is not what this script measures")
    groin = "_groin" if ("_groin" in os.path.basename(str(arm_a)) and "nogroin" not in os.path.basename(str(arm_a))) else ""
    return OUTPUT_ROOT / va / f"{preset}{groin}"


# Provenance block written above the captured output
def _report_header(arm_a, arm_b, preset):
    return (
        "=" * 74 + "\n"
        f"generated   {datetime.datetime.now():%Y-%m-%d %H:%M:%S} by "
        f"{os.path.basename(__file__)}\n"
        f"period      {START_YEAR}-{END_YEAR}\n"
        f"preset      {preset}\n"
        f"arm A       {_arm_provenance(arm_a)}\n"
        f"arm B       {_arm_provenance(arm_b)}\n"
        "\n"
        "This file is the console output of the run that wrote the CSVs and\n"
        "GIFs beside it. If the arm identities above do not match the runs on\n"
        "disk, this report is stale -- re-run the comparison.\n"
        + "=" * 74 + "\n\n")


# Run: resolve the arms, compare them, write the report
def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    raw_runs = PROJECT_BASE_DIR / "output" / "raw_runs"
    parser.add_argument("--period", type=int, default=DEFAULT_PERIOD,
                        help="hindcast start year, a HATTERAS_PERIODS key "
                             f"(default {DEFAULT_PERIOD}); sets the window, "
                             "the arm names and the output root")
    parser.add_argument("--preset", default=DEFAULT_PRESET,
                        help="source/sink preset; names both arms and the "
                             f"output subdirectory (default {DEFAULT_PRESET})")
    parser.add_argument("--arm-a", default=None,
                        help="run directory with relocations OFF "
                             "(overrides --preset)")
    parser.add_argument("--arm-b", default=None,
                        help="run directory with relocations ON "
                             "(overrides --preset)")
    parser.add_argument("--out", default=None,
                        help="output directory "
                             "(default OUTPUT_ROOT/<topo version of the arms>/<preset>[_groin])")
    args = parser.parse_args()
    set_period(args.period)

    preset, _ = resolve_be_preset(args.preset)
    name_a, name_b = arm_names(preset)
    # Both arms share a preset, so their folder is resolved once
    period_dir = preset_dir_for(raw_runs, (START_YEAR, END_YEAR), preset)
    args.arm_a = args.arm_a or str(period_dir / name_a)
    args.arm_b = args.arm_b or str(period_dir / name_b)

    out_dir = Path(args.out).resolve() if args.out else default_out_dir(args.arm_a, args.arm_b, preset)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Delete the old report first, so a failed run leaves none rather than a stale one
    report_path = out_dir / "report.txt"
    if report_path.exists():
        report_path.unlink()

    tee = _Tee()
    failure = None
    try:
        with contextlib.redirect_stdout(tee):
            _compare(args, preset, out_dir)
    except BaseException as exc:                      # noqa: BLE001
        failure = exc
        raise
    finally:
        body = tee.getvalue()
        if failure is not None:
            body += ("\n" + "!" * 74 + "\n"
                     f"RUN FAILED before completing: "
                     f"{type(failure).__name__}: {failure}\n"
                     "This report is PARTIAL. The artifacts in this folder are "
                     "incomplete or absent.\n" + "!" * 74 + "\n")
        report_path.write_text(
            _report_header(args.arm_a, args.arm_b, preset) + body,
            encoding="utf-8")
        print(f"report -> {report_path}")


# The comparison itself; everything it prints becomes report.txt
def _compare(args, preset, out_dir):
    tables_dir = out_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    print("=" * 74)
    print("NC-12 RELOCATION: emergent vs prescribed")
    print("=" * 74)
    print(f"  preset              {preset}")
    print(f"  arm A (free)        {args.arm_a}")
    print(f"  arm B (prescribed)  {args.arm_b}")

    span = (HATTERAS_FIRST_ROAD_DOMAIN, HATTERAS_LAST_ROAD_DOMAIN)
    cascade_a = load_cascade(args.arm_a)
    cascade_b = load_cascade(args.arm_b)
    series_a = road_series(cascade_a, HATTERAS_DOMAINS, *span)
    series_b = road_series(cascade_b, HATTERAS_DOMAINS, *span)
    targets = historical_targets(START_YEAR, END_YEAR)
    if not targets:
        raise SystemExit(f"no relocation event falls in {START_YEAR}-{END_YEAR}; "
                         "there is nothing for this comparison to score")

    print(f"  period              {START_YEAR}-{END_YEAR}")
    print(f"\n  managed road domains  A: {len(series_a)}   B: {len(series_b)}")
    print(f"  historical domains    {sorted(targets)}")
    for year in sorted(set(targets.values())):
        moved = sorted(g for g, y in targets.items() if y == year)
        print(f"    {year}  GIS {moved}")

    # 0. the checks that decide whether the rest means anything
    print("\n" + "-" * 74)
    print("0. VALIDITY CHECKS")
    print("-" * 74)

    first_event = min(targets.values())
    ok, offenders = check_determinism(series_a, series_b, first_event,
                                      START_YEAR)
    print(f"  arms identical {START_YEAR}-{first_event - 1}: "
          f"{'YES' if ok else 'NO'}")
    if not ok:
        print(f"    !! diverged before the first event in GIS {offenders}")
        print("    !! that is a bug, not a result -- everything below is "
              "uninterpretable")

    idx_df = check_event_indexing(series_b, targets, START_YEAR)
    if not idx_df.empty:
        agree = int(idx_df["matches"].sum())
        print(f"  prescribed jumps land on the event year: "
              f"{agree}/{len(idx_df)} domains")
        if agree != len(idx_df):
            print(idx_df[~idx_df["matches"]].to_string(index=False))
        idx_df.to_csv(tables_dir / "indexing_check.csv", index=False)

    # 1. first relocation year
    print("\n" + "-" * 74)
    print("1. FIRST MODELLED RELOCATION, free-running arm")
    print("-" * 74)
    first_df = score_first_year(series_a, targets, START_YEAR)
    print(first_df.to_string(index=False))
    first_df.to_csv(tables_dir / "first_relocation_year.csv", index=False)

    dated = first_df.dropna(subset=["error_years"])
    if not dated.empty:
        err = dated["error_years"].astype(float)
        print(f"\n  domains that relocated at all   {len(dated)}/{len(first_df)}")
        print(f"  median signed error             {err.median():+.1f} yr")
        print(f"  mean absolute error             {err.abs().mean():.1f} yr")
    else:
        print("\n  no historical domain relocated in the free-running arm")

    # 1b. how close did the misses come?
    print("\n" + "-" * 74)
    print("1b. NEAR MISSES: how much more dune migration would have fired it")
    print("-" * 74)
    near = near_miss_table(series_a, targets, START_YEAR)
    near.to_csv(tables_dir / "near_miss_margin.csv", index=False)
    if near.empty:
        print("  every managed domain relocated at least once")
    else:
        print("  The trigger is `setback < 0`, strict, and the setback moves in")
        print("  whole 10 m cells (ShorelineChangeTS counts cells). So a road at")
        print("  setback 0 has NOT fired -- it needs one more full cell. "
              "'needed'")
        print("  below is that extra landward migration, = closest + 10 m.\n")
        hist = near[near["kind"] == "historical"]
        if not hist.empty:
            print("  HISTORICAL domains that never relocated:")
            print("   gis  hist  start_m  closest_m  at_year  cells_left  "
                  "needed_m")
            for _, r in hist.iterrows():
                print(f"  {r['gis']:>4} {int(r['historical_year']):>5} "
                      f"{r['start_setback_m']:>8.0f} {r['min_setback_m']:>10.0f} "
                      f"{int(r['closest_year']):>8} {r['cells_remaining']:>11} "
                      f"{r['migration_needed_m']:>9.0f}")
        ctrl = near[near["kind"] == "control"]
        if not ctrl.empty:
            print(f"\n  CONTROL domains (history never relocated these), "
                  f"closest 10 of {len(ctrl)}:")
            print("   gis  start_m  closest_m  at_year  cells_left  needed_m")
            for _, r in ctrl.head(10).iterrows():
                print(f"  {r['gis']:>4} {r['start_setback_m']:>8.0f} "
                      f"{r['min_setback_m']:>10.0f} {int(r['closest_year']):>8} "
                      f"{r['cells_remaining']:>11} "
                      f"{r['migration_needed_m']:>9.0f}")
            print(f"\n  control margin: median {ctrl['min_setback_m'].median():.0f} m, "
                  f"min {ctrl['min_setback_m'].min():.0f} m")
            tight = ctrl[ctrl["cells_remaining"] <= 1]
            print(f"  control domains within ONE cell of firing: {len(tight)}"
                  + (f"  GIS {tight['gis'].tolist()}" if len(tight) else ""))
            print("  -- a 0/45 false-positive count is only reassuring if these")
            print("     margins are wide; this is where that gets checked.")
        print(f"\n  saved -> near_miss_margin.csv")

    # 2. hit / miss
    print("\n" + "-" * 74)
    print("2. HIT / MISS, with false positives")
    print("-" * 74)
    conf_rows = []
    for tol in TOLERANCE_YEARS:
        conf = score_confusion(series_a, targets, START_YEAR, tol)
        conf_rows.append({k: v for k, v in conf.items()
                          if not k.endswith("_domains") or k.endswith("s")})
        print(f"\n  +/-{tol} yr   hits {conf['hits']}/"
              f"{conf['historical_domains']}  "
              f"(recall {conf['recall']:.2f})")
        print(f"           misses      GIS {conf['miss_domains']}")
        print(f"           false pos   {conf['false_positives']}/"
              f"{conf['control_domains']} control domains "
              f"(rate {conf['false_positive_rate']:.2f})")
        if conf["false_positives"]:
            shown = conf["false_positive_domains"][:20]
            more = "" if len(shown) == conf["false_positives"] else " ..."
            print(f"                       GIS {shown}{more}")
    pd.DataFrame([{k: (v if not isinstance(v, list) else " ".join(map(str, v)))
                   for k, v in r.items()} for r in conf_rows]).to_csv(
        tables_dir / "confusion.csv", index=False)

    # 3. setback trajectories
    print("\n" + "-" * 74)
    print("3. SETBACK TRAJECTORIES")
    print("-" * 74)
    traj_df, long_df = score_trajectories(
        series_a, series_b, targets, START_YEAR, HATTERAS_RELOCATION_CHECK_2004)
    hist = traj_df[traj_df["historical"]]
    where = ("the last year of the window" if CHECK_YEAR == END_YEAR else
             f"year {CHECK_YEAR - START_YEAR} of {END_YEAR - START_YEAR}, mid-window")
    print(f"  position cross-check at {CHECK_YEAR}, {where}: measured_2004_m is")
    print("  RoadOffset_2004 (2008 line vs 2004-start row 0); free/prescribed_at_check_m")
    print("  are each arm's modelled setback that year.\n")
    print(hist.to_string(index=False))
    traj_df.to_csv(tables_dir / "setback_summary.csv", index=False)
    long_df.to_csv(tables_dir / "setback_by_year.csv", index=False)
    print(f"\n  saved per-year setbacks for {long_df['gis'].nunique()} domains "
          f"-> setback_by_year.csv")

    # 4. outcomes
    print("\n" + "-" * 74)
    print("4. ROAD OUTCOMES: did prescribing the relocations change the fate?")
    print("-" * 74)
    out_rows = []
    for label, cas in (("free", cascade_a), ("prescribed", cascade_b)):
        for row in roadway_module.summarise_road_management(
                cas, HATTERAS_DOMAINS, *span):
            out_rows.append(dict(arm=label, **row))
    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(tables_dir / "road_outcomes.csv", index=False)

    for label in ("free", "prescribed"):
        arm = out_df[out_df["arm"] == label]
        print(f"  {label:<11} drowned {int(arm['drowned'].sum()):>3}   "
              f"blocked {int(arm['relocation_blocked'].sum()):>3}   "
              f"relocations {int(arm['relocations'].sum()):>5}   "
              f"of {len(arm)} managed domains")

    wide = out_df.pivot(index="gis", columns="arm",
                        values=["drowned", "relocation_blocked", "reason"])
    changed = wide[wide[("reason", "free")] != wide[("reason", "prescribed")]]
    if changed.empty:
        print("\n  no domain changed outcome between the two arms")
    else:
        print(f"\n  {len(changed)} domain(s) changed outcome:")
        for gis in changed.index:
            print(f"    GIS {gis:>3}  free: {wide.loc[gis, ('reason', 'free')]:<20}"
                  f"  prescribed: {wide.loc[gis, ('reason', 'prescribed')]}")

    # 5. animation
    print()
    print("-" * 74)
    print("5. ANIMATION")
    print("-" * 74)
    shore_a = load_shoreline_matrix(args.arm_a)
    shore_b = load_shoreline_matrix(args.arm_b)
    if shore_a is None or shore_b is None:
        print("  skipped: a run has no *_shoreline_matrix.npy to animate")
    else:
        info_a = RunInfo(run_name=os.path.basename(args.arm_a),
                         run_dir=args.arm_a, start_year=START_YEAR,
                         end_year=END_YEAR)
        info_b = RunInfo(run_name=os.path.basename(args.arm_b),
                         run_dir=args.arm_b, start_year=START_YEAR,
                         end_year=END_YEAR)
        back_a = back_barrier_matrix(cascade_a)
        back_b = back_barrier_matrix(cascade_b)
        for name, lo, hi in windows_in_period(GIF_WINDOWS, targets):
            place = out_dir / PLACE_DIR[name]
            place.mkdir(parents=True, exist_ok=True)
            make_road_relocation_gif(
                (shore_a, info_a), (shore_b, info_b), series_a, series_b, lo, hi,
                str(place / GIF_FILES["lines"]), back_a=back_a, back_b=back_b,
                event_years=targets, gif_config=GIF_CONFIG,
                title=f"NC-12 and the dune line \u2014 {name}")
        for name, lo, hi in windows_in_period(TOPO_WINDOWS, targets):
            place = out_dir / PLACE_DIR[name]
            place.mkdir(parents=True, exist_ok=True)
            make_topography_gif(
                cascade_a, cascade_b, series_a, series_b, lo, hi,
                str(place / GIF_FILES["topography"]),
                START_YEAR, event_years=targets, gif_config=GIF_CONFIG,
                title=f"Hatteras topography and NC-12 \u2014 {name}",
                planform_note=PLANFORM_NOTE)

    print("\n" + "=" * 74)
    print(f"artifacts -> {out_dir}")
    print("=" * 74)


if __name__ == "__main__":
    main()
