#!/usr/bin/env python3
"""Does CASCADE relocate NC-12 where and when history did?

THE QUESTION
    `roadway_manager` relocates the road on its own when the dune line
    overruns it -- `road_relocation_checks` fires the moment the setback goes
    negative. Separately, the pipeline can PRESCRIBE the two historical
    relocations (1989, GIS 84-87; 1999, GIS 9-14) as measured displacements.
    This script runs the two against each other:

        arm A  relocations OFF -- the module decides on its own
        arm B  relocations ON  -- the measured displacements are applied

    and asks whether A reproduces B's timing and footprint unaided.

WHAT IS DELIBERATELY NOT DONE
    Arm A runs the module EXACTLY as built. `_road_relocation_setback` is left
    at its initialised value -- the run's starting setback -- rather than being
    set to a design distance or to the measured post-relocation alignment.
    That matters because six of the ten historical domains (10-13, 85-87)
    start at a setback of 0.0 m, so in those domains a relocation puts the road
    back where it already was and the trigger re-fires the next year the dune
    line moves landward. That ratcheting is a property of the module, and
    reporting it is the point; designing around it would hide it.

WHAT THE COMPARISON CAN AND CANNOT SAY
    CASCADE's trigger is purely geometric: dune line overruns road. There is no
    storm damage, no cost, and no maintenance decision in it. So a match here
    means "the modelled physics would have overrun NC-12 near that year", NOT
    "NCDOT would have moved the road then". The second question is outside what
    this module represents, and no configuration of it gets there.

BACKGROUND EROSION
    Both arms run under whichever source/sink preset the runs were driven with.
    Note that `edgeBE` carries rates on GIS 1 and 90 ONLY, so at every domain
    under test here edgeBE and zeroBE are the same forcing: the dune-line
    retreat that fires the trigger comes entirely from Barrier3D/BRIE dynamics.
    Only `calibBE` puts a background-erosion term on the relocation domains,
    which makes a calibBE re-run the natural sensitivity test once those
    source/sink terms are updated.

TIME INDEXING
    `RoadwayManager` writes every time series at `time_index - 1`, and the
    loop applies a historical event before the update for `start_year +
    time_step`. So index i in `_road_setback_TS` is calendar year
    START_YEAR + i. The prescribed arm is used to CHECK that rather than
    assume it: arm B's setbacks must jump at exactly 1989 and 1999.

USAGE
    python scripts/hatteras_ms/HAT_relocation_comparison.py
    python scripts/hatteras_ms/HAT_relocation_comparison.py --arm-a DIR --arm-b DIR

Author: Hannah A. Henry, UNC CECL
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
PROJECT_BASE_DIR = _HERE.parents[2]
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live in scripts/hatteras_ms/.")
for _path in (SCRIPTS_DIR, _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from cascade_pipeline import roadway as roadway_module          # noqa: E402
from hatteras_site_config import (                              # noqa: E402
    HATTERAS_DOMAINS,
    HATTERAS_FIRST_ROAD_DOMAIN,
    HATTERAS_LAST_ROAD_DOMAIN,
    HATTERAS_RELOCATION_CHECK_2004,
    HATTERAS_ROAD_EVENTS,
    resolve_be_preset,
)
from cascade_pipeline.roadway import RelocationEvent            # noqa: E402
from cascade_pipeline.run_info import RunInfo                   # noqa: E402
from cascade_pipeline.run_registry import preset_dir_for        # noqa: E402
from cascade_pipeline.plotting.road_relocation_gif import (     # noqa: E402
    make_all_road_gifs,
    make_topography_gif,
)
from cascade_pipeline.plotting.shoreline_gif import (           # noqa: E402
    DEFAULT_GIF_CONFIG,
)

START_YEAR = 1984
END_YEAR = 2004

# The scenario both arms run. `full_management` is the status-quo hindcast,
# and it is the only scenario where a relocation arm is meaningful and the
# village management is also present.
ARM_SCENARIO_TOKENS = ("road", "bdm", "nogroin")

DEFAULT_PRESET = "zeroBE"


def arm_names(preset):
    """The two run directory names this comparison reads, for one preset.

    Built with the same token rule the hindcast derives RUN_NAME with in its
    section 7.5: the arms are identical in every token except `reloc`. Both
    are derived from ONE preset argument on purpose -- an arm A read against
    an arm B from a different preset would produce a clean-looking comparison
    of two different forcings, and nothing downstream would notice, because
    every check in this file is a difference between the arms.

    Args:
        preset: A source/sink preset key or deprecated alias.

    Returns:
        (arm_a_name, arm_b_name), relocations off and on.
    """
    canonical, _ = resolve_be_preset(preset)
    head, tail = ARM_SCENARIO_TOKENS[0], ARM_SCENARIO_TOKENS[1:]
    stem = f"HAT_{START_YEAR}_{END_YEAR}_{canonical}_{head}"
    return (f"{stem}_" + "_".join(tail),
            f"{stem}_reloc_" + "_".join(tail))

def _slug(text):
    """Filename-safe form of a window name."""
    return re.sub(r"[^A-Za-z0-9]+", "", str(text))


# One directory per preset underneath it. The background-erosion preset is
# the axis this comparison is repeated over -- edgeBE and zeroBE differ only
# at GIS 1 and 90, while calibBE is the only one carrying a source/sink term
# on the relocation domains themselves -- so the artifacts have to be kept
# apart or the second run silently overwrites the first.
OUTPUT_ROOT = PROJECT_BASE_DIR / "output" / "comparisons" / "relocation_1984_2004"

# Tolerance windows for the hit/miss matrix. Two are reported rather than one
# because the answer is sensitive to it and a single number would hide that.
TOLERANCE_YEARS = (2, 5)

# These animations are READ, not watched. The question they answer -- in which
# year does this domain's road step landward, and does the model do it when
# history did -- needs the viewer to hold one frame long enough to find the
# domain, read the year clock, and compare the two panels. At the shared
# default of 3 fps that is 333 ms per frame and the eye cannot do it; a 20-year
# run is over in seven seconds.
#
# 1 fps gives a second per model year and a ~21 s loop. Slower than a general
# shoreline animation wants, which is why this overrides rather than changing
# GifConfig's default: the per-run shoreline GIFs show a smooth trend where 3
# fps reads fine, and it is only the discrete, dated relocation events that
# need dwell time.
GIF_CONFIG = dataclasses.replace(DEFAULT_GIF_CONFIG, fps=1)

# Alongshore windows for the animation. The two event windows are padded a
# few domains beyond the event footprint so the relocating stretch is seen
# against road that is NOT relocating -- an unpadded window shows every
# domain stepping at once and reads as a global effect.
GIF_WINDOWS = (
    ("full island", HATTERAS_FIRST_ROAD_DOMAIN, HATTERAS_LAST_ROAD_DOMAIN),
    ("1999 event (GIS 9-14)", 9, 20),
    ("1989 event (GIS 84-87)", 80, 90),
)

# Windows for the topographic raster. The full island is included despite
# being 4100 cells wide and a few hundred tall -- rendered wide and short,
# that IS the shape of Hatteras, and it is the view that reads as the island
# rather than as a chart. The event windows carry the cross-shore detail.
TOPO_WINDOWS = (
    ("full island", HATTERAS_FIRST_ROAD_DOMAIN, HATTERAS_LAST_ROAD_DOMAIN),
    ("1999 event (GIS 9-14)", 6, 20),
    ("1989 event (GIS 84-87)", 82, 90),
)

# Recorded on every topographic frame. The measured island planform spans
# 6.3 km of cross-shore offset across the real domains, but Cascade's
# shoreline_offset reaches BRIE in decameters where BRIE reads metres, so
# the run carries a tenth of it. Stated on the figure rather than silently
# corrected: the animation shows the island the model actually ran.
PLANFORM_NOTE = ("planform note: this run's shoreline offset is 1/10 of the "
                 "measured island curvature (shoreline_offset unit mismatch, under review)")


# =============================================================================
# Loading
# =============================================================================

def load_shoreline_matrix(run_dir):
    """Loads the plan-view shoreline matrix a run saved beside its model.

    Args:
        run_dir: The run directory.

    Returns:
        A 2-D [n_years, total_domains] array in metres, raw x_s_TS
        convention, or None if the run did not save one.
    """
    hits = sorted(glob.glob(os.path.join(run_dir, "*_shoreline_matrix.npy")))
    return np.load(hits[0]) if hits else None


def back_barrier_matrix(cascade):
    """Back-barrier shoreline per domain per year, metres, raw convention.

    The shoreline matrix the run saves carries only x_s. Drawing the island
    with any width needs x_b as well, on the same sign convention so the two
    can share an axis.

    Args:
        cascade: A finished Cascade instance.

    Returns:
        A [n_years, n_domains] array in metres, or None if x_b_TS is absent.
    """
    b3d = getattr(cascade, "barrier3d", None)
    if not b3d or not hasattr(b3d[0], "x_b_TS"):
        return None
    n = min(len(bb.x_b_TS) for bb in b3d)
    return np.array([[float(bb.x_b_TS[t]) * 10.0 for bb in b3d]
                     for t in range(n)])


def load_cascade(run_dir):
    """Loads the pickled Cascade a run wrote with `cascade.save(run_dir)`.

    Args:
        run_dir: Directory holding exactly one .npz model state.

    Returns:
        The Cascade instance.

    Raises:
        FileNotFoundError: If the directory holds no .npz, which means the run
            was driven with HAT_SAVE_MODEL_STATE=false and cannot be compared.
    """
    matches = sorted(glob.glob(os.path.join(run_dir, "*.npz")))
    if not matches:
        raise FileNotFoundError(
            f"no .npz model state in {run_dir}. Re-run with "
            f"HAT_SAVE_MODEL_STATE=1 -- the roadway managers' time series only "
            f"exist inside the saved model.")
    # np.load changes nothing on disk, but Cascade.save() os.chdir()s into the
    # run directory, so anything downstream that assumes cwd is the repo root
    # must not rely on it. Absolute paths are used throughout this file.
    with np.load(matches[0], allow_pickle=True) as handle:
        return handle["cascade"][0]


def road_series(cascade, geometry, first_gis, last_gis):
    """Pulls the per-domain roadway time series out of a finished run.

    Reads the managers CASCADE already holds rather than re-deriving anything.
    Only domains the run actually managed are returned: a domain outside
    `roadway_management_module` has a RoadwayManager object that was never
    called, and its all-zero series would read as "never relocated" rather
    than "never asked".

    Args:
        cascade: A Cascade instance after its run.
        geometry: DomainGeometry describing the padded array.
        first_gis: First GIS domain carrying road.
        last_gis: Last GIS domain carrying road.

    Returns:
        A {gis: dict} mapping, each dict holding setback, relocated, elevation
        (numpy arrays indexed by years-since-START_YEAR) plus the drowned and
        relocation_blocked flags.
    """
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


def historical_targets(start_year, end_year):
    """Maps each domain a relocation event moves to that event's year.

    Read from HATTERAS_ROAD_EVENTS rather than retyped, so dropping a domain
    from an event (as GIS 15 was) drops it from the scoring too.

    Args:
        start_year: First calendar year of the period.
        end_year: Last calendar year of the period.

    Returns:
        A {gis: event_year} dict for relocation events inside the period.
    """
    return {gis: event.year
            for event in HATTERAS_ROAD_EVENTS
            if isinstance(event, RelocationEvent) and event.enabled
            and start_year <= event.year <= end_year
            for gis in event.displacement_m}


# =============================================================================
# Scoring
# =============================================================================

def first_relocation_year(relocated_ts, start_year):
    """Calendar year of the first modelled relocation, or None.

    Args:
        relocated_ts: `_road_relocated_TS`, indexed by years since start_year.
        start_year: First calendar year of the run.

    Returns:
        The calendar year, or None if the domain never relocated.
    """
    fired = np.flatnonzero(relocated_ts > 0)
    return int(start_year + fired[0]) if fired.size else None


def relocation_years(relocated_ts, start_year):
    """Every calendar year the domain relocated."""
    return [int(start_year + i) for i in np.flatnonzero(relocated_ts > 0)]


def relocation_margin(entry, start_year):
    """How close a road came to firing the relocation trigger, and never did.

    THE TRIGGER, EXACTLY. `roadway_manager.road_relocation_checks` does

        road_setback = road_setback + dune_migrated      # dune_migrated < 0 landward
        if road_setback < 0:  relocate

    and the caller supplies

        dune_migration = barrier3d.ShorelineChangeTS[t-1] * 10       # m

    `ShorelineChangeTS` counts WHOLE dam cells, so the setback only ever moves
    in 10 m steps and the test is STRICT. Both facts matter for a margin:

      * a setback sitting at exactly 0.0 has NOT fired. The dune line has
        reached the road and stopped there. It needs one more full cell.
      * so the extra landward migration needed is `min_setback + 10 m`, not
        `min_setback`. Reporting the setback alone would say GIS 84 needed 0 m
        more, which is wrong -- it needed one more cell.

    Args:
        entry: One road_series() value, with its "setback" array.
        start_year: First calendar year of the run.

    Returns:
        dict with the closest approach, the year it happened, the cells still
        between dune and road there, and the extra migration that would fire.
    """
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


def near_miss_table(series, targets, start_year):
    """Every managed domain that never relocated, ranked by how close it came.

    Covers the CONTROL domains as well as the historical ones on purpose. The
    comparison's headline false-positive count is 0/45, which reads as "the
    module is appropriately conservative" -- but a control domain sitting one
    cell from firing is a different statement about robustness than one sitting
    thirty cells away, and only this table separates them.

    Args:
        series: road_series() output for the free-running arm.
        targets: {gis: historical_event_year}.
        start_year: First calendar year of the run.

    Returns:
        DataFrame of non-relocating domains, closest first.
    """
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


def score_first_year(series, targets, start_year):
    """Builds the first-relocation-year table -- the primary result.

    Robust to the ratcheting described in the module docstring: however many
    times a domain relocates afterwards, the FIRST firing is the model's
    answer to "when did the dune line reach the road".

    Args:
        series: road_series() output for the free-running arm.
        targets: {gis: historical_event_year}.
        start_year: First calendar year of the run.

    Returns:
        A DataFrame, one row per historical domain.
    """
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
            # Only meaningful where the trigger never fired. On a domain that
            # DID relocate the setback is reset by _apply_relocation, so its
            # minimum is an artifact of the reset rather than a near miss.
            **({k: None for k in ("min_setback_m", "closest_year",
                                  "cells_remaining", "migration_needed_m",
                                  "years_at_min", "end_setback_m")}
               if first is not None else margin),
        ))
    return pd.DataFrame(rows)


def score_confusion(series, targets, start_year, tolerance):
    """Hit/miss matrix at one tolerance window.

    A HIT is a historical domain that relocated within +/-tolerance years of
    its event. A FALSE POSITIVE is a managed road domain history never
    relocated that relocated anyway, at any time. The false-positive count is
    the half of this that a per-domain error table cannot show: a module that
    relocates everywhere scores perfectly on the historical domains while
    saying nothing.

    Args:
        series: road_series() output for the free-running arm.
        targets: {gis: historical_event_year}.
        start_year: First calendar year of the run.
        tolerance: Half-width of the match window, in years.

    Returns:
        A dict of counts and rates.
    """
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


def last_managed_index(entry):
    """Index of the last year the manager actually ran for this domain.

    Dated from `_road_ele_TS`, NOT from the setback series. A setback of
    exactly 0.0 m is legitimate here -- it means the road sits on the dune
    line, which is where six of the ten historical domains start -- so a
    nonzero test on the setback would report those domains as never managed.
    Road ELEVATION has no such ambiguity: the module stops managing the moment
    it drops below 0 m MHW, so its last non-zero entry dates the last managed
    year. This is the same signal `summarise_road_management` dates from.

    Args:
        entry: One road_series() value, carrying "elevation".

    Returns:
        The index, or -1 if the manager never ran.
    """
    written = np.flatnonzero(np.asarray(entry["elevation"], dtype=float) != 0)
    return int(written[-1]) if written.size else -1


def score_trajectories(series_a, series_b, targets, start_year, check_2004):
    """Per-domain setback trajectory comparison between the two arms.

    Compared over the window BOTH arms were still managing the domain. If one
    arm's road drowns, its series stops being written; extending the
    comparison past that point would score an unwritten zero against a live
    setback and report a difference that is really the end of the record.

    Args:
        series_a: road_series() for the free-running arm.
        series_b: road_series() for the prescribed arm.
        targets: {gis: historical_event_year}.
        start_year: First calendar year of the run.
        check_2004: {gis: measured_setback_m}, the independent cross-check.

    Returns:
        A (summary_df, long_df) tuple. long_df is one row per domain-year and
        is what the GIF and any trajectory plot read; it carries a `managed`
        flag per arm so a plot can stop the line where the record stops.
    """
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
            measured_2004_m=check_2004.get(gis),
        ))
    return pd.DataFrame(summary), pd.DataFrame(long)


def check_determinism(series_a, series_b, first_event_year, start_year):
    """Confirms the two arms are identical before the first prescribed event.

    Both arms are the same configuration up to the first event, and Barrier3D
    runs on a seeded RNG with a prescribed storm file, so every series must
    agree exactly through `first_event_year - 1`. A divergence there is a bug,
    not a result, and everything downstream would be uninterpretable -- so
    this is checked rather than assumed.

    Args:
        series_a: road_series() for the free-running arm.
        series_b: road_series() for the prescribed arm.
        first_event_year: Calendar year of the earliest prescribed event.
        start_year: First calendar year of the run.

    Returns:
        A (ok, offenders) tuple; offenders lists the GIS domains that diverged.
    """
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


def check_event_indexing(series_b, targets, start_year):
    """Confirms index i really is calendar year start_year + i.

    Arm B's setback is displaced by a measured amount at the event year and by
    dune migration alone in every other year, so the largest single-year jump
    in the prescribed arm must land on the event year. If it does not, the
    time indexing in this file is wrong and every year reported here is off.

    Args:
        series_b: road_series() for the prescribed arm.
        targets: {gis: historical_event_year}.
        start_year: First calendar year of the run.

    Returns:
        A DataFrame with the located jump year beside the expected one.
    """
    rows = []
    for gis, event_year in sorted(targets.items()):
        entry = series_b.get(gis)
        if entry is None:
            continue
        setback = entry["setback"]
        # Bound the search by the last MANAGED year, dated from the elevation
        # series -- see last_managed_index. Unwritten trailing years would
        # otherwise contribute a spurious jump back to zero.
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


# =============================================================================
# Report
# =============================================================================

class _Tee:
    """Mirrors everything printed to the console into a buffer.

    The report is the run's own console output rather than a second rendering
    of the same numbers, which is the point: a separately-composed summary can
    disagree with the CSVs beside it, and this cannot.
    """

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


def _arm_provenance(run_dir):
    """One line describing which run a comparison arm actually read."""
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


def _report_header(arm_a, arm_b, preset):
    """Provenance block written above the captured output.

    WHY THIS EXISTS. On 2026-08-25 a re-run of this comparison rewrote every
    CSV and GIF in output/comparisons/relocation_1984_2004/<preset>/ and left
    the report.txt from 2026-08-22 sitting beside them -- the script had lost
    its report-writing step, so nothing overwrote it. For three days that
    folder held a report describing DIFFERENT runs from the CSVs next to it,
    with nothing on its face to say so. It even named the pre-restructure flat
    run paths, which by then did not exist.

    So every report now carries the identity of the two runs it was built
    from. A report whose arms do not match the runs on disk is visible at a
    glance instead of having to be inferred from file mtimes.
    """
    return (
        "=" * 74 + "\n"
        f"generated   {datetime.datetime.now():%Y-%m-%d %H:%M:%S} by "
        f"{os.path.basename(__file__)}\n"
        f"preset      {preset}\n"
        f"arm A       {_arm_provenance(arm_a)}\n"
        f"arm B       {_arm_provenance(arm_b)}\n"
        "\n"
        "This file is the console output of the run that wrote the CSVs and\n"
        "GIFs beside it. If the arm identities above do not match the runs on\n"
        "disk, this report is stale -- re-run the comparison.\n"
        + "=" * 74 + "\n\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    raw_runs = PROJECT_BASE_DIR / "output" / "raw_runs"
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
                             "(default OUTPUT_ROOT/<preset>)")
    args = parser.parse_args()

    preset, _ = resolve_be_preset(args.preset)
    name_a, name_b = arm_names(preset)
    # Runs are filed [<forcing arm>/]<period>/<preset>/. Both relocation arms
    # share a preset by construction -- arm_names derives both from the one
    # token -- so the directory is resolved once and used for both. Resolved
    # rather than joined: the join had no slot for the forcing-arm component,
    # and "arm" here means the relocation switch, not that one.
    period_dir = preset_dir_for(raw_runs, (START_YEAR, END_YEAR), preset)
    args.arm_a = args.arm_a or str(period_dir / name_a)
    args.arm_b = args.arm_b or str(period_dir / name_b)

    out_dir = Path(args.out).resolve() if args.out else (OUTPUT_ROOT / preset)
    out_dir.mkdir(parents=True, exist_ok=True)

    # A STALE report is worse than a missing one -- see _report_header. Delete
    # any previous report BEFORE doing the work, so a run that dies partway
    # leaves this folder with no report rather than with the previous run's.
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


def _compare(args, preset, out_dir):
    """The comparison itself. Everything it prints becomes report.txt."""
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

    print(f"\n  managed road domains  A: {len(series_a)}   B: {len(series_b)}")
    print(f"  historical domains    {sorted(targets)}")
    for year in sorted(set(targets.values())):
        moved = sorted(g for g, y in targets.items() if y == year)
        print(f"    {year}  GIS {moved}")

    # --- 0. the checks that decide whether the rest means anything ----------
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
        idx_df.to_csv(out_dir / "indexing_check.csv", index=False)

    # --- 1. first relocation year -------------------------------------------
    print("\n" + "-" * 74)
    print("1. FIRST MODELLED RELOCATION, free-running arm")
    print("-" * 74)
    first_df = score_first_year(series_a, targets, START_YEAR)
    print(first_df.to_string(index=False))
    first_df.to_csv(out_dir / "first_relocation_year.csv", index=False)

    dated = first_df.dropna(subset=["error_years"])
    if not dated.empty:
        err = dated["error_years"].astype(float)
        print(f"\n  domains that relocated at all   {len(dated)}/{len(first_df)}")
        print(f"  median signed error             {err.median():+.1f} yr")
        print(f"  mean absolute error             {err.abs().mean():.1f} yr")
    else:
        print("\n  no historical domain relocated in the free-running arm")

    # --- 1b. how close did the misses come? ---------------------------------
    print("\n" + "-" * 74)
    print("1b. NEAR MISSES: how much more dune migration would have fired it")
    print("-" * 74)
    near = near_miss_table(series_a, targets, START_YEAR)
    near.to_csv(out_dir / "near_miss_margin.csv", index=False)
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

    # --- 2. hit / miss --------------------------------------------------------
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
        out_dir / "confusion.csv", index=False)

    # --- 3. setback trajectories ----------------------------------------------
    print("\n" + "-" * 74)
    print("3. SETBACK TRAJECTORIES")
    print("-" * 74)
    traj_df, long_df = score_trajectories(
        series_a, series_b, targets, START_YEAR, HATTERAS_RELOCATION_CHECK_2004)
    hist = traj_df[traj_df["historical"]]
    print(hist.to_string(index=False))
    traj_df.to_csv(out_dir / "setback_summary.csv", index=False)
    long_df.to_csv(out_dir / "setback_by_year.csv", index=False)
    print(f"\n  saved per-year setbacks for {long_df['gis'].nunique()} domains "
          f"-> setback_by_year.csv")

    # --- 4. outcomes -----------------------------------------------------------
    print("\n" + "-" * 74)
    print("4. ROAD OUTCOMES: did prescribing the relocations change the fate?")
    print("-" * 74)
    out_rows = []
    for label, cas in (("free", cascade_a), ("prescribed", cascade_b)):
        for row in roadway_module.summarise_road_management(
                cas, HATTERAS_DOMAINS, *span):
            out_rows.append(dict(arm=label, **row))
    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(out_dir / "road_outcomes.csv", index=False)

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

    # --- 5. animation --------------------------------------------------------
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
        make_all_road_gifs(
            (shore_a, info_a), (shore_b, info_b), series_a, series_b,
            GIF_WINDOWS, str(out_dir), event_years=targets,
            back_a=back_a, back_b=back_b, gif_config=GIF_CONFIG)
        for name, lo, hi in TOPO_WINDOWS:
            make_topography_gif(
                cascade_a, cascade_b, series_a, series_b, lo, hi,
                os.path.join(str(out_dir), f"road_topography_{_slug(name)}.gif"),
                START_YEAR, event_years=targets, gif_config=GIF_CONFIG,
                title=f"Hatteras topography and NC-12 \u2014 {name}",
                planform_note=PLANFORM_NOTE)

    print("\n" + "=" * 74)
    print(f"artifacts -> {out_dir}")
    print("=" * 74)


if __name__ == "__main__":
    main()
