# ==============================================================================
# HAT_period_input_check.py
#
# Does every hindcast period resolve everything it needs, and does each of
# those things exist on disk?
#
# WHY THIS EXISTS
#   HATTERAS_PERIODS names six forcing inputs per period as RELATIVE PATHS.
#   Nothing checks them until a run is most of the way through section 5, and
#   two of the six are read later still -- so a period wired against a file
#   that was never built fails partway into a run that has already spent
#   minutes building an island. Four periods since 2026-09-11, two of them
#   wired ahead of inputs still being digitised, made that a certainty rather
#   than a risk.
#
#   It is also the answer to "is this period runnable yet", which is otherwise
#   answered by starting a run and waiting.
#
# WHAT IT DOES NOT DO
#   It does not validate VALUES. A setback of the right shape measured against
#   the wrong topography is a real failure this cannot see; that is what
#   HAT_road_setback_audit.py and the extractor audits are for. This checks
#   resolution and existence, the class of failure that costs a run rather
#   than a result.
#
# EXIT CODE
#   0 when every period is runnable, 1 when any is blocked, so a batch driver
#   can gate on it.
#
#     python HAT_period_input_check.py
#     python HAT_period_input_check.py --period 2010
#
# Author: Hannah A. Henry, UNC CECL
# ==============================================================================

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
PROJECT_ROOT = _HERE.parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from hatteras_site_config import (                      # noqa: E402
    HATTERAS_PERIODS, HATTERAS_BE_PRESETS, HATTERAS_DOMAINS,
    HATTERAS_ROAD_EVENTS, HATTERAS_NOURISHMENT_PROJECTS,
    HATTERAS_ROAD_ELEVATION_FILE, HATTERAS_FIRST_ROAD_DOMAIN,
    HATTERAS_LAST_ROAD_DOMAIN)
from hat_topo_version import domain_arrays, topo_dirs   # noqa: E402
from hat_observed_rates import (                        # noqa: E402
    COASTSAT_LRR_ROOT, TRANSECT_FILE)

INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"

# Resolved, not built: hat_observed_rates owns where the rate fits live, so
# this check cannot look somewhere the runner does not.

OK = "ok"
MISSING = "MISSING"


class Finding:
    """One checked input: what it is, whether it is there, what it holds."""

    def __init__(self, label, status, detail="", blocking=True):
        self.label = label
        self.status = status
        self.detail = detail
        self.blocking = blocking

    @property
    def failed(self):
        return self.status != OK

    def line(self):
        mark = "  " if not self.failed else ("!!" if self.blocking else " ~")
        return "{0} {1:<22} {2:<9} {3}".format(
            mark, self.label, self.status, self.detail)


def check_file(label, path, blocking=True):
    if not path.is_file():
        return Finding(label, MISSING, str(path), blocking)
    return Finding(label, OK, path.name, blocking)


def check_storms(path, run_years):
    """Exists, loads, and covers every model step the run will spend."""
    found = check_file("storm series", path)
    if found.failed:
        return found
    storms = np.load(path)
    steps = storms[:, 0].astype(int)
    found.detail = "{0} storms, steps {1}-{2}".format(
        len(storms), steps.min(), steps.max())
    if steps.max() < run_years:
        found.status = "SHORT"
        found.detail += " but the run spends {0} steps".format(run_years)
    return found


def check_offsets(path, geometry):
    """Exists, and is the padded length CASCADE expects."""
    found = check_file("island offsets", path)
    if found.failed:
        return found
    values = np.loadtxt(path, delimiter=",", skiprows=1).reshape(-1)
    found.detail = "{0} domains".format(values.size)
    if values.size != geometry.total_domains:
        found.status = "SHAPE"
        found.detail += ", expected {0}".format(geometry.total_domains)
    return found


def check_setbacks(path):
    """Exists, is the 2-row format, and covers the road span."""
    found = check_file("road setback", path)
    if found.failed:
        return found
    raw = np.loadtxt(path, delimiter=",")
    expected = HATTERAS_LAST_ROAD_DOMAIN - HATTERAS_FIRST_ROAD_DOMAIN + 1
    found.detail = "{0} x {1}".format(*raw.shape)
    if raw.ndim != 2 or raw.shape[0] != 2 or raw.shape[1] != expected:
        found.status = "SHAPE"
        found.detail += ", expected 2 x {0}".format(expected)
        return found
    # Not a failure, but worth seeing: a zero setback puts the road on the dune
    # line, where one cell of retreat re-fires a relocation.
    found.detail += ", {0} at zero".format(int((raw[1] == 0).sum()))
    return found


def check_topography(product, geometry):
    """The product resolves to a version, and its domain arrays are on disk."""
    try:
        _, _, version = topo_dirs(product)
    except Exception as exc:                 # the resolver raises loudly
        return Finding("topography", "UNRESOLVED",
                       "{0}: {1}".format(product, exc))

    elev, dunes = domain_arrays(
        product, first_gis=geometry.first_gis_id, last_gis=geometry.last_gis_id)
    absent = [p for p in elev + dunes if not Path(p).is_file()]
    found = Finding("topography", OK, "{0} / {1}, {2} arrays".format(
        product, version, len(elev) + len(dunes)))
    if absent:
        found.status = MISSING
        found.detail = "{0} / {1}, {2} of {3} arrays absent".format(
            product, version, len(absent), len(elev) + len(dunes))
    return found


def check_presets(start_year):
    """Which source/sink presets are solved for this period.

    Not blocking: zeroBE alone is enough to run. A period missing edgeBE simply
    cannot be run under it, which is a fact about the calibration rather than
    about the inputs.
    """
    solved = sorted(name for name, rates in HATTERAS_BE_PRESETS.items()
                    if start_year in rates)
    missing = sorted(set(HATTERAS_BE_PRESETS) - set(solved))
    found = Finding("source/sink", OK, ", ".join(solved), blocking=False)
    if missing:
        found.status = "PARTIAL"
        found.detail = "{0}  (unsolved: {1})".format(
            ", ".join(solved), ", ".join(missing))
    return found


def events_in_window(start_year, end_year):
    """What the record fires inside this window.

    The window is start..end-1, matching run_cascade_simulation: an event dated
    exactly on the end year belongs to the next period and never fires here.
    """
    lines = []
    for event in HATTERAS_ROAD_EVENTS:
        year = getattr(event, "year", None)
        if year is None or not (start_year <= year < end_year):
            continue
        moved = getattr(event, "displacement_m", None)
        if moved:
            lines.append("{0}  road relocation, GIS {1}".format(
                year, sorted(moved)))
        else:
            lines.append("{0}  {1}".format(
                year, getattr(event, "note", "road event")))
    for project in HATTERAS_NOURISHMENT_PROJECTS:
        if start_year <= project.year < end_year:
            lines.append("{0}  {1}".format(project.year, project.name))
    return sorted(lines)


def check_period(start_year):
    period = HATTERAS_PERIODS[start_year]
    end_year = period["end_year"]
    run_years = end_year - start_year

    print("\n" + "=" * 78)
    print("{0}-{1}   {2} model years, simulating {0}-{3}".format(
        start_year, end_year, run_years, end_year - 1))
    print("=" * 78)

    findings = [
        check_topography(period["topo_product"], HATTERAS_DOMAINS),
        check_storms(INIT_ROOT / period["storm_file"], run_years),
        check_offsets(INIT_ROOT / period["island_offset_file"],
                      HATTERAS_DOMAINS),
        check_setbacks(INIT_ROOT / period["road_setback_file"]),
        check_file("road elevation", INIT_ROOT / HATTERAS_ROAD_ELEVATION_FILE),
        check_file("observed rates",
                   COASTSAT_LRR_ROOT / "{0}_{1}".format(start_year, end_year)
                   / TRANSECT_FILE),
        check_presets(start_year),
    ]
    for found in findings:
        print(found.line())

    fills = "on, {0} m3/m".format(period["nourishment_volume"]) \
        if period["enable_nourishment"] else "off"
    print("   {0:<22} {1}".format("nourishment", fills))
    print("   {0:<22} {1} m/yr".format("sea level",
                                       period["sea_level_rise_rate"]))

    events = events_in_window(start_year, end_year)
    print("   {0:<22} {1}".format(
        "events in window", events[0] if events else "none"))
    for line in events[1:]:
        print("   {0:<22} {1}".format("", line))

    blocked = [f for f in findings if f.failed and f.blocking]
    if blocked:
        print("\n   BLOCKED: " + ", ".join(f.label for f in blocked))
    else:
        print("\n   runnable")
    return not blocked


def main():
    parser = argparse.ArgumentParser(
        description="resolve and check the inputs of every hindcast period")
    parser.add_argument("--period", type=int, action="append",
                        help="start year; repeatable. Default: all of them.")
    args = parser.parse_args()

    periods = args.period or sorted(HATTERAS_PERIODS)
    unknown = [p for p in periods if p not in HATTERAS_PERIODS]
    if unknown:
        parser.error("no such period {0}; have {1}".format(
            unknown, sorted(HATTERAS_PERIODS)))

    results = {p: check_period(p) for p in periods}

    print("\n" + "=" * 78)
    ready = [p for p, good in results.items() if good]
    blocked = [p for p, good in results.items() if not good]
    print("runnable: {0}".format(ready if ready else "none"))
    if blocked:
        print("blocked : {0}".format(blocked))
    print("=" * 78)
    return 1 if blocked else 0


if __name__ == "__main__":
    sys.exit(main())
