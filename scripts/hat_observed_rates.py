# ==============================================================================
# hat_observed_rates.py
#
# WHERE IS THE OBSERVED SHORELINE DATA - the CoastSat chainage, the per-window
# rate fits the model is graded against, and the transect-to-domain lookup that
# ties them to Barrier3D domains?
#
# WHY THIS EXISTS
#   Until 2026-09-12 these lived under `scripts/input_prep/5-scr/CoastSat/`.
#   They are DATA, and the rate CSV is a MODEL INPUT: section 8 of the hindcast
#   runner reads it on every run and the calibrated source/sink preset is fitted
#   against it. Everything else the model ingests comes from
#   `data/hatteras_init/`, so one input lived somewhere no one would look for
#   it, and anyone archiving or sharing the data tree shipped a model that
#   could not run.
#
#   Twenty-two files built that path by hand. Moving it once meant editing all
#   of them, which is the failure `hat_topo_version.py` was written to end for
#   topography. Same answer here: the location is resolved ONCE, and a window
#   that is not on disk is a loud error naming the ones that are.
#
# THE LAYOUT, which mirrors the producer folders so provenance is obvious
#
#     data/hatteras_init/5-scr/
#         coastsat_timeseries/        raw per-transect chainage, by CoastSat site
#         transect_domains/           the lookup, the transect layer, the domain
#                                     polygons, and the verification set
#         coastsat_lrr/               rate fits, one folder per WINDOW
#             1984_2004/  1996_2010/  2004_2024/  2010_2024/
#             old_time_periods/       retired windows, not for use
#             two_period_comparison/  rodanthe_plots/  old_dsas_comparisons/
#         coastsat_timeseries_lrr/    the 5-year-bin fits
#         shoreline_inventory/        study-area and reference shorelines
#         shoreline_change_patterns/  trajectory classification output
#         scr-dsas-1978-2019/         the DSAS rates, a different source
#
# A WINDOW IS <start>_<end>, NOT A START YEAR. A rate fit spans an interval, so
# it is named for one - unlike a dune line or a road alignment, which is a
# survey at a moment and is named for its year. That split is the naming rule
# the data tree follows (Hannah, 2026-09-12).
#
# USAGE
#     from hat_observed_rates import lrr_csv, transect_lookup
#     path = lrr_csv(1996, 2010)          # raises, listing windows, if absent
# ==============================================================================

from __future__ import annotations

from pathlib import Path

_HERE = Path(__file__).resolve()
PROJECT_ROOT = _HERE.parents[1]
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"

SCR_ROOT = INIT_ROOT / "5-scr"

COASTSAT_TIMESERIES = SCR_ROOT / "coastsat_timeseries"
COASTSAT_LRR_ROOT = SCR_ROOT / "coastsat_lrr"
TRANSECT_DOMAINS = SCR_ROOT / "transect_domains"
TIMESERIES_LRR = SCR_ROOT / "coastsat_timeseries_lrr"
SHORELINE_INVENTORY = SCR_ROOT / "shoreline_inventory"
SHORELINE_PATTERNS = SCR_ROOT / "shoreline_change_patterns"
DSAS_ROOT = SCR_ROOT / "scr-dsas-1978-2019"

# The per-window products a rate fit writes.
TRANSECT_FILE = "transect_lrr_full.csv"
DOMAIN_FILE = "domain_lrr_summary.csv"

# Folders under coastsat_lrr/ that are not windows. Listed so windows() can
# report what IS available without having to parse every directory name.
_NOT_WINDOWS = {"old_time_periods", "two_period_comparison",
                "rodanthe_plots", "old_dsas_comparisons"}


def windows():
    """Every rate window on disk, as (start, end) pairs, oldest first."""
    if not COASTSAT_LRR_ROOT.is_dir():
        return []
    found = []
    for path in COASTSAT_LRR_ROOT.iterdir():
        if not path.is_dir() or path.name in _NOT_WINDOWS:
            continue
        start, _, end = path.name.partition("_")
        if start.isdigit() and end.isdigit():
            found.append((int(start), int(end)))
    return sorted(found)


def window_dir(start_year, end_year):
    """The folder holding one window's rate fit.

    Raises:
        FileNotFoundError: If that window has not been built, naming the ones
            that have. A window is built by
            scripts/input_prep/5-scr/CoastSat/coastsat_domain_lrr_fixed.py.
    """
    path = COASTSAT_LRR_ROOT / "{0}_{1}".format(start_year, end_year)
    if not path.is_dir():
        raise FileNotFoundError(
            "no observed rates for {0}-{1}. On disk: {2}. Build it with "
            "coastsat_domain_lrr_fixed.py --start-year {0} --end-year {1}"
            .format(start_year, end_year,
                    ", ".join("{0}-{1}".format(*w) for w in windows())
                    or "none"))
    return path


def lrr_csv(start_year, end_year):
    """Per-transect rate fit for one window -- the model's grading target."""
    path = window_dir(start_year, end_year) / TRANSECT_FILE
    if not path.is_file():
        raise FileNotFoundError(
            "{0} exists but holds no {1}".format(path.parent, TRANSECT_FILE))
    return path


def domain_csv(start_year, end_year):
    """Per-domain summary of one window's rate fit."""
    path = window_dir(start_year, end_year) / DOMAIN_FILE
    if not path.is_file():
        raise FileNotFoundError(
            "{0} exists but holds no {1}".format(path.parent, DOMAIN_FILE))
    return path


def transect_lookup():
    """Transect id -> Barrier3D domain, written by coastsat_domain_mapping.py."""
    path = TRANSECT_DOMAINS / "transect_domain_lookup.csv"
    if not path.is_file():
        raise FileNotFoundError(
            "no transect-to-domain lookup at {0}".format(path))
    return path


if __name__ == "__main__":
    print("5-scr root       {0}".format(SCR_ROOT))
    print("rate windows     {0}".format(
        ", ".join("{0}-{1}".format(*w) for w in windows()) or "none"))
    print("transect lookup  {0}".format(
        "ok" if (TRANSECT_DOMAINS / "transect_domain_lookup.csv").is_file()
        else "MISSING"))
