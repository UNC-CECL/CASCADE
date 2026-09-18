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
#         duneline_vs_coastsat/       dune-line change vs the CoastSat
#                                     shoreline, one folder per window
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
#     from site_layer.hat_observed_rates import lrr_csv, transect_lookup
#     path = lrr_csv(1996, 2010)          # raises, listing windows, if absent
# ==============================================================================

from __future__ import annotations

from pathlib import Path

# RUN AS A FILE, NOT IMPORTED. `python scripts/site_layer/hat_observed_rates.py` puts this
# file's OWN folder on sys.path, not scripts/, so a `site_layer.` import cannot
# resolve and the __main__ block below would die on it. Importing the module
# the normal way never takes this branch -- __package__ is "site_layer" then.
if __package__ in (None, ""):
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents
                    if (_p / "pyproject.toml").exists())
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"

SCR_ROOT = INIT_ROOT / "5-scr"

COASTSAT_TIMESERIES = SCR_ROOT / "coastsat_timeseries"
COASTSAT_LRR_ROOT = SCR_ROOT / "coastsat_lrr"
TRANSECT_DOMAINS = SCR_ROOT / "transect_domains"
TIMESERIES_LRR = SCR_ROOT / "coastsat_timeseries_lrr"
DUNELINE_VS_COASTSAT = SCR_ROOT / "duneline_vs_coastsat"
# The dune-line LRR product (2026-09-16): the same layout as coastsat_lrr/,
# one folder per window, an OLS through every island-wide dune line inside
# the window per transect. Written by
# scripts/input_prep/5-scr/duneline_lrr/duneline_lrr.py.
DUNELINE_LRR_ROOT = SCR_ROOT / "duneline_lrr"
SHORELINE_INVENTORY = SCR_ROOT / "shoreline_inventory"
SHORELINE_PATTERNS = SCR_ROOT / "shoreline_change_patterns"
DSAS_ROOT = SCR_ROOT / "scr-dsas-1978-2019"
# The four windows drawn on one y axis (coastsat_lrr_windows.py).
COASTSAT_LRR_WINDOWS = SCR_ROOT / "coastsat_lrr_windows"
# Two-window comparison figures (coastsat_two_period_comparison.py).
TWO_PERIOD_COMPARISON = COASTSAT_LRR_ROOT / "two_period_comparison"
# Retired windows (1978-1997, 1997-2019 and their specific-dates variants),
# kept for the DSAS comparison in 6-scr-smooth, never a grading target.
COASTSAT_LRR_SUPERSEDED = COASTSAT_LRR_ROOT / "superseded_20260810"

# Single files in transect_domains/ that scripts outside 5-scr read by name.
# HAT_domains.json holds the 90 real domain boxes; the map_elements polygons
# in 9-figures are NOT the model domains.
DOMAIN_BOXES = TRANSECT_DOMAINS / "HAT_domains.json"
TRANSECT_LAYER = TRANSECT_DOMAINS / "CoastSat_transect_layer.geojson"

# The per-window products a rate fit writes.
TRANSECT_FILE = "transect_lrr_full.csv"
DOMAIN_FILE = "domain_lrr_summary.csv"

# Folders under coastsat_lrr/ that are not windows. Listed so windows() can
# report what IS available without having to parse every directory name.
_NOT_WINDOWS = {"old_time_periods", "two_period_comparison",
                "rodanthe_plots", "old_dsas_comparisons",
                "superseded_20260810", "custom"}


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


def dune_lrr_csv(start_year, end_year):
    """Per-transect OLS through the island-wide dune lines of one window,
    in the coastsat_lrr layout. Raises, naming the producer, if absent."""
    path = DUNELINE_LRR_ROOT / "{0}_{1}".format(start_year, end_year) / TRANSECT_FILE
    if not path.is_file():
        have = sorted(p.name for p in DUNELINE_LRR_ROOT.glob("*_*")
                      if (p / TRANSECT_FILE).is_file()) if DUNELINE_LRR_ROOT.is_dir() else []
        raise FileNotFoundError(
            "no dune-line LRR for {0}-{1}; have {2}. Build it with "
            "scripts/input_prep/5-scr/duneline_lrr/duneline_lrr.py".format(
                start_year, end_year, have or "none"))
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


# THE EXTENSION (2026-09-16, the Pea Island extension experiment). The
# transects beyond GIS 1-90, numbered by their whole-island domain polygon
# (hat_extension_domains.join_origins), get their own lookup and their own rate
# fit per window, beside the surveyed ones and never merged into them:
#
#     transect_domains/transect_domain_lookup_ext.csv
#     coastsat_lrr/<window>/ext/transect_lrr_full.csv       extension only
#     coastsat_lrr/<window>/ext/transect_lrr_with_base.csv  surveyed + extension
#
# The with_base file is what an extended-geometry run loads as its active
# dataset (one LOESS over the whole reach); the surveyed file stays the
# scoring table for GIS 2-89 so extended and base runs are graded alike.
# Built by scripts/input_prep/5-scr/CoastSat/coastsat_extension_lrr.py.
EXT_DIR = "ext"
WITH_BASE_FILE = "transect_lrr_with_base.csv"


def transect_lookup_ext():
    path = TRANSECT_DOMAINS / "transect_domain_lookup_ext.csv"
    if not path.is_file():
        raise FileNotFoundError(
            "no extension lookup at {0}; build it with "
            "coastsat_extension_lrr.py".format(path))
    return path


def lrr_csv_ext(start_year, end_year, with_base=True):
    """The extension's rate fit for one window (with the surveyed reach by
    default), raising if the extension has not been built for it."""
    path = (window_dir(start_year, end_year) / EXT_DIR
            / (WITH_BASE_FILE if with_base else TRANSECT_FILE))
    if not path.is_file():
        raise FileNotFoundError(
            "no extension rates for {0}-{1} at {2}; build them with "
            "coastsat_extension_lrr.py --start-year {0} --end-year {1}"
            .format(start_year, end_year, path))
    return path
