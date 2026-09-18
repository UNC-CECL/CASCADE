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
# THE LAYOUT, grouped by job (2026-09-18). Each group is numbered in the
# order the work runs: observations, then the frame that ties them to
# domains, then the rate fits built on both, then the comparisons.
#
#     data/hatteras_init/5-scr/
#         1-observations/             measured or digitized, not fitted by us
#             coastsat_timeseries/    raw per-transect chainage, by CoastSat site
#             dsas_1978_2019/         the DSAS rates, a different source
#             shoreline_inventory/    study-area and reference shorelines
#         2-transect-frame/
#             transect_domains/       the lookup, the transect layer, the domain
#                                     polygons, and the verification set
#         3-rates/                    the MODEL TARGETS: tables + one figure per
#                                     window (rates_figures.py).
#                                     Grouped by source since 2026-09-18.
#             coastsat/
#                 lrr/<window>/       the OLS rate fits, one folder per WINDOW:
#                                     1984_2004 1996_2010 2004_2024 2010_2024,
#                                     and 1996_2024 (CONTEXT, not graded)
#                 endpoint/<window>/  net change between +/-6-month means at
#                                     the dune-line dates (m and m/yr)
#                 5yr_bins/<window>/  the OLS in successive 5-year bins,
#                                     1996_2010 2010_2024 1996_2024
#             duneline/
#                 endpoint/<window>/  net change between the two dune lines
#                                     (m and m/yr; replaced duneline_lrr/)
#         4-comparisons/
#             coastsat_windows/       the four windows on one y axis
#             duneline_vs_coastsat/   dune-line change vs the CoastSat
#                                     shoreline, one folder per window
#             duneline_windows/       net dune-line change in metres, a
#                                     long window and its halves
#             net_change_1996_2024/   CoastSat vs dune-line net change,
#                                     1996-2024 and its halves (09-18)
#             trajectory_patterns/    trajectory classification output
#             two_period_comparison/  1984-2004 against 2004-2024
#         archive/                    retired windows, old 5-year-bin runs,
#                                     the Rodanthe poster figures; not for use
#
# Before 2026-09-18 all of these sat side by side at the top of 5-scr/, under
# their old names (scr-dsas-1978-2019, coastsat_timeseries_lrr,
# coastsat_lrr_windows, shoreline_change_patterns).
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
OBSERVATIONS = SCR_ROOT / "1-observations"
TRANSECT_FRAME = SCR_ROOT / "2-transect-frame"
RATES = SCR_ROOT / "3-rates"
COMPARISONS = SCR_ROOT / "4-comparisons"
ARCHIVE = SCR_ROOT / "archive"

COASTSAT_TIMESERIES = OBSERVATIONS / "coastsat_timeseries"
# 3-rates is grouped BY SOURCE since 2026-09-18 (Hannah): coastsat/{lrr,
# endpoint,5yr_bins} and duneline/endpoint. Until then the four sat flat as
# coastsat_lrr/, coastsat_endpoint/, coastsat_5yr_bins/, duneline_endpoint/.
# 3-rates holds the tables and ONE house-style figure per window beside them
# (scripts/input_prep/5-scr/rates_figures.py); comparisons are in 4-comparisons.
COASTSAT_RATES = RATES / "coastsat"
DUNELINE_RATES = RATES / "duneline"
COASTSAT_LRR_ROOT = COASTSAT_RATES / "lrr"
TRANSECT_DOMAINS = TRANSECT_FRAME / "transect_domains"
TIMESERIES_LRR = COASTSAT_RATES / "5yr_bins"
DUNELINE_VS_COASTSAT = COMPARISONS / "duneline_vs_coastsat"
# The stored dune-line observation (2026-09-18): the NET CHANGE between the
# two lines that bound a window, per transect and per domain, in m and m/yr,
# one folder per window. Written by
# scripts/input_prep/5-scr/duneline_endpoint/duneline_endpoint.py. It replaced
# duneline_lrr/ (an OLS through every line in the window, 2026-09-16; Hannah:
# "we are tracking net change"), now under archive/duneline_lrr_retired_20260918/.
DUNELINE_ENDPOINT_ROOT = DUNELINE_RATES / "endpoint"
# The CoastSat counterpart (2026-09-18): net change in the CoastSat shoreline
# between +/-6-month window means centred on the SAME dune-line survey dates,
# so the two products difference like for like. Written by
# scripts/input_prep/5-scr/coastsat_endpoint/coastsat_endpoint.py.
COASTSAT_ENDPOINT_ROOT = COASTSAT_RATES / "endpoint"
# Both endpoint products use the same two file names.
ENDPOINT_TRANSECT_FILE = "transect_endpoint.csv"
ENDPOINT_DOMAIN_FILE = "domain_endpoint_summary.csv"
DUNE_ENDPOINT_TRANSECT_FILE = ENDPOINT_TRANSECT_FILE
DUNE_ENDPOINT_DOMAIN_FILE = ENDPOINT_DOMAIN_FILE
# The CoastSat and dune-line net change side by side, 1996-2024 and its
# halves (net_change_1996_2024.py, 2026-09-18).
NET_CHANGE_1996_2024 = COMPARISONS / "net_change_1996_2024"
SHORELINE_INVENTORY = OBSERVATIONS / "shoreline_inventory"
SHORELINE_PATTERNS = COMPARISONS / "trajectory_patterns"
DSAS_ROOT = OBSERVATIONS / "dsas_1978_2019"
# The four windows drawn on one y axis (coastsat_lrr_windows.py).
COASTSAT_LRR_WINDOWS = COMPARISONS / "coastsat_windows"
# Net dune-line change over a long window and its two halves, in metres
# (duneline_windows.py, 2026-09-18): the endpoint mirror of the CoastSat
# halves figure, one folder per window.
DUNELINE_WINDOWS = COMPARISONS / "duneline_windows"
# Two-window comparison figures (coastsat_two_period_comparison.py).
TWO_PERIOD_COMPARISON = COMPARISONS / "two_period_comparison"
# Retired windows (1978-1997, 1997-2019 and their specific-dates variants),
# kept for the DSAS comparison in 6-scr-smooth, never a grading target.
# Archived, but still read, so still resolved.
COASTSAT_LRR_SUPERSEDED = ARCHIVE / "coastsat_lrr_superseded_20260810"

# 6-scr-smooth: what the LOESS smoothing does to the observed rates. Resolved
# here too because its outputs are read outside their producers (2026-09-18,
# when the two folders lost their HAT_*_output names).
#     method_comparison/   transect-based against domain-averaged smoothing;
#                          03_cascade_inputs/ is read by the overwash work
#     dsas_vs_coastsat/    the two rate sources, both smoothed, on the
#                          retired 1978-1997 / 1997-2019 windows
SMOOTH_ROOT = INIT_ROOT / "6-scr-smooth"
SMOOTH_METHOD_COMPARISON = SMOOTH_ROOT / "method_comparison"
SMOOTH_DSAS_VS_COASTSAT = SMOOTH_ROOT / "dsas_vs_coastsat"

# Single files in transect_domains/ that scripts outside 5-scr read by name.
# HAT_domains.json holds the 90 real domain boxes; the older 1000 m polygons,
# now in map_elements/archive/, are NOT the model domains.
DOMAIN_BOXES = TRANSECT_DOMAINS / "HAT_domains.json"
TRANSECT_LAYER = TRANSECT_DOMAINS / "CoastSat_transect_layer.geojson"

# The per-window products a rate fit writes.
TRANSECT_FILE = "transect_lrr_full.csv"
DOMAIN_FILE = "domain_lrr_summary.csv"

# Folders under coastsat_lrr/ that are not windows. Listed so windows() can
# report what IS available without having to parse every directory name.
# Everything but "custom" moved out on 2026-09-18; the names stay listed so a
# stray copy restored from an old checkout is not mistaken for a window.
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


def coastsat_endpoint_csv(start_year, end_year, level="transect"):
    """Net CoastSat shoreline change for one window, between +/-6-month means
    about the dune-line survey dates: `level` "transect" or "domain". Raises,
    naming the producer, if absent."""
    name = {"transect": ENDPOINT_TRANSECT_FILE, "domain": ENDPOINT_DOMAIN_FILE}[level]
    path = COASTSAT_ENDPOINT_ROOT / "{0}_{1}".format(start_year, end_year) / name
    if not path.is_file():
        raise FileNotFoundError(
            "no CoastSat endpoint for {0}-{1}. Build it with "
            "scripts/input_prep/5-scr/coastsat_endpoint/coastsat_endpoint.py".format(
                start_year, end_year))
    return path


def dune_endpoint_csv(start_year, end_year, level="transect"):
    """Net dune-line change for one window: `level` "transect" (per 100 m
    transect) or "domain" (per GIS domain). Raises, naming the producer, if
    absent."""
    name = {"transect": DUNE_ENDPOINT_TRANSECT_FILE,
            "domain": DUNE_ENDPOINT_DOMAIN_FILE}[level]
    path = DUNELINE_ENDPOINT_ROOT / "{0}_{1}".format(start_year, end_year) / name
    if not path.is_file():
        have = sorted(p.name for p in DUNELINE_ENDPOINT_ROOT.glob("*_*")
                      if (p / name).is_file()) if DUNELINE_ENDPOINT_ROOT.is_dir() else []
        raise FileNotFoundError(
            "no dune-line endpoint for {0}-{1}; have {2}. Build it with "
            "scripts/input_prep/5-scr/duneline_endpoint/duneline_endpoint.py".format(
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
