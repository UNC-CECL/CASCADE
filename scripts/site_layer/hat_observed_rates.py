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
#                 total_change/       the same rate as a DISTANCE: LRR(W) x the
#                     <window>/       years of W, beside the observed change.
#                                     1996_2024 1996_2010 2010_2024
#                                     (was lrr_projected/ until 2026-09-21 --
#                                     nothing in it was ever a projection)
#                 projected/          the 1996-2024 LRR carried onto a window it
#                     <window>/       was NOT fitted on: x 14 yr over 1996_2010
#                                     and 2010_2024. The only projections here.
#             duneline/
#                 endpoint/<window>/  net change between the two dune lines
#                                     (m and m/yr; replaced duneline_lrr/)
#         4-comparisons/
#             shoreline_vs_duneline/  (one question folder since 2026-09-19)
#                 coastsat_endpoint_vs_duneline_endpoint/<window>/
#                     (+ all_windows_stacked/)
#                 coastsat_total_change_vs_duneline_endpoint/<window>/
#                     (+ all_windows_stacked/, change_between_periods/)
#             duneline_positions/     where the 1997/2009/2023 dune lines sat:
#                                     maps, zooms, dune-road, beach width
#             (coastsat_windows/ and duneline_windows/ archived 2026-09-19 as
#             duplicates of 3-rates; trajectory_patterns/ and
#             two_period_comparison/ outputs deleted as stale)
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
# (scripts/input_prep/5-scr/3-rates/rates_figures.py); comparisons are in 4-comparisons.
COASTSAT_RATES = RATES / "coastsat"
DUNELINE_RATES = RATES / "duneline"
COASTSAT_LRR_ROOT = COASTSAT_RATES / "lrr"

# ---------------------------------------------------------------------------
# WHAT EACH WINDOW IS. Five window folders sit as peers under most products
# and they are NOT equivalent: two chains, and one context window that is not
# graded at all. Nothing in a folder name says so, which is how 1984_2004 gets
# read as a current result (Hannah, 2026-09-21). This dict is the ONE
# definition -- data/hatteras_init/5-scr/WINDOWS.md is written from it and
# coastsat_total_change.py imports it rather than keeping its own copy.
#
# All four periods are live in hatteras_site_config.HATTERAS_PERIODS; "legacy"
# here means superseded as the MAIN chain, not deleted. See
# [[cascade-canonical-periods]]: 1996 -> 2010 -> 2024 is the main chain.
CURRENT_CHAIN = (1996, 2010, 2024)
LEGACY_CHAIN = (1984, 2004, 2024)
WINDOW_ROLE = {
    (1996, 2010): "Calibration period",
    (2010, 2024): "Test period (held out)",
    (1996, 2024): "Full period (context, not graded)",
    (1984, 2004): "Calibration period, legacy chain",
    (2004, 2024): "Test period, legacy chain",
}
# The one-line gloss under each role, for READMEs and captions.
WINDOW_NOTE = {
    (1996, 2010): "the model is fitted here",
    (2010, 2024): "held out; nothing is fitted to it",
    (1996, 2024): "spans the whole current chain; CONTEXT ONLY, no run is "
                  "graded against it",
    (1984, 2004): "the 1984-start chain, superseded as the main chain 2026-09",
    (2004, 2024): "the 1984-start chain, superseded as the main chain 2026-09",
}


def window_role(window, with_note=False):
    """"Calibration period" for (1996, 2010); "" for a window with no role."""
    role = WINDOW_ROLE.get(tuple(window), "")
    if with_note and role:
        return f"{role} - {WINDOW_NOTE[tuple(window)]}"
    return role


def is_current_chain(window):
    """True for the 1996 -> 2010 -> 2024 chain and its full-period context."""
    return tuple(window) in {(1996, 2010), (2010, 2024), (1996, 2024)}

TRANSECT_DOMAINS = TRANSECT_FRAME / "transect_domains"
TIMESERIES_LRR = COASTSAT_RATES / "5yr_bins"
# Shoreline vs dune line, ONE question folder since 2026-09-19 (Hannah):
# endpoint_net_change/<window>/ (coastsat_vs_duneline.py, was 4-comparisons/
# coastsat_vs_duneline/), endpoint_net_change/chains/ (net_change_vs_duneline.py, was
# 4-comparisons/net_change_1996_2024/), total_change/<window>/
# (total_change_vs_duneline.py).
SHORELINE_VS_DUNELINE = COMPARISONS / "shoreline_vs_duneline"
# Named for BOTH sides since 2026-09-22 (Hannah): the folder alone says what
# was measured and how, with no lookup. The dune side is the same measured
# endpoint in every comparison here -- the 2009 digitized line minus the
# 1997 one -- so what the old names failed to say is how the SHORELINE was
# read, which is the only thing that differs.
#   was endpoint_net_change/  -> coastsat_endpoint_vs_duneline_endpoint/
#   was total_change/         -> coastsat_total_change_vs_duneline_endpoint/
COASTSAT_ENDPOINT_VS_DUNELINE = (SHORELINE_VS_DUNELINE
                                 / "coastsat_endpoint_vs_duneline_endpoint")
# The old name, kept so nothing that still imports it breaks silently.
DUNELINE_VS_COASTSAT = COASTSAT_ENDPOINT_VS_DUNELINE
# The stored dune-line observation (2026-09-18): the NET CHANGE between the
# two lines that bound a window, per transect and per domain, in m and m/yr,
# one folder per window. Written by
# scripts/input_prep/5-scr/3-rates/duneline/duneline_endpoint.py. It replaced
# duneline_lrr/ (an OLS through every line in the window, 2026-09-16; Hannah:
# "we are tracking net change"), now under archive/duneline_lrr_retired_20260918/.
DUNELINE_ENDPOINT_ROOT = DUNELINE_RATES / "endpoint"
# The CoastSat counterpart (2026-09-18): net change in the CoastSat shoreline
# between +/-6-month window means centred on the SAME dune-line survey dates,
# so the two products difference like for like. Written by
# scripts/input_prep/5-scr/3-rates/coastsat/endpoint/coastsat_endpoint.py.
COASTSAT_ENDPOINT_ROOT = COASTSAT_RATES / "endpoint"
# ---------------------------------------------------------------------------
# THE VOCABULARY (Hannah, by interview, 2026-09-21). An LRR turned into a
# DISTANCE is named by the window it was FITTED on, not by the arithmetic:
#
#   TOTAL SHORELINE CHANGE   the rate is evaluated over the SAME window it was
#                            fitted on. LRR(1996-2010) x 14 yr is total change,
#                            so is LRR(1996-2024) x 28 yr. A shorter span
#                            INSIDE the fit window (the 25.72 yr dune-line
#                            interval) is still total change; the span is named
#                            in the title, not in the folder.
#   PROJECTED SHORELINE      the rate is carried onto a window it was NOT
#   CHANGE                   fitted on. LRR(1996-2024) x 14 yr over 1996-2010
#                            or 2010-2024 is the only case in the project.
#   OBSERVED CHANGE          no rate anywhere: the mean position over the end
#                            calendar year minus the mean over the start one.
#
# Before 2026-09-21 the folder `lrr_projected/` held all three windows built
# from their OWN rate -- i.e. no projections at all -- which is what the rename
# below fixes. Old paths: lrr_projected/ -> total_change/.
# ---------------------------------------------------------------------------
# TOTAL shoreline change (2026-09-19, Hannah's advisor; renamed 2026-09-21):
# per transect lrr_m_yr x (end - start) years, the rate fitted on that same
# window, beside the OBSERVED change between calendar-year mean positions at
# the two ends. Windows 1996_2024, 1996_2010, 2010_2024. Written by
# scripts/input_prep/5-scr/3-rates/coastsat/total_change/coastsat_total_change.py.
COASTSAT_TOTAL_CHANGE_ROOT = COASTSAT_RATES / "total_change"
# PROJECTED shoreline change (2026-09-21, Hannah by interview): the 1996-2024
# LRR carried over each 14-yr half, against the same observed change. Only
# 1996_2010 and 2010_2024 exist -- 1996_2024 would BE the total change above.
# Written by the same script, --product projected.
COASTSAT_PROJECTED_ROOT = COASTSAT_RATES / "projected"
PROJECTED_RATE_WINDOW = (1996, 2024)
# Both endpoint products use the same two file names.
ENDPOINT_TRANSECT_FILE = "transect_endpoint.csv"
ENDPOINT_DOMAIN_FILE = "domain_endpoint_summary.csv"
DUNE_ENDPOINT_TRANSECT_FILE = ENDPOINT_TRANSECT_FILE
DUNE_ENDPOINT_DOMAIN_FILE = ENDPOINT_DOMAIN_FILE
# The CoastSat and dune-line net change side by side, 1996-2024 and its
# halves (net_change_vs_duneline.py, 2026-09-18): the net-change chain figure.
NET_CHANGE_1996_2024 = COASTSAT_ENDPOINT_VS_DUNELINE / "all_windows_stacked"
# TOTAL shoreline change against the dune line's measured net change, per
# domain, one folder per window: 1996_2010, 2010_2024 and 1996_2024, each
# window's rate fitted on that same window, x the CALENDAR span (14, 14, 28
# yr). Plus chains/ and difference/. Written by
# scripts/input_prep/5-scr/4-comparisons/shoreline_vs_duneline/total_change_vs_duneline.py.
#
# ONE tree since 2026-09-21 (Hannah, by interview). It absorbed two folders
# that were the same quantity under two names:
#   lrr_net_change/      the halves, each on its own rate -- renamed, it was
#                        never a projection and its docstring said so.
#   projected/1996_2024/ the 1996-2024 rate x the 25.72 yr dune-line interval.
#                        Not a projection either: same fit window, shorter
#                        span. Its numbers were already carried here as the
#                        *_dune_interval_m columns (verified identical to
#                        0 m before the merge), so the headline figure is the
#                        28 yr calendar span and the dune interval is a column
#                        and a caption line. Retired to superseded_20260921/.
TOTAL_CHANGE_VS_DUNELINE = (SHORELINE_VS_DUNELINE
                            / "coastsat_total_change_vs_duneline_endpoint")
# The same comparison with the shoreline side PROJECTED instead: the
# 1996-2024 LRR carried onto each 14-yr half, against the dune line measured
# over that half (2026-09-22, Hannah). 1996_2010 and 2010_2024 only -- over
# the full period the rate window IS the change window, which is
# TOTAL_CHANGE_VS_DUNELINE/1996_2024. Written by the same script,
# total_change_vs_duneline.py --product projected.
PROJECTED_VS_DUNELINE_ENDPOINT = (SHORELINE_VS_DUNELINE
                                  / "coastsat_projected_vs_duneline_endpoint")
# Where the dune line sat in 1997, 2009 and 2023: maps, imagery zooms, and
# its distance to NC-12 and to the CoastSat shoreline (duneline_positions.py,
# 2026-09-18).
DUNELINE_POSITIONS = COMPARISONS / "duneline_positions"
SHORELINE_INVENTORY = OBSERVATIONS / "shoreline_inventory"
# Outputs deleted 2026-09-19 (stale, on the old 1984/2004 periods); the
# scripts would regenerate here.
SHORELINE_PATTERNS = COMPARISONS / "trajectory_patterns"
DSAS_ROOT = OBSERVATIONS / "dsas_1978_2019"
# The four windows drawn on one y axis (coastsat_lrr_windows.py). Since
# 2026-09-19 only its 2 x 2 is drawn, into 3-rates/coastsat/lrr/; the
# per-window and halves figures duplicated 3-rates and were archived with
# the old folder (archive/2026-09-19_4-comparisons_duplicates/).
COASTSAT_LRR_WINDOWS = COASTSAT_LRR_ROOT
# RETIRED 2026-09-19: duneline_windows.py duplicated 3-rates/duneline/
# endpoint (window + chain figures); its output is archived.
DUNELINE_WINDOWS = ARCHIVE / "2026-09-19_4-comparisons_duplicates" / "duneline_windows"
# Two-window comparison figures (coastsat_two_period_comparison.py). Outputs
# deleted 2026-09-19 (stale, old periods); the script would regenerate here.
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
            scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_domain_lrr.py.
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
            "scripts/input_prep/5-scr/3-rates/coastsat/endpoint/coastsat_endpoint.py".format(
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
            "scripts/input_prep/5-scr/3-rates/duneline/duneline_endpoint.py".format(
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
# Built by scripts/input_prep/5-scr/3-rates/coastsat/extension/coastsat_extension_lrr.py.
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
