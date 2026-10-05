"""
Where is the observed shoreline data: the CoastSat chainage, the rate fits and the transect lookup?

    from site_layer.hat_observed_rates import lrr_csv, transect_lookup

Resolves data/hatteras_init/5-scr/ once; the rate CSV is a model input, and a
window not on disk raises, naming the ones that are. Details: scripts/site_layer/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

from pathlib import Path

# Run as a file, scripts/ is not on sys.path: add it so the site_layer imports resolve
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
# 3-rates is grouped by source: coastsat/{lrr,endpoint,5yr_bins} and duneline/endpoint
COASTSAT_RATES = RATES / "coastsat"
DUNELINE_RATES = RATES / "duneline"
COASTSAT_LRR_ROOT = COASTSAT_RATES / "lrr"
# The lrr field smoothed alongshore at two LOWESS widths, per window (coastsat_lrr_smoothed.py)
COASTSAT_LRR_SMOOTHED_ROOT = COASTSAT_RATES / "lrr_smoothed"

# What each window is: two chains and one context window; the one definition (WINDOWS.md is written from it)
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
# Shoreline vs dune line: one question folder
SHORELINE_VS_DUNELINE = COMPARISONS / "shoreline_vs_duneline"
# Named for both sides, so the folder says how the shoreline was read
COASTSAT_ENDPOINT_VS_DUNELINE = (SHORELINE_VS_DUNELINE
                                 / "coastsat_endpoint_vs_duneline_endpoint")
# The old name, kept so nothing that still imports it breaks silently.
DUNELINE_VS_COASTSAT = COASTSAT_ENDPOINT_VS_DUNELINE
# The dune-line net change between the two lines bounding a window (duneline_endpoint.py)
DUNELINE_ENDPOINT_ROOT = DUNELINE_RATES / "endpoint"
# The CoastSat counterpart, on the same survey dates (coastsat_endpoint.py)
COASTSAT_ENDPOINT_ROOT = COASTSAT_RATES / "endpoint"
# Vocabulary: total change (rate x its own window), projected (onto another window), observed (no rate)

# Total shoreline change: each window's rate x its own span, beside the observed change
COASTSAT_TOTAL_CHANGE_ROOT = COASTSAT_RATES / "total_change"
# Projected shoreline change: the 1996-2024 rate carried over each 14-yr half
COASTSAT_PROJECTED_ROOT = COASTSAT_RATES / "projected"
PROJECTED_RATE_WINDOW = (1996, 2024)
# Window convergence: nested windows walking toward 1996-2024 from both sides
# Renamed from window_convergence/ 2026-10-02, when the 1996-2026 sibling arrived
COASTSAT_WINDOW_CONVERGENCE_ROOT = COASTSAT_RATES / "window_convergence_1996_2024"


# Question first: 1-rate_profiles, 2-r_bias_rmse, 3-settling_window, experiments
# (2-r_bias_rmse added and settling renumbered 2 -> 3, Hannah 2026-10-01)
WINDOW_PROFILES_DIR = "1-rate_profiles"
WINDOW_SCORES_DIR = "2-r_bias_rmse"
SETTLING_WINDOW_DIR = "3-settling_window"
# Per-transect positions with the long-term fit and the two halves (Hannah, 2026-10-02)
SPLIT_WINDOWS_DIR = "4-split_windows"
WINDOW_EXPERIMENTS_DIR = "experiments"

# The settling sweep's three scales, keyed as its --scale choices are
SETTLING_SCALE_DIRS = {
    "sites": "a-eight_sites",
    "all": "b-every_transect",
    "domains": "c-domain_means",
}


def window_convergence_root(ref_start=1996, ref_end=2024) -> Path:
    """The window-convergence tree for one reference record. 1996-2024 is
    `window_convergence_1996_2024/`; a record that runs PAST 2024 gets its own sibling
    tree with the same 1-/2-/3- layout, `window_convergence_<start>_<end>/`
    (1996-2026, Hannah 2026-10-02). A record cut short stays an experiment
    under the 1996-2024 tree (see window_convergence_dir)."""
    ref_start, ref_end = int(ref_start), int(ref_end)
    if (ref_start, ref_end) == (1996, 2024):
        return COASTSAT_WINDOW_CONVERGENCE_ROOT
    return COASTSAT_RATES / "window_convergence_{0}_{1}".format(ref_start, ref_end)


def _extends_record(ref_start, ref_end):
    return int(ref_start) == 1996 and int(ref_end) > 2024


def window_convergence_dir(direction, anchor_year,
                          ref_start=1996, ref_end=2024) -> Path:
    """One settling sweep's folder: `forward_from_1996` or `backward_from_2024`.

    NOT `<start>_<end>`, although rule 2 would ask for it: a span name claims
    ONE interval and each of these folders holds twenty-five of them. What
    they share is the end that is PINNED, so that is what the name gives.

    The full 1996-2024 record files under `3-settling_window/`. Any other
    record span is a DIFFERENT EXPERIMENT, not a version of the same product,
    because every window in it is fitted against a different reference, so it
    files under `experiments/record_cut_<end>/` (or `record_<start>_<end>/` if
    the start moved too). Filing them apart also stops a truncated forward run
    overwriting the full one, which shares its pinned year and so its folder
    name (Hannah, 2026-09-23).

    Two directions, because the pair brackets the answer (Hannah, 2026-09-23).
    Forward pins 1996 and walks the end year outward: how much record do you
    need from the start of the chain? Backward pins 2024 and walks the start
    year back: how late can a window begin and still recover the long-term
    rate? Both converge on the same 1996-2024 reference from opposite sides.
    """
    _check_direction(direction)
    ref_start, ref_end = int(ref_start), int(ref_end)
    if (ref_start, ref_end) == (1996, 2024) or _extends_record(ref_start, ref_end):
        base = window_convergence_root(ref_start, ref_end) / SETTLING_WINDOW_DIR
    else:
        record = ("record_cut_{0}".format(ref_end) if ref_start == 1996
                  else "record_{0}_{1}".format(ref_start, ref_end))
        base = COASTSAT_WINDOW_CONVERGENCE_ROOT / WINDOW_EXPERIMENTS_DIR / record
    return base / "{0}_from_{1}".format(direction, int(anchor_year))


def window_profiles_dir(direction, anchor_year, ref_start=1996, ref_end=2024) -> Path:
    """One rate-profile family's folder, under `1-rate_profiles/` of the
    record's tree; the profiles have no truncated experiment."""
    _check_direction(direction)
    return (window_convergence_root(ref_start, ref_end) / WINDOW_PROFILES_DIR
            / "{0}_from_{1}".format(direction, int(anchor_year)))


def window_scores_dir(direction=None, anchor_year=None,
                      ref_start=1996, ref_end=2024) -> Path:
    """The r / bias / RMSE scores of the rate profiles, under `2-r_bias_rmse/`
    of the record's tree: the folder itself (both directions together) or,
    given a direction, that direction's folder."""
    base = window_convergence_root(ref_start, ref_end) / WINDOW_SCORES_DIR
    if direction is None:
        return base
    _check_direction(direction)
    return base / "{0}_from_{1}".format(direction, int(anchor_year))


def split_windows_dir() -> Path:
    """The per-transect split-window figures and picks, under `4-split_windows/`."""
    return COASTSAT_WINDOW_CONVERGENCE_ROOT / SPLIT_WINDOWS_DIR


def _check_direction(direction):
    if direction not in ("forward", "backward"):
        raise ValueError(
            "direction is 'forward' (pinned start) or 'backward' (pinned end), "
            "not {0!r}".format(direction))


WINDOW_CONVERGENCE_SWEEP_FILE = "window_convergence_transects.csv"
WINDOW_CONVERGENCE_SUMMARY_FILE = "convergence_summary.csv"
# Both endpoint products use the same two file names.
ENDPOINT_TRANSECT_FILE = "transect_endpoint.csv"
ENDPOINT_DOMAIN_FILE = "domain_endpoint_summary.csv"
DUNE_ENDPOINT_TRANSECT_FILE = ENDPOINT_TRANSECT_FILE
DUNE_ENDPOINT_DOMAIN_FILE = ENDPOINT_DOMAIN_FILE
# The CoastSat and dune-line net change side by side, 1996-2024 and its halves
NET_CHANGE_1996_2024 = COASTSAT_ENDPOINT_VS_DUNELINE / "all_windows_stacked"
# Total shoreline change against the dune line's net change, per domain, per window
TOTAL_CHANGE_VS_DUNELINE = (SHORELINE_VS_DUNELINE
                            / "coastsat_total_change_vs_duneline_endpoint")
# The same comparison with the shoreline side projected (1996_2010, 2010_2024 only)
PROJECTED_VS_DUNELINE_ENDPOINT = (SHORELINE_VS_DUNELINE
                                  / "coastsat_projected_vs_duneline_endpoint")
# Where the dune line sat in 1997, 2009 and 2023 (duneline_positions.py)
DUNELINE_POSITIONS = COMPARISONS / "duneline_positions"
# One period's mean shoreline over two windows, differenced (coastsat_mean_shoreline_windows.py)
MEAN_SHORELINE_WINDOWS = COMPARISONS / "mean_shoreline_windows"
SHORELINE_INVENTORY = OBSERVATIONS / "shoreline_inventory"

# The detrended position: the signal the whole island shares (an observation, not a rate)
DETRENDED_POSITION = OBSERVATIONS / "detrended_position"
DETRENDED_POSITION_MATRIX = "annual_medians_detrended.csv"
# Outputs deleted as stale; the scripts would regenerate here
SHORELINE_PATTERNS = COMPARISONS / "trajectory_patterns"
DSAS_ROOT = OBSERVATIONS / "dsas_1978_2019"

# The mean shoreline as a line on the ground: a position, so the geolocation is real work
MEAN_SHORELINE_ROOT = OBSERVATIONS / "mean_shoreline"


def mean_shoreline_label(start, end) -> str:
    """The token a window's folder and files are named with.

    Calendar years give `1995_1997`. Since 2026-09-29 a window can also be
    given by ISO dates -- the +/-1 yr windows centred on a start DEM's flights
    -- and then the token is the dates, `1995-10-12_1997-10-12` (Hannah: name
    them by exact dates, beside the calendar folders)."""
    def token(v):
        s = str(v)
        return s if "-" in s else str(int(s))
    return "{0}_{1}".format(token(start), token(end))


def mean_shoreline_dir(start, end) -> Path:
    """One averaging window's folder: the line, the per-transect means, the
    figure and PROVENANCE.md. Named for the span, per rule 2. `start`/`end`
    are calendar years or ISO dates (see mean_shoreline_label)."""
    return MEAN_SHORELINE_ROOT / mean_shoreline_label(start, end)


def mean_shoreline_geojson(start, end) -> Path:
    """The mean shoreline as ONE LineString, EPSG:26918, carrying the same
    metadata properties a digitised dune line carries so the 2-brie-offset
    intersection step reads it unchanged."""
    label = mean_shoreline_label(start, end)
    return mean_shoreline_dir(start, end) / "shoreline_mean_{0}.geojson".format(label)


def mean_shoreline_csv(start, end) -> Path:
    """Per CoastSat transect: the window mean, its scatter and count, the
    dates it spans, and the geolocated mean point."""
    label = mean_shoreline_label(start, end)
    return mean_shoreline_dir(start, end) / "transect_means_{0}.csv".format(label)
# The four windows on one y axis (coastsat_lrr_windows.py)
COASTSAT_LRR_WINDOWS = COASTSAT_LRR_ROOT
# Two-window comparison figures; outputs deleted as stale
TWO_PERIOD_COMPARISON = COMPARISONS / "two_period_comparison"
# Retired windows, kept for the DSAS comparison in 6-scr-smooth
COASTSAT_LRR_SUPERSEDED = ARCHIVE / "coastsat_lrr_superseded_20260810"

# 6-scr-smooth: what the LOWESS smoothing does to the observed rates
SMOOTH_ROOT = INIT_ROOT / "6-scr-smooth"
SMOOTH_METHOD_COMPARISON = SMOOTH_ROOT / "method_comparison"
SMOOTH_DSAS_VS_COASTSAT = SMOOTH_ROOT / "dsas_vs_coastsat"

# Single files in transect_domains/ read by name outside 5-scr (HAT_domains.json = the 90 real boxes)
DOMAIN_BOXES = TRANSECT_DOMAINS / "HAT_domains.json"
TRANSECT_LAYER = TRANSECT_DOMAINS / "CoastSat_transect_layer.geojson"

# The per-window products a rate fit writes.
TRANSECT_FILE = "transect_lrr_full.csv"
DOMAIN_FILE = "domain_lrr_summary.csv"

# Folders under coastsat_lrr/ that are not windows
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


# The extension: its own lookup and rate fit per window, never merged into the surveyed ones
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
