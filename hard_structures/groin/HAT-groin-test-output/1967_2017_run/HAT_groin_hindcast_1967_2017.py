"""
HAT_groin_threeway_hindcast_1967_2017.py
================================
Extended groin-DETERIORATION-test hindcast (1967-2017, GIS D2-D12). Built
from HAT_groin_hindcast_1967_1997.py -- same inputs/conventions -- extended
from a 30-year (1967-1997) groin-only test to a 51-year (1967-2017) test
that also exercises the 1995->2003 deterioration ramp, which needs the
longer window to actually play out (both dates fall after 1997, so the
original 30-year test never reached them).

  RUN_MATRIX = ["no_groin", "groin"]
    -> HAT_1967_2018_M60_deterioration_no_groin   (groin OFF -- erosive baseline)
    -> HAT_1967_2018_M60_deterioration_groin      (groin ON, WITH deterioration)

The "no_groin" run needs nothing beyond your working base setup. The "groin"
run additionally requires:
  1. HAT_groin_module.py importable (same folder / on PYTHONPATH), and
  2. the inert pre-AST hook in cascade.py (sets cascade._groin_callback).
If the hook is missing, the groin run warns loudly (diagnostics stay empty)
so you never mistake a no-op for a real groin run.

IMPORTANT -- END_YEAR CONVENTION: RUN_YEARS = END_YEAR - START_YEAR is
EXCLUSIVE of END_YEAR (the original 1967-1997 run, with END_YEAR=1997, only
ever simulated through 1996). To actually reach 2017 as the final modeled
year, END_YEAR is set to 2018 below (51 model years, 1967-2017 inclusive) --
matching the storm file built by HAT_build_1967_2017_storms.py.

DETERIORATION (Section 2, GROIN_DETERIORATION_*): last repair 1995 (after
Hurricane Gordon damage in 1994) -> linear decline -> Hurricane Isabel 2003
locks in the new deteriorated state. Floor fraction M/3 (~0.333), from
Katherine's "one groin still functional out of three" framing -- see chat
history with Laura/Katherine on the Coastal Sediments 2027 abstract for the
full reasoning behind these specific numbers.

Historical beach nourishment (Section 2c: 1971, 1973) is applied identically
to EVERY run in RUN_MATRIX -- same treatment as the edge BE correction
(Section 2b) -- since these are documented real-world events, not an
experimental variable. See Section 2c for full sourcing.

Plot with HAT_groin_effect_comparison.py:
    RUN_NO_GROIN = "HAT_1967_2018_M60_deterioration_no_groin"
    RUN_GROIN    = "HAT_1967_2018_M60_deterioration_groin"
(no_groin first = baseline for the difference figure).
"""

import os
import sys
import numpy as np
import pandas as pd

# ── Which Cascade to use ──────────────────────────────────────────────────────
# While TESTING the groin, use a SANDBOX copy of cascade.py (named cascade_groin.py)
# that lives INSIDE the cascade package folder, next to the real cascade.py, with
# the 3-line groin hook added. Your real cascade.py stays untouched. Once the groin
# is proven, fold the hook into the real cascade.py and set USE_SANDBOX_CASCADE=False.
USE_SANDBOX_CASCADE = True

if USE_SANDBOX_CASCADE:
    from cascade.cascade_groin import Cascade   # hooked sandbox copy
else:
    from cascade import Cascade                 # real package (hook folded in)

# =============================================================================
# SECTION 1: DOMAIN CONFIGURATION  (D2-D12)
# =============================================================================
NUM_REAL_DOMAINS   = 11
NUM_BUFFER_DOMAINS = 15
FIRST_FILE_NUMBER  = 2
LAST_FILE_NUMBER   = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1     # 12
TOTAL_DOMAINS      = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS  # 41
START_REAL_INDEX   = NUM_BUFFER_DOMAINS                            # 15
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS           # 26


def _gis_to_pad(gis_id):
    """D2->15, D5->18, D6->19, D12->25."""
    return START_REAL_INDEX + (gis_id - FIRST_FILE_NUMBER)


# =============================================================================
# SECTION 2: WHICH RUNS + GROIN CONFIG
# =============================================================================
# Pick which run(s) to do. Put ONE in the list to run a single experiment at a
# time; list several to run them in sequence. Names match the comparison folders and
# the plotter's RUNS list.
#   "no_groin"  -> HAT_1967_2018_M60_deterioration_no_groin   (erosive baseline)
#   "groin"     -> HAT_1967_2018_M60_deterioration_groin      (dipole at D5/D6, WITH
#                  deterioration; needs hook)
#   "groin_be"  -> HAT_1967_2018_M60_deterioration_groin_be   (groin + background erosion)
RUN_MATRIX = ["groin"]              # <- edit this line to choose the run(s)

# Buxton groin: sits in D6 (source/accretion), starves D5 (sink/erosion).
GROIN_UPDRIFT_GIS   = 6
GROIN_DOWNDRIFT_GIS = 5
GROIN_TRAPPING_RATE_M_YR = 60.0    # M -- the single knob; tune to observed updrift
GROIN_INSTALL_YEAR  = 1970         # inert before this (free 1967-69 control window)

# Deterioration: last repair 1995 (after Hurricane Gordon, 1994) -> linear
# decline -> Hurricane Isabel 2003 locks in the new deteriorated state.
# Delay is relative to GROIN_INSTALL_YEAR (this script's own 1970, not the
# true historical 1969), so the module resolves the right calendar years
# regardless of which install-year convention this test uses.
GROIN_DETERIORATION_DELAY_YEARS = 1995 - GROIN_INSTALL_YEAR   # = 25
GROIN_DETERIORATION_MODE        = "linear_ramp"
GROIN_DETERIORATION_RAMP_YEARS  = 2003 - 1995                  # = 8
# Floor fraction: originally M/3 (~0.333, Katherine's "one groin still
# functional out of three" framing), but that removed too much of the groin
# signal -- testing 0.5 (50% trapping retained) instead as a less severe
# deterioration assumption. Still folded into the ramp itself (not run as a
# separate standalone experiment) per the plan worked out with
# Laura/Katherine for the Coastal Sediments 2027 abstract.
GROIN_DETERIORATION_FRACTION    = 0.50

# Optional regional background erosion for a "groin_be" run (only used if that
# key is in RUN_MATRIX). m/yr, negative = erosive.
REGIONAL_BE_RATE_M_YR = 0

# =============================================================================
# SECTION 2b: EDGE SOURCE/SINK CORRECTION (buffer-orientation boundary fix)
# =============================================================================
# Mirrors the main 1984-2024 hindcast's edge correction at GIS 1 / GIS 90: the
# outermost REAL domains (here D2 and D12) sit directly against buffer padding,
# and the buffer's flat/repeated orientation can introduce an artificial
# alongshore-transport signal right at that boundary. Set a background_erosion
# value (m/yr, same sign convention as REGIONAL_BE_RATE_M_YR: negative =
# erosive) at the edge domain(s) to correct for it.
#
# This is a STRUCTURAL fix (same role as GIS 1 / 90 in the main script), not a
# scientific choice like REGIONAL_BE_RATE_M_YR below -- so it's applied to
# EVERY run in RUN_MATRIX (no_groin, groin, groin_be alike), and it STACKS
# additively with REGIONAL_BE_RATE_M_YR on groin_be runs rather than being
# overwritten by it. Set both edge values to 0.0 to disable.
APPLY_EDGE_BE_CORRECTION = True
EDGE_BE_RATES_GIS = {
    2:  5.0,   # D2  -- south edge, against buffer (solve for this)
    12: 10.0,   # D12 -- north edge, against buffer (solve for this)
}

# =============================================================================
# SECTION 2c: HISTORICAL BEACH NOURISHMENT (1971, 1973)
# =============================================================================
# Two documented NPS nourishment projects fall within this 1967-2017 window,
# both targeting the erosion embayment adjacent to the (Navy, 1969) groin
# field this script already models. Applied to EVERY run in RUN_MATRIX
# (no_groin, groin, groin_be alike) -- same treatment as the edge BE
# correction above -- since these are real historical events, not a
# scientific/experimental variable.
#
# SOURCES:
#   - Machemehl (1973): reports the 1971 volume as 300,000 cy, sourced from a
#     man-made lake at Cape Point, pumped via 14-in cutterhead dredge
#     (JA LaPort Dredging Co.) ~3.5 mi to a discharge point near the Hatteras
#     Court Motel, left to migrate south under normal littoral drift.
#   - NPS (1980), pg 48: reports the 1971 volume as 200,000 cy (borrow
#     material "proved insufficient to have any significant impact"); reports
#     the 1973 volume as 1,300,000 cy from an interior Cape Point borrow area
#     (basin still visible today as altered vegetation), discharged via a
#     16-in dredge + 3 boosters ~4 mi north near the Hatteras Court Motel,
#     widening the beach ~500 ft over a cited 5,000-ft reach.
#   - Dolan, Hayden, Riddel & Ponton (1974), "1973 Buxton Beach Nourishment
#     Project: An Annotated Photographic Atlas," NPS Contract No.
#     CX5000031059 (the primary, station-surveyed source behind NPS 1980's
#     1973 figures): independently confirms 200,000 cy for 1971 (agreeing
#     with NPS 1980, not Machemehl); gives exact engineering stations for the
#     1973 project (south limit STA 2235+00, explicitly excluding the groin
#     cells -- "no material was pumped into the groin cells, southern
#     sediment drift caused a build-up" -- north limit near STA 2164+80,
#     MP41.0) and confirms the Navy's 1969 groin field is this script's
#     modeled groin.
#
# VOLUME CHOICE: using NPS(1980)'s 200,000 cy for 1971 rather than
# Machemehl's 300,000 cy -- two independent sources (NPS 1980, Dolan et al.
# 1974) agree on 200,000 cy against one source for 300,000 cy.
#
# DOMAIN RANGE (derived, NOT yet confirmed against the project's own GIS
# domain shapefile -- verify before treating these as final):
#   1973: templated engineered fill, ~D6-D10. Anchored on the atlas's own
#     station data (south limit STA 2235 ~= 0.24 domains north of the groin;
#     north limit STA 2164+80 ~= 4.5 domains north of the groin), converted
#     at 500 m/domain and 100 ft/station, using this script's GROIN_UPDRIFT/
#     DOWNDRIFT_GIS (D5/D6) as the shared anchor point with the atlas's
#     "Navy Groins" landmark.
#   1971: NOT a templated fill -- the atlas gives no station range for it,
#     only that it targeted the same "north embayment" and was discharged as
#     a point source near the motels, then left to migrate south under
#     natural littoral drift. Modeled here as a SINGLE-domain point
#     injection at D8 (the motels, per the derived station crosswalk),
#     deliberately NOT spread across a template -- letting BRIE's own
#     alongshore transport carry it south is the same "extent is emergent,
#     not prescribed" logic already used for the groin dipole, and matches
#     what the historical record says actually happened.
ENABLE_HISTORICAL_NOURISHMENT = True

_CY_TO_M3       = 0.764555   # cubic yards -> cubic metres
DOMAIN_LENGTH_M = 500        # alongshore width of one CASCADE domain (m)

# REMINDER: every volume below is stored/printed in m^3, NOT cubic yards.
# 1 cy = 0.764555 m^3, so the m^3 figures look ~24% SMALLER than the cy
# numbers you'll find quoted in the source reports (e.g. 1973's "1,300,000
# cy" becomes 993,921.5 m^3 total). If a number here looks low compared to
# what you remember reading in a report, check units before assuming an
# error -- it's very likely just cy vs m^3, not a mistake.

# Calendar years of nourishment events - SORTED, one entry per distinct year
HAT_BN_YEARS = [1971, 1973]
#                 ^       ^
#              point-source  templated fill
#              near D8        D6-D10

# HAT_BN_VOLUME_BY_DOMAIN
# Key   : GIS domain ID (2-12)
# Value : [volume_for_1971_m3, volume_for_1973_m3]
#         0 = not nourished that year; non-zero = total m^3 for that domain
HAT_BN_VOLUME_BY_DOMAIN = {
    # --- 1971: single-domain point injection at D8 (Hatteras Court Motel) ---
    #   200,000 cy x 0.764555 = 152,911.0 m^3, placed entirely in D8
     8: [round(200_000 * _CY_TO_M3, 1), 0],
    # --- 1973: templated fill, D6-D10 (5 domains, 1,300,000 cy) ---
    #   1,300,000 / 5 x 0.764555 = 198,784.3 m^3/domain
     6: [0, round(1_300_000 / 5 * _CY_TO_M3, 1)],
     7: [0, round(1_300_000 / 5 * _CY_TO_M3, 1)],
     9: [0, round(1_300_000 / 5 * _CY_TO_M3, 1)],
    10: [0, round(1_300_000 / 5 * _CY_TO_M3, 1)],
}
# NOTE: D8 needs entries for BOTH years (1971 point injection AND its share of
# the 1973 template), so its list is built by combining the two rather than
# appearing twice as a dict key.
HAT_BN_VOLUME_BY_DOMAIN[8][1] = round(1_300_000 / 5 * _CY_TO_M3, 1)

if not ENABLE_HISTORICAL_NOURISHMENT:
    HAT_BN_YEARS = []
    HAT_BN_VOLUME_BY_DOMAIN = {}

# Per-domain flag passed to Cascade's beach_nourishment_module: True only for
# domains that receive a nourishment event at some point in this window.
NOURISHMENT_MANAGEMENT_ON = [False] * TOTAL_DOMAINS
for _gis_id in HAT_BN_VOLUME_BY_DOMAIN:
    _pad_idx = _gis_to_pad(_gis_id)
    if 0 <= _pad_idx < TOTAL_DOMAINS:
        NOURISHMENT_MANAGEMENT_ON[_pad_idx] = True

# =============================================================================
# SECTION 3: PERIOD / FILE PATHS   (identical to base run)
# =============================================================================
PROJECT_BASE_DIR   = r"/"
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, "data", "hatteras_init")
OUTPUT_BASE_DIR    = os.path.join(PROJECT_BASE_DIR, "output", "raw_runs")
PARAMETER_FILE     = "Hatteras-CASCADE-parameters.yaml"

START_YEAR = 1967
END_YEAR   = 2018   # EXCLUSIVE convention (RUN_YEARS = END_YEAR - START_YEAR)
                     # -> 51 model years, final simulated year = 2017
RUN_YEARS  = END_YEAR - START_YEAR          # 51

SEA_LEVEL_RISE_RATE = 0.004
SEA_LEVEL_CONSTANT  = True

GROIN_INIT_DIR = os.path.join(
    PROJECT_BASE_DIR, "scripts", "groin", "HAT-hindcast-groin-test", "groin_init",
)
STORM_FILE = os.path.join(
    GROIN_INIT_DIR, "storms", "1967_2017",
    "1967_2017_groin_storms.npy",
)
ISLAND_OFFSET_FILE = os.path.join(
    GROIN_INIT_DIR, "2-brie-offset",
    "Island_Dune_Offsets_1967_D2_D12_PADDED_41.csv",
)

TOPO_DUNE_INIT_YEAR = "2009"
TOPO_DUNE_SUBFOLDER = "2009"

# =============================================================================
# SECTION 4: SIMULATION PARAMETERS   (identical to base run)
# =============================================================================
BERM_ELEVATION = 1.7
MHW_ELEVATION  = 0.36
NUM_CORES      = 1

WAVE_HEIGHT_M                  = 1.0
FIXED_WAVE_PERIOD              = 8
FIXED_WAVE_ASYMMETRY           = 0.7
FIXED_WAVE_ANGLE_HIGH_FRACTION = 0.1

DUNE_REBUILD_HEIGHT     = 3.0
REBUILD_ELEV_THRESHOLD  = 0.01   # dam
OVERWASH_TO_DUNE        = 0.0
OVERWASH_FILTER_DEFAULT = 0.0

RUN_NAME_SUFFIX = "edge_calibrated"    # -> HAT_1967_2018_M60_deterioration_{run_key}

# Auto-generate + save figures into each run's folder when the run finishes.
MAKE_FIGURES = True
MAKE_RUN_GIF = True     # animate the modeled shoreline evolving over the run

os.chdir(PROJECT_BASE_DIR)
os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)


# =============================================================================
# SECTION 5: LOADERS   (verbatim from base run)
# =============================================================================
def check_inputs_exist():
    for label, path in [("STORM_FILE", STORM_FILE),
                        ("ISLAND_OFFSET_FILE", ISLAND_OFFSET_FILE),
                        ("PARAMETER_FILE", os.path.join(HATTERAS_DATA_BASE, PARAMETER_FILE))]:
        if not os.path.isfile(path):
            print(f"CRITICAL ERROR: Missing data file ({label})")
            print(f"  as given: {path}")
            print(f"  abspath:  {os.path.abspath(path)}")
            sys.exit(1)
    print("  All required input files found.")


def load_island_offset_dam():
    offset_all = np.loadtxt(ISLAND_OFFSET_FILE, skiprows=1, delimiter=",")
    offset_dam = offset_all / 10.0
    if offset_dam.size != TOTAL_DOMAINS:
        sys.exit(f"ERROR: offset has {offset_dam.size} values, expected {TOTAL_DOMAINS}.")
    print(f"  Loaded island offsets: {offset_dam.size} domains (dam)")
    return list(offset_dam)


def build_file_lists():
    elev, dune = [], []
    for _ in range(START_REAL_INDEX):
        dune.append(os.path.join(HATTERAS_DATA_BASE, "buffer", "sample_1_dune.npy"))
        elev.append(os.path.join(HATTERAS_DATA_BASE, "buffer", "sample_1_topography.npy"))
    for i_list in range(START_REAL_INDEX, END_REAL_INDEX):
        file_num = FIRST_FILE_NUMBER + (i_list - START_REAL_INDEX)
        dune.append(os.path.join(HATTERAS_DATA_BASE, "dunes", TOPO_DUNE_SUBFOLDER,
                                 f"domain_{file_num}_dune_{TOPO_DUNE_INIT_YEAR}.npy"))
        elev.append(os.path.join(HATTERAS_DATA_BASE, "topography", TOPO_DUNE_SUBFOLDER,
                                 f"domain_{file_num}_topography_{TOPO_DUNE_INIT_YEAR}.npy"))
    for _ in range(END_REAL_INDEX, TOTAL_DOMAINS):
        dune.append(os.path.join(HATTERAS_DATA_BASE, "buffer", "sample_1_dune.npy"))
        elev.append(os.path.join(HATTERAS_DATA_BASE, "buffer", "sample_1_topography.npy"))
    print(f"  Generated {len(elev)} elevation + {len(dune)} dune file paths")

    missing = [p for p in set(elev + dune) if not os.path.isfile(p)]
    if missing:
        print("CRITICAL ERROR: missing init file(s):")
        for p in missing[:10]:
            print("  ", p)
        if len(missing) > 10:
            print(f"   ... and {len(missing) - 10} more")
        sys.exit(1)
    return elev, dune


def get_x_s_TS(b3d):
    if hasattr(b3d, "x_s_TS"):
        return np.asarray(b3d.x_s_TS, dtype=float)
    if hasattr(b3d, "_x_s_TS"):
        return np.asarray(b3d._x_s_TS, dtype=float)
    raise AttributeError("No x_s_TS / _x_s_TS on Barrier3D object.")


def build_shoreline_matrix(cascade, to_meters=True):
    b3d_list = cascade.barrier3d
    ndom = len(b3d_list)
    nt   = len(get_x_s_TS(b3d_list[0]))
    shoreline = np.zeros((nt, ndom), dtype=float)
    for j in range(ndom):
        shoreline[:, j] = get_x_s_TS(b3d_list[j])
    if to_meters:
        shoreline *= 10.0
    return shoreline


def build_nourishment_arrays_from_manual_inputs():
    """
    Build per-year nourishment-on and volume arrays for the CASCADE time loop,
    from HAT_BN_YEARS + HAT_BN_VOLUME_BY_DOMAIN (Section 2c). Mirrors the main
    1984-2024 hindcast's build_nourishment_arrays_from_manual_inputs() exactly.
    Years outside [START_YEAR, END_YEAR] are silently skipped, so
    ENABLE_HISTORICAL_NOURISHMENT=False returns all-zero arrays with no other
    code change needed.

    Returns
    -------
    nourishment_on_by_year     : dict {year: np.ndarray[TOTAL_DOMAINS]}, 1/0
    nourishment_volume_by_year : dict {year: list[TOTAL_DOMAINS]}, m^3/m
    """
    nourishment_on_by_year     = {}
    nourishment_volume_by_year = {}

    for year in range(START_YEAR, END_YEAR + 1):
        nourishment_on_by_year[year]     = np.zeros(TOTAL_DOMAINS)
        nourishment_volume_by_year[year] = [0.0] * TOTAL_DOMAINS

    for gis_id, volumes_m3 in HAT_BN_VOLUME_BY_DOMAIN.items():
        if len(HAT_BN_YEARS) != len(volumes_m3):
            raise ValueError(
                f"GIS domain {gis_id}: HAT_BN_YEARS and volume list must have "
                f"the same length ({len(HAT_BN_YEARS)} vs {len(volumes_m3)})."
            )

        pad_idx = _gis_to_pad(gis_id)
        if not (0 <= pad_idx < TOTAL_DOMAINS):
            print(f"  WARNING: GIS {gis_id} -> pad {pad_idx} out of range - skipped.")
            continue

        for year, total_m3 in zip(HAT_BN_YEARS, volumes_m3):
            if year < START_YEAR or year > END_YEAR:
                continue   # event outside this period - skip silently

            volume_m3_per_m = float(total_m3) / DOMAIN_LENGTH_M
            nourishment_on_by_year[year][pad_idx]     = 1
            nourishment_volume_by_year[year][pad_idx] = volume_m3_per_m

    # Print schedule summary
    has_events = False
    for year in range(START_YEAR, END_YEAR + 1):
        active_pad = np.where(nourishment_on_by_year[year] == 1)[0]
        if len(active_pad) > 0:
            has_events = True
            active_gis = [
                FIRST_FILE_NUMBER + (idx - START_REAL_INDEX)
                for idx in active_pad
                if START_REAL_INDEX <= idx < END_REAL_INDEX
            ]
            total_vol = (
                np.sum(np.asarray(nourishment_volume_by_year[year], dtype=float))
                * DOMAIN_LENGTH_M
            )
            print(f"  {year}: GIS domains {active_gis}  |  total = {total_vol:,.0f} m^3")

    if not has_events:
        print("  (no nourishment events in this period's date range)")

    return nourishment_on_by_year, nourishment_volume_by_year


# =============================================================================
# SECTION 6: RUN ONE  (base-run body + optional groin attachment)
# =============================================================================
def _save_run_figures(run_name, run_dir, shoreline_m):
    """Generate + save this run's figures into its own folder, reusing the
    plotting functions from HAT_plot_groin_runs.py (single source of truth --
    no duplicated plot code). Single-run figures: position change, rate,
    trajectories. Never blocks the run: any plotting error is caught + reported."""
    try:
        import matplotlib
        matplotlib.use("Agg")   # headless save, no popups from the run script
        import HAT_plot_groin_runs as P
    except Exception as e:
        print(f"  [figures skipped] could not import plotter: {e}")
        return

    runs_data = {run_name: shoreline_m}
    fig_makers = [P.fig_position_change, P.fig_change_rate, P.fig_trajectories,
                  P.fig_model_vs_observed, P.fig_position_planform]
    for maker in fig_makers:
        try:
            fig, tag = maker(runs_data)
            out = os.path.join(run_dir, f"{run_name}_PLOT_{tag}.png")
            fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
            import matplotlib.pyplot as plt
            plt.close(fig)
            print(f"  Saved figure: {os.path.basename(out)}")
        except Exception as e:
            print(f"  [figure '{maker.__name__}' skipped] {e}")


def _save_run_gif(run_name, run_dir, shoreline_m):
    """Animate the modeled shoreline over the run (real domains D2-D12), year by
    year, in POSITION mode matching the main hindcast: real planform relative to
    the year-0 alongshore mean, ocean at bottom (seaward downward). Shows the
    fillet growing at the groin against the real island orientation."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation, PillowWriter
    except Exception as e:
        print(f"  [run GIF skipped] animation unavailable: {e}")
        return

    OCEAN_AT_BOTTOM = True
    nt = shoreline_m.shape[0]
    gis = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)
    flip = -1.0   # x_s increases landward -> flip so + = seaward
    pos = flip * shoreline_m[:, START_REAL_INDEX:END_REAL_INDEX]
    ref_mean = np.nanmean(pos[0])                 # year-0 alongshore mean = the 0
    series = pos - ref_mean                        # real planform, position mode
    ymin, ymax = series.min(), series.max()
    pad = 0.1 * (ymax - ymin if ymax > ymin else 1)

    fig, ax = plt.subplots(figsize=(10, 5))

    def draw(t):
        ax.clear()
        # year-0 reference planform (dashed grey)
        ax.plot(gis, series[0], color="0.6", ls="--", lw=1.4, zorder=2)
        # current year
        ax.plot(gis, series[t], marker="o", ms=5, lw=2.2, color="#FF8C00", zorder=4)
        # shade seaward/landward of reference
        ax.fill_between(gis, series[0], series[t],
                        where=(series[t] >= series[0]), color="#4a90d9",
                        alpha=0.25, zorder=1)
        ax.fill_between(gis, series[0], series[t],
                        where=(series[t] < series[0]), color="#d95f5f",
                        alpha=0.25, zorder=1)
        ax.axvline(5.5, color="#B71C1C", ls="--", lw=1.5, alpha=0.9, zorder=3)
        if OCEAN_AT_BOTTOM:
            ax.set_ylim(ymax + pad, ymin - pad)   # seaward downward
        else:
            ax.set_ylim(ymin - pad, ymax + pad)
        ax.set_xlabel(f"GIS Domain ID (D{FIRST_FILE_NUMBER}-D{LAST_FILE_NUMBER})")
        up_word = "landward" if OCEAN_AT_BOTTOM else "seaward"
        ax.set_ylabel(f"Cross-shore position (m, rel. {START_YEAR} mean)  {up_word} \u25b2")
        ax.set_title(f"{run_name}  \u2014  {START_YEAR + t}")
        ax.text(0.02, 0.06, str(START_YEAR + t), transform=ax.transAxes,
                fontsize=22, fontweight="bold", color="#FF8C00", alpha=0.8)
        ax.grid(alpha=0.3)

    try:
        anim = FuncAnimation(fig, draw, frames=nt, interval=250)
        out = os.path.join(run_dir, f"{run_name}_PLOT_shoreline_evolution.gif")
        anim.save(out, writer=PillowWriter(fps=4))
        plt.close(fig)
        print(f"  Saved GIF: {os.path.basename(out)}")
    except Exception as e:
        plt.close(fig)
        print(f"  [run GIF skipped] {e}")


def run_one(run_key, island_offset_dam, elevation_files, dune_files,
            historical_nourishment_on_by_year, historical_nourishment_volume_by_year):
    groin_on = run_key in ("groin", "groin_be")
    be_on    = run_key == "groin_be"

    run_name = f"HAT_{START_YEAR}_{END_YEAR}_{RUN_NAME_SUFFIX}_{run_key}"
    print("\n" + "=" * 78)
    print(f"RUN: {run_name}   (groin={'ON' if groin_on else 'off'}, "
          f"BE={'ON' if be_on else 'off'}, "
          f"edge_correction={'ON' if APPLY_EDGE_BE_CORRECTION else 'off'})")
    if APPLY_EDGE_BE_CORRECTION:
        print(f"  Edge BE correction: "
              + ", ".join(f"D{g}={r:+.1f} m/yr" for g, r in EDGE_BE_RATES_GIS.items()))
    print("=" * 78)

    be = [0.0] * TOTAL_DOMAINS
    if APPLY_EDGE_BE_CORRECTION:
        for gis, rate in EDGE_BE_RATES_GIS.items():
            be[_gis_to_pad(gis)] += rate
    if be_on:
        for gis in range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1):
            be[_gis_to_pad(gis)] += REGIONAL_BE_RATE_M_YR

    cascade = Cascade(
        HATTERAS_DATA_BASE,
        run_name,
        storm_file=STORM_FILE,
        elevation_file=elevation_files,
        dune_file=dune_files,
        parameter_file=PARAMETER_FILE,

        berm_elevation=BERM_ELEVATION,
        MHW=MHW_ELEVATION,

        wave_height=WAVE_HEIGHT_M,
        wave_period=FIXED_WAVE_PERIOD,
        wave_asymmetry=FIXED_WAVE_ASYMMETRY,
        wave_angle_high_fraction=FIXED_WAVE_ANGLE_HIGH_FRACTION,

        sea_level_rise_rate=SEA_LEVEL_RISE_RATE,
        sea_level_rise_constant=SEA_LEVEL_CONSTANT,

        background_erosion=be,
        alongshore_section_count=TOTAL_DOMAINS,
        time_step_count=RUN_YEARS,

        min_dune_growth_rate=[0.55] * TOTAL_DOMAINS,
        max_dune_growth_rate=[0.95] * TOTAL_DOMAINS,
        num_cores=NUM_CORES,

        roadway_management_module=[False] * TOTAL_DOMAINS,
        beach_nourishment_module=NOURISHMENT_MANAGEMENT_ON,
        sandbag_management_on=[False] * TOTAL_DOMAINS,
        alongshore_transport_module=True,
        community_economics_module=False,

        dune_design_elevation=[DUNE_REBUILD_HEIGHT] * TOTAL_DOMAINS,
        dune_minimum_elevation=[REBUILD_ELEV_THRESHOLD] * TOTAL_DOMAINS,

        overwash_filter=[OVERWASH_FILTER_DEFAULT] * TOTAL_DOMAINS,
        overwash_to_dune=[OVERWASH_TO_DUNE] * TOTAL_DOMAINS,

        enable_shoreline_offset=True,
        shoreline_offset=island_offset_dam,

        nourishment_volume=0,
        nourishment_interval=None,
    )
    print("  Cascade built OK.")

    groin_cb = None
    if groin_on:
        try:
            from scripts.groin_module.hindcast_groin_test.version_control.HAT_groin_module import GroinCallback
        except ImportError as e:
            sys.exit(f"ERROR: groin run needs HAT_groin_module.py importable: {e}")

        groin_cb = GroinCallback(
            updrift_pad=_gis_to_pad(GROIN_UPDRIFT_GIS),
            downdrift_pad=_gis_to_pad(GROIN_DOWNDRIFT_GIS),
            trapping_rate_m_yr=GROIN_TRAPPING_RATE_M_YR,
            start_year=START_YEAR,
            install_year=GROIN_INSTALL_YEAR,
            n_domains=TOTAL_DOMAINS,
            deterioration_delay_years=GROIN_DETERIORATION_DELAY_YEARS,
            deterioration_mode=GROIN_DETERIORATION_MODE,
            deterioration_ramp_years=GROIN_DETERIORATION_RAMP_YEARS,
            deterioration_fraction=GROIN_DETERIORATION_FRACTION,
        )
        cascade._groin_callback = groin_cb
        print(f"  Groin attached: updrift D{GROIN_UPDRIFT_GIS} "
              f"(pad {_gis_to_pad(GROIN_UPDRIFT_GIS)}), "
              f"downdrift D{GROIN_DOWNDRIFT_GIS} (pad {_gis_to_pad(GROIN_DOWNDRIFT_GIS)}), "
              f"M={GROIN_TRAPPING_RATE_M_YR} m/yr, install {GROIN_INSTALL_YEAR}, "
              f"deterioration @ +{GROIN_DETERIORATION_DELAY_YEARS}yr "
              f"({GROIN_DETERIORATION_MODE}, floor={GROIN_DETERIORATION_FRACTION:.3f})")

    print(f"  Stepping {RUN_YEARS} years...")
    historical_nourishment_log = []
    for time_step in range(RUN_YEARS - 1):
        current_year = START_YEAR + time_step
        print(f"\r    Year {time_step + 1}/{RUN_YEARS}", end="", flush=True)

        # Reset per-step nourishment flags so nothing carries over from prior year.
        cascade.nourish_now = np.zeros(TOTAL_DOMAINS)

        # ── HISTORICAL BEACH NOURISHMENT (1971, 1973 -- see Section 2c) ────
        if current_year in historical_nourishment_on_by_year:
            nourish_now = np.asarray(
                historical_nourishment_on_by_year[current_year], dtype=float
            )
            nourish_vol = np.asarray(
                historical_nourishment_volume_by_year[current_year], dtype=float
            )

            if np.any(nourish_now == 1):
                print(f"\n  -> Applying historical BN in {current_year}:")
                cascade.nourish_now = nourish_now

                for iB3D in range(TOTAL_DOMAINS):
                    if nourish_now[iB3D] != 1:
                        continue

                    # IMPORTANT: must set Cascade's OWN nourishment_volume list,
                    # NOT cascade.nourishments[iB3D]._nourishment_volume directly.
                    # cascade.update() overwrites nourishments[iB3D].nourishment_volume
                    # from cascade._nourishment_volume[iB3D] every single call (see
                    # cascade_groin.py, right before nourishments[iB3D].update() is
                    # called) -- setting the object's attribute directly gets
                    # silently discarded the moment update() runs, which is exactly
                    # why the requested volume never reached x_s.
                    cascade.nourishment_volume[iB3D] = float(nourish_vol[iB3D])

                    gis_id = FIRST_FILE_NUMBER + (iB3D - START_REAL_INDEX)
                    vol_m3_per_m = float(nourish_vol[iB3D])
                    print(
                        f"    GIS {gis_id:3d} (pad {iB3D:3d}): "
                        f"{vol_m3_per_m:.1f} m^3/m  |  "
                        f"{vol_m3_per_m * DOMAIN_LENGTH_M:,.0f} m^3 total"
                    )
                    historical_nourishment_log.append(dict(
                        run_name                    = run_name,
                        model_year                  = current_year,
                        time_step                   = time_step,
                        padded_index                = iB3D,
                        gis_domain                  = gis_id,
                        nourishment_volume_m3_per_m = vol_m3_per_m,
                        nourishment_volume_m3_total = vol_m3_per_m * DOMAIN_LENGTH_M,
                    ))

        cascade.update()
        if getattr(cascade, "b3d_break", False):
            print(f"\n    Model stopped early at year {time_step + 1} (b3d_break).")
            break
    print("\n  Stepping complete.")

    if groin_cb is not None and len(groin_cb.year_TS) == 0:
        print("\n" + "!" * 78)
        print("WARNING: groin callback was never called. The pre-AST hook in")
        print("cascade.py is missing, so this 'groin' run is identical to no_groin.")
        print("Add the 3-line hook to cascade.py, then re-run.")
        print("!" * 78)

    run_dir = os.path.join(OUTPUT_BASE_DIR, run_name)
    os.makedirs(run_dir, exist_ok=True)
    cascade.save(run_dir)
    shoreline_m = build_shoreline_matrix(cascade, to_meters=True)
    np.save(os.path.join(run_dir, f"{run_name}_shoreline_matrix.npy"), shoreline_m)
    print(f"  Saved run to: {run_dir}   (matrix {shoreline_m.shape})")

    if MAKE_FIGURES:
        _save_run_figures(run_name, run_dir, shoreline_m)
    if MAKE_RUN_GIF:
        _save_run_gif(run_name, run_dir, shoreline_m)

    if groin_cb is not None and len(groin_cb.year_TS) > 0:
        pd.DataFrame(groin_cb.diagnostics_frame()).to_csv(
            os.path.join(run_dir, f"{run_name}_groin_diagnostics.csv"), index=False)
        print(f"  Saved groin diagnostics ({len(groin_cb.year_TS)} yrs)")

    if len(historical_nourishment_log) > 0:
        bn_df = pd.DataFrame(historical_nourishment_log)
        bn_csv = os.path.join(run_dir, f"{run_name}_historical_BN_log.csv")
        bn_df.to_csv(bn_csv, index=False)
        print(f"  Saved BN log ({len(bn_df)} events): {bn_csv}")

    delta = shoreline_m[-1, START_REAL_INDEX:END_REAL_INDEX] - \
            shoreline_m[0, START_REAL_INDEX:END_REAL_INDEX]
    print(f"  Real-domain shoreline change D{FIRST_FILE_NUMBER}-D{LAST_FILE_NUMBER} (m, raw end-start):")
    for i, d in enumerate(delta):
        print(f"    D{FIRST_FILE_NUMBER + i:<3d} {d:+.1f} m")

    return run_name


# =============================================================================
# SECTION 7: MAIN
# =============================================================================
def main():
    print("=" * 78)
    print(f"GROIN-TEST HINDCAST  {START_YEAR}-{END_YEAR}  "
          f"D{FIRST_FILE_NUMBER}-D{LAST_FILE_NUMBER}  ({TOTAL_DOMAINS} padded)")
    print(f"Run matrix: {RUN_MATRIX}")
    print("=" * 78)

    print("\nChecking inputs...")
    check_inputs_exist()
    island_offset_dam = load_island_offset_dam()
    elevation_files, dune_files = build_file_lists()

    print("\nBuilding historical nourishment schedule (1971, 1973)...")
    hist_nourish_on, hist_nourish_vol = build_nourishment_arrays_from_manual_inputs()

    produced = []
    for run_key in RUN_MATRIX:
        produced.append(run_one(run_key, island_offset_dam, elevation_files, dune_files,
                                 hist_nourish_on, hist_nourish_vol))

    print("\n" + "=" * 78)
    print("DONE. Runs produced:")
    for r in produced:
        print(f"   {r}")
    print("\nPlot with HAT_plot_groin_runs.py:")
    print(f"   RUNS = {produced}")
    print("=" * 78)


if __name__ == "__main__":
    main()
