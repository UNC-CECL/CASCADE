"""
The 1967-2017 groin rig (GIS D2-D12, 41 padded domains): no-groin and groin runs with the 1995-2003 deterioration ramp.

    python HAT_groin_hindcast_1967_2017.py

Runs every key in RUN_MATRIX on the 1967-2017 storm file and island offset
under HAT-groin-buxton-input/groin_init/, with the 1971 and 1973 fills on
every run; writes each run to output/calibration/groin_rig/<run>/. The sweep,
its worker and the edge solve import this file. Needs CASCADE (cascade_groin).
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""

import os
import pathlib
import sys
import numpy as np
import pandas as pd


# --- CONFIG ------------------------------------------------------------------
# Section 1: which Cascade (the sandbox cascade_groin.py carries the groin hook; README)

USE_SANDBOX_CASCADE = True

# Section 2: domains (D2-D12)

NUM_REAL_DOMAINS   = 11
NUM_BUFFER_DOMAINS = 15
FIRST_FILE_NUMBER  = 2
LAST_FILE_NUMBER   = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1     # 12
TOTAL_DOMAINS      = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS  # 41
START_REAL_INDEX   = NUM_BUFFER_DOMAINS                            # 15
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS           # 26

# Section 3: which runs, and the groin

# Run keys, in order: "no_groin", "groin" (needs the hook), "groin_be" (groin + regional BE)
RUN_MATRIX = ["no_groin", "groin"]  # <- edit this line to choose the run(s)

# Buxton groin: sits in D6 (source/accretion), starves D5 (sink/erosion).
GROIN_UPDRIFT_GIS   = 6
GROIN_DOWNDRIFT_GIS = 5
GROIN_TRAPPING_RATE_M_YR = 60.0    # M -- the single knob; tune to observed updrift
GROIN_INSTALL_YEAR  = 1970         # inert before this (free 1967-69 control window)

# Deterioration: last repair 1995, linear decline, Isabel 2003 locks it in; delay counts from install (README)
GROIN_DETERIORATION_DELAY_YEARS = 1995 - GROIN_INSTALL_YEAR   # = 25
GROIN_DETERIORATION_MODE        = "linear_ramp"
GROIN_DETERIORATION_RAMP_YEARS  = 2003 - 1995                  # = 8
# Floor fraction: the decided pair (2026-08-30); the sweep overrides it per cell (README)
GROIN_DETERIORATION_FRACTION    = 0.60   # decided pair, 2026-08-30

# Regional background erosion for a "groin_be" run, m/yr, negative = erosive
REGIONAL_BE_RATE_M_YR = 0

# Section 3b: edge source/sink correction at D2 and D12 (structural, every run; README)

APPLY_EDGE_BE_CORRECTION = True
EDGE_BE_RATES_GIS = {
    2:  5.0,   # D2  -- south edge, against buffer (solve for this)
    12: 10.0,   # D12 -- north edge, against buffer (solve for this)
}

# Section 3c: historical nourishment, 1971 and 1973, on every run (sources and domains: README)

ENABLE_HISTORICAL_NOURISHMENT = True

_CY_TO_M3       = 0.764555   # cubic yards -> cubic metres
DOMAIN_LENGTH_M = 500        # alongshore width of one CASCADE domain (m)

# Fill years, sorted: 1971 point source near D8, 1973 template D6-D10 (volumes are m^3, not cy)
HAT_BN_YEARS = [1971, 1973]

# GIS domain -> [m^3 in 1971, m^3 in 1973]; 0 = not nourished that year
HAT_BN_VOLUME_BY_DOMAIN = {
    # 1971: 200,000 cy, all in D8 (Hatteras Court Motel)
     8: [round(200_000 * _CY_TO_M3, 1), 0],
    # 1973: 1,300,000 cy over D6-D10, 198,784.3 m^3 per domain
     6: [0, round(1_300_000 / 5 * _CY_TO_M3, 1)],
     7: [0, round(1_300_000 / 5 * _CY_TO_M3, 1)],
     9: [0, round(1_300_000 / 5 * _CY_TO_M3, 1)],
    10: [0, round(1_300_000 / 5 * _CY_TO_M3, 1)],
}

# Section 4: period and file paths

# Repo root, found by searching upward (ORGANIZATION.md rule 5)
PROJECT_BASE_DIR = str(next(
    p for p in pathlib.Path(__file__).resolve().parents
    if (p / "pyproject.toml").exists()))
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, "data", "hatteras_init")
# Rig runs file under calibration/groin_rig, never raw_runs (README)
OUTPUT_BASE_DIR    = os.path.join(PROJECT_BASE_DIR, "output", "calibration", "groin_rig")
PARAMETER_FILE     = "Hatteras-CASCADE-parameters.yaml"

START_YEAR = 1967
END_YEAR   = 2018   # EXCLUSIVE convention (RUN_YEARS = END_YEAR - START_YEAR)
                     # -> 51 model years, final simulated year = 2017
RUN_YEARS  = END_YEAR - START_YEAR          # 51

SEA_LEVEL_RISE_RATE = 0.004
SEA_LEVEL_CONSTANT  = True

GROIN_INIT_DIR = os.path.join(
    PROJECT_BASE_DIR, "hard-structures", "groin", "HAT-groin-buxton-input", "groin_init",
)
STORM_FILE = os.path.join(
    GROIN_INIT_DIR, "storms", "1967_2017",
    "1967_2017_groin_storms.npy",
)
ISLAND_OFFSET_FILE = os.path.join(
    GROIN_INIT_DIR, "island_offset",
    "Island_Dune_Offsets_1967_D2_D12_PADDED_41.csv",
)

# Topography product: 1984-start, the one the production period-1 groin fit reads (README)
RIG_TOPO_PRODUCT = "1984-start"

TOPO_DUNE_INIT_YEAR = "2009"   # legacy label, no longer used to build filenames
TOPO_DUNE_SUBFOLDER = "2009"

# Section 5: simulation parameters

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

# Figures and a shoreline GIF saved into each run's folder when it finishes
MAKE_FIGURES = True
MAKE_RUN_GIF = True     # animate the modeled shoreline evolving over the run
# -----------------------------------------------------------------------------


# Pick the Cascade build
if USE_SANDBOX_CASCADE:
    from cascade.cascade_groin import Cascade   # hooked sandbox copy
else:
    from cascade import Cascade                 # real package (hook folded in)


# GIS domain id -> padded index (D2->15, D5->18, D6->19, D12->25)
def _gis_to_pad(gis_id):
    return START_REAL_INDEX + (gis_id - FIRST_FILE_NUMBER)


# D8 takes both years: the 1971 point source and its share of the 1973 template
HAT_BN_VOLUME_BY_DOMAIN[8][1] = round(1_300_000 / 5 * _CY_TO_M3, 1)

if not ENABLE_HISTORICAL_NOURISHMENT:
    HAT_BN_YEARS = []
    HAT_BN_VOLUME_BY_DOMAIN = {}

# Per-domain beach_nourishment_module flag: True where a fill lands in this window
NOURISHMENT_MANAGEMENT_ON = [False] * TOTAL_DOMAINS
for _gis_id in HAT_BN_VOLUME_BY_DOMAIN:
    _pad_idx = _gis_to_pad(_gis_id)
    if 0 <= _pad_idx < TOTAL_DOMAINS:
        NOURISHMENT_MANAGEMENT_ON[_pad_idx] = True

os.chdir(PROJECT_BASE_DIR)
os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)


# Stop if the storm, offset or parameter file is missing
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


# Island offset per padded domain, metres in the file, decameters out
def load_island_offset_dam():
    offset_all = np.loadtxt(ISLAND_OFFSET_FILE, skiprows=1, delimiter=",")
    offset_dam = offset_all / 10.0
    if offset_dam.size != TOTAL_DOMAINS:
        sys.exit(f"ERROR: offset has {offset_dam.size} values, expected {TOTAL_DOMAINS}.")
    print(f"  Loaded island offsets: {offset_dam.size} domains (dam)")
    return list(offset_dam)


# Padded elevation and dune file lists, resolved through hat_topo_version, never pinned (README)
def build_file_lists():
    import sys as _sys
    _scripts = os.path.join(PROJECT_BASE_DIR, "scripts")
    if _scripts not in _sys.path:
        _sys.path.insert(0, _scripts)
    from site_layer.hat_topo_version import topo_dirs, array_name, BUFFER_DIR
    topo_dir, dune_dir, _topo_run = topo_dirs(RIG_TOPO_PRODUCT)
    print(f"  topography: {RIG_TOPO_PRODUCT}/{_topo_run}")
    buf_dune = os.path.join(str(BUFFER_DIR), "sample_1_dune.npy")
    buf_elev = os.path.join(str(BUFFER_DIR), "sample_1_topography.npy")

    elev, dune = [], []
    for _ in range(START_REAL_INDEX):
        dune.append(buf_dune)
        elev.append(buf_elev)
    for i_list in range(START_REAL_INDEX, END_REAL_INDEX):
        file_num = FIRST_FILE_NUMBER + (i_list - START_REAL_INDEX)
        dune.append(os.path.join(str(dune_dir),
                                 array_name("dune", file_num)))
        elev.append(os.path.join(str(topo_dir),
                                 array_name("topography", file_num)))
    for _ in range(END_REAL_INDEX, TOTAL_DOMAINS):
        dune.append(buf_dune)
        elev.append(buf_elev)
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


# Per-year nourishment on/off and m^3/m arrays from the Section 3c schedule
def build_nourishment_arrays_from_manual_inputs():
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


# A Barrier3D object's shoreline time series, under either attribute name
def get_x_s_TS(b3d):
    if hasattr(b3d, "x_s_TS"):
        return np.asarray(b3d.x_s_TS, dtype=float)
    if hasattr(b3d, "_x_s_TS"):
        return np.asarray(b3d._x_s_TS, dtype=float)
    raise AttributeError("No x_s_TS / _x_s_TS on Barrier3D object.")


# Shoreline matrix [year x padded domain], metres by default
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


# Save the run's figures with HAT_plot_groin_runs.py's functions; never blocks the run
def _save_run_figures(run_name, run_dir, shoreline_m):
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


# Animate D2-D12 year by year against the year-0 planform, ocean at bottom
def _save_run_gif(run_name, run_dir, shoreline_m):
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

    # One frame: year-0 reference, this year, shading between, the groin line
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
        ax.set_ylabel(f"Cross-shore position (m, rel. {START_YEAR} mean)  {up_word} ▲")
        ax.set_title(f"{run_name}  —  {START_YEAR + t}")
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


# One run: build Cascade, attach the groin if asked, step with the fills, save matrix, figures, logs
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
        # The shared groin model in the package, the same one production fits (README)
        try:
            from cascade.groin import GroinCallback
        except ImportError as e:
            sys.exit(f"ERROR: groin run needs cascade.groin importable: {e}")

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

        # Historical fills, 1971 and 1973
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

                    # Cascade's own nourishment_volume list; update() overwrites the object's attribute (README)
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


# Run: check inputs, build offsets, file lists and the fill schedule, then every key in RUN_MATRIX
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
