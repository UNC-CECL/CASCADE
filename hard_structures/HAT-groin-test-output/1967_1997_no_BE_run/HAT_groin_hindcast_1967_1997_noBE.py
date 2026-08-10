"""
HAT_groin_hindcast_1967_1997_noBE.py
===============================
Groin-test hindcast (1967-1997, GIS D2-D12). Built from the proven base-run
script (HAT_base_run_1967_1997.py) -- same inputs, same conventions -- extended
to run TWO experiments and save each to its own folder:

  RUN_MATRIX = ["no_groin", "groin"]
    -> HAT_1967_1997_groinTest_no_groin   (groin OFF -- the erosive baseline)
    -> HAT_1967_1997_groinTest_groin      (groin ON  -- dipole at D5/D6)

The "no_groin" run needs nothing beyond your working base setup. The "groin" run
additionally requires:
  1. HAT_groin_module.py importable (same folder / on PYTHONPATH), and
  2. the inert pre-AST hook in cascade.py (sets cascade._groin_callback).
If the hook is missing, the groin run warns loudly (diagnostics stay empty) so
you never mistake a no-op for a real groin run.

Plot with HAT_plot_groin_runs.py:
    RUNS = ["HAT_1967_1997_groinTest_no_groin", "HAT_1967_1997_groinTest_groin"]
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
#   "no_groin"  -> HAT_1967_1997_groinTest_no_groin   (erosive baseline)
#   "groin"     -> HAT_1967_1997_groinTest_groin      (dipole at D5/D6; needs hook)
#   "groin_be"  -> HAT_1967_1997_groinTest_groin_be   (groin + background erosion)
RUN_MATRIX = ["groin"]              # <- edit this line to choose the run(s)

# Buxton groin: sits in D6 (source/accretion), starves D5 (sink/erosion).
GROIN_UPDRIFT_GIS   = 6
GROIN_DOWNDRIFT_GIS = 5
GROIN_TRAPPING_RATE_M_YR = 70.0    # M -- the single knob; tune to observed updrift
GROIN_INSTALL_YEAR  = 1970         # inert before this (free 1967-69 control window)

# Optional regional background erosion for a "groin_be" run (only used if that
# key is in RUN_MATRIX). m/yr, negative = erosive.
REGIONAL_BE_RATE_M_YR = 0

# =============================================================================
# SECTION 3: PERIOD / FILE PATHS   (identical to base run)
# =============================================================================
PROJECT_BASE_DIR   = r"/"
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, "data", "hatteras_init")
OUTPUT_BASE_DIR    = os.path.join(PROJECT_BASE_DIR, "comparison", "raw_runs")
PARAMETER_FILE     = "Hatteras-CASCADE-parameters.yaml"

START_YEAR = 1967
END_YEAR   = 1997
RUN_YEARS  = END_YEAR - START_YEAR          # 30

SEA_LEVEL_RISE_RATE = 0.004
SEA_LEVEL_CONSTANT  = True

GROIN_INIT_DIR = os.path.join(
    PROJECT_BASE_DIR, "scripts", "groin_module_noBE", "hindcast_groin_test", "groin_init",
)
STORM_FILE = os.path.join(
    GROIN_INIT_DIR, "storms", "1967_1997", "1967_1997_grointest_storms.npy",
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

RUN_NAME_SUFFIX = "noBE_70M"    # -> HAT_1967_1997_groinTest_{run_key}

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


def run_one(run_key, island_offset_dam, elevation_files, dune_files):
    groin_on = run_key in ("groin", "groin_be")
    be_on    = run_key == "groin_be"

    run_name = f"HAT_{START_YEAR}_{END_YEAR}_{RUN_NAME_SUFFIX}_{run_key}"
    print("\n" + "=" * 78)
    print(f"RUN: {run_name}   (groin={'ON' if groin_on else 'off'}, "
          f"BE={'ON' if be_on else 'off'})")
    print("=" * 78)

    be = [0.0] * TOTAL_DOMAINS
    if be_on:
        for gis in range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1):
            be[_gis_to_pad(gis)] = REGIONAL_BE_RATE_M_YR

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
        beach_nourishment_module=[False] * TOTAL_DOMAINS,
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
        )
        cascade._groin_callback = groin_cb
        print(f"  Groin attached: updrift D{GROIN_UPDRIFT_GIS} "
              f"(pad {_gis_to_pad(GROIN_UPDRIFT_GIS)}), "
              f"downdrift D{GROIN_DOWNDRIFT_GIS} (pad {_gis_to_pad(GROIN_DOWNDRIFT_GIS)}), "
              f"M={GROIN_TRAPPING_RATE_M_YR} m/yr, install {GROIN_INSTALL_YEAR}")

    print(f"  Stepping {RUN_YEARS} years...")
    for time_step in range(RUN_YEARS - 1):
        print(f"\r    Year {time_step + 1}/{RUN_YEARS}", end="", flush=True)
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

    produced = []
    for run_key in RUN_MATRIX:
        produced.append(run_one(run_key, island_offset_dam, elevation_files, dune_files))

    print("\n" + "=" * 78)
    print("DONE. Runs produced:")
    for r in produced:
        print(f"   {r}")
    print("\nPlot with HAT_plot_groin_runs.py:")
    print(f"   RUNS = {produced}")
    print("=" * 78)


if __name__ == "__main__":
    main()
