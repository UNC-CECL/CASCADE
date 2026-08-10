#!/usr/bin/env python3
"""
PEA ISLAND: RUN CASCADE + PLOT OVERWASH FLUX (Qow)
Period: 1992–2010

Modified to:
  1. Read/check BermEl and MHW from YAML
  2. Pass forced BermEl and MHW values to Cascade
  3. Print actual BermEl and MHW values inside Barrier3D after initialization
"""

import os
import sys
import yaml
import numpy as np
import matplotlib.pyplot as plt
from cascade.cascade import Cascade

# =============================================================================
# SECTION 1: DOMAIN CONFIGURATION
# =============================================================================

NUM_REAL_DOMAINS = 41
NUM_BUFFER_DOMAINS = 15

FIRST_FILE_NUMBER = 80
LAST_FILE_NUMBER = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1

TOTAL_DOMAINS = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS
START_REAL_INDEX = NUM_BUFFER_DOMAINS
END_REAL_INDEX = START_REAL_INDEX + NUM_REAL_DOMAINS

FIRST_ROAD_DOMAIN = 9
LAST_ROAD_DOMAIN = 90
START_ROAD_INDEX = (FIRST_ROAD_DOMAIN - 1) + NUM_BUFFER_DOMAINS
END_ROAD_INDEX = (LAST_ROAD_DOMAIN - 1) + NUM_BUFFER_DOMAINS + 1

print("=" * 80)
print("PEA ISLAND OVERWASH / Qow CONFIGURATION")
print("=" * 80)
print(f"Real Domains:   {NUM_REAL_DOMAINS} IDs {FIRST_FILE_NUMBER}..{LAST_FILE_NUMBER}")
print(f"Buffers:        {NUM_BUFFER_DOMAINS} each side")
print(f"TOTAL_DOMAINS:  {TOTAL_DOMAINS}")
print(f"Real index span: [{START_REAL_INDEX}..{END_REAL_INDEX - 1}]")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 2: FILE PATHS
# =============================================================================

PROJECT_BASE_DIR = r"/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE"
PEA_DATA_BASE = os.path.join(PROJECT_BASE_DIR, "data", "PeaIsland_init")
OUTPUT_BASE_DIR = os.path.join(PROJECT_BASE_DIR, "comparison", "raw_runs", "overwash")

DUNE_OFFSET_FILE = "/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE/data/PeaIsland_init/Island_offset/org/Island_Dune_Offsets_1992_2010_PADDED_71_old.csv"

STORM_FILE = os.path.join(
    PEA_DATA_BASE,
    "Storms",
    "storms_1992_2010_old.npy",
)

ROAD_SETBACK_FILE = os.path.join(
    PEA_DATA_BASE,
    "Roads",
    "raw_offset",
    "1978_road_setback_CLEAN.csv",
)

PARAMETER_FILE = os.path.join(
    PEA_DATA_BASE,
    "barrier3d-NoMngment-parameters.yaml",
)

os.chdir(PROJECT_BASE_DIR)
os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)

# =============================================================================
# SECTION 3: SIMULATION PARAMETERS — 1992–2010 ONLY
# =============================================================================

START_YEAR = 1992
END_YEAR = 2010
RUN_YEARS = END_YEAR - START_YEAR

YEAR_COLUMN_INDEX = 0
RUN_NAME = f"PEA_{START_YEAR}_{END_YEAR}_natural_Qow_1992_2010"

NUM_CORES = 4

SEA_LEVEL_RISE_RATE = 0.005
SEA_LEVEL_CONSTANT = True

# =============================================================================
# FORCE THESE VALUES INTO CASCADE
# =============================================================================
# These are in meters when passed to CASCADE.
# CASCADE/Barrier3D may internally store them in decameters.
# Example:
#   1.7 m may print as 0.17 inside Barrier3D
#   0.36 m may print as 0.036 inside Barrier3D

BERM_ELEVATION_M = 1.7
MHW_M = 0.36

ENABLE_ROADWAY_MANAGEMENT = False
ENABLE_NOURISHMENT = False
ENABLE_SANDBAG_PLACEMENT = False

DUNE_REBUILD_HEIGHT = 3.0
REBUILD_ELEV_THRESHOLD = 0.01

FIXED_WAVE_HEIGHT = 2.0
FIXED_WAVE_PERIOD = 7
FIXED_WAVE_ASYMMETRY = 0.8
FIXED_WAVE_ANGLE_HIGH_FRACTION = 0.1

ROAD_ELEVATION = 1.45
ROAD_WIDTH = 20.0

OVERWASH_FILTER = 0
OVERWASH_TO_DUNE = 9
NOURISHMENT_VOLUME = 0

LABEL_EVERY_N_DOMAINS = 2
PLOT_EVERY_N_YEARS = 1

print("Simulation Configuration:")
print(f"  Start Year: {START_YEAR}")
print(f"  End Year:   {END_YEAR}")
print(f"  RUN_YEARS:  {RUN_YEARS}")
print(f"  Storm file: {STORM_FILE}")
print(f"  Run Name:   {RUN_NAME}")
print(f"  Forced BermEl argument: {BERM_ELEVATION_M} m")
print(f"  Forced MHW argument:    {MHW_M} m")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 4: LOAD INPUT DATA
# =============================================================================

print("Loading input data...")

try:
    dune_offset_all = np.loadtxt(DUNE_OFFSET_FILE, skiprows=1, delimiter=",")
    dune_offset_dam = dune_offset_all[:, YEAR_COLUMN_INDEX] / 10.0
    print(f"✓ Loaded dune offsets: {dune_offset_dam.size} values")

    road_setbacks_raw = np.loadtxt(ROAD_SETBACK_FILE, skiprows=1, delimiter=",")
    print(f"✓ Loaded road setbacks: {road_setbacks_raw.size} values")

except FileNotFoundError as e:
    print(f"❌ Missing file: {e.filename}")
    sys.exit(1)
except Exception as e:
    print(f"❌ Error loading input data: {e}")
    sys.exit(1)

BACKGROUND_EROSION_RATES = [0.0] * TOTAL_DOMAINS

road_setbacks_full = np.zeros(TOTAL_DOMAINS)
num_road_values = min(len(road_setbacks_raw), END_ROAD_INDEX - START_ROAD_INDEX)
road_setbacks_full[
    START_ROAD_INDEX:START_ROAD_INDEX + num_road_values
] = road_setbacks_raw[:num_road_values]

ROADWAY_MANAGEMENT_ON = [ENABLE_ROADWAY_MANAGEMENT] * TOTAL_DOMAINS
NOURISHMENT_MANAGEMENT_ON = [ENABLE_NOURISHMENT] * TOTAL_DOMAINS

print(f"✓ Road setback array prepared for {TOTAL_DOMAINS} domains")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 5: ELEVATION + DUNE FILE LISTS
# =============================================================================

print("Generating elevation + dune file lists...")

ELEVATION_FILE_PATHS = []
DUNE_FILE_PATHS = []

for _ in range(START_REAL_INDEX):
    DUNE_FILE_PATHS.append(
        os.path.join(PEA_DATA_BASE, "Buffer", "sample_1_dune_Pea.npy")
    )
    ELEVATION_FILE_PATHS.append(
        os.path.join(PEA_DATA_BASE, "Buffer", "sample_1_topography_Pea.npy")
    )

for i_list in range(START_REAL_INDEX, END_REAL_INDEX):
    file_num = FIRST_FILE_NUMBER + (i_list - START_REAL_INDEX)

    DUNE_FILE_PATHS.append(
        os.path.join(
            PEA_DATA_BASE,
            "Dunes",
            "1996_dunes",
            f"resampled_1996_domains_{file_num}_dune.npy",
        )
    )

    ELEVATION_FILE_PATHS.append(
        os.path.join(
            PEA_DATA_BASE,
            "Topo",
            "2011_modified",
            f"resampled_2011_domains_{file_num}_topography.npy",
        )
    )

for _ in range(END_REAL_INDEX, TOTAL_DOMAINS):
    DUNE_FILE_PATHS.append(
        os.path.join(PEA_DATA_BASE, "Buffer", "sample_1_dune_Pea.npy")
    )
    ELEVATION_FILE_PATHS.append(
        os.path.join(PEA_DATA_BASE, "Buffer", "sample_1_topography_Pea.npy")
    )

print(f"✓ Elevation files: {len(ELEVATION_FILE_PATHS)}")
print(f"✓ Dune files:      {len(DUNE_FILE_PATHS)}")
print("=" * 80 + "\n")

# =============================================================================
# DOMAIN HELPERS
# =============================================================================

def real_domain_to_padded_index(real_domain):
    return START_REAL_INDEX + (real_domain - FIRST_FILE_NUMBER)


def padded_index_to_real_domain(padded_index):
    if START_REAL_INDEX <= padded_index < END_REAL_INDEX:
        return FIRST_FILE_NUMBER + (padded_index - START_REAL_INDEX)
    return None


# =============================================================================
# BERM / MHW CHECK HELPER
# =============================================================================

def print_berm_mhw_check(cascade=None, note=""):
    """
    Print BermEl and MHW from:
      1. YAML file
      2. Script-forced constants
      3. Actual Barrier3D objects after Cascade initialization
    """

    print("\n" + "=" * 80)
    print(f"BERM / MHW CHECK {note}")
    print("=" * 80)

    try:
        with open(PARAMETER_FILE, "r") as f:
            yaml_params = yaml.safe_load(f)

        print(f"YAML file: {PARAMETER_FILE}")
        print(f"YAML BermEl: {yaml_params.get('BermEl')}")
        print(f"YAML MHW:    {yaml_params.get('MHW')}")

    except Exception as e:
        print(f"Could not read YAML file: {e}")

    print("\nScript-forced values passed to Cascade:")
    print(f"  BERM_ELEVATION_M = {BERM_ELEVATION_M}")
    print(f"  MHW_M            = {MHW_M}")

    if cascade is not None:
        print("\nActual values inside Barrier3D objects:")

        check_indices = [
            0,
            START_REAL_INDEX,
            real_domain_to_padded_index(100),
            real_domain_to_padded_index(110),
            real_domain_to_padded_index(120),
            TOTAL_DOMAINS - 1,
        ]

        for idx in check_indices:
            b3d = cascade.barrier3d[idx]

            berm_val = getattr(b3d, "BermEl", None)
            mhw_val = getattr(b3d, "MHW", None)

            real_domain = padded_index_to_real_domain(idx)

            print(
                f"  padded index {idx:02d}, "
                f"real domain {real_domain}: "
                f"BermEl = {berm_val}, "
                f"MHW = {mhw_val}"
            )

        print("\nUnit note:")
        print("  If BermEl prints as 0.17, that usually means 1.7 m stored as decameters.")
        print("  If MHW prints as 0.036, that usually means 0.36 m stored as decameters.")
        print("  If BermEl prints as 1.7 and MHW as 0.36, then this CASCADE version stores them in meters.")

    print("=" * 80 + "\n")


# =============================================================================
# CASCADE RUNNER
# =============================================================================

def run_cascade_simulation(
    nt,
    name,
    storm_file,
    alongshore_section_count,
    num_cores,
    rmin,
    rmax,
    elevation_file,
    dune_file,
    dune_design_elevation,
    dune_minimum_elevation,
    road_ele,
    road_width,
    road_setback,
    overwash_filter,
    overwash_to_dune,
    nourishment_volume,
    background_erosion,
    roadway_management_on,
    beach_dune_manager_on,
    sea_level_rise_rate,
    sea_level_constant,
    enable_shoreline_offset,
    shoreline_offset,
    wave_height,
    wave_period,
    wave_asymmetry,
    wave_angle_high_fraction,
):
    datadir = PEA_DATA_BASE

    print_berm_mhw_check(cascade=None, note="BEFORE CASCADE INIT")

    cascade = Cascade(
        datadir,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file=PARAMETER_FILE,

        # ---------------------------------------------------------------------
        # IMPORTANT:
        # These force the values used by CASCADE, matching your other run.
        # ---------------------------------------------------------------------
        berm_elevation=BERM_ELEVATION_M,
        MHW=MHW_M,

        wave_height=wave_height,
        wave_period=wave_period,
        wave_asymmetry=wave_asymmetry,
        wave_angle_high_fraction=wave_angle_high_fraction,

        sea_level_rise_rate=sea_level_rise_rate,
        sea_level_rise_constant=sea_level_constant,

        background_erosion=background_erosion,
        alongshore_section_count=alongshore_section_count,
        time_step_count=nt,

        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,

        num_cores=num_cores,

        roadway_management_module=roadway_management_on,
        beach_nourishment_module=beach_dune_manager_on,
        alongshore_transport_module=True,
        community_economics_module=False,

        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,

        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,

        overwash_filter=overwash_filter,
        overwash_to_dune=overwash_to_dune,

        enable_shoreline_offset=enable_shoreline_offset,
        shoreline_offset=shoreline_offset,

        nourishment_volume=nourishment_volume,
        nourishment_interval=None,
    )

    print_berm_mhw_check(cascade=cascade, note="AFTER CASCADE INIT")

    # Required for your CASCADE version:
    # update() expects sandbag settings as lists, even when sandbags are off.
    cascade._sandbag_management_on = [False] * alongshore_section_count
    cascade._sandbag_elevation = [0.0] * alongshore_section_count

    for time_step in range(nt - 1):
        print(f"\rYear {time_step + 1}/{nt}", end="", flush=True)
        cascade.update()

        if getattr(cascade, "b3d_break", False):
            print(f"\n⚠️ stopped at year {time_step + 1} because b3d_break=True")
            break

    print_berm_mhw_check(cascade=cascade, note="AFTER FULL RUN")

    save_path = os.path.join(OUTPUT_BASE_DIR, name)
    os.makedirs(save_path, exist_ok=True)
    cascade.save(save_path)
    print(f"\n✓ saved run: {save_path}")

    return cascade


# =============================================================================
# Qow HELPERS
# =============================================================================

def get_Qow_ts(b3d):
    for name in ("QowTS", "Qow_TS", "_QowTS", "_Qow_TS"):
        if hasattr(b3d, name):
            arr = np.asarray(getattr(b3d, name), dtype=float).squeeze()
            if arr.ndim == 1 and arr.size > 1:
                return arr, name

    if hasattr(b3d, "Qow"):
        arr = np.asarray(getattr(b3d, "Qow"), dtype=float).squeeze()
        if arr.ndim == 1 and arr.size > 1:
            return arr, "Qow"

    raise AttributeError(
        "No Qow time series found on this Barrier3D object. "
        "Try: print([k for k in dir(cascade.barrier3d[15]) "
        "if 'qow' in k.lower() or 'ow' in k.lower()])"
    )


def build_Qow_matrix(cascade, real_only=True, use_model_domain_numbers=False):
    b3d_list = cascade.barrier3d

    if real_only:
        dom_ids = list(range(START_REAL_INDEX, END_REAL_INDEX))
    else:
        dom_ids = list(range(len(b3d_list)))

    ts0, attr = get_Qow_ts(b3d_list[dom_ids[0]])
    nt = ts0.size

    Q = np.zeros((nt, len(dom_ids)), dtype=float)
    Q[:, 0] = ts0

    for j, d in enumerate(dom_ids[1:], start=1):
        ts, _ = get_Qow_ts(b3d_list[d])
        if ts.size != nt:
            raise ValueError(f"Domain {d}: Qow time series length {ts.size} != {nt}")
        Q[:, j] = ts

    years = START_YEAR + np.arange(nt)

    if real_only:
        if use_model_domain_numbers:
            domain_labels = list(range(1, len(dom_ids) + 1))
        else:
            domain_labels = list(range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1))
    else:
        domain_labels = dom_ids

    return years, Q, attr, domain_labels


# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def plot_Qow_lines_by_year(
    cascade,
    real_only=True,
    use_model_domain_numbers=False,
    label_every_n_domains=2,
    plot_every_n_years=1,
    save_fig=False,
    out_path=None,
    cmap_name="Greens",   # options: "viridis", "plasma", "Greens", "Blues", "coolwarm"
):
    years, Q, attr, domain_labels = build_Qow_matrix(
        cascade,
        real_only=real_only,
        use_model_domain_numbers=use_model_domain_numbers,
    )

    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)

    # Years that will actually be plotted
    plot_indices = list(range(0, len(years), plot_every_n_years))

    # Create gradient colors
    cmap = plt.get_cmap(cmap_name)
    colors = cmap(np.linspace(0.15, 0.95, len(plot_indices)))

    for color, i in zip(colors, plot_indices):
        ax.plot(
            domain_labels,
            Q[i, :],
            lw=1.8,
            alpha=0.9,
            color=color,
            label=str(years[i]),
        )

    ax.set_title(f"Pea Island Overwash Flux by Domain, {START_YEAR}–{END_YEAR} ({attr})")

    ax.set_xlabel(
        "Model Domain #" if (use_model_domain_numbers and real_only)
        else "Real Domain ID" if real_only
        else "Alongshore Domain Index"
    )

    ax.set_ylabel("Overwash Flux, Qow")
    ax.grid(alpha=0.3)

    xticks = domain_labels[::label_every_n_domains]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, rotation=45)

    ax.legend(
        title="Year",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=8,
    )

    if save_fig and out_path is not None:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"✓ saved figure: {out_path}")

    plt.show()

    print(f"Used Qow attribute: {attr}")
    print(f"Qow matrix shape: {Q.shape}")
    print(f"Qow range: {np.nanmin(Q):.4f} .. {np.nanmax(Q):.4f}")

def plot_Qow_mean_and_sum(cascade, real_only=True, use_model_domain_numbers=False):
    years, Q, attr, domain_labels = build_Qow_matrix(
        cascade,
        real_only=real_only,
        use_model_domain_numbers=use_model_domain_numbers,
    )

    q_mean = np.nanmean(Q, axis=1)
    q_sum = np.nansum(Q, axis=1)

    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)

    ax.plot(years, q_mean, lw=2, label="Mean Qow")
    ax.plot(years, q_sum, lw=2, label="Sum Qow")

    ax.set_title(f"Pea Island Overwash Flux Over Time, {START_YEAR}–{END_YEAR} ({attr})")
    ax.set_xlabel("Year")
    ax.set_ylabel("Overwash Flux, Qow")
    ax.grid(alpha=0.3)
    ax.legend()

    plt.show()


def plot_single_year_Qow(
    cascade,
    year_to_plot,
    real_only=True,
    use_model_domain_numbers=False,
):
    years, Q, attr, domain_labels = build_Qow_matrix(
        cascade,
        real_only=real_only,
        use_model_domain_numbers=use_model_domain_numbers,
    )

    if year_to_plot not in years:
        raise ValueError(
            f"Year {year_to_plot} not found. "
            f"Available years: {years[0]} to {years[-1]}"
        )

    idx = np.where(years == year_to_plot)[0][0]

    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)

    ax.plot(domain_labels, Q[idx, :], lw=2, marker="o")
    ax.set_title(f"Pea Island Overwash Flux by Domain in {year_to_plot} ({attr})")

    ax.set_xlabel(
        "Model Domain #" if (use_model_domain_numbers and real_only)
        else "Real Domain ID" if real_only
        else "Alongshore Domain Index"
    )

    ax.set_ylabel("Overwash Flux, Qow")
    ax.grid(alpha=0.3)

    xticks = domain_labels[::2]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, rotation=45)

    plt.show()


def save_Qow_csv(cascade, out_csv, real_only=True, use_model_domain_numbers=False):
    years, Q, attr, domain_labels = build_Qow_matrix(
        cascade,
        real_only=real_only,
        use_model_domain_numbers=use_model_domain_numbers,
    )

    header = ["Year"] + [f"Domain_{d}" for d in domain_labels]
    data = np.column_stack([years, Q])

    os.makedirs(os.path.dirname(out_csv), exist_ok=True)

    with open(out_csv, "w") as f:
        f.write(",".join(header) + "\n")
        for row in data:
            f.write(",".join(map(str, row)) + "\n")

    print(f"✓ Qow CSV saved to: {out_csv}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    RMIN = [0.55] * TOTAL_DOMAINS
    RMAX = [0.95] * TOTAL_DOMAINS

    DUNE_DESIGN_ELEVATION = [DUNE_REBUILD_HEIGHT] * TOTAL_DOMAINS
    DUNE_MINIMUM_ELEVATION = [REBUILD_ELEV_THRESHOLD] * TOTAL_DOMAINS

    cascade = run_cascade_simulation(
        nt=RUN_YEARS,
        name=RUN_NAME,
        storm_file=STORM_FILE,
        alongshore_section_count=TOTAL_DOMAINS,
        num_cores=NUM_CORES,

        rmin=RMIN,
        rmax=RMAX,

        elevation_file=ELEVATION_FILE_PATHS,
        dune_file=DUNE_FILE_PATHS,

        dune_design_elevation=DUNE_DESIGN_ELEVATION,
        dune_minimum_elevation=DUNE_MINIMUM_ELEVATION,

        road_ele=ROAD_ELEVATION,
        road_width=ROAD_WIDTH,
        road_setback=road_setbacks_full,

        overwash_filter=OVERWASH_FILTER,
        overwash_to_dune=OVERWASH_TO_DUNE,

        nourishment_volume=NOURISHMENT_VOLUME,
        background_erosion=BACKGROUND_EROSION_RATES,

        roadway_management_on=ROADWAY_MANAGEMENT_ON,
        beach_dune_manager_on=NOURISHMENT_MANAGEMENT_ON,

        sea_level_rise_rate=SEA_LEVEL_RISE_RATE,
        sea_level_constant=SEA_LEVEL_CONSTANT,

        enable_shoreline_offset=True,
        shoreline_offset=dune_offset_dam,

        wave_height=FIXED_WAVE_HEIGHT,
        wave_period=FIXED_WAVE_PERIOD,
        wave_asymmetry=FIXED_WAVE_ASYMMETRY,
        wave_angle_high_fraction=FIXED_WAVE_ANGLE_HIGH_FRACTION,
    )

    plot_Qow_mean_and_sum(
        cascade,
        real_only=True,
        use_model_domain_numbers=False,
    )

    out_fig = os.path.join(
        OUTPUT_BASE_DIR,
        f"PEA_Qow_lines_{START_YEAR}_{END_YEAR}.png",
    )

    plot_Qow_lines_by_year(
        cascade,
        real_only=True,
        use_model_domain_numbers=False,
        label_every_n_domains=LABEL_EVERY_N_DOMAINS,
        plot_every_n_years=PLOT_EVERY_N_YEARS,
        save_fig=True,
        out_path=out_fig,
    )

    # Optional: plot only one year
    # plot_single_year_Qow(
    #     cascade,
    #     year_to_plot=1992,
    #     real_only=True,
    #     use_model_domain_numbers=False,
    # )

    out_csv = os.path.join(
        OUTPUT_BASE_DIR,
        f"PEA_Qow_by_year_domain_{START_YEAR}_{END_YEAR}.csv",
    )

    save_Qow_csv(
        cascade,
        out_csv=out_csv,
        real_only=True,
        use_model_domain_numbers=False,
    )


if __name__ == "__main__":
    main()