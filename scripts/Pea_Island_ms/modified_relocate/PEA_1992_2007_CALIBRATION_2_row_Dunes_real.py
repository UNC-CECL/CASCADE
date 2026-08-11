#!/usr/bin/env python3
"""
PEA ISLAND: 1992–2007 CALIBRATION
Natural vs Roadway Management vs Roadway + Historical BN
Historical road relocation uses a one-step request queued on RoadwayManager.

This version adds diagnostics for threshold-based beach nourishment:

1. Historical BN is still applied from manual year/domain/volume arrays.
2. Automatic beach-width threshold BN can be turned ON/OFF with:
       APPLY_AUTOMATIC_BN_THRESHOLD
3. The code logs:
       - where/when beach width falls below threshold
       - whether threshold BN was actually applied to CASCADE
       - how many managed domains triggered each year
4. The code also plots:
       - yearly relative shoreline change + BN
       - yearly absolute shoreline position + BN
"""

import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml
import inspect
import cascade.roadway_manager as roadway_module
import cascade.cascade as cascade_module
import cascade.brie_coupler as brie_coupler



def get_barrier3d_array(model, public_name, private_name):
    """Read a Barrier3D array using either its public or private name."""
    if hasattr(model, public_name):
        return np.asarray(getattr(model, public_name))

    if hasattr(model, private_name):
        return np.asarray(getattr(model, private_name))

    raise AttributeError(
        f"Could not find {public_name!r} or {private_name!r} "
        f"inside {type(model).__name__}"
    )


def check_barrier3d_orientation(
        cascade_model,
        barrier3d_index,
        saved_topography_path,
):
    """
    Compare a saved topography .npy file with the array loaded by Barrier3D.

    Expected saved shape:
        (cross_shore, alongshore), such as (200, 50)
    """

    saved_path = Path(saved_topography_path)

    if not saved_path.is_file():
        raise FileNotFoundError(
            f"Topography file does not exist:\n{saved_path}"
        )

    saved = np.load(saved_path)

    # Some CASCADE versions use barrier3d; others may use _barrier3d.
    if hasattr(cascade_model, "barrier3d"):
        barrier_models = cascade_model.barrier3d
    elif hasattr(cascade_model, "_barrier3d"):
        barrier_models = cascade_model._barrier3d
    else:
        raise AttributeError(
            "Could not find 'barrier3d' or '_barrier3d' in the CASCADE model."
        )

    b3d = barrier_models[barrier3d_index]

    loaded = get_barrier3d_array(
        b3d,
        public_name="InteriorDomain",
        private_name="_InteriorDomain",
    )

    print("\n" + "=" * 76)
    print("BARRIER3D ORIENTATION VERIFICATION")
    print("=" * 76)
    print(f"Input file:       {saved_path.name}")
    print(f"Barrier3D index:  {barrier3d_index}")
    print(f"Saved shape:      {saved.shape}")
    print(f"Loaded shape:     {loaded.shape}")

    candidates = {
        "unchanged": saved,
        "alongshore reversed": saved[:, ::-1],
        "cross-shore reversed": saved[::-1, :],
        "both axes reversed": saved[::-1, ::-1],
        "transposed": saved.T,
    }

    matched = []

    print("\nComparison with the Barrier3D array:")

    for description, candidate in candidates.items():
        is_match = (
            candidate.shape == loaded.shape
            and np.allclose(
                candidate,
                loaded,
                equal_nan=True,
                rtol=1e-8,
                atol=1e-10,
            )
        )

        print(f"  {description:26s}: {is_match}")

        if is_match:
            matched.append(description)

    if matched == ["unchanged"]:
        print("\nRESULT: Barrier3D preserved the saved array orientation.")

    elif "alongshore reversed" in matched:
        print(
            "\nRESULT: The array loaded by Barrier3D is reversed "
            "alongshore relative to the saved file."
        )

    elif not matched:
        print(
            "\nRESULT: The loaded Barrier3D array does not exactly match any "
            "simple flip or transpose of the saved input."
        )

    else:
        print(f"\nRESULT: Matching transformation(s): {matched}")

    print("=" * 76)

# Ensure Cascade constructs manager objects from the modified roadway module.
cascade_module.RoadwayManager = roadway_module.RoadwayManager
Cascade = cascade_module.Cascade

print("\nCASCADE ROADWAY IMPORT CHECK")
print("cascade.py:")
print(cascade_module.__file__)
print("roadway_manager.py:")
print(inspect.getfile(roadway_module.RoadwayManager))
print(
    "Cascade uses modified RoadwayManager:",
    cascade_module.RoadwayManager is roadway_module.RoadwayManager,
)
print(
    "request_relocation available:",
    hasattr(cascade_module.RoadwayManager, "request_relocation"),
)
print(
    "request_relocation signature:",
    inspect.signature(cascade_module.RoadwayManager.request_relocation)
    if hasattr(cascade_module.RoadwayManager, "request_relocation")
    else "MISSING",
)
print()


# =============================================================================
# CASCADE / BARRIER3D YAML-PATH COMPATIBILITY FIX
# =============================================================================

def set_yaml_compatible(variable, value, file_name):
    """Update a Barrier3D YAML file when CASCADE supplies either a path or prefix.

    Some CASCADE versions pass ``parameter_file`` as a Barrier3D prefix, while
    ``brie_coupler.set_yaml`` attempts to open that prefix as a literal file.
    Barrier3D itself later appends ``-parameters.yaml``. This replacement makes
    both components resolve the same file.
    """
    candidate = os.fspath(file_name)

    if os.path.isfile(candidate):
        resolved_file = candidate
    elif os.path.isfile(candidate + "-parameters.yaml"):
        resolved_file = candidate + "-parameters.yaml"
    else:
        raise FileNotFoundError(
            "Could not locate the Barrier3D parameter file.\n"
            f"Checked literal path: {candidate}\n"
            f"Checked prefix path:  {candidate}-parameters.yaml"
        )

    with open(resolved_file, "r") as f:
        yaml_data = yaml.safe_load(f)

    if yaml_data is None:
        yaml_data = {}

    yaml_data[variable] = value

    with open(resolved_file, "w") as f:
        yaml.safe_dump(yaml_data, f, sort_keys=False)


# Patch the function used by initialize_equal() before any Cascade is created.
brie_coupler.set_yaml = set_yaml_compatible


# =============================================================================
# SECTION 1: DOMAIN CONFIGURATION
# =============================================================================

BERM_ELEVATION_M = 1.7
MHW_M = 0.36

NUM_REAL_DOMAINS = 40
NUM_BUFFER_DOMAINS_LEFT = 71
NUM_BUFFER_DOMAINS_RIGHT = 0

FIRST_FILE_NUMBER = 80
LAST_FILE_NUMBER = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1

TOTAL_DOMAINS = NUM_BUFFER_DOMAINS_LEFT + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS_RIGHT

START_REAL_INDEX = NUM_BUFFER_DOMAINS_LEFT
END_REAL_INDEX = START_REAL_INDEX + NUM_REAL_DOMAINS

print("=" * 80)
print("PEA ISLAND CASCADE DOMAIN CONFIGURATION")
print("=" * 80)
print(f"Real Domains:   {NUM_REAL_DOMAINS} IDs {FIRST_FILE_NUMBER}..{LAST_FILE_NUMBER}")
print(f"Buffers:        left={NUM_BUFFER_DOMAINS_LEFT}, right={NUM_BUFFER_DOMAINS_RIGHT}")
print(f"TOTAL_DOMAINS:  {TOTAL_DOMAINS}")
print(f"Real index span: [{START_REAL_INDEX}..{END_REAL_INDEX - 1}]")
print("=" * 80 + "\n")


# =============================================================================
# SECTION 2: FILE PATHS
# =============================================================================

PROJECT_BASE_DIR = r"/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE"
PEA_DATA_BASE = os.path.join(PROJECT_BASE_DIR, "data", "PeaIsland_init")
OUTPUT_BASE_DIR = os.path.join(PROJECT_BASE_DIR, "output", "raw_runs", "HannahParam_storm_Dunemaxel7", "1992_2007","CMP_VS_CASCADE", "2real_row_dune")

DUNE_OFFSET_FILE = os.path.join(
    PEA_DATA_BASE,
    "Island_offset",
    "New_offsets_hannah_base",
    "new",
    "Whole_Island_Dune_Offsets_1992_CASCADE_Input_unpadded.csv",
)

STORM_FILE = os.path.join(
    PEA_DATA_BASE,
    "Storms",
    "1992_2007",
    "storms_1992_2007.npy",
)


PARAMETER_FILE = r"/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE/data/PeaIsland_init/barrier3d-NoMngment-parameters.yaml"

# Domain-averaged CMP shoreline-change rates for 1992–2024.
# Update only this path if your CMP file has a different filename/location.
# CSV and Excel files are both supported.
CMP_FILE = r"/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE/scripts/input_preperation/Shoreline_changes/CMP_CASCADE_Domain_LRR_1992_2024_Real_Domains_80_119.csv"

# For Excel input, use a sheet name or zero-based sheet index. Ignored for CSV.
CMP_SHEET_NAME = 0

# Preferred column names. The loader also checks common aliases.
CMP_DOMAIN_COLUMN = "domain_id"
CMP_RATE_COLUMN = "mean_LRR_original_CMP_sign_m_yr"
CMP_STD_COLUMN = "std_dev_m_yr"
CMP_COUNT_COLUMN = "CMP_transect_count"

# The original CMP sign convention is opposite to the CASCADE plotting convention.
# Flip it once so accretion is positive and erosion is negative.
CMP_FLIP_SIGN = True

os.chdir(PROJECT_BASE_DIR)
os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)


def create_run_parameter_file(datadir, run_name, dmaxel_m):
    """Create a protected Barrier3D YAML file for one CASCADE run.

    Barrier3D expects the actual file to be named::

        <parameter_prefix>-parameters.yaml

    CASCADE must receive only ``parameter_prefix``. It must not receive the
    complete YAML path or a string ending in ``.yaml``.
    """
    if not os.path.isfile(PARAMETER_FILE):
        raise FileNotFoundError(f"Master YAML file not found: {PARAMETER_FILE}")

    with open(PARAMETER_FILE, "r") as f:
        params = yaml.safe_load(f)

    if params is None:
        params = {}

    # Reset the physical input values before every run.
    params["MHW"] = float(MHW_M)
    params["BermEl"] = float(BERM_ELEVATION_M)
    params["Dmaxel"] = float(dmaxel_m)

    # This is a SHORT PREFIX only. Do not add a directory, '.yaml', or
    # '-parameters' here. Barrier3D adds '-parameters.yaml' internally.
    run_parameter_prefix = f"{run_name}_barrier3d"

    # The actual YAML must be located inside datadir because Barrier3D joins
    # datadir and the parameter prefix when loading inputs.
    run_parameter_file = os.path.join(
        datadir,
        run_parameter_prefix + "-parameters.yaml",
    )

    with open(run_parameter_file, "w") as f:
        yaml.safe_dump(
            params,
            f,
            sort_keys=False,
            default_flow_style=False,
        )

    if run_parameter_prefix.endswith(".yaml"):
        raise RuntimeError(
            "The Barrier3D parameter prefix must not contain '.yaml': "
            f"{run_parameter_prefix}"
        )

    if run_parameter_prefix.endswith("-parameters"):
        raise RuntimeError(
            "The Barrier3D parameter prefix must not end with '-parameters': "
            f"{run_parameter_prefix}"
        )

    if not os.path.isfile(run_parameter_file):
        raise FileNotFoundError(
            f"Run-specific YAML was not created: {run_parameter_file}"
        )

    return run_parameter_prefix, run_parameter_file


# =============================================================================
# SECTION 3: SIMULATION PARAMETERS
# =============================================================================

START_YEAR = 1992
END_YEAR = 2007
ANNUAL_UPDATE_COUNT = END_YEAR - START_YEAR + 1
TIME_STEP_COUNT = ANNUAL_UPDATE_COUNT + 1
RUN_YEARS = ANNUAL_UPDATE_COUNT

TO_METERS = True
SEA_LEVEL_RISE_RATE = 0.0048
SEA_LEVEL_CONSTANT = True

NUM_CORES = 4

# Barrier3D maximum dune elevation in the original YAML input units (meters).
# A fresh run-specific YAML file is created for every simulation, so the master
# YAML is never changed by CASCADE/Barrier3D.
DMAXEL_M = 3.4

DUNE_REBUILD_HEIGHT = 3.0
REBUILD_ELEV_THRESHOLD = 0.01
REBUILD_DUNE_THRESHOLD = 1.0

BEACH_WIDTH_THRESHOLD_VALUE = 30.0

# =============================================================================
# IMPORTANT SWITCH
# =============================================================================
# False = historical BN only.
# True  = historical BN + threshold-based BN diagnostic/application.
#
# To see whether threshold BN was effective, set this to True.
# To keep historical BN only, set this to False.
# =============================================================================

APPLY_AUTOMATIC_BN_THRESHOLD = False

FIXED_WAVE_PERIOD = 7
FIXED_WAVE_ASYMMETRY = 0.8
FIXED_WAVE_ANGLE_HIGH_FRACTION = 0.1
WAVE_HEIGHTS_TO_TEST = [2]

DOMAIN_TICK_STEP = 5
FLIP_SIGN_MODEL = True

# Use linear-regression rate for a like-for-like comparison with CMP LRR.
# Allowed values: "lrr" or "endpoint".
CASCADE_RATE_METHOD = "lrr"

# Offset CSV column positions.
# Column 0 = Domain_ID; column 1 = 1992 shoreline/dune offset in meters.
DOMAIN_ID_COLUMN_INDEX = 0
OFFSET_COLUMN_INDEX = 1
FIRST_OFFSET_DOMAIN_ID = 9
LAST_OFFSET_DOMAIN_ID = 119

DOMAIN_LENGTH_M = 500.0

MAKE_YEARLY_BN_GIF = True
GIF_DURATION_SECONDS = 4

print("Simulation Configuration:")
print(
    f"  Start Year: {START_YEAR} | End Year: {END_YEAR} | "
    f"annual updates: {ANNUAL_UPDATE_COUNT} | saved states: {TIME_STEP_COUNT}"
)
print(f"  Storm file: {STORM_FILE}")
print(f"  WAVE_HEIGHTS: {WAVE_HEIGHTS_TO_TEST}")
print(f"  DUNE_REBUILD_HEIGHT: {DUNE_REBUILD_HEIGHT} m")
print(f"  REBUILD_ELEV_THRESHOLD: {REBUILD_ELEV_THRESHOLD} m")
print(f"  REBUILD_DUNE_THRESHOLD: {REBUILD_DUNE_THRESHOLD} m")
print(f"  BEACH_WIDTH_THRESHOLD_VALUE: {BEACH_WIDTH_THRESHOLD_VALUE} m")
print(f"  APPLY_AUTOMATIC_BN_THRESHOLD: {APPLY_AUTOMATIC_BN_THRESHOLD}")
print(f"  FLIP_SIGN_MODEL={FLIP_SIGN_MODEL}")
print(f"  CASCADE_RATE_METHOD={CASCADE_RATE_METHOD}")
print("=" * 80 + "\n")


# =============================================================================
# SECTION 4: MANUAL HISTORICAL NOURISHMENT ARRAYS
# =============================================================================

BN_YEARS_BY_DOMAIN = {
    111: [2003],
    112: [2003],
    113: [2002, 2003, 2009],
    114: [1992, 2002, 2003, 2009],
    115: [1992, 2002, 2003, 2004],
    116: [1992, 1993, 2001, 2002, 2003, 2004, 2008, 2009],
    117: [1992, 1993, 1995, 2001, 2002, 2003, 2004, 2008, 2009],
    118: [1992, 2001, 2003, 2004, 2008, 2009],
    119: [2001, 2004, 2008, 2009],
}

BN_VOLUME_M3_BY_DOMAIN = {
    111: [102954.25],
    112: [102954.25],

    113: [96763.6, 102954.25, 236628.75],

    114: [431200.0, 96763.6, 102954.25, 236628.75],

    115: [431200.0, 96763.6, 102954.25, 74972.4],

    116: [
        49146.6667,
        173294.0,
        102741.25,
        96763.6,
        102954.25,
        74972.4,
        316731.5,
        236628.75,
    ],

    117: [
        49146.6667,
        173294.0,
        52184.8,
        102741.25,
        96763.6,
        102954.25,
        74972.4,
        316731.5,
        236628.75,
    ],

    118: [
        49146.6667,
        102741.25,
        102954.25,
        74972.4,
        187431.0,
        316731.5,
    ],

    119: [
        102741.25,
        74972.4,
        187431.0,
        316731.5,
    ],
}


# =============================================================================
# SECTION 5: LOAD INPUT DATA
# =============================================================================

print("Loading input data...")

try:
    # -------------------------------------------------------------------------
    # Load the two-column shoreline/dune-offset CSV using explicit column
    # positions:
    #     column 0 = Domain_ID
    #     column 1 = offset in meters
    #
    # Sorting by Domain_ID is critical. It prevents the CSV row order from
    # changing the alongshore placement of domains in CASCADE.
    # -------------------------------------------------------------------------
    offset_df = pd.read_csv(DUNE_OFFSET_FILE)

    if offset_df.shape[1] != 2:
        raise ValueError(
            "The shoreline-offset CSV must contain exactly two columns:\n"
            "  column 0: Domain_ID\n"
            "  column 1: offset in meters\n"
            f"Found {offset_df.shape[1]} columns: {list(offset_df.columns)}\n"
            f"File: {DUNE_OFFSET_FILE}"
        )

    original_offset_headers = list(offset_df.columns)

    # Rename strictly by position; the second header may be Offset, offset_m,
    # 1992, or another valid label.
    offset_df = offset_df.iloc[:, [
        DOMAIN_ID_COLUMN_INDEX,
        OFFSET_COLUMN_INDEX,
    ]].copy()
    offset_df.columns = ["Domain_ID", "Offset_m"]

    offset_df["Domain_ID"] = pd.to_numeric(
        offset_df["Domain_ID"],
        errors="coerce",
    )
    offset_df["Offset_m"] = pd.to_numeric(
        offset_df["Offset_m"],
        errors="coerce",
    )

    invalid_rows = offset_df[
        offset_df["Domain_ID"].isna()
        | offset_df["Offset_m"].isna()
    ]

    if not invalid_rows.empty:
        raise ValueError(
            "The shoreline-offset CSV contains missing or nonnumeric values:\n"
            + invalid_rows.to_string(index=False)
        )

    # Domain IDs must be whole numbers.
    noninteger_ids = (
        offset_df["Domain_ID"]
        != np.floor(offset_df["Domain_ID"])
    )

    if noninteger_ids.any():
        raise ValueError(
            "Domain_ID values must be integers. Invalid rows:\n"
            + offset_df.loc[noninteger_ids].to_string(index=False)
        )

    offset_df["Domain_ID"] = offset_df["Domain_ID"].astype(int)

    duplicate_ids = offset_df.loc[
        offset_df["Domain_ID"].duplicated(keep=False),
        "Domain_ID",
    ].tolist()

    if duplicate_ids:
        raise ValueError(
            "Duplicate Domain_ID values found in the shoreline-offset CSV: "
            f"{sorted(set(duplicate_ids))}"
        )

    # This run uses D9-D119 only. A legacy D120 row is removed by ID, not by
    # blindly dropping the final CSV row.
    if 120 in offset_df["Domain_ID"].values:
        print("Offset file includes legacy domain D120; removing D120 by ID.")
        offset_df = offset_df.loc[
            offset_df["Domain_ID"] != 120
        ].copy()

    # Critical correction: sort the two-column data by Domain_ID before
    # creating the CASCADE offset vector.
    offset_df = (
        offset_df
        .sort_values("Domain_ID")
        .reset_index(drop=True)
    )

    expected_domain_ids = np.arange(
        FIRST_OFFSET_DOMAIN_ID,
        LAST_OFFSET_DOMAIN_ID + 1,
        dtype=int,
    )
    actual_domain_ids = offset_df["Domain_ID"].to_numpy(dtype=int)

    if not np.array_equal(actual_domain_ids, expected_domain_ids):
        missing_ids = sorted(
            set(expected_domain_ids) - set(actual_domain_ids)
        )
        unexpected_ids = sorted(
            set(actual_domain_ids) - set(expected_domain_ids)
        )

        raise ValueError(
            "The shoreline-offset Domain_ID mapping does not match the current "
            "CASCADE layout.\n"
            f"Expected exactly D{FIRST_OFFSET_DOMAIN_ID}-"
            f"D{LAST_OFFSET_DOMAIN_ID} ({len(expected_domain_ids)} rows).\n"
            f"Found {len(actual_domain_ids)} rows.\n"
            f"Missing IDs: {missing_ids}\n"
            f"Unexpected IDs: {unexpected_ids}"
        )

    # Column index 1 is the correct offset column.
    dune_offset_m = offset_df["Offset_m"].to_numpy(dtype=float)

    if dune_offset_m.size != TOTAL_DOMAINS:
        raise ValueError(
            f"Expected {TOTAL_DOMAINS} shoreline-offset values "
            f"({NUM_BUFFER_DOMAINS_LEFT} left buffers + "
            f"{NUM_REAL_DOMAINS} real domains), but found "
            f"{dune_offset_m.size}."
        )

    # CASCADE/Barrier3D shoreline coordinates are in decameters.
    dune_offset_dam = dune_offset_m / 10.0

    print("\nValidated shoreline-offset input:")
    print(f"  File:             {DUNE_OFFSET_FILE}")
    print(f"  Original headers: {original_offset_headers}")
    print(
        f"  Domain IDs:       D{actual_domain_ids[0]}-"
        f"D{actual_domain_ids[-1]}"
    )
    print(f"  Offset column:    index {OFFSET_COLUMN_INDEX} (meters)")
    print(f"  Number of rows:   {dune_offset_dam.size}")
    print(
        f"  Offset range:     {dune_offset_m.min():.3f} to "
        f"{dune_offset_m.max():.3f} m"
    )

    # Print the real-domain mapping so a shifted or reversed vector is visible
    # before the simulation begins.
    print("\nReal-domain offset mapping:")
    previous_offset_m = None

    for domain_id in range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1):
        cascade_index = START_REAL_INDEX + (
            domain_id - FIRST_FILE_NUMBER
        )
        offset_m = dune_offset_m[cascade_index]

        if previous_offset_m is None:
            jump_text = "first real domain"
        else:
            jump_text = (
                f"jump from previous = "
                f"{offset_m - previous_offset_m:+.3f} m"
            )

        print(
            f"  D{domain_id} -> CASCADE index {cascade_index:3d} -> "
            f"{offset_m:10.3f} m ({dune_offset_dam[cascade_index]:9.3f} dam) | "
            f"{jump_text}"
        )

        previous_offset_m = offset_m

except FileNotFoundError as error:
    print(f"CRITICAL ERROR: Missing data file: {error.filename}")
    sys.exit(1)
except Exception as error:
    print(f"CRITICAL ERROR loading shoreline-offset data: {error}")
    raise

# Background shoreline-change rates in m/yr.
# CASCADE uses a negative value for erosion.
BACKGROUND_EROSION_RATES = [0.0] * TOTAL_DOMAINS

DOMAIN_80_INDEX = START_REAL_INDEX
DOMAIN_119_INDEX = END_REAL_INDEX - 1

# Domain 80: 1 m/yr erosion.
BACKGROUND_EROSION_RATES[DOMAIN_80_INDEX] = 0.0

# Domain 119: no background erosion for now.
BACKGROUND_EROSION_RATES[DOMAIN_119_INDEX] = 0.0

print("\nBackground erosion configuration:")
print(
    f"  Domain 80 | CASCADE index {DOMAIN_80_INDEX}: "
    f"{BACKGROUND_EROSION_RATES[DOMAIN_80_INDEX]:.1f} m/yr"
)
print(
    f"  Domain 119 | CASCADE index {DOMAIN_119_INDEX}: "
    f"{BACKGROUND_EROSION_RATES[DOMAIN_119_INDEX]:.1f} m/yr"
)

# BACKGROUND_EROSION_RATES = [ #buffers
#                            0,0,0,0,0,
#                            0,0,0,0,0,
#                            0,0,0,0,0,
#
#                            0,0,0,0,0,
#                            -0,0,0,0,0,
#                            -0,-0,0,0,0,
#                            0,0,-0,0,-0,
#                            0,0,0,0,0,
#                            0,0,0,0,0,
#                            0,0,0,0,0,
#                            0,0,0,0,-9,-18,
#
#                            #buffers:
#                            0,0,0,0,0,
#                            0,0,0,0,0,
#                            0,0,0,0,0]


# =============================================================================
# ROADWAY MANAGEMENT + HISTORICAL ROAD RELOCATION
# =============================================================================

# Original 1992 roadway setback for real Pea Island domains 80-119 [m].
# Buffer domains remain unmanaged and retain a zero placeholder setback.
ROAD_SETBACK_REAL_M = [
    450, 360, 180, 94, 85, 125, 100, 85, 159, 143,
    168, 159, 208, 211, 250, 270, 240, 185, 107, 119,
    70, 33, 12, 10, 33, 55, 52, 60, 109, 113,
    104, 98, 86, 48, 65, 66, 55, 44, 62, 110,
]

if len(ROAD_SETBACK_REAL_M) != NUM_REAL_DOMAINS:
    raise ValueError(
        f"ROAD_SETBACK_REAL_M has {len(ROAD_SETBACK_REAL_M)} values; "
        f"expected {NUM_REAL_DOMAINS}."
    )

road_setbacks_full = np.zeros(TOTAL_DOMAINS, dtype=float)
road_setbacks_full[START_REAL_INDEX:END_REAL_INDEX] = ROAD_SETBACK_REAL_M

# Roadway management is active only for real domains 80-119.
ROADWAY_MANAGEMENT_ON = [False] * TOTAL_DOMAINS
ROADWAY_MANAGEMENT_ON[START_REAL_INDEX:END_REAL_INDEX] = [True] * NUM_REAL_DOMAINS

BRIDGE_DOMAINS = [119]

for _bridge_domain in BRIDGE_DOMAINS:
    _bridge_index = START_REAL_INDEX + (
        _bridge_domain - FIRST_FILE_NUMBER
    )
    ROADWAY_MANAGEMENT_ON[_bridge_index] = False
# -----------------------------------------------------------------------------
# HISTORICAL ROAD-RELOCATION EVENTS — ROADWAY MANAGER TRIGGER
# -----------------------------------------------------------------------------
# Each event supplies the calendar year, affected domains, and the target
# POST-RELOCATION setback. Immediately before the normal cascade.update(), the
# run script calls RoadwayManager.request_relocation(target_setback) for each
# affected domain. RoadwayManager consumes that request during its next update
# and performs the relocation, elevation, dune, bulldozing, and time-series
# updates. Every queued request is an observed hindcast boundary condition, so
# relocate_now automatically bypasses event-year width and drowning rejection.
#
# cascade.py remains unchanged, and the run script never directly assigns
# roadway._road_setback or its time-series values.
#
# For 1996, the target setbacks below equal the supplied 1992 initial setback
# plus the documented landward movement distance:
#   D100: 70 + 100 = 170 m
#   D101: 33 + 200 = 233 m
#   D102: 12 +  90 = 102 m
#   D103: 10 +  90 = 100 m
#   D104: 33 +  90 = 123 m
#   D105: 55 +  90 = 145 m
#   D106: 52 +  90 = 142 m
#   D107: 60 +  90 = 150 m
#
# Do not activate an event until a target setback or a landward movement is
# supplied for every affected domain. The helper below raises an error rather
# than silently placing the historical road at an assumed location.
HISTORICAL_ROAD_EVENTS = [
    dict(
        year=1996,
        domain_start=100,
        domain_end=107,
        type="relocate",
        note="NC-12 relocated landward in 1996, Pea Island NWR",
        post_relocation_setback_m={
            100: 170.0,
            101: 233.0,
            102: 102.0,
            103: 100.0,
            104: 123.0,
            105: 145.0,
            106: 142.0,
            107: 150.0,
        },
    ),
    # dict(
    #     year=1998,
    #     domain_start=82,
    #     domain_end=88,
    #     type="relocate",
    #     note="NC-12 S-curves relocation in 1998",
    #     # Fill every affected domain before activating this event:
    #     post_relocation_setback_m={
    #         # 82: ...,
    #         # 83: ...,
    #         # 84: ...,
    #         # 85: ...,
    #         # 86: ...,
    #         # 87: ...,
    #         # 88: ...,
    #     },
    # ),
]

print("\nRoadway management configuration:")
print("  Buffers: OFF")
print(f"  Real domains {FIRST_FILE_NUMBER}-{LAST_FILE_NUMBER}: ON")
for _event in HISTORICAL_ROAD_EVENTS:
    _mapping = _event.get("post_relocation_setback_m", {})
    if _mapping:
        print(
            f"  Relocation year {_event['year']}: prescribed post-relocation "
            f"setbacks={_mapping} m"
        )
    else:
        print(
            f"  Relocation year {_event['year']}: domains "
            f"{_event['domain_start']}-{_event['domain_end']}, but no target "
            "setbacks are defined"
        )
print("=" * 80 + "\n")


# =============================================================================
# SECTION 6: ELEVATION + DUNE FILE LISTS
# =============================================================================

print("Generating elevation + dune file paths...")

ELEVATION_FILE_PATHS = []
DUNE_FILE_PATHS = []

LEFT_BUFFER_SUBDIR = os.path.join(
    PEA_DATA_BASE,
    "Buffer",
    "last_one_repeat",
    "leftside",
)
LEFT_BUFFER_DUNE_FILE = os.path.join(
    LEFT_BUFFER_SUBDIR,
    "buffer_80_actual_dune_2row.npy",
)
LEFT_BUFFER_TOPO_FILE = os.path.join(
    LEFT_BUFFER_SUBDIR,
    "buffer_80_actual_topo.npy",
)

# 71 left-side buffer domains corresponding to logical IDs 9-79.
for _ in range(NUM_BUFFER_DOMAINS_LEFT):
    DUNE_FILE_PATHS.append(LEFT_BUFFER_DUNE_FILE)
    ELEVATION_FILE_PATHS.append(LEFT_BUFFER_TOPO_FILE)

# Real domains 80-119 only.
for file_num in range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1):
    DUNE_FILE_PATHS.append(
        os.path.join(
            PEA_DATA_BASE,
            "1996",
            "dune_topo_new",
            "2011_Roya_1996_cross_ROAD_2row_alongshore_corrected",
            "dunes",
            f"domain_{file_num}_dune_2011.npy",
        )
    )
    ELEVATION_FILE_PATHS.append(
        os.path.join(
            PEA_DATA_BASE,
            "1996",
            "dune_topo_new",
            "2011_Roya_1996_cross_ROAD_2row_alongshore_corrected",
            "topography",
            f"domain_{file_num}_topography_2011.npy",
        )
    )

if len(ELEVATION_FILE_PATHS) != TOTAL_DOMAINS:
    raise ValueError(
        f"Generated {len(ELEVATION_FILE_PATHS)} elevation files; "
        f"expected {TOTAL_DOMAINS}."
    )

if len(DUNE_FILE_PATHS) != TOTAL_DOMAINS:
    raise ValueError(
        f"Generated {len(DUNE_FILE_PATHS)} dune files; expected {TOTAL_DOMAINS}."
    )

missing_files = [
    path
    for path in DUNE_FILE_PATHS + ELEVATION_FILE_PATHS
    if not os.path.isfile(path)
]

if missing_files:
    print("\nMissing dune/topography files:")
    for path in sorted(set(missing_files))[:20]:
        print(f"  {path}")
    raise FileNotFoundError("One or more dune/topography files are missing.")

print(f"Generated {len(ELEVATION_FILE_PATHS)} elevation files")
print(f"Generated {len(DUNE_FILE_PATHS)} dune files")
print("=" * 80 + "\n")


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def real_domain_to_padded_index(real_domain):
    return START_REAL_INDEX + (real_domain - FIRST_FILE_NUMBER)


def padded_index_to_real_domain(padded_index):
    if START_REAL_INDEX <= padded_index < END_REAL_INDEX:
        return FIRST_FILE_NUMBER + (padded_index - START_REAL_INDEX)
    return None


def real_domain_to_x(real_domain_id):
    return START_REAL_INDEX + (real_domain_id - FIRST_FILE_NUMBER)


def get_x_s_TS(b3d):
    if hasattr(b3d, "x_s_TS"):
        return np.asarray(b3d.x_s_TS, dtype=float)
    if hasattr(b3d, "_x_s_TS"):
        return np.asarray(b3d._x_s_TS, dtype=float)
    raise AttributeError("No shoreline time series found on Barrier3D object.")


def build_shoreline_matrix(cascade, to_meters=True):
    b3d_list = cascade.barrier3d
    ndom = len(b3d_list)
    nt = len(get_x_s_TS(b3d_list[0]))

    shoreline = np.zeros((nt, ndom), dtype=float)

    for j in range(ndom):
        shoreline[:, j] = get_x_s_TS(b3d_list[j])

    if to_meters:
        shoreline *= 10.0

    return shoreline


def build_relative_shoreline_change_matrix(cascade, to_meters=True, flip_sign=True):
    shoreline_m = build_shoreline_matrix(cascade, to_meters=to_meters)
    shoreline_change_m = shoreline_m - shoreline_m[0, :]

    if flip_sign:
        shoreline_change_m *= -1.0

    return shoreline_change_m


def _resolve_column(df, preferred, aliases, required=True):
    """Return the first matching dataframe column, case-insensitively."""
    available = {str(col).strip().lower(): col for col in df.columns}

    candidates = []
    if preferred is not None:
        candidates.append(preferred)
    candidates.extend(aliases)

    for candidate in candidates:
        key = str(candidate).strip().lower()
        if key in available:
            return available[key]

    if required:
        raise ValueError(
            f"Could not identify required column. Tried: {candidates}. "
            f"Available columns: {list(df.columns)}"
        )

    return None


def load_cmp_domain_rates(
    file_path,
    sheet_name=0,
    domain_col=None,
    rate_col=None,
    std_col=None,
    count_col=None,
    flip_sign=False,
):
    """Load domain-averaged CMP 1992–2024 shoreline-change rates.

    Required fields are a CASCADE real-domain ID (80–119) and a rate in m/yr.
    Optional standard-deviation and transect-count fields are retained when
    present. CSV and Excel inputs are supported.
    """
    if file_path is None or not os.path.isfile(file_path):
        raise FileNotFoundError(
            "CMP shoreline-rate file was not found. Update CMP_FILE near the "
            f"top of the script. Current value: {file_path}"
        )

    extension = os.path.splitext(file_path)[1].lower()

    if extension in {".xlsx", ".xls", ".xlsm"}:
        df = pd.read_excel(file_path, sheet_name=sheet_name)
    elif extension in {".csv", ".txt"}:
        df = pd.read_csv(file_path)
    else:
        raise ValueError(
            f"Unsupported CMP file type '{extension}'. Use CSV or Excel."
        )

    domain_name = _resolve_column(
        df,
        domain_col,
        aliases=[
            "domain", "domain_id", "real_domain", "cascade_domain",
            "CASCADE domain", "CASCADE_domain_ID",
        ],
        required=True,
    )

    rate_name = _resolve_column(
        df,
        rate_col,
        aliases=[
            "mean_LRR_m_yr",
            "mean_LRR_CoastSat_sign_m_yr",
            "CMP_LRR_1992_2024_m_yr",
            "cmp_lrr_1992_2024",
            "CMP 1992–2024 (m/yr)",
            "shoreline_change_rate_m_yr",
            "rate_m_yr",
        ],
        required=True,
    )

    std_name = _resolve_column(
        df,
        std_col,
        aliases=[
            "std_dev_m_yr", "std_m_yr", "CMP_std_m_yr",
            "cmp_std_1992_2024", "standard_deviation_m_yr",
        ],
        required=False,
    )

    count_name = _resolve_column(
        df,
        count_col,
        aliases=[
            "CMP_transect_count", "transect_count", "n_transects", "count",
        ],
        required=False,
    )

    out = pd.DataFrame(
        {
            "domain_id": pd.to_numeric(df[domain_name], errors="coerce"),
            "cmp_rate_m_yr": pd.to_numeric(df[rate_name], errors="coerce"),
        }
    )

    if std_name is not None:
        out["cmp_std_m_yr"] = pd.to_numeric(df[std_name], errors="coerce")
    else:
        out["cmp_std_m_yr"] = np.nan

    if count_name is not None:
        out["cmp_transect_count"] = pd.to_numeric(
            df[count_name], errors="coerce"
        )
    else:
        out["cmp_transect_count"] = np.nan

    out = out.dropna(subset=["domain_id", "cmp_rate_m_yr"]).copy()
    out["domain_id"] = out["domain_id"].astype(int)

    out = out.loc[
        out["domain_id"].between(FIRST_FILE_NUMBER, LAST_FILE_NUMBER)
    ].copy()

    if out.empty:
        raise ValueError(
            f"CMP file contains no valid domains in the required range "
            f"{FIRST_FILE_NUMBER}–{LAST_FILE_NUMBER}."
        )

    # If multiple rows exist per domain, combine them safely.
    aggregation = {
        "cmp_rate_m_yr": "mean",
        "cmp_std_m_yr": "mean",
        "cmp_transect_count": "sum",
    }
    out = out.groupby("domain_id", as_index=False).agg(aggregation)

    if flip_sign:
        out["cmp_rate_m_yr"] *= -1.0

    out = out.sort_values("domain_id").reset_index(drop=True)

    print("\nCMP shoreline-change-rate input:")
    print(f"  File: {file_path}")
    if extension in {".xlsx", ".xls", ".xlsm"}:
        print(f"  Sheet: {sheet_name}")
    print(f"  Domain column: {domain_name}")
    print(f"  Rate column: {rate_name}")
    print(f"  Standard-deviation column: {std_name}")
    print(f"  Transect-count column: {count_name}")
    print(f"  Loaded domains: {len(out)}")
    print(
        f"  Domain range: {out['domain_id'].min()}–{out['domain_id'].max()}"
    )
    print(
        f"  CMP rate range: {out['cmp_rate_m_yr'].min():.3f} to "
        f"{out['cmp_rate_m_yr'].max():.3f} m/yr"
    )

    return out


def calculate_cascade_shoreline_change_rate(
    shoreline_m,
    method="lrr",
    flip_sign=True,
):
    """Calculate one shoreline-change rate per CASCADE domain in m/yr."""
    shoreline_m = np.asarray(shoreline_m, dtype=float)

    if shoreline_m.ndim != 2:
        raise ValueError(
            f"shoreline_m must be 2-D (time, domain); got {shoreline_m.shape}."
        )

    nt = shoreline_m.shape[0]
    if nt < 2:
        raise ValueError("At least two shoreline time steps are required.")

    method = str(method).strip().lower()

    if method == "endpoint":
        elapsed_years = float(nt - 1)
        rates = (shoreline_m[-1, :] - shoreline_m[0, :]) / elapsed_years

    elif method == "lrr":
        years = np.arange(nt, dtype=float)
        rates = np.full(shoreline_m.shape[1], np.nan, dtype=float)

        for domain_index in range(shoreline_m.shape[1]):
            values = shoreline_m[:, domain_index]
            valid = np.isfinite(years) & np.isfinite(values)

            if np.count_nonzero(valid) >= 2:
                rates[domain_index] = np.polyfit(
                    years[valid], values[valid], 1
                )[0]

    else:
        raise ValueError(
            "CASCADE_RATE_METHOD must be either 'lrr' or 'endpoint'."
        )

    if flip_sign:
        rates *= -1.0

    return rates


def calculate_comparison_statistics(observed, modeled):
    """Return calibration statistics for matched finite CMP/CASCADE rates."""
    observed = np.asarray(observed, dtype=float)
    modeled = np.asarray(modeled, dtype=float)

    valid = np.isfinite(observed) & np.isfinite(modeled)
    observed = observed[valid]
    modeled = modeled[valid]

    if len(observed) == 0:
        return dict(n=0, bias=np.nan, mae=np.nan, rmse=np.nan, r2=np.nan)

    residual = modeled - observed
    bias = np.mean(residual)
    mae = np.mean(np.abs(residual))
    rmse = np.sqrt(np.mean(residual ** 2))

    if len(observed) >= 2 and np.std(observed) > 0 and np.std(modeled) > 0:
        r = np.corrcoef(observed, modeled)[0, 1]
        r2 = r ** 2
    else:
        r2 = np.nan

    return dict(n=len(observed), bias=bias, mae=mae, rmse=rmse, r2=r2)


def build_nourishment_arrays_from_manual_inputs():
    nourishment_on_by_year = {}
    nourishment_volume_by_year = {}

    for year in range(START_YEAR, END_YEAR + 1):
        nourishment_on_by_year[year] = np.zeros(TOTAL_DOMAINS)
        nourishment_volume_by_year[year] = [0.0] * TOTAL_DOMAINS

    for real_domain, years in BN_YEARS_BY_DOMAIN.items():

        if real_domain not in BN_VOLUME_M3_BY_DOMAIN:
            raise ValueError(f"Missing nourishment volume list for domain {real_domain}")

        volumes = BN_VOLUME_M3_BY_DOMAIN[real_domain]

        if len(years) != len(volumes):
            raise ValueError(
                f"Domain {real_domain}: years and volumes must have same length."
            )

        padded_index = real_domain_to_padded_index(real_domain)

        if padded_index < 0 or padded_index >= TOTAL_DOMAINS:
            print(f"Skipping domain {real_domain}: outside CASCADE range")
            continue

        for year, total_volume_m3 in zip(years, volumes):

            if year < START_YEAR or year > END_YEAR:
                continue

            volume_m3_per_m = float(total_volume_m3) / DOMAIN_LENGTH_M

            nourishment_on_by_year[year][padded_index] = 1
            nourishment_volume_by_year[year][padded_index] = volume_m3_per_m

    print("\nHistorical nourishment schedule from manual arrays:")
    for year in range(START_YEAR, END_YEAR + 1):
        active = np.where(nourishment_on_by_year[year] == 1)[0]
        if len(active) > 0:
            real_domains = [
                FIRST_FILE_NUMBER + (idx - START_REAL_INDEX)
                for idx in active
                if START_REAL_INDEX <= idx < END_REAL_INDEX
            ]

            total_volume = (
                np.sum(np.asarray(nourishment_volume_by_year[year], dtype=float))
                * DOMAIN_LENGTH_M
            )

            print(
                f"  {year}: domains {real_domains}, "
                f"total volume = {total_volume:,.0f} m3"
            )

    return nourishment_on_by_year, nourishment_volume_by_year


def add_site_annotations(ax):
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()

    yrange = ymax - ymin
    xrange = xmax - xmin

    x_dom_80 = real_domain_to_x(80)
    x_dom_81 = real_domain_to_x(81)
    x_dom_110 = real_domain_to_x(110)

    x_rodanthe = 0.5 * (x_dom_80 + x_dom_81)

    site_line_kwargs = dict(
        color="0.55",
        linestyle="--",
        linewidth=1.0,
        alpha=0.60,
        zorder=0,
    )

    ax.axvline(x_rodanthe, **site_line_kwargs)
    ax.axvline(x_dom_110, **site_line_kwargs)

    ax.text(
        x_rodanthe,
        ymax - 0.035 * yrange,
        "Rodanthe",
        fontsize=8.5,
        color="0.35",
        ha="center",
        va="top",
    )

    ax.text(
        x_dom_110,
        ymax - 0.035 * yrange,
        "Pea Island\nVisitor Center",
        fontsize=8.5,
        color="0.35",
        ha="center",
        va="top",
    )

    # ax.annotate(
    #     "Oregon Inlet\nTerminal Groin\nDomain 120",
    #     xy=(x_dom_120, ymax - 0.20 * yrange),
    #     xytext=(x_dom_120 - 10, ymax - 0.02 * yrange),
    #     fontsize=9,
    #     color="0.22",
    #     ha="center",
    #     va="top",
    #     arrowprops=dict(
    #         arrowstyle="-|>",
    #         color="0.35",
    #         lw=1.2,
    #         shrinkA=4,
    #         shrinkB=4,
    #         connectionstyle="arc3,rad=-0.18",
    #     ),
    #     bbox=dict(
    #         boxstyle="round,pad=0.25",
    #         fc="white",
    #         ec="0.70",
    #         alpha=0.88,
    #     ),
    # )

    # ax.annotate(
    #     "Rodanthe\nDomains 80–81",
    #     xy=(x_rodanthe, ymin + 0.35 * yrange),
    #     xytext=(x_rodanthe + 8, ymin + 0.15 * yrange),
    #     fontsize=9,
    #     color="0.22",
    #     ha="center",
    #     va="center",
    #     arrowprops=dict(
    #         arrowstyle="-|>",
    #         color="0.35",
    #         lw=1.2,
    #         shrinkA=4,
    #         shrinkB=4,
    #         connectionstyle="arc3,rad=0.22",
    #     ),
    #     bbox=dict(
    #         boxstyle="round,pad=0.25",
    #         fc="white",
    #         ec="0.70",
    #         alpha=0.88,
    #     ),
    # )
    #
    # ax.annotate(
    #     "Pea Island\nVisitor Center\nDomain 110",
    #     xy=(x_dom_110, ymax - 0.35 * yrange),
    #     xytext=(x_dom_110 - 8, ymax - 0.13 * yrange),
    #     fontsize=9,
    #     color="0.22",
    #     ha="center",
    #     va="top",
    #     arrowprops=dict(
    #         arrowstyle="-|>",
    #         color="0.35",
    #         lw=1.2,
    #         shrinkA=4,
    #         shrinkB=4,
    #         connectionstyle="arc3,rad=0.18",
    #     ),
    #     bbox=dict(
    #         boxstyle="round,pad=0.25",
    #         fc="white",
    #         ec="0.70",
    #         alpha=0.88,
    #     ),
    # )

    x_arrow = xmax - 0.035 * xrange

    ax.annotate(
        "Accretion  +",
        xy=(x_arrow, ymax - 0.20 * yrange),
        xytext=(x_arrow, ymax - 0.34 * yrange),
        fontsize=9,
        color="0.25",
        ha="center",
        va="center",
        arrowprops=dict(
            arrowstyle="-|>",
            color="0.35",
            lw=1.2,
        ),
    )

    ax.annotate(
        "Erosion  −",
        xy=(x_arrow, ymin + 0.20 * yrange),
        xytext=(x_arrow, ymin + 0.34 * yrange),
        fontsize=9,
        color="0.25",
        ha="center",
        va="center",
        arrowprops=dict(
            arrowstyle="-|>",
            color="0.35",
            lw=1.2,
        ),
    )


def add_real_domain_top_axis(ax):
    top_ax = ax.secondary_xaxis("top")
    top_ax.set_xlabel(
        "Real Pea Island domain ID",
        fontsize=10,
        fontweight="bold",
        labelpad=8,
    )

    top_positions = []
    top_labels = []

    for dom_id in range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP):
        j = START_REAL_INDEX + (dom_id - FIRST_FILE_NUMBER)
        top_positions.append(j)
        top_labels.append(str(dom_id))

    top_ax.set_xticks(top_positions)
    top_ax.set_xticklabels(top_labels, fontsize=9)

    return top_ax


def get_roadway_manager_list(cascade):
    """Return the actual RoadwayManager instance list."""
    if hasattr(cascade, "_roadways"):
        roadways = cascade._roadways
    elif hasattr(cascade, "roadways"):
        roadways = cascade.roadways
    else:
        raise AttributeError(
            "Could not find RoadwayManager objects on the Cascade instance. "
            "Expected cascade._roadways or cascade.roadways."
        )

    if roadways is None:
        raise RuntimeError("CASCADE roadway-manager list is None.")

    return roadways



def queue_historical_road_relocations(
    cascade,
    current_year,
    historical_road_events,
    alongshore_section_count,
    roadway_management_on,
):
    """Queue prescribed relocations on RoadwayManager objects for one update.

    The RoadwayManager request is consumed and automatically cleared during the
    following ``cascade.update()``. Therefore, ``cascade.py`` does not need to
    accept or forward any new arguments.

    Event records may define either:

    * ``post_relocation_setback_m``: absolute target setback from the modeled
      dune/interior edge; or
    * ``landward_movement_m``: distance added to the road's current modeled
      setback at the start of the event year.

    Returns
    -------
    list of dict
        Diagnostics collected before the model update.
    """
    active_event_records = []

    if not historical_road_events:
        return active_event_records

    roadways = get_roadway_manager_list(cascade)

    for event in historical_road_events:
        if current_year != int(event["year"]):
            continue
        if event.get("type", "relocate") != "relocate":
            continue

        domain_start = int(event["domain_start"])
        domain_end = int(event["domain_end"])
        post_setbacks = event.get("post_relocation_setback_m", {})
        movement_distances = event.get("landward_movement_m", {})
        prescribed_road_elevations = event.get("road_elevation_m", {})

        print(
            f"\nApplying HISTORICAL ROAD RELOCATION in {current_year}: "
            f"{event.get('note', '')}"
        )

        for real_domain in range(domain_start, domain_end + 1):
            iB3D = real_domain_to_padded_index(real_domain)

            if not (0 <= iB3D < alongshore_section_count):
                raise IndexError(
                    f"Historical relocation domain {real_domain} maps to "
                    f"invalid CASCADE index {iB3D}."
                )

            if not roadway_management_on[iB3D]:
                raise RuntimeError(
                    f"Historical relocation requested for unmanaged domain "
                    f"{real_domain} (CASCADE index {iB3D})."
                )

            roadway = roadways[iB3D]
            if roadway is None or not hasattr(roadway, "request_relocation"):
                raise RuntimeError(
                    "The imported RoadwayManager does not provide "
                    "request_relocation(). Replace cascade/roadway_manager.py "
                    "with the no-cascade-change hindcast version. "
                    f"Affected domain: {real_domain}."
                )

            if hasattr(roadway, "road_setback"):
                old_setback = float(roadway.road_setback)
            else:
                old_setback = float(roadway._road_setback)

            if real_domain in post_setbacks:
                target_setback = float(post_setbacks[real_domain])
                target_source = "prescribed post-relocation setback"
            elif real_domain in movement_distances:
                movement_m = float(movement_distances[real_domain])
                target_setback = old_setback + movement_m
                target_source = "current modeled setback + prescribed movement"
            else:
                raise ValueError(
                    f"Historical relocation in {current_year} has no target "
                    f"setback or movement distance for domain {real_domain}."
                )

            if not np.isfinite(target_setback) or target_setback < 0.0:
                raise ValueError(
                    f"Invalid relocation target for domain {real_domain} in "
                    f"{current_year}: {target_setback}"
                )

            # Queue the observed hindcast event. Do not directly change the
            # current road position; RoadwayManager enforces it during update().
            prescribed_road_elevation = prescribed_road_elevations.get(real_domain)
            roadway.request_relocation(
                target_setback,
                road_elevation=prescribed_road_elevation,
            )

            active_event_records.append(
                dict(
                    model_year=current_year,
                    padded_index=iB3D,
                    real_domain=real_domain,
                    road_setback_before_m=old_setback,
                    requested_relocation_setback_m=target_setback,
                    prescribed_road_elevation_m=prescribed_road_elevation,
                    setback_source=target_source,
                    note=event.get("note", ""),
                )
            )

            print(
                f"  Domain {real_domain} (index {iB3D}): queued enforced "
                f"hindcast relocation from current setback {old_setback:.1f} m "
                f"to target {target_setback:.1f} m [{target_source}]"
            )

    return active_event_records



# =============================================================================
# YEARLY PLOTTING FUNCTIONS
# =============================================================================

def plot_yearly_relative_shoreline_and_bn(
    cascades_by_label,
    hist_nourish_volume_by_year,
    output_dir,
    filename_prefix,
    flip_sign=True,
    make_gif=True,
    gif_duration_seconds=4,
):
    yearly_dir = os.path.join(output_dir, filename_prefix)
    os.makedirs(yearly_dir, exist_ok=True)

    domain_numbers = np.arange(TOTAL_DOMAINS)

    shoreline_change_by_label = {}
    max_nt = None

    for label, cascade in cascades_by_label.items():
        shoreline_change_m = build_relative_shoreline_change_matrix(
            cascade,
            to_meters=True,
            flip_sign=flip_sign,
        )

        shoreline_change_by_label[label] = shoreline_change_m

        if max_nt is None:
            max_nt = shoreline_change_m.shape[0]
        else:
            max_nt = min(max_nt, shoreline_change_m.shape[0])

    if max_nt is None:
        print("No cascades available for yearly plotting.")
        return []

    all_vals = []

    for shoreline_change_m in shoreline_change_by_label.values():
        vals = shoreline_change_m[:max_nt, :]
        vals = vals[np.isfinite(vals)]
        all_vals.extend(vals)

    if len(all_vals) > 0:
        y_min = np.nanmin(all_vals)
        y_max = np.nanmax(all_vals)

        if y_max > y_min:
            y_pad = 0.15 * (y_max - y_min)
        else:
            y_pad = 5.0

        y_min -= y_pad
        y_max += y_pad
    else:
        y_min, y_max = -10.0, 10.0

    all_bn_total_m3 = []

    for year in range(START_YEAR, END_YEAR + 1):
        if year in hist_nourish_volume_by_year:
            vols_total_m3 = (
                np.asarray(hist_nourish_volume_by_year[year], dtype=float)
                * DOMAIN_LENGTH_M
            )
            all_bn_total_m3.extend(vols_total_m3[vols_total_m3 > 0])

    if len(all_bn_total_m3) > 0:
        max_bn_volume = max(all_bn_total_m3)
    else:
        max_bn_volume = 1.0

    preferred_bn_label = "Roadway management + historical BN"

    if preferred_bn_label in shoreline_change_by_label:
        bn_curve_label = preferred_bn_label
    else:
        bn_curve_label = list(shoreline_change_by_label.keys())[-1]

    png_files = []

    for t in range(max_nt):
        model_year = START_YEAR + t

        if model_year in hist_nourish_volume_by_year:
            bn_total_m3 = (
                np.asarray(hist_nourish_volume_by_year[model_year], dtype=float)
                * DOMAIN_LENGTH_M
            )
        else:
            bn_total_m3 = np.zeros(TOTAL_DOMAINS)

        active_bn = np.where(bn_total_m3 > 0)[0]

        fig, (ax, ax_bn) = plt.subplots(
            2,
            1,
            figsize=(22, 9),
            sharex=True,
            gridspec_kw={"height_ratios": [3.2, 1.15]},
            constrained_layout=True,
        )

        for label, shoreline_change_m in shoreline_change_by_label.items():
            ax.plot(
                domain_numbers,
                shoreline_change_m[t, :],
                linewidth=2.2,
                label=label,
            )

        ax.axvline(START_REAL_INDEX, linestyle="--", linewidth=1.0, color="0.45", alpha=0.70)
        ax.axvline(END_REAL_INDEX - 1, linestyle="--", linewidth=1.0, color="0.45", alpha=0.70)
        ax.axhline(0.0, linestyle="--", linewidth=1.2, color="0.25", alpha=0.75)

        if len(active_bn) > 0:
            x_bn = []
            y_bn = []

            for idx in active_bn:
                idx = int(idx)

                x_bn.append(idx)
                y_bn.append(shoreline_change_by_label[bn_curve_label][t, idx])

                ax.axvline(
                    idx,
                    linestyle=":",
                    linewidth=1.0,
                    color="red",
                    alpha=0.45,
                    zorder=2,
                )

            ax.scatter(
                x_bn,
                y_bn,
                s=35,
                marker="o",
                color="red",
                edgecolor="black",
                linewidth=0.6,
                zorder=10,
                label="Historical BN domain",
            )

        ax.set_xlim(0, TOTAL_DOMAINS - 1)
        ax.set_ylim(y_min, y_max)

        if flip_sign:
            ax.set_ylabel(
                f"Relative shoreline change from {START_YEAR} (m)\n"
                "Accretion = positive, erosion = negative",
                fontsize=11,
                fontweight="bold",
            )
        else:
            ax.set_ylabel(
                f"Raw shoreline-position change from {START_YEAR} (m)",
                fontsize=11,
                fontweight="bold",
            )

        ax.set_title(
            f"Pea Island relative shoreline change with historical BN | {model_year}",
            fontsize=14,
            fontweight="bold",
            pad=12,
        )

        ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.35)
        ax.legend(loc="best", fontsize=9, frameon=True, framealpha=0.92)

        add_site_annotations(ax)
        add_real_domain_top_axis(ax)

        ax_bn.bar(
            domain_numbers,
            bn_total_m3,
            width=0.8,
            color="tab:blue",
            alpha=0.85,
        )

        total_bn_this_year = np.sum(bn_total_m3)

        if len(active_bn) > 0:
            bn_title = (
                f"Historical beach nourishment in {model_year}: YES | "
                f"Total volume = {total_bn_this_year:,.0f} m³"
            )
        else:
            bn_title = f"Historical beach nourishment in {model_year}: NO"

        ax_bn.set_title(bn_title, fontsize=12, fontweight="bold")
        ax_bn.set_ylabel("Historical BN\nvolume\n(m³/domain)", fontsize=10, fontweight="bold")
        ax_bn.set_xlabel("CASCADE domain index including buffers", fontsize=11, fontweight="bold")
        ax_bn.set_ylim(0, max_bn_volume * 1.25)
        ax_bn.grid(True, axis="y", linestyle=":", linewidth=0.8, alpha=0.35)

        xticks = np.arange(0, TOTAL_DOMAINS, DOMAIN_TICK_STEP)
        ax_bn.set_xticks(xticks)
        ax_bn.set_xticklabels([str(i) for i in xticks], rotation=45, ha="right", fontsize=9)

        for idx in active_bn:
            idx = int(idx)

            real_domain = padded_index_to_real_domain(idx)
            volume = bn_total_m3[idx]

            ax_bn.text(
                idx,
                volume + max_bn_volume * 0.025,
                f"D{real_domain}\n{volume/1000:.0f}k",
                ha="center",
                va="bottom",
                fontsize=8,
                color="0.25",
            )

        fig_out = os.path.join(yearly_dir, f"{filename_prefix}_{model_year}.png")
        fig.savefig(fig_out, dpi=220, bbox_inches="tight")
        plt.close(fig)

        png_files.append(fig_out)

    print("\nSaved yearly relative shoreline + historical BN plots to:")
    print(yearly_dir)

    if make_gif:
        try:
            import imageio.v2 as imageio

            images = [imageio.imread(f) for f in png_files]
            gif_out = os.path.join(output_dir, f"{filename_prefix}.gif")

            imageio.mimsave(
                gif_out,
                images,
                duration=gif_duration_seconds,
                loop=0,
            )

            print("\nSaved yearly relative shoreline + historical BN GIF:")
            print(gif_out)

        except ImportError:
            print("\nimageio is not installed. GIF was not created.")
            print("Install it with:")
            print("pip install imageio")

    return png_files


def plot_yearly_absolute_shoreline_and_bn(
    cascades_by_label,
    hist_nourish_volume_by_year,
    output_dir,
    filename_prefix,
    make_gif=True,
    gif_duration_seconds=4,
):
    yearly_dir = os.path.join(output_dir, filename_prefix)
    os.makedirs(yearly_dir, exist_ok=True)

    domain_numbers = np.arange(TOTAL_DOMAINS)

    shoreline_position_by_label = {}
    max_nt = None

    for label, cascade in cascades_by_label.items():
        shoreline_m = build_shoreline_matrix(
            cascade,
            to_meters=True,
        )

        shoreline_position_by_label[label] = shoreline_m

        if max_nt is None:
            max_nt = shoreline_m.shape[0]
        else:
            max_nt = min(max_nt, shoreline_m.shape[0])

    if max_nt is None:
        print("No cascades available for absolute shoreline plotting.")
        return []

    all_vals = []

    for shoreline_m in shoreline_position_by_label.values():
        vals = shoreline_m[:max_nt, :]
        vals = vals[np.isfinite(vals)]
        all_vals.extend(vals)

    if len(all_vals) > 0:
        y_min = np.nanmin(all_vals)
        y_max = np.nanmax(all_vals)

        if y_max > y_min:
            y_pad = 0.08 * (y_max - y_min)
        else:
            y_pad = 50.0

        y_min -= y_pad
        y_max += y_pad
    else:
        y_min, y_max = -6000.0, -5000.0

    all_bn_total_m3 = []

    for year in range(START_YEAR, END_YEAR + 1):
        if year in hist_nourish_volume_by_year:
            vols_total_m3 = (
                np.asarray(hist_nourish_volume_by_year[year], dtype=float)
                * DOMAIN_LENGTH_M
            )
            all_bn_total_m3.extend(vols_total_m3[vols_total_m3 > 0])

    if len(all_bn_total_m3) > 0:
        max_bn_volume = max(all_bn_total_m3)
    else:
        max_bn_volume = 1.0

    preferred_bn_label = "Roadway management + historical BN"

    if preferred_bn_label in shoreline_position_by_label:
        bn_curve_label = preferred_bn_label
    else:
        bn_curve_label = list(shoreline_position_by_label.keys())[-1]

    png_files = []

    for t in range(max_nt):
        model_year = START_YEAR + t

        if model_year in hist_nourish_volume_by_year:
            bn_total_m3 = (
                np.asarray(hist_nourish_volume_by_year[model_year], dtype=float)
                * DOMAIN_LENGTH_M
            )
        else:
            bn_total_m3 = np.zeros(TOTAL_DOMAINS)

        active_bn = np.where(bn_total_m3 > 0)[0]

        fig, (ax, ax_bn) = plt.subplots(
            2,
            1,
            figsize=(22, 9),
            sharex=True,
            gridspec_kw={"height_ratios": [3.2, 1.15]},
            constrained_layout=True,
        )

        for label, shoreline_m in shoreline_position_by_label.items():
            ax.plot(
                domain_numbers,
                shoreline_m[t, :],
                linewidth=2.2,
                label=label,
            )

        ax.axvline(START_REAL_INDEX, linestyle="--", linewidth=1.0, color="0.45", alpha=0.70)
        ax.axvline(END_REAL_INDEX - 1, linestyle="--", linewidth=1.0, color="0.45", alpha=0.70)

        if len(active_bn) > 0:
            x_bn = []
            y_bn = []

            for idx in active_bn:
                idx = int(idx)

                x_bn.append(idx)
                y_bn.append(shoreline_position_by_label[bn_curve_label][t, idx])

                ax.axvline(
                    idx,
                    linestyle=":",
                    linewidth=1.0,
                    color="red",
                    alpha=0.45,
                    zorder=2,
                )

            ax.scatter(
                x_bn,
                y_bn,
                s=35,
                marker="o",
                color="red",
                edgecolor="black",
                linewidth=0.6,
                zorder=10,
                label="Historical BN domain",
            )

        ax.set_xlim(0, TOTAL_DOMAINS - 1)
        ax.set_ylim(y_min, y_max)

        ax.set_ylabel(
            "Absolute shoreline position, x_s (m)",
            fontsize=11,
            fontweight="bold",
        )

        ax.set_title(
            f"Pea Island absolute shoreline position with historical BN | {model_year}",
            fontsize=14,
            fontweight="bold",
            pad=12,
        )

        ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.35)
        ax.legend(loc="best", fontsize=9, frameon=True, framealpha=0.92)

        add_site_annotations(ax)
        add_real_domain_top_axis(ax)

        ax_bn.bar(
            domain_numbers,
            bn_total_m3,
            width=0.8,
            color="tab:blue",
            alpha=0.85,
        )

        total_bn_this_year = np.sum(bn_total_m3)

        if len(active_bn) > 0:
            bn_title = (
                f"Historical beach nourishment in {model_year}: YES | "
                f"Total volume = {total_bn_this_year:,.0f} m³"
            )
        else:
            bn_title = f"Historical beach nourishment in {model_year}: NO"

        ax_bn.set_title(bn_title, fontsize=12, fontweight="bold")
        ax_bn.set_ylabel("Historical BN\nvolume\n(m³/domain)", fontsize=10, fontweight="bold")
        ax_bn.set_xlabel("CASCADE domain index including buffers", fontsize=11, fontweight="bold")
        ax_bn.set_ylim(0, max_bn_volume * 1.25)
        ax_bn.grid(True, axis="y", linestyle=":", linewidth=0.8, alpha=0.35)

        xticks = np.arange(0, TOTAL_DOMAINS, DOMAIN_TICK_STEP)
        ax_bn.set_xticks(xticks)
        ax_bn.set_xticklabels([str(i) for i in xticks], rotation=45, ha="right", fontsize=9)

        for idx in active_bn:
            idx = int(idx)

            real_domain = padded_index_to_real_domain(idx)
            volume = bn_total_m3[idx]

            ax_bn.text(
                idx,
                volume + max_bn_volume * 0.025,
                f"D{real_domain}\n{volume/1000:.0f}k",
                ha="center",
                va="bottom",
                fontsize=8,
                color="0.25",
            )

        fig_out = os.path.join(yearly_dir, f"{filename_prefix}_{model_year}.png")
        fig.savefig(fig_out, dpi=220, bbox_inches="tight")
        plt.close(fig)

        png_files.append(fig_out)

    print("\nSaved yearly absolute shoreline position + historical BN plots to:")
    print(yearly_dir)

    if make_gif:
        try:
            import imageio.v2 as imageio

            images = [imageio.imread(f) for f in png_files]
            gif_out = os.path.join(output_dir, f"{filename_prefix}.gif")

            imageio.mimsave(
                gif_out,
                images,
                duration=gif_duration_seconds,
                loop=0,
            )

            print("\nSaved yearly absolute shoreline position + historical BN GIF:")
            print(gif_out)

        except ImportError:
            print("\nimageio is not installed. Absolute shoreline GIF was not created.")
            print("Install it with:")
            print("pip install imageio")

    return png_files


def plot_threshold_bn_effectiveness(
    threshold_bn_effective_log,
    output_dir,
    run_name,
):
    if len(threshold_bn_effective_log) == 0:
        print(f"\n{run_name}: No threshold BN records to plot.")
        return

    df = pd.DataFrame(threshold_bn_effective_log)

    real_df = df.dropna(subset=["real_domain"]).copy()

    if real_df.empty:
        print(f"\n{run_name}: Threshold BN log has no real-domain records to plot.")
        return

    real_df["real_domain"] = real_df["real_domain"].astype(int)

    years = np.arange(START_YEAR, END_YEAR + 1)
    domains = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)

    crossed_matrix = np.zeros((len(domains), len(years)))
    applied_matrix = np.zeros((len(domains), len(years)))

    for _, row in real_df.iterrows():
        dom = int(row["real_domain"])
        yr = int(row["model_year"])

        if dom in domains and yr in years:
            i = np.where(domains == dom)[0][0]
            j = np.where(years == yr)[0][0]

            crossed_matrix[i, j] = 1

            if bool(row["threshold_bn_applied_to_model"]):
                applied_matrix[i, j] = 1

    fig, ax = plt.subplots(figsize=(14, 8), constrained_layout=True)

    im = ax.imshow(
        crossed_matrix,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )

    ax.set_xticks(np.arange(len(years)))
    ax.set_xticklabels(years, rotation=45, ha="right")

    ax.set_yticks(np.arange(len(domains)))
    ax.set_yticklabels(domains)

    ax.set_xlabel("Model year", fontsize=12, fontweight="bold")
    ax.set_ylabel("Real Pea Island domain", fontsize=12, fontweight="bold")

    ax.set_title(
        f"Beach-width threshold crossings by year/domain\n{run_name}",
        fontsize=14,
        fontweight="bold",
    )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Threshold crossed: 1 = yes, 0 = no")

    out_png = os.path.join(
        output_dir,
        f"{run_name}_threshold_BN_crossing_heatmap.png",
    )

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"\nSaved threshold crossing heatmap: {out_png}")

    fig, ax = plt.subplots(figsize=(14, 8), constrained_layout=True)

    im = ax.imshow(
        applied_matrix,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )

    ax.set_xticks(np.arange(len(years)))
    ax.set_xticklabels(years, rotation=45, ha="right")

    ax.set_yticks(np.arange(len(domains)))
    ax.set_yticklabels(domains)

    ax.set_xlabel("Model year", fontsize=12, fontweight="bold")
    ax.set_ylabel("Real Pea Island domain", fontsize=12, fontweight="bold")

    ax.set_title(
        f"Threshold BN actually applied to CASCADE\n{run_name}",
        fontsize=14,
        fontweight="bold",
    )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Threshold BN applied: 1 = yes, 0 = no")

    out_png = os.path.join(
        output_dir,
        f"{run_name}_threshold_BN_applied_heatmap.png",
    )

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved threshold applied heatmap: {out_png}")


# =============================================================================
# CASCADE RUNNER
# =============================================================================

def run_cascade_simulation(
    nt,
    name,
    storm_file,
    alongshore_section_count,
    num_cores,
    beach_width_threshold,
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
    rebuild_dune_threshold,
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
    historical_nourishment_on_by_year=None,
    historical_nourishment_volume_by_year=None,
    historical_road_events=None,
):
    datadir = os.path.join(PROJECT_BASE_DIR, "data", "PeaIsland_init")

    dune_rebuild_log = []
    historical_nourishment_log = []
    historical_road_relocation_log = []
    threshold_bn_crossing_log = []
    threshold_bn_effective_log = []

    # Create a clean YAML copy for this run. The master YAML is never passed
    # directly to CASCADE, so Dmaxel cannot be repeatedly rewritten there.
    run_parameter_prefix, run_parameter_file = create_run_parameter_file(
        datadir=datadir,
        run_name=name,
        dmaxel_m=DMAXEL_M,
    )

    with open(run_parameter_file, "r") as f:
        yaml_params = yaml.safe_load(f)

    print("\nPARAMETER CHECK:")
    print(f"  Master YAML file: {PARAMETER_FILE}")
    print(f"  Barrier3D datadir: {datadir}")
    print(f"  Run parameter prefix: {run_parameter_prefix!r}")
    print(f"  Run YAML file:    {run_parameter_file}")
    print(f"  Run YAML exists:  {os.path.isfile(run_parameter_file)}")
    print(f"  Run YAML MHW:     {yaml_params.get('MHW')} m")
    print(f"  Run YAML BermEl:  {yaml_params.get('BermEl')} m")
    print(f"  Run YAML Dmaxel:  {yaml_params.get('Dmaxel')} m")
    print(f"  CASCADE MHW argument: {MHW_M} m")
    print(f"  CASCADE berm_elevation argument: {BERM_ELEVATION_M} m")
    print(f"  Requested Dmaxel: {DMAXEL_M} m")
    print(f"  APPLY_AUTOMATIC_BN_THRESHOLD: {APPLY_AUTOMATIC_BN_THRESHOLD}")

    expected_parameter_file = os.path.join(
        datadir,
        run_parameter_prefix + "-parameters.yaml",
    )
    print(f"  File Barrier3D will request: {expected_parameter_file}")

    if expected_parameter_file != run_parameter_file:
        raise RuntimeError(
            "Barrier3D parameter filename mismatch:\n"
            f"Created:  {run_parameter_file}\n"
            f"Expected: {expected_parameter_file}"
        )
    domain80_index = 71  # 71 buffers before real domain 80

    print("\nDOMAIN 80 FILE ACTUALLY PASSED TO CASCADE:")
    print(ELEVATION_FILE_PATHS[domain80_index])

    topo80 = np.load(ELEVATION_FILE_PATHS[domain80_index])

    print("Shape:", topo80.shape)
    print("Min:", np.nanmin(topo80))
    print("Max:", np.nanmax(topo80))

    print("\nCHECKING ALL CASCADE INPUT SHAPES...")

    for i, (dune_path, topo_path) in enumerate(
            zip(DUNE_FILE_PATHS, ELEVATION_FILE_PATHS)
    ):
        dune = np.load(dune_path)
        topo = np.load(topo_path)

        if dune.shape != (100,):
            raise ValueError(
                f"CASCADE index {i}: bad dune shape {dune.shape}\n"
                f"File: {dune_path}\n"
                "Expected (100,) for two dune rows."
            )

        if topo.shape != (200, 50):
            raise ValueError(
                f"CASCADE index {i}: bad topo shape {topo.shape}\n"
                f"File: {topo_path}"
            )

    print(
        f"OK: all {len(DUNE_FILE_PATHS)} domains have "
        "dune=(100,) and topo=(200,50)"
    )

    cascade = Cascade(
        datadir,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        # PREFIX ONLY: Barrier3D appends "-parameters.yaml" itself.
        parameter_file=run_parameter_prefix,

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
        alongshore_transport_module=True,

        beach_nourishment_module=beach_dune_manager_on,

        community_economics_module=False,

        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,

        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,

        nourishment_interval=None,
        nourishment_volume=nourishment_volume,

        overwash_filter=overwash_filter,
        overwash_to_dune=overwash_to_dune,

        enable_shoreline_offset=enable_shoreline_offset,
        shoreline_offset=shoreline_offset,
    )

    # =========================================================================
    # BARRIER3D ORIENTATION VERIFICATION
    # =========================================================================
    # Run immediately after CASCADE/Barrier3D initialization and before any
    # model updates. The check compares the exact topography files passed into
    # this run against the arrays stored inside the corresponding Barrier3D
    # objects. It is diagnostic only and does not modify any array.
    #
    # D80 and D81 are checked because their shared boundary is the first real-
    # domain connection. D119 is also checked as the final real domain.
    orientation_check_domains = [80, 81, 119]

    print("\nRUNNING BARRIER3D ORIENTATION VERIFICATION...")

    for real_domain_id in orientation_check_domains:
        barrier3d_index = START_REAL_INDEX + (
            real_domain_id - FIRST_FILE_NUMBER
        )

        if not (0 <= barrier3d_index < len(elevation_file)):
            raise IndexError(
                f"Domain {real_domain_id} maps to invalid Barrier3D index "
                f"{barrier3d_index}."
            )

        saved_topography_path = elevation_file[barrier3d_index]

        print(
            f"\nDomain D{real_domain_id} -> Barrier3D index "
            f"{barrier3d_index}"
        )
        print(f"Actual file passed to CASCADE: {saved_topography_path}")

        check_barrier3d_orientation(
            cascade_model=cascade,
            barrier3d_index=barrier3d_index,
            saved_topography_path=saved_topography_path,
        )

    print("BARRIER3D ORIENTATION VERIFICATION COMPLETE.\n")

    roadways = get_roadway_manager_list(cascade)
    first_managed_index = next(
        i
        for i, management_on in enumerate(roadway_management_on)
        if management_on
    )
    test_roadway = roadways[first_managed_index]

    print("\nROADWAY INSTANCE CHECK")
    print("CASCADE index:", first_managed_index)
    print("Object type:", type(test_roadway))
    print("Object module:", type(test_roadway).__module__)
    print(
        "request_relocation available:",
        hasattr(test_roadway, "request_relocation"),
    )
    if hasattr(test_roadway, "request_relocation"):
        print(
            "request_relocation signature:",
            inspect.signature(test_roadway.request_relocation),
        )
    print()

    print("\nBarrier3D berm elevation check:")
    print("cascade.barrier3d[0].BermEl raw =", cascade.barrier3d[0].BermEl)
    print("cascade.barrier3d[0].BermEl meters =", cascade.barrier3d[0].BermEl * 10)

    for i in range(START_REAL_INDEX, END_REAL_INDEX):
        real_id = FIRST_FILE_NUMBER + (i - START_REAL_INDEX)
        print(
            f"Domain {real_id}: "
            f"BermEl_raw={cascade.barrier3d[i].BermEl}, "
            f"BermEl_m={cascade.barrier3d[i].BermEl * 10:.3f}"
        )

    print("\nBarrier3D Dmaxel check:")
    for i in range(START_REAL_INDEX, END_REAL_INDEX):
        real_id = FIRST_FILE_NUMBER + (i - START_REAL_INDEX)
        b3d = cascade.barrier3d[i]

        if hasattr(b3d, "Dmaxel"):
            dmaxel_internal = b3d.Dmaxel
        elif hasattr(b3d, "_Dmaxel"):
            dmaxel_internal = b3d._Dmaxel
        else:
            dmaxel_internal = None

        print(
            f"Domain {real_id}: "
            f"Dmaxel_input={DMAXEL_M:.3f} m, "
            f"Dmaxel_internal={dmaxel_internal}"
        )

    # Confirm that the run-specific YAML still contains the requested input.
    # Even if Barrier3D changes this copy, the master YAML remains protected.
    with open(run_parameter_file, "r") as f:
        yaml_after_initialization = yaml.safe_load(f)

    dmaxel_after = yaml_after_initialization.get("Dmaxel")
    print(
        "\nRun YAML after CASCADE initialization: "
        f"Dmaxel={dmaxel_after}"
    )
    if float(dmaxel_after) != float(DMAXEL_M):
        print(
            "WARNING: Barrier3D changed Dmaxel in the run-specific YAML copy. "
            "The master YAML is still safe and the value will be reset for "
            "the next run."
        )

    cascade._sandbag_management_on = [False] * alongshore_section_count
    cascade._sandbag_elevation = [0.0] * alongshore_section_count

    dune_rebuild_threshold_abs = rebuild_dune_threshold + (
        cascade.barrier3d[0].BermEl * 10
    )

    print(f"\n{name}: absolute dune rebuild threshold = {dune_rebuild_threshold_abs:.3f} m MHW")

    managed_mask = np.asarray(beach_dune_manager_on, dtype=bool)
    n_managed = int(np.sum(managed_mask))



    for time_step in range(nt - 1):
        current_year = START_YEAR + time_step

        cascade.nourish_now = np.zeros(alongshore_section_count)
        cascade.rebuild_dune_now = np.zeros(alongshore_section_count)

        # ---------------------------------------------------------------------
        # HISTORICAL BN
        # ---------------------------------------------------------------------
        active_historical_bn = np.array([], dtype=int)
        pre_historical_bn_state = {}

        if historical_nourishment_on_by_year is not None:
            if current_year in historical_nourishment_on_by_year:

                nourish_now = np.asarray(
                    historical_nourishment_on_by_year[current_year],
                    dtype=float,
                )

                nourish_volume = np.asarray(
                    historical_nourishment_volume_by_year[current_year],
                    dtype=float,
                )

                active_historical_bn = np.where(nourish_now == 1)[0]

                if active_historical_bn.size > 0:
                    print(f"\n{name}: applying HISTORICAL BN in {current_year}")
                    cascade.nourish_now = nourish_now

                    for iB3D in active_historical_bn:
                        iB3D = int(iB3D)

                        if not managed_mask[iB3D]:
                            raise RuntimeError(
                                f"Historical BN requested for unmanaged CASCADE "
                                f"index {iB3D}, domain "
                                f"{padded_index_to_real_domain(iB3D)}."
                            )

                        if hasattr(
                            cascade.nourishments[iB3D],
                            "_nourishment_volume",
                        ):
                            cascade.nourishments[
                                iB3D
                            ]._nourishment_volume = float(
                                nourish_volume[iB3D]
                            )
                        elif hasattr(
                            cascade.nourishments[iB3D],
                            "nourishment_volume",
                        ):
                            cascade.nourishments[
                                iB3D
                            ].nourishment_volume = float(
                                nourish_volume[iB3D]
                            )
                        else:
                            raise AttributeError(
                                "Could not find nourishment volume attribute on "
                                f"cascade.nourishments[{iB3D}]"
                            )

                        b3d = cascade.barrier3d[iB3D]
                        nourishment = cascade.nourishments[iB3D]

                        current_index = max(int(b3d.time_index) - 1, 0)
                        current_index = min(
                            current_index,
                            len(nourishment.beach_width) - 1,
                            len(get_x_s_TS(b3d)) - 1,
                        )

                        pre_historical_bn_state[iB3D] = {
                            "beach_width_before_m": float(
                                nourishment.beach_width[current_index]
                            ),
                            "shoreline_before_m": float(
                                get_x_s_TS(b3d)[current_index] * 10.0
                            ),
                            "volume_m3_per_m": float(
                                nourish_volume[iB3D]
                            ),
                        }

        # ---------------------------------------------------------------------
        # HISTORICAL ROAD RELOCATION — QUEUE ONE ROADWAY-MANAGER REQUEST
        # ---------------------------------------------------------------------
        active_historical_road_events = queue_historical_road_relocations(
            cascade=cascade,
            current_year=current_year,
            historical_road_events=historical_road_events,
            alongshore_section_count=alongshore_section_count,
            roadway_management_on=roadway_management_on,
        )

        for event_record in active_historical_road_events:
            event_record["run_name"] = name
            event_record["time_step"] = time_step

        print(f"\r{name}: Time Step {time_step}", end="", flush=True)

        # =============================================================================
        # CHECK HOW BARRIER3D ACTUALLY LOADED THE TWO DUNE ROWS
        # =============================================================================

        print("\n" + "=" * 80)
        print("INTERNAL BARRIER3D DUNE CHECK")
        print("=" * 80)

        for domain_id in [80, 81, 90, 100, 110, 119]:

            idx = START_REAL_INDEX + (domain_id - FIRST_FILE_NUMBER)
            b3d = cascade.barrier3d[idx]

            dune0 = np.asarray(b3d.DuneDomain[0], dtype=float)

            print(f"\nDomain {domain_id} | CASCADE index {idx}")
            print("DuneDomain[0] shape:", dune0.shape)

            print(
                "Internal DuneWidth:",
                getattr(b3d, "_DuneWidth",
                        getattr(b3d, "DuneWidth", "NOT FOUND"))
            )

            print(
                "DuneParamMultipleRows:",
                getattr(
                    b3d,
                    "_DuneParamMultipleRows",
                    getattr(b3d, "DuneParamMultipleRows", "NOT FOUND")
                )
            )

            if dune0.ndim == 2:

                for row in range(dune0.shape[1]):
                    print(
                        f"  Row {row}: "
                        f"min={np.nanmin(dune0[:, row]):.4f}, "
                        f"mean={np.nanmean(dune0[:, row]):.4f}, "
                        f"max={np.nanmax(dune0[:, row]):.4f}"
                    )

                if dune0.shape[1] >= 2:
                    print(
                        "  Rows identical:",
                        np.allclose(
                            dune0[:, 0],
                            dune0[:, 1],
                            equal_nan=True
                        )
                    )

                    print(
                        "  Mean absolute row difference:",
                        np.nanmean(
                            np.abs(dune0[:, 0] - dune0[:, 1])
                        )
                    )

                crest = np.nanmax(dune0, axis=1)

                print(
                    f"  Crest: min={np.nanmin(crest):.4f}, "
                    f"mean={np.nanmean(crest):.4f}, "
                    f"max={np.nanmax(crest):.4f}"
                )

        print("=" * 80)



        # Keep the original CASCADE API. RoadwayManager consumes its queued
        # historical relocation request inside this normal update.
        cascade.update()

        # =====================================================================
        # RUNTIME / EARLY-STOP DIAGNOSTIC
        # =====================================================================

        time_indices = np.array(
            [b3d.time_index for b3d in cascade.barrier3d],
            dtype=int
        )

        drowned_indices = [
            i
            for i, b3d in enumerate(cascade.barrier3d)
            if getattr(b3d, "drown_break", 0)
        ]

        print(
            f"\nYEAR {current_year}: "
            f"b3d_break={cascade.b3d_break} | "
            f"time_index min={time_indices.min()} "
            f"max={time_indices.max()} | "
            f"drowned={drowned_indices}"
        )

        if cascade.b3d_break:
            print("\n!!! CASCADE STOPPED EARLY !!!")

            for i in drowned_indices:
                real_domain = padded_index_to_real_domain(i)

                print(
                    f"  CASCADE index {i}, "
                    f"real domain={real_domain}, "
                    f"InteriorDomain shape="
                    f"{cascade.barrier3d[i].InteriorDomain.shape}"
                )

            break

        if cascade.b3d_break:
            break

        # ---------------------------------------------------------------------
        # VERIFY QUEUED ROADWAY-MANAGER HISTORICAL RELOCATION AFTER THE UPDATE
        # ---------------------------------------------------------------------
        if active_historical_road_events:
            roadways = get_roadway_manager_list(cascade)

            for event_record in active_historical_road_events:
                iB3D = int(event_record["padded_index"])
                roadway = roadways[iB3D]

                ts_index = max(int(roadway._time_index) - 1, 0)
                road_relocated_ts = getattr(roadway, "_road_relocated_TS", None)
                if road_relocated_ts is not None and ts_index < len(road_relocated_ts):
                    relocation_applied = bool(road_relocated_ts[ts_index])
                else:
                    relocation_applied = False

                record = dict(event_record)
                record.update(
                    road_setback_after_update_m=float(roadway._road_setback),
                    road_width_after_update_m=float(roadway._road_width),
                    road_elevation_after_update_m=float(roadway._road_ele),
                    relocation_break=bool(
                        getattr(roadway, "_relocation_break", 0)
                    ),
                    drown_break=bool(getattr(roadway, "_drown_break", 0)),
                    prescribed_event_requested=True,
                    prescribed_event_applied=relocation_applied,
                )
                historical_road_relocation_log.append(record)

                print(
                    f"\nHistorical relocation verification | "
                    f"year={current_year} | domain={record['real_domain']} | "
                    f"target={record['requested_relocation_setback_m']:.2f} m | "
                    f"after update={record['road_setback_after_update_m']:.2f} m | "
                    f"applied={record['prescribed_event_applied']} | "
                    f"relocation_break={record['relocation_break']} | "
                    f"drown_break={record['drown_break']}"
                )

        # ---------------------------------------------------------------------
        # VERIFY HISTORICAL BN AFTER THE UPDATE
        # ---------------------------------------------------------------------
        for iB3D in active_historical_bn:
            iB3D = int(iB3D)

            b3d = cascade.barrier3d[iB3D]
            nourishment = cascade.nourishments[iB3D]

            new_index = max(int(b3d.time_index) - 1, 0)
            new_index = min(
                new_index,
                len(nourishment.beach_width) - 1,
                len(get_x_s_TS(b3d)) - 1,
            )

            beach_width_after_m = float(
                nourishment.beach_width[new_index]
            )
            shoreline_after_m = float(
                get_x_s_TS(b3d)[new_index] * 10.0
            )

            before = pre_historical_bn_state[iB3D]
            beach_width_change_m = (
                beach_width_after_m
                - before["beach_width_before_m"]
            )
            shoreline_position_change_m = (
                shoreline_after_m
                - before["shoreline_before_m"]
            )

            effect_detected = bool(
                abs(beach_width_change_m) > 1e-6
                or abs(shoreline_position_change_m) > 1e-6
            )

            real_domain = padded_index_to_real_domain(iB3D)

            historical_nourishment_log.append(
                dict(
                    run_name=name,
                    model_year=current_year,
                    time_step=time_step,
                    padded_index=iB3D,
                    real_domain=real_domain,
                    nourishment_type="historical",
                    nourishment_volume_m3_per_m=before[
                        "volume_m3_per_m"
                    ],
                    nourishment_volume_total_m3=before[
                        "volume_m3_per_m"
                    ] * DOMAIN_LENGTH_M,
                    beach_width_before_m=before[
                        "beach_width_before_m"
                    ],
                    beach_width_after_m=beach_width_after_m,
                    beach_width_change_m=beach_width_change_m,
                    shoreline_before_m=before[
                        "shoreline_before_m"
                    ],
                    shoreline_after_m=shoreline_after_m,
                    shoreline_position_change_m=(
                        shoreline_position_change_m
                    ),
                    nourishment_effect_detected=effect_detected,
                )
            )

            print(
                f"\nBN verification | year={current_year} "
                f"| domain={real_domain} "
                f"| beach-width change={beach_width_change_m:+.3f} m "
                f"| shoreline change="
                f"{shoreline_position_change_m:+.3f} m "
                f"| effect detected={effect_detected}"
            )

        t = cascade.barrier3d[0].time_index

        tmp_rebuild_dune = np.zeros(alongshore_section_count)
        tmp_threshold_nourish_now = np.zeros(alongshore_section_count)

        for iB3D in range(alongshore_section_count):

            if cascade.community_break[iB3D]:
                continue

            if not beach_dune_manager_on[iB3D]:
                continue

            # -----------------------------------------------------------------
            # THRESHOLD BN DIAGNOSTIC
            # -----------------------------------------------------------------
            if APPLY_AUTOMATIC_BN_THRESHOLD:

                beach_width_now = float(cascade.nourishments[iB3D].beach_width[t - 1])
                threshold_now = float(beach_width_threshold[iB3D])

                if beach_width_now < threshold_now:

                    tmp_threshold_nourish_now[iB3D] = 1

                    real_domain = padded_index_to_real_domain(iB3D)

                    threshold_bn_crossing_log.append(
                        dict(
                            run_name=name,
                            model_year=START_YEAR + time_step + 1,
                            time_step=time_step + 1,
                            padded_index=iB3D,
                            real_domain=real_domain,
                            beach_width_m=beach_width_now,
                            beach_width_threshold_m=threshold_now,
                            threshold_crossed=True,
                        )
                    )

            # -----------------------------------------------------------------
            # DUNE REBUILD LOGIC
            # -----------------------------------------------------------------
            DuneDomainCrest = (
                cascade.barrier3d[iB3D].DuneDomain[t - 1, :, :].max(axis=1)
            )

            DuneCrestMin = (
                np.min(DuneDomainCrest) + cascade.barrier3d[iB3D].BermEl
            ) * 10

            if DuneCrestMin < dune_rebuild_threshold_abs:
                tmp_rebuild_dune[iB3D] = 1

                real_domain = padded_index_to_real_domain(iB3D)

                dune_rebuild_log.append(
                    dict(
                        run_name=name,
                        model_year=START_YEAR + time_step + 1,
                        time_step=time_step + 1,
                        padded_index=iB3D,
                        real_domain=real_domain,
                        dune_crest_min_m=float(DuneCrestMin),
                        rebuild_threshold_m=float(dune_rebuild_threshold_abs),
                        roadway_management_on=bool(roadway_management_on[iB3D]),
                        beach_dune_manager_on=bool(beach_dune_manager_on[iB3D]),
                        road_setback=float(road_setback[iB3D]),
                    )
                )

        # ---------------------------------------------------------------------
        # APPLY THRESHOLD BN USING CURRENT OCRACOKE-STYLE np.all LOGIC
        # ---------------------------------------------------------------------
        if APPLY_AUTOMATIC_BN_THRESHOLD and n_managed > 0:

            n_triggered = int(np.sum(tmp_threshold_nourish_now[managed_mask] == 1))

            threshold_bn_applied = False

            if np.all(tmp_threshold_nourish_now[managed_mask] == 1):
                cascade.nourish_now = tmp_threshold_nourish_now
                threshold_bn_applied = True

            active_threshold_bn = np.where(tmp_threshold_nourish_now == 1)[0]

            for idx in active_threshold_bn:
                idx = int(idx)

                real_domain = padded_index_to_real_domain(idx)

                if hasattr(cascade.nourishments[idx], "beach_width"):
                    beach_width_logged = float(cascade.nourishments[idx].beach_width[t - 1])
                else:
                    beach_width_logged = np.nan

                threshold_bn_effective_log.append(
                    dict(
                        run_name=name,
                        model_year=START_YEAR + time_step + 1,
                        time_step=time_step + 1,
                        padded_index=idx,
                        real_domain=real_domain,
                        beach_width_m=beach_width_logged,
                        beach_width_threshold_m=float(beach_width_threshold[idx]),
                        threshold_crossed=True,
                        n_managed_domains=n_managed,
                        n_triggered_domains=n_triggered,
                        threshold_bn_applied_to_model=threshold_bn_applied,
                        applied_rule="np.all managed domains below threshold",
                    )
                )

            if n_triggered > 0:
                print(
                    f"\n{name}: Threshold BN diagnostic year {START_YEAR + time_step + 1}: "
                    f"{n_triggered}/{n_managed} managed domains below threshold | "
                    f"applied_to_model={threshold_bn_applied}"
                )

        # ---------------------------------------------------------------------
        # APPLY DUNE REBUILD USING OCRACOKE-STYLE np.all LOGIC
        # ---------------------------------------------------------------------
        if n_managed > 0 and np.any(tmp_rebuild_dune == 1):
            cascade.rebuild_dune_now = tmp_rebuild_dune

    save_path = os.path.join(OUTPUT_BASE_DIR, name)
    os.makedirs(save_path, exist_ok=True)

    cascade.save(save_path)
    print(f"\nSaved: {save_path}")

    cascade.dune_rebuild_log = dune_rebuild_log
    cascade.historical_nourishment_log = historical_nourishment_log
    cascade.historical_road_relocation_log = historical_road_relocation_log
    cascade.threshold_bn_crossing_log = threshold_bn_crossing_log
    cascade.threshold_bn_effective_log = threshold_bn_effective_log

    # -------------------------------------------------------------------------
    # SAVE HISTORICAL ROAD-RELOCATION LOG
    # -------------------------------------------------------------------------
    if len(historical_road_relocation_log) > 0:
        relocation_df = pd.DataFrame(historical_road_relocation_log)
        relocation_csv = os.path.join(
            save_path,
            "historical_road_relocation_log.csv",
        )
        relocation_df.to_csv(relocation_csv, index=False)
        print(f"Saved historical road-relocation log: {relocation_csv}")

    # -------------------------------------------------------------------------
    # SAVE DUNE REBUILD LOG
    # -------------------------------------------------------------------------
    if len(dune_rebuild_log) > 0:
        rebuild_df = pd.DataFrame(dune_rebuild_log)
        rebuild_csv = os.path.join(
            save_path,
            "dune_rebuild_log.csv",
        )
        rebuild_df.to_csv(rebuild_csv, index=False)

        print("\n" + "=" * 80)
        print(f"DUNE REBUILD CONDITION SUMMARY FOR RUN: {name}")
        print("=" * 80)
        print(f"Number of domain-year records below threshold: {len(rebuild_df)}")
        print("\nBy year:")
        print(rebuild_df.groupby("model_year").size())
        print("\nBy real domain:")
        print(rebuild_df.groupby("real_domain").size())
        print(f"\nSaved dune rebuild log: {rebuild_csv}")
    else:
        print(f"\n{name}: No dune rebuild condition was detected.")

    # -------------------------------------------------------------------------
    # SAVE HISTORICAL BN LOG
    # -------------------------------------------------------------------------
    if len(historical_nourishment_log) > 0:
        hist_bn_df = pd.DataFrame(historical_nourishment_log)
        hist_bn_csv = os.path.join(
            save_path,
            "historical_BN_log.csv",
        )
        hist_bn_df.to_csv(hist_bn_csv, index=False)

        print("\n" + "=" * 80)
        print(f"HISTORICAL BN SUMMARY FOR RUN: {name}")
        print("=" * 80)
        print(f"Number of domain-year historical BN records: {len(hist_bn_df)}")
        print("\nBy year:")
        print(hist_bn_df.groupby("model_year").size())
        print("\nBy real domain:")
        print(hist_bn_df.groupby("real_domain").size())
        if "nourishment_effect_detected" in hist_bn_df.columns:
            print("\nEffect detection:")
            print(
                hist_bn_df[
                    "nourishment_effect_detected"
                ].value_counts(dropna=False)
            )

        print(f"\nSaved historical BN log: {hist_bn_csv}")
    else:
        print(f"\n{name}: No historical BN was applied.")

    # -------------------------------------------------------------------------
    # SAVE THRESHOLD CROSSING LOG
    # -------------------------------------------------------------------------
    if len(threshold_bn_crossing_log) > 0:
        threshold_cross_df = pd.DataFrame(threshold_bn_crossing_log)
        threshold_cross_csv = os.path.join(
            save_path,
            "threshold_BN_crossing_log.csv",
        )
        threshold_cross_df.to_csv(threshold_cross_csv, index=False)

        print("\n" + "=" * 80)
        print(f"THRESHOLD BN CROSSING SUMMARY FOR RUN: {name}")
        print("=" * 80)
        print(f"Number of threshold-crossing domain-year records: {len(threshold_cross_df)}")
        print("\nBy year:")
        print(threshold_cross_df.groupby("model_year").size())
        print("\nBy real domain:")
        print(threshold_cross_df.groupby("real_domain").size())
        print(f"\nSaved threshold BN crossing log: {threshold_cross_csv}")
    else:
        print(f"\n{name}: No beach-width threshold crossings were detected.")

    # -------------------------------------------------------------------------
    # SAVE THRESHOLD EFFECTIVE LOG
    # -------------------------------------------------------------------------
    if len(threshold_bn_effective_log) > 0:
        threshold_effective_df = pd.DataFrame(threshold_bn_effective_log)

        threshold_effective_csv = os.path.join(
            save_path,
            "threshold_BN_effective_log.csv",
        )

        threshold_effective_df.to_csv(threshold_effective_csv, index=False)

        print("\n" + "=" * 80)
        print(f"THRESHOLD BN EFFECTIVENESS SUMMARY FOR RUN: {name}")
        print("=" * 80)

        print(f"Number of threshold-crossing domain-year records: {len(threshold_effective_df)}")

        print("\nBy year:")
        print(
            threshold_effective_df
            .groupby("model_year")
            .agg(
                triggered_domains=("padded_index", "count"),
                managed_domains=("n_managed_domains", "max"),
                applied_to_model=("threshold_bn_applied_to_model", "max"),
            )
        )

        print("\nBy real domain:")
        print(
            threshold_effective_df
            .groupby("real_domain")
            .agg(
                times_triggered=("model_year", "count"),
                ever_applied_to_model=("threshold_bn_applied_to_model", "max"),
            )
        )

        print(f"\nSaved threshold BN effectiveness log: {threshold_effective_csv}")

        plot_threshold_bn_effectiveness(
            threshold_bn_effective_log=threshold_bn_effective_log,
            output_dir=OUTPUT_BASE_DIR,
            run_name=name,
        )

    else:
        print(f"\n{name}: No threshold BN effectiveness records were created.")
        print(f"APPLY_AUTOMATIC_BN_THRESHOLD={APPLY_AUTOMATIC_BN_THRESHOLD}")

    return cascade


# =============================================================================
# MAIN
# =============================================================================

def main():
    RMIN = [0.55] * TOTAL_DOMAINS
    RMAX = [0.95] * TOTAL_DOMAINS

    DUNE_DESIGN_ELEVATION = [DUNE_REBUILD_HEIGHT] * TOTAL_DOMAINS
    DUNE_MINIMUM_ELEVATION = [REBUILD_ELEV_THRESHOLD] * TOTAL_DOMAINS

    BEACH_WIDTH_THRESHOLD = [BEACH_WIDTH_THRESHOLD_VALUE] * TOTAL_DOMAINS

    ROAD_ELEVATION_REAL = [
        1.0776, 1.2634, 1.4387, 1.1686, 1.3991,
        1.2738, 1.2564, 1.3914, 1.8516, 1.3563,
        0.8953, 0.7497, 0.9807, 0.7779, 0.8456,
        0.8496, 0.9289, 0.9280, 1.0394, 1.3601,
        1.2931, 1.4435, 1.5030, 1.4476, 1.3255,
        1.3430, 1.1604, 1.6557, 1.0107, 1.1600,
        1.2845, 0.9283, 1.8228, 1.4769, 1.2812,
        1.2199, 1.1808, 0.9350, 1.3043, 1.6772,
    ]

    ROAD_ELEVATION = [0.0] * TOTAL_DOMAINS
    ROAD_ELEVATION[START_REAL_INDEX:END_REAL_INDEX] = ROAD_ELEVATION_REAL

    ROAD_WIDTH = [10.0] * TOTAL_DOMAINS

    OVERWASH_FILTER = 0
    OVERWASH_TO_DUNE = 9

    NOURISHMENT_VOLUME = 100

    HIST_NOURISH_ON, HIST_NOURISH_VOLUME = build_nourishment_arrays_from_manual_inputs()

    cmp_df = load_cmp_domain_rates(
        file_path=CMP_FILE,
        sheet_name=CMP_SHEET_NAME,
        domain_col=CMP_DOMAIN_COLUMN,
        rate_col=CMP_RATE_COLUMN,
        std_col=CMP_STD_COLUMN,
        count_col=CMP_COUNT_COLUMN,
        flip_sign=CMP_FLIP_SIGN,
    )

    NATURAL_ROADWAY_ON = [False] * TOTAL_DOMAINS
    NATURAL_BEACH_DUNE_ON = [False] * TOTAL_DOMAINS

    ROADWAY_BEACH_DUNE_ON = [False] * TOTAL_DOMAINS
    ROADWAY_BEACH_DUNE_ON[START_REAL_INDEX:END_REAL_INDEX] = [True] * NUM_REAL_DOMAINS

    scenarios = [
        # dict(
        #     scenario_name="natural",
        #     label="Natural",
        #     roadway_on=NATURAL_ROADWAY_ON,
        #     beach_dune_on=NATURAL_BEACH_DUNE_ON,
        #     use_historical_nourishment=False,
        # ),
        # dict(
        #     scenario_name="roadwayMngmt",
        #     label="Roadway management",
        #     roadway_on=ROADWAY_MANAGEMENT_ON,
        #     beach_dune_on=ROADWAY_BEACH_DUNE_ON,
        #     use_historical_nourishment=False,
        # ),
        dict(
            scenario_name="duneRebuild_HIST_BN",
            label="Dune rebuilding + historical BN",
            roadway_on=ROADWAY_MANAGEMENT_ON,
            beach_dune_on=ROADWAY_BEACH_DUNE_ON,
            use_historical_nourishment=True,
            use_historical_road_relocation=True,
        ),
    ]

    rate_profiles = {}
    cascades_by_label = {}

    for scenario in scenarios:
        for Hs in WAVE_HEIGHTS_TO_TEST:
            run_name = (
                f"PEA_{START_YEAR}_{END_YEAR}_"
                f"DuneRebuild_HistBN_Hs{Hs:.1f}"
            ).replace(".", "p")

            print("\n" + "=" * 80)
            print(f"Starting run: {run_name}")
            print("=" * 80)

            cascade = run_cascade_simulation(
                nt=TIME_STEP_COUNT,
                name=run_name,
                storm_file=STORM_FILE,
                alongshore_section_count=TOTAL_DOMAINS,
                num_cores=NUM_CORES,

                beach_width_threshold=BEACH_WIDTH_THRESHOLD,

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

                rebuild_dune_threshold=REBUILD_DUNE_THRESHOLD,

                roadway_management_on=scenario["roadway_on"],
                beach_dune_manager_on=scenario["beach_dune_on"],

                sea_level_rise_rate=SEA_LEVEL_RISE_RATE,
                sea_level_constant=SEA_LEVEL_CONSTANT,

                enable_shoreline_offset=True,
                shoreline_offset=dune_offset_dam,

                wave_height=Hs,
                wave_period=FIXED_WAVE_PERIOD,
                wave_asymmetry=FIXED_WAVE_ASYMMETRY,
                wave_angle_high_fraction=FIXED_WAVE_ANGLE_HIGH_FRACTION,

                historical_nourishment_on_by_year=(
                    HIST_NOURISH_ON if scenario["use_historical_nourishment"] else None
                ),
                historical_nourishment_volume_by_year=(
                    HIST_NOURISH_VOLUME if scenario["use_historical_nourishment"] else None
                ),
                historical_road_events=(
                    HISTORICAL_ROAD_EVENTS
                    if scenario.get("use_historical_road_relocation", False)
                    else None
                ),
            )

            # =====================================================================
            # CHECK ACTUAL BARRIER3D OVERWASH
            # =====================================================================

            print("\n" + "=" * 80)
            print("OVERWASH CHECK — REAL DOMAINS 80–119")
            print("=" * 80)

            for domain_id in range(80, 120):

                idx = START_REAL_INDEX + (domain_id - FIRST_FILE_NUMBER)

                b3d = cascade.barrier3d[idx]

                qow = np.asarray(b3d.QowTS, dtype=float)

                print(
                    f"Domain {domain_id}: "
                    f"QowTS length={len(qow)}, "
                    f"max={np.nanmax(qow):.3f} m3/m"
                )

                # QowTS[0] = initial condition
                annual_qow = qow[1:]

                years = np.arange(
                    START_YEAR,
                    START_YEAR + len(annual_qow)
                )

                active = annual_qow > 0

                if np.any(active):

                    for year, value in zip(
                            years[active],
                            annual_qow[active]
                    ):
                        print(
                            f"    OVERWASH: {year} "
                            f"Qow={value:.3f} m3/m"
                        )

            print("=" * 80)

            # Continue with your existing code
            cascades_by_label[scenario["label"]] = cascade

            shoreline_m = build_shoreline_matrix(
                cascade,
                to_meters=TO_METERS
            )

            cascades_by_label[scenario["label"]] = cascade

            shoreline_m = build_shoreline_matrix(cascade, to_meters=TO_METERS)

            change_rate = calculate_cascade_shoreline_change_rate(
                shoreline_m=shoreline_m,
                method=CASCADE_RATE_METHOD,
                flip_sign=FLIP_SIGN_MODEL,
            )

            rate_profiles[(scenario["label"], Hs)] = change_rate

    # =============================================================================
    # YEARLY RELATIVE SHORELINE CHANGE + HISTORICAL BN PLOTS
    # =============================================================================

    threshold_tag = "thresholdBN_ON" if APPLY_AUTOMATIC_BN_THRESHOLD else "thresholdBN_OFF"

    yearly_prefix = (
        f"PEA_{START_YEAR}_{END_YEAR}_yearly_relative_shoreline_HIST_BN_{threshold_tag}"
    )

    # plot_yearly_relative_shoreline_and_bn(
    #     cascades_by_label=cascades_by_label,
    #     hist_nourish_volume_by_year=HIST_NOURISH_VOLUME,
    #     output_dir=OUTPUT_BASE_DIR,
    #     filename_prefix=yearly_prefix,
    #     flip_sign=True,
    #     make_gif=MAKE_YEARLY_BN_GIF,
    #     gif_duration_seconds=GIF_DURATION_SECONDS,
    # )

    # =============================================================================
    # YEARLY ABSOLUTE SHORELINE POSITION + HISTORICAL BN PLOTS
    # =============================================================================

    absolute_yearly_prefix = (
        f"PEA_{START_YEAR}_{END_YEAR}_yearly_absolute_shoreline_HIST_BN_{threshold_tag}"
    )

    # plot_yearly_absolute_shoreline_and_bn(
    #     cascades_by_label=cascades_by_label,
    #     hist_nourish_volume_by_year=HIST_NOURISH_VOLUME,
    #     output_dir=OUTPUT_BASE_DIR,
    #     filename_prefix=absolute_yearly_prefix,
    #     make_gif=MAKE_YEARLY_BN_GIF,
    #     gif_duration_seconds=GIF_DURATION_SECONDS,
    # )

    # =============================================================================
    # CMP 1992–2024 VS CASCADE 1992–2007 CALIBRATION COMPARISON
    # =============================================================================

    real_domain_ids = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)

    comparison_df = cmp_df.set_index("domain_id").reindex(real_domain_ids).reset_index()
    comparison_df = comparison_df.rename(columns={"index": "domain_id"})

    fig, ax = plt.subplots(figsize=(18, 7), constrained_layout=True)

    cmp_has_std = comparison_df["cmp_std_m_yr"].notna().any()

    if cmp_has_std:
        ax.errorbar(
            comparison_df["domain_id"],
            comparison_df["cmp_rate_m_yr"],
            yerr=comparison_df["cmp_std_m_yr"],
            linestyle="--",
            marker="o",
            linewidth=1.8,
            markersize=5.5,
            capsize=3,
            label="CMP LRR 1992–2024 (mean ± 1 SD)",
            zorder=10,
        )
    else:
        ax.plot(
            comparison_df["domain_id"],
            comparison_df["cmp_rate_m_yr"],
            linestyle="--",
            marker="o",
            linewidth=1.8,
            markersize=5.5,
            label="CMP LRR 1992–2024",
            zorder=10,
        )

    statistics_rows = []

    for (scenario_label, Hs), full_rate in rate_profiles.items():
        cascade_real_rate = np.asarray(full_rate, dtype=float)[
            START_REAL_INDEX:END_REAL_INDEX
        ]

        safe_scenario = (
            scenario_label.lower()
            .replace(" ", "_")
            .replace("+", "plus")
            .replace("–", "-")
        )
        model_column = f"cascade_{safe_scenario}_Hs{Hs:.1f}_m_yr"
        comparison_df[model_column] = cascade_real_rate

        stats = calculate_comparison_statistics(
            observed=comparison_df["cmp_rate_m_yr"].to_numpy(),
            modeled=cascade_real_rate,
        )

        statistics_rows.append(
            {
                "scenario": scenario_label,
                "wave_height_m": Hs,
                "rate_method": CASCADE_RATE_METHOD,
                **stats,
            }
        )

        stat_label = (
            f"{scenario_label}, Hs={Hs:g} m, "
            f"CASCADE {START_YEAR}–{END_YEAR} {CASCADE_RATE_METHOD.upper()} | "
            f"RMSE={stats['rmse']:.2f}, bias={stats['bias']:.2f} m/yr"
        )

        ax.plot(
            real_domain_ids,
            cascade_real_rate,
            linewidth=2.4,
            label=stat_label,
        )

        comparison_df[f"residual_{safe_scenario}_Hs{Hs:.1f}_m_yr"] = (
            cascade_real_rate - comparison_df["cmp_rate_m_yr"].to_numpy()
        )

    ax.axhline(
        0.0,
        linestyle="--",
        linewidth=1.2,
        color="0.25",
        alpha=0.75,
    )

    ax.set_xticks(
        np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP)
    )
    ax.set_xlim(FIRST_FILE_NUMBER - 0.5, LAST_FILE_NUMBER + 0.5)

    ax.set_xlabel(
        "Real Pea Island CASCADE domain ID",
        fontsize=12,
        fontweight="bold",
    )

    ax.set_ylabel(
        "Shoreline change rate (m/yr)\nAccretion = positive; erosion = negative",
        fontsize=12,
        fontweight="bold",
    )

    ax.set_title(
        f"CMP 1992–2024 vs CASCADE {START_YEAR}–{END_YEAR} calibration\n"
        f"CASCADE rate method: {CASCADE_RATE_METHOD.upper()} | "
        f"SLR={SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr",
        fontsize=14,
        fontweight="bold",
        pad=14,
    )

    ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.35)
    ax.minorticks_on()
    ax.grid(True, which="minor", linestyle=":", linewidth=0.5, alpha=0.18)

    all_y_values = [comparison_df["cmp_rate_m_yr"].to_numpy()]
    if cmp_has_std:
        all_y_values.extend(
            [
                (
                    comparison_df["cmp_rate_m_yr"]
                    - comparison_df["cmp_std_m_yr"]
                ).to_numpy(),
                (
                    comparison_df["cmp_rate_m_yr"]
                    + comparison_df["cmp_std_m_yr"]
                ).to_numpy(),
            ]
        )

    for full_rate in rate_profiles.values():
        all_y_values.append(
            np.asarray(full_rate, dtype=float)[START_REAL_INDEX:END_REAL_INDEX]
        )

    finite_y = np.concatenate(
        [np.asarray(values, dtype=float).ravel() for values in all_y_values]
    )
    finite_y = finite_y[np.isfinite(finite_y)]

    if finite_y.size > 0:
        y_min = np.min(finite_y)
        y_max = np.max(finite_y)
        y_pad = 0.18 * (y_max - y_min) if y_max > y_min else 1.0
        ax.set_ylim(y_min - y_pad, y_max + y_pad)

    ax.legend(
        loc="best",
        fontsize=9,
        frameon=True,
        framealpha=0.92,
    )

    fig.text(
        0.5,
        -0.02,
        "CMP uses the 1992–2024 observational LRR. CASCADE uses the "
        f"{START_YEAR}–{END_YEAR} calibration simulation; signs are aligned so "
        "positive is accretion and negative is erosion.",
        ha="center",
        va="center",
        fontsize=9,
        style="italic",
        color="0.35",
    )

    final_run_dir = os.path.join(OUTPUT_BASE_DIR, run_name)
    os.makedirs(final_run_dir, exist_ok=True)

    figure_out = os.path.join(
        final_run_dir,
        "CMP_1992_2024_vs_CASCADE_1992_2007_shoreline_change_rate.png",
    )
    comparison_csv_out = os.path.join(
        final_run_dir,
        "CMP_1992_2024_vs_CASCADE_1992_2007_domain_comparison.csv",
    )
    statistics_csv_out = os.path.join(
        final_run_dir,
        "CMP_1992_2024_vs_CASCADE_1992_2007_statistics.csv",
    )

    fig.savefig(figure_out, dpi=300, bbox_inches="tight")
    comparison_df.to_csv(comparison_csv_out, index=False)
    pd.DataFrame(statistics_rows).to_csv(statistics_csv_out, index=False)

    print(f"\nSaved CMP vs CASCADE shoreline-rate plot: {figure_out}")
    print(f"Saved domain comparison table: {comparison_csv_out}")
    print(f"Saved calibration statistics: {statistics_csv_out}")

    plt.show()


if __name__ == "__main__":
    main()