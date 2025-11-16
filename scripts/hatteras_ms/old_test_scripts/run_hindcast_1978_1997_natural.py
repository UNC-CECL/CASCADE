# ==============================================================================
# HINDCAST SIMULATION SCRIPT – HATTERAS ISLAND (1978–1997) — NATURAL ONLY
# ==============================================================================

"""
Natural-processes-only CASCADE hindcast for Hatteras Island:
- No roadway effects (setbacks/elevation), no nourishment, no dune rebuilds, no sandbags.
- Alongshore transport ON, constant SLR as configured.
- Buffers included on both ends for stable alongshore fluxes.

Author: Hannah Henry
Affiliation: University of North Carolina at Chapel Hill
Date Created: 2024
Last Modified: 2025 (natural-only cleanup)
"""

# ------------------------------
# Imports
# ------------------------------
import os
import copy
import logging
import numpy as np
import pandas as pd
from cascade.cascade import Cascade  # CASCADE core model class

# ------------------------------
# Run identity & logging
# ------------------------------
run_name = 'Hindcast_1978_1997_NATURAL'
logging.basicConfig(
    filename='../cascade_hindcast.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logging.info(f"Starting simulation run: {run_name}")

# ------------------------------
# Working directory (adjust if needed)
# ------------------------------
os.chdir('/')

# ------------------------------
# Core inputs
# ------------------------------
s_file = r'/data/hatteras_init/storms/hindcast_storms/HAT_1978_2022_Final_Hindcast_Storms.npy'

# Model geometry
REAL_DOMAINS = 92
LEFT_BUFFERS = 15
RIGHT_BUFFERS = 15
Total_B3D_Number = LEFT_BUFFERS + REAL_DOMAINS + RIGHT_BUFFERS  # = 122

# Dune offsets (shoreline offset feature) — choose one start year
start_year = 1978  # natural hindcast start
dune_load_name = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\1978\Relative_Dune_1978_Offset.csv'
dune_df = pd.read_csv(dune_load_name)
# Expect a column named "1978" in meters; CASCADE shoreline_offset is commonly in decimeters
dune_offset_real = dune_df[str(start_year)].to_numpy() * 10.0  # m → dm (domain-length = 92)

# ------------------------------
# Time configuration
# ------------------------------
# 1978 → 1997 inclusive = 20 time steps (years).
run_years = 20

# ------------------------------
# Background shoreline change (m/yr)
# ------------------------------
# Your original list (length 92) was all zeros; keep it but make it explicit and pad for buffers.
background_erosion_real = [0.0] * REAL_DOMAINS
background_erosion = [0.0] * LEFT_BUFFERS + background_erosion_real + [0.0] * RIGHT_BUFFERS  # length 122

# ------------------------------
# Natural-only road configuration (all zero / False)
# ------------------------------
road_setbacks = np.zeros(Total_B3D_Number, dtype=float)   # cm if used by internals; keeping zeros safe
road_ele = [0.0] * Total_B3D_Number                       # m
road_cells = [False] * Total_B3D_Number                   # not used, but explicit

# ------------------------------
# Elevation & dune files (with buffers)
# ------------------------------
e_file = []  # elevation files (length 122)
d_file = []  # dune files (length 122)

# Left buffers (placeholder arrays)
for _ in range(LEFT_BUFFERS):
    d_file.append(r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\buffer\Sample_1_dune.npy')
    e_file.append(r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\buffer\Sample_1_topography.npy')

# Real domains (1..92)
for i in range(1, REAL_DOMAINS + 1):
    d_file.append(fr'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\dunes\domain_{i}_resampled_dune.npy')
    e_file.append(fr'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\topography\domain_{i}_resampled_topography.npy')

# Right buffers (placeholder arrays)
for _ in range(RIGHT_BUFFERS):
    d_file.append(r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\buffer\Sample_1_dune.npy')
    e_file.append(r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\buffer\Sample_1_topography.npy')

# ------------------------------
# Helper: run alongshore-connected model
# ------------------------------
def alongshore_connected(
    nt,
    name,
    storm_file,
    alongshore_section_count,  # should be 122 when using buffers
    num_cores,
    # Physical parameters / thresholds (natural-only values where applicable)
    rmin,                      # min dune growth rate per REAL domain (len 92) — CASCADE param meaning
    rmax,                      # max dune growth rate per REAL domain (len 92)
    elevation_file,            # len 122
    dune_file,                 # len 122
    dune_design_elevation,     # len 92
    dune_minimum_elevation,    # scalar or len 92
    # Roads & management (off / zeroed)
    road_ele,                  # len 122 (zeros)
    road_width,                # scalar (unused when roads off)
    road_setback,              # len 122 (zeros)
    # Natural-only: no nourishment / rebuilds
    background_erosion,        # len 122
    # Sea level
    sea_level_rise_rate=0.0056,
    sea_level_constant=True,
    # Shoreline offset
    enable_shoreline_offset=True,
    shoreline_offset_dm=None   # len 92 in decimeters if enabled
):
    """
    Natural-only wrapper that instantiates and runs CASCADE with management modules OFF.
    """
    datadir = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init'

    cascade = Cascade(
        datadir=datadir,
        name=name,
        # Core storm & terrain
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,

        # Physical parameter file (keep yours)
        parameter_file="Hatteras-CASCADE-parameters.yaml",

        # Offshore wave climate (as you had; adjust as needed)
        wave_height=1.0,
        wave_period=7.0,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,

        # Sea level
        sea_level_rise_rate=sea_level_rise_rate,
        sea_level_rise_constant=sea_level_constant,

        # Background shoreline change
        background_erosion=background_erosion,

        # Grid & time
        alongshore_section_count=alongshore_section_count,  # 122 with buffers
        time_step_count=nt,
        min_dune_growth_rate=rmin,   # len 92 (applies within real domains)
        max_dune_growth_rate=rmax,   # len 92
        num_cores=num_cores,

        # MODULE TOGGLES — NATURAL ONLY
        roadway_management_module=False,     # roads OFF
        alongshore_transport_module=True,    # keep transport ON
        beach_nourishment_module=False,      # nourishment OFF
        community_economics_module=False,    # OFF

        # Road config (ignored when module is off, but pass zeros to be safe)
        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,

        # Dune thresholds (used physically, not management)
        dune_design_elevation=dune_design_elevation,     # len 92
        dune_minimum_elevation=dune_minimum_elevation,   # scalar or len 92

        # Sandbags OFF
        sandbag_management_on=[False] * alongshore_section_count,
        sandbag_elevation=False,

        # Shoreline offset
        enable_shoreline_offset=enable_shoreline_offset,
        shoreline_offset=shoreline_offset_dm if shoreline_offset_dm is not None else [0]*REAL_DOMAINS
    )

    # ---- Run loop ----
    for _ in range(nt - 1):
        cascade.update()
        if cascade.b3d_break:
            break

    # Save outputs
    base_output_dir = "output"
    run_output_dir = os.path.join(base_output_dir, name)
    os.makedirs(run_output_dir, exist_ok=True)
    cascade.save(run_output_dir)
    logging.info("Natural-only simulation completed successfully.")

    return cascade

# ==============================================================================
# Configure & Execute — Hatteras (92 real + 30 buffer = 122)
# ==============================================================================
def alongshore_uniform_natural():
    """
    Natural-only, alongshore-connected simulation with buffers.
    """
    number_barrier3d_models = REAL_DOMAINS  # for 92-length parameter arrays
    alongshore_with_buffers = Total_B3D_Number  # 122

    # Domain-wise phys parameters (len 92)
    beach_width_threshold = [30.0] * number_barrier3d_models  # retained param; not used by manager
    rmin = [0.55] * number_barrier3d_models
    rmax = [0.95] * number_barrier3d_models
    dune_design_elevation = [1.5] * number_barrier3d_models
    dune_minimum_elevation = 1.0  # scalar OK

    # Shoreline offset (enable if you want to respect your 1978 offsets)
    shoreline_offset_enabled = True
    shoreline_offset_dm = dune_offset_real  # len 92, in decimeters

    # Compute a 122-length shoreline offset array by padding buffers with zeros (applies only to real domains internally)
    # CASCADE expects len 92 here (per-domain), so we pass the 92-length array directly above.

    # NATURAL-ONLY: management OFF, but pass zero-valued arrays
    num_cores = 4
    road_width = 20.0  # unused when roads OFF
    nourishment_volume = 0.0  # not used
    overwash_to_dune = 0.0    # not used

    sea_level_rise_rate = 0.0056
    sea_level_constant = True

    alongshore_connected(
        nt=run_years,
        name=run_name,
        storm_file=s_file,
        alongshore_section_count=alongshore_with_buffers,  # 122
        num_cores=num_cores,
        rmin=rmin,
        rmax=rmax,
        elevation_file=e_file,
        dune_file=d_file,
        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,
        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setbacks,
        background_erosion=background_erosion,
        sea_level_rise_rate=sea_level_rise_rate,
        sea_level_constant=sea_level_constant,
        enable_shoreline_offset=shoreline_offset_enabled,
        shoreline_offset_dm=shoreline_offset_dm
    )

# Execute
if __name__ == "__main__":
    alongshore_uniform_natural()
