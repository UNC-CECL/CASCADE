"""
Build CASCADE background_erosion array (m/yr) from DSAS total shoreline change.

- Input:  DSAS CSV with one row per GIS domain (1–90),
          containing total shoreline change over a period (m).
- Output: CSV with 120-domain background erosion rates (m/yr),
          aligned with CASCADE domain indexing:
          15 left buffers + 90 real domains + 15 right buffers.

Mapping:
    GIS domain 1–90  → CASCADE indices 15–104  (0-based)
    buffers (0–14, 105–119) are filled using nearest real-domain rate.
"""

import os
import numpy as np
import pandas as pd

# ============================================================
# USER SETTINGS – EDIT THESE FOR EACH PERIOD
# ============================================================

# 1) Clean DSAS file with total shoreline change per domain
DSAS_FILE = r"/data/hatteras_init/shoreline_change/dsas_1978_1997_CLEAN.csv"

# 2) Column names in that DSAS file
DOMAIN_COL = "domain_id"           # GIS domain ID field (1–90)
TOTAL_CHANGE_COL = "obs_total_change_m"  # total shoreline change over period (m)

# 3) Length of the time period for this DSAS dataset (years)
YEARS = 19  # 1978–1997 → 19 years; for 1997–2019 you’d set YEARS = 22

# 4) CASCADE domain configuration (should match your run script)
NUM_REAL_DOMAINS = 90
NUM_BUFFER_DOMAINS = 15
TOTAL_DOMAINS = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS
LEFT_BUFFER_END = NUM_BUFFER_DOMAINS - 1    # = 14, last left buffer index

# 5) Output file path for the 120-domain background erosion (m/yr)
OUT_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\background_erosion\background_erosion_1978_1997_120domains.csv"

# 6) Optional cleaning parameters (can leave as-is for now)
ZERO_SMALL_RATES = True
SMALL_RATE_THRESHOLD = 0.05  # m/yr → set |rate| < this to zero
CLIP_EXTREME_RATES = True
MAX_ABS_RATE = 5.0           # m/yr → cap |rate| at this value

# ============================================================
# MAIN LOGIC
# ============================================================

def main():
    # --- Load DSAS data ---
    if not os.path.exists(DSAS_FILE):
        raise FileNotFoundError(f"DSAS file not found: {DSAS_FILE}")

    dsas = pd.read_csv(DSAS_FILE)

    if DOMAIN_COL not in dsas.columns or TOTAL_CHANGE_COL not in dsas.columns:
        raise ValueError(
            f"Expected columns '{DOMAIN_COL}' and '{TOTAL_CHANGE_COL}' "
            f"in {DSAS_FILE}, but got: {list(dsas.columns)}"
        )

    print(f"Loaded DSAS file: {DSAS_FILE}")
    print(f"Rows: {len(dsas)}")
    print(dsas[[DOMAIN_COL, TOTAL_CHANGE_COL]].head())

    # --- Compute rate (m/yr) from total change (m) over YEARS ---
    dsas["rate_m_per_yr"] = dsas[TOTAL_CHANGE_COL] / YEARS

    # --- Initialize full 120-domain array with zeros ---
    background_erosion = np.zeros(TOTAL_DOMAINS, dtype=float)

    # --- Map GIS 1–90 → CASCADE indices 15–104 ---
    for _, row in dsas.iterrows():
        gis_dom = int(row[DOMAIN_COL])            # 1–90
        rate = float(row["rate_m_per_yr"])        # m/yr

        cascade_idx = gis_dom + LEFT_BUFFER_END   # 1→15, ..., 90→104

        if 0 <= cascade_idx < TOTAL_DOMAINS:
            background_erosion[cascade_idx] = rate
        else:
            print(f"⚠️ Skipping GIS domain {gis_dom}: mapped index {cascade_idx} out of range")

    # --- Extend rates into buffers (so pattern is continuous) ---
    # Left buffers 0–14 get the first real-domain rate (index 15)
    first_real_idx = NUM_BUFFER_DOMAINS              # 15
    last_real_idx = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS - 1  # 104

    first_real_rate = background_erosion[first_real_idx]
    last_real_rate = background_erosion[last_real_idx]

    background_erosion[0:NUM_BUFFER_DOMAINS] = first_real_rate
    background_erosion[NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS:] = last_real_rate

    # --- Optional cleaning: zero tiny rates ---
    if ZERO_SMALL_RATES:
        small_mask = np.abs(background_erosion) < SMALL_RATE_THRESHOLD
        background_erosion[small_mask] = 0.0

    # --- Optional cleaning: clip extreme rates ---
    if CLIP_EXTREME_RATES:
        background_erosion = np.clip(background_erosion, -MAX_ABS_RATE, MAX_ABS_RATE)

    # --- Build output dataframe ---
    out_df = pd.DataFrame({
        "cascade_index": np.arange(TOTAL_DOMAINS),
        "background_erosion_m_per_yr": background_erosion
    })

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    out_df.to_csv(OUT_CSV, index=False)

    print("\nSaved background erosion array to:")
    print("  ", OUT_CSV)
    print("\nPreview:")
    print(out_df.head(20))
    print("\nStats:")
    print("  Min rate:", background_erosion.min(), "m/yr")
    print("  Max rate:", background_erosion.max(), "m/yr")


if __name__ == "__main__":
    main()
