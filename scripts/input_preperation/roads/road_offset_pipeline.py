"""
Hatteras CASCADE Road Setback Generator
=======================================

Purpose
-------
Compute CASCADE-ready road–dune setback distances for Hatteras Island.
This script reads raw road and dune offset CSVs (exported from ArcGIS),
computes the minimum offshore distance per domain for both road and dune,
calculates the road–dune setback in meters, and writes a single 2-row
CASCADE-ready CSV file.

Final Output
------------
CASCADE_RoadSetback_<YEAR>.csv

Row 1: Domain IDs (1–105)
Row 2: Road setback distances (meters)

Notes
-----
- CASCADE expects 105 real domains (1–105).
- No alongshore padding required for road arrays.
- Units are meters; convert to centimeters/decimeters externally if needed.

Author: Hannah A. Henry (refactored pipeline)
"""

import os
import unicodedata
import numpy as np
import pandas as pd

# =============================================================================
# 1. USER CONFIGURATION
# =============================================================================

YEAR = 1978

RAW_ROAD_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\roads\offset\1978\1978_road_offset_raw.csv"
RAW_DUNE_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\1978\1978_duneline_offset_raw.csv"

OUTPUT_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\roads\offset"
FINAL_OUTPUT_NAME = f"CASCADE_RoadSetback_{YEAR}.csv"

# CASCADE domain window
START_DOMAIN = 1
END_DOMAIN   = 105

# Column names from ArcGIS exports
ROAD_ID_COLUMN = "domain_id_1"
DUNE_ID_COLUMN = "domain_id"

ROAD_DIST_COLUMN = "ORIG_LEN"
DUNE_DIST_COLUMN = "ORIG_LEN"


# =============================================================================
# 2. FUNCTIONS
# =============================================================================

def load_and_standardize_csv(path):
    """Load CSV, clean column names, drop unnamed columns."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

    df = pd.read_csv(path)
    df.columns = [unicodedata.normalize("NFKD", col).strip() for col in df.columns]
    df = df.loc[:, ~df.columns.str.contains("^Unnamed")]
    return df


def compute_min_offsets(road_df, dune_df):
    """Compute (min road distance) – (min dune distance) per domain."""
    road_grouped = road_df.groupby(ROAD_ID_COLUMN)[ROAD_DIST_COLUMN].min().reset_index()
    road_grouped.rename(columns={ROAD_ID_COLUMN: "domain_id",
                                 ROAD_DIST_COLUMN: "road_offset"}, inplace=True)

    dune_grouped = dune_df.groupby(DUNE_ID_COLUMN)[DUNE_DIST_COLUMN].min().reset_index()
    dune_grouped.rename(columns={DUNE_ID_COLUMN: "domain_id",
                                 DUNE_DIST_COLUMN: "dune_offset"}, inplace=True)

    merged = pd.merge(road_grouped, dune_grouped, on="domain_id", how="inner")
    merged["Road_Setback_m"] = merged["road_offset"] - merged["dune_offset"]
    merged = merged.sort_values("domain_id").reset_index(drop=True)
    return merged


def filter_domain_range(df, start_id, end_id):
    """Keep only domains in CASCADE window (1–105)."""
    mask = (df["domain_id"] >= start_id) & (df["domain_id"] <= end_id)
    filtered = df.loc[mask].copy()
    filtered = filtered.sort_values("domain_id").reset_index(drop=True)
    return filtered


def save_two_row_cascade(df, output_path):
    """Write final 2-row CASCADE format."""
    df = df[["domain_id", "Road_Setback_m"]].apply(pd.to_numeric, errors="coerce").dropna()

    domain_ids = df["domain_id"].to_numpy()
    setbacks = df["Road_Setback_m"].to_numpy()

    array_out = np.vstack([domain_ids, setbacks])

    np.savetxt(
        output_path,
        array_out,
        delimiter=",",
        fmt="%.3f"
    )

    print(f"\nCASCADE road setback file saved:")
    print(f"  {output_path}")
    print(f"  Shape: {array_out.shape[0]} rows x {array_out.shape[1]} columns")


# =============================================================================
# 3. MAIN PIPELINE
# =============================================================================

def main():
    print(f"\n=== Generating CASCADE Road Setbacks for {YEAR} ===")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load raw datasets
    road_raw = load_and_standardize_csv(RAW_ROAD_CSV)
    dune_raw = load_and_standardize_csv(RAW_DUNE_CSV)

    # Compute setbacks
    offsets = compute_min_offsets(road_raw, dune_raw)

    # Filter to CASCADE domain window
    offsets_filtered = filter_domain_range(offsets, START_DOMAIN, END_DOMAIN)

    # Save final CASCADE 2-row file
    final_path = os.path.join(OUTPUT_DIR, FINAL_OUTPUT_NAME)
    save_two_row_cascade(offsets_filtered, final_path)

    print("\nPipeline complete.\n")


if __name__ == "__main__":
    main()
