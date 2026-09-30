"""
Raw per-transect offsets -> the relative island offset per domain, padded to 120 domains for the model.

    python scripts/input_prep/2-brie-offset/1-produce/island_offset_hybrid.py --year 1996 --version v1

Relative offset per domain (metres, baseline = minimum), buffer domains
closed with a Hermite curve, written unpadded and padded with a diagnostic
figure under 2-brie-offset/<year>/. Details: scripts/input_prep/2-brie-offset/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 1. user configuration

# Paths anchored on the repo root and YEAR, so either start can be built
import argparse as _argparse
import sys as _sys

_PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                     if (_p / "pyproject.toml").exists())
_sys.path.insert(0, str(_PROJECT_ROOT / "scripts"))
from site_layer.hat_topo_version import (DUNE_LINE_FOR_YEAR, DEFAULT_OFFSET_SOURCE,  # noqa: E402
                                         OFFSET_SOURCES, dune_raw_file_for_year,
                                         offset_basename, offset_file,
                                         offset_start_dir,
                                         shoreline_raw_file_for_year)
from site_layer.hat_extension_domains import BASE_GEOMETRY, GEOMETRIES, SURVEYED_GIS, gis_bounds  # noqa: E402
from site_layer.hat_topo_version import offset_version  # noqa: E402
from cascade_pipeline.hindcast import pad_offset_ring  # noqa: E402

from site_layer.hat_topo_version import BRIE_ROOT as _BRIE_ROOT  # noqa: E402

_ap = _argparse.ArgumentParser(description="island offsets for one hindcast start")
# Any start year DUNE_LINE_FOR_YEAR pairs with a dune-line vintage
_ap.add_argument("--year", type=int, default=2004,
                 choices=tuple(sorted(DUNE_LINE_FOR_YEAR)))
# Versioned output: --version writes <year>/<version>/; CURRENT picks the one read
_ap.add_argument("--version", default=None,
                 help="write to <year>/<version>/ instead of <year>/ (e.g. v2)")
# The raw file to build from; defaults to the vintage's current one
_ap.add_argument("--raw-file", default=None,
                 help="raw CSV to read instead of the vintage's file in raw_offsets/")
# EXTENDED GEOMETRY (2026-09-16, the Pea Island extension experiment)
_ap.add_argument("--geometry", default=None,
                 choices=[g for g in GEOMETRIES if g != BASE_GEOMETRY],
                 help="an extended reach, to <year>/ext/<geometry>/")
# WHICH FEATURE THE OFFSET IS MEASURED FROM (2026-09-22)
_ap.add_argument("--source", default=DEFAULT_OFFSET_SOURCE, choices=OFFSET_SOURCES,
                 help="the feature the offset is measured from "
                      f"(default {DEFAULT_OFFSET_SOURCE})")
_args = _ap.parse_args()
YEAR = _args.year
VERSION = _args.version
GEOMETRY = _args.geometry
SOURCE = _args.source
if GEOMETRY and VERSION:
    _ap.error("--geometry builds are filed under ext/, not as a version")
if GEOMETRY and SOURCE != DEFAULT_OFFSET_SOURCE:
    _ap.error("the extended geometries are only built from the dune line")

# A source names its own raw file
RAW_FILE = (str(Path(_args.raw_file).resolve()) if _args.raw_file
            else str(dune_raw_file_for_year(YEAR)) if SOURCE == DEFAULT_OFFSET_SOURCE
            else str(shoreline_raw_file_for_year(YEAR)))
RAW_EXT_FILE = (str(Path(RAW_FILE).parent / "ext"
                    / (Path(RAW_FILE).stem + "_ext.csv")) if GEOMETRY else None)

_START_DIR = offset_start_dir(YEAR, SOURCE)
# ext/ sits under the source too (2026-09-22), and is checked against its CURRENT
OUTPUT_DIR    = str(_START_DIR / "ext" / GEOMETRY if GEOMETRY
                    else _START_DIR / VERSION if VERSION
                    else _START_DIR)
OUTPUT_BASENAME = offset_basename(YEAR, SOURCE)
# --- CONFIG ------------------------------------------------------------------
# Feature names for titles and warnings, so a shoreline build is not called 'Dune'
FEATURE_LABEL = {"duneline": "Dune", "shoreline": "Shoreline"}[SOURCE]
FEATURE_NOUN = {"duneline": "dune line", "shoreline": "mean shoreline"}[SOURCE]

START_DOMAIN, END_DOMAIN = gis_bounds(GEOMETRY or BASE_GEOMETRY)
B3D_GRIDS    = list(range(START_DOMAIN, END_DOMAIN + 1))

PADDING_ZEROS = 15
TARGET_LENGTH = (END_DOMAIN - START_DOMAIN + 1) + 2 * PADDING_ZEROS  # 120 for GIS 1-90

DOMAIN_SPACING_M = 500.0          # BRIE's dy; hatteras_site_config.HATTERAS_DOMAINS
UNSTABLE_ANGLE_DEG = 42.0         # (1.2 sin^2 - cos^2) changes sign here

COL_MAP = {
    "Domain_ID": "domain_id",
    "Distance":  "ORIG_LEN",
    "Transect":  "LineID",
}
# -----------------------------------------------------------------------------

# 2. functions

# Compute mean relative dune raw_offset per domain from raw CSV
def calculate_relative_offset(file_path, year, col_map, grids):
    print(f"\n--- Processing {year} ---")
    print(f"Input file: {file_path}")

    try:
        raw_df = pd.concat([pd.read_csv(f) for f in
                            ([file_path] if isinstance(file_path, str) else file_path)],
                           ignore_index=True)
    except FileNotFoundError as e:
        print(f"ERROR: File not found: {e.filename}")
        return None

    for key, col in col_map.items():
        if col not in raw_df.columns:
            print(f"ERROR: Column '{col}' not found.")
            return None

    data_df = pd.DataFrame({
        "B3D_Grid": raw_df[col_map["Domain_ID"]],
        "Distance": raw_df[col_map["Distance"]],
        "Transect": raw_df[col_map["Transect"]],
    })

    mean_distances = []
    seen_domains   = []

    for grid_id in grids:
        subset = data_df[data_df["B3D_Grid"] == int(grid_id)]
        if subset.empty:
            print(f"  Warning: No data for domain {grid_id}.")
            continue

        try:
            min_t = int(subset["Transect"].min())
            max_t = int(subset["Transect"].max())
        except ValueError:
            print(f"  Warning: Invalid transect data for domain {grid_id}. Skipping.")
            continue

        distances = []
        for t_id in range(min_t, max_t + 1):
            t_vals = subset[subset["Transect"] == t_id]
            if not t_vals.empty:
                distances.append(t_vals.iloc[0]["Distance"])

        if not distances:
            continue

        mean_distances.append(float(np.mean(distances)))
        seen_domains.append(grid_id)

    if not mean_distances:
        print(f"FATAL: No valid domains processed for year {year}.")
        return None

    baseline         = min(mean_distances)
    relative_offsets = np.subtract(mean_distances, baseline)

    print(f"  {len(mean_distances)} domains processed.")
    print(f"  Baseline distance = {baseline:.3f} m (min mean).")

    return pd.DataFrame({"Domain_ID": seen_domains, str(year): relative_offsets})


# The angle BRIE reads between each padded domain and the next, wrapping from the last back to the ...
def shoreline_angles_deg(padded):
    return np.degrees(np.arctan2(np.diff(np.r_[padded, padded[0]]), DOMAIN_SPACING_M))


# The padded profile and the shoreline angle BRIE reads, one figure
def plot_buffer_diagnostic(padded, year, padding, output_dir, output_basename):
    from site_layer.hat_figure_style import (C, INK_MUTED, _title, apply_style,
                                             figsize, open_frame, record_caption,
                                             save)
    apply_style()
    n_real = len(padded) - 2 * padding
    x = np.arange(len(padded)) - padding + START_DOMAIN     # GIS numbering, buffers outside
    theta = shoreline_angles_deg(padded)
    real = slice(padding, padding + n_real)
    fig, (ax_o, ax_t) = plt.subplots(2, 1, figsize=figsize("double", height=5.0),
                                     sharex=True, constrained_layout=True)
    for ax in (ax_o, ax_t):
        for lo, hi in ((x[0] - 0.5, START_DOMAIN - 0.5), (END_DOMAIN + 0.5, x[-1] + 0.5)):
            ax.axvspan(lo, hi, color=C["BASE_FILL"], alpha=0.6, lw=0)
        ax.grid(axis="y")
        open_frame(ax)
    ax_o.plot(x, padded / 1000.0, color=C["INK"], lw=1.4)
    ax_o.plot(x[real], padded[real] / 1000.0, color=C["LATE"], lw=1.8)
    ax_o.set_ylabel("Cross-shore offset (km)")
    # No build version in the title: the folder carries it
    _title(ax_o, 0, f"Initial {FEATURE_NOUN} offset, {year}")
    # the wrap from the last domain back to the first is drawn at the right end
    ax_t.plot(x, theta, color=C["INK"], lw=1.2)
    for sign in (1, -1):
        ax_t.axhline(sign * UNSTABLE_ANGLE_DEG, color=C["ACCENT"], lw=0.9, ls="--")
    ax_t.axhline(0, color=INK_MUTED, lw=0.6)
    ax_t.set_ylabel("Shoreline angle (°)")
    ax_t.set_xlabel("Alongshore domain (GIS numbering)")
    _title(ax_t, 1, "Shoreline orientation between adjacent domains")
    path = Path(output_dir) / f"{output_basename}_buffer_diagnostic.png"
    save(fig, path, vector=False, dpi=200, close=True)
    record_caption(path, (
        f"The padded {FEATURE_NOUN} offset for the {year} start, in metres as "
        f"written to {output_basename}_PADDED_{len(padded)}.csv. (a) Offset along "
        f"the padded domains; the real reach GIS {START_DOMAIN}-{END_DOMAIN} in "
        f"blue, the {padding} buffer domains each side shaded. The buffers close "
        "BRIE's periodic domain from the last real domain back round to the "
        "first along a cubic Hermite matched to the island's end slopes "
        "(cascade_pipeline.hindcast.pad_offset_ring), the same array offset_mode "
        "'metres' hands Cascade. (b) The shoreline angle between each padded "
        "domain and the next (atan2 of the offset step over 500 m), the last "
        f"point being the wrap; dashed at +/-{UNSTABLE_ANGLE_DEG:g} degrees, past "
        "which BRIE's alongshore diffusivity changes sign."))
    print(f"\nDiagnostic figure saved to:\n  {path}")
    return path


# 3. main

# Run: relative offsets, the padded file, the diagnostic figure
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 3.1. Compute relative offsets
    result = calculate_relative_offset(
        file_path=[RAW_FILE, RAW_EXT_FILE] if GEOMETRY else RAW_FILE,
        year=YEAR,
        col_map=COL_MAP,
        grids=B3D_GRIDS,
    )

    if result is None:
        print("Offset calculation failed. Exiting.")
        return

    if GEOMETRY:
        # The extension must not move the surveyed reach
        _base = offset_file(YEAR, "unpadded", source=SOURCE)
        if not _base.is_file():
            raise SystemExit(f"no surveyed build to check against: {_base}")
        base = pd.read_csv(_base).set_index("Domain_ID")[str(YEAR)]
        mine = result.set_index("Domain_ID")[str(YEAR)]
        lo, hi = SURVEYED_GIS
        missing = [g for g in B3D_GRIDS if g not in mine.index]
        diff = (mine.loc[lo:hi] - base.loc[lo:hi]).abs().max()
        print(f"\nGeometry {GEOMETRY}: GIS {START_DOMAIN}..{END_DOMAIN}, "
              f"{len(mine)} domains, {len(missing)} without a {FEATURE_NOUN} {missing}")
        _ver = offset_version(YEAR, SOURCE) or "(flat)"
        print(f"  surveyed slice vs {YEAR}/{_ver}: max |diff| {diff:.6f} m")
        if missing or diff > 1e-6:
            raise SystemExit("extension build changed or lost surveyed domains; refusing")

    # Save unpadded file with Domain_ID
    unpadded_path = os.path.join(OUTPUT_DIR, f"{OUTPUT_BASENAME}_CASCADE_Input_unpadded.csv")
    result.to_csv(unpadded_path, index=False)
    print(f"\nUnpadded offsets saved to:\n  {unpadded_path}")

    # CASCADE-format: raw_offset column only
    cascade_df = result[[str(YEAR)]].copy()
    cascade_unpadded_path = os.path.join(OUTPUT_DIR, f"{OUTPUT_BASENAME}_CASCADE_Input.csv")
    cascade_df.to_csv(cascade_unpadded_path, index=False)
    print(f"Unpadded CASCADE-format file saved to:\n  {cascade_unpadded_path}")

    # 3.2. Pad with the smooth wrap-around the model uses
    real_m = cascade_df[str(YEAR)].to_numpy(dtype=float)
    padded = pad_offset_ring(real_m, PADDING_ZEROS)
    if len(padded) != TARGET_LENGTH:
        raise SystemExit(f"padded length {len(padded)}, expected {TARGET_LENGTH}")
    theta = shoreline_angles_deg(padded)
    real = slice(PADDING_ZEROS, PADDING_ZEROS + len(real_m))
    buf = np.r_[theta[:PADDING_ZEROS], theta[PADDING_ZEROS + len(real_m) - 1:]]
    print("\nPadding summary (smooth wrap-around, pad_offset_ring):")
    print(f"  Buffer domains per side : {PADDING_ZEROS}")
    print(f"  Real span               : {real_m.min():.1f} - {real_m.max():.1f} m")
    print(f"  Buffer span             : {np.r_[padded[:PADDING_ZEROS], padded[-PADDING_ZEROS:]].min():.1f}"
          f" - {np.r_[padded[:PADDING_ZEROS], padded[-PADDING_ZEROS:]].max():.1f} m")
    print(f"  Largest angle, real     : {np.abs(theta[real][:-1]).max():.1f} deg")
    print(f"  Largest angle, buffer   : {np.abs(buf).max():.1f} deg "
          f"(BRIE goes anti-diffusive past ~{UNSTABLE_ANGLE_DEG:g})")
    if np.abs(theta).max() > UNSTABLE_ANGLE_DEG:
        print("  WARNING: an angle exceeds the anti-diffusive limit")

    padded_path = os.path.join(OUTPUT_DIR, f"{OUTPUT_BASENAME}_PADDED_{TARGET_LENGTH}.csv")
    pd.DataFrame({str(YEAR): padded}).to_csv(padded_path, index=False)
    print(f"\nSUCCESS: Padded CASCADE input saved to:\n  {padded_path}")

    # 3.3. Diagnostic figure
    plot_buffer_diagnostic(padded, YEAR, PADDING_ZEROS, OUTPUT_DIR, OUTPUT_BASENAME)


if __name__ == "__main__":
    main()
