"""
HAT_geometric_distance_sanity_check.py
========================================
General-purpose tool: computes per-domain distance-to-datum for ANY set of
shoreline geojsons (dune lines, wet/dry lines, or a mix -- doesn't matter
which), using the offshore datum line + domain polygons, and produces two
diagnostic figures so the result can be checked BY EYE, not just trusted
numerically.

Consolidates and generalizes HAT_duneline_geometric_distance.py and
HAT_wetdry_geometric_distance.py into one reusable tool. Handles both file
styles automatically:
  - one feature per file, no 'year' property (e.g. duneline_1967.geojson)
    -> labeled from the filename
  - many dated features in one file (e.g. wet_dry_shorelines_groin.geojson)
    -> labeled from each feature's own 'year' (+ month, if a year repeats)

FIGURES
-------
1. Spatial map (real coordinates): domain polygons, the datum line, and
   every shoreline overlaid, colored by year (older=blue, newer=red via
   colormap) -- lets you SEE whether shorelines actually fall inside the
   expected domains and stay roughly parallel to the datum line. This is
   the figure most likely to catch a wrong CRS, wrong file, or reversed
   geometry immediately, before ever looking at a number.
2. Distance profile: distance-to-datum (or change-from-REFERENCE_KEY, if
   set) per domain, one line per shoreline, same color scheme as the map.

Also runs the same validation/sanity checks already established:
  - if VALIDATION_KNOWN_GOOD is set, checks a named shoreline's min-
    subtracted values against it (as done for 1967 dune line).
  - prints a pairwise seaward/landward check between any two named
    shorelines you list in SEAWARD_CHECK_PAIRS (e.g. wet/dry vs dune line
    at the same year -- wet/dry should be seaward everywhere).

Author: Hannah A. Henry, UNC CECL
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pyproj import Transformer
from shapely.geometry import shape
from shapely.ops import transform as shp_transform

# =============================================================================
# CONFIG
# =============================================================================
INPUT_DIR  = r"/hard-structures/groin\hindcast_groin_test\input_prep\shoreline_position\input"
OUTPUT_DIR = r"/scripts/groin/HAT-hindcast-groin-test/input_prep/shoreline_position/shoreline_position_output"

TARGET_CRS = "EPSG:26918"

DOMAIN_POLYGONS_FILE = os.path.join(INPUT_DIR, "domains_subset_2_12.geojson")
DATUM_LINE_FILE       = os.path.join(INPUT_DIR, "offshore_datum_line.geojson")

# Any mix of single-feature and multi-feature (dated) shoreline files.
SHORELINE_FILES = [
    os.path.join(INPUT_DIR, "duneline_1967.geojson"),
    os.path.join(INPUT_DIR, "duneline_1978.geojson"),
    os.path.join(INPUT_DIR, "duneline_1997.geojson"),
    os.path.join(INPUT_DIR, "duneline_2017.geojson"),
    os.path.join(INPUT_DIR, "wet_dry_shorelines_groin.geojson"),
]

N_SAMPLES = 25   # points sampled along each domain's clipped shoreline segment

# Change-from-reference CSVs are saved for each key in this list (one CSV per
# key), ready to feed straight into HAT_groin_effect_comparison.py's observed-
# change logic. Add both a dune-line and a wet/dry reference if you want to
# compare either baseline. Leave empty to skip (raw-distance CSV only).
REFERENCE_KEYS = ["duneline_1967", "wetdry_1967"]

# Which reference (if any) the *figure's* second panel shows change relative
# to -- independent from REFERENCE_KEYS above (that controls saved CSVs).
FIGURE_REFERENCE_KEY = "wetdry_1967"

# Optional validation against a known-good min-subtracted series (see
# HAT_duneline_geometric_distance.py). Set KEY to None to skip.
VALIDATION_KEY = "duneline_1967"
VALIDATION_KNOWN_GOOD = {
    2: 895.0, 3: 772.2, 4: 603.2, 5: 471.6, 6: 368.0,
    7: 300.6, 8: 233.6, 9: 166.4, 10: 102.4, 11: 44.6, 12: 0.0,
}
VALIDATION_TOLERANCE_M = 15.0

# Optional pairwise seaward/landward checks: (should_be_seaward, reference).
# Prints whether the first is seaward of (smaller distance than) the second
# at every domain both have data for.
SEAWARD_CHECK_PAIRS = [
    ("wetdry_1967", "duneline_1967"),
]

FIGURE_DPI = 200

# =============================================================================
# GEOMETRY HELPERS
# =============================================================================
def _load_geojson(path):
    with open(path) as f:
        data = json.load(f)
    src_crs = data.get("crs", {}).get("properties", {}).get("name", TARGET_CRS)
    return data, src_crs


def _reproject(geom, src_crs, target_crs=TARGET_CRS):
    if src_crs == target_crs:
        return geom
    transformer = Transformer.from_crs(src_crs, target_crs, always_xy=True)
    return shp_transform(lambda x, y, z=None: transformer.transform(x, y), geom)


def load_domain_polygons(path):
    data, src_crs = _load_geojson(path)
    polys = {}
    for feat in data["features"]:
        geom = _reproject(shape(feat["geometry"]), src_crs)
        polys[feat["properties"]["domain_id"]] = geom
    return polys


def load_single_geom(path):
    data, src_crs = _load_geojson(path)
    return _reproject(shape(data["features"][0]["geometry"]), src_crs)


def load_shorelines(path):
    """Load every shoreline feature in a file, auto-labeling each one.

    - Multi-feature files with a 'year' property (e.g. wet/dry): labeled
      'wetdry_<year>' or 'wetdry_<year>_<month>' if the year repeats.
    - Single-feature files with no 'year' property (e.g. duneline_YYYY.geojson):
      labeled from the filename stem.

    Returns dict {label: (geometry, sort_key)}; sort_key is a float year
    (with a small month-based fraction added for duplicate years) so
    figures can be colored/ordered chronologically regardless of label text.
    """
    data, src_crs = _load_geojson(path)
    stem = os.path.splitext(os.path.basename(path))[0]
    out = {}
    seen_years = {}
    is_multi_feature = len(data["features"]) > 1

    month_frac = {
        "january": 0.04, "february": 0.12, "march": 0.20, "april": 0.29,
        "may": 0.37, "june": 0.45, "july": 0.54, "august": 0.62,
        "september": 0.70, "october": 0.79, "november": 0.87, "december": 0.95,
    }

    for feat in data["features"]:
        props = feat.get("properties", {})
        year = props.get("year") if is_multi_feature else None
        geom = _reproject(shape(feat["geometry"]), src_crs)

        if year is not None:
            prefix = "wetdry" if "wet" in stem.lower() or "dry" in stem.lower() else stem
            key = f"{prefix}_{year}"
            month = str(props.get("month", "")).lower()
            sort_key = float(year) + month_frac.get(month, 0.0)
            if key in out:
                existing_geom, existing_sort = out.pop(key)
                out[f"{key}_{list(seen_years.get(key, ['prior']))[-1]}"] = (existing_geom, existing_sort)
                key = f"{key}_{month}" if month else f"{key}_2"
            seen_years.setdefault(key, []).append(month)
            out[key] = (geom, sort_key)
        else:
            # Single-feature file, no year -- derive label + sort key from filename.
            digits = "".join(c for c in stem if c.isdigit())
            sort_key = float(digits) if digits else 0.0
            out[stem] = (geom, sort_key)

    return out


def per_domain_mean_distance(line, domain_polys, datum_line, n_samples=N_SAMPLES):
    results = {}
    for d, poly in domain_polys.items():
        clipped = line.intersection(poly)
        if clipped.is_empty:
            results[d] = None
            continue
        if clipped.geom_type == "LineString":
            lines = [clipped]
        elif clipped.geom_type == "MultiLineString":
            lines = list(clipped.geoms)
        elif clipped.geom_type == "Point":
            results[d] = datum_line.distance(clipped)
            continue
        elif clipped.geom_type == "MultiPoint":
            results[d] = float(np.mean([datum_line.distance(p) for p in clipped.geoms]))
            continue
        else:
            lines = []
        dists = []
        for ln in lines:
            for i in range(n_samples):
                pt = ln.interpolate(i / (n_samples - 1), normalized=True)
                dists.append(datum_line.distance(pt))
        results[d] = float(np.mean(dists)) if dists else None
    return results


# =============================================================================
# FIGURES
# =============================================================================
def fig_spatial_map(domain_polys, datum_line, shorelines, out_path):
    """Real-coordinate map: domain polygons, datum line, every shoreline
    colored chronologically. The figure most likely to catch a CRS/file/
    geometry mistake by eye before trusting any number."""
    fig, ax = plt.subplots(figsize=(10, 12))

    for d, poly in domain_polys.items():
        xs, ys = poly.exterior.xy
        ax.fill(xs, ys, color="#dddddd", edgecolor="#999999", alpha=0.5, zorder=1)
        cx, cy = poly.centroid.x, poly.centroid.y
        ax.text(cx, cy, f"D{d}", fontsize=8, ha="center", va="center",
                color="#555555", zorder=2)

    dx, dy = datum_line.xy
    ax.plot(dx, dy, color="black", ls="--", lw=1.5, label="Offshore datum line", zorder=3)

    sort_keys = [v[1] for v in shorelines.values()]
    vmin, vmax = min(sort_keys), max(sort_keys)
    norm = plt.Normalize(vmin=vmin, vmax=vmax if vmax > vmin else vmin + 1)
    cmap = cm.get_cmap("coolwarm")

    for label, (geom, sort_key) in sorted(shorelines.items(), key=lambda kv: kv[1][1]):
        color = cmap(norm(sort_key))
        lines = [geom] if geom.geom_type == "LineString" else list(geom.geoms)
        for ln in lines:
            xs, ys = ln.xy
            ax.plot(xs, ys, color=color, lw=1.3, alpha=0.85, zorder=4)

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Year (approx.)", fontsize=9)

    # Zoom to the domain polygons' own extent, not the full (much longer)
    # shoreline/datum-line length -- that's the point of this figure.
    all_bounds = [poly.bounds for poly in domain_polys.values()]
    xmin = min(b[0] for b in all_bounds)
    ymin = min(b[1] for b in all_bounds)
    xmax = max(b[2] for b in all_bounds)
    ymax = max(b[3] for b in all_bounds)
    pad_x = 0.10 * (xmax - xmin)
    pad_y = 0.05 * (ymax - ymin)
    ax.set_xlim(xmin - pad_x, xmax + pad_x)
    ax.set_ylim(ymin - pad_y, ymax + pad_y)

    ax.set_aspect("equal")
    ax.set_xlabel("Easting (m, EPSG:26918)")
    ax.set_ylabel("Northing (m, EPSG:26918)")
    ax.set_title("Spatial sanity check \u2014 domains, datum line, shorelines "
                 "(zoomed to domain range)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=FIGURE_DPI, facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def fig_distance_profile(all_raw, domains, out_path, reference_key=None):
    """Distance-to-datum (or change-from-reference) per domain, one line
    per shoreline, colored chronologically to match fig_spatial_map."""
    sort_keys = {k: v[1] for k, v in all_raw.items()}
    vmin, vmax = min(sort_keys.values()), max(sort_keys.values())
    norm = plt.Normalize(vmin=vmin, vmax=vmax if vmax > vmin else vmin + 1)
    cmap = cm.get_cmap("coolwarm")

    fig, ax = plt.subplots(figsize=(11, 6))

    ref_vals = all_raw[reference_key][0] if reference_key in all_raw else None

    for label, (result, sort_key) in sorted(all_raw.items(), key=lambda kv: kv[1][1]):
        vals = [result.get(d) for d in domains]
        if ref_vals is not None:
            vals = [
                (v - ref_vals[d]) if (v is not None and ref_vals.get(d) is not None) else np.nan
                for d, v in zip(domains, vals)
            ]
        color = cmap(norm(sort_key))
        ax.plot(domains, vals, "o-", color=color, lw=1.4, ms=4, alpha=0.85, label=label)

    ax.set_xlabel("GIS Domain ID")
    ylabel = (f"Change from {reference_key} (m)" if reference_key in all_raw
              else "Distance to datum (m)")
    ax.set_ylabel(ylabel)
    ax.set_title("Distance profile sanity check", fontsize=12, fontweight="bold")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6.5, loc="center left", bbox_to_anchor=(1.01, 0.5), ncol=1)
    fig.tight_layout()
    fig.savefig(out_path, dpi=FIGURE_DPI, facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def save_change_table(all_raw, domains, reference_key, output_dir):
    """Change-from-reference_key CSV for every shoreline, ready to feed
    HAT_groin_effect_comparison.py's observed-change logic directly. Mirrors
    HAT_duneline_geometric_distance.py / HAT_wetdry_geometric_distance.py's
    combined-change-table format."""
    if reference_key not in all_raw:
        print(f"  REFERENCE_KEYS: '{reference_key}' not found among loaded "
              f"shorelines -- skipping its change table.")
        return None

    ref_vals, _ = all_raw[reference_key]
    df = pd.DataFrame({"Domain_ID": domains})
    for label, (result, _) in sorted(all_raw.items(), key=lambda kv: kv[1][1]):
        change = [
            (result[d] - ref_vals[d])
            if (result.get(d) is not None and ref_vals.get(d) is not None)
            else np.nan
            for d in domains
        ]
        df[f"change_from_{reference_key}_{label}_m"] = change

    out_path = os.path.join(output_dir, f"Change_from_{reference_key}_D2_D12.csv")
    df.to_csv(out_path, index=False)
    print(f"  Saved: {out_path}")
    return df


# =============================================================================
# MAIN
# =============================================================================
def main():
    print("=" * 78)
    print("GEOMETRIC DISTANCE + SANITY-CHECK FIGURES")
    print("=" * 78)

    for label, path in [("DOMAIN_POLYGONS_FILE", DOMAIN_POLYGONS_FILE),
                        ("DATUM_LINE_FILE", DATUM_LINE_FILE)] + \
                       [(f"SHORELINE_FILES[{i}]", p) for i, p in enumerate(SHORELINE_FILES)]:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Missing {label}: {os.path.abspath(path)}")
    print("  All input files found.")

    domain_polys = load_domain_polygons(DOMAIN_POLYGONS_FILE)
    datum_line   = load_single_geom(DATUM_LINE_FILE)
    domains = sorted(domain_polys.keys())
    print(f"  {len(domain_polys)} domain polygons, datum line length={datum_line.length:.0f} m")

    shorelines = {}
    for path in SHORELINE_FILES:
        loaded = load_shorelines(path)
        shorelines.update(loaded)
        print(f"  {os.path.basename(path)}: {len(loaded)} shoreline(s) -> "
              f"{sorted(loaded.keys())}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("\nComputing per-domain distances...")
    all_raw = {}
    for label, (geom, sort_key) in shorelines.items():
        result = per_domain_mean_distance(geom, domain_polys, datum_line)
        n_ok = sum(1 for v in result.values() if v is not None)
        print(f"  {label}: {n_ok}/{len(domains)} domains matched")
        all_raw[label] = (result, sort_key)

    # --- Validation check ---
    if VALIDATION_KEY and VALIDATION_KEY in all_raw and VALIDATION_KNOWN_GOOD:
        print(f"\nValidating '{VALIDATION_KEY}' against known-good values...")
        result, _ = all_raw[VALIDATION_KEY]
        vd = sorted(VALIDATION_KNOWN_GOOD.keys())
        raw = np.array([result[d] for d in vd])
        rel = raw - np.nanmin(raw)
        known = np.array([VALIDATION_KNOWN_GOOD[d] for d in vd])
        max_diff = np.max(np.abs(rel - known))
        print(f"  Max difference: {max_diff:.1f} m (tolerance {VALIDATION_TOLERANCE_M} m)  "
              f"-- {'OK' if max_diff <= VALIDATION_TOLERANCE_M else 'CHECK THIS'}")
        print(f"  Correlation: {np.corrcoef(rel, known)[0, 1]:.6f}")

    # --- Seaward/landward pairwise checks ---
    for seaward_key, ref_key in SEAWARD_CHECK_PAIRS:
        if seaward_key in all_raw and ref_key in all_raw:
            r1, _ = all_raw[seaward_key]
            r2, _ = all_raw[ref_key]
            shared = [d for d in domains if r1.get(d) is not None and r2.get(d) is not None]
            ok = all(r1[d] < r2[d] for d in shared)
            print(f"\n'{seaward_key}' seaward of '{ref_key}' at all {len(shared)} shared domains: {ok}")

    # --- Save combined raw-distance table ---
    df = pd.DataFrame({"Domain_ID": domains})
    for label, (result, _) in all_raw.items():
        df[label] = [result.get(d) for d in domains]
    csv_out = os.path.join(OUTPUT_DIR, "geometric_distances_all_shorelines.csv")
    df.to_csv(csv_out, index=False)
    print(f"\nSaved: {csv_out}")

    # --- Save change-from-reference tables (the actual plotting-ready shoreline_position_output) ---
    print("\nBuilding change-from-reference tables...")
    for ref_key in REFERENCE_KEYS:
        save_change_table(all_raw, domains, ref_key, OUTPUT_DIR)

    # --- Figures ---
    print("\nBuilding sanity-check figures...")
    fig_spatial_map(domain_polys, datum_line,
                    {k: v for k, v in shorelines.items()},
                    os.path.join(OUTPUT_DIR, "sanity_check_spatial_map.png"))
    fig_distance_profile(all_raw, domains,
                          os.path.join(OUTPUT_DIR, "sanity_check_distance_profile.png"),
                          reference_key=FIGURE_REFERENCE_KEY)

    print("\nDone.")


if __name__ == "__main__":
    main()
