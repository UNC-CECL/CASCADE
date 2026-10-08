"""
Per-domain distance to the offshore datum for any set of shoreline geojsons, with two figures to check it by eye.

    python HAT_geometric_distance_sanity_check.py

Reads the domain polygons, the datum line and SHORELINE_FILES from INPUT_DIR;
writes the raw-distance CSV, one change-from-reference CSV per REFERENCE_KEYS
and two sanity-check figures to OUTPUT_DIR. Needs numpy, pandas, matplotlib,
pyproj, shapely. Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""

import os
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pyproj import Transformer
from shapely.geometry import shape
from shapely.ops import transform as shp_transform

# --- CONFIG ------------------------------------------------------------------
REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
GROIN = REPO / "hard-structures" / "groin"
INPUT_DIR  = str(GROIN / "3-hindcast" / "1-dipole-1967-2017" / "inputs" / "groin_init" / "island_offset" / "input")
OUTPUT_DIR = str(GROIN / "1-observations" / "wetdry_photo_positions")

TARGET_CRS = "EPSG:26918"

DOMAIN_POLYGONS_FILE = os.path.join(INPUT_DIR, "domains_subset_2_12.geojson")
DATUM_LINE_FILE       = os.path.join(INPUT_DIR, "offshore_datum_line.geojson")

# Any mix of single-feature and multi-feature (dated) shoreline files
SHORELINE_FILES = [
    os.path.join(INPUT_DIR, "duneline_1967.geojson"),
    os.path.join(INPUT_DIR, "duneline_1978.geojson"),
    os.path.join(INPUT_DIR, "duneline_1997.geojson"),
    os.path.join(INPUT_DIR, "duneline_2017.geojson"),
    os.path.join(INPUT_DIR, "wet_dry_shorelines_groin.geojson"),
]

N_SAMPLES = 25   # points sampled along each domain's clipped shoreline segment

REFERENCE_KEYS = ["duneline_1967", "wetdry_1967"]   # one change CSV per key; [] = raw CSV only (README)

FIGURE_REFERENCE_KEY = "wetdry_1967"   # the profile figure's reference, separate from REFERENCE_KEYS

# Known-good min-subtracted series to validate against; KEY None skips it (README)
VALIDATION_KEY = "duneline_1967"
VALIDATION_KNOWN_GOOD = {
    2: 895.0, 3: 772.2, 4: 603.2, 5: 471.6, 6: 368.0,
    7: 300.6, 8: 233.6, 9: 166.4, 10: 102.4, 11: 44.6, 12: 0.0,
}
VALIDATION_TOLERANCE_M = 15.0

# Pairwise checks: (should_be_seaward, reference), smaller distance at every shared domain
SEAWARD_CHECK_PAIRS = [
    ("wetdry_1967", "duneline_1967"),
]

FIGURE_DPI = 200
# -----------------------------------------------------------------------------


# A geojson and its CRS name (TARGET_CRS if it names none)
def _load_geojson(path):
    with open(path) as f:
        data = json.load(f)
    src_crs = data.get("crs", {}).get("properties", {}).get("name", TARGET_CRS)
    return data, src_crs


# A shapely geometry reprojected into target_crs
def _reproject(geom, src_crs, target_crs=TARGET_CRS):
    if src_crs == target_crs:
        return geom
    transformer = Transformer.from_crs(src_crs, target_crs, always_xy=True)
    return shp_transform(lambda x, y, z=None: transformer.transform(x, y), geom)


# Domain polygons keyed by domain_id
def load_domain_polygons(path):
    data, src_crs = _load_geojson(path)
    polys = {}
    for feat in data["features"]:
        geom = _reproject(shape(feat["geometry"]), src_crs)
        polys[feat["properties"]["domain_id"]] = geom
    return polys


# The first feature of a file, reprojected
def load_single_geom(path):
    data, src_crs = _load_geojson(path)
    return _reproject(shape(data["features"][0]["geometry"]), src_crs)


# Every shoreline in a file as {label: (geometry, sort_key)}, labelled from year or filename (README)
def load_shorelines(path):
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
            # Single-feature file, no year: label and sort key from the filename
            digits = "".join(c for c in stem if c.isdigit())
            sort_key = float(digits) if digits else 0.0
            out[stem] = (geom, sort_key)

    return out


# Mean distance to the datum line of the shoreline clipped to each domain
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


# Change-from-reference CSV for every shoreline, in the groin comparison's format
def save_change_table(all_raw, domains, reference_key, output_dir):
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


# Map in real coordinates: domains, datum line, every shoreline coloured by year
def fig_spatial_map(domain_polys, datum_line, shorelines, out_path):
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

    # Zoom to the domain polygons, not the much longer shoreline and datum line
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
    ax.set_title("Spatial sanity check — domains, datum line, shorelines "
                 "(zoomed to domain range)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=FIGURE_DPI, facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# Distance to datum (or change from a reference) per domain, coloured as the map
def fig_distance_profile(all_raw, domains, out_path, reference_key=None):
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


# Run: load, measure, validate, check seaward pairs, write tables and figures
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

    # Validate against the known-good series
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

    # Seaward/landward pairwise checks
    for seaward_key, ref_key in SEAWARD_CHECK_PAIRS:
        if seaward_key in all_raw and ref_key in all_raw:
            r1, _ = all_raw[seaward_key]
            r2, _ = all_raw[ref_key]
            shared = [d for d in domains if r1.get(d) is not None and r2.get(d) is not None]
            ok = all(r1[d] < r2[d] for d in shared)
            print(f"\n'{seaward_key}' seaward of '{ref_key}' at all {len(shared)} shared domains: {ok}")

    # Combined raw-distance table
    df = pd.DataFrame({"Domain_ID": domains})
    for label, (result, _) in all_raw.items():
        df[label] = [result.get(d) for d in domains]
    csv_out = os.path.join(OUTPUT_DIR, "geometric_distances_all_shorelines.csv")
    df.to_csv(csv_out, index=False)
    print(f"\nSaved: {csv_out}")

    # Change-from-reference tables, the plotting-ready output
    print("\nBuilding change-from-reference tables...")
    for ref_key in REFERENCE_KEYS:
        save_change_table(all_raw, domains, ref_key, OUTPUT_DIR)

    # Figures
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
