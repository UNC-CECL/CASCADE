"""
HAT_trajectory_map_v1.py
======================
Produces a geographic map of Hatteras Island with CASCADE domains coloured
by shoreline trajectory classification (from domain_trajectory_metrics.csv).

Uses:
  - HAT_domains.json        : CASCADE domain polygons (EPSG:3725 / UTM)
  - CoastSat_transect_layer.geojson : transect lines (WGS84) for context
  - domain_trajectory_metrics.csv   : comparison from HAT_shoreline_trajectory_classification.py

Outputs
-------
  fig_map_trajectory.png   — plan-view map coloured by trajectory class
  fig_map_lrr.png          — plan-view map coloured by full-record LRR magnitude

Dependencies
------------
  pip install geopandas pyproj shapely matplotlib
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap, Normalize
import geopandas as gpd
from shapely.geometry import shape, LineString, Point
from pyproj import Transformer

# ============================================================
# CONFIG
# ============================================================

DOMAINS_GEOJSON   = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\HAT_domains.json"
TRANSECTS_GEOJSON = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\CoastSat_transect_layer.geojson"
METRICS_CSV       = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\shoreline_change_patterns\output\domain_trajectory_metrics.csv"
OUTLINE_SHP      = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\figure_making\map_elements\hatteras_outline\HAT_island_outline.shp"
OUTPUT_DIR        = r"/scripts/input_prep/shoreline_change_patterns/old_drafts/trajectory_map"

# Hatteras bounding box (WGS84) — used to filter transects
HAT_LON_MIN, HAT_LON_MAX = -75.82, -75.38
HAT_LAT_MIN, HAT_LAT_MAX =  35.15,  35.82

# Domain CRS (from HAT_domains.json)
DOMAIN_CRS  = "EPSG:3725"
TRANSECT_CRS = "EPSG:4326"   # WGS84

# Trajectory classification colour scheme (matches classification script)
CLASS_COLORS = {
    "Persistent Erosion":    "#b2182b",
    "Persistent Accretion":  "#2166ac",
    "Persistently Stable":   "#4dac26",
    "Switching: Acc→Ero":    "#f4a582",
    "Switching: Ero→Acc":    "#92c5de",
    "Transitioning":         "#d9d9d9",
    "Insufficient Data":     "#eeeeee",
}

# Notable domain labels to annotate on map
LABEL_DOMAINS = {
    "Cape Hatteras\n(Buxton)": 7,
    "Avon": 26,
    "Wimble\nShoals": 67,
    "Tri-Village\n(Salvo/Waves/\nRodanthe)": 75,
}

# ============================================================
# HELPERS
# ============================================================

def load_domain_gdf(path):
    """Load domain polygons, reproject to WGS84 for plotting."""
    gdf = gpd.read_file(path)
    gdf = gdf.set_crs(DOMAIN_CRS, allow_override=True)
    gdf = gdf.to_crs("EPSG:4326")
    gdf = gdf.rename(columns={"domain_id": "domain"})
    return gdf


def load_hatteras_transects(path):
    """Load CoastSat transects filtered to Hatteras bounding box."""
    with open(path) as f:
        data = json.load(f)
    feats = []
    for feat in data["features"]:
        coords = feat["geometry"]["coordinates"]
        lon0 = coords[0][0]
        lat0 = coords[0][1]
        if (HAT_LON_MIN < lon0 < HAT_LON_MAX and
                HAT_LAT_MIN < lat0 < HAT_LAT_MAX):
            feats.append(feat)
    gdf = gpd.GeoDataFrame.from_features(feats, crs="EPSG:4326")
    return gdf


def load_island_outline(path):
    """Load Hatteras Island outline shapefile, reproject to WGS84 if needed."""
    gdf = gpd.read_file(path)
    if gdf.crs is None or str(gdf.crs) != "EPSG:4326":
        gdf = gdf.to_crs("EPSG:4326")
    return gdf


def make_lrr_cmap():
    return LinearSegmentedColormap.from_list(
        "lrr_map",
        [(0.00, "#4a0010"),
         (0.15, "#b2182b"),
         (0.35, "#ef8a62"),
         (0.50, "#f7f7f7"),
         (0.65, "#67a9cf"),
         (0.85, "#2166ac"),
         (1.00, "#08306b")]
    )


def add_domain_labels(ax, gdf, label_domains):
    """Annotate selected domains by name."""
    for label, dom_id in label_domains.items():
        row = gdf[gdf["domain"] == dom_id]
        if row.empty:
            continue
        centroid = row.geometry.centroid.iloc[0]
        ax.annotate(
            label,
            xy=(centroid.x, centroid.y),
            xytext=(centroid.x + 0.08, centroid.y + 0.04),
            fontsize=6.5, ha="left", va="center",
            color="#222222",
            arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5),
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#cccccc",
                      alpha=0.85, lw=0.5),
            zorder=6
        )


# ============================================================
# FIGURE 1 — Trajectory class map
# ============================================================

def plot_trajectory_map(domain_gdf, transect_gdf, outline_gdf, metrics, out_path):
    fig, ax = plt.subplots(figsize=(7, 14))
    ax.set_facecolor("#d0e8f5")   # ocean blue
    fig.patch.set_facecolor("white")

    # ── Draw domain polygons coloured by trajectory class ────────────────────
    for _, row in domain_gdf.iterrows():
        dom = row["domain"]
        if dom in metrics.index:
            traj  = metrics.loc[dom, "trajectory"]
            lrr_f = metrics.loc[dom, "dom_lrr_full"]
            color = CLASS_COLORS.get(traj, "#eeeeee")
            # Opacity encodes full-record LRR magnitude (stronger signal = more opaque)
            mag   = min(abs(lrr_f) / 4.0, 1.0) if not np.isnan(lrr_f) else 0.4
            alpha = 0.45 + 0.55 * mag
        else:
            color = CLASS_COLORS["Insufficient Data"]
            alpha = 0.5
            traj  = "Insufficient Data"

        gpd.GeoSeries([row.geometry]).plot(
            ax=ax, color=color, alpha=alpha,
            edgecolor="white", linewidth=0.4, zorder=2
        )

        # High variability hatching overlay
        if dom in metrics.index and metrics.loc[dom, "high_var"]:
            gpd.GeoSeries([row.geometry]).plot(
                ax=ax, color="none",
                edgecolor="#ff7f00", linewidth=1.2,
                linestyle="--", zorder=3
            )

    # ── Island outline ───────────────────────────────────────────────────────
    outline_gdf.plot(ax=ax, facecolor="none", edgecolor="#1a1a1a",
                     linewidth=0.8, zorder=4)

    # ── Draw transect lines as thin grey lines (island structure) ─────────────
    transect_gdf.plot(ax=ax, color="#555555", linewidth=0.15,
                      alpha=0.25, zorder=1)

    # ── Domain number labels for every 5th domain ────────────────────────────
    for _, row in domain_gdf.iterrows():
        dom = row["domain"]
        if dom % 5 == 0:
            c = row.geometry.centroid
            ax.text(c.x, c.y, str(dom),
                    fontsize=5, ha="center", va="center",
                    color="#333333", fontweight="bold", zorder=5)

    add_domain_labels(ax, domain_gdf, LABEL_DOMAINS)

    # ── Formatting ────────────────────────────────────────────────────────────
    ax.set_xlim(HAT_LON_MIN, HAT_LON_MAX)
    ax.set_ylim(HAT_LAT_MIN, HAT_LAT_MAX)
    ax.set_xlabel("Longitude", fontsize=8)
    ax.set_ylabel("Latitude", fontsize=8)
    ax.tick_params(labelsize=7)

    # North arrow
    ax.annotate("N", xy=(0.95, 0.97), xytext=(0.95, 0.94),
                xycoords="axes fraction",
                fontsize=10, ha="center", fontweight="bold",
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.5))

    # Scale bar (rough — 0.1° lon ≈ ~8.5 km at 35°N)
    sb_x0, sb_y = HAT_LON_MIN + 0.05, HAT_LAT_MIN + 0.03
    ax.plot([sb_x0, sb_x0 + 0.1], [sb_y, sb_y],
            color="black", lw=2, zorder=6)
    ax.text(sb_x0 + 0.05, sb_y + 0.01, "~8.5 km",
            ha="center", fontsize=6, zorder=6)

    # Legend
    patches = []
    for cls, col in CLASS_COLORS.items():
        patches.append(mpatches.Patch(color=col, label=cls, alpha=0.85))
    patches.append(mpatches.Patch(facecolor="none", edgecolor="#ff7f00",
                                   linestyle="--", linewidth=1.2,
                                   label="High variability"))
    ax.legend(handles=patches, fontsize=6.5, loc="lower left",
              framealpha=0.92, title="Trajectory class",
              title_fontsize=7, edgecolor="#cccccc")

    # Opacity note
    ax.text(0.02, 0.01,
            "Opacity encodes |LRR| magnitude\n(darker = faster rate of change)",
            transform=ax.transAxes, fontsize=5.5, color="#555555",
            va="bottom")

    ax.set_title(
        "Hatteras Island\nShoreline trajectory classification\n"
        "CoastSat 1984–2024  ·  threshold ±1 m/yr",
        fontsize=10, fontweight="bold", pad=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"  Trajectory map saved → {out_path}")


# ============================================================
# FIGURE 2 — LRR magnitude map (continuous colour scale)
# ============================================================

def plot_lrr_map(domain_gdf, transect_gdf, outline_gdf, metrics, out_path):
    lrr_full = metrics["dom_lrr_full"].values
    valid     = lrr_full[~np.isnan(lrr_full)]
    vabs      = np.ceil(max(abs(valid.min()), abs(valid.max())) / 0.5) * 0.5

    lrr_cmap = make_lrr_cmap()
    norm      = Normalize(vmin=-vabs, vmax=vabs)

    fig, axes = plt.subplots(1, 2, figsize=(14, 14),
                             gridspec_kw={"wspace": 0.05})

    for ax, period, lrr_col, title_sfx in [
        (axes[0], "Period 1\n(1984–2004)", "dom_lrr_p1", "Period 1  (1984–2004)"),
        (axes[1], "Period 2\n(2004–2024)", "dom_lrr_p2", "Period 2  (2004–2024)"),
    ]:
        ax.set_facecolor("#d0e8f5")

        for _, row in domain_gdf.iterrows():
            dom = row["domain"]
            lrr = metrics.loc[dom, lrr_col] if dom in metrics.index else np.nan
            if np.isnan(lrr):
                color, alpha = "#cccccc", 0.5
            else:
                color = lrr_cmap(norm(lrr))
                alpha = 0.90

            gpd.GeoSeries([row.geometry]).plot(
                ax=ax, color=color, alpha=alpha,
                edgecolor="white", linewidth=0.3, zorder=2
            )

        transect_gdf.plot(ax=ax, color="#333333", linewidth=0.12,
                          alpha=0.20, zorder=1)
        outline_gdf.plot(ax=ax, facecolor="none", edgecolor="#1a1a1a",
                         linewidth=0.8, zorder=4)

        # Domain number labels every 5
        for _, row in domain_gdf.iterrows():
            if row["domain"] % 5 == 0:
                c = row.geometry.centroid
                ax.text(c.x, c.y, str(row["domain"]),
                        fontsize=4.5, ha="center", va="center",
                        color="#222222", fontweight="bold", zorder=5)

        ax.set_xlim(HAT_LON_MIN, HAT_LON_MAX)
        ax.set_ylim(HAT_LAT_MIN, HAT_LAT_MAX)
        ax.set_xlabel("Longitude", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_title(title_sfx, fontsize=10, fontweight="bold", pad=6)

    axes[0].set_ylabel("Latitude", fontsize=8)
    axes[1].set_yticklabels([])

    # Shared colorbar
    sm = cm.ScalarMappable(cmap=lrr_cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation="horizontal",
                        pad=0.04, shrink=0.6, aspect=40)
    cbar.set_label("Linear rate of change  (m/yr)\nblue = accretion  ·  red = erosion",
                   fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    fig.suptitle(
        "Hatteras Island  |  Shoreline linear rate of change by period\n"
        "CoastSat domain-averaged LRR  ·  1984–2024",
        fontsize=11, fontweight="bold", y=0.98)

    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"  LRR map saved → {out_path}")


# ============================================================
# MAIN
# ============================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Loading island outline …")
    outline_gdf = load_island_outline(OUTLINE_SHP)
    print(f"  {len(outline_gdf)} outline features loaded.\n")

    print("Loading domain polygons …")
    domain_gdf = load_domain_gdf(DOMAINS_GEOJSON)
    print(f"  {len(domain_gdf)} domain polygons loaded, reprojected to WGS84.")

    print("Loading Hatteras transects …")
    transect_gdf = load_hatteras_transects(TRANSECTS_GEOJSON)
    print(f"  {len(transect_gdf)} Hatteras transects loaded.")

    print("Loading trajectory metrics …")
    metrics = pd.read_csv(METRICS_CSV, index_col="domain")
    print(f"  {len(metrics)} domains in metrics table.")

    print("\nGenerating figures …")

    plot_trajectory_map(
        domain_gdf, transect_gdf, outline_gdf, metrics,
        os.path.join(OUTPUT_DIR, "fig_map_trajectory.png")
    )

    plot_lrr_map(
        domain_gdf, transect_gdf, outline_gdf, metrics,
        os.path.join(OUTPUT_DIR, "fig_map_lrr.png")
    )

    print("\nDone.")


if __name__ == "__main__":
    main()
