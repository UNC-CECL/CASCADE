"""
HAT_trajectory_map.py
======================
Geographic maps of Hatteras Island with CASCADE domains coloured by
shoreline trajectory classification or LRR magnitude.

CONFIG options:
  USE_SATELLITE = True   → Esri WorldImagery satellite basemap (requires internet)
  USE_SATELLITE = False  → plain ocean-blue background (no internet needed)

Dependencies
------------
  pip install geopandas contextily shapely matplotlib pyproj
  pip install matplotlib-scalebar   (optional — for accurate scale bar)
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap, Normalize
import geopandas as gpd
import matplotlib.patheffects as pe

# ============================================================
# CONFIG  — edit here
# ============================================================

DOMAINS_GEOJSON   = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\HAT_domains.json"
TRANSECTS_GEOJSON = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\CoastSat_transect_layer.geojson"
METRICS_CSV       = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\shoreline_change_patterns\output\domain_trajectory_metrics.csv"
OUTLINE_SHP       = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\figure_making\map_elements\hatteras_outline\HAT_island_outline.shp"
OUTPUT_DIR        = r"/scripts/input_prep/shoreline_change_patterns/old_drafts/trajectory_map_v2"

# ── Basemap choice ────────────────────────────────────────────────────────────
USE_SATELLITE = True   # True = Esri satellite tiles; False = plain ocean blue

OCEAN_COLOR   = "#d0e8f5"   # used when USE_SATELLITE = False
TILE_ZOOM     = 13           # satellite tile zoom level (12=regional, 14=detailed)

# ── Bounding box (WGS84) ──────────────────────────────────────────────────────
HAT_LON_MIN, HAT_LON_MAX = -75.82, -75.38
HAT_LAT_MIN, HAT_LAT_MAX =  35.15,  35.82

# Web Mercator (EPSG:3857) bounds — used to set axis limits for satellite tiles
# (computed from WGS84 bounds above via pyproj)
HAT_X_MIN, HAT_X_MAX = -8_440_244, -8_391_263
HAT_Y_MIN, HAT_Y_MAX =  4_184_284,  4_275_882

# ── CRS ───────────────────────────────────────────────────────────────────────
DOMAIN_CRS = "EPSG:3725"   # UTM — HAT_domains.json native
# Satellite tiles require Web Mercator; plain mode uses WGS84
PLOT_CRS = "EPSG:3857" if USE_SATELLITE else "EPSG:4326"

# ── Trajectory colours (matches classification script) ────────────────────────
CLASS_COLORS = {
    "Persistent Erosion":   "#b2182b",
    "Persistent Accretion": "#2166ac",
    "Persistently Stable":  "#4dac26",
    "Switching: Acc→Ero":   "#f4a582",
    "Switching: Ero→Acc":   "#92c5de",
    "Transitioning":        "#d9d9d9",
    "Insufficient Data":    "#eeeeee",
}

# ── Place-name annotations ────────────────────────────────────────────────────
LABEL_DOMAINS = {
    "Buxton":                               7,
    "Avon":                                26,
    "Wimble Shoals":                       67,
    "Tri-Village\n(Salvo/Waves/Rodanthe)": 75,
}

# ============================================================
# BASEMAP HELPERS
# ============================================================

def apply_basemap(ax):
    """Add satellite tiles or plain background depending on USE_SATELLITE."""
    if USE_SATELLITE:
        # Must set axis limits in Web Mercator BEFORE adding tiles
        ax.set_xlim(HAT_X_MIN, HAT_X_MAX)
        ax.set_ylim(HAT_Y_MIN, HAT_Y_MAX)
        import contextily as ctx
        ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery,
                        crs=PLOT_CRS, zoom=TILE_ZOOM, attribution_size=5)
    else:
        ax.set_facecolor(OCEAN_COLOR)


def label_color():
    """White labels on satellite, dark labels on plain background."""
    return "white" if USE_SATELLITE else "#1a1a1a"


def annotation_box_style():
    return dict(boxstyle="round,pad=0.25",
                fc="#222222" if USE_SATELLITE else "white",
                ec=label_color(), alpha=0.80, lw=0.5)


# ============================================================
# DATA LOADERS
# ============================================================

def load_island_outline(path):
    gdf = gpd.read_file(path)
    if gdf.crs is None or str(gdf.crs) != "EPSG:4326":
        gdf = gdf.to_crs("EPSG:4326")
    return gdf.to_crs(PLOT_CRS)


def load_domain_gdf(path, outline_gdf_plot_crs):
    """Load domain polygons, clip to island outline, reproject to PLOT_CRS."""
    gdf = gpd.read_file(path)
    gdf = gdf.set_crs(DOMAIN_CRS, allow_override=True).to_crs("EPSG:4326")
    gdf = gdf.rename(columns={"domain_id": "domain"})
    # Clip rectangular domains to actual island shape
    island_wgs = outline_gdf_plot_crs.to_crs("EPSG:4326").union_all()
    gdf["geometry"] = gdf.geometry.intersection(island_wgs)
    gdf = gdf[~gdf.geometry.is_empty]
    return gdf.to_crs(PLOT_CRS)


def load_hatteras_transects(path):
    with open(path) as f:
        data = json.load(f)
    feats = [f for f in data["features"]
             if (HAT_LON_MIN < f["geometry"]["coordinates"][0][0] < HAT_LON_MAX
                 and HAT_LAT_MIN < f["geometry"]["coordinates"][0][1] < HAT_LAT_MAX)]
    return gpd.GeoDataFrame.from_features(feats, crs="EPSG:4326").to_crs(PLOT_CRS)


def make_lrr_cmap():
    return LinearSegmentedColormap.from_list(
        "lrr_map",
        [(0.00, "#4a0010"), (0.15, "#b2182b"), (0.35, "#ef8a62"),
         (0.50, "#f7f7f7"), (0.65, "#67a9cf"), (0.85, "#2166ac"),
         (1.00, "#08306b")]
    )


# ============================================================
# ANNOTATION HELPERS
# ============================================================

def add_north_arrow(ax):
    lc = label_color()
    ax.annotate("", xy=(0.95, 0.97), xytext=(0.95, 0.91),
                xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=lc,
                                lw=2.0, mutation_scale=14))
    ax.text(0.95, 0.98, "N", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=9,
            color=lc, fontweight="bold")


def add_scalebar(ax):
    try:
        from matplotlib_scalebar.scalebar import ScaleBar
        sb = ScaleBar(1, units="m", length_fraction=0.15,
                      location="lower right", color=label_color(),
                      box_alpha=0.0, font_properties={"size": 7})
        ax.add_artist(sb)
    except ImportError:
        # Manual fallback: draw a ~10 km bar in map units
        xlim = ax.get_xlim(); ylim = ax.get_ylim()
        bar_len = (xlim[1] - xlim[0]) * 0.12
        x0 = xlim[0] + 0.74 * (xlim[1] - xlim[0])
        y0 = ylim[0] + 0.03 * (ylim[1] - ylim[0])
        ax.plot([x0, x0 + bar_len], [y0, y0],
                color=label_color(), lw=2.5, zorder=10)
        ax.text(x0 + bar_len / 2, y0 + (ylim[1]-ylim[0])*0.01,
                "~10 km", ha="center", color=label_color(),
                fontsize=6, zorder=10)


def draw_outline(ax, outline_gdf):
    outline_gdf.plot(ax=ax, facecolor="none",
                     edgecolor=label_color(), linewidth=0.9, zorder=5)


def draw_domain_numbers(ax, domain_gdf, every=5, fontsize=6):
    for _, row in domain_gdf.iterrows():
        if row["domain"] % every == 0:
            c = row.geometry.centroid
            ax.text(c.x, c.y, str(row["domain"]),
                    fontsize=fontsize, ha="center", va="center",
                    color="white", fontweight="bold", zorder=6,
                    path_effects=[
                        pe.withStroke(linewidth=2, foreground="black")])


def draw_place_labels(ax, domain_gdf):
    lc = label_color()
    for label, dom_id in LABEL_DOMAINS.items():
        row = domain_gdf[domain_gdf["domain"] == dom_id]
        if row.empty:
            continue
        c = row.geometry.centroid.iloc[0]
        ax.annotate(
            label,
            xy=(c.x, c.y),
            xytext=(c.x + (ax.get_xlim()[1]-ax.get_xlim()[0])*0.06, c.y),
            fontsize=6.5, ha="left", va="center", color=lc,
            arrowprops=dict(arrowstyle="-", color=lc, lw=0.6),
            bbox=annotation_box_style(),
            zorder=8
        )


# ============================================================
# FIGURE 1 — Trajectory class map
# ============================================================

def plot_trajectory_map(domain_gdf, outline_gdf, metrics, out_path):
    fig = plt.figure(figsize=(8, 17))
    # Map panel + thin legend panel below
    gs  = plt.GridSpec(2, 1, height_ratios=[15, 1], hspace=0.03)
    ax     = fig.add_subplot(gs[0])
    ax_leg = fig.add_subplot(gs[1])
    ax_leg.axis("off")

    # Basemap
    apply_basemap(ax)

    # Domain polygons
    for _, row in domain_gdf.iterrows():
        dom = row["domain"]
        if dom in metrics.index:
            traj  = metrics.loc[dom, "trajectory"]
            lrr_f = metrics.loc[dom, "dom_lrr_full"]
            color = CLASS_COLORS.get(traj, "#eeeeee")
            if USE_SATELLITE:
                # Satellite: high base opacity, slight magnitude modulation
                mag   = min(abs(lrr_f) / 4.0, 1.0) if not np.isnan(lrr_f) else 0.5
                alpha = 0.72 + 0.25 * mag   # range 0.72–0.97
            else:
                mag   = min(abs(lrr_f) / 4.0, 1.0) if not np.isnan(lrr_f) else 0.3
                alpha = 0.45 + 0.50 * mag
        else:
            color, alpha = CLASS_COLORS["Insufficient Data"], 0.4

        gpd.GeoSeries([row.geometry], crs=PLOT_CRS).plot(
            ax=ax, color=color, alpha=alpha,
            edgecolor="white", linewidth=0.5, zorder=3)

        # High-variability: orange dashed border
        if dom in metrics.index and metrics.loc[dom, "high_var"]:
            gpd.GeoSeries([row.geometry], crs=PLOT_CRS).plot(
                ax=ax, color="none", edgecolor="#ff7f00",
                linewidth=1.4, linestyle="--", zorder=4)

    draw_outline(ax, outline_gdf)
    draw_domain_numbers(ax, domain_gdf)
    draw_place_labels(ax, domain_gdf)
    add_north_arrow(ax)
    add_scalebar(ax)

    ax.set_axis_off()
    ax.set_title(
        "Hatteras Island  |  Shoreline trajectory classification\n"
        "CoastSat 1984–2024  ·  threshold ±1 m/yr",
        fontsize=10, fontweight="bold", pad=6)

    # Legend — outside map in its own panel
    patches = [mpatches.Patch(color=v, label=k, alpha=0.85)
               for k, v in CLASS_COLORS.items()]
    patches.append(mpatches.Patch(facecolor="none", edgecolor="#ff7f00",
                                  linestyle="--", linewidth=1.5,
                                  label="High variability"))
    ax_leg.legend(handles=patches, fontsize=7, loc="center",
                  ncol=4, framealpha=0.0,
                  title="Trajectory class  (opacity encodes |LRR| magnitude — darker = faster change)",
                  title_fontsize=7, edgecolor="none", columnspacing=1.0)

    fig.savefig(out_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Trajectory map saved → {out_path}")


# ============================================================
# FIGURE 2 — LRR magnitude map, Period 1 vs Period 2
# ============================================================

def plot_lrr_map(domain_gdf, outline_gdf, metrics, out_path):
    lrr_vals = np.concatenate([
        metrics["dom_lrr_p1"].dropna().values,
        metrics["dom_lrr_p2"].dropna().values
    ])
    vabs     = np.ceil(max(abs(lrr_vals.min()), abs(lrr_vals.max())) / 0.5) * 0.5
    lrr_cmap = make_lrr_cmap()
    norm      = Normalize(vmin=-vabs, vmax=vabs)

    # Two map panels side by side + shared colorbar below
    fig = plt.figure(figsize=(16, 17))
    gs  = plt.GridSpec(2, 2, height_ratios=[15, 0.7],
                       hspace=0.03, wspace=0.04,
                       top=0.97, bottom=0.04, left=0.02, right=0.98)
    axes  = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
    ax_cb = fig.add_subplot(gs[1, :])
    ax_cb.axis("off")

    for ax, lrr_col, title_sfx in [
        (axes[0], "dom_lrr_p1", "Period 1  (1984–2004)"),
        (axes[1], "dom_lrr_p2", "Period 2  (2004–2024)"),
    ]:
        apply_basemap(ax)

        for _, row in domain_gdf.iterrows():
            dom = row["domain"]
            lrr = metrics.loc[dom, lrr_col] if dom in metrics.index else np.nan
            color = lrr_cmap(norm(lrr)) if not np.isnan(lrr) else "#cccccc"
            alpha = 0.85 if not np.isnan(lrr) else 0.4
            gpd.GeoSeries([row.geometry], crs=PLOT_CRS).plot(
                ax=ax, color=color, alpha=alpha,
                edgecolor="white", linewidth=0.4, zorder=3)

        draw_outline(ax, outline_gdf)
        draw_domain_numbers(ax, domain_gdf, fontsize=4)
        add_north_arrow(ax)
        add_scalebar(ax)
        ax.set_axis_off()
        prefix = "Hatteras Island  |  CoastSat LRR  ·  "
        ax.set_title(prefix + title_sfx, fontsize=9, fontweight="bold", pad=4)

    # Colorbar in its own panel below both maps
    sm = cm.ScalarMappable(cmap=lrr_cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax_cb, orientation="horizontal",
                        fraction=0.6, shrink=0.50, aspect=40, pad=0.0)
    cbar.set_label(
        "Linear rate of change  (m/yr)  —  "
        "red = erosion  ·  grey = stable  ·  blue = accretion",
        fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    fig.savefig(out_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  LRR map saved → {out_path}")


# ============================================================
# MAIN
# ============================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print(f"Basemap mode: {'satellite' if USE_SATELLITE else 'plain blue'}\n")

    print("Loading island outline …")
    outline_gdf = load_island_outline(OUTLINE_SHP)
    print(f"  {len(outline_gdf)} outline features.\n")

    print("Loading domain polygons (clipping to island outline) …")
    domain_gdf = load_domain_gdf(DOMAINS_GEOJSON, outline_gdf)
    print(f"  {len(domain_gdf)} domain polygons loaded and clipped.\n")

    print("Loading trajectory metrics …")
    metrics = pd.read_csv(METRICS_CSV, index_col="domain")
    print(f"  {len(metrics)} domains.\n")

    print("Generating figures …")

    plot_trajectory_map(
        domain_gdf, outline_gdf, metrics,
        os.path.join(OUTPUT_DIR, "fig_map_trajectory.png"))

    plot_lrr_map(
        domain_gdf, outline_gdf, metrics,
        os.path.join(OUTPUT_DIR, "fig_map_lrr.png"))

    print("\nDone.")


if __name__ == "__main__":
    main()
