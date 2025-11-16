#!/usr/bin/env python3
"""
CASCADE: Domain Heatmaps for Road Elevation Validation — Publication-Ready Split Panels
-------------------------------------------------------------------------------------
Enhanced version with improved layout, styling, and readability for publication.
"""

import os
import argparse
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

# Optional scalebar
try:
    from matplotlib_scalebar.scalebar import ScaleBar

    HAVE_SCALEBAR = True
except Exception:
    HAVE_SCALEBAR = False

# ------------------------
# ENHANCED CONFIG
# ------------------------
CONFIG = {
    # Core data
    "domains_path": r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\domains\HAT_domains.shp",
    "domains_layer": None,
    "domains_id_field": "ID",
    "merged_csv": r"outputs\road_validation\HAT_1978_table_merged.csv",

    # Context layers
    "island_outline": r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\maps\HAT_island_outline.shp",
    "road_centerline": r"",

    # Split settings
    "split_enabled": True,
    "split_threshold": 46,

    # Output and styling
    "outdir": r"outputs",
    "tag": "HAT_1978",
    "dpi": 600,
    "figsize": (10, 8),  # Wider to prevent overlap
    "color_cmap": "RdBu_r",
    "use_symmetric_limits": True,
    "fixed_lim_m": 0.40,

    # Enhanced styling
    "linewidth_domains": 0.3,
    "linewidth_context": 0.8,  # Thicker for visibility
    "context_gray": 0.15,  # Darker for better contrast
    "ocean_color": "#e6f2ff",  # Light blue background
    "land_color": "#f5f5f0",  # Light tan background

    # Labels
    "domain_label": True,
    "domain_label_field": "ID",
    "label_every_n": 1,  # Label ALL domains
    "label_fontsize": 5.5,  # Adjusted for all labels
    "label_halo_width": 1.8,  # Slightly reduced halo

    # Publication style
    "use_panel_labels": True,  # Add (a), (b) labels
    "panel_label_size": 12,
    "title_size": 14,
    "subtitle_size": 11,
    "cbar_label_size": 11,
    "cbar_tick_size": 9,
}

METRICS = [
    ("Diff_Buffer_vs_Transect", "Elevation difference (Buffer − Transect) [m]"),
    ("Diff_Road10m_vs_Domain10m", "Elevation difference (Road 10 m − Domain 10 m) [m]"),
]


# ------------------------
# Helpers
# ------------------------
def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def read_gdf(path, layer=None):
    if not path:
        return None
    if not os.path.exists(path):
        print(f"[warn] Context layer not found: {path}")
        return None
    if path.lower().endswith(".gpkg"):
        return gpd.read_file(path, layer=layer)
    return gpd.read_file(path)


def to_match_crs(gdf_reference, gdf_other):
    if gdf_other is None:
        return None
    try:
        if gdf_other.crs != gdf_reference.crs:
            return gdf_other.to_crs(gdf_reference.crs)
        return gdf_other
    except Exception:
        return gdf_other


def add_north_arrow(ax, xy=(0.5, 0.02), size=0.025, cfg=None):
    """Cleaner north arrow centered at bottom"""
    from matplotlib.patches import FancyArrowPatch, Polygon
    # Create a filled triangle arrow
    arrow_width = 0.012
    arrow_head = Polygon([
        (xy[0], xy[1] + size),  # tip
        (xy[0] - arrow_width, xy[1]),  # left base
        (xy[0] + arrow_width, xy[1])  # right base
    ], transform=ax.transAxes, facecolor='black', edgecolor='black',
        linewidth=0.5, zorder=10)
    ax.add_patch(arrow_head)

    # N label
    ax.text(xy[0], xy[1] + size + 0.012, "N", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=9, fontweight="bold",
            color="black", zorder=10)


def add_scale_bar(ax, cfg):
    """Enhanced scale bar positioning"""
    if not HAVE_SCALEBAR:
        return
    try:
        scalebar = ScaleBar(dx=1, units="m", location="lower left",
                            box_alpha=0.7, box_color="white",
                            color="black", length_fraction=0.25,
                            font_properties={'size': 8})
        ax.add_artist(scalebar)
    except Exception as e:
        print(f"[warn] Could not add scale bar: {e}")


def compute_norm(series, fixed_lim, symmetric=True):
    s = series.dropna()
    if s.empty:
        return None, None
    if fixed_lim is not None:
        lim = float(fixed_lim)
    else:
        vmin, vmax = s.min(), s.max()
        lim = max(abs(vmin), abs(vmax)) if symmetric else None
    if symmetric:
        return TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim), (-lim, lim)
    else:
        vmin, vmax = s.min(), s.max()
        return TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax), (vmin, vmax)


def panel_extent(gdf_aoi, pad_frac=0.08):
    """Increased padding to prevent panel label overlap"""
    xmin, ymin, xmax, ymax = gdf_aoi.total_bounds
    dx = (xmax - xmin) * pad_frac
    dy = (ymax - ymin) * pad_frac
    return (xmin - dx, xmax + dx, ymin - dy, ymax + dy)


def add_background(ax, gdf_aoi, cfg, pad_frac=0.08):
    """Add ocean background for context"""
    xmin, ymin, xmax, ymax = gdf_aoi.total_bounds
    dx = (xmax - xmin) * pad_frac
    dy = (ymax - ymin) * pad_frac
    extent = (xmin - dx, xmax + dx, ymin - dy, ymax + dy)
    rect = Rectangle((extent[0], extent[2]),
                     extent[1] - extent[0], extent[3] - extent[2],
                     facecolor=cfg["ocean_color"], edgecolor='none',
                     zorder=0, transform=ax.transData)
    ax.add_patch(rect)


def plot_panel(ax, gdf_subset, ctx_island, ctx_road, field, norm, cfg,
               label_offset=0, panel_letter=None):
    """Enhanced panel plotting with better styling"""
    ax.set_axis_off()

    # Set extent first
    extent = panel_extent(gdf_subset if not gdf_subset.empty else gdf_subset)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])

    # Background
    add_background(ax, gdf_subset, cfg)

    # Island outline (filled for context)
    if ctx_island is not None and not ctx_island.empty:
        ctx_island.plot(ax=ax, facecolor=cfg["land_color"],
                        edgecolor=str(cfg["context_gray"]),
                        linewidth=cfg["linewidth_context"], alpha=1.0, zorder=1)

    # Domains (main data)
    gdf_subset.plot(column=field, ax=ax, cmap=cfg["color_cmap"], norm=norm,
                    linewidth=cfg["linewidth_domains"], edgecolor="k",
                    alpha=0.9, zorder=2)

    # Road centerline
    if ctx_road is not None and not ctx_road.empty:
        ctx_road.plot(ax=ax, color="black", linewidth=cfg["linewidth_context"] * 0.8,
                      linestyle="-", alpha=0.7, zorder=3)

    # Domain labels with enhanced halo
    if cfg["domain_label"] and cfg["domain_label_field"] in gdf_subset.columns:
        reps = gdf_subset.representative_point()
        for i, ((x, y), label) in enumerate(
                zip(reps.geometry.apply(lambda p: (p.x, p.y)),
                    gdf_subset[cfg["domain_label_field"]])
        ):
            if (i + label_offset) % cfg["label_every_n"] == 0:
                txt = ax.text(x, y, str(label), ha="center", va="center",
                              fontsize=cfg["label_fontsize"],
                              color="#333333", alpha=1.0, fontweight="normal",
                              zorder=4)
                txt.set_path_effects([
                    path_effects.withStroke(linewidth=cfg["label_halo_width"],
                                            foreground="white", alpha=0.9)
                ])

    # Panel label (a, b, etc.)
    if panel_letter and cfg["use_panel_labels"]:
        ax.text(0.02, 0.98, f"({panel_letter})", transform=ax.transAxes,
                fontsize=cfg["panel_label_size"], fontweight="bold",
                va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="none", alpha=0.8))

    # North arrow at bottom right
    add_north_arrow(ax, xy=(0.92, 0.05), size=0.035, cfg=cfg)
    if panel_letter == "a":  # Only on first panel
        add_scale_bar(ax, cfg)


# ------------------------
# Main
# ------------------------
def main():
    ap = argparse.ArgumentParser(description="Publication-ready split-panel domain heatmaps.")
    ap.add_argument("--domains_path")
    ap.add_argument("--domains_layer")
    ap.add_argument("--domains_id")
    ap.add_argument("--merged_csv")
    ap.add_argument("--outdir")
    ap.add_argument("--tag")
    ap.add_argument("--island_outline")
    ap.add_argument("--road_centerline")
    args = ap.parse_args()

    cfg = CONFIG.copy()
    if args.domains_path:      cfg["domains_path"] = args.domains_path
    if args.domains_layer is not None: cfg["domains_layer"] = args.domains_layer
    if args.domains_id:        cfg["domains_id_field"] = args.domains_id
    if args.merged_csv:        cfg["merged_csv"] = args.merged_csv
    if args.outdir:            cfg["outdir"] = args.outdir
    if args.tag:               cfg["tag"] = args.tag
    if args.island_outline:    cfg["island_outline"] = args.island_outline
    if args.road_centerline:   cfg["road_centerline"] = args.road_centerline

    # Load data
    gdf = read_gdf(cfg["domains_path"], cfg["domains_layer"])
    if gdf is None:
        raise SystemExit(f"Domains not found: {cfg['domains_path']}")
    df = pd.read_csv(cfg["merged_csv"])

    # Find ID columns
    csv_id = next((c for c in ["Domain_ID", "domain_id", "ID", "id"] if c in df.columns), None)
    if csv_id is None:
        raise SystemExit("CSV needs a domain id column ('Domain_ID' or 'ID').")

    dom_id = cfg["domains_id_field"]
    if dom_id not in gdf.columns:
        match = next((c for c in gdf.columns if c.lower() == dom_id.lower()), None)
        if not match:
            raise SystemExit(f"'{dom_id}' not in polygons. Fields: {list(gdf.columns)}")
        dom_id = match

    gdf_join = gdf.merge(df, left_on=dom_id, right_on=csv_id, how="left", validate="1:1")

    # Context layers
    island = to_match_crs(gdf_join, read_gdf(cfg["island_outline"]))
    road = to_match_crs(gdf_join, read_gdf(cfg["road_centerline"]))

    # Split subsets
    if cfg["split_enabled"]:
        south = gdf_join[gdf_join[dom_id] < cfg["split_threshold"]].copy()
        north = gdf_join[gdf_join[dom_id] >= cfg["split_threshold"]].copy()
    else:
        south = gdf_join.copy()
        north = None

    # Process each metric
    for field, cbar_label in METRICS:
        if field not in gdf_join.columns:
            print(f"[info] Skipping '{field}' (not in merged CSV).")
            continue

        # Compute shared normalization
        norm, (vmin, vmax) = compute_norm(gdf_join[field], cfg["fixed_lim_m"],
                                          cfg["use_symmetric_limits"])
        if norm is None:
            print(f"[info] Skipping '{field}' (no finite values).")
            continue

        # Create figure with gridspec for precise control
        if cfg["split_enabled"] and north is not None and not north.empty and not south.empty:
            fig = plt.figure(figsize=cfg["figsize"])
            gs = gridspec.GridSpec(1, 3, figure=fig, width_ratios=[1, 1, 0.04],
                                   wspace=0.01, hspace=0.05,
                                   left=0.02, right=0.93, top=0.93, bottom=0.08)

            axL = fig.add_subplot(gs[0, 0])
            axR = fig.add_subplot(gs[0, 1])
            cax = fig.add_subplot(gs[0, 2])

            # Shorter, more concise title
            year = cfg['tag'].split('_')[-1]
            fig.suptitle(f"Elevation Difference: Buffer vs. Transect ({year})",
                         fontsize=cfg["title_size"], fontweight="bold", y=0.97)

            # Plot panels
            plot_panel(axL, south, island, road, field, norm, cfg,
                       label_offset=0, panel_letter="a")
            plot_panel(axR, north, island, road, field, norm, cfg,
                       label_offset=1, panel_letter="b")

            # Subtitle for each panel
            axL.text(0.5, -0.02, "Southern Section", transform=axL.transAxes,
                     ha="center", va="top", fontsize=cfg["subtitle_size"],
                     fontweight="semibold")
            axR.text(0.5, -0.02, "Northern Section", transform=axR.transAxes,
                     ha="center", va="top", fontsize=cfg["subtitle_size"],
                     fontweight="semibold")

            # Colorbar
            sm = plt.cm.ScalarMappable(cmap=cfg["color_cmap"], norm=norm)
            sm._A = []
            cbar = plt.colorbar(sm, cax=cax)
            cbar.set_label(cbar_label, fontsize=cfg["cbar_label_size"],
                           labelpad=12)
            cbar.ax.tick_params(labelsize=cfg["cbar_tick_size"])

            # Emphasize zero line
            cbar.ax.axhline(0, color='black', linewidth=1.2, alpha=0.7)

            # Better tick formatting
            cbar.locator = mticker.MaxNLocator(nbins=9, symmetric=True)
            cbar.update_ticks()

        else:
            # Single panel fallback
            fig, ax = plt.subplots(figsize=cfg["figsize"])
            year = cfg['tag'].split('_')[-1]
            ax.set_title(f"Road Elevation Difference: Buffer versus Transect Method ({year})",
                         fontsize=cfg["title_size"], fontweight="bold", pad=12)
            plot_panel(ax, gdf_join, island, road, field, norm, cfg)

            sm = plt.cm.ScalarMappable(cmap=cfg["color_cmap"], norm=norm)
            sm._A = []
            cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
            cbar.set_label(cbar_label, fontsize=cfg["cbar_label_size"])
            cbar.ax.tick_params(labelsize=cfg["cbar_tick_size"])
            cbar.ax.axhline(0, color='black', linewidth=1.2, alpha=0.7)
            cbar.locator = mticker.MaxNLocator(nbins=9, symmetric=True)
            cbar.update_ticks()

        # Save outputs
        out_dir = ensure_dir(os.path.join(cfg["outdir"], "road_validation", "hatteras_outline"))
        base = os.path.join(out_dir, f"{cfg['tag']}_{field}_pubready")

        plt.savefig(base + ".png", dpi=cfg["dpi"], bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        plt.savefig(base + ".pdf", bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        plt.close(fig)
        print(f"[ok] Wrote publication-ready figure: {base}.png")

    # Export GeoPackage
    out_gpkg = os.path.join(cfg["outdir"], "road_validation", "hatteras_outline",
                            f"{cfg['tag']}_domain_heatmap.gpkg")
    try:
        gdf_join.to_file(out_gpkg, layer="domains_joined", driver="GPKG")
        print(f"[ok] Wrote joined GeoPackage: {out_gpkg}")
    except Exception as e:
        print(f"[warn] Could not write GeoPackage ({e}).")


if __name__ == "__main__":
    main()