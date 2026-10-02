"""
Observed CoastSat shoreline position through time either side of the Buxton groin field.

    python scripts/figure_making/shoreline/buxton_groin_shoreline_position.py

Inputs: the CoastSat transect series and transect layer under 5-scr, the domain boxes,
the groin position from hatteras_site_config and the Buxton fill from nourishment_projects.csv.
Basemap tiles need the network the first time (cached after). Needs geopandas, contextily, rasterio.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-02
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "5-scr" / "1-observations" / "mean_shoreline"))
import scr_paths  # noqa: E402,F401
from coastsat_lrr import compute_lrr, filter_dates, load_timeseries  # noqa: E402
from coastsat_mean_shoreline import timeseries_file, transect_geometry  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    INK, INK_MUTED, _north_arrow, _scalebar, apply_style, caption,
    figsize, figure_dir, open_frame, save, spines_for_image, support_dir, title)
from site_layer.hat_observed_rates import DOMAIN_BOXES, transect_lookup  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
CRS = "EPSG:26918"
GIS_RANGE = (3, 8)                 # domains either side of the groin field (GIS 5.5)
YEARS = (1984, 2025)               # whole record bar the January 2026 passes
BASE_YEARS = (1984, 1986)          # the reference position: median of these years' passes
END_YEARS = (2023, 2025)           # the latest position drawn on the map
MIN_OBS_PER_YEAR = 3               # fewer passes in a year -> that year is left blank
NOURISHMENT = REPO / "data" / "hatteras_init" / "4-mgmt-forcing" / "nourishment" / "nourishment_projects.csv"
FILL_NAME = "Buxton shore protection"
TILE_ZOOM = 16
TILE_CACHE = Path(tempfile.gettempdir()) / "hat_tile_cache"
MAP_PAD_X = (350.0, 850.0)         # metres landward / seaward of the shoreline on the map
C_SOUTH, C_NORTH = "#4d9221", "#e69f00"   # green south, amber north: clear of the red-blue scales
CMAP = "RdBu"                      # red landward, blue seaward
C_BASE_LINE, C_END_LINE = "0.35", "0.05"   # neutral: red and blue mean landward/seaward here
LS_BASE = (0, (3, 1.5))
LS_GROIN = (0, (5, 2, 1, 2))
FILL_STYLE = dict(color=INK, lw=1.2, ls=":")
CLIM = 80.0                        # colour limit, m
RATE_CLIM = 2.5                    # rate colour limit, m/yr; every transect falls inside it
OUT = figure_dir("observations", "shoreline")
STEM = "buxton_groin_shoreline_position"
# -----------------------------------------------------------------------------


# Northing of the groin field: GIS position 5.5 between the domain-box centroids
def groin_northing() -> tuple[float, float]:
    import geopandas as gpd
    boxes = gpd.read_file(DOMAIN_BOXES).to_crs(CRS)
    cy = boxes.geometry.centroid.y.to_numpy()
    pos = HATTERAS_ANNOTATIONS.groins["Buxton Groin"]
    gis = np.arange(1, len(boxes) + 1)          # file order is south -> north, GIS 1-90
    return float(np.interp(pos, gis, cy)), pos


# Transects of the chosen domains with origin and seaward unit vector
def transects() -> pd.DataFrame:
    lk = pd.read_csv(transect_lookup())
    lk = lk[lk["domain_number"].between(*GIS_RANGE)].copy()
    geom = transect_geometry(set(lk["transect_id"]))
    lk = lk[lk["transect_id"].isin(geom)]
    for i, k in enumerate(("x0", "y0", "ux", "uy")):
        lk[k] = lk["transect_id"].map(lambda t: geom[t][i])
    lk["number"] = lk["transect_id"].str.rsplit("_", n=1).str[1].astype(int)
    return lk.sort_values("number").reset_index(drop=True)


# Annual median chainage per transect, and the reference and latest medians
def annual_positions(tr: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for t in tr.itertuples():
        obs = load_timeseries(str(timeseries_file(t.transect_id)))
        obs["year"] = obs["date"].dt.year
        obs = obs[obs["year"].between(*YEARS)]
        g = obs.groupby("year")["chainage_m"].agg(["median", "count"])
        g = g[g["count"] >= MIN_OBS_PER_YEAR]
        base = obs.loc[obs["year"].between(*BASE_YEARS), "chainage_m"].median()
        end = obs.loc[obs["year"].between(*END_YEARS), "chainage_m"].median()
        lrr = compute_lrr(filter_dates(obs, f"{YEARS[0]}-01-01", f"{YEARS[1]}-12-31"))
        for y, r in g.iterrows():
            rows.append(dict(transect_id=t.transect_id, year=int(y), n=int(r["count"]),
                             median_chainage_m=r["median"], base_chainage_m=base,
                             end_chainage_m=end, lrr_m_yr=lrr["lrr_m_yr"],
                             lrr_unc95_m_yr=lrr["unc_m_yr"], lrr_n_obs=lrr["n_obs"]))
    out = pd.DataFrame(rows)
    out["change_m"] = out["median_chainage_m"] - out["base_chainage_m"]
    return out


# Alongshore distance of each transect from the groin field, along the reference shoreline
def alongshore(tr: pd.DataFrame, pos: pd.DataFrame, y_groin: float) -> pd.DataFrame:
    base = pos.groupby("transect_id")["base_chainage_m"].first()
    end = pos.groupby("transect_id")["end_chainage_m"].first()
    tr = tr.copy()
    tr["base"] = tr["transect_id"].map(base)
    tr["end"] = tr["transect_id"].map(end)
    for k in ("lrr_m_yr", "lrr_unc95_m_yr", "lrr_n_obs"):
        tr[k] = tr["transect_id"].map(pos.groupby("transect_id")[k].first())
    tr["bx"], tr["by"] = tr.x0 + tr.base * tr.ux, tr.y0 + tr.base * tr.uy
    tr["ex"], tr["ey"] = tr.x0 + tr.end * tr.ux, tr.y0 + tr.end * tr.uy
    s = np.r_[0.0, np.cumsum(np.hypot(np.diff(tr.bx), np.diff(tr.by)))]
    s0 = np.interp(y_groin, tr.by.to_numpy(), s)   # northing rises with transect number here
    tr["dist_m"] = s - s0
    tr["side"] = np.where(tr["dist_m"] < 0, "south", "north")
    return tr


# The Buxton fill as it is entered as a model input
def buxton_fill() -> dict:
    df = pd.read_csv(NOURISHMENT)
    r = df[(df["name"] == FILL_NAME) & df["source"].str.startswith("model input")].iloc[0]
    return dict(year=int(r["year"]), first_gis=int(r["first_gis"]), last_gis=int(r["last_gis"]))


# Esri World Imagery for a UTM window
def basemap(b):
    import contextily as cx
    from rasterio.warp import transform_bounds
    cx.set_cache_dir(str(TILE_CACHE))
    w, s, e, n = transform_bounds(CRS, "EPSG:3857", *b)
    img, ext = cx.bounds2img(w, s, e, n, zoom=TILE_ZOOM, source=cx.providers.Esri.WorldImagery, ll=False)
    return cx.warp_tiles(img, ext, t_crs=CRS)


# (a) the reach on imagery: transects by rate, the two shorelines, the groin field
def draw_map(ax, tr, y_groin):
    x0 = tr[["bx", "ex"]].min().min() - MAP_PAD_X[0]
    x1 = tr[["bx", "ex"]].max().max() + MAP_PAD_X[1]
    y0, y1 = tr.by.min() - 60, tr.by.max() + 60
    img, ext = basemap((x0, y0, x1, y1))
    ax.imshow(img, extent=ext, zorder=0, interpolation="bilinear")
    norm = plt.Normalize(-RATE_CLIM, RATE_CLIM)
    cmap = plt.get_cmap(CMAP)
    for t in tr.itertuples():
        xs = [t.bx - 160 * t.ux, t.bx + 220 * t.ux]
        ys = [t.by - 160 * t.uy, t.by + 220 * t.uy]
        ax.plot(xs, ys, color=cmap(norm(t.lrr_m_yr)), lw=1.0, zorder=2,
                path_effects=[pe.withStroke(linewidth=1.6, foreground="0.1")])
    ax.plot(tr.bx, tr.by, color=C_BASE_LINE, lw=1.3, ls=LS_BASE, zorder=4)
    ax.plot(tr.ex, tr.ey, color=C_END_LINE, lw=1.3, zorder=4)
    ax.plot([np.interp(y_groin, tr.by, tr.bx) + 60, x1], [y_groin, y_groin], color="white", lw=1.0,
            ls=LS_GROIN, zorder=5)
    ax.plot([x1 - 0.04 * (x1 - x0)], [y_groin], marker="<", ms=6, color="white", mec=INK,
            mew=0.6, zorder=6, ls="none")
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    spines_for_image(ax)
    _scalebar(ax, 500.0, show_cells=False)
    _north_arrow(ax, x=0.12, y=0.86)
    ax.text(0.985, 0.008, "Imagery: Esri World Imagery", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=6.5, color=INK_MUTED, zorder=25,
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", boxstyle="square,pad=0.15"))
    return plt.cm.ScalarMappable(norm=norm, cmap=cmap)


# The rate scale for (a), inset over open water in the lower right
def rate_bar(ax, sm):
    ax.add_patch(Rectangle((0.505, 0.115), 0.48, 0.205, transform=ax.transAxes, facecolor="white",
                           alpha=0.88, edgecolor="none", zorder=20))
    cax = ax.inset_axes([0.56, 0.275, 0.37, 0.026], zorder=21)
    cb = plt.colorbar(sm, cax=cax, orientation="horizontal", ticks=[-2, -1, 0, 1, 2])
    cb.ax.tick_params(labelsize=7, length=2, pad=1)
    cb.set_label("Shoreline change\nrate (m/yr)\n← landward | seaward →", fontsize=6.5, labelpad=2)
    cb.outline.set_linewidth(0.5)


# (b) change from the reference position, transect by transect, through time
def draw_hovmoller(ax, tr, pos, fill, y_fill_lo):
    grid = pos.pivot(index="transect_id", columns="year", values="change_m").reindex(tr.transect_id)
    years = np.arange(YEARS[0], YEARS[1] + 1)
    grid = grid.reindex(columns=years)
    d = tr["dist_m"].to_numpy() / 1000.0
    edges_d = np.r_[d[0] - (d[1] - d[0]) / 2, (d[:-1] + d[1:]) / 2, d[-1] + (d[-1] - d[-2]) / 2]
    edges_y = np.r_[years - 0.5, years[-1] + 0.5]
    m = ax.pcolormesh(edges_y, edges_d, np.ma.masked_invalid(grid.to_numpy()), cmap=CMAP,
                      vmin=-CLIM, vmax=CLIM, shading="flat", rasterized=True)
    ax.axhline(0.0, color=INK, lw=0.9, ls=LS_GROIN)
    ax.plot([fill["year"], fill["year"]], [y_fill_lo, edges_d[-1]], **FILL_STYLE)
    ax.set_ylim(edges_d[0], edges_d[-1])
    ax.set_xlim(edges_y[0], edges_y[-1])
    ax.set_ylabel("Alongshore distance\nfrom groin field (km)")
    for side, va, yy in (("north ↑", "bottom", 0.04), ("south ↓", "top", -0.04)):
        ax.text(YEARS[0] + 0.3, yy, side, ha="left", va=va, fontsize=7.5, color=INK,
                bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", boxstyle="square,pad=0.12"))
    return m


# (c) the mean change of each side, with the spread across its transects
def draw_sides(ax, tr, pos, fill):
    p = pos.merge(tr[["transect_id", "side"]], on="transect_id")
    ax.axhline(0.0, color=INK_MUTED, lw=0.6)
    ax.axvline(fill["year"], zorder=0, **FILL_STYLE)
    for side, c in (("south", C_SOUTH), ("north", C_NORTH)):
        g = p[p.side == side].groupby("year")["change_m"]
        med, q1, q3 = g.median(), g.quantile(0.25), g.quantile(0.75)
        ax.fill_between(med.index, q1, q3, color=c, alpha=0.2, lw=0)
        ax.plot(med.index, med.values, color=c, lw=1.4)
    ax.set_xlim(YEARS[0] - 0.5, YEARS[1] + 0.5)
    ax.set_ylabel("Shoreline change since\n1984–1986 (m)")
    ax.set_xlabel("Year")
    ax.grid(axis="y")
    open_frame(ax)


# Run: data, three panels, caption, save
def main() -> None:
    apply_style()
    y_groin, gis_pos = groin_northing()
    tr = transects()
    pos = annual_positions(tr)
    tr = alongshore(tr, pos, y_groin)
    fill = buxton_fill()
    n_s, n_n = int((tr.side == "south").sum()), int((tr.side == "north").sum())

    # Layout: the map on the left, the two time panels stacked on the right
    fig = plt.figure(figsize=figsize("double", 0.62), layout="constrained")
    gs = fig.add_gridspec(2, 2, width_ratios=[0.62, 1.0], height_ratios=[1.15, 1.0])
    ax_map = fig.add_subplot(gs[:, 0])
    ax_h = fig.add_subplot(gs[0, 1])
    ax_s = fig.add_subplot(gs[1, 1], sharex=ax_h)

    # Draw
    sm = draw_map(ax_map, tr, y_groin)
    m = draw_hovmoller(ax_h, tr, pos, fill, 0.0)
    draw_sides(ax_s, tr, pos, fill)
    plt.setp(ax_h.get_xticklabels(), visible=False)
    cb = fig.colorbar(m, ax=ax_h, location="right", shrink=0.95, aspect=18, pad=0.02, extend="both")
    cb.set_label("Change (m)\n← landward | seaward →", fontsize=8)
    cb.outline.set_linewidth(0.5)
    rate_bar(ax_map, sm)
    title(ax_map, 0, "Buxton, GIS 3–8")
    title(ax_h, 1, "Each transect")
    title(ax_s, 2, "Each side")

    # One legend for the whole figure
    handles = [Line2D([], [], color=C_SOUTH, lw=1.4), Line2D([], [], color=C_NORTH, lw=1.4),
               Line2D([], [], color=C_BASE_LINE, lw=1.3, ls=LS_BASE),
               Line2D([], [], color=C_END_LINE, lw=1.3),
               Line2D([], [], color=INK, lw=0.9, ls=LS_GROIN),
               Line2D([], [], **FILL_STYLE)]
    labels = ["South of groin field", "North of groin field",
              f"Shoreline {BASE_YEARS[0]}–{BASE_YEARS[1]}", f"Shoreline {END_YEARS[0]}–{END_YEARS[1]}",
              "Groin field", f"Beach fill {fill['year']} (GIS {fill['first_gis']}–{fill['last_gis']})"]
    fig.legend(handles, labels, loc="outside lower center", ncol=3, frameon=False)

    # Table of what is drawn
    table = pos.merge(tr[["transect_id", "dist_m", "side"]], on="transect_id")
    table.round(2).to_csv(support_dir(OUT) / f"{STEM}_annual.csv", index=False)
    rates = tr[["transect_id", "domain_number", "side", "dist_m", "lrr_m_yr", "lrr_unc95_m_yr",
                "lrr_n_obs"]]
    rates.round(4).to_csv(support_dir(OUT) / f"{STEM}_lrr.csv", index=False)
    med_s = tr.loc[tr.side == "south", "lrr_m_yr"].median()
    med_n = tr.loc[tr.side == "north", "lrr_m_yr"].median()
    print(f"LRR {tr.lrr_m_yr.min():.2f} to {tr.lrr_m_yr.max():.2f} m/yr; "
          f"median south {med_s:.2f}, north {med_n:.2f}")

    caption(fig, (
        f"Observed shoreline position near the Buxton groin field, Hatteras Island, from the CoastSat "
        f"satellite-derived shoreline record ({YEARS[0]}–{YEARS[1]}). (a) The reach (GIS domains "
        f"{GIS_RANGE[0]}–{GIS_RANGE[1]}) on Esri World Imagery, north up: the {len(tr)} CoastSat "
        f"transects used, {n_s} south and {n_n} north of the groin field, each coloured by its "
        f"shoreline change rate: the linear regression rate (ordinary least-squares slope of all "
        f"CoastSat positions against time, 1 January {YEARS[0]} to 31 December {YEARS[1]}, m/yr, "
        f"positive seaward), red landward and blue seaward, the scale clipped at ±{RATE_CLIM:g} m/yr "
        f"(no transect reaches it); the "
        f"shoreline as the median of each transect's positions over {BASE_YEARS[0]}–{BASE_YEARS[1]} "
        f"(dashed charcoal) and {END_YEARS[0]}–{END_YEARS[1]} (solid near-black); the dash-dot "
        f"white line and arrowhead mark "
        f"the groin field at the model's groin position (GIS {gis_pos}, the boundary between "
        f"domains 5 and 6). (b) Change in each transect's annual median position from its "
        f"{BASE_YEARS[0]}–{BASE_YEARS[1]} median, along the transect, positive seaward, by alongshore "
        f"distance from the groin field measured along the {BASE_YEARS[0]}–{BASE_YEARS[1]} "
        f"shoreline (positive north), red landward and blue seaward, white near zero; the dash-dot "
        f"line is the groin field; years with fewer than {MIN_OBS_PER_YEAR} satellite positions "
        f"are blank; the colour scale is clipped at ±{CLIM:.0f} m. The dotted vertical line marks the "
        f"{fill['year']} Buxton beach fill (GIS {fill['first_gis']}–{fill['last_gis']}, from the "
        f"management record), which begins at the groin field and extends north past this view. "
        f"(c) The median of (b) over the transects on each side, south of the groin field "
        f"green and north amber, with the interquartile range across them shaded; the dotted "
        f"line is the fill year. Positions are CoastSat's chainages "
        f"along each transect, not referenced to a vertical or survey datum beyond CoastSat's own "
        f"processing; annual medians still carry tide, wave and season noise of several metres."))
    save(fig, OUT / STEM, close=True)


if __name__ == "__main__":
    main()
