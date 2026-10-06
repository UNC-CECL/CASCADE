"""
Where each nourishment report says the sand went, drawn against the model domains it was mapped to.

    python scripts/input_prep/4-mgmt-forcings/nourishment_reported_extent.py

Places every limit a source names (a groin, a street, a refuge boundary, a distance
from a pier) on the 2009 period's mean CoastSat shoreline (2008-08-17 to 2010-08-17),
reads off the GIS domain it
falls in, and draws the reported stretch beside the model footprint on imagery.
Writes reported_limits.csv and one six-panel figure. A check on the forcing; nothing
here feeds a run.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-06
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(_HERE.parent))

import geopandas as gpd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from shapely.geometry import Point  # noqa: E402

from beach_nourishment import COMMUNITY_COLOUR, GROIN_FILE, RECORD_ONLY_PROJECTS  # noqa: E402
from site_layer.hat_figure_style import (C, INK, MAP_TEXT, apply_style, figsize, letter_at,  # noqa: E402
                                         north_dart, record_caption, save, scale_bar_km)
from site_layer.hat_observed_rates import DOMAIN_BOXES, mean_shoreline_geojson  # noqa: E402
from site_layer.hat_topo_version import MGMT_ROOT, shoreline_window_for_year  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS, HATTERAS_NOURISHMENT_PROJECTS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT_DIR = MGMT_ROOT / "nourishment" / "reported_extent"
# The 2009 period's shoreline window (+/-1 yr of the 2009 USACE lidar), the one its island offset uses
SHORELINE_WINDOW = shoreline_window_for_year(2009)
SHORELINE = mean_shoreline_geojson(*SHORELINE_WINDOW)
MAP_CRS = "EPSG:26918"
MI = 1609.344
FT = 0.3048
PAD_DOMAINS = 1.5            # domains of context beyond the footprint and the reported limits
PANEL_ASPECT = 2.0           # panel height / width
TILE_ZOOM = 15
TILE_CACHE = Path(tempfile.gettempdir()) / "hat_tile_cache"
BRACKET_GAP_M = 120          # the reported-extent bracket sits this far seaward of the domain boxes

# Located points (lat, lon) and where each came from
REFUGE_SOUTH = (35.6064617, -75.4706281)          # Pea Island NWR, southernmost boundary vertex, OpenStreetMap
HAULOVER = (35.30044634138652, -75.51508959896135)  # Haulover Day Use Area start, picked by Hannah 2026-10-04
DUE_EAST = (35.35694, -75.50126)                  # Due East (road), OpenStreetMap way centre
ASKINS_N = (35.32564, -75.50802)                  # Askins Creek North Drive, OpenStreetMap way centre
GREENWOOD = (35.32875, -75.50710)                 # Greenwood Place, OpenStreetMap way centre
PAMPAS = (35.33797, -75.50517)                    # Pampas Street, OpenStreetMap way centre

SRC = {
    "rodanthe": "https://beachapedia.org/State_of_the_Beach/State_Reports/NC/Beach_Fill (quoting the 2013 USACE public notice, via Island Free Press)",
    "buxton17": "https://outerbanksvoice.com/2018/03/01/delayed-buxton-beach-nourishment-project-is-finally-done/ ; https://islandfreepress.org/blog/decades-of-shoreline-engineering-the-long-history-of-a-changing-buxton-beach/",
    "buxton22": "https://content.govdelivery.com/accounts/NCDARECOUNTY/bulletins/3285a2f",
    "avon22": "https://www.darenc.gov/government/beach-nourishment/avon-beach-nourishement ; https://dredgewire.com/2022-avon-beach-nourishment/",
    "avon26": "https://islandfreepress.org/outer-banks-news/more-than-75-of-avon-nourishment-project-completed-as-work-moves-south/ ; https://islandfreepress.org/blog/avon-and-buxton-beach-nourishment-faqs-2026-edition/",
    "buxton26": "https://islandfreepress.org/blog/avon-and-buxton-beach-nourishment-faqs-2026-edition/",
}
# -----------------------------------------------------------------------------


class Shore:
    """The mean shoreline as one north-running line, with arc length and domain lookups."""

    def __init__(self):
        g = gpd.read_file(SHORELINE).to_crs(MAP_CRS)
        line = max(g.geometry.explode(index_parts=False), key=lambda x: x.length)
        if line.coords[0][1] > line.coords[-1][1]:
            line = type(line)(list(line.coords)[::-1])
        self.line = line
        self.dom = gpd.read_file(DOMAIN_BOXES).to_crs(MAP_CRS).set_index("domain_id").sort_index()

    # Point on the shoreline nearest a lat/lon, as arc length
    def s_of(self, latlon):
        p = gpd.GeoSeries([Point(latlon[1], latlon[0])], crs=4326).to_crs(MAP_CRS).iloc[0]
        return self.line.project(p)

    def xy(self, s):
        q = self.line.interpolate(s)
        return q.x, q.y

    # GIS domain containing a shoreline point, by northing, and metres into it from its south edge
    def domain_of(self, s):
        x, y = self.xy(s)
        for d, row in self.dom.iterrows():
            b = row.geometry.bounds
            if b[1] <= y < b[3]:
                return int(d), y - b[1]
        return None, np.nan

    # Shoreline point level with a lat/lon (same northing), for a limit picked off the island rather than the beach
    def s_at_northing(self, latlon):
        y = gpd.GeoSeries([Point(latlon[1], latlon[0])], crs=4326).to_crs(MAP_CRS).iloc[0].y
        ss = np.arange(0, self.line.length, 2.0)
        ys = np.array([self.xy(v)[1] for v in ss])
        return float(ss[np.argmin(np.abs(ys - y))])

    # Shoreline point at the pier: domain south edge plus the configured fraction, by northing
    def s_of_pier(self, name):
        d, frac = HATTERAS_ANNOTATIONS.piers[name]
        b = self.dom.loc[d].geometry.bounds
        y = b[1] + frac * (b[3] - b[1])
        ys = np.array([self.xy(s)[1] for s in np.arange(0, self.line.length, 10.0)])
        return float(np.argmin(np.abs(ys - y)) * 10.0)


# Each project's reported limits: (label on the map, reported wording, how located, arc length)
def reported(shore):
    groin = gpd.read_file(GROIN_FILE).to_crs(MAP_CRS)
    south_groin = groin.geometry.iloc[0]
    s_groin = shore.line.project(south_groin.centroid)
    s_haul = shore.s_at_northing(HAULOVER)
    s_refuge = shore.s_of(REFUGE_SOUTH)
    s_rod_n = s_refuge + 1.5 * MI
    s_avon_pier = shore.s_of_pier("Avon Pier")
    buxton_limits = [
        ("southernmost groin", "the groin at the old Cape Hatteras Lighthouse site", "digitized groin (groins_hatteras.geojson)", s_groin),
        ("Haulover Day Use Area", "Haulover Day Use Area", "point picked by Hannah on imagery; the limit is at its northing", s_haul)]
    return {
        ("Rodanthe emergency fill", 2014): dict(
            key="rodanthe", length_mi=2.13, limits=[
                ("2.13 mi south of the north end", "2.13 miles of beach ... into the Mirlo Beach community to just north of the Rodanthe pier",
                 "measured: north end minus 2.13 mi along the shoreline", s_rod_n - 2.13 * MI),
                ("1.5 mi north of the refuge border", "from 1.5 miles north of the Pea Island National Wildlife Refuge border",
                 "measured: OpenStreetMap refuge boundary plus 1.5 mi along the shoreline", s_rod_n)],
            marks=[("refuge border", s_refuge), ("Rodanthe Pier", shore.s_of_pier("Rodanthe Pier"))]),
        ("Buxton beach nourishment", 2017): dict(key="buxton17", length_mi=2.94, limits=buxton_limits, marks=[]),
        ("Buxton shore protection", 2022): dict(key="buxton22", length_mi=2.9, limits=buxton_limits, marks=[]),
        ("Avon shore protection", 2022): dict(
            key="avon22", length_mi=2.5, limits=[
                ("NPS / Avon boundary", "to the National Park Service / Avon boundary (south village limit, just north of ORV Ramp 38)",
                 "estimated: Askins Creek North Drive (OpenStreetMap), the config's south limit", shore.s_of(ASKINS_N)),
                ("Due East Rd", "3,000 feet north of Avon Pier at Due East Road", "geocoded: Due East (OpenStreetMap)",
                 shore.s_of(DUE_EAST))],
            marks=[("Avon Pier", s_avon_pier)]),
        ("Avon 2026", 2026): dict(
            key="avon26", length_mi=1.0, limits=[
                ("Greenwood Place", "southward toward Greenwood Place", "geocoded: Greenwood Place (OpenStreetMap)",
                 shore.s_of(GREENWOOD)),
                ("just south of Avon Pier", "from just south of the Avon Fishing Pier (progress report, mid-project: crews \"shifted operations south from the Pampas Drive area\" once the northern section was done)",
                 "measured: Avon Pier position (site config) less 100 m", s_avon_pier - 100.0)],
            marks=[("Pampas St", shore.s_of(PAMPAS)), ("Avon Pier", s_avon_pier)]),
        ("Buxton 2026", 2026): dict(
            key="buxton26", length_mi=2.9,
            limits=[(a, "the southernmost groin in Buxton" if i == 0 else "the Haulover area", c, s) for i, (a, _, c, s) in
                    enumerate(buxton_limits)], marks=[]),
    }


def main():
    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    shore = Shore()
    rep = reported(shore)
    model = {(p.name, p.year): (p, True) for p in HATTERAS_NOURISHMENT_PROJECTS}
    model.update({(p.name, p.year): (p, False) for p in RECORD_ONLY_PROJECTS})

    rows = []
    for k, r in rep.items():
        p, in_model = model[k]
        first, last = min(p.gis_domains), max(p.gis_domains)
        s0, s1 = sorted(x[3] for x in r["limits"])
        for (label, words, how, s), end in zip(sorted(r["limits"], key=lambda t: t[3]), ("south", "north")):
            d, into = shore.domain_of(s)
            lon, lat = gpd.GeoSeries([Point(*shore.xy(s))], crs=MAP_CRS).to_crs(4326).iloc[0].coords[0]
            rows.append(dict(project=k[0], year=k[1], in_model=in_model, limit=end, label=label, reported_wording=words,
                             source=SRC[r["key"]], lat=round(lat, 6), lon=round(lon, 6), how_located=how,
                             gis_domain=d, m_into_domain_from_south=round(into),
                             reported_length_mi=r["length_mi"], located_length_km=round((s1 - s0) / 1000, 2),
                             model_first_gis=first, model_last_gis=last,
                             model_length_km=len(p.gis_domains) * 0.5))
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT_DIR / "reported_limits.csv", index=False)

    import contextily as cx
    from rasterio.warp import transform_bounds
    cx.set_cache_dir(str(TILE_CACHE))
    keys = sorted(rep, key=lambda k: (k[1], k[0]))
    keys = [keys[i] for i in (0, 1, 2, 3, 4, 5)]
    fig, axes = plt.subplots(2, 3, figsize=figsize("double", height=9.6))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.07, wspace=0.05, hspace=0.17)
    dom = shore.dom
    cen = dom.geometry.centroid
    for i, (ax, k) in enumerate(zip(axes.flat, keys)):
        r = rep[k]
        p, in_model = model[k]
        town = k[0].split()[0]
        dark, light = COMMUNITY_COLOUR[town]
        first, last = min(p.gis_domains), max(p.gis_domains)
        sub = tab[(tab.project == k[0]) & (tab.year == k[1])]
        lo_d = int(min(first, sub.gis_domain.min()) - PAD_DOMAINS)
        hi_d = int(max(last, sub.gis_domain.max()) + PAD_DOMAINS)
        lo_d, hi_d = max(1, lo_d), min(90, hi_d)
        x0, y0, x1, y1 = dom.loc[lo_d:hi_d].total_bounds
        x1 += 2 * BRACKET_GAP_M
        cxm, cym = (x0 + x1) / 2, (y0 + y1) / 2
        hh = max(y1 - y0, (x1 - x0) * PANEL_ASPECT) / 2 * 1.02
        bx = (cxm - hh / PANEL_ASPECT, cym - hh, cxm + hh / PANEL_ASPECT, cym + hh)
        w, s_, e, n = transform_bounds(MAP_CRS, "EPSG:3857", *bx)
        img, ext = cx.bounds2img(w, s_, e, n, zoom=TILE_ZOOM, source=cx.providers.Esri.WorldImagery, ll=False)
        img, ext = cx.warp_tiles(img, ext, t_crs=MAP_CRS)
        ax.imshow(img, extent=ext, zorder=0)
        dom.loc[lo_d:hi_d].plot(ax=ax, facecolor="none", edgecolor="white", lw=0.4, alpha=0.75, zorder=2)
        dom.loc[first:last].plot(ax=ax, facecolor=light, edgecolor=dark, lw=0.9, alpha=0.4, zorder=3)
        for d in range(lo_d, hi_d + 1):
            if d in dom.index:
                b = dom.loc[d].geometry.bounds
                ax.text(b[0] + 70, cen.loc[d].y, str(d), ha="left", va="center",
                        **{**MAP_TEXT, "fontsize": 6.5, "fontweight": "bold" if first <= d <= last else "normal"})
        # Limits are alongshore positions: each is a line across its domain box at its northing,
        # and the reported stretch a bracket just seaward of the boxes
        s0, s1 = sorted(x[3] for x in r["limits"])
        y_lo, y_hi = shore.xy(s0)[1], shore.xy(s1)[1]
        xb = dom.loc[lo_d:hi_d].total_bounds[2] + BRACKET_GAP_M
        for xx, lw, col in ((xb, 5.0, "white"), (xb, 3.0, INK)):
            ax.plot([xx, xx], [y_lo, y_hi], color=col, lw=lw, solid_capstyle="butt", zorder=6 if col == "white" else 7)
        for (label, _, how, sv) in r["limits"]:
            y = shore.xy(sv)[1]
            d, _ = shore.domain_of(sv)
            b = dom.loc[d].geometry.bounds
            ax.plot([b[0], xb], [y, y], color="white", lw=2.6, zorder=6)
            ax.plot([b[0], xb], [y, y], color=INK, lw=1.2, ls=(0, (4, 2)), zorder=7)
            # Labels sit inside the bracket: below the north limit, above the south one
            north = sv == s1
            ax.text(b[0] + 70, y + (-1 if north else 1) * 0.010 * (bx[3] - bx[1]), label, ha="left",
                    va="top" if north else "bottom", **{**MAP_TEXT, "fontsize": 7, "fontweight": "bold"})
        for label, sv in r["marks"]:
            y = shore.xy(sv)[1]
            d, _ = shore.domain_of(sv)
            px = dom.loc[d].geometry.bounds[2] - 200
            ax.plot(px, y, marker="v" if "Pier" in label else "s", ms=5, color="white", mec=INK, mew=0.6, zorder=8)
            ax.text(px - 90, y, label, ha="right", va="center",
                    **{**MAP_TEXT, "fontsize": 6.5, "fontstyle": "italic"})
        if town == "Buxton":
            hp = gpd.GeoSeries([Point(HAULOVER[1], HAULOVER[0])], crs=4326).to_crs(MAP_CRS).iloc[0]
            ax.plot(hp.x, hp.y, marker="o", ms=5, color="white", mec=INK, mew=0.8, zorder=9)
            g = gpd.read_file(GROIN_FILE).to_crs(MAP_CRS)
            g.plot(ax=ax, color="white", lw=3.0, zorder=7)
            g.plot(ax=ax, color=C["GROIN"], lw=1.6, zorder=7)
        ax.set_xlim(bx[0], bx[2])
        ax.set_ylim(bx[1], bx[3])
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.6)
        letter_at(ax, i, 0.03, 0.985)
        ax.text(0.03, 0.935, f"{town}, {k[1]}" + ("" if in_model else "\nnot in the model"), transform=ax.transAxes,
                ha="left", va="top", **{**MAP_TEXT, "fontsize": 8.5, "fontstyle": "italic"})
        scale_bar_km(ax, length_m=1000, segments=2, x=0.06, y=0.05)
        north_dart(ax, (bx[2] - 0.13 * (bx[2] - bx[0]), bx[1] + 0.08 * (bx[3] - bx[1])),
                   arrow_m=0.05 * (bx[3] - bx[1]))
        so, no = sub.iloc[0], sub.iloc[1]
        ax.text(0.5, -0.02,
                f"Reported {r['length_mi']:g} mi ({r['length_mi'] * MI / 1000:.1f} km); limits as located "
                f"{so.located_length_km:.1f} km\nlimits in GIS {so.gis_domain} and {no.gis_domain}; "
                f"model GIS {first}–{last} ({len(p.gis_domains) * 0.5:.1f} km)",
                transform=ax.transAxes, ha="center", va="top", fontsize=7, color=INK, linespacing=1.35)
    handles = [Line2D([], [], color=INK, lw=3.0, label="Reported extent"),
               Line2D([], [], color=INK, lw=1.2, ls=(0, (4, 2)), label="Reported limit"),
               Patch(facecolor="0.85", edgecolor="0.4", lw=0.8, alpha=0.8, label="Model footprint (community colour)"),
               Line2D([], [], color=C["GROIN"], lw=1.6, label="Buxton groins"),
               Line2D([], [], color="white", marker="v", mec=INK, ls="none", ms=5, label="Pier")]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=7.5, bbox_to_anchor=(0.5, 0.0))
    out = save(fig, OUT_DIR / "nourishment_reported_extent", close=True)[0]
    record_caption(out, (
        "Reported nourishment extents against the model footprints, on Esri World Imagery (current), north up. "
        f"Each limit a source names is placed on the {SHORELINE_WINDOW[0]} to {SHORELINE_WINDOW[1]} mean CoastSat shoreline and drawn as a dashed line "
        "across its domain at that northing; the black bracket beside the domains spans the reported stretch. "
        "Domains are assigned by northing, so only the alongshore position matters. Shaded boxes are the 500 m model "
        "domains the fill is applied to (2026 panels: record only, not in the model). Limits are geocoded "
        "(OpenStreetMap street or boundary), measured (a reported distance along the shoreline from a located "
        "point), or estimated; reported_limits.csv gives the wording, source and method for each, and the "
        "domain it falls in. Rodanthe 2014's limits come from the 2013 permit notice (2.13 mi from 1.5 mi "
        "north of the Pea Island refuge border); the as-built project was reported as about 2 mi. Avon 2026's "
        "north limit is the pier less 100 m; the progress report also names Pampas Street, marked."))
    with pd.option_context("display.width", 220, "display.max_columns", 30, "display.max_colwidth", 40):
        print(tab[["project", "year", "limit", "label", "gis_domain", "m_into_domain_from_south",
                   "reported_length_mi", "located_length_km", "model_first_gis", "model_last_gis"]].to_string(index=False))
    print(out)


if __name__ == "__main__":
    main()
