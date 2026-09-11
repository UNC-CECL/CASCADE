"""
HAT_overwash_vs_footprint.py
==============================================================================
Does the observed overwash of 1984–1997 line up with the rows the 1984
reconstruction adds and removes?

    python HAT_overwash_vs_footprint.py

THE TWO RECORDS
    The footprint (2-domain-reconstruction-1984/2-extent/footprint_1984_by_domain.csv)
    gives every domain the median shift between the digitised 1984 and 1997
    dune lines, positive where the 1984 line lay seaward, and the row count the
    10 m rule keeps of it: rows ADDED where the dune retreated over the window,
    rows REMOVED where it advanced.

    The overwash record (observations/Hatteras_Overwash_Data.xlsx) gives every
    domain, per image, whether washover was visible.

THE WINDOW (Hannah, 2026-09-10: strictly between the line dates)
    Both dune lines were digitised from the same USGS photographs the overwash
    record uses: the 1984 line from the 19 Sep 1984 frame, the 1997 line from
    the 12 Oct 1997 frame. So the images that can speak to the shift are the
    ones taken AFTER the 1984 frame and UP TO AND INCLUDING the 1997 frame:
    nine images, Aug 1985 to Oct 1997. The Diana image is out (the 1984 line
    was drawn on it, so that overwash predates the line); the Bonnie image is
    out (it postdates the 1997 line).

WHAT COUNTS AS AGREEMENT (Hannah: show both readings, do not blame absence)
    Physically, overwash flattens the dune and pushes the vegetation break
    landward, so overwash should go with rows ADDED. Overwash in a rows-REMOVED
    domain is the disagreement worth a look. A rows-added domain with NO
    overwash in the window is listed as unexplained, not as a mismatch: the
    nine images are two to four years apart and washover fades from imagery
    within a few years, so absence is weak evidence. Both readings are
    reported: "given overwash, which action?" and "given the action, was
    there overwash?".

FLAGS (kept, not dropped)
    SPREAD_STRADDLES_ZERO from the footprint (the p10–p90 of the shift
    crosses zero); the erosion hotspots and jetties from domains.geojson;
    NC-12 relocated 1984–2004 (road_relocation_1984_2004.csv); and shoreline
    erosion faster than ERODE_THRESH by the CoastSat 1984–2004 LRR, with the
    DSAS 1978–1997 mean rate carried beside it because it sits closer to the
    window. Erosion retreats the dune line without any overwash, which is
    why the last flag exists.

OUTPUT   data/hatteras_init/8-overwash-analysis/vs-footprint/
    overwash_vs_footprint_by_domain.csv    the joined table, one row per domain
    overwash_vs_footprint_contingency.csv  overwashed x action, all and unflagged
    overwash_vs_footprint_summary.txt      the readings in words, with the lists
    ../figures/vs-footprint/
        overwash_vs_footprint_alongshore.png   images, footprint bars, flags, by domain
        overwash_vs_footprint_summary.png      shift by overwash status; share overwashed per action
        overwash_vs_footprint_map.png          three alongshore sections, zoomed, each with
                                               the three layers: overwash, footprint, reading
        overwash_vs_footprint_map_island.png   the whole island, the same three layers
    The map needs D:/Hatteras_GIS (domain boxes, coastline), through
    overwash_map_periods.load_geometry.
==============================================================================
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(HERE))

from hat_figure_style import (C, DOMAIN_AXIS_LABEL, INK, _title,     # noqa: E402
                              apply_style, figsize, open_frame, save,
                              spines_for_image, town_bands)
from overwash_data import (OUT_DIR, SECTIONS, load_observations,     # noqa: E402
                           remove_caption, upsert_caption)
from overwash_heatmap_multiperiod import CLR_PART_EDGE, CLR_PART_FACE     # noqa: E402
from overwash_map_periods import (CLR_BOX, CLR_LAND, CLR_ROAD, PAD_E,     # noqa: E402
                                  PAD_N, PAD_S, PAD_W, count_cmap,
                                  draw_island, load_geometry,
                                  scalebar_and_north)

INIT = REPO / "data" / "hatteras_init"
FOOTPRINT = (INIT / "1-barrier3d-domains/1984-start/2-domain-reconstruction-1984"
             / "2-extent/footprint_1984_by_domain.csv")
RELOCATION = INIT / "4-mgmt-forcing/road_relocation/1984_2004/road_relocation_1984_2004.csv"
COASTSAT = (INIT / "6-scr-smooth/HAT_loess_method_comparison_output/03_cascade_inputs"
            / "cascade_lrr_inputs_transect_based.csv")
DSAS = INIT / "5-scr/scr-dsas-1978-2019/dsas_1978_1997_domain_means.csv"
DOMAIN_FILE = Path("D:/Hatteras_GIS/domains.geojson")

VS_DIR = OUT_DIR / "vs-footprint"
FIG_DIR = OUT_DIR / "figures" / "vs-footprint"

LINE_1984 = pd.Timestamp("1984-09-19")     # the frame the 1984 dune line was drawn on
LINE_1997 = pd.Timestamp("1997-10-12")     # the frame the 1997 dune line was drawn on
ERODE_THRESH = -2.0                        # m/yr, CoastSat LRR 1984–2004, "eroding fast"

# Fallback if the drive is absent: what domains.geojson says (read 2026-09-10).
HOTSPOT_FALLBACK = {10: "Buxton", 11: "Buxton", 12: "Buxton", 13: "Buxton",
                    14: "Buxton", 15: "Buxton", 16: "Buxton",
                    32: "Avon", 33: "Avon", 34: "Avon", 35: "Avon", 36: "Avon"}
# (the geojson also names "Rodanthe 'S' Curves"; read from the drive when present)
ARMOR_FALLBACK = {6: "3 jetties, from 1930s and 1970s"}

# The footprint bars take the LIGHTER RdBu pair so they cannot be read as the
# overwash accent red of panel (a); the overwash marker in (b) is the accent.
CLR_ADD = "#d6604d"       # rows added (1984 line seaward)
CLR_REMOVE = "#4393c3"    # rows removed
CLR_NONE = "#9a9a9a"
CLR_OW = C["ACCENT"]
CLR_INK = INK
ACTION_CLR = {"add": CLR_ADD, "remove": CLR_REMOVE, "none": CLR_NONE}


# =================================================================== inputs
def load_inputs():
    obs, domains, matrix = load_observations()
    win = (obs["Imagery_Date"] > LINE_1984) & (obs["Imagery_Date"] <= LINE_1997)
    idx = [int(i) for i in obs.index[win]]
    print(f"  window: {len(idx)} images, "
          f"{obs.loc[idx[0], 'Imagery_Date']:%Y-%m-%d} to "
          f"{obs.loc[idx[-1], 'Imagery_Date']:%Y-%m-%d}")

    fp = pd.read_csv(FOOTPRINT)
    fp = fp[["domain", "topo_version", "shift_m_median", "shift_m_p10",
             "shift_m_p90", "spread_straddles_zero", "n_cells", "action", "flags"]]

    rel = pd.read_csv(RELOCATION)[["domain", "classification",
                                   "median_signed_landward_m"]]
    cs = pd.read_csv(COASTSAT)[["domain", "cs_lrr_1984_2004"]]
    ds = pd.read_csv(DSAS)[["domain_id", "annual_rate_m_per_yr"]].rename(
        columns={"domain_id": "domain", "annual_rate_m_per_yr": "dsas_rate_1978_1997_m_yr"})

    hotspot, armor = HOTSPOT_FALLBACK, ARMOR_FALLBACK
    if DOMAIN_FILE.exists():
        import geopandas as gpd
        g = gpd.read_file(DOMAIN_FILE)
        # `hotspot` is the 0/1 column; `eros_hot` carries the name ("5 - Buxton")
        hot = g[g["hotspot"].astype(int) == 1]
        hotspot = {int(r.domain_id): str(r.eros_hot).split("-", 1)[-1].strip()
                   for r in hot.itertuples()}
        arm = g[g["armor"].astype(str).str.strip().ne("")]
        armor = {int(r.domain_id): str(r.armor) for r in arm.itertuples()}
    else:
        print(f"  WARNING: {DOMAIN_FILE} not reachable; hotspot and armor "
              f"flags from the fallback lists in this script")
    return obs, domains, matrix, idx, fp, rel, cs, ds, hotspot, armor


def build_table(obs, domains, matrix, idx, fp, rel, cs, ds, hotspot, armor):
    sub = matrix[idx]
    labels = [f"{obs.loc[i, 'Imagery_Date']:%Y-%m}" for i in idx]
    rows = []
    for j, d in enumerate(domains):
        col = sub[:, j]
        seen = [labels[k] for k in range(len(idx)) if col[k] >= 0.5]
        rows.append(dict(domain=int(d),
                         n_images_assessed=int(np.sum(~np.isnan(col))),
                         n_images_overwash=len(seen),
                         overwashed=len(seen) > 0,
                         overwash_images="; ".join(seen)))
    t = pd.DataFrame(rows).merge(fp, on="domain", how="left")
    t = t.merge(rel, on="domain", how="left").merge(cs, on="domain", how="left")
    t = t.merge(ds, on="domain", how="left")
    t["road_relocated_1984_2004"] = t["classification"].eq("relocated")
    t["erosion_hotspot"] = t["domain"].map(hotspot).fillna("")
    t["armor"] = t["domain"].map(armor).fillna("")
    t["eroding_fast"] = t["cs_lrr_1984_2004"] < ERODE_THRESH
    t["spread_straddles_zero"] = t["spread_straddles_zero"].astype(bool)
    t["unflagged"] = ~(t["spread_straddles_zero"] | t["road_relocated_1984_2004"]
                       | t["erosion_hotspot"].ne("") | t["armor"].ne(""))

    def reading(r):
        if r["overwashed"] and r["action"] == "add":
            return "agrees (overwash, rows added)"
        if r["overwashed"] and r["action"] == "remove":
            return "DISAGREES (overwash, rows removed)"
        if r["overwashed"]:
            return "overwash, no rows changed"
        if r["action"] == "add":
            return "unexplained retreat (rows added, no overwash seen)"
        return "no overwash seen"
    t["reading"] = t.apply(reading, axis=1)
    order = ["domain", "action", "n_cells", "shift_m_median", "shift_m_p10",
             "shift_m_p90", "spread_straddles_zero", "overwashed",
             "n_images_overwash", "n_images_assessed", "overwash_images",
             "reading", "road_relocated_1984_2004", "median_signed_landward_m",
             "erosion_hotspot", "armor", "cs_lrr_1984_2004", "eroding_fast",
             "dsas_rate_1978_1997_m_yr", "unflagged", "topo_version", "flags"]
    return t[order]


# ================================================================= readings
def contingency(t, label):
    ct = pd.crosstab(t["overwashed"].map({True: "overwashed", False: "not seen"}),
                     t["action"]).reindex(index=["overwashed", "not seen"],
                                          columns=["add", "none", "remove"]).fillna(0).astype(int)
    ct.index.name = f"{label} (n={len(t)})"
    return ct


def fisher(t):
    """overwashed x (add vs remove), 'none' left out. Small counts: report, do not lean on."""
    try:
        from scipy.stats import fisher_exact
    except ImportError:
        return None
    s = t[t["action"].isin(["add", "remove"])]
    a = int(((s.action == "add") & s.overwashed).sum())
    b = int(((s.action == "add") & ~s.overwashed).sum())
    c = int(((s.action == "remove") & s.overwashed).sum())
    d = int(((s.action == "remove") & ~s.overwashed).sum())
    odds, p = fisher_exact([[a, b], [c, d]])
    return dict(a=a, b=b, c=c, d=d, odds=odds, p=p)


def summary_text(t, idx, obs):
    lines = []
    dates = [f"{obs.loc[i, 'Imagery_Date']:%d %b %Y}" for i in idx]
    lines.append("Overwash 1984-1997 against the 1984 footprint (rows added / removed)")
    lines.append("=" * 72)
    lines.append(f"Window: {len(idx)} images strictly between the dune-line frames "
                 f"({LINE_1984:%d %b %Y} excluded, {LINE_1997:%d %b %Y} included):")
    lines.append("  " + ", ".join(dates))
    lines.append("")
    for label, s in (("all domains", t), ("unflagged only", t[t.unflagged])):
        ct = contingency(s, label)
        lines.append(ct.to_string())
        lines.append("")
        ow = s[s.overwashed]
        lines.append(f"Reading 1, given overwash ({len(ow)} domains): "
                     + ", ".join(f"{k} {int((ow.action == k).sum())}"
                                 for k in ("add", "none", "remove")))
        parts = []
        for k in ("add", "none", "remove"):
            g = s[s.action == k]
            parts.append(f"{k} {int(g.overwashed.sum())}/{len(g)}")
        lines.append("Reading 2, given the action, overwashed: " + ", ".join(parts))
        f = fisher(s)
        if f:
            lines.append(f"Fisher exact, overwashed x (add vs remove): "
                         f"[[{f['a']}, {f['b']}], [{f['c']}, {f['d']}]]  "
                         f"odds {f['odds']:.2f}, p = {f['p']:.3f}  (small counts)")
        lines.append("")
    lines.append("Lists")
    lines.append("-" * 72)

    def show(sub):
        out = []
        for r in sub.itertuples():
            fl = []
            if r.spread_straddles_zero:
                fl.append("shift straddles 0")
            if r.road_relocated_1984_2004:
                fl.append("NC-12 relocated")
            if r.erosion_hotspot:
                fl.append(f"hotspot {r.erosion_hotspot}")
            if r.armor:
                fl.append("armor")
            if r.eroding_fast:
                fl.append(f"CoastSat {r.cs_lrr_1984_2004:+.1f} m/yr")
            out.append(f"  D{r.domain:>2}  shift {r.shift_m_median:+6.1f} m, "
                       f"{r.n_cells:+d} rows;  overwash in {r.overwash_images or '-'}"
                       + (f";  [{', '.join(fl)}]" if fl else ""))
        return out or ["  (none)"]
    for title, mask in (
        ("Agrees: overwash seen, rows added", t.overwashed & (t.action == "add")),
        ("DISAGREES: overwash seen, rows removed", t.overwashed & (t.action == "remove")),
        ("Overwash seen, no rows changed", t.overwashed & (t.action == "none")),
        ("Unexplained retreat: rows added, no overwash seen (absence is weak evidence)",
         ~t.overwashed & (t.action == "add")),
    ):
        lines.append(title)
        lines += show(t[mask])
        lines.append("")
    lines.append(f"Flags: shift straddles zero {int(t.spread_straddles_zero.sum())}, "
                 f"NC-12 relocated {int(t.road_relocated_1984_2004.sum())}, "
                 f"erosion hotspot {int(t.erosion_hotspot.ne('').sum())}, "
                 f"armor {int(t.armor.ne('').sum())}, "
                 f"eroding faster than {ERODE_THRESH:+.0f} m/yr {int(t.eroding_fast.sum())}; "
                 f"unflagged {int(t.unflagged.sum())} of {len(t)}.")
    return "\n".join(lines) + "\n"


# =================================================================== figures
def write(fig, out):
    """PNG for slides and a PDF beside it for the manuscript (house `save`)."""
    save(fig, out, close=True)
    print(f"  wrote {out.relative_to(REPO)} (+ .pdf)")


FIG_ALONG = "overwash_vs_footprint_alongshore.png"
FIG_SUMMARY = "overwash_vs_footprint_summary.png"
FIG_MAP = "overwash_vs_footprint_map.png"              # three sections, zoomed
FIG_MAP_ISLAND = "overwash_vs_footprint_map_island.png"  # the whole island, three layers
SECTION_SPLIT = [(1, 30), (31, 60), (61, 90)]

# Diverging classes for the footprint on the map: |rows| 1, 2–3, 4+ (ColorBrewer RdBu)
ADD_CLASS = ["#f4a582", "#d6604d", "#b2182b"]
REM_CLASS = ["#92c5de", "#4393c3", "#2166ac"]
READING_CLR = {
    "agrees (overwash, rows added)": "#b2182b",
    "DISAGREES (overwash, rows removed)": "#2166ac",
    "overwash, no rows changed": "#e08214",
    "unexplained retreat (rows added, no overwash seen)": "#fddbc7",
}
READING_LABEL = {
    "agrees (overwash, rows added)": "overwash, rows added",
    "DISAGREES (overwash, rows removed)": "overwash, rows removed",
    "overwash, no rows changed": "overwash, no rows changed",
    "unexplained retreat (rows added, no overwash seen)": "rows added, no overwash",
}


def _row_class(n):
    a = abs(int(n))
    return 0 if a == 1 else (1 if a <= 3 else 2)


def _frameless_legend(ax, handles, ncol, y=-0.012, fontsize=7.0):
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, y), ncol=ncol,
              fontsize=fontsize, frameon=False, handlelength=1.4, handleheight=1.0,
              columnspacing=1.3, borderaxespad=0, labelspacing=0.5)


def fig_alongshore(t, obs, domains, matrix, idx, out):
    """Three stacked panels on one domain axis, at the double-column width:
    the images, the footprint, the flags."""
    n_img = len(idx)
    x0, x1 = domains[0] - 0.5, domains[-1] + 0.5
    ticks = [d for d in domains if d == 1 or d % 10 == 0]

    fig_w = figsize("double")[0]
    LM, RM = 1.66, 0.40                     # flag names on the left, "rows" right
    W = fig_w - LM - RM
    TM, G = 0.34, 0.38
    H_OW, H_BAR, H_FL = 0.185 * n_img, 2.05, 0.185 * 4
    XL, LG, BM = 0.34, 0.72, 0.06
    fig_h = TM + H_OW + G + H_BAR + G + H_FL + XL + LG + BM
    fig = plt.figure(figsize=figsize("double", height=fig_h))

    def ax_at(y_from_bottom, h):
        return fig.add_axes([LM / fig_w, y_from_bottom / fig_h, W / fig_w, h / fig_h])

    y = fig_h - TM - H_OW
    ax_ow = ax_at(y, H_OW)
    y -= G + H_BAR
    ax_bar = ax_at(y, H_BAR)
    y -= G + H_FL
    ax_fl = ax_at(y, H_FL)

    # (a) the images between the frames
    for k, i in enumerate(idx):
        v = matrix[i]
        for j, d in enumerate(domains):
            if np.isnan(v[j]):
                ax_ow.add_patch(mpatches.Rectangle((d - 0.5, k), 1, 1, fc=CLR_PART_FACE,
                                                   ec=CLR_PART_EDGE, hatch="....", lw=0))
            elif v[j] >= 0.5:
                ax_ow.add_patch(mpatches.Rectangle((d - 0.5, k), 1, 1, fc=CLR_OW, ec="none"))
        ax_ow.axhline(k, color="white", lw=0.5)
    for _, lo, _, _ in SECTIONS[1:]:
        ax_ow.axvline(lo - 0.5, color="#bfbfbf", lw=0.7)
    town_bands(ax_ow)                      # the villages, named here only
    ax_ow.set_xlim(x0, x1)
    ax_ow.set_ylim(n_img, 0)
    ax_ow.set_yticks(np.arange(n_img) + 0.5)
    ax_ow.set_yticklabels([f"{obs.loc[i, 'Imagery_Date']:%Y  %d %b}"
                           + ("*" if obs.loc[i, "poor"] else "") for i in idx], fontsize=7.5)
    ax_ow.tick_params(axis="y", length=0)
    ax_ow.set_xticks(ticks)
    ax_ow.set_xticklabels([])
    ax_ow.tick_params(axis="x", length=2.5, pad=2)
    spines_for_image(ax_ow)
    _title(ax_ow, 0, "overwash seen, Aug 1985 to Oct 1997")

    # (b) the footprint
    for r in t.itertuples():
        c = ACTION_CLR[r.action]
        ax_bar.bar(r.domain, r.shift_m_median, width=1.0, color=c, zorder=2,
                   hatch="////" if r.spread_straddles_zero else None,
                   edgecolor="white" if r.spread_straddles_zero else "none", lw=0)
        ax_bar.plot([r.domain, r.domain], [r.shift_m_p10, r.shift_m_p90],
                    color="#555555", lw=0.5, alpha=0.7, zorder=3)
        if r.overwashed:
            ytxt = (max(r.shift_m_p90, 0) + 4) if r.shift_m_median >= 0 else 4
            ax_bar.plot(r.domain, ytxt, marker="v", ms=3.2, color=CLR_OW, mec="none", zorder=4)
    ax_bar.axhline(0, color=CLR_INK, lw=0.7)
    for _, lo, _, _ in SECTIONS[1:]:
        ax_bar.axvline(lo - 0.5, color="#bfbfbf", lw=0.7, zorder=1)
    town_bands(ax_bar, label=False)
    ax_bar.set_xlim(x0, x1)
    ymin = min(-30, float(t.shift_m_p10.min()) - 5)
    ymax = max(40, float(t.shift_m_p90.max()) + 10) + 6
    ax_bar.set_ylim(ymin, ymax)
    ax_bar.set_xticks(ticks)
    ax_bar.set_xticklabels([str(x) for x in ticks])
    ax_bar.tick_params(axis="both", length=2.5, pad=2)
    ax_bar.set_ylabel("dune-line shift 1984 to 1997 (m)\n+ where the 1984 line lay seaward",
                      fontsize=8)
    ax_r = ax_bar.twinx()
    ax_r.set_ylim(ymin / 10, ymax / 10)
    ax_r.set_yticks(np.arange(np.ceil(ymin / 10), np.floor(ymax / 10) + 1, 2))
    ax_r.tick_params(axis="y", length=2.5)
    ax_r.set_ylabel("rows", fontsize=8)
    ax_r.spines["top"].set_visible(False)
    ax_r.spines["right"].set_visible(True)
    _title(ax_bar, 1, "dune-line shift 1984 to 1997")

    # (c) flags
    # The straddle is not a row here: it is already the hatching in (b).
    flag_rows = [
        ("NC-12 relocated, 1984 to 2004", t.road_relocated_1984_2004, C["ROAD"]),
        # erosion in warm colours, infrastructure in black and grey
        ("erosion hotspot", t.erosion_hotspot.ne(""), "#e08214"),
        ("jetties", t.armor.ne(""), "#8c8c8c"),
        (f"shoreline eroding > {abs(ERODE_THRESH):.0f} m/yr", t.eroding_fast, "#d6604d"),
    ]
    for k, (lab, mask, col) in enumerate(flag_rows):
        for d in t.loc[mask, "domain"]:
            ax_fl.add_patch(mpatches.Rectangle((d - 0.5, k + 0.15), 1, 0.7, fc=col, ec="none"))
    for _, lo, _, _ in SECTIONS[1:]:
        ax_fl.axvline(lo - 0.5, color="#bfbfbf", lw=0.7)
    town_bands(ax_fl, label=False)
    ax_fl.set_xlim(x0, x1)
    ax_fl.set_ylim(len(flag_rows), 0)
    ax_fl.set_yticks(np.arange(len(flag_rows)) + 0.5)
    ax_fl.set_yticklabels([f[0] for f in flag_rows], fontsize=7.5)
    ax_fl.tick_params(axis="y", length=0)
    ax_fl.set_xticks(ticks)
    ax_fl.set_xticklabels([str(x) for x in ticks])
    ax_fl.tick_params(axis="x", length=2.5, pad=2)
    for sp in ax_fl.spines.values():
        sp.set_visible(True)
        sp.set_color("#bbbbbb")
    _title(ax_fl, 2, "management and erosion context")
    ax_fl.set_xlabel(DOMAIN_AXIS_LABEL, labelpad=3)

    # one key under the figure: the panels are full of data and the village
    # names have the top of (a)
    fig.legend(handles=[
        mpatches.Patch(fc=CLR_OW, label="overwash present (a)"),
        mpatches.Patch(fc=CLR_PART_FACE, ec=CLR_PART_EDGE, hatch="....", lw=0,
                       label="image does not reach the domain (a)"),
        mpatches.Patch(fc=CLR_ADD, label="rows added (1984 line seaward)"),
        mpatches.Patch(fc=CLR_REMOVE, label="rows removed (1984 line landward)"),
        mpatches.Patch(fc=CLR_NONE, label="unchanged"),
        mpatches.Patch(fc="#dddddd", ec="white", hatch="////",
                       label="shift spread straddles zero"),
        plt.Line2D([0], [0], marker="v", color="none", markerfacecolor=CLR_OW,
                   markeredgecolor="none", markersize=5, label="overwash seen in (a)")],
        loc="lower center", bbox_to_anchor=(0.5, BM / fig_h), ncol=3, fontsize=7.5,
        frameon=False, handlelength=1.4, handleheight=1.0, columnspacing=1.6,
        borderaxespad=0)
    write(fig, out)


def fig_summary(t, out):
    fig = plt.figure(figsize=figsize("double", aspect=0.46))
    ax_a = fig.add_axes([0.095, 0.21, 0.35, 0.66])
    ax_b = fig.add_axes([0.60, 0.21, 0.37, 0.66])

    # (a) shift by overwash status, with median and IQR
    rng = np.random.default_rng(0)
    ymin = float(t.shift_m_median.min()) - 8
    ymax = float(t.shift_m_median.max()) + 14
    for xi, (lab, mask) in enumerate((("no overwash seen", ~t.overwashed),
                                      ("overwash seen", t.overwashed))):
        g = t[mask]
        q1, med, q3 = np.percentile(g.shift_m_median, [25, 50, 75])
        ax_a.add_patch(mpatches.Rectangle((xi - 0.3, q1), 0.6, q3 - q1, fc="#f0f0f0",
                                          ec="none", zorder=1))
        ax_a.plot([xi - 0.3, xi + 0.3], [med, med], color=CLR_INK, lw=1.3, zorder=4)
        xs = xi + rng.uniform(-0.2, 0.2, len(g))
        ax_a.scatter(xs, g.shift_m_median, s=13, c=[ACTION_CLR[a] for a in g.action],
                     edgecolors="white", linewidths=0.4, zorder=3)
        ax_a.text(xi, ymax - 1, f"n = {len(g)}\nmedian {med:+.0f} m", ha="center", va="top",
                  fontsize=7.5, color=CLR_INK)
    ax_a.axhline(0, color=CLR_INK, lw=0.6, zorder=2)
    ax_a.set_xlim(-0.6, 1.6)
    ax_a.set_ylim(ymin, ymax)
    ax_a.set_xticks([0, 1])
    ax_a.set_xticklabels(["no overwash seen\n1985 to 1997", "overwash seen\n1985 to 1997"],
                         fontsize=7.5)
    ax_a.tick_params(axis="both", length=2.5)
    ax_a.set_ylabel("dune-line shift 1984 to 1997 (m)", fontsize=8)
    open_frame(ax_a)
    _title(ax_a, 0, "shift by overwash status")

    # (b) share of each action that was overwashed, all domains and unflagged
    acts = [("add", "rows added"), ("none", "unchanged"), ("remove", "rows removed")]
    for k, (a, lab) in enumerate(acts):
        for off, (sub, col, alpha) in enumerate(((t, ACTION_CLR[a], 1.0),
                                                 (t[t.unflagged], ACTION_CLR[a], 0.45))):
            g = sub[sub.action == a]
            n, m = int(g.overwashed.sum()), len(g)
            frac = n / m if m else 0
            yb = k + (0.18 if off == 0 else -0.18)
            ax_b.barh(yb, frac, height=0.32, color=col, alpha=alpha, zorder=2)
            ax_b.text(frac + 0.015, yb, f"{n} / {m}", va="center", ha="left",
                      fontsize=7.5, color=CLR_INK)
    ax_b.set_xlim(0, 1.0)
    ax_b.set_ylim(-0.6, len(acts) - 0.4)
    ax_b.set_yticks(range(len(acts)))
    ax_b.set_yticklabels([lab for _, lab in acts], fontsize=7.5)
    ax_b.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax_b.set_xticklabels(["0", "25", "50", "75", "100 %"])
    ax_b.tick_params(axis="both", length=2.5)
    ax_b.invert_yaxis()
    ax_b.set_xlabel("domains with overwash seen", fontsize=8)
    open_frame(ax_b)
    _title(ax_b, 1, "overwash by footprint action")
    handles = [mpatches.Patch(fc="#777777", label="all 90 domains"),
               mpatches.Patch(fc="#777777", alpha=0.45, label="unflagged domains only")]
    ax_b.legend(handles=handles, loc="lower right", fontsize=7.5, frameon=False,
                handlelength=1.3, borderaxespad=0.2)
    write(fig, out)


def fig_map_island(t, geo, out):
    dom, land, road, bounds = geo
    yspan = (bounds[3] - bounds[1]) + PAD_S + PAD_N
    r_a = yspan / ((bounds[2] - bounds[0]) + PAD_W + PAD_E)
    r_bc = yspan / ((bounds[2] - bounds[0]) + 1500.0 + PAD_E)
    fig_w = figsize("double")[0]
    LM, G, RM, TM, LG, BM = 0.16, 0.22, 0.10, 0.36, 1.32, 0.08
    # three islands across the printed width; the height follows from them
    H = (fig_w - LM - 2 * G - RM) / (1 / r_a + 2 / r_bc)
    w_a, w_bc = H / r_a, H / r_bc
    fig_h = TM + H + LG + BM
    fig = plt.figure(figsize=figsize("double", height=fig_h))

    def ax_at(x, w):
        return fig.add_axes([x / fig_w, (BM + LG) / fig_h, w / fig_w, H / fig_h])

    ax_a = ax_at(LM, w_a)
    ax_b = ax_at(LM + w_a + G, w_bc)
    ax_c = ax_at(LM + w_a + G + w_bc + G, w_bc)
    by_dom = t.set_index("domain")

    # (a) overwash seen, number of images
    vmax = max(int(t.n_images_overwash.max()), 1)
    cols = count_cmap(vmax)
    fill = {int(d): cols[int(n) - 1] for d, n in zip(t.domain, t.n_images_overwash) if n > 0}
    draw_island(ax_a, dom, land, road, bounds, fill, "a", "overwash seen")
    scalebar_and_north(ax_a)
    h = [mpatches.Patch(fc=CLR_LAND, ec=CLR_BOX, lw=0.4, label="none seen")]
    h += [mpatches.Patch(fc=cols[i], label=f"{i + 1} image{'s' if i else ''}") for i in range(vmax)]
    h += [plt.Line2D([0], [0], color=CLR_ROAD, lw=0.9, label="NC-12")]
    _frameless_legend(ax_a, h, ncol=2)

    # (b) the footprint, rows added / removed in three classes
    fill, alpha = {}, {}
    for r in t.itertuples():
        if r.n_cells > 0:
            fill[r.domain] = ADD_CLASS[_row_class(r.n_cells)]
        elif r.n_cells < 0:
            fill[r.domain] = REM_CLASS[_row_class(r.n_cells)]
        if r.spread_straddles_zero and r.domain in fill:
            alpha[r.domain] = 0.45
    draw_island(ax_b, dom, land, road, bounds, fill, "b", "the 1984 footprint",
                reach_labels=False, alpha_of=alpha)
    h = [mpatches.Patch(fc=ADD_CLASS[i], label=f"added, {lab}")
         for i, lab in enumerate(("1 row", "2 to 3", "4 or more"))]
    h += [mpatches.Patch(fc=REM_CLASS[i], label=f"removed, {lab}")
          for i, lab in enumerate(("1 row", "2 to 3", "4 or more"))]
    h += [mpatches.Patch(fc=CLR_LAND, ec=CLR_BOX, lw=0.4, label="unchanged"),
          mpatches.Patch(fc=ADD_CLASS[1], alpha=0.45, label="paler: spread straddles 0")]
    _frameless_legend(ax_b, h, ncol=1)

    # (c) the reading
    fill = {int(r.domain): READING_CLR[r.reading] for r in t.itertuples()
            if r.reading in READING_CLR}
    draw_island(ax_c, dom, land, road, bounds, fill, "c", "the reading, per domain",
                reach_labels=False)
    present = [k for k in READING_CLR if (t.reading == k).any()]
    h = [mpatches.Patch(fc=READING_CLR[k], label=f"{READING_LABEL[k]}  ({int((t.reading == k).sum())})")
         for k in present]
    h += [mpatches.Patch(fc=CLR_LAND, ec=CLR_BOX, lw=0.4,
                         label="no overwash,\nno rows added  "
                               f"({int((~t.overwashed & (t.action != 'add')).sum())})")]
    _frameless_legend(ax_c, h, ncol=1)

    write(fig, out)


def _layers(t):
    """The three fills the maps share: (title, fill_of, alpha_of, legend handles)."""
    vmax = max(int(t.n_images_overwash.max()), 1)
    cols = count_cmap(vmax)
    fill_a = {int(d): cols[int(n) - 1] for d, n in zip(t.domain, t.n_images_overwash) if n > 0}
    h_a = [mpatches.Patch(fc=CLR_LAND, ec=CLR_BOX, lw=0.4, label="none seen")]
    h_a += [mpatches.Patch(fc=cols[i], label=f"{i + 1} image{'s' if i else ''}") for i in range(vmax)]
    h_a += [plt.Line2D([0], [0], color=CLR_ROAD, lw=0.9, label="NC-12")]

    fill_b, alpha_b = {}, {}
    for r in t.itertuples():
        if r.n_cells > 0:
            fill_b[r.domain] = ADD_CLASS[_row_class(r.n_cells)]
        elif r.n_cells < 0:
            fill_b[r.domain] = REM_CLASS[_row_class(r.n_cells)]
        if r.spread_straddles_zero and r.domain in fill_b:
            alpha_b[r.domain] = 0.45
    h_b = [mpatches.Patch(fc=ADD_CLASS[i], label=f"added, {lab}")
           for i, lab in enumerate(("1 row", "2 to 3", "4 or more"))]
    h_b += [mpatches.Patch(fc=REM_CLASS[i], label=f"removed, {lab}")
            for i, lab in enumerate(("1 row", "2 to 3", "4 or more"))]
    h_b += [mpatches.Patch(fc=CLR_LAND, ec=CLR_BOX, lw=0.4, label="unchanged"),
            mpatches.Patch(fc=ADD_CLASS[1], alpha=0.45, label="paler: spread straddles 0")]

    fill_c = {int(r.domain): READING_CLR[r.reading] for r in t.itertuples()
              if r.reading in READING_CLR}
    present = [k for k in READING_CLR if (t.reading == k).any()]
    h_c = [mpatches.Patch(fc=READING_CLR[k],
                          label=f"{READING_LABEL[k]}  ({int((t.reading == k).sum())})")
           for k in present]
    h_c += [mpatches.Patch(fc=CLR_LAND, ec=CLR_BOX, lw=0.4,
                           label=f"no overwash, no rows added  "
                                 f"({int((~t.overwashed & (t.action != 'add')).sum())})")]
    return [("overwash seen", fill_a, {}, h_a),
            ("dune-line shift 1984 to 1997", fill_b, alpha_b, h_b),
            ("reading", fill_c, {}, h_c)]


def fig_map_sections(t, geo, out):
    """Three alongshore sections, each zoomed, each with the three layers."""
    dom, land, road, bounds = geo
    layers = _layers(t)[:2]          # overwash seen, dune-line shift; the reading is in the tables
    PAD_NS, pad_w, pad_e = 500.0, 900.0, 1000.0
    groups = []
    for lo, hi in SECTION_SPLIT:
        sub = dom[dom["domain_id"].between(lo, hi)].reset_index(drop=True)
        b = sub.total_bounds
        yspan = (b[3] - b[1]) + 2 * PAD_NS
        xspan = (b[2] - b[0]) + pad_w + pad_e
        groups.append(dict(lo=lo, hi=hi, dom=sub, b=b, ratio=yspan / xspan,
                           title=f"domains {lo}–{hi}"))

    LM, G_IN, G_GROUP, RM = 0.16, 0.10, 0.34, 0.10
    TM, LG, BM = 0.52, 0.95, 0.08
    n_lay = len(layers)
    fig_w = figsize("double")[0]
    # every panel the same height; the width each one needs follows its span
    H = ((fig_w - LM - RM - 2 * G_GROUP - len(groups) * (n_lay - 1) * G_IN)
         / sum(n_lay / g["ratio"] for g in groups))
    for g in groups:
        g["w"] = H / g["ratio"]
    fig_h = TM + H + LG + BM
    fig = plt.figure(figsize=figsize("double", height=fig_h))
    letters = iter("abcdefghi")
    x = LM
    legend_anchor_x = []
    for k, g in enumerate(groups):
        x_group0 = x
        for li, (lname, fill, alpha, handles) in enumerate(layers):
            ax = fig.add_axes([x / fig_w, (BM + LG) / fig_h, g["w"] / fig_w, H / fig_h])
            first = li == 0
            # the panel is one column wide: the letter alone goes above it and
            # the legend beneath names the layer
            draw_island(ax, g["dom"], land, road, g["b"], fill, next(letters), "",
                        reach_labels=first, pad_w=pad_w, alpha_of=alpha, pad_e=pad_e,
                        pad_s=PAD_NS, pad_n=PAD_NS, label_every=5,
                        reach_rotation=90)
            if k == 0 and first:
                scalebar_and_north(ax, length_m=2000.0)
            if k == 0:
                legend_anchor_x.append(x / fig_w)
            x += g["w"] + G_IN
        x += G_GROUP - G_IN
        fig.text((x_group0 + 0.5 * (n_lay * g["w"] + (n_lay - 1) * G_IN)) / fig_w,
                 (BM + LG + H + 0.30) / fig_h,
                 g["title"], ha="center", va="bottom", fontsize=9, color=CLR_INK)

    # one legend per layer, side by side under the figure: the panels are a
    # single column wide, too narrow to carry a legend under each of them
    for (lname, _, _, handles), lx in zip(layers, (LM / fig_w, 0.42)):
        fig.legend(handles=handles, loc="upper left",
                   bbox_to_anchor=(lx, (BM + LG - 0.06) / fig_h),
                   ncol=2 if len(handles) > 6 else 1,
                   fontsize=7, frameon=False, handlelength=1.4, handleheight=1.0,
                   borderaxespad=0, labelspacing=0.4, title=lname, title_fontsize=7.5,
                   alignment="left")
    write(fig, out)


# ================================================================== captions
def captions(t, idx, obs):
    ow = t[t.overwashed]
    n_add, n_rem = int((t.action == "add").sum()), int((t.action == "remove").sum())
    f = fisher(t)
    window = (f"the {len(idx)} images taken after the frame the 1984 dune line was "
              f"digitised on ({LINE_1984:%d %b %Y}) and up to the frame the 1997 line "
              f"was digitised on ({LINE_1997:%d %b %Y})")
    readings = (f"Given overwash ({len(ow)} domains): {int((ow.action == 'add').sum())} rows "
                f"added, {int((ow.action == 'none').sum())} unchanged, "
                f"{int((ow.action == 'remove').sum())} rows removed. Given the action: "
                f"{int(t[t.action == 'add'].overwashed.sum())} of {n_add} rows-added and "
                f"{int(t[t.action == 'remove'].overwashed.sum())} of {n_rem} rows-removed "
                f"domains show overwash.")
    absence = ("A rows-added domain without overwash is not counted against the "
               "footprint: the images are two to four years apart and washover fades "
               "from imagery, so absence is weak evidence.")
    along = (
        f"Observed overwash 1985 to 1997 against the 1984 footprint, by domain (1 at "
        f"Cape Point, 90 at Pea Island; the light bands are the villages, named across "
        f"the top of (a)). (a) {window}; purple where overwash was present, "
        f"dotted where the image does not reach the domain, * poor image. (b) The "
        f"footprint's median shift between the two dune lines with the p10 to p90 over "
        f"50 profiles, positive where the 1984 line lay seaward (rows added, red), "
        f"negative where it lay landward (rows removed, blue), grey where under one "
        f"10 m cell; hatched where the spread straddles zero; right axis in 10 m rows; ▼ "
        f"marks the domains with overwash in (a). (c) Flags carried, not dropped: "
        f"NC-12 relocated over 1984 to 2004, the Buxton, Avon and Rodanthe S-curves "
        f"erosion hotspots, the domain-6 jetties, and CoastSat 1984 to 2004 shoreline "
        f"erosion faster than {abs(ERODE_THRESH):.0f} m/yr; the straddle flag is the "
        f"hatching in (b). {readings} {absence}")
    summ = (
        f"The comparison in two numbers. (a) The dune-line shift of every domain by "
        f"whether overwash was seen in {window}, coloured by footprint action as in the "
        f"alongshore figure; bar at the median, box over the interquartile range. "
        f"(b) The share of each action's domains with overwash seen, all 90 domains "
        f"(solid) and the {int(t.unflagged.sum())} domains carrying no flag (pale), "
        f"counts beside the bars. {readings}"
        + (f" Fisher exact on overwashed against rows added versus removed, unchanged "
           f"left out: p = {f['p']:.3f}; the table is the evidence, not the p-value." if f else "")
        + f" {absence}")
    mp = (
        f"The whole island at one scale, for orientation (NC 1:80k coastline, UTM 18N), the "
        f"90 domain boxes clipped to land, NC-12 as the dark line, domain numbers every "
        f"ten on the ocean side, reaches on the sound side of (a). (a) Number of images "
        f"with overwash in {window}, in the purple of the overwash record. (b) The "
        f"1984 footprint: rows added where the 1984 dune line lay seaward of the 1997 "
        f"line, rows removed where it lay landward, in three classes; paler where the "
        f"shift's p10 to p90 straddles zero. (c) The reading per domain: overwash seen "
        f"and rows added is the pairing the physics predicts; overwash seen with rows "
        f"removed would be the contradiction (none this run); rows added without overwash "
        f"seen is listed, not counted against the footprint. {readings} {absence}")
    secs = (
        f"The comparison on the island outline in three alongshore sections, "
        f"domains 1 to 30, 31 to 60 and 61 to 90 (NC 1:80k coastline, UTM 18N; "
        f"the domain boxes clipped to land, NC-12 as the dark line, domain numbers "
        f"every five on the ocean side, reaches along the sound side of the first "
        f"panel of each section). Within each section, left: the number of images "
        f"with overwash in {window}; right: the dune-line shift 1984 to 1997 as the "
        f"rows the 1984 reconstruction adds where the 1984 dune line lay seaward of "
        f"the 1997 line and removes where it lay landward, in three classes, paler "
        f"where the shift's p10 to p90 straddles zero. {readings} {absence}")
    return [(FIG_ALONG, along), (FIG_SUMMARY, summ), (FIG_MAP, secs), (FIG_MAP_ISLAND, mp)]


# ===================================================================== main
def main():
    apply_style()
    VS_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    obs, domains, matrix, idx, fp, rel, cs, ds, hotspot, armor = load_inputs()
    t = build_table(obs, domains, matrix, idx, fp, rel, cs, ds, hotspot, armor)

    t.to_csv(VS_DIR / "overwash_vs_footprint_by_domain.csv", index=False)
    pd.concat([contingency(t, "all domains"),
               contingency(t[t.unflagged], "unflagged only")]).to_csv(
        VS_DIR / "overwash_vs_footprint_contingency.csv")
    text = summary_text(t, idx, obs)
    (VS_DIR / "overwash_vs_footprint_summary.txt").write_text(text, encoding="utf-8")
    print(text)

    fig_alongshore(t, obs, domains, matrix, idx, FIG_DIR / FIG_ALONG)
    fig_summary(t, FIG_DIR / FIG_SUMMARY)
    geo = load_geometry()
    fig_map_sections(t, geo, FIG_DIR / FIG_MAP)
    fig_map_island(t, geo, FIG_DIR / FIG_MAP_ISLAND)
    remove_caption("overwash_vs_footprint.png")
    for name, text in captions(t, idx, obs):
        p = upsert_caption(name, "vs-footprint", text)
    print(f"  wrote {p.relative_to(REPO)}")


if __name__ == "__main__":
    main()
