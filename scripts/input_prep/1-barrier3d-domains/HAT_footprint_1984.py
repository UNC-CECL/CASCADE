#!/usr/bin/env python3
r"""
HAT_footprint_1984.py
==============================================================================
Where the Barrier3D interior has to grow or shrink for its seaward row 0 to
stand where the 1984 dune line stood - per domain, in whole 10 m cells, in
BOTH directions - and what the 1984 road setback becomes when it does.

SCOPE ONLY. No array is written, no elevation is fabricated, no model-facing
CSV is touched. The rows are drawn BLANK: this says where they land and how
many, so the fill can be argued separately against a known footprint.

DECIDED WITH HANNAH, 2026-09-07 (each one changes what the model would ingest)
    symmetric     rows are ADDED where the 1984 line lies seaward of the 1997
                  line and REMOVED where it lies landward. The earlier work
                  (layers v3-v8, deleted 2026-09-07) floored negatives to 0.
    paired median the per-domain shift is the MEDIAN OF THE 50 PAIRED
                  PER-PROFILE DIFFERENCES, line_1997 - line_1984 on the same
                  raster row. duneline_retreat_1984_1997.csv stores the
                  difference of two medians instead; the two disagree by one
                  cell at 13 domains (GIS 80: 7 vs 6). The paired form also
                  gives a p10-p90 spread, which the stored file cannot.
    10 m rule     N = trunc(shift / 10 m): a row only when a FULL cell of
                  change is measured. Understates by 0-9 m at every changed
                  domain, always toward less change; the residual is a column.
    row-0 setback the new setback keeps the model's reference (metres landward
                  of interior row 0, which is one cell behind the picked crest):
                      setback_new(p) = (road(p) - row0(p)) + (line97(p) - line84(p))
                  per profile, UNROUNDED, median per domain. NOT road minus the
                  1984 dune line: the digitized lines trace the toe, ~19 m
                  seaward of row 0 (IQR 15-28 m over the road domains), so that
                  reading is ~20 m larger everywhere and would silently change
                  the convention every run so far has used. It is kept as a
                  record column, `setback_raw84_m`, and nothing reads it.

WHAT IS ASSUMED (and cannot be checked from these files)
    * the 1984 and 1997 lines trace the SAME feature. 1997 carries metadata
      saying "light/dark elevation break"; 1984 carries none.
    * one integer per domain. The interior is rectangular, so the alongshore
      median stands for 50 profiles and the spread inside a domain is lost.
    * 1997 stands for 1996. The surface is 1996 ALACE; the line is a year later.
    * the 1984 dune crest equalled the 1996 one - the dune array is not
      re-estimated.
    * removal deletes SURVEYED 1996 beach/foredune cells that were ocean in 1984.
    * only the island width moves. The ocean shoreline is the shoreline-offset
      input; rows at the dune move the bay edge.
    * the 1984 road line is the 1978 export (deliberate, recorded elsewhere).

INPUTS (all already on disk; nothing is re-measured against GIS here)
    row-insert-scope/duneline-shift/duneline_shift_{1984,1997}_profiles.csv
        per (domain, profile): the line's crossing as a cross-shore cell and
        interior row 0, both in the extractor's own c0/shear frame
    4-mgmt-forcing/road_offset/dunestart_offset/1984/RoadOffset_1984_profiles.csv
        per (domain, profile): the road's seaward cell and interior row 0, same
        frame. Row 0 is asserted identical across the three files.
    4-mgmt-forcing/road_offset/dunestart_offset/1984/RoadOffset_1984_domains.csv
        the setback the model currently receives (setback_model_m) and flags
    1984-start/dune-topo/<CURRENT>/topography   rows_now, and the grid figure

OUTPUTS  row-insert-scope/
    footprint_1984_by_domain.csv     one row per domain (the audit table)
    footprint_1984_profiles.csv      the per-profile join the medians come from
    HAT_footprint_1984.txt           the report
    figures/1-scope/HAT_footprint_1984_grid.png   the Barrier3D grid, both signs
    figures/1-scope/HAT_footprint_1984_plan.png   plan view on the DEM, both lines, NC-12
    figures/1-scope/HAT_footprint_1984_rows.png     rows per domain, signed
    figures/1-scope/HAT_footprint_1984_shift.png    the paired shift with spread and the rows kept
    figures/1-scope/HAT_footprint_1984_setback.png  the road setback now and from the new row 0
    figures/CAPTIONS.md              the words that go under the figures

USAGE
    python HAT_footprint_1984.py            # everything
    python HAT_footprint_1984.py --no-plan  # skip the slow DEM panel
==============================================================================
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
INIT = REPO / "data" / "hatteras_init"
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "0-elevation" / "3-figures"))
from hat_topo_version import (  # noqa: E402
    array_name, duneline_shift_dir, insert_figures_dir, topo_dirs)
import HAT_plot_duneline_offset as off  # noqa: E402  the house style + map loaders
from hat_figure_style import elevation_cmap  # noqa: E402  the elevation classes

PRODUCT = "1984-start"
CELL_M = 10.0
SCOPE_DIR = INIT / "1-barrier3d-domains" / PRODUCT / "row-insert-scope"
SHIFT_DIR = duneline_shift_dir(PRODUCT)
ROAD_DIR = INIT / "4-mgmt-forcing" / "road_offset" / "dunestart_offset" / "1984"
FIG_DIR = insert_figures_dir(PRODUCT, "1-scope")
CAPTIONS = SCOPE_DIR / "figures" / "CAPTIONS.md"

# Colours: the RdBu pair the dune-line figures use. Red is 1984 / seaward /
# ground ADDED; blue is 1997 / landward / ground REMOVED. Same meaning on every
# panel of every figure here.
C_ADD, C_ADD_FILL = off.C_1984, off.C_1984_FILL
C_REM, C_REM_FILL = off.C_1997, off.C_1997_FILL
C_LAND, C_WATER = "#f0e6c8", "#a8c8e0"
C_ROAD = "#1a1a1a"
INK = off.INK

ROWS_SHOWN = 190          # every interior row (the deepest domain has 189); was 40 until 2026-09-07, when Hannah asked for the full domains
DOMAINS_PER_STRIP = 30
NEAR_ZERO_M = 10.0        # a new setback under one cell is flagged


# =============================================================================
# THE NUMBERS
# =============================================================================

def load_profiles() -> pd.DataFrame:
    """The per-profile join: both dune-line crossings and the road, one frame."""
    p84 = pd.read_csv(SHIFT_DIR / "duneline_shift_1984_profiles.csv")
    p97 = pd.read_csv(SHIFT_DIR / "duneline_shift_1997_profiles.csv")
    road = pd.read_csv(ROAD_DIR / "RoadOffset_1984_profiles.csv")

    m = p84.merge(p97, on=["domain", "profile"], suffixes=("_84", "_97"))
    if not (m["interior_row0_cell_84"] == m["interior_row0_cell_97"]).all():
        raise SystemExit("row 0 differs between the 1984 and 1997 profile files "
                         "- they were measured on different extractions")
    m = m.rename(columns={"interior_row0_cell_84": "row0_cell",
                          "duneline_cell_84": "line84_cell",
                          "duneline_cell_97": "line97_cell"})
    m = m[["domain", "profile", "row0_cell", "line84_cell", "line97_cell"]]
    # + = the 1984 line lies SEAWARD of the 1997 line = the island retreated
    #     = rows to ADD. Cells grow landward, so seaward is the smaller index.
    m["shift_m"] = (m["line97_cell"] - m["line84_cell"]) * CELL_M

    r = road[["domain", "profile", "interior_row0_cell", "road_seaward_cell"]]
    m = m.merge(r, on=["domain", "profile"], how="left")
    has = m["road_seaward_cell"].notna()
    if not (m.loc[has, "interior_row0_cell"] == m.loc[has, "row0_cell"]).all():
        raise SystemExit("row 0 differs between the road and dune-line profile "
                         "files - re-run HAT_road_offset_from_dune_start.py on "
                         "the current extraction")
    m = m.drop(columns=["interior_row0_cell"])
    m["setback_v2_m"] = (m["road_seaward_cell"] - m["row0_cell"]) * CELL_M
    m["setback_new_m"] = m["setback_v2_m"] + m["shift_m"]          # row-0 convention
    m["setback_raw84_m"] = (m["road_seaward_cell"] - m["line84_cell"]) * CELL_M
    return m


def n_cells(shift_m: float) -> int:
    """trunc(shift / 10): a row only once a FULL cell of change is measured."""
    return int(np.trunc(shift_m / CELL_M))


def by_domain(prof: pd.DataFrame, topo_dir: Path, topo_name: str) -> pd.DataFrame:
    dom_csv = pd.read_csv(ROAD_DIR / "RoadOffset_1984_domains.csv").set_index("domain")
    rows = []
    for d, g in prof.groupby("domain"):
        s = g["shift_m"].to_numpy()
        med = float(np.median(s))
        n = n_cells(med)
        p10, p90 = float(np.percentile(s, 10)), float(np.percentile(s, 90))
        f = topo_dir / array_name("topography", d)
        rows_now = int(np.load(f).shape[0]) if f.is_file() else -1
        rec = {
            "domain": int(d), "topo_version": topo_name, "n_profiles": len(s),
            "shift_m_median": round(med, 1), "shift_m_p10": round(p10, 1),
            "shift_m_p90": round(p90, 1),
            "spread_straddles_zero": int(p10 < 0 < p90),
            "n_cells": n, "action": "add" if n > 0 else ("remove" if n < 0 else "none"),
            "residual_m": round(med - n * CELL_M, 1),
            "rows_now": rows_now, "rows_after": rows_now + n if rows_now > 0 else -1,
        }
        rd = g.dropna(subset=["road_seaward_cell"])
        rec["n_road_profiles"] = len(rd)
        flags = []
        if len(rd):
            sb_new = rd["setback_new_m"].to_numpy()
            rec["setback_v2_m"] = round(float(np.median(rd["setback_v2_m"])), 1)
            rec["setback_new_m"] = round(float(np.median(sb_new)), 1)
            rec["setback_new_p10_m"] = round(float(np.percentile(sb_new, 10)), 1)
            rec["setback_new_p90_m"] = round(float(np.percentile(sb_new, 90)), 1)
            rec["setback_derived_m"] = round(rec["setback_v2_m"] + n * CELL_M, 1)
            rec["setback_raw84_m"] = round(float(np.median(rd["setback_raw84_m"])), 1)
            if rec["setback_new_m"] < 0:
                flags.append("NEGATIVE_NEW")
            elif rec["setback_new_m"] < NEAR_ZERO_M:
                flags.append("NEAR_ZERO_NEW(<%.0fm)" % NEAR_ZERO_M)
        else:
            for k in ("setback_v2_m", "setback_new_m", "setback_new_p10_m",
                      "setback_new_p90_m", "setback_derived_m", "setback_raw84_m"):
                rec[k] = np.nan
        # what the model receives TODAY (floored, drowning-relocated), for the
        # before/after panel; NaN where the road is outside the managed span
        rec["setback_model_now_m"] = (float(dom_csv.loc[d, "setback_model_m"])
                                      if d in dom_csv.index else np.nan)
        if d in dom_csv.index and isinstance(dom_csv.loc[d, "flags"], str) \
                and "EXCLUDED_FROM_SPAN" in dom_csv.loc[d, "flags"]:
            flags.append("ROAD_EXCLUDED_FROM_SPAN")
        if rec["spread_straddles_zero"] and n != 0:
            flags.append("SPREAD_STRADDLES_ZERO")
        rec["flags"] = ",".join(flags)
        rows.append(rec)
    return pd.DataFrame(rows).set_index("domain").sort_index()


# =============================================================================
# FIGURE 1 - THE BARRIER3D GRID, BOTH SIGNS
# =============================================================================

def _elev_rgba(topo_dam: np.ndarray) -> np.ndarray:
    """Elevation classes (hat_figure_style) for one interior, as RGBA."""
    cmap, norm, _ = elevation_cmap()
    return cmap(norm(topo_dam * CELL_M))


def _community_bar(ax, y_bar: float, x_lo: float, x_hi: float) -> None:
    """The communities as a bar above the strip, names above it, from
    HATTERAS_ANNOTATIONS - the same object every other island figure uses."""
    ann = off.HATTERAS_ANNOTATIONS
    for name, (lo, hi) in ann.town_spans.items():
        a, b = max(lo - 0.5, x_lo), min(hi + 0.5, x_hi)
        if b <= a:
            continue
        ax.plot([a, b], [y_bar, y_bar], color=ann.color_town_span, lw=4.0,
                solid_capstyle="butt", zorder=8, clip_on=False)
        if (b - a) >= 0.5 * (hi - lo + 1):
            ax.text((a + b) / 2, y_bar - 0.9, name, ha="center", va="bottom",
                    fontsize=8, fontweight="bold", color=INK, clip_on=False)
    for name, gid in ann.village_lines.items():
        if x_lo <= gid <= x_hi:
            ax.plot([gid, gid], [y_bar - 0.7, y_bar + 0.7], color=ann.color_village_line,
                    lw=1.0, zorder=9, clip_on=False)
            ax.text(gid, y_bar - 0.9, name, ha="center", va="bottom", fontsize=7,
                    color=ann.color_village_line, clip_on=False)


def fig_grid(tab: pd.DataFrame, topo_dir: Path) -> Path:
    """
    Every interior, every row, in ONE frame: current interior row 0 is y = 0 on
    every domain. Existing cells are drawn in the project's elevation classes,
    so the dune ridge, the backbarrier flat and the sound-side marsh read as
    what they are; added rows sit ABOVE 0 (between the new row 0 and the old),
    blank because no fill has been chosen; removed rows are the existing rows
    0..|N|-1, hatched. The black tick is where interior row 0 ends up - the
    1984 dune line, one cell behind the crest - and the dark bar is NC-12 at
    its measured position, which does not move.
    """
    off.apply_style()
    doms = tab.index.to_numpy()
    n_by = tab["n_cells"].to_dict()
    top = int(max(1, tab["n_cells"].max()))
    ymin = -top - 9.0                                    # room for labels + communities
    groups = [doms[i:i + DOMAINS_PER_STRIP] for i in range(0, len(doms), DOMAINS_PER_STRIP)]
    add_rgba = np.array(matplotlib.colors.to_rgba(C_ADD_FILL))

    fig, axes = plt.subplots(len(groups), 1, figsize=(15.0, 5.2 * len(groups) + 1.4),
                             constrained_layout=True)
    axes = np.atleast_1d(axes)
    for k, (ax, g) in enumerate(zip(axes, groups)):
        cols = []
        for d in g:
            topo = np.load(topo_dir / array_name("topography", d))
            n = n_by[d]
            H = top + ROWS_SHOWN
            img = np.ones((H, topo.shape[1], 4))                 # white = off the array
            take = min(topo.shape[0], ROWS_SHOWN)
            img[top:top + take] = _elev_rgba(topo[:take])
            if n > 0:
                img[top - n:top] = add_rgba
            cols.append(img)
        ax.imshow(np.concatenate(cols, axis=1), aspect="auto", interpolation="nearest",
                  origin="upper", extent=[g[0] - 0.5, g[-1] + 0.5, ROWS_SHOWN, -top])
        for d in g:
            n = n_by[d]
            ax.axvline(d + 0.5, color="0.55", linewidth=0.35, zorder=3)
            if n < 0:
                ax.add_patch(Rectangle((d - 0.5, 0.0), 1.0, -n, facecolor="none",
                                       edgecolor=C_REM, hatch="//////", linewidth=0.0, zorder=4))
                ax.add_patch(Rectangle((d - 0.5, 0.0), 1.0, -n, facecolor="none",
                                       edgecolor=C_REM, linewidth=0.6, zorder=4))
            ax.plot([d - 0.5, d + 0.5], [-n, -n], color=INK if n else "0.5",
                    linewidth=1.7 if n else 0.6, solid_capstyle="butt", zorder=6)
            if n:
                ax.text(d, -top - 0.6, f"{n:+d}", ha="center", va="bottom", fontsize=7.2,
                        color=C_ADD if n > 0 else C_REM, fontweight="bold")
            sb = tab.loc[d, "setback_v2_m"]
            if np.isfinite(sb):
                ax.add_patch(Rectangle((d - 0.32, sb / CELL_M), 0.64, 2.0,
                                       facecolor=C_ROAD, edgecolor="none", zorder=7))
        _community_bar(ax, -top - 5.0, g[0] - 0.5, g[-1] + 0.5)
        ax.set_xlim(g[0] - 0.5, g[-1] + 0.5)
        ax.set_ylim(ROWS_SHOWN, ymin)
        ax.set_xticks([d for d in g if d % 5 == 0])
        ax.set_xticks(list(g), minor=True)
        ax.set_yticks(range(0, ROWS_SHOWN, 25))
        ax.set_ylabel("cross-shore cell\n(0 = current interior row 0)")
        ax.grid(axis="y", color="0.9", linewidth=0.4)
        ax.set_axisbelow(True)
        sec = ax.secondary_yaxis("right", functions=(lambda c: c * CELL_M, lambda m: m / CELL_M))
        sec.set_ylabel("m landward of current row 0")
        sec.set_yticks(range(0, int(ROWS_SHOWN * CELL_M), 500))
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(True)
        off._title(ax, k, f"Domains {g[0]}\u2013{g[-1]}")
    axes[-1].set_xlabel("domain (1 = south, Cape Hatteras)")

    cmap, _norm, bounds = elevation_cmap()
    labels = ["below 0 (water)"] + [f"{lo:g}\u2013{hi:g}" for lo, hi in zip(bounds[1:-2], bounds[2:-1])] \
        + [f"above {bounds[-2]:g}"]
    handles = [Patch(facecolor=cmap(i), edgecolor="0.4", linewidth=0.4, label=lab)
               for i, lab in enumerate(labels)]
    handles += [Patch(facecolor="white", edgecolor="0.4", linewidth=0.4, label="off the array"),
                Patch(facecolor=C_ADD_FILL, edgecolor="none", label="rows added (blank, no fill yet)"),
                Patch(facecolor="none", edgecolor=C_REM, hatch="//////", label="existing rows removed"),
                Line2D([0], [0], color=INK, linewidth=1.7, label="new interior row 0 (1984 dune line)"),
                Patch(facecolor=C_ROAD, edgecolor="none", label="NC-12 (1984), measured position"),
                Line2D([0], [0], color=off.HATTERAS_ANNOTATIONS.color_town_span, lw=4.0, label="community")]
    fig.legend(handles=handles, loc="outside lower center", ncol=7, fontsize=8,
               title="existing interior, elevation classes (m MHW)", title_fontsize=8)
    p = FIG_DIR / "HAT_footprint_1984_grid.png"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


def fig_rows(tab: pd.DataFrame) -> Path:
    """Rows per domain, signed, on its own - the communities banded along the
    axis so a domain can be placed without the map."""
    off.apply_style()
    doms = tab.index.to_numpy()
    nn = tab["n_cells"].to_numpy()
    fig, ax = plt.subplots(figsize=(13.0, 3.8), constrained_layout=True)
    ann = off.HATTERAS_ANNOTATIONS
    for name, (lo, hi) in ann.town_spans.items():
        ax.axvspan(lo - 0.5, hi + 0.5, color="0.93", zorder=0)
        ax.text((lo + hi) / 2, 0.985, name, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=7.5, color=off.INK_MUTED)
    ax.bar(doms, nn, width=0.82, color=np.where(nn > 0, C_ADD, C_REM), edgecolor="none", zorder=3)
    ax.axhline(0, color=INK, linewidth=0.7, zorder=2)
    ax.set_xlabel("domain (1 = south, Cape Hatteras)")
    ax.set_ylabel("rows added (+) / removed (\u2212)")
    ax.set_xlim(doms[0] - 0.8, doms[-1] + 0.8)
    ax.set_xticks([1] + list(range(10, int(doms.max()) + 1, 10)))
    ax.set_xticks(list(doms), minor=True)
    ax.set_yticks(range(int(nn.min()), int(nn.max()) + 1))
    ax.grid(axis="y", color="0.92", linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    add, rem = nn[nn > 0], nn[nn < 0]
    ax.text(0.2, 0.03,
            f"{len(add)} domains gain {int(add.sum())} rows;  {len(rem)} domains lose "
            f"{int(-rem.sum())};  {int((nn == 0).sum())} unchanged.  "
            f"N = trunc(median paired shift / 10 m)",
            transform=ax.transAxes, fontsize=8, ha="left", va="bottom", color=off.INK_MUTED)
    ax.legend(handles=[Patch(facecolor=C_ADD, label="rows added"),
                       Patch(facecolor=C_REM, label="existing rows removed")],
              loc="upper center", ncol=2)
    p = FIG_DIR / "HAT_footprint_1984_rows.png"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


# =============================================================================
# FIGURE 2 - PLAN VIEW ON THE DEM
# =============================================================================

def fig_plan(tab: pd.DataFrame) -> Path:
    """
    The footprint in plan view, in the layout of the dune-line figure
    HAT_duneline_offset_lines_island_3panel.png (2026-09-07, Hannah's request):
    three panels of thirty domains, each cropped to the strip the two dune lines
    and NC-12 occupy, a grey hillshade backdrop with no readable elevation, no
    coordinate ticks (scale bar and north arrow instead), the communities as a
    bracket in the ocean margin, domain numbers on the landward edge.

    Two encodings, because one cannot work alone at island scale: each domain
    box is shaded by N (red added, blue removed), and the TRUE-SCALE band is
    drawn on top - the 1997 dune line offset by N x 10 m seaward (add) or
    landward (remove). The 1997 line is the map-space proxy for the existing
    array's seaward edge; the digitized line itself sits ~19 m seaward of row 0.
    """
    from shapely.geometry import box as _box
    off.apply_style()
    print("  plan view: loading the island mosaic and re-sampling the lines ...")
    elev, _surv, extent, _n = off.m.load_mosaic()
    gdf = off.load_domains()
    lines = off.load_lines(gdf.crs)
    geoms = {yr: g.union_all() for yr, g in lines.items()}
    _rows, samples = off.measure(gdf, geoms)
    drawn = off.clip_for_drawing(lines, gdf.union_all())
    roads_all = off.m.load_roads(gdf.crs, clip_to=gdf.union_all())
    roads = {1984: roads_all[1984]} if 1984 in roads_all else {}
    by_dom = {r["domain"]: r for r in off.read_rows(off.OUT_DIR / off.CSV_NAME)}

    n_by = tab["n_cells"].to_dict()
    nmax = int(max(abs(v) for v in n_by.values()))
    reds = plt.get_cmap("Reds")(np.linspace(0.30, 0.90, nmax))
    blues = plt.get_cmap("Blues")(np.linspace(0.30, 0.90, nmax))

    ids = np.sort(gdf["domain_id"].astype(int).to_numpy())
    groups = np.array_split(ids, off.N_PANELS)
    pad_m = off.LINES_ISLAND_PAD_M

    # Each panel's window is the full extent of its thirty domain BOXES plus
    # a pad - not the strip the lines occupy, as the dune-line figure crops.
    # The boxes are the subject here (they carry N), so every one is shown
    # whole, sound side included (Hannah, 2026-09-07).
    wins = []
    for group in groups:
        sub = gdf[gdf["domain_id"].astype(int).isin(group)]
        bx = sub.total_bounds
        wins.append((bx[0] - pad_m, bx[2] + pad_m,
                     bx[1] - 0.2 * pad_m, bx[3] + 0.2 * pad_m))
    ratios = [(w[1] - w[0]) / (w[3] - w[2]) for w in wins]

    panel_h = 13.0
    fig, axes = plt.subplots(1, len(groups),
                             figsize=(panel_h * sum(ratios) + 0.35 * len(groups) + 0.8,
                                      panel_h + 1.9),
                             constrained_layout=True,
                             gridspec_kw=dict(width_ratios=ratios))
    axes = np.atleast_1d(axes)

    for i, (ax, group, win) in enumerate(zip(axes, groups, wins)):
        x0, x1, y0, y1 = win
        sub = gdf[gdf["domain_id"].astype(int).isin(group)]
        off._hillshade(ax, elev, extent, res=off.GRID_10M)
        for _, r in sub.iterrows():
            d = int(r["domain_id"])
            n = n_by.get(d, 0)
            b = r.geometry.bounds
            if n:
                ax.add_patch(Rectangle((b[0], b[1]), b[2] - b[0], b[3] - b[1],
                                       facecolor=(reds if n > 0 else blues)[abs(n) - 1],
                                       edgecolor="none", alpha=0.45, zorder=2))
                m_ = samples["domain"] == d
                y = samples["y"][m_]
                x97 = samples["x1997"][m_]
                ok = np.isfinite(x97)
                if ok.sum() > 2:
                    # cross-shore grows landward = west, so seaward is +x
                    x_edge = x97[ok] + n * CELL_M
                    poly = np.concatenate([np.column_stack([x97[ok], y[ok]]),
                                           np.column_stack([x_edge[::-1], y[ok][::-1]])])
                    ax.add_patch(plt.Polygon(poly, closed=True,
                                             facecolor=C_ADD if n > 0 else C_REM,
                                             edgecolor="none", alpha=0.95, zorder=5))
        sub.boundary.plot(ax=ax, color="0.45", linewidth=0.4, zorder=4)
        off.draw_lines(ax, drawn, scale=off.LINES_ISLAND_LINE_SCALE,
                       style=off.SIMPLE_LINE_STYLE)
        if roads:
            off.m.draw_roads(ax, roads, scale=0.55)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        off._title(ax, i, f"Domains {group.min()}\u2013{group.max()}")
        off._places(ax, gdf, y0, y1)
        off._structures(ax, gdf, by_dom, y0, y1)
        if i == 0:
            off._end_label(ax, off.HATTERAS_ANNOTATIONS.low_end_label, False, y0, y1)
            off._north_arrow(ax, x=0.86, y=0.09)
        if i == len(groups) - 1:
            off._end_label(ax, off.HATTERAS_ANNOTATIONS.high_end_label, True, y0, y1)
        off._scalebar(ax, length_m=500.0)
        # Labels inside each box's landward (sound-side) edge: every changed
        # domain with its N, in its colour; every fifth unchanged domain as a
        # plain number, as the reference does. Inside the box rather than at
        # the panel edge because the boxes are staggered across the panel.
        for d in group:
            n = n_by.get(d, 0)
            b = gdf[gdf["domain_id"].astype(int) == d].total_bounds
            if n:
                ax.text(b[0] + 45.0, (b[1] + b[3]) / 2, f"{d}  {n:+d}",
                        fontsize=7.5, fontweight="bold", ha="left", va="center",
                        color=C_ADD if n > 0 else C_REM, zorder=6,
                        bbox=dict(facecolor="white", alpha=0.75, edgecolor="none",
                                  boxstyle="square,pad=0.1"))
            elif d % 5 == 0 or d in (ids.min(), ids.max()):
                ax.text(b[0] + 45.0, (b[1] + b[3]) / 2, str(d),
                        fontsize=7.5, ha="left", va="center", color=off.INK_MUTED,
                        zorder=6,
                        bbox=dict(facecolor="white", alpha=0.7, edgecolor="none",
                                  boxstyle="square,pad=0.1"))

    present = sorted(set(v for v in n_by.values() if v))
    handles = ([Patch(facecolor=reds[v - 1], alpha=0.45, edgecolor="0.4",
                      label=f"+{v} row{'s' if v > 1 else ''}") for v in present if v > 0]
               + [Patch(facecolor=blues[-v - 1], alpha=0.45, edgecolor="0.4",
                        label=f"\u2212{-v} row{'s' if v < -1 else ''}")
                  for v in sorted(present, key=abs) if v < 0]
               + [Patch(facecolor=C_ADD, label="added, true scale"),
                  Patch(facecolor=C_REM, label="removed, true scale")]
               + [Line2D([0], [0], label=f"{yr} dune line",
                         **dict(off.SIMPLE_LINE_STYLE[yr], linewidth=2.0)) for yr in (1984, 1997)]
               + (off.m.road_legend_handles(roads) if roads else [])
               + [Line2D([0], [0], color=off.HATTERAS_ANNOTATIONS.color_town_span, lw=5.0,
                         label="community")])
    fig.legend(handles=handles, loc="outside lower center",
               ncol=min(len(handles), 9), fontsize=8.5)

    p = FIG_DIR / "HAT_footprint_1984_plan.png"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


# =============================================================================
# FIGURE 3 - SIGNED N WITH SPREAD; SETBACK BEFORE AND AFTER
# =============================================================================

BLOCKS = (((9, 14), 1999), ((84, 87), 1989))      # the NC-12 relocation blocks


def _community_bands(ax) -> None:
    ann = off.HATTERAS_ANNOTATIONS
    for name, (lo, hi) in ann.town_spans.items():
        ax.axvspan(lo - 0.5, hi + 0.5, color="0.93", zorder=0)
        ax.text((lo + hi) / 2, 0.985, name, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=7.5, color=off.INK_MUTED)
    for (lo, hi), yr in BLOCKS:
        ax.axvspan(lo - .5, hi + .5, facecolor="none", edgecolor="#2c6e49", lw=0.9,
                   ls=(0, (3, 2)), zorder=1)
        ax.text((lo + hi) / 2, 0.94, f"relocated {yr}", transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=7, color="#2c6e49")


def fig_shift(tab: pd.DataFrame) -> Path:
    """The measured shift, its spread, and what the 10 m rule keeps of it."""
    off.apply_style()
    doms = tab.index.to_numpy()
    n = tab["n_cells"].to_numpy()
    med = tab["shift_m_median"].to_numpy()
    lo, hi = tab["shift_m_p10"].to_numpy(), tab["shift_m_p90"].to_numpy()
    fig, ax = plt.subplots(figsize=(13.0, 4.4), constrained_layout=True)
    _community_bands(ax)
    ax.bar(doms, n * CELL_M, width=0.82, color=np.where(n > 0, C_ADD_FILL, C_REM_FILL),
           edgecolor="none", zorder=2)
    ax.errorbar(doms, med, yerr=[med - lo, hi - med], fmt="o", ms=3.0, color=INK,
                ecolor="0.45", elinewidth=0.7, capsize=0, zorder=4)
    for yv in (-CELL_M, CELL_M):
        ax.axhline(yv, color="0.6", linewidth=0.6, linestyle="--", zorder=1)
    ax.axhline(0, color=INK, linewidth=0.6, zorder=1)
    ax.text(89.5, CELL_M + 1.5, "\u00b11 cell", ha="right", va="bottom", fontsize=7, color="0.45")
    ax.set_xlim(0.2, 90.8)
    ax.set_xticks([1] + list(range(10, 91, 10)))
    ax.set_xticks(list(doms), minor=True)
    ax.set_ylim(float(np.nanmin(lo)) - 6, float(np.nanmax(hi)) + 16)
    ax.set_xlabel("domain (1 = south, Cape Hatteras)")
    ax.set_ylabel("1984 line \u2212 1997 line (m)\n+ seaward (rows added)")
    ax.grid(axis="y", color="0.92", linewidth=0.5)
    ax.set_axisbelow(True)
    fig.legend(handles=[Line2D([0], [0], marker="o", ms=3.0, color=INK, linestyle="none",
                               label="median of 50 paired profiles, p10\u2013p90"),
                        Patch(facecolor=C_ADD_FILL, label="kept by the 10 m rule, N \u00d7 10 m (rows added)"),
                        Patch(facecolor=C_REM_FILL, label="kept by the 10 m rule (rows removed)"),
                        Line2D([0], [0], color="0.6", linestyle="--", label="\u00b11 cell"),
                        Line2D([0], [0], color="#2c6e49", lw=0.9, ls=(0, (3, 2)), label="NC-12 relocation block")],
               loc="outside lower center", ncol=5, fontsize=8)
    p = FIG_DIR / "HAT_footprint_1984_shift.png"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


def fig_setback(tab: pd.DataFrame) -> Path:
    """The road setback the model receives now, and under the new row 0."""
    off.apply_style()
    rd = tab.dropna(subset=["setback_new_m"])
    x = rd.index.to_numpy()
    now = rd["setback_model_now_m"].to_numpy()
    new = rd["setback_new_m"].to_numpy()
    p10, p90 = rd["setback_new_p10_m"].to_numpy(), rd["setback_new_p90_m"].to_numpy()
    nn = rd["n_cells"].to_numpy()
    fig, ax = plt.subplots(figsize=(13.0, 4.6), constrained_layout=True)
    _community_bands(ax)
    ax.vlines(x, now, new, color="0.75", linewidth=1.0, zorder=2)
    ax.plot(x, now, "o", ms=3.4, color="0.5", mfc="white", mew=1.0, zorder=3)
    ax.errorbar(x, new, yerr=[np.clip(new - p10, 0, None), np.clip(p90 - new, 0, None)],
                fmt="none", ecolor="0.6", elinewidth=0.8, capsize=0, zorder=3)
    for mask, c in (((nn > 0), C_ADD), ((nn < 0), C_REM), ((nn == 0), "0.3")):
        ax.plot(x[mask], new[mask], "o", ms=3.6, color=c, zorder=4)
    ax.axhline(0, color=INK, linewidth=0.6, zorder=1)
    ax.set_yscale("symlog", linthresh=50, linscale=1.2)
    ax.set_yticks([-20, 0, 10, 20, 50, 100, 200, 500])
    ax.set_yticklabels(["\u221220", "0", "10", "20", "50", "100", "200", "500"])
    ax.set_ylim(-32, 900)
    ax.set_ylabel("NC-12 setback (m landward of row 0)")
    ax.set_xlabel("domain (1 = south, Cape Hatteras)")
    ax.set_xlim(0.2, 90.8)
    ax.set_xticks([1] + list(range(10, 91, 10)))
    ax.set_xticks(list(rd.index), minor=True)
    ax.grid(axis="y", color="0.92", linewidth=0.5)
    ax.set_axisbelow(True)
    # the three domains worth reading off: the two floored at 0 today and the
    # one that loses most of its setback. Labels pushed sideways, clear of the
    # neighbours' spread bars.
    for d, dx, ha in ((16, 10, "left"), (85, -10, "right"), (86, 10, "left")):
        if d in rd.index:
            ax.annotate(f"{d}: {rd.loc[d, 'setback_model_now_m']:.0f} \u2192 {rd.loc[d, 'setback_new_m']:.0f} m",
                        xy=(d, rd.loc[d, "setback_new_m"]), xytext=(dx, 0),
                        textcoords="offset points", ha=ha, va="center", fontsize=7, color=INK,
                        bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.15"))
    fig.legend(handles=[Line2D([0], [0], marker="o", ms=3.4, color="0.5", mfc="white", mew=1.0, linestyle="none",
                               label="setback the model receives now (v2, floored at 0)"),
                        Line2D([0], [0], marker="o", ms=3.6, color=C_ADD, linestyle="none",
                               label="from the new row 0, rows added"),
                        Line2D([0], [0], marker="o", ms=3.6, color=C_REM, linestyle="none",
                               label="from the new row 0, rows removed"),
                        Line2D([0], [0], marker="o", ms=3.6, color="0.3", linestyle="none",
                               label="from the new row 0, unchanged"),
                        Line2D([0], [0], color="0.6", lw=0.8, label="p10\u2013p90 over the profiles"),
                        Line2D([0], [0], color="#2c6e49", lw=0.9, ls=(0, (3, 2)), label="NC-12 relocation block")],
               loc="outside lower center", ncol=6, fontsize=8)
    p = FIG_DIR / "HAT_footprint_1984_setback.png"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


# =============================================================================
# WORDS
# =============================================================================

def write_report(tab: pd.DataFrame, topo_name: str, figs: list[Path]) -> Path:
    add = tab[tab.n_cells > 0]
    rem = tab[tab.n_cells < 0]
    rd = tab.dropna(subset=["setback_new_m"])
    L = []
    w = L.append
    w("=" * 78)
    w("1984 FOOTPRINT - where the interior grows or shrinks to stand at the 1984 dune line")
    w("=" * 78)
    w(f"written        {datetime.now():%Y-%m-%d %H:%M}")
    w(f"topography     {PRODUCT} / {topo_name}   (read, never modified)")
    w(f"shift from     {SHIFT_DIR.relative_to(REPO)}  (per-profile files, PAIRED)")
    w(f"road from      {ROAD_DIR.relative_to(REPO)}  (per-profile file, same frame)")
    for f in figs:
        w(f"figure         {f.relative_to(REPO)}")
    w("")
    w("SCOPE ONLY. No array written, no fill chosen, no model-facing CSV touched.")
    w("")
    w("RULES (decided 2026-09-07)")
    w("  shift(p)   = (line_1997(p) - line_1984(p)) * 10 m, per profile; + = 1984 line seaward")
    w("  shift      = median over the 50 paired profiles (NOT the difference of two medians)")
    w("  N          = trunc(shift / 10 m)   symmetric: + rows added, - existing rows removed")
    w("  setback'   = median_p[ (road(p) - row0(p)) + shift(p) ]   row-0 convention, unrounded")
    w("")
    w("-" * 78)
    w("SUMMARY")
    w("-" * 78)
    w(f"  domains                          {len(tab)}")
    w(f"  rows ADDED     {len(add):2d} domains, {int(add.n_cells.sum()):3d} rows, largest +{int(add.n_cells.max())} (domain {int(add.n_cells.idxmax())})")
    w(f"  rows REMOVED   {len(rem):2d} domains, {int(-rem.n_cells.sum()):3d} rows, largest {int(rem.n_cells.min())} (domain {int(rem.n_cells.idxmin())})")
    w(f"  unchanged      {int((tab.n_cells == 0).sum()):2d} domains (|shift| < 10 m)")
    w(f"  residual |shift - 10 N|: median {tab.residual_m.abs().median():.1f} m, max {tab.residual_m.abs().max():.1f} m, always toward less change")
    w(f"  domains whose p10-p90 spread straddles zero: {int(tab.spread_straddles_zero.sum())} "
      f"(of which {int(((tab.spread_straddles_zero == 1) & (tab.n_cells != 0)).sum())} still get rows)")
    w("")
    w("  N        domains")
    for v in sorted(tab.n_cells.unique()):
        ds = tab.index[tab.n_cells == v].tolist()
        w(f"  {v:+3d}  {len(ds):3d}   {ds if v != 0 else ''}")
    w("")
    w("-" * 78)
    w("ROAD SETBACK, row-0 convention (metres landward of interior row 0)")
    w("-" * 78)
    w(f"  road domains measured            {len(rd)}")
    w(f"  new setback negative             {rd.index[rd.setback_new_m < 0].tolist()}")
    w(f"  new setback under one cell       {rd.index[(rd.setback_new_m >= 0) & (rd.setback_new_m < NEAR_ZERO_M)].tolist()}")
    w(f"  currently floored at 0           {rd.index[rd.setback_model_now_m == 0].tolist()}")
    w("  the raw-lines reading (road - 1984 line) is ~%.0f m larger, median, and is NOT used: the digitized"
      % float((rd.setback_raw84_m - rd.setback_new_m).median()))
    w("  lines trace the toe, ~19 m seaward of row 0, so it would change the model's reference frame.")
    w("")
    cols = ["n_cells", "shift_m_median", "shift_m_p10", "shift_m_p90", "residual_m",
            "setback_model_now_m", "setback_new_m", "setback_new_p10_m", "setback_new_p90_m", "flags"]
    w("  domain  N   shift   p10    p90   resid  sb_now  sb_new  p10    p90   flags")
    for d, r in rd.iterrows():
        w(f"  {d:6d} {r.n_cells:+3d} {r.shift_m_median:7.1f} {r.shift_m_p10:6.1f} {r.shift_m_p90:6.1f} "
          f"{r.residual_m:6.1f} {r.setback_model_now_m:7.1f} {r.setback_new_m:7.1f} "
          f"{r.setback_new_p10_m:6.1f} {r.setback_new_p90_m:6.1f}  {r['flags']}")
    w("")
    w("-" * 78)
    w("ASSUMPTIONS THIS CANNOT CHECK")
    w("-" * 78)
    for s in ("the 1984 and 1997 lines trace the same feature (1997 says 'light/dark break'; 1984 has no metadata)",
              "one integer per domain: the alongshore median stands for 50 profiles",
              "1997 stands for 1996: the surface is 1996 ALACE, the line a year later",
              "the 1984 dune crest equalled the 1996 one; the dune array is not re-estimated",
              "removal deletes SURVEYED 1996 beach/foredune cells that were ocean in 1984",
              "only island width moves; the ocean shoreline is the shoreline-offset input",
              "the 1984 road line is the 1978 export"):
        w(f"  * {s}")
    w("")
    w("At GIS 63-67 the 1997 line lies 90-130 m seaward of row 0: the pick sits well behind the")
    w("modern foredune there. N does not depend on it (line minus line); the interiors do.")
    p = SCOPE_DIR / "HAT_footprint_1984.txt"
    p.write_text("\n".join(L) + "\n", encoding="utf-8")
    return p


def write_captions(tab: pd.DataFrame, topo_name: str) -> None:
    add = tab[tab.n_cells > 0]
    rem = tab[tab.n_cells < 0]
    rd = tab.dropna(subset=["setback_new_m"])
    stats = (f"{len(add)} domains gain {int(add.n_cells.sum())} rows and {len(rem)} lose "
             f"{int(-rem.n_cells.sum())}; {int((tab.n_cells == 0).sum())} are unchanged. "
             f"N = trunc(median paired shift / 10 m), so a row appears only once a full cell of "
             f"change is measured; the largest are +{int(add.n_cells.max())} at GIS {int(add.n_cells.idxmax())} "
             f"and {int(rem.n_cells.min())} at GIS {int(rem.n_cells.idxmin())}.")
    sections = {
        "HAT_footprint_1984_grid.png":
            f"(a\u2013c) The 90 Barrier3D domains as the model indexes them, 50 alongshore cells each, every "
            f"cross-shore row down the page from the CURRENT interior row 0 at 0 on every domain "
            f"(topography {topo_name}; interiors run to 189 rows, white is off the array, cells are not "
            f"square; right axis in metres). Existing cells are shaded by elevation class (m above MHW), "
            f"so the dune ridge, the backbarrier flat and the sound-side marsh are distinguishable. Red rows "
            f"above 0 are the rows that would be added to bring row 0 to the 1984 dune line, blank because "
            f"no fill has been chosen; blue hatching marks existing rows 0 to |N|\u22121 that would be "
            f"removed where the island has prograded since 1984; the signed count is printed above each "
            f"changed domain. The black tick is where interior row 0 ends up; the dark bar is NC-12 at its "
            f"measured 1984 position (seaward edge, 20 m wide), which does not move, so its distance to the "
            f"tick is the new setback. Communities and villages along the top from the site "
            f"configuration. {stats}",
        "HAT_footprint_1984_rows.png":
            f"Rows per domain under the 10 m rule, positive where rows are added (the 1984 dune line lay "
            f"seaward of the 1997 line) and negative where existing rows are removed, with the communities "
            f"banded along the axis. {stats}",
        "HAT_footprint_1984_plan.png":
            f"(a–c) The same footprint in plan view, in three equal-aspect panels of thirty domains, "
            f"south at left, each showing its domain boxes (2000 × 500 m) whole. Grey relief is the "
            f"2009-2014-1996 DEM hillshaded at 10 m and carries no readable elevation. The 1984 (red) and "
            f"1997 (blue) dune lines, NC-12 in 1984 (dashed), and each domain box shaded by N, red where "
            f"rows are added and blue where they are removed, labelled on the landward edge. The solid "
            f"band is the same N at true scale, the 1997 dune line offset seaward (add) or landward "
            f"(remove) by N × 10 m, so one row is a 10 m sliver that needs a zoom; the band is anchored on "
            f"the 1997 line as the map proxy for the existing array's seaward edge, which itself sits "
            f"about 19 m seaward of interior row 0. The gap between the band's outer edge and the 1984 "
            f"line is the truncation to whole cells. Communities as a bracket in the ocean margin, the "
            f"pier and groin as seaward marks; scale bar 500 m = 50 cells. {stats}",
        "HAT_footprint_1984_shift.png":
            f"Per domain, the median of the 50 paired per-profile differences between the 1997 and 1984 "
            f"dune-line crossings (points, p10\u2013p90 bars), positive where the 1984 line lies seaward, "
            f"and the part of it the 10 m rule keeps as whole rows (filled bars, N \u00d7 10 m; red added, "
            f"blue removed). Dashed guides at \u00b11 cell; points inside them become no rows, and the "
            f"distance from a point to its bar is the truncation residual, always toward less change. "
            f"Communities banded along the axis; the two NC-12 relocation blocks outlined. {stats}",
        "HAT_footprint_1984_setback.png":
            f"For the {len(rd)} road domains, the NC-12 setback the model receives today (open circles, the "
            f"v2 measurement floored at 0) and the setback from the new row 0 (filled, coloured by what "
            f"happens to the domain, with the p10\u2013p90 over the profiles), both in metres landward of "
            f"interior row 0 on a symmetric-log axis. The new value is (road \u2212 row 0) + shift per "
            f"profile, then the median, unrounded. GIS 85 and 86, floored at 0 today, become "
            f"{rd.loc[85, 'setback_new_m']:.0f} and {rd.loc[86, 'setback_new_m']:.0f} m; GIS 16 falls to "
            f"{rd.loc[16, 'setback_new_m']:.0f} m after losing four rows. No new setback is negative.",
    }
    head = ("# Figure captions\n\nWritten by `HAT_footprint_1984.py` (`write_captions`). The figures "
            "carry no in-image titles or footnotes on purpose; use these under them.\n")
    existing = CAPTIONS.read_text(encoding="utf-8") if CAPTIONS.is_file() else head
    parts = existing.split("\n## ")
    keep = [parts[0]] + [p for p in parts[1:] if not any(p.startswith(f"`{k}`") for k in sections)]
    text = "\n## ".join(keep).rstrip() + "\n"
    for k, v in sections.items():
        text += f"\n## `{k}`\n\n{v}\n"
    CAPTIONS.write_text(text, encoding="utf-8")


# =============================================================================

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-plan", action="store_true", help="skip the DEM plan view")
    args = ap.parse_args()

    topo_dir, _dune_dir, topo_name = topo_dirs(PRODUCT)
    print(f"topography : {PRODUCT} / {topo_name}")
    prof = load_profiles()
    tab = by_domain(prof, topo_dir, topo_name)

    SCOPE_DIR.mkdir(parents=True, exist_ok=True)
    prof.round(2).to_csv(SCOPE_DIR / "footprint_1984_profiles.csv", index=False)
    tab.to_csv(SCOPE_DIR / "footprint_1984_by_domain.csv")
    print(f"wrote {SCOPE_DIR / 'footprint_1984_by_domain.csv'}")

    figs = [fig_grid(tab, topo_dir), fig_rows(tab), fig_shift(tab), fig_setback(tab)]
    if not args.no_plan:
        figs.append(fig_plan(tab))
    for f in figs:
        print(f"wrote {f}")
    rep = write_report(tab, topo_name, figs)
    write_captions(tab, topo_name)
    print(f"wrote {rep}\nwrote {CAPTIONS}")
    print(open(rep, encoding="utf-8").read().split("ROAD SETBACK")[0])


if __name__ == "__main__":
    main()
