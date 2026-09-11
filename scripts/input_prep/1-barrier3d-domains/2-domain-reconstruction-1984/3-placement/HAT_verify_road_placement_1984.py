#!/usr/bin/env python3
r"""
HAT_verify_road_placement_1984.py
==============================================================================
Is NC-12 placed where the 1984 road and dune lines say it was?  (Hannah,
2026-09-08: "help me ensure that the roadway is being placed correctly and
that the offset matches reality")

The check runs in three frames and they have to agree:

  MAP        the 1984 dune line and the 1984 NC-12 centreline on the 1 m
             lidar, in map metres, no Barrier3D processing. Two distances
             per domain: ALONG the extractor's profiles (raster rows, the
             frame the model indexes) and PERPENDICULAR (frame-free, nearest
             point on the road from samples along the dune line). They differ
             only by the obliquity of the profiles to the island.
  ROW 0      the same distance expressed as the model needs it: metres
             landward of interior row 0 (one cell behind the picked 1996
             crest). The 1984 dune line is a toe, ~19 m seaward of row 0,
             so the row-0 setback is the toe-to-road distance minus that
             per-profile feature term:
                 setback_new = (road - line84) - (row0 - line97)
                             = (road - row0) + (line97 - line84)      [identity]
             It is the SAME measurement; only the reference changes.
  MODEL      what v3 receives: RoadSetback_1984_dunestart.csv must equal
             setback_new_m; the road rows are int(setback/10) and +1; the
             truncation loses 0-10 m; the array must hold rows_before + N
             rows and the road must sit on it. With rows added the road as
             placed lies on the v2 cells N rows behind today's pavement
             (the block goes in behind it); with rows removed in front of
             the road it lies on the old pavement's first row or 0-2 rows
             seaward of it (`road_cells_offset`: the row count comes from
             the shift median, the setback from the median of the sum).
             bulldoze() overwrites the road rows with the road elevation
             every year, so which measured cells lie under the pavement
             does not change the run; the distance from the crest does.

OUTPUTS  2-domain-reconstruction-1984/3-placement/road_placement_check_1984.csv     per road domain
         2-domain-reconstruction-1984/3-placement/HAT_road_placement_check_1984.txt  the report
         figures/3-placement/behind-road/
             rows-{added,removed}/HAT_road_placement_check_GIS<N>.png
                 map (1 m lidar, the two
                 dune lines, NC-12, row 0, the model's road rows and the
                 block / removed rows in map space) beside the v3 grid
             HAT_road_placement_check_island.png   every road domain: the
                 map distance along profiles against perpendicular, and
                 the setback today / the 1984 one / the one the model holds

USAGE
    python HAT_verify_road_placement_1984.py                 # 85, 63, 49, 16
    python HAT_verify_road_placement_1984.py --domains 85,84
==============================================================================
"""
from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from shapely.geometry import Point


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "0-elevation" / "3-figures"))
from hat_topo_version import (insert_figures_dir, insert_figures_dir_for_domain,  # noqa: E402
                              array_name, dune_topo_root, duneline_shift_dir, topo_dirs, insert_scope_step)
from hat_figure_style import (C, DOMAIN_AXIS_LABEL, caption, elevation_cmap,     # noqa: E402
                              figsize, open_frame, save, town_bands)
import HAT_plot_duneline_offset as off                                          # noqa: E402

INIT = REPO / "data" / "hatteras_init"
PRODUCT = "1984-start"
SCOPE_DIR = INIT / "1-barrier3d-domains" / PRODUCT / "2-domain-reconstruction-1984"
FOOTPRINT_CSV = insert_scope_step(PRODUCT, "2-extent") / "footprint_1984_by_domain.csv"
STEP_DIR = insert_scope_step(PRODUCT, "3-placement")           # the check's table and report (2026-09-09)
ROAD_DIR = INIT / "4-mgmt-forcing" / "road_offset" / "dunestart_offset" / "1984"
SHIFT_DIR = duneline_shift_dir(PRODUCT)
FIG_DIR = insert_figures_dir(PRODUCT, "3-placement", "road-check")
BUILT = "v3"
CELL_M = 10.0
ROAD_ROWS = 2
HALF_M = 10.0                 # the geojson is a centreline; the model road is 20 m
BERM_EL_M = 1.7
SAMPLE_M = 10.0               # spacing of samples along the 1984 dune line
C_ADD, C_REM, C_ROAD, C_ROAD_OLD, INK = off.C_1984, off.C_1997, C["ROAD"], C["BASE"], off.INK


# =============================================================================
# THE MEASUREMENTS, PER PROFILE
# =============================================================================

def load_profiles() -> pd.DataFrame:
    road = pd.read_csv(ROAD_DIR / "RoadOffset_1984_profiles.csv")
    l84 = pd.read_csv(SHIFT_DIR / "duneline_shift_1984_profiles.csv")
    l97 = pd.read_csv(SHIFT_DIR / "duneline_shift_1997_profiles.csv")
    m = l84.merge(l97, on=["domain", "profile"], suffixes=("_84", "_97"))
    if not (m["interior_row0_cell_84"] == m["interior_row0_cell_97"]).all():
        raise SystemExit("row 0 differs between the 1984 and 1997 profile files")
    m = m.rename(columns={"interior_row0_cell_84": "row0"})
    m = m[["domain", "profile", "row0", "duneline_cell_84", "duneline_cell_97"]]
    r = road[["domain", "profile", "interior_row0_cell", "road_seaward_cell", "road_landward_cell",
              "interior_x", "interior_y", "road_x"]]
    m = m.merge(r, on=["domain", "profile"], how="left")
    has = m["road_seaward_cell"].notna()
    if not (m.loc[has, "interior_row0_cell"] == m.loc[has, "row0"]).all():
        raise SystemExit("row 0 differs between the road and dune-line profile files")
    # rows RELATIVE to interior row 0 (negative = seaward); map x = interior_x - row * 10
    m["r_line84"] = m["duneline_cell_84"] - m["row0"]
    m["r_line97"] = m["duneline_cell_97"] - m["row0"]
    m["r_road"] = m["road_seaward_cell"] - m["row0"]
    m["raw84_m"] = (m["r_road"] - m["r_line84"]) * CELL_M          # toe -> road, on the map, along the profile
    m["feat97_m"] = -m["r_line97"] * CELL_M                         # the 1997 line seaward of row 0 (+)
    m["shift_m"] = (m["r_line97"] - m["r_line84"]) * CELL_M          # + = 1984 line seaward = rows added
    m["sb_v2_m"] = m["r_road"] * CELL_M
    m["sb_new_m"] = m["sb_v2_m"] + m["shift_m"]
    # identity: sb_new + feat97 == raw84, per profile, exactly
    ok = m["road_seaward_cell"].notna()
    assert np.allclose(m.loc[ok, "sb_new_m"] + m.loc[ok, "feat97_m"], m.loc[ok, "raw84_m"])
    return m


def perpendicular(gdf, line84, road84, d: int) -> dict:
    """Nearest distance from samples along the 1984 dune line (inside the
    domain box) to the 1984 NC-12 centreline, minus the half width. Frame-free."""
    box = gdf[gdf["domain_id"].astype(int) == d].geometry.iloc[0]
    seg = line84.intersection(box)
    road = road84.intersection(box.buffer(400.0))
    if seg.is_empty or road.is_empty:
        return dict(perp_med_m=np.nan, perp_p10_m=np.nan, perp_p90_m=np.nan, n_perp=0)
    parts = list(getattr(seg, "geoms", [seg]))
    dist = []
    for g in parts:
        if g.geom_type != "LineString" or g.length < SAMPLE_M:
            continue
        for s in np.arange(0.0, g.length, SAMPLE_M):
            p = g.interpolate(s)
            dist.append(p.distance(road) - HALF_M)
    if not dist:
        return dict(perp_med_m=np.nan, perp_p10_m=np.nan, perp_p90_m=np.nan, n_perp=0)
    a = np.asarray(dist)
    return dict(perp_med_m=round(float(np.median(a)), 1), perp_p10_m=round(float(np.percentile(a, 10)), 1),
                perp_p90_m=round(float(np.percentile(a, 90)), 1), n_perp=int(a.size))


def read_setback_csv(p: Path) -> dict:
    rows = list(csv.reader(open(p, newline="")))
    ids = [int(float(x)) for x in rows[0] if x.strip()]
    vals = [float(x) for x in rows[1] if x.strip()]
    return dict(zip(ids, vals))


# =============================================================================
# THE CHECK, PER DOMAIN
# =============================================================================

def check(tab: pd.DataFrame, prof: pd.DataFrame, gdf, line84, road84, topo_v2: Path, topo_v3: Path,
          sb_csv: dict) -> pd.DataFrame:
    rows = []
    for d, g in prof.groupby("domain"):
        rd = g.dropna(subset=["road_seaward_cell"])
        if not len(rd) or d not in tab.index:
            continue
        t = tab.loc[d]
        n = int(t["n_cells"])
        rec = {"domain": int(d), "n_cells": n, "n_road_profiles": len(rd)}
        rec["raw84_along_m"] = round(float(np.median(rd["raw84_m"])), 1)
        rec["feat97_m"] = round(float(np.median(rd["feat97_m"])), 1)
        rec["shift_m"] = round(float(np.median(g["shift_m"])), 1)
        rec["setback_v2_m"] = round(float(np.median(rd["sb_v2_m"])), 1)
        rec["setback_new_m"] = round(float(np.median(rd["sb_new_m"])), 1)
        assert abs(rec["setback_new_m"] - float(t["setback_new_m"])) < 0.06, (d, "footprint table disagrees")
        # the two medians do not add: the row count follows the shift median, the
        # setback the median of the per-profile sum
        rec["median_nonadditivity_m"] = round(rec["setback_new_m"] - (rec["setback_v2_m"] + rec["shift_m"]), 1)
        rec.update(perpendicular(gdf, line84, road84, int(d)))
        rec["perp_minus_along_m"] = (round(rec["perp_med_m"] - rec["raw84_along_m"], 1)
                                     if np.isfinite(rec["perp_med_m"]) else np.nan)
        # the model
        r_v2 = int(rec["setback_v2_m"] // CELL_M)
        r_new = int(rec["setback_new_m"] // CELL_M)
        rec["road_row_v2"] = r_v2
        rec["road_row_new"] = r_new
        rec["model_setback_m"] = r_new * CELL_M
        rec["truncation_m"] = round(rec["setback_new_m"] - r_new * CELL_M, 1)
        rec["csv_setback_m"] = sb_csv.get(int(d), np.nan)
        rec["csv_matches"] = int(abs(rec["csv_setback_m"] - rec["setback_new_m"]) < 0.06)
        z2 = np.load(topo_v2 / array_name("topography", d))
        z3 = np.load(topo_v3 / array_name("topography", d))
        rec["rows_v2"], rec["rows_v3"] = int(z2.shape[0]), int(z3.shape[0])
        rec["rows_match"] = int(z3.shape[0] == z2.shape[0] + n)
        rec["road_on_array"] = int(r_new + ROAD_ROWS + 1 < z3.shape[0])   # bulldoze reads road_end + 1
        if n > 0:
            ins = int(t["insert_row_behind_road"])
            rec["v2_cells_under_road"] = f"{r_new}..{r_new + 1}"
            rec["pavement_offset_cells"] = r_new - r_v2       # the placed road is this many v2 rows behind today's pavement
            rec["block_rows_v3"] = f"{ins}..{ins + n - 1}"
        elif n < 0:
            ins = int(t["insert_row_behind_road"])
            rec["v2_cells_under_road"] = f"{r_new - n}..{r_new - n + 1}"
            rec["pavement_offset_cells"] = int(t["road_cells_offset"])   # old pavement's first row minus the placed road row
            rec["block_rows_v3"] = f"{ins}..{ins - n - 1} removed"
        else:
            rec["v2_cells_under_road"] = f"{r_new}..{r_new + 1}"
            rec["pavement_offset_cells"] = r_new - r_v2
            rec["block_rows_v3"] = ""
        flags = []
        if not rec["csv_matches"]:
            flags.append("CSV_NE_SETBACK")
        if not rec["rows_match"]:
            flags.append("ROWS_NE_V2_PLUS_N")
        if not rec["road_on_array"]:
            flags.append("ROAD_OFF_ARRAY")
        if np.isfinite(rec["perp_minus_along_m"]) and abs(rec["perp_minus_along_m"]) > CELL_M:
            flags.append("PERP_VS_ALONG>10m")
        if abs(rec["median_nonadditivity_m"]) >= CELL_M:
            flags.append("MEDIANS_DISAGREE>=1cell")
        if n < 0 and rec["pavement_offset_cells"] >= 2:
            flags.append("ROAD_2ROWS_SEAWARD_OF_PAVEMENT")
        rec["flags"] = ",".join(flags)
        rows.append(rec)
    return pd.DataFrame(rows).set_index("domain").sort_index()


# =============================================================================
# FIGURE, ONE DOMAIN: MAP BESIDE GRID
# =============================================================================

def _band(ax, pr: pd.DataFrame, r0: float, r1: float, **kw):
    """A cross-shore band over rows r0..r1 (relative to row 0, inclusive) on
    every profile of the domain, as one polygon in map space."""
    pr = pr.sort_values("interior_y")
    y = pr["interior_y"].to_numpy()
    xa = pr["interior_x"].to_numpy() - r0 * CELL_M + CELL_M / 2
    xb = pr["interior_x"].to_numpy() - r1 * CELL_M - CELL_M / 2
    # profiles are raster rows 10 m tall: extend the polygon to the cell edges
    yy = np.concatenate([[y[0] - CELL_M / 2], y, [y[-1] + CELL_M / 2]])
    xa = np.concatenate([[xa[0]], xa, [xa[-1]]])
    xb = np.concatenate([[xb[0]], xb, [xb[-1]]])
    poly = np.concatenate([np.column_stack([xa, yy]), np.column_stack([xb[::-1], yy[::-1]])])
    ax.add_patch(plt.Polygon(poly, closed=True, **kw))


def fig_domain(d: int, chk: pd.DataFrame, tab: pd.DataFrame, prof: pd.DataFrame, gdf, lines, roads,
               topo_v3: Path, dune_v3: Path) -> Path:
    off.apply_style()
    c = chk.loc[d]
    t = tab.loc[d]
    n = int(c["n_cells"])
    pr = prof[(prof["domain"] == d) & prof["road_seaward_cell"].notna()].copy()
    r_v2, r_new = int(c["road_row_v2"]), int(c["road_row_new"])
    ins = int(t["insert_row_behind_road"]) if n != 0 else None
    box = gdf[gdf["domain_id"].astype(int) == d].geometry.iloc[0]

    fig = plt.figure(figsize=figsize("double", height=5.4), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0])
    ax = fig.add_subplot(gs[0, 0])
    ag = fig.add_subplot(gs[0, 1])

    # ---- (a) the map --------------------------------------------------------
    arr, extent = off.load_1m(gdf, [d])
    if arr is not None:
        off._hillshade(ax, arr, extent, res=1.0)
    gdf[gdf["domain_id"].astype(int) == d].boundary.plot(ax=ax, color="0.3", linewidth=0.8, zorder=4)
    drawn = off.clip_for_drawing(lines, box.buffer(30.0))
    off.draw_lines(ax, drawn, scale=1.0, style=off.LINE_STYLE)
    off.m.draw_roads(ax, roads, scale=0.9)
    # interior row 0 on every profile
    pr_s = pr.sort_values("interior_y")
    ax.plot(pr_s["interior_x"], pr_s["interior_y"], color=INK, lw=1.2, ls=(0, (4, 2)), zorder=6)
    # today's pavement rows (v2 frame = the map)
    _band(ax, pr, r_v2, r_v2 + 1, facecolor="none", edgecolor=C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)), zorder=7)
    # the model's road AS PLACED, on the v2 cells it actually covers: v3 row r
    # is v2 row r in front of the seam and r + |N| behind it, so with rows
    # removed the two road rows can straddle the seam and sit apart on the map
    def v2_row(r):
        return r if (n >= 0 or r < ins) else r - n
    for rr in (r_new, r_new + 1):
        _band(ax, pr, v2_row(rr), v2_row(rr), facecolor=C_ROAD, edgecolor="none", alpha=0.45, zorder=7)
    r_map = v2_row(r_new)
    if n > 0:
        _band(ax, pr, r_new + ROAD_ROWS, r_new + ROAD_ROWS + n - 1, facecolor=C_ADD, edgecolor=C_ADD,
              alpha=0.35, hatch="////", lw=0.8, zorder=6)
    elif n < 0:
        _band(ax, pr, ins, ins - n - 1, facecolor=C_REM, edgecolor=C_REM, alpha=0.35, hatch="////",
              lw=0.8, zorder=6)
    # name the band, at the top of the panel above its cross-shore position
    top = pr["interior_y"].max() + CELL_M / 2
    if n > 0:
        xb = float(pr["interior_x"].median()) - (r_new + ROAD_ROWS + n / 2 - 0.5) * CELL_M
        ax.text(xb, top - 12, f"+{n} rows", ha="center", va="top",
                fontsize=7, color=C_ADD, fontweight="bold", zorder=9,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.15"))
    elif n < 0:
        xb = float(pr["interior_x"].median()) - (ins + (-n) / 2 - 0.5) * CELL_M
        ax.text(xb, top - 12, f"−{-n} rows", ha="center", va="top",
                fontsize=7, color=C_REM, fontweight="bold", zorder=9,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.15"))
    # distance arrows on the middle profile
    mid = pr_s.iloc[len(pr_s) // 2]
    y_mid = float(mid["interior_y"])
    x0 = float(mid["interior_x"])
    x_l84 = x0 - float(mid["r_line84"]) * CELL_M
    x_road_v2 = x0 - r_v2 * CELL_M + CELL_M / 2          # seaward edge of today's pavement
    x_road_new = x0 - r_map * CELL_M + CELL_M / 2        # seaward edge of the road as placed
    dy = 0.0
    x_l97 = x0 - float(mid["r_line97"]) * CELL_M
    # window: the seaward ~700 m of the box around the road (set before the
    # labels, which need the edges)
    b = box.bounds
    x_hi = min(b[2], max(x_l84, x_l97, float(pr["interior_x"].max())) + 120.0)
    x_lo = max(b[0], min(x_road_new, x_road_v2) - max(abs(n), 2) * CELL_M - 220.0)
    shift = float(c["shift_m"])                       # + = the 1984 line seaward of the 1997 line = rows added
    sb_v2, sb_new = float(c["setback_v2_m"]), float(c["setback_new_m"])
    arith = (f"{sb_v2:.0f} {'+' if shift >= 0 else '\u2212'} {abs(shift):.0f} "
             f"{'=' if abs(sb_v2 + shift - sb_new) < 0.6 else '\u2248'} {sb_new:.0f} m")
    c_shift = C_ADD if n > 0 else C_REM
    for (xa, xb, lab, col, k) in (
            (x_l84, x_road_v2, f"1984 dune line \u2192 NC-12: {c['raw84_along_m']:.0f} m", C_ADD, 0),
            (x_l97, x_l84, f"dune-line shift: {abs(shift):.0f} m "
                           f"{'seaward' if shift > 0 else 'landward'} \u2192 {n:+d} rows", c_shift, 1),
            (x0, x_road_v2, f"setback measured on the 1996 surface: {sb_v2:.0f} m", C_ROAD_OLD, 2),
            (x0, x_road_new, f"1984 setback (model input): {sb_new:.0f} m", C_ROAD, 3)):
        yy = y_mid + (1.5 - k) * 58.0
        ax.annotate("", xy=(xb, yy), xytext=(xa, yy),
                    arrowprops=dict(arrowstyle="<->", color=col, lw=1.4, shrinkA=0, shrinkB=0), zorder=9)
        # the label sits over its arrow's middle unless that is near a panel
        # edge, where it hangs off the arrow's inner end instead of spilling out
        xm = (xa + xb) / 2
        if xm > x_lo + 0.62 * (x_hi - x_lo):
            xt, ha = max(xa, xb), "right"
        elif xm < x_lo + 0.38 * (x_hi - x_lo):
            xt, ha = min(xa, xb), "left"
        else:
            xt, ha = xm, "center"
        ax.text(xt, yy + 6, lab, ha=ha, va="bottom", fontsize=7, color=col, fontweight="bold",
                zorder=9, bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.15"))
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(b[1] - 10, b[3] + 10)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    off._scalebar(ax, length_m=100.0)
    off._north_arrow(ax, x=0.94, y=0.12)
    ax.text(0.99, 0.5, "ocean", transform=ax.transAxes, ha="right", va="center", fontsize=8,
            color=off.INK_MUTED, rotation=90)
    off._title(ax, 0, f"GIS {d}: the lines on the 1 m lidar")

    # ---- (b) the v3 grid ----------------------------------------------------
    cmap, norm, bounds = elevation_cmap()
    z = np.load(topo_v3 / array_name("topography", d)) * CELL_M
    dune = np.load(dune_v3 / array_name("dune", d)) * CELL_M + BERM_EL_M
    img = np.concatenate([np.tile(dune[None, :], (ROAD_ROWS, 1)), z], axis=0)
    R = min(img.shape[0], max(40, r_new + abs(n) + 16) + ROAD_ROWS)
    ncol = z.shape[1]
    ag.imshow(img[:R], cmap=cmap, norm=norm, aspect="auto", interpolation="nearest", origin="upper",
              extent=[-0.5, ncol - 0.5, R - 0.5, -0.5])
    ag.axhline(ROAD_ROWS - 0.5, color=INK, lw=1.0, zorder=5)
    ag.text(ncol - 1.0, ROAD_ROWS / 2 - 0.5, "dune", ha="right", va="center", fontsize=7, fontweight="bold",
            color="white", zorder=6)
    y_road = r_new + ROAD_ROWS
    ag.add_patch(Rectangle((-0.5, y_road - 0.5), ncol, ROAD_ROWS, facecolor=C_ROAD, edgecolor=C_ROAD, lw=1.2,
                           alpha=0.45, zorder=5))
    ag.text(ncol - 1.0, y_road + ROAD_ROWS / 2 - 0.5, f"NC-12, {c['setback_new_m']:.0f} m",
            ha="right", va="center", fontsize=7, fontweight="bold", color="white", zorder=6)
    if n > 0:
        ag.add_patch(Rectangle((-0.5, ins + ROAD_ROWS - 0.5), ncol, n, facecolor="none", edgecolor=C_ADD, lw=1.6,
                               zorder=5))
        ag.text(0.5, ins + ROAD_ROWS + n - 0.3, f"+{n} rows, copied from landward",
                ha="left", va="top", fontsize=7, color=C_ADD, fontweight="bold", zorder=6,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.15"))
        y_pav = r_v2 + ROAD_ROWS
    elif n < 0:
        y_seam = ins + ROAD_ROWS - 0.5
        ag.axhline(y_seam, xmax=0.60, color=C_REM, lw=2.0, ls=(0, (3, 1.5)), zorder=8)
        seam_txt = (f"seam: rows {ins}–{ins - n - 1} removed" if ins > 0 else
                    f"seam at row 0: rows 0\u2013{-n - 1} removed")
        ag.annotate(seam_txt,
                    xy=(ncol * 0.35, y_seam), xytext=(ncol * 0.35, y_road + ROAD_ROWS + 4.0),
                    ha="center", va="top", fontsize=7, color=C_REM, fontweight="bold", zorder=9,
                    arrowprops=dict(arrowstyle="-|>", color=C_REM, lw=1.2, mutation_scale=12, shrinkB=0),
                    bbox=dict(facecolor="white", alpha=0.9, edgecolor=C_REM, lw=0.8, boxstyle="square,pad=0.25"))
        y_pav = r_v2 + n + ROAD_ROWS
    else:
        y_pav = r_v2 + ROAD_ROWS
    ag.add_patch(Rectangle((-0.5, y_pav - 0.5), ncol, ROAD_ROWS, facecolor="none", edgecolor=C_ROAD_OLD, lw=1.0,
                           ls=(0, (2, 2)), zorder=5))
    ag.set_xlabel("alongshore cell")
    ag.set_ylabel("cross-shore cell (0 = the dune)")
    ag.set_yticks(range(0, R, 10))
    sec = ag.secondary_yaxis("right", functions=(lambda cc: (cc - ROAD_ROWS) * CELL_M,
                                                  lambda mm: mm / CELL_M + ROAD_ROWS))
    sec.set_ylabel("m landward of interior row 0")
    for sp in ("top", "right"):
        ag.spines[sp].set_visible(True)
    off._title(ag, 1, f"GIS {d}: the model domain")

    labels = ["< 0 m (water)"] + [f"{lo:g}–{hi:g} m" for lo, hi in zip(bounds[1:-2], bounds[2:-1])] \
        + [f"> {bounds[-2]:g} m"]
    handles = [Patch(facecolor=cmap(i), edgecolor="0.4", lw=0.4, label=lab) for i, lab in enumerate(labels)]
    handles += [Line2D([0], [0], **dict(off.LINE_STYLE[1984], linewidth=2.0), label="1984 dune line"),
                Line2D([0], [0], **dict(off.LINE_STYLE[1997], linewidth=2.0), label="1997 dune line"),
                Line2D([0], [0], color="black", lw=1.6, ls=(0, (4, 2)), label="NC-12 (1984 line)"),
                Line2D([0], [0], color=INK, lw=1.2, ls=(0, (4, 2)), label="interior row 0 (crest + 1)"),
                Patch(facecolor=C_ROAD, alpha=0.45, edgecolor=C_ROAD, label="NC-12 road rows at the 1984 setback (model input)"),
                Patch(facecolor="none", edgecolor=C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)),
                      label="NC-12 road rows at the setback measured on the 1996 surface"),
                Patch(facecolor=C_ADD, alpha=0.35, edgecolor=C_ADD, hatch="////", label="rows inserted landward of NC-12"),
                Patch(facecolor=C_REM, alpha=0.35, edgecolor=C_REM, hatch="////", label="rows removed seaward of NC-12"),
                Line2D([0], [0], color=C_REM, lw=2.0, ls=(0, (3, 1.5)), label="seam left by the removal")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, fontsize=7, frameon=False,
               title="interior elevation (m above MHW)", title_fontsize=7)
    caption(fig, f"GIS {d}. (a) The 1984 dune line and the 1984 NC-12 centreline on the 1 m lidar, "
                 f"south at the bottom, ocean to the right; no coordinate ticks, scale bar 100 m. Interior "
                 f"row 0 (one cell behind the picked 1996 crest) is the dashed black line on every profile. "
                 f"The arrows are the measurements the check compares: the dune line to the road "
                 f"({c['raw84_along_m']:.0f} m along the extractor's profiles, {c['perp_med_m']:.0f} m "
                 f"measured perpendicular), the 1984-to-1997 dune-line shift ({shift:+.0f} m, which the 10 m "
                 f"rule turns into {n:+d} rows), the setback measured on the 1996 surface "
                 f"({sb_v2:.0f} m from row 0) and the 1984 setback the model receives "
                 f"({sb_v2:.0f} {'+' if shift >= 0 else '−'} {abs(shift):.0f} = {sb_new:.0f} m, "
                 f"which the model resolves to whole rows {r_new}–{r_new + 1}). "
                 f"(b) The same domain as the model holds it: two dune rows on top drawn at berm + dune "
                 f"height, then every interior row, in elevation classes; NC-12 at the 1984 setback filled "
                 f"and at the setback measured on the 1996 surface outlined; the rows the footprint inserts "
                 f"or removes marked. bulldoze() overwrites the road rows every year, so which cells lie "
                 f"under the pavement does not change a run — the distance from the crest does.")
    p = insert_figures_dir_for_domain(PRODUCT, "3-placement", d, under="road-check") / f"HAT_road_placement_check_GIS{d}.png"
    save(fig, p, vector=False, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


# =============================================================================
# FIGURE, EVERY ROAD DOMAIN
# =============================================================================

def fig_island(chk: pd.DataFrame) -> Path:
    off.apply_style()
    doms = chk.index.to_numpy()
    fig, (a, b) = plt.subplots(2, 1, figsize=figsize("double", height=5.2), sharex=True,
                               constrained_layout=True,
                               gridspec_kw=dict(height_ratios=[0.8, 1.2]))
    town_bands(a)
    town_bands(b, label=False)
    # (a) the map distance two ways
    a.axhspan(-CELL_M, CELL_M, color="0.96", zorder=0)
    a.axhline(0, color=INK, lw=0.6)
    a.plot(doms, chk["perp_minus_along_m"], "o", ms=3.8, color=C_ADD, mec="white", mew=0.4, zorder=4,
           label="perpendicular minus along-profile distance, 1984 dune line to NC-12")
    a.plot(doms, chk["median_nonadditivity_m"], "s", ms=3.4, color=C_REM, mec="white", mew=0.4, zorder=4,
           label="setback median minus (measured median + shift median)")
    a.set_ylabel("m")
    a.set_ylim(-32, 32)
    a.grid(axis="y")
    a.set_axisbelow(True)
    open_frame(a)
    a.legend(loc="lower left", ncol=2, fontsize=7, frameon=False)
    off._title(a, 0, "two residuals against one cell (±10 m, shaded)")
    # (b) the setback in three forms
    for col, lab, mk, cc, dx in (("setback_v2_m", "measured on the 1996 surface (floored at 0)", "o", C["BASE"], -0.25),
                                 ("setback_new_m", "1984: measured + dune-line shift", "D", C_ADD, 0.0),
                                 ("model_setback_m", "1984, resolved to whole 10 m rows", "s", C_ROAD, 0.25)):
        v = chk[col].to_numpy().astype(float)
        if col == "setback_v2_m":
            v = np.maximum(v, 0.0)
        b.plot(doms + dx, v, mk, ms=3.6, color=cc, mec="white", mew=0.4, linestyle="none", zorder=4, label=lab)
    b.set_yscale("symlog", linthresh=50, linscale=1.0)
    b.set_yticks([0, 10, 20, 50, 100, 200, 500])
    b.set_yticklabels(["0", "10", "20", "50", "100", "200", "500"])
    b.set_ylabel("NC-12 setback\n(m landward of row 0)")
    b.set_xlabel(DOMAIN_AXIS_LABEL)
    b.set_xlim(0.2, 90.8)
    b.set_xticks([1] + list(range(10, 91, 10)))
    b.set_xticks(doms, minor=True)
    b.grid(axis="y")
    b.set_axisbelow(True)
    open_frame(b)
    fig.legend(handles=b.get_legend_handles_labels()[0], loc="outside lower center", ncol=3,
               fontsize=7, frameon=False, title="NC-12 setback", title_fontsize=7)
    off._title(b, 1, "NC-12 setback, three ways")
    caption(fig, "Every domain with a model road, south at left; villages banded. (a) Two residuals that "
                 "would show a frame problem, against one Barrier3D cell (±10 m, shaded): the 1984 dune "
                 "line to NC-12 measured perpendicular minus the same distance along the extractor's "
                 "profiles, and the median of the per-profile setbacks minus the sum of the two medians "
                 "(the row count follows the shift median, the setback the median of the sum, so they need "
                 "not agree). (b) The NC-12 setback in metres landward of interior row 0 on a "
                 "symmetric-log axis: as measured on the 1996 surface, the 1984 value the model receives "
                 "(measured + dune-line shift), and that value resolved to whole 10 m rows, which "
                 "truncates toward the crest by 0–10 m.")
    p = FIG_DIR / "HAT_road_placement_check_island.png"
    save(fig, p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


# =============================================================================

def write_report(chk: pd.DataFrame, figs: list[Path]) -> Path:
    L = []
    w = L.append
    w("=" * 78)
    w("NC-12 PLACEMENT CHECK, 1984 - the lines on the lidar against what v3 holds")
    w("=" * 78)
    w(f"written   {datetime.now():%Y-%m-%d %H:%M}")
    w(f"inputs    {ROAD_DIR.relative_to(REPO)}/RoadOffset_1984_profiles.csv")
    w(f"          {SHIFT_DIR.relative_to(REPO)}/duneline_shift_{{1984,1997}}_profiles.csv")
    w(f"          {FOOTPRINT_CSV.relative_to(REPO)}")
    w(f"          dune-topo/{BUILT}/RoadSetback_1984_dunestart.csv, dune-topo/{BUILT}/topography")
    w("          raw-duneline-geojson/duneline_1984.geojson, raw_offset/1984/nc12_1984.geojson")
    for f in figs:
        w(f"figure    {f.relative_to(REPO)}")
    w("")
    w("WHAT IS CHECKED")
    w("  map      1984 dune line -> 1984 NC-12, in metres, ALONG the extractor's profiles and")
    w("           PERPENDICULAR (frame-free). Same number to within a cell = the profiles are")
    w("           not oblique enough to matter.")
    w("  row 0    setback_new = (road - line84) - (row0 - line97): the map distance minus the")
    w("           toe-to-row-0 feature term, per profile (identity, asserted).")
    w(f"  model    {BUILT}'s setback CSV == setback_new_m; rows == v2 rows + N; road rows on the")
    w("           array; truncation to the cell; which v2 cells lie under the road as placed.")
    w("")
    w("SUMMARY")
    w(f"  road domains checked                 {len(chk)}")
    w(f"  CSV setback == setback_new_m         {int(chk.csv_matches.sum())} of {len(chk)}")
    w(f"  rows == v2 + N                       {int(chk.rows_match.sum())} of {len(chk)}")
    w(f"  road rows on the array               {int(chk.road_on_array.sum())} of {len(chk)}")
    pa = chk["perp_minus_along_m"].dropna()
    w(f"  perpendicular - along, m             median {pa.median():+.1f}, p10 {pa.quantile(.1):+.1f}, "
      f"p90 {pa.quantile(.9):+.1f}, |max| {pa.abs().max():.1f} (GIS {int(pa.abs().idxmax())})")
    w(f"  feature term (1997 line vs row 0), m median {chk.feat97_m.median():.1f}, "
      f"IQR {chk.feat97_m.quantile(.25):.1f}..{chk.feat97_m.quantile(.75):.1f}")
    w(f"  truncation to the cell, m            median {chk.truncation_m.median():.1f}, max {chk.truncation_m.max():.1f} "
      f"(always toward the crest)")
    na = chk["median_nonadditivity_m"]
    w(f"  setback vs today+shift medians, m    median {na.median():+.1f}, |max| {na.abs().max():.1f} "
      f"(GIS {int(na.abs().idxmax())}); >= 1 cell at {chk.index[na.abs() >= CELL_M].tolist()}")
    add = chk[chk.n_cells > 0]
    rem = chk[chk.n_cells < 0]
    w(f"  rows added ({len(add)} domains)      road as placed sits N v2 rows behind today's pavement at "
      f"{int((add.pavement_offset_cells == add.n_cells).sum())} of {len(add)}; the block behind it")
    w(f"  rows removed ({len(rem)} domains)    road as placed vs the old pavement's first row: "
      + ", ".join(f"{int((rem.pavement_offset_cells == k).sum())} at {k} row{'s' if k != 1 else ''} seaward"
                  for k in sorted(rem.pavement_offset_cells.unique())))
    w(f"  flagged                              {chk.index[chk['flags'] != ''].tolist()}")
    w("")
    w("PER DOMAIN")
    w("  domain   N  line84->road  perp  feat97  sb_v2  sb_new  csv    row  model  trunc  v2 cells under  offset  flags")
    for d, r in chk.iterrows():
        w(f"  {d:6d} {r.n_cells:+3d}  {r.raw84_along_m:11.1f}  {r.perp_med_m:5.0f}  {r.feat97_m:6.1f}  "
          f"{r.setback_v2_m:5.0f}  {r.setback_new_m:6.1f}  {r.csv_setback_m:5.1f}  {r.road_row_new:3d}  "
          f"{r.model_setback_m:5.0f}  {r.truncation_m:5.1f}  {r.v2_cells_under_road:>14s}  {r.pavement_offset_cells:6d}  {r['flags']}")
    w("")
    w("READ WITH")
    w("  * the road rows are overwritten with the road elevation by roadway_manager.bulldoze()")
    w("    every year, so the cells under the pavement do not change the run; the distance")
    w("    from the crest (the setback) and the cells in front of and behind the road do.")
    w("  * where rows are added the road as placed is N v2 rows behind today's pavement on the")
    w("    MAP (the block goes in behind it), but at its 1984 distance from the crest in the")
    w("    model - the advisor's placement books the lost width on the sound side.")
    w("  * where rows are removed in front of the road the road keeps its map position to")
    w("    within the truncation: 0-2 rows seaward of the old pavement.")
    p = STEP_DIR / "HAT_road_placement_check_1984.txt"
    p.write_text("\n".join(L) + "\n", encoding="utf-8")
    return p


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--domains", default="85,63,49,16")
    args = ap.parse_args()
    tab = pd.read_csv(FOOTPRINT_CSV).set_index("domain")
    prof = load_profiles()
    topo_v2, _d2, name_v2 = topo_dirs(PRODUCT, override="v2")
    root = dune_topo_root(PRODUCT)
    topo_v3, dune_v3 = root / BUILT / "topography", root / BUILT / "dunes"
    if not topo_v3.is_dir():
        raise SystemExit(f"\n{topo_v3} not built - run HAT_build_footprint_version.py first\n")
    sb_csv = read_setback_csv(root / BUILT / "RoadSetback_1984_dunestart.csv")
    print("  loading the domain boxes, the dune lines and NC-12 ...")
    gdf = off.load_domains()
    lines = off.load_lines(gdf.crs)
    roads_all = off.m.load_roads(gdf.crs, clip_to=gdf.union_all().buffer(400.0))
    roads = {1984: roads_all[1984]}
    line84 = lines[1984].union_all()
    road84 = roads[1984].union_all()
    chk = check(tab, prof, gdf, line84, road84, topo_v2, topo_v3, sb_csv)
    chk.to_csv(STEP_DIR / "road_placement_check_1984.csv")
    print(f"wrote {SCOPE_DIR / 'road_placement_check_1984.csv'}")
    figs = []
    for d in [int(x) for x in args.domains.split(",") if x.strip()]:
        if d not in chk.index:
            print(f"  GIS {d}: no model road - skipped")
            continue
        figs.append(fig_domain(d, chk, tab, prof, gdf, lines, roads, topo_v3, dune_v3))
        print(f"wrote {figs[-1]}")
    figs.append(fig_island(chk))
    print(f"wrote {figs[-1]}")
    rep = write_report(chk, figs)
    print(f"wrote {rep}\n")
    print(open(rep, encoding="utf-8").read().split("PER DOMAIN")[0])


if __name__ == "__main__":
    main()
