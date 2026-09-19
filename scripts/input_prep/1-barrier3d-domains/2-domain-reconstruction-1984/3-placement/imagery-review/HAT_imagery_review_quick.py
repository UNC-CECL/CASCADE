#!/usr/bin/env python3
r"""
HAT_imagery_review_quick.py
==============================================================================
The basic review (Hannah's advisor, 2026-09-10): one domain at a time, the
1984 and 1997 photographs side by side, the 1984 ROAD OFFSET drawn on them
raw and as modeled, and two questions. Nothing to pick, nothing to toggle.

    offset_ok   yes | no | unclear   does the modeled 1984 road offset look
                                     right against the photographs?
    rows_ok     yes | no | unclear   does N, the rows the footprint adds or
                                     removes, look right?
    notes       free text

THE OFFSETS (decided with Hannah, 2026-09-10, revised the same afternoon)
    on the photographs, one arrow per panel in the year's colour: that year's
             dune line to the 1984 NC-12 line, both from the geojsons,
             measured along the 50 profiles (the line's crossing nearest row
             0, the road's crossing landward of it, moved HALF_M seaward to
             the pavement edge), median over the profiles. BOTH years
             reference the 1984 roadway (Hannah, 2026-09-10: "none of the
             offset ... should be in reference to the 2004 road position");
             the 2004 line is drawn for context only. What you would measure
             by hand on the map, year by year, against the road as it was.
    table    for the record in the text: `setback_raw84_m` from the footprint
             table, the same 1984 distance but to the road MASK in whole
             cells; it runs ~19 m larger than the model's reference because
             the digitized line traces the toe, a cell or two seaward of
             interior row 0.
    modeled  the setback CASCADE receives for v3 (`setback_new_m`, row-0
             convention, no floor) cut to whole cells: the road's first row is
             int(setback / 10) from interior row 0, and the pavement is the two
             rows from there. Drawn on the MODEL panel only, as the dark rows
             and an arrow from row 0 in rows.
    today    for reference only, in the text: the setback measured on the 1996
             surface (`setback_v2_m`, signed) and what the model holds for it
             (`setback_model_now_m`, floored at 0 where negative).
    Both arrows are drawn on the profile nearest the middle of the domain, so
    they start on the line they belong to; their LENGTH is the domain median,
    which is the number printed. The ~1.2 m stated accuracy of the
    georeferencing and the 10 m cell are the yardsticks.

WHAT IS ON THE PHOTOGRAPHS (settled in interview, Hannah, 2026-09-10, after
the overlays had grown to nine and "made it more confusing")
    Two verdicts, two comparisons, nothing else:
    * OFFSET. The 1984 NC-12 line and, on the middle profile, ONE white tick:
      interior row 0 measured seaward from the pavement edge - on (a) as the
      model places it for 1984 (the setback after its cut to cells, the same
      number panel (c) shows), on (b) as the DEM has it. The reviewer asks
      whether the crest visible in that year's photograph sits on the tick.
    * ROWS N. Both digitized dune lines (1984 red, 1997 blue) on both
      photographs, and one bracket between them on the middle profile
      labelled with the median shift and the rows it became. The reviewer
      asks whether each line follows the vegetation edge of its year and
      whether the retreat looks like the number.
    Plus the domain box. Removed from the photographs: interior row 0 as a
    line, the implied 1984 crest, the dune-line-to-road arrows, the row-0-to-
    road arrows, the dune search window, the 2004 road line. Their numbers
    stay in the side panel; the model-side marks stay on panel (c).

THE MODEL PANEL (upper right; the photographs stack down the left, 1984 over
1997, and the legend sits lower right - Hannah, 2026-09-10)
    The PROCESSED domain, the model input as the model holds it: the
    straightened Barrier3D grid of dune-topo/v3 (two dune rows at berm + dune
    height, then the interior, in the elevation classes of the placement
    check), cross-shore rows across with the ocean on the right and alongshore
    cells up the page, south at the bottom, so it faces the same way as the
    photographs. Nothing is mapped back through the shear: the dune is a
    straight band because that is what the model gets. With the model-side
    measurements on it: interior row 0, the road rows at the 1984 setback
    (dark), the rows the footprint inserted (red outline) or the seam it left
    (blue dashes), and the modeled offset as an arrow in rows. (The outline
    of today's pavement rows was dropped 2026-09-10: it is not part of the
    model input and read as a second road.) If dune-topo/v3 is not on disk the CURRENT version is
    drawn instead and the panel says so.

WHAT IT WRITES
    On Save: offset_ok, rows_ok, notes, reviewed_by, reviewed_at for the
    domain on screen into 2-domain-reconstruction-1984/3-placement/imagery-review/
    imagery_review_1984.csv (the sheet HAT_imagery_review_1984.py writes; that
    script keeps these columns when it re-runs). Nothing else in the sheet is
    touched, and the six older verdict columns are left as they are.
    "Save figure" writes the view on screen at 200 dpi to
        figures/3-placement/imagery-review/{rows-added,rows-removed,unchanged}/
            HAT_imagery_review_quick_GIS<N>.png
    with its caption in figures/CAPTIONS.md (one section for the set). The
    figure carries a header (domain, island section, N, the 1984 setback and
    its row), split legends for the photographs and the model input, and a
    source line (USGS release and DOI, stated accuracy, which road line each
    year is measured against).
    "Summarize" (or --summary) tallies the sheet:
        HAT_imagery_review_quick.txt                 counts, the "no" domains
        figures/3-placement/imagery-review/island/HAT_imagery_review_quick.png
            (a) N per domain along the island, coloured by rows_ok
            (b) the raw and modeled 1984 offset per domain, with offset_ok

KEYS   Right / Left  next / previous domain    Ctrl+S  save
       Space or B    flip the year in blink mode
       (keys are ignored while the cursor is in a text box)

Photographs are read through the same cache as the full window
(~/.cascade/imagery_review_cache/), so a domain seen in either is instant in
the other.

USAGE
    python HAT_imagery_review_quick.py                  # every domain in the sheet
    python HAT_imagery_review_quick.py --domains 85,63  # a subset
    python HAT_imagery_review_quick.py --summary        # tally only, no window
    python HAT_imagery_review_quick.py --smoke          # open, draw one, screenshot, close
==============================================================================
"""
from __future__ import annotations

import argparse
import getpass
import sys
import threading
import tkinter as tk
from tkinter import ttk
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import HAT_imagery_review_1984 as R  # noqa: E402   data, overlays, sheet
from HAT_imagery_review_gui import Data, CACHE_DIR, YEAR_COLOUR  # noqa: E402   the loader and its cache

TITLE_SIZE, TITLE_PAD = 10.5, 6                           # one title style for every panel


def _square_window(pr, pl, bounds):
    """The photographs' window made SQUARE about its centre, so that the panel
    boxes are equal squares (Hannah, 2026-09-10). Wraps the batch script's
    domain_window; the loader reads the tiles for this window and caches them
    under it, so the padding is real photograph, not blank."""
    x_lo, x_hi, y_lo, y_hi = _domain_window_orig(pr, pl, bounds)
    side = max(x_hi - x_lo, y_hi - y_lo)
    cx, cy = (x_lo + x_hi) / 2, (y_lo + y_hi) / 2
    return cx - side / 2, cx + side / 2, cy - side / 2, cy + side / 2


_domain_window_orig = R.domain_window
R.domain_window = _square_window                          # the loader looks it up on R at call time
from site_layer.hat_topo_version import insert_figures_dir, topo_dirs, array_name  # noqa: E402
from site_layer.hat_figure_style import elevation_cmap  # noqa: E402

off = R.off
BUILT = "v3"                      # the version the model panel shows; falls back to CURRENT
ROAD_ROWS = R.ROAD_ROWS
BERM_EL_M = 1.7                   # dune rows drawn at berm + dune height, as the placement check
VOCAB = ["", "yes", "no", "unclear"]
QUESTION = {"offset_ok": "1984 road offset OK?", "rows_ok": "rows N OK?"}
C_MODEL = "white"                 # the model's own measurements (row 0 → road)
VERDICT_C = {"yes": "#4d4d4d", "no": R.C_1984 if hasattr(R, "C_1984") else "#b2182b",
             "unclear": "#e08214"}
SUMMARY_PNG = "HAT_imagery_review_quick.png"
SUMMARY_TXT = R.STEP_DIR / "HAT_imagery_review_quick.txt"


def _num(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def _s(v) -> str:
    return v.strip() if isinstance(v, str) else ""


# =============================================================================
# THE OFFSETS OF ONE DOMAIN
# =============================================================================

def offsets(t: pd.Series) -> dict:
    """raw, modeled (row and metres) and today's, from the footprint table."""
    raw = _num(t.get("setback_raw84_m", np.nan))
    new = _num(t.get("setback_new_m", np.nan))
    v2 = _num(t.get("setback_v2_m", np.nan))
    now = _num(t.get("setback_model_now_m", np.nan))
    row = _num(t.get("road_row_new", np.nan))
    if not np.isfinite(row) and np.isfinite(new):
        row = float(int(new // R.CELL_M))            # the model's own cut, as HAT_footprint_1984.py
    row = int(row) if np.isfinite(row) else None
    return {"raw": raw, "new": new, "row": row,
            "model": row * R.CELL_M if row is not None else np.nan, "v2": v2, "now": now}


ROAD_REF = 1984                                          # every offset is to the 1984 NC-12 line (Hannah, 2026-09-10)
LINE_PAIR = {1984: (1984, ROAD_REF), 1997: (1997, ROAD_REF)}   # photo year -> (dune line year, NC-12 line year)


def line_offsets(pr: pd.DataFrame, lines: dict, roads: dict, box) -> dict:
    """The offset between that year's dune line and that year's NC-12, from the
    geojsons, along the 50 profiles: per profile the line's crossing nearest
    interior row 0 and the road's crossing nearest landward of it, the road
    centreline moved HALF_M seaward to the pavement's edge; median over the
    profiles. No cells, no mask: what the photograph's overlays show.
    Returns {photo year: dict(m, n, x_line, x_road, y)} with the middle
    profile's crossings for the arrow. Both years are measured to the 1984
    NC-12 line (ROAD_REF): the question is the 1984 roadway position, and the
    1997 arrow shows how far the dune had moved from it by then."""
    from shapely.geometry import LineString
    area = box.buffer(80.0)
    clipped = off.clip_for_drawing(lines, area)
    prs = pr.sort_values("interior_y")
    i_mid = int(np.argmin(np.abs(prs["interior_y"].to_numpy() - prs["interior_y"].median())))
    out = {}
    for y, (ly, ry) in LINE_PAIR.items():
        if ly not in clipped or ry not in roads or not len(clipped[ly]):
            continue
        lg = clipped[ly].union_all()
        rg = roads[ry].clip(area).union_all() if len(roads[ry]) else None
        if rg is None or rg.is_empty:
            continue
        vals, vals_r0, mid, mid_r0 = [], [], None, None
        for k, (_, p) in enumerate(prs.iterrows()):
            x0, yy = float(p["interior_x"]), float(p["interior_y"])
            seg = LineString([(x0 + 600.0, yy), (x0 - 1500.0, yy)])
            xs_l = [g.x for g in _points(lg.intersection(seg))]
            xs_r = [g.x for g in _points(rg.intersection(seg))]
            if not xs_l or not xs_r:
                continue
            x_l = min(xs_l, key=lambda v: abs(v - x0))
            land = [v for v in xs_r if v < x_l]
            x_r = max(land) if land else min(xs_r, key=lambda v: abs(v - x_l))
            x_r += R.HALF_M                                # centreline -> seaward pavement edge
            vals.append(x_l - x_r)
            vals_r0.append(x0 - x_r)                       # the model's reference to the same road edge
            if k == i_mid:
                mid, mid_r0 = (x_l, x_r, yy), (x0, x_r, yy)
        if vals:
            m = float(np.median(vals))
            if mid is None:                                # middle profile had no crossing: draw at the median row
                mid = (float(prs["interior_x"].median()), float(prs["interior_x"].median()) - m,
                       float(prs["interior_y"].median()))
            out[y] = dict(m=m, n=len(vals), x_line=mid[0], x_road=mid[0] - m, y=mid[2])
        if vals_r0 and "row0" not in out and ry == ROAD_REF:
            m0 = float(np.median(vals_r0))
            if mid_r0 is None:
                mid_r0 = (float(prs["interior_x"].median()), float(prs["interior_x"].median()) - m0,
                          float(prs["interior_y"].median()))
            out["row0"] = dict(m=m0, n=len(vals_r0), x_line=mid_r0[0], x_road=mid_r0[0] - m0, y=mid_r0[2])
    return out


def _points(g) -> list:
    if g.is_empty:
        return []
    if g.geom_type == "Point":
        return [g]
    if hasattr(g, "geoms"):
        return [q for gg in g.geoms for q in _points(gg)]
    return [g.representative_point()]                      # a line lying along the profile: one crossing


C_WINDOW = "#ffd92f"              # the dune search window, on the DEM panel only


def search_window(frames, d: int) -> tuple[int, int] | None:
    """The extractor's picked dune search window for one domain, (i0, i1) as
    cross-shore cells of the straightened profile array, i1 exclusive - the
    same cell frame as `row0` in the profile table. None if not picked."""
    try:
        _ro, _ext, windows = frames._extractor()
    except SystemExit as e:
        print(f"  search windows unavailable: {e}")
        return None
    w = windows.get(f"domain_{d}")
    if not w or "i0" not in w:
        return None
    return int(w["i0"]), int(w["i1"])


def draw_search_window(ax, pr: pd.DataFrame, win: tuple[int, int] | None) -> None:
    """The search window on the map: per profile the cells i0..i1-1 relative
    to that profile's row 0 (map x = interior_x - row x 10)."""
    if win is None:
        return
    i0, i1 = win
    prs = pr.sort_values("interior_y")
    y = prs["interior_y"].to_numpy()
    r0 = i0 - prs["row0"].to_numpy()                      # relative rows, per profile
    r1 = i1 - 1 - prs["row0"].to_numpy()
    xa = prs["interior_x"].to_numpy() - r0 * R.CELL_M + R.CELL_M / 2
    xb = prs["interior_x"].to_numpy() - r1 * R.CELL_M - R.CELL_M / 2
    yy = np.concatenate([[y[0] - R.CELL_M / 2], y, [y[-1] + R.CELL_M / 2]])
    xa = np.concatenate([[xa[0]], xa, [xa[-1]]])
    xb = np.concatenate([[xb[0]], xb, [xb[-1]]])
    poly = np.concatenate([np.column_stack([xa, yy]), np.column_stack([xb[::-1], yy[::-1]])])
    ax.add_patch(plt.Polygon(poly, closed=True, facecolor=C_WINDOW, alpha=0.22, edgecolor="none", zorder=3))
    ax.add_patch(plt.Polygon(poly, closed=True, facecolor="none", edgecolor=C_WINDOW, alpha=0.95, lw=1.0,
                             ls=(0, (2, 2)), zorder=3))


def draw_crest_tick(ax, year: int, meas: dict, o: dict, slot: int = 0) -> None:
    """ONE tick per photograph on the middle profile: interior row 0 measured
    seaward from the 1984 pavement edge. On the 1984 photograph it is the row
    the model places for 1984 (the setback after its cut to cells, the same
    number panel (c) shows); on the 1997 photograph it is the DEM's own row 0.
    The reviewer judges whether the crest visible in that photograph sits on
    the tick. A dotted connector from the road says what the metres refer to."""
    if "row0" not in meas:
        return
    e = meas["row0"]
    x_road, yy = e["x_road"], e["y"] - 20.0 - 40.0 * slot
    if year == 1984:
        if o["row"] is None:
            return
        dist, lab = o["model"], f"row 0, 1984 model: {o['model']:.0f} m"
    else:
        dist, lab = e["m"], f"row 0, DEM: {e['m']:.0f} m"
    x_tick = x_road + dist
    halo = [pe.withStroke(linewidth=3.0, foreground="black")]
    ax.plot([x_road, x_tick], [yy, yy], color=C_MODEL, lw=1.2, ls=(0, (1, 2)), zorder=11, path_effects=halo)
    ax.plot([x_tick, x_tick], [yy - 25.0, yy + 25.0], color=C_MODEL, lw=2.6, zorder=12, path_effects=halo)
    ax.text(x_tick, yy - 29.0, lab, ha="center", va="top", fontsize=9, color=C_MODEL, fontweight="bold",
            zorder=13, path_effects=halo)                    # below the tick; the bracket's label sits above


def draw_shift_bracket(ax, pr: pd.DataFrame, meas: dict, shift: float, n: int) -> None:
    """One bracket between the two dune lines on the middle profile, labelled
    with the median shift and the rows it became. Ends at each line's
    crossing of that profile; the LENGTH drawn is the domain median, which
    is the number printed, so the bracket and the label agree."""
    if not np.isfinite(shift):
        return
    prs = pr.sort_values("interior_y")
    i = int(np.argmin(np.abs(prs["interior_y"].to_numpy() - prs["interior_y"].median())))
    p = prs.iloc[i]
    x84 = meas[1984]["x_line"] if 1984 in meas else float(p["interior_x"] - p["r_line84"] * R.CELL_M)
    yy = float(p["interior_y"]) + 45.0
    x97 = x84 - shift
    halo = [pe.withStroke(linewidth=3.0, foreground="black")]
    ax.plot([x84, x97], [yy, yy], color=C_MODEL, lw=1.4, zorder=11, path_effects=halo)
    for xx in (x84, x97):
        ax.plot([xx, xx], [yy - 10.0, yy + 10.0], color=C_MODEL, lw=1.4, zorder=11, path_effects=halo)
    ax.text((x84 + x97) / 2, yy + 13.0, f"shift {shift:+.0f} m = {n:+d} rows", ha="center", va="bottom",
            fontsize=9, color=C_MODEL, fontweight="bold", zorder=13, path_effects=halo)


def offset_lines(d: int, t: pd.Series, o: dict, cover: dict, meas: dict | None = None) -> list[str]:
    n = int(t["n_cells"])
    meas = meas or {}
    L = [f"GIS {d}    N = {n:+d} rows    ({t['action']})",
         f"shift 1984->1997: {t['shift_m_median']:+.0f} m   (p10 {t['shift_m_p10']:+.0f}, p90 {t['shift_m_p90']:+.0f})",
         "",
         "ROAD OFFSET, dune line -> 1984 pavement edge (the photos)",
         "  negative = the dune line lies landward of the 1984 road"]
    for y, (ly, ry) in LINE_PAIR.items():
        if y in meas:
            L.append(f"  {y}     {meas[y]['m']:6.0f} m   {ly} dune line -> NC-12 {ry} line, "
                     f"{meas[y]['n']} profiles")
        else:
            L.append(f"  {y}         -     no crossing of both lines in this domain")
    if 1984 in meas and 1997 in meas:
        L.append(f"  change   {meas[1997]['m'] - meas[1984]['m']:+6.0f} m   1997 minus 1984")
    if "row0" in meas:
        L.append(f"  row 0    {meas['row0']['m']:6.0f} m   interior row 0 (1996 crest) -> NC-12 1984 line, "
                 f"on the map (panel b)")
        sh = _num(t.get("shift_m_median"))
        if np.isfinite(sh):
            L.append(f"  crest'84 {meas['row0']['m'] + sh:6.0f} m   implied 1984 crest (row 0 + shift) -> "
                     f"NC-12 (panel a; v3's assumption)")
    L += ["", "THE MODEL'S 1984 SETBACK, in three steps"]
    if o["row"] is not None and np.isfinite(o["v2"]):
        shift = _num(t.get("shift_m_median"))
        L += [f"  1 measured  {o['v2']:6.0f} m   1984 road mask -> interior row 0, on the DEM",
              f"  2 + shift   {shift:+6.0f} m   dune-line shift 1984->1997  =>  {o['new']:.1f} m",
              f"  3 cut       {o['model']:6.0f} m   int({o['new']:.1f} / 10) = row {o['row']}: the road rows",
              f"  v2 ran with {o['now']:.0f} m" + (" (step 1 floored at 0)" if o["v2"] < 0 else "")
              + f"; v3 runs with step 3.  Table raw: {o['raw']:.0f} m (1984 dune line -> mask)"]
    elif o["row"] is not None:
        L.append(f"  modeled  {o['model']:6.0f} m   row {o['row']} from interior row 0 (setback {o['new']:.1f} m)")
    else:
        L.append("  no model road in this domain (GIS 1-8)")
    where = {"road": "behind the road as placed" if n > 0 else "in front of today's pavement",
             "crest": "behind the crest row"}.get(str(t.get("insert_anchor", "")), "-")
    L += ["",
          f"v3: {t.get('rows_behind_road', '') or 'no rows change'}" + (f"   ({where})" if n else ""),
          f"flags: {t['flags'] if isinstance(t['flags'], str) else '-'}",
          "photo cover: " + ", ".join(f"{y} {c:.0%}" for y, c in cover.items()),
          "",
          "photos: white tick = row 0 from the 1984 road, (a) as placed",
          "for 1984, (b) the DEM's; bracket = the dune-line shift.",
          "model panel: dark rows the road at the 1984 setback,",
          "white arrow the modeled offset, red outline rows inserted,",
          "blue dashes the seam of rows removed."]
    return L


# =============================================================================
# THE MODEL DOMAIN
# =============================================================================

_GRID: dict = {}


def grid_dirs() -> tuple[Path, Path, str]:
    """dune-topo/v3 if it is on disk, else the resolved CURRENT version."""
    if "dirs" not in _GRID:
        try:
            _GRID["dirs"] = topo_dirs(R.PRODUCT, override=BUILT)
        except SystemExit:
            _GRID["dirs"] = topo_dirs(R.PRODUCT)
            print(f"  dune-topo/{BUILT} not on disk; the model panel draws {_GRID['dirs'][2]}")
    return _GRID["dirs"]


def load_grid(d: int) -> tuple[np.ndarray, np.ndarray, str]:
    """(topography m, dune m, version name) for one domain, cached."""
    if d not in _GRID:
        topo, dune, name = grid_dirs()
        z = np.load(topo / array_name("topography", d)) * R.CELL_M
        dn = np.load(dune / array_name("dune", d)) * R.CELL_M + BERM_EL_M
        _GRID[d] = (z, dn, name)
    return _GRID[d]


def model_r_max(t: pd.Series, o: dict, pr: pd.DataFrame, win, n_rows: int) -> int:
    """The landward-most interior row the model panel shows: at least what the
    photographs show, and past the road and the footprint's rows."""
    n = int(t["n_cells"])
    ins = _num(t.get("insert_row_behind_road"))
    ins = int(ins) if (n != 0 and np.isfinite(ins)) else 0
    r_win = int(np.ceil((float(pr["interior_x"].max()) - win[0]) / R.CELL_M)) + 1
    return min(n_rows - ROAD_ROWS - 1, max(r_win, (o["row"] or 0) + abs(n) + 12, ins + abs(n) + 8))


def draw_model(ax, d: int, t: pd.Series, o: dict, pr: pd.DataFrame, win, gd) -> None:
    """The PROCESSED domain, the model input as the model holds it: the
    straightened grid, dune rows a straight band, cross-shore rows across with
    the ocean on the right and alongshore cells up the page (south at the
    bottom, as the photographs). Nothing is mapped back through the shear."""
    z, dune, name = load_grid(d)
    n = int(t["n_cells"])
    cmap, norm, _bounds = elevation_cmap()
    img = np.concatenate([np.tile(dune[None, :], (ROAD_ROWS, 1)), z], axis=0)   # img row k = interior row k - 2
    C = R.CELL_M
    r_new = o["row"]
    ins = int(t["insert_row_behind_road"]) if n != 0 and np.isfinite(_num(t.get("insert_row_behind_road"))) else None
    ncol = img.shape[1]
    r_max = max(model_r_max(t, o, pr, win, img.shape[0]), ncol - ROAD_ROWS - 1)   # at least a square of cells
    # x = interior row (dune rows at -2, -1), y = alongshore cell; ocean right
    grid = img[:r_max + ROAD_ROWS + 1].T                                       # (alongshore, rows)
    ax.imshow(grid, cmap=cmap, norm=norm, aspect="equal", interpolation="nearest", origin="lower",
              extent=[-ROAD_ROWS - 0.5, r_max + 0.5, -0.5, ncol - 0.5], zorder=1)
    ax.axvline(-0.5, color="white", lw=0.8, ls=(0, (2, 3)), alpha=0.7, zorder=6)   # seaward edge of row 0
    halo = [pe.withStroke(linewidth=3.0, foreground="black")]
    lab = dict(fontsize=8, fontweight="bold", zorder=9, path_effects=halo)
    ax.text(-ROAD_ROWS + 0.5, ncol - 1.5, "dune", ha="center", va="top", color="white", rotation=90, **lab)
    y_ar = ncol * 0.72                               # the arrow high, the road label at the foot
    if r_new is not None:
        ax.axvspan(r_new - 0.5, r_new + 1.5, facecolor=R.C_ROAD, edgecolor="none", alpha=0.5, zorder=5)
        # the label hangs off the ocean edge when the road is close to it, so it stays inside the panel
        x_edge = -ROAD_ROWS - 0.3
        near = r_new < 12
        if r_new > 0:
            ax.annotate("", xy=(r_new - 0.5, y_ar), xytext=(-0.5, y_ar), zorder=12,
                        arrowprops=dict(arrowstyle="<->", color="white", lw=2.0, shrinkA=0, shrinkB=0,
                                        path_effects=halo))
        ax.text(x_edge if near else (r_new - 1) / 2, y_ar + 1.2, f"{o['model']:.0f} m ({r_new} rows)",
                ha="right" if near else "center", va="bottom", color="white", **lab)
    if n > 0 and ins is not None:
        ax.add_patch(plt.Rectangle((ins - 0.5, -0.5), n, ncol, facecolor="none", edgecolor=R.C_ADD, lw=2.0,
                                   zorder=7))                                  # named in the legend, not here
    elif n < 0 and ins is not None:
        ax.axvline(ins - 0.5, color=R.C_REM, lw=2.4, ls=(0, (3, 1.5)), zorder=7, path_effects=halo)
    # a SQUARE data window (the four panel boxes are equal squares): the rows
    # shown set the side; the alongshore range is centred on the domain and
    # padded with blank where the side exceeds 50 cells
    side = r_max + ROAD_ROWS + 1
    cy = (ncol - 1) / 2
    ax.set_xlim(r_max + 0.5, -ROAD_ROWS - 0.5)                                  # ocean right
    ax.set_ylim(cy - side / 2, cy + side / 2)
    ax.set_aspect("equal")
    ax.set_xlabel("cross-shore row from interior row 0 (dune rows seaward of it), and m", fontsize=8)
    ax.set_ylabel("alongshore cell (south at bottom), and m", fontsize=8)
    xt = [r for r in range(0, r_max + 1, 10)]
    ax.set_xticks(xt)
    ax.set_xticklabels([f"{r}\n{r * C:.0f} m" for r in xt])
    yt = [c for c in range(0, ncol, 10)]
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{c}  {c * C:.0f} m" for c in yt])
    ax.tick_params(labelsize=7)
    ax.text(0.99, 0.5, "ocean", transform=ax.transAxes, ha="right", va="center", fontsize=9, color="white",
            rotation=90, zorder=9, path_effects=halo)
    ax.set_title(f"(c) Model input: the Barrier3D grid, dune-topo/{name}"
                 + ("" if name == BUILT else f"  ({BUILT} not on disk)"), loc="left", pad=TITLE_PAD, color=R.INK,
                 fontweight="bold", fontsize=TITLE_SIZE)


def _long_date(s: str) -> str:
    """'1984-09-19' -> '19 September 1984'; a bare year stays a year."""
    try:
        return datetime.strptime(s, "%Y-%m-%d").strftime("%d %B %Y").lstrip("0")
    except ValueError:
        return s


def elevation_handles() -> list:
    cmap, _norm, bounds = elevation_cmap()
    labels = ["below 0 (water)"] + [f"{lo:g}-{hi:g}" for lo, hi in zip(bounds[1:-2], bounds[2:-1])] \
        + [f"above {bounds[-2]:g}"]
    return [Patch(facecolor=cmap(i), edgecolor="0.4", lw=0.4, label=f"{lab} m") for i, lab in enumerate(labels)]


# =============================================================================
# THE WINDOW
# =============================================================================

class App:
    def __init__(self, data: Data, ids: list[int], smoke: bool = False):
        self.data, self.ids, self.i, self.smoke = data, ids, 0, smoke
        self.blink_year: int | None = None
        self.blink_images: dict[int, object] = {}
        self.axes: list = []
        self.sections: dict[int, str] = {}                # island section per domain, for the header
        try:
            sec = pd.read_csv(R.ROAD_DIR / "RoadOffset_1984_domains.csv", usecols=["domain", "section"])
            self.sections = {int(r.domain): str(r.section) for r in sec.itertuples()}
        except Exception as e:
            print(f"  sections not read ({e}); the header goes without them")

        self.root = tk.Tk()
        self.blink = tk.BooleanVar(value=False)
        self.show_dem = tk.BooleanVar(value=False)        # the 1 m DEM as a fourth panel, on request
        self.root.title("HAT quick imagery review - the 1984 road offset and rows")
        self.root.geometry("1700x1000")
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)

        left = ttk.Frame(self.root)
        left.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        right = ttk.Frame(self.root, width=400, padding=8)
        right.pack(side=tk.RIGHT, fill=tk.Y)
        right.pack_propagate(False)

        off.apply_style()
        self.fig = Figure(figsize=(12.5, 9.0), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=left)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, left, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        self.canvas.mpl_connect("scroll_event", self.on_scroll)

        self.build_right(right)
        self.root.bind("<Right>", lambda e: self._key(lambda: self.step(+1)))
        self.root.bind("<Left>", lambda e: self._key(lambda: self.step(-1)))
        self.root.bind("<space>", lambda e: self._key(self.flip))
        self.root.bind("<b>", lambda e: self._key(self.flip))
        self.root.bind("<Control-s>", lambda e: self.save())
        self.status("ready")
        self.show(0)

    # ---- the right-hand panel ------------------------------------------------
    def build_right(self, f):
        nav = ttk.Frame(f)
        nav.pack(fill=tk.X)
        ttk.Button(nav, text="< Prev", command=lambda: self.step(-1)).pack(side=tk.LEFT)
        self.combo = ttk.Combobox(nav, state="readonly", width=30)
        self.combo.pack(side=tk.LEFT, padx=4, fill=tk.X, expand=True)
        self.combo.bind("<<ComboboxSelected>>", lambda e: self.show(self.combo.current()))
        ttk.Button(nav, text="Next >", command=lambda: self.step(+1)).pack(side=tk.LEFT)
        self.refresh_combo()

        self.info = tk.Text(f, height=24, width=50, wrap="word", relief="flat",
                            font=("Consolas", 9), background="#f4f4f4")
        self.info.pack(fill=tk.X, pady=(8, 4))
        self.info.configure(state="disabled")

        ttk.Checkbutton(f, text="blink mode: one panel, Space flips the year",
                        variable=self.blink, command=lambda: self.show(self.i)).pack(anchor="w", pady=(4, 0))
        ttk.Checkbutton(f, text="DEM panel: the 1 m surface row 0 was picked on, with the search window",
                        variable=self.show_dem, command=lambda: self.show(self.i)).pack(anchor="w")

        ttk.Label(f, text="Verdict", font=("Segoe UI", 10, "bold")).pack(anchor="w", pady=(10, 0))
        self.fields: dict[str, ttk.Combobox] = {}
        for c in R.QUICK_COLS:
            row = ttk.Frame(f)
            row.pack(fill=tk.X, pady=2)
            ttk.Label(row, text=QUESTION[c], width=22).pack(side=tk.LEFT)
            cb = ttk.Combobox(row, values=VOCAB, state="readonly", width=12)
            cb.pack(side=tk.LEFT)
            self.fields[c] = cb
        ttk.Label(f, text="notes").pack(anchor="w", pady=(6, 0))
        self.notes = tk.Text(f, height=4, width=50, wrap="word", font=("Segoe UI", 9))
        self.notes.pack(fill=tk.X)

        btns = ttk.Frame(f)
        btns.pack(fill=tk.X, pady=(10, 0))
        ttk.Button(btns, text="Save  (Ctrl+S)", command=self.save).pack(side=tk.LEFT)
        ttk.Button(btns, text="Save & next", command=lambda: (self.save(), self.step(+1))).pack(side=tk.LEFT, padx=6)
        ttk.Button(btns, text="Summarize", command=self.summarize).pack(side=tk.RIGHT)
        ttk.Button(btns, text="Save figure", command=self.export).pack(side=tk.RIGHT, padx=6)
        self.status_var = tk.StringVar()
        ttk.Label(f, textvariable=self.status_var, foreground="#444", wraplength=380).pack(anchor="w", pady=(8, 0))
        ttk.Label(f, text=f"sheet: {R.SHEET.relative_to(R.SCOPE_DIR.parent)}", foreground="#888",
                  wraplength=380).pack(anchor="w", side=tk.BOTTOM)

    def sheet(self) -> pd.DataFrame:
        if not R.SHEET.is_file():
            return pd.DataFrame()
        df = pd.read_csv(R.SHEET, dtype=str).fillna("")
        df["domain"] = df["domain"].astype(int)
        return df.set_index("domain")

    def refresh_combo(self):
        sh = self.sheet()
        labels = []
        for d in self.ids:
            n = int(self.data.tab.loc[d, "n_cells"])
            role = sh.loc[d, "role"] if d in sh.index and "role" in sh else ("control" if n == 0 else "changed")
            done = d in sh.index and any(_s(sh.loc[d, c]) for c in R.QUICK_COLS if c in sh)
            labels.append(f"GIS {d:>2}   {n:+d} rows   {role:8s} {'done' if done else ''}")
        self.combo["values"] = labels
        if labels:
            self.combo.current(min(self.i, len(labels) - 1))

    def status(self, msg: str):
        self.status_var.set(f"{datetime.now():%H:%M:%S}  {msg}")
        self.root.update_idletasks()

    def _key(self, fn):
        if isinstance(self.root.focus_get(), (tk.Text, ttk.Combobox)):
            return
        fn()

    # ---- navigation ---------------------------------------------------------
    def step(self, k: int):
        j = self.i + k
        if 0 <= j < len(self.ids):
            self.show(j)

    def show(self, j: int):
        self.i = j
        d = self.ids[j]
        self.combo.current(j)
        self.status(f"loading GIS {d} ...")
        try:
            rec = self.data.load(d)
        except Exception as e:
            self.status(f"GIS {d}: {e}")
            return
        self.fill_form(d)
        self.draw(d, rec)
        self.fill_info(d, rec)
        self.status(f"GIS {d} shown")
        if j + 1 < len(self.ids):
            threading.Thread(target=self._prefetch, args=(self.ids[j + 1],), daemon=True).start()
        if self.smoke:
            self.root.after(1500, self._smoke_done)

    def _prefetch(self, d: int):
        try:
            self.data.load(d)
        except Exception:
            pass

    # ---- drawing ------------------------------------------------------------
    def draw(self, d: int, rec: dict):
        self.fig.clear()
        self.blink_images = {}
        pr, t = rec["pr"], rec["t"]
        n = int(t["n_cells"])
        o = offsets(t)
        x_lo, x_hi, y_lo, y_hi = rec["win"]
        years = list(self.data.years)
        blink = self.blink.get()
        n_photo = 1 if blink else len(years)
        gd = self.data.gdf[self.data.gdf["domain_id"].astype(int) == d]
        box = gd.geometry.iloc[0]
        drawn = off.clip_for_drawing(self.data.lines, box.buffer(30.0))
        try:
            meas = line_offsets(pr, self.data.lines, self.data.roads, box)
        except Exception as e:
            print(f"  GIS {d}: line offsets not measured ({e})")
            meas = {}
        self._meas = meas
        # photographs down the left (1984 over 1997), the model input upper right,
        # the legend lower right; in blink mode one photograph left, the model right.
        # Four equal cells; every panel is a square data window with equal aspect,
        # so the four boxes come out the same size and line up (Hannah, 2026-09-10).
        # The header has its own row above; titles sit at one pad; the strip
        # below holds the legend when the DEM panel takes the legend's cell.
        dem = self.show_dem.get()
        dem_slot = None
        gkw = dict(wspace=0.16, hspace=0.24, left=0.03, right=0.99)
        if n_photo == 1 and not dem:
            gs = self.fig.add_gridspec(1, 2, top=0.90, bottom=0.10, **gkw)
            photo_slots, model_slot, leg_slot = [gs[0, 0]], gs[0, 1], None
        elif n_photo == 1:                                   # blink + DEM: photo | model over DEM | legend
            gs = self.fig.add_gridspec(2, 2, top=0.92, bottom=0.05, **gkw)
            photo_slots, model_slot, dem_slot, leg_slot = [gs[0, 0]], gs[0, 1], gs[1, 0], gs[1, 1]
        elif dem:                                            # photos | model over DEM | a narrow legend column
            gs = self.fig.add_gridspec(2, 3, width_ratios=[1.0, 1.0, 0.56], top=0.92, bottom=0.05,
                                       wspace=gkw["wspace"], hspace=gkw["hspace"], left=gkw["left"], right=0.995)
            photo_slots, model_slot, dem_slot, leg_slot = [gs[0, 0], gs[1, 0]], gs[0, 1], gs[1, 1], gs[:, 2]
        else:
            gs = self.fig.add_gridspec(2, 2, top=0.92, bottom=0.05, **gkw)
            photo_slots, model_slot, leg_slot = [gs[0, 0], gs[1, 0]], gs[0, 1], gs[1, 1]
        axes = []
        for s in photo_slots:
            axes.append(self.fig.add_subplot(s, sharex=axes[0] if axes else None, sharey=axes[0] if axes else None))
        ag = self.fig.add_subplot(model_slot)
        try:
            draw_model(ag, d, t, o, pr, rec["win"], gd)
        except Exception as e:                       # a missing array must not kill the review
            ag.text(0.5, 0.5, f"model panel unavailable: {e}", transform=ag.transAxes, ha="center", va="center",
                    color=off.INK_MUTED, fontsize=9, wrap=True)
            ag.set_xticks([])
            ag.set_yticks([])

        shift = _num(t.get("shift_m_median"))               # + = the 1984 line seaward of the 1997 line

        roads_ref = {ROAD_REF: self.data.roads[ROAD_REF]} if ROAD_REF in self.data.roads else {}

        def overlays(ax, y=None):
            # the settled set (interview 2026-09-10): the 1984 road, both dune
            # lines, the domain box, one crest tick per year, one shift bracket
            gd.boundary.plot(ax=ax, color="0.3", linewidth=0.8, zorder=4)
            off.draw_lines(ax, drawn, scale=1.0, style=off.LINE_STYLE)
            off.m.draw_roads(ax, roads_ref, scale=0.9)
            draw_shift_bracket(ax, pr, meas, shift, n)
            for k, yy in enumerate([y] if y is not None else years):   # one per panel; blink mode stacks both
                draw_crest_tick(ax, yy, meas, o, slot=k)
            ax.set_xlim(x_lo, x_hi)
            ax.set_ylim(y_lo, y_hi)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.text(0.99, 0.5, "ocean", transform=ax.transAxes, ha="right", va="center", fontsize=9,
                    color="white", rotation=90, zorder=9)

        def photo(ax, y, visible=True):
            shown = rec["photos"].get(y)
            if shown is not None:
                return ax.imshow(shown, extent=(x_lo, x_hi, y_lo, y_hi), origin="upper", zorder=1,
                                 interpolation="bilinear", visible=visible)
            ax.set_facecolor("0.93")
            return ax.text(0.5, 0.5, f"no {y} photograph here", transform=ax.transAxes,
                           ha="center", va="center", visible=visible)

        if blink:
            ax = axes[0]
            for y in years:
                self.blink_images[y] = photo(ax, y, visible=False)
            if self.blink_year not in years:
                self.blink_year = years[0]
            overlays(ax)
            off._scalebar(ax, length_m=100.0)
            off._north_arrow(ax, x=0.94, y=0.12)
        else:
            for k, y in enumerate(years):
                ax = axes[k]
                photo(ax, y)
                overlays(ax, y)
                ax.set_title(f"({chr(97 + k)}) {_long_date(self.data.imagery[y].date)}, aerial photograph",
                             loc="left", color=YEAR_COLOUR.get(y, R.INK), fontweight="bold", pad=TITLE_PAD,
                             fontsize=TITLE_SIZE)
                off._scalebar(ax, length_m=100.0)
                if k == 0:
                    off._north_arrow(ax, x=0.94, y=0.12)

        # ---- (d) the DEM, on request: the surface row 0 was picked on, with the
        # diagnostic layers that belong to it (row 0 line, search window)
        ad = None
        if dem_slot is not None:
            ad = self.fig.add_subplot(dem_slot, sharex=axes[0], sharey=axes[0])
            try:
                arr, extent = self.data.lidar(d)
                if arr is not None:
                    off._hillshade(ad, arr, extent, res=1.0)
            except Exception as e:
                ad.text(0.5, 0.5, f"DEM unavailable: {e}", transform=ad.transAxes, ha="center", va="center",
                        color=off.INK_MUTED, fontsize=9, wrap=True)
            draw_search_window(ad, pr, search_window(self.data.frames, d))
            prs_ = pr.sort_values("interior_y")
            ad.plot(prs_["interior_x"], prs_["interior_y"], color="white", lw=1.0, ls=(0, (2, 3)), zorder=6)
            overlays(ad, 1997)                             # the DEM's own tick is the 1997 one (row 0 of the DEM)
            ad.set_title("(d) 1 m DEM, hillshade: the 1996 foredune on the 2009 backdune", loc="left",
                         color=R.INK, fontweight="bold", fontsize=TITLE_SIZE, pad=TITLE_PAD)

        # legends: one for the photographs, one for the model input
        h_photo = [Line2D([0], [0], **dict(off.LINE_STYLE[1984], linewidth=2.0), label="1984 dune line (digitized)"),
                   Line2D([0], [0], **dict(off.LINE_STYLE[1997], linewidth=2.0), label="1997 dune line (digitized)")]
        h_photo += off.m.road_legend_handles(roads_ref)
        h_photo += [Line2D([0], [0], color="0.3", lw=0.8, label="domain box"),
                    Line2D([0], [0], color="0.45", lw=2.4, marker="|", ms=9, mew=2.4, ls=(0, (1, 2)),
                           label="row 0 tick (a model, b DEM)"),
                    Line2D([0], [0], color="0.45", lw=1.4, marker="|", ms=7, mew=1.4,
                           label="dune-line shift = rows N")]
        h_model = [Line2D([0], [0], color="0.5", lw=2.0, marker="|", ms=7, mew=2.0,
                          label="modeled offset: row 0 → road"),
                   Patch(facecolor=R.C_ROAD, alpha=0.5, label="NC-12 rows, 1984 setback")]
        if n > 0:
            h_model.append(Patch(facecolor="none", edgecolor=R.C_ADD, lw=2.0,
                                 label=f"{n} rows inserted (copied)"))
        elif n < 0:
            h_model.append(Line2D([0], [0], color=R.C_REM, lw=2.4, ls=(0, (3, 1.5)),
                                  label=f"seam: {-n} rows removed"))
        h_model += elevation_handles()
        if ad is not None:
            h_photo += [Line2D([0], [0], color="0.4", lw=1.0, ls=(0, (2, 3)), label="(d) interior row 0, per profile"),
                        Patch(facecolor=C_WINDOW, alpha=0.45, edgecolor=C_WINDOW, lw=1.0, ls=(0, (2, 2)),
                              label="(d) dune search window")]
        if leg_slot is not None:
            al = self.fig.add_subplot(leg_slot)
            al.axis("off")
            # side by side in the legend cell; stacked when the legend is the narrow column beside the DEM
            stacked = ad is not None and n_photo == 2
            l1 = al.legend(handles=h_photo, loc="upper left", bbox_to_anchor=(0.0, 1.0), fontsize=9.5,
                           frameon=False, title="Photographs (a, b)" + (", DEM (d)" if ad is not None else ""),
                           title_fontsize=10.5, alignment="left", handlelength=2.0, labelspacing=0.4,
                           borderaxespad=0.0)
            al.add_artist(l1)
            al.legend(handles=h_model, loc="upper left", bbox_to_anchor=(0.0, 0.52) if stacked else (0.55, 1.0),
                      fontsize=9.5, frameon=False, title="Model input (c), m MHW", title_fontsize=10.5,
                      alignment="left", handlelength=2.0, labelspacing=0.4, borderaxespad=0.0)
        else:
            self.fig.legend(handles=h_photo + h_model[:4], loc="lower center", ncol=5, fontsize=9, frameon=False)
        # the header and the source line
        o_txt = (f"1984 setback {o['new']:.1f} m → row {o['row']}" if o["row"] is not None
                 else "no model road")
        self.fig.suptitle(f"GIS {d}  ·  {self.sections.get(d, '')}  ·  N = {n:+d} rows  ·  {o_txt}",
                          x=0.01, y=0.99, ha="left", va="top", fontsize=10.5, fontweight="bold", color=R.INK)
        self.fig.text(0.02, 0.006,
                      "Photographs: USGS Henderson release (doi 10.5066/P1CXBCDW), georeferenced to the 2007 orthophotos, "
                      "stated horizontal accuracy 1.2 m. NC-12 is the 1984 line (the 1978 export); the tick is interior "
                      "row 0 measured seaward from the pavement edge on the middle profile, (a) as the model places it for "
                      "1984, (b) as the DEM has it; the bracket is the median dune-line shift over the 50 profiles. "
                      "Model: Barrier3D, 10 m cells, elevations m MHW.",
                      fontsize=7, color=off.INK_MUTED, ha="left", va="bottom", wrap=True)
        self.axes = axes + [ag] + ([ad] if ad is not None else [])
        if blink:
            self._show_blink_year()
        self.canvas.draw_idle()

    def on_scroll(self, event):
        ax = event.inaxes
        if ax is None or ax not in self.axes or event.xdata is None:
            return
        f = 0.8 if event.button == "up" else 1.25
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        ax.set_xlim(event.xdata - (event.xdata - x0) * f, event.xdata + (x1 - event.xdata) * f)
        ax.set_ylim(event.ydata - (event.ydata - y0) * f, event.ydata + (y1 - event.ydata) * f)
        self.canvas.draw_idle()

    # ---- blink ----------------------------------------------------------------
    def _show_blink_year(self):
        for y, im in self.blink_images.items():
            im.set_visible(y == self.blink_year)
        if self.blink.get() and self.axes and self.blink_year is not None:
            rec = self.data.load(self.ids[self.i])
            self.axes[0].set_title(f"{self.data.imagery[self.blink_year].date}   cover "
                                   f"{rec['cover'].get(self.blink_year, 0):.0%}    (Space flips)",
                                   loc="left", color=YEAR_COLOUR.get(self.blink_year, "0.2"))
        self.canvas.draw_idle()

    def flip(self):
        if not self.blink.get() or not self.blink_images:
            return
        years = list(self.blink_images)
        self.blink_year = years[(years.index(self.blink_year) + 1) % len(years)]
        self._show_blink_year()

    # ---- the form -------------------------------------------------------------
    def fill_info(self, d: int, rec: dict):
        self.info.configure(state="normal")
        self.info.delete("1.0", tk.END)
        self.info.insert("1.0", "\n".join(offset_lines(d, rec["t"], offsets(rec["t"]), rec["cover"],
                                                        getattr(self, "_meas", {}))))
        self.info.configure(state="disabled")

    def fill_form(self, d: int):
        sh = self.sheet()
        for c, cb in self.fields.items():
            cb.set(_s(sh.loc[d, c]) if d in sh.index and c in sh else "")
        self.notes.delete("1.0", tk.END)
        self.notes.insert("1.0", _s(sh.loc[d, "notes"]) if d in sh.index and "notes" in sh else "")

    def save(self):
        d = self.ids[self.i]
        if not R.SHEET.is_file():
            self.status("no sheet on disk - run HAT_imagery_review_1984.py first")
            return
        df = pd.read_csv(R.SHEET, dtype=str).fillna("")
        df["domain"] = df["domain"].astype(int)
        df = df.set_index("domain")
        if d not in df.index:
            self.status(f"GIS {d} is not in the sheet - run the batch script on it first")
            return
        for c in R.QUICK_COLS + ["notes", "reviewed_by", "reviewed_at"]:
            if c not in df:
                df[c] = ""
        for c, cb in self.fields.items():
            df.loc[d, c] = cb.get()
        df.loc[d, "notes"] = self.notes.get("1.0", tk.END).strip().replace("\n", " / ")
        df.loc[d, "reviewed_by"] = getpass.getuser()
        df.loc[d, "reviewed_at"] = f"{datetime.now():%Y-%m-%d %H:%M}"
        df.to_csv(R.SHEET)
        self.refresh_combo()
        self.status(f"GIS {d} saved to {R.SHEET.name}")

    def export(self):
        """The view on screen as a figure, into the figures tree by the sign of N."""
        d = self.ids[self.i]
        try:
            out = (R.insert_figures_dir_for_domain(R.PRODUCT, "3-placement", d, under="imagery-review")
                   / f"HAT_imagery_review_quick_GIS{d}.png")
            self.fig.savefig(out, dpi=200, facecolor="white")
            R.upsert_caption(
                "## `HAT_imagery_review_quick_GIS<N>.png` (3-placement/imagery-review/rows-added, rows-removed, unchanged)",
                "The quick imagery review of one domain, saved from the review window. (a) The 19 September 1984 "
                "and (b) the 12 October 1997 aerial photographs of the same window (USGS Henderson release, doi "
                "10.5066/P1CXBCDW, georeferenced to the 2007 orthophotos, stated accuracy 1.2 m, 0.5 m drawing "
                "resolution with a 2–98 % stretch), ocean to the right, with the 1984 (red) and 1997 (blue) dune "
                "lines as digitized, NC-12 in 1984 (dashed, the 1978 export) and 2004 (solid), interior row 0 "
                "NC-12 as digitized in 1984 (the 1978 export), and two marks on the middle profile: a white tick at "
                "interior row 0 measured seaward from the pavement edge, on (a) as the model places it for 1984 "
                "(the 1984 setback cut to cells, the number panel (c) shows) and on (b) as the DEM has it, so the "
                "crest visible in each photograph can be judged against it; and a bracket between the two dune "
                "lines labelled with the median shift over the 50 profiles and the rows it became. "
                "(c) The model input for "
                "the same domain: the straightened Barrier3D grid of dune-topo/v3, cross-shore rows across with "
                "the ocean on the right and alongshore cells up the page (south at the bottom), the two dune rows "
                "drawn at berm + dune height, elevations in m MHW; on it interior row 0 (white dashes), the road "
                "rows at the 1984 setback (dark), "
                "the rows the footprint inserted (red outline) or the seam where it removed rows (blue dashes), "
                "and the modeled offset from row 0 to the road as an arrow in rows. Header: domain, island "
                "section, N and the 1984 setback with the row it cuts to. Panels (a) and (b) share one view; (c) "
                "is in cells. The verdicts for the domain are in `imagery_review_1984.csv`.")
            self.status(f"figure saved: {out.relative_to(R.SCOPE_DIR)}")
        except Exception as e:
            self.status(f"figure not saved: {e}")

    def summarize(self):
        try:
            fig, rep = summarize()
            self.status(f"summary written: {fig.name}, {rep.name}")
        except Exception as e:
            self.status(f"summary failed: {e}")

    # ---- lifecycle ------------------------------------------------------------
    def _smoke_done(self):
        out = CACHE_DIR.parent / "imagery_review_quick_smoke.png"       # outside the repo
        self.fig.savefig(out, dpi=100)
        self.export()                                                    # exercises the figure export too
        print(f"smoke: drew GIS {self.ids[self.i]}, saved {out}; {self.status_var.get()}")
        self.on_close()

    def on_close(self):
        self.root.quit()
        self.root.destroy()

    def run(self):
        self.root.mainloop()


# =============================================================================
# THE TALLY
# =============================================================================

def summarize(sheet: Path | None = None) -> tuple[Path, Path]:
    sheet = sheet or R.SHEET
    if not sheet.is_file():
        raise SystemExit(f"{sheet} not found - run HAT_imagery_review_1984.py first")
    df = pd.read_csv(sheet, dtype=str).fillna("")
    df["domain"] = df["domain"].astype(int)
    df = df.set_index("domain").sort_index()
    for c in R.QUICK_COLS + ["notes"]:
        if c not in df:
            df[c] = ""
    df["n"] = df["n_cells"].map(_num).astype(int)
    tab = pd.read_csv(R.FOOTPRINT_CSV).set_index("domain")
    o = pd.DataFrame({d: offsets(tab.loc[d]) for d in df.index}).T
    df["raw"] = o["raw"].astype(float)
    df["model"] = o["model"].astype(float)
    ch, ct = df[df["n"] != 0], df[df["n"] == 0]

    off.apply_style()
    fig, (ax, ab) = plt.subplots(2, 1, figsize=(14.0, 8.0), constrained_layout=True, sharex=True)
    for d, r in ch.iterrows():
        col = VERDICT_C.get(_s(r["rows_ok"]))
        ax.bar(d, r["n"], width=0.8, facecolor=col or "white", edgecolor=col or "0.5", lw=1.0, zorder=3)
    for d, r in ct.iterrows():
        col = VERDICT_C.get(_s(r["rows_ok"]))
        ax.plot(d, 0, marker="o", ms=6, mfc=col or "white", mec="0.3", mew=1.0, zorder=4, ls="none")
    ax.axhline(0, color=off.INK, lw=0.8, zorder=2)
    ax.set_ylabel("rows in the footprint, N (signed)")
    ax.grid(True, axis="y", alpha=0.4)
    n_j = int(ch["rows_ok"].map(_s).astype(bool).sum())
    off._title(ax, 0, f"rows N by the reviewer: {n_j} of {len(ch)} changed domains judged, "
                      f"{int(ct['rows_ok'].map(_s).astype(bool).sum())} of {len(ct)} controls")
    handles = [Patch(facecolor=c, edgecolor=c, label=f"rows_ok = {k}") for k, c in VERDICT_C.items()]
    handles.append(Patch(facecolor="white", edgecolor="0.5", label="not yet judged"))
    handles.append(Line2D([0], [0], marker="o", mfc="white", mec="0.3", ls="none", label="control (N = 0)"))
    ax.legend(handles=handles, loc="upper left", fontsize=8, ncol=5, frameon=False)

    rd = df[np.isfinite(df["model"])]
    ab.plot(rd.index, rd["raw"], marker="o", ms=5, color="#c9a800", ls="none", label="raw: 1984 dune line to NC-12, as digitized")
    ab.plot(rd.index, rd["model"], marker="s", ms=5, color="0.25", ls="none", label="modeled: row cut of the 1984 setback")
    for d, r in rd.iterrows():
        ab.plot([d, d], [r["raw"], r["model"]], color="0.7", lw=0.8, zorder=1)
        v = _s(r["offset_ok"])
        if v in ("no", "unclear"):
            ab.plot(d, max(r["raw"], r["model"]) + 15, marker="x" if v == "no" else "?", color=VERDICT_C[v],
                    ms=7, mew=1.8, ls="none", zorder=5)
    ab.set_xlim(0, 91)
    ab.set_xlabel("GIS domain (south at left)")
    ab.set_ylabel("1984 road offset, m")
    ab.set_xticks(range(5, 91, 5))
    ab.grid(True, axis="y", alpha=0.4)
    n_o = int(rd["offset_ok"].map(_s).astype(bool).sum())
    off._title(ab, 1, f"the 1984 road offset, raw and modeled: {n_o} of {len(rd)} road domains judged; "
                      f"x = offset_ok 'no', ? = 'unclear'")
    ab.legend(loc="upper left", fontsize=8, frameon=False)
    out = insert_figures_dir(R.PRODUCT, "3-placement", "imagery-review/island") / SUMMARY_PNG
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    def counts(sub, col):
        v = sub[col].map(_s)
        return ", ".join(f"{k} {int((v == k).sum())}" for k in ("yes", "no", "unclear")) + f", blank {int((v == '').sum())}"

    L = [f"HAT_imagery_review_quick.txt - the quick review, tallied ({datetime.now():%Y-%m-%d %H:%M})",
         f"sheet: {sheet}", "",
         "offset_ok: does the modeled 1984 road offset look right against the photographs?",
         "rows_ok:   does N, the rows the footprint adds or removes, look right?", ""]
    for name, sub in (("rows added", ch[ch["n"] > 0]), ("rows removed", ch[ch["n"] < 0]), ("controls", ct)):
        L.append(f"{name.upper()} ({len(sub)})")
        for c in R.QUICK_COLS:
            L.append(f"    {c:10s} {counts(sub, c)}")
    L.append("")
    flagged = df[df["offset_ok"].map(_s).isin(["no", "unclear"]) | df["rows_ok"].map(_s).isin(["no", "unclear"])]
    L.append(f"DOMAINS ANSWERED 'no' OR 'unclear' ({len(flagged)})")
    L.append("  domain    N   raw    model   offset_ok  rows_ok   notes")
    for d, r in flagged.iterrows():
        L.append(f"  {d:>6}  {r['n']:+3d}  {r['raw']:5.0f}  {r['model']:5.0f}   {_s(r['offset_ok']) or '-':9s}  "
                 f"{_s(r['rows_ok']) or '-':8s}  {_s(r['notes'])}")
    L.append("")
    L.append(f"raw minus modeled offset over the {len(rd)} road domains: median "
             f"{np.median(rd['raw'] - rd['model']):+.0f} m (the toe-to-row-0 gap plus the cell cut)")
    L.append(f"figure: {out}")
    SUMMARY_TXT.write_text("\n".join(L) + "\n", encoding="utf-8")
    R.upsert_caption(
        f"## `{SUMMARY_PNG}` (3-placement/imagery-review/island)",
        "The quick imagery review, tallied. (a) Every reviewed domain along the island, south at left: "
        "the signed rows of the 1984 footprint as a bar, coloured by the reviewer's answer to whether N "
        "looks right against the 1984 and 1997 photographs (dark yes, red no, orange unclear; hollow not "
        "yet judged); the controls, unchanged neighbours of the changed runs, as circles at zero. (b) The "
        "1984 road offset per road domain, raw (yellow circles: the 1984 dune line to the 1984 NC-12 "
        "line, both as digitized, median over the 50 profiles) and as modeled (dark squares: the 1984 "
        "setback in the row-0 convention cut to whole 10 m cells, the road's first row in v3), joined by "
        "a grey tie; a red cross marks a domain whose modeled offset the reviewer judged wrong, a "
        "question mark one judged unclear. Counts and the flagged domains in "
        f"`HAT_imagery_review_quick.txt`; the answers in `imagery_review_1984.csv`.")
    print(f"wrote {out}\nwrote {SUMMARY_TXT}")
    return out, SUMMARY_TXT


# =============================================================================

def main() -> None:
    ap = argparse.ArgumentParser(description="the quick imagery review: offsets raw and modeled, two questions")
    ap.add_argument("--domains", default="", help="comma-separated GIS ids (default: the sheet's domains)")
    ap.add_argument("--years", default=",".join(str(y) for y in R.DEFAULT_YEARS))
    ap.add_argument("--summary", action="store_true", help="tally the sheet and exit; no window")
    ap.add_argument("--smoke", action="store_true", help="open, draw the first domain, screenshot, close")
    a = ap.parse_args()
    if a.summary:
        summarize()
        return
    years = [int(y) for y in a.years.split(",") if y.strip()]
    data = Data(years)
    if a.domains:
        ids = [int(x) for x in a.domains.split(",")]
    elif R.SHEET.is_file():
        ids = sorted(int(d) for d in pd.read_csv(R.SHEET)["domain"])
    else:
        ids = sorted(int(d) for d in data.tab.index if int(data.tab.loc[d, "n_cells"]) != 0)
    App(data, ids, smoke=a.smoke).run()


if __name__ == "__main__":
    main()
