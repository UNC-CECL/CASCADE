#!/usr/bin/env python3
r"""
HAT_report_row_insert_scope.py
==============================================================================
Which domains would have interior rows ADDED behind the dune, which would have
rows REMOVED, and how many - drawn on the Barrier3D grid the way the model
would hold it, with a report and a per-domain table.

SCOPE ONLY. Nothing is written into a topography version, no array is
modified, no elevation is fabricated. The added rows are drawn BLANK.

REWRITTEN 2026-09-07 for the symmetric footprint. Until then this script
computed N itself as round(shift / 10) FLOORED AT ZERO and reported what the
negatives "would have implied". Hannah's decision that day made the footprint
symmetric (rows removed where the island prograded), changed the per-domain
statistic to the median of PAIRED per-profile differences, and changed the
rounding to a 10 m threshold (trunc). Those rules are applied in ONE place,
HAT_footprint_1984.py, which writes `footprint_1984_by_domain.csv`; this script
READS that table rather than re-deriving N, so the two cannot disagree. The
rules, and the assumptions behind them, are in that script's docstring.

WHAT THIS ADDS TO THE FOOTPRINT SCRIPT
    * the grid drawn AS THE MODEL WOULD HOLD IT: dune rows on top, then the
      added rows, then the existing interior pushed down the page by N; where
      rows are removed, the existing rows 0..|N|-1 are hatched between the
      dune and the rows that survive. HAT_footprint_1984_grid.png draws the
      same footprint in the CURRENT frame (row 0 fixed); this one shows the
      stack. NC-12 is drawn at its measured position in both.
    * the cross-check against the independent easting-frame measurement
      (`0-elevation/2009-2014-1996-duneline/duneline_offset_by_domain.csv`:
      raw easting in the axis-aligned box, 1 m sampling, no extractor frame,
      no row 0). Same two geojsons, so agreement bounds the FRAME, not the
      lines. The same trunc rule is applied to it.

THE PLAN VIEW moved. `HAT_row_insert_plan.png` (add-only) was retired the same
day; the plan view of the symmetric footprint is
`figures/2-footprint-1984/seaward/HAT_footprint_1984_plan.png`, drawn by HAT_footprint_1984.py.
Drawing it twice under two names would be one figure with two provenances.

OUTPUTS (1-barrier3d-domains/1984-start/row-insert-scope/)
    HAT_row_insert_scope.txt          the report
    row_insert_scope_by_domain.csv    per domain: signed N, shift, the easting
                                      cross-check, rows now/after
    figures/2-footprint-1984/seaward/HAT_row_insert_grid.png             the stacked grid, both signs
    figures/2-footprint-1984/behind-road/HAT_row_insert_grid_behindroad.png  the same rows behind NC-12 (--anchor road)
    figures/2-footprint-1984/HAT_row_insert_rows.png                     rows per domain, signed

USAGE
    python HAT_footprint_1984.py          # first - writes the footprint table
    python HAT_report_row_insert_scope.py
==============================================================================
"""

from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle


def _find_project_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit(f"cannot find data/hatteras_init above {start}")


REPO = _find_project_root(Path(__file__).resolve())
INIT = REPO / "data" / "hatteras_init"
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "0-elevation" / "3-figures"))
from hat_topo_version import (  # noqa: E402
    array_name, insert_figures_dir, topo_dirs)
import HAT_plot_duneline_offset as off  # noqa: E402  the house style
from hat_figure_style import elevation_cmap  # noqa: E402  the elevation classes

PRODUCT = "1984-start"
START_DIR = INIT / "1-barrier3d-domains" / PRODUCT
OUT_DIR = START_DIR / "row-insert-scope"
FIG_DIR = insert_figures_dir(PRODUCT, "2-footprint-1984")            # rows per domain
FIG_SEAWARD = insert_figures_dir(PRODUCT, "2-footprint-1984", "seaward")
FIG_ROAD = insert_figures_dir(PRODUCT, "2-footprint-1984", "behind-road")

# THE footprint: the one table every consumer of N reads (2026-09-07).
FOOTPRINT_CSV = OUT_DIR / "footprint_1984_by_domain.csv"
# The independent easting-frame measurement, for the cross-check column only.
EASTING_CSV = (INIT / "0-elevation" / "2009-2014-1996-duneline"
               / "duneline_offset_by_domain.csv")

CELL_M = 10.0
DUNE_ROWS = 2            # DuneWidth; row 1 is a copy of row 0
SENTINEL_DAM = -0.30     # the extractor's water sentinel, in decametres
ANCHOR = "dune"          # set from --anchor; "road" draws the behind-the-road placement
ROWS_SHOWN = 200         # every row: dune (2) + up to 7 added + the deepest interior (189); was 46 until 2026-09-07, when Hannah asked for the full domains
DOMAINS_PER_STRIP = 30

# Colours: the RdBu pair every 1984 dune-line figure uses. Red = 1984 line
# seaward = ground ADDED; blue = 1984 line landward = ground REMOVED.
C_ADD, C_ADD_FILL = off.C_1984, off.C_1984_FILL
C_REM = off.C_1997
C_DUNE, C_LAND, C_WATER, C_ROAD = "#c8a165", "#f0e6c8", "#a8c8e0", "#1a1a1a"
INK = off.INK


def n_cells(metres: float) -> int:
    """trunc(shift / 10): a row only once a FULL cell of change is measured.
    The rule HAT_footprint_1984.py applies; repeated here ONLY for the
    easting-frame cross-check, which that script does not carry."""
    return int(np.trunc(metres / CELL_M))


def read_footprint(path: Path) -> dict[int, dict]:
    if not path.exists():
        raise SystemExit(
            f"\nfootprint table not found:\n    {path}\n"
            f"  Write it first:  python HAT_footprint_1984.py\n"
            f"  (this report reads N from it rather than re-deriving it)\n")
    out = {}
    for r in csv.DictReader(open(path, encoding="utf-8")):
        out[int(r["domain"])] = r
    return out


def read_easting(path: Path) -> dict[int, float]:
    if not path.exists():
        print(f"  NOTE: {path.name} absent - cross-check column left blank")
        return {}
    out = {}
    for r in csv.DictReader(open(path)):
        v = r["offset_med_m"]
        if v not in ("", "nan"):
            out[int(r["domain"])] = float(v)
    return out


def _f(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


# =============================================================================
# THE FIGURE - the grid as the model would hold it
# =============================================================================

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


def build_strip(domains, topo_dir, n_by):
    """
    One alongshore strip as RGBA, rows DOWN the page from the dune. The dune
    stays put. Added rows go in between the dune and the existing interior,
    which is pushed down by N. Removed rows are the existing rows 0..|N|-1,
    left in place here and hatched by the caller: what the model would hold is
    the interior starting at row |N|. Existing cells carry the project's
    elevation classes; white is off the array.
    """
    cmap, norm, _ = elevation_cmap()
    dune_rgba = np.array(matplotlib.colors.to_rgba(C_DUNE))
    add_rgba = np.array(matplotlib.colors.to_rgba(C_ADD_FILL))
    cols = []
    for d in domains:
        topo = np.load(topo_dir / array_name("topography", d))
        n = n_by[d]
        ins = INSERT_ROW.get(d, 0) if ANCHOR == "road" else 0   # cells landward of row 0
        img = np.ones((ROWS_SHOWN, topo.shape[1], 4))
        img[:DUNE_ROWS] = dune_rgba
        z = cmap(norm(topo * CELL_M))
        # rows 0..ins-1 stay where they are; the block goes in at `ins`; the
        # rest of the interior is pushed down by N (add) or stays (remove -
        # the caller hatches the rows that go)
        head = min(ins, topo.shape[0])
        img[DUNE_ROWS:DUNE_ROWS + head] = z[:head]
        r0 = DUNE_ROWS + head + max(n, 0)
        img[DUNE_ROWS + head:r0] = add_rgba
        take = min(topo.shape[0] - head, ROWS_SHOWN - r0)
        if take > 0:
            img[r0:r0 + take] = z[head:head + take]
        cols.append(img)
    return np.concatenate(cols, axis=1)


INSERT_ROW: dict = {}     # domain -> insert_row_behind_road, filled by main()


def fig_grid(rows, topo_dir):
    off.apply_style()
    n_by = {r["domain"]: r["n_rows"] for r in rows}
    sb_by = {r["domain"]: r["setback_v2_m"] for r in rows}
    doms = sorted(n_by)
    groups = [doms[i:i + DOMAINS_PER_STRIP]
              for i in range(0, len(doms), DOMAINS_PER_STRIP)]
    ymin = -9.5                                   # room for labels + communities

    fig, axes = plt.subplots(len(groups), 1, figsize=(15.0, 5.2 * len(groups) + 1.4),
                             constrained_layout=True)
    axes = np.atleast_1d(axes)
    for k, (ax, g) in enumerate(zip(axes, groups)):
        ax.imshow(build_strip(g, topo_dir, n_by), aspect="auto", interpolation="nearest",
                  origin="upper", extent=[g[0] - 0.5, g[-1] + 0.5, ROWS_SHOWN, 0])
        for d in g:
            n = n_by[d]
            ax.axvline(d + 0.5, color="0.55", linewidth=0.35, zorder=3)
            ins = INSERT_ROW.get(d, 0) if ANCHOR == "road" else 0
            if n < 0:
                # the existing rows ins..ins+|n|-1, which the model would not hold
                ax.add_patch(Rectangle((d - 0.5, DUNE_ROWS + ins), 1.0, -n, facecolor="none",
                                       edgecolor=C_REM, hatch="//////", linewidth=0.0, zorder=4))
                ax.add_patch(Rectangle((d - 0.5, DUNE_ROWS + ins), 1.0, -n, facecolor="none",
                                       edgecolor=C_REM, linewidth=0.6, zorder=4))
            if n:
                ax.text(d, -0.9, f"{n:+d}", fontsize=7.2, ha="center", va="bottom",
                        color=C_ADD if n > 0 else C_REM, fontweight="bold")
            # NC-12 at its measured position (seaward edge, 20 m). It does not
            # move: with rows added it is pushed down with the interior; with
            # rows removed it stays where the surviving rows put it.
            sb = sb_by[d]
            if np.isfinite(sb):
                # dune anchor: the road is pushed down with the interior. Road
                # anchor: the rows go in BEHIND it, so it stays put.
                y = DUNE_ROWS + (max(n, 0) if ANCHOR == "dune" else 0) + sb / CELL_M
                ax.add_patch(Rectangle((d - 0.32, y), 0.64, 2.0, facecolor=C_ROAD,
                                       edgecolor="none", zorder=7))
        _community_bar(ax, -5.5, g[0] - 0.5, g[-1] + 0.5)
        ax.set_xlim(g[0] - 0.5, g[-1] + 0.5)
        ax.set_ylim(ROWS_SHOWN, ymin)
        ax.set_xticks([d for d in g if d % 5 == 0])
        ax.set_xticks(list(g), minor=True)
        ax.set_yticks(range(0, ROWS_SHOWN, 25))
        ax.set_ylabel("cross-shore cell\n(0 = the dune)")
        ax.grid(axis="y", color="0.9", linewidth=0.4)
        ax.set_axisbelow(True)
        sec = ax.secondary_yaxis("right", functions=(lambda c: c * CELL_M, lambda m: m / CELL_M))
        sec.set_ylabel("m landward of the dune")
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
                Patch(facecolor=C_DUNE, edgecolor="none", label=f"dune rows ({DUNE_ROWS})"),
                Patch(facecolor=C_ADD_FILL, edgecolor="none", label="rows added (blank, no fill yet)"),
                Patch(facecolor="none", edgecolor=C_REM, hatch="//////", label="existing rows removed"),
                Patch(facecolor=C_ROAD, edgecolor="none", label="NC-12 (1984), measured position"),
                Line2D([0], [0], color=off.HATTERAS_ANNOTATIONS.color_town_span, lw=4.0, label="community")]
    fig.legend(handles=handles, loc="outside lower center", ncol=7, fontsize=8,
               title="existing interior, elevation classes (m MHW)", title_fontsize=8)
    p = (FIG_SEAWARD / "HAT_row_insert_grid.png" if ANCHOR == "dune"
         else FIG_ROAD / "HAT_row_insert_grid_behindroad.png")
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


def fig_rows(rows):
    """Rows per domain, signed, on its own - the communities banded along the
    axis so a domain can be placed without the map."""
    off.apply_style()
    doms = np.array([r["domain"] for r in rows])
    nn = np.array([r["n_rows"] for r in rows])
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
    p = FIG_DIR / "HAT_row_insert_rows.png"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


# =============================================================================
# THE REPORT
# =============================================================================

def write_report(rows, topo_name, path, fig_path):
    n = np.array([r["n_rows"] for r in rows])
    add = [r for r in rows if r["n_rows"] > 0]
    rem = [r for r in rows if r["n_rows"] < 0]
    has_e = [r for r in rows if r["shift_m_easting"] == r["shift_m_easting"]]
    agree = [r for r in has_e if r["n_rows_easting"] == r["n_rows"]]
    differ = [r for r in has_e if r["n_rows_easting"] != r["n_rows"]]

    L = []
    w = L.append
    w("=" * 78)
    w("1984 ROW FOOTPRINT - SCOPE REPORT (symmetric: rows added AND removed)")
    w("=" * 78)
    w(f"written        {datetime.now():%Y-%m-%d %H:%M}")
    w(f"topography     {PRODUCT} / {topo_name}   (read, never modified)")
    w(f"N from         {FOOTPRINT_CSV.relative_to(REPO)}   (HAT_footprint_1984.py)")
    w(f"cross-check    {EASTING_CSV.relative_to(REPO)}")
    w(f"figure         {fig_path.relative_to(REPO)}")
    w(f"figure         {(FIG_DIR / 'HAT_row_insert_rows.png').relative_to(REPO)}")
    w(f"plan view      figures/2-footprint-1984/seaward/HAT_footprint_1984_plan.png  (HAT_footprint_1984.py)")
    w("")
    w("SCOPE ONLY. No array is written, no elevation is fabricated. This says")
    w("WHERE the interior grows or shrinks and by HOW MANY cells, so the fill")
    w("for the added rows can be argued separately against a known footprint.")
    w("")
    w("shift    = median over 50 paired profiles of (line_1997 - line_1984), m")
    w("N_rows   = trunc(shift / 10 m)          decided 2026-09-07")
    w("")
    w("  N > 0  the 1984 line is SEAWARD of the 1997 line; the island has")
    w("         retreated, so N rows are ADDED behind the dune.")
    w("  N < 0  the 1984 line is LANDWARD; the island has prograded, so the")
    w("         existing rows 0..|N|-1 are REMOVED. (Until 2026-09-07 this was")
    w("         floored to 0; the layers built that way, v3-v8, are gone.)")
    w("  |shift| < 10 m   no full cell of change is measured; left alone.")
    w("")
    w("-" * 78)
    w("SUMMARY")
    w("-" * 78)
    w(f"  domains examined                {len(rows)}")
    w(f"  rows ADDED     {len(add):2d} domains, {int(sum(r['n_rows'] for r in add)):3d} rows, "
      f"largest +{max(r['n_rows'] for r in add)} (domain {max(add, key=lambda r: r['n_rows'])['domain']})")
    w(f"  rows REMOVED   {len(rem):2d} domains, {int(-sum(r['n_rows'] for r in rem)):3d} rows, "
      f"largest {min(r['n_rows'] for r in rem)} (domain {min(rem, key=lambda r: r['n_rows'])['domain']})")
    w(f"  unchanged      {int((n == 0).sum()):2d} domains")
    w("")
    w("  N        domains")
    for v in sorted(set(n.tolist())):
        ids = [r["domain"] for r in rows if r["n_rows"] == v]
        w(f"  {v:+3d}  {len(ids):>4}   {ids if v else ''}")
    w("")
    w("-" * 78)
    w("CROSS-CHECK against the independent easting-frame measurement")
    w("-" * 78)
    w("  Same two geojsons, completely different frame: raw easting inside the")
    w("  axis-aligned domain box, 1 m alongshore sampling, no extractor c0, no")
    w("  shear, no row 0. It bounds the FRAME, not the lines. The same trunc")
    w("  rule is applied to it.")
    w("")
    d = np.array([r["shift_m"] - r["shift_m_easting"] for r in has_e])
    w(f"  identical N in                  {len(agree)} of {len(has_e)} domains")
    w(f"  differ                          {len(differ)}")
    w(f"  median difference in metres     {np.median(d):+.2f}")
    w(f"  mean absolute difference        {np.abs(d).mean():.2f} m "
      f"({np.abs(d).mean() / CELL_M:.2f} cells)")
    w(f"  largest difference              {np.abs(d).max():.1f} m")
    if differ:
        w("")
        w("  the domains where N differs:")
        for r in differ:
            w(f"    domain {r['domain']:>3}   extractor frame "
              f"{r['shift_m']:+7.1f} m -> {r['n_rows']:+d}   "
              f"easting frame {r['shift_m_easting']:+7.1f} m -> "
              f"{r['n_rows_easting']:+d}")
    w("")
    w("-" * 78)
    w("CAVEATS - read before using the footprint")
    w("-" * 78)
    w("  1. THE ADDED ROWS ARE FABRICATED. No survey covers land that had eroded")
    w("     away by 1996. This report sizes the block; it does not fill it.")
    w("  2. THE REMOVED ROWS ARE SURVEYED. Where the island prograded, rows")
    w("     0..|N|-1 are 1996 beach and foredune that were ocean in 1984;")
    w("     removing them deletes measured cells on the strength of the line pair.")
    w("  3. THE LINE IS 1997, THE DEM's BEACH IS 1996. One year of change is")
    w("     unaccounted for. Recorded, not corrected.")
    w("  4. THE TWO LINES MAY NOT BE THE SAME FEATURE. duneline_1997 records")
    w("     'light/dark elevation break'; duneline_1984 carries no metadata.")
    w("  5. THE 10 m RULE UNDERSTATES. trunc keeps 0-9 m less than measured at")
    w("     every changed domain, always toward less change. The residual is a")
    w("     column of the footprint table.")
    w("  6. N IS NOT A RATE. Row 0 is a 1996 feature only where ALACE reached")
    w("     the dune and 2009 elsewhere.")
    w("")
    w("-" * 78)
    w("PER DOMAIN")
    w("-" * 78)
    w(f"  {'dom':>4} {'N':>3} {'shift_m':>9} {'p10':>7} {'p90':>7} "
      f"{'easting_m':>10} {'N_e':>4} {'rows_now':>9} {'rows_after':>11}  note")
    for r in rows:
        note = {1: "ADD", -1: "REMOVE"}.get(int(np.sign(r["n_rows"])), "")
        e = (f"{r['shift_m_easting']:+10.1f}"
             if r["shift_m_easting"] == r["shift_m_easting"] else " " * 10)
        ne = (f"{r['n_rows_easting']:>+4d}"
              if r["shift_m_easting"] == r["shift_m_easting"] else "   -")
        w(f"  {r['domain']:>4} {r['n_rows']:>+3d} {r['shift_m']:>+9.1f} "
          f"{r['shift_p10_m']:>+7.1f} {r['shift_p90_m']:>+7.1f} {e} {ne} "
          f"{r['rows_now']:>9} {r['rows_after']:>11}  {note}")
    w("")
    w("=" * 78)
    path.write_text("\n".join(L) + "\n", encoding="utf-8")
    return "\n".join(L)


# =============================================================================
# MAIN
# =============================================================================

def main(base=None, anchor="dune"):
    global ANCHOR
    ANCHOR = anchor
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    topo_dir, _dune_dir, topo_name = topo_dirs(PRODUCT, override=base)
    print(f"topography : {PRODUCT} / {topo_name}")
    print(f"N from     : {FOOTPRINT_CSV.relative_to(REPO)}")

    fp = read_footprint(FOOTPRINT_CSV)
    easting = read_easting(EASTING_CSV)

    rows = []
    for d in sorted(fp):
        r = fp[d]
        f = topo_dir / array_name("topography", d)
        if not f.exists():
            print(f"  WARNING: no topography array for domain {d} - skipped")
            continue
        rows_now = int(np.load(f).shape[0])
        n = int(r["n_cells"])
        e = easting.get(d, float("nan"))
        rows.append({
            "domain": d,
            "n_rows": n,
            "shift_m": _f(r["shift_m_median"]),
            "shift_p10_m": _f(r["shift_m_p10"]),
            "shift_p90_m": _f(r["shift_m_p90"]),
            "shift_m_easting": e,
            "n_rows_easting": n_cells(e) if e == e else -99,
            "rows_now": rows_now,
            "rows_after": rows_now + n,
            "setback_v2_m": _f(r.get("setback_v2_m", "")),
            "insert_anchor": r.get("insert_anchor", ""),
            "insert_row_behind_road": _f(r.get("insert_row_behind_road", "")),
            "modified": int(n != 0),
        })
        if np.isfinite(_f(r.get("insert_row_behind_road", ""))):
            INSERT_ROW[d] = int(_f(r["insert_row_behind_road"]))

    cp = OUT_DIR / "row_insert_scope_by_domain.csv"
    with open(cp, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)

    fig = fig_grid(rows, topo_dir)
    if anchor == "road":
        # the table, the rows figure and the report are the same for both
        # placements (N does not change); only the grid is redrawn
        print(f"\n  figure : {fig}   (rows behind the road; table/report unchanged)")
        return
    fig_bars = fig_rows(rows)
    tp = OUT_DIR / "HAT_row_insert_scope.txt"
    text = write_report(rows, topo_name, tp, fig)
    print("\n" + "\n".join(text.splitlines()[:34]))
    print(f"\n  report : {tp}\n  table  : {cp}\n  figure : {fig}\n  figure : {fig_bars}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Scope the symmetric 1984 row footprint on the grid. "
                    "Writes nothing into a topography version.")
    ap.add_argument("--base", default=None,
                    help="topography version to read as the existing grid. "
                         "Default: whatever hat_topo_version resolves. N does "
                         "not depend on it; rows_now and the figure do.")
    ap.add_argument("--anchor", choices=("dune", "road"), default="dune",
                    help="dune: rows at the seaward edge (row 0 -> the 1984 line). "
                         "road: the same rows BEHIND the road, the crest-to-road "
                         "strip kept as measured (advisor's placement, 2026-09-07); "
                         "writes HAT_row_insert_grid_behindroad.png only")
    a = ap.parse_args()
    main(a.base, a.anchor)
