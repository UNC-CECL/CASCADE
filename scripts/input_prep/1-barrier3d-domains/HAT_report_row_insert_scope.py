#!/usr/bin/env python3
r"""
HAT_report_row_insert_scope.py
==============================================================================
Which domains would have interior rows inserted behind the dune, and how many.

SCOPE ONLY. Nothing is written into a topography version, no array is modified,
no elevation is fabricated. This answers one question - WHERE would the 1984
seaward-row insert touch the island, and by how many Barrier3D cells - so that
the fill rule can be argued about separately, against a known footprint.

The thing that actually performs the insert is HAT_insert_seaward_rows.py. This
script deliberately duplicates none of its fill logic; it reads the same
measurement that script's `--shift-source date` reads, and applies the same
N rule.

WHERE N COMES FROM
------------------
    per profile p, for each dune-line vintage y:
        shift_y(p) = ( interior_row0_cell(p) - duneline_cell_y(p) ) * 10 m
    domain value   = median over the 50 profiles
    N_metres       = shift_1984 - shift_1997  =  line_1997 - line_1984
    N_rows         = round( N_metres / 10 ), floored at 0

Read from `1984-start/row-insert-scope/duneline-shift/`, written by
HAT_measure_duneline_shift.py in the EXTRACTOR's own frame - its `c0` and its
per-profile shear - so N is counted in the cells CASCADE indexes rather than in
raw easting.

Two things cancel in that subtraction and both matter. Interior row 0 cancels
algebraically, so no assumption about where row 0 sits survives into N. And
since both lines were digitized to the same feature from imagery, the
DEFINITIONAL offset between "a digitized dune line" and "the model's row 0"
cancels with it. Island-wide the uncorrected 1984-vs-row-0 number is +18.9 m,
of which the 1997 line alone accounts for +16.2 m - about 85% was definitional.

SIGN, AND WHY NEGATIVES INSERT NOTHING
--------------------------------------
    N_metres > 0   the 1984 line sits SEAWARD of the 1997 line. The island has
                   retreated since 1984, so rows are ADDED behind the dune to
                   put interior row 0 back at the 1984 position.
    N_metres < 0   the 1984 line sits LANDWARD of the 1997 line.

A negative is floored to zero and the domain is left alone. It is NOT turned
into a removal. Deleting measured interior rows to represent apparent
progradation would destroy surveyed ground on the strength of a line pair whose
own comparability is unresolved, and the asymmetry is a choice, not an
oversight - so this report counts the negative domains, lists them, and prints
what they would have implied had the rule been symmetric.

THE CROSS-CHECK
---------------
`0-elevation/2009-2014-1996-duneline/duneline_offset_by_domain.csv` measures the
same 1984-minus-1997 quantity a completely different way: raw easting inside the
axis-aligned domain box, 1 m alongshore sampling, no extractor frame, no row 0,
no shear. Both are reported per domain and the disagreement is printed. They are
not independent of the same two geojsons, so agreement bounds the FRAME, not the
lines.

THIS IS NOT HYPOTHETICAL
------------------------
`dune-topo/v5` is already `v3` plus exactly these rows, built by
HAT_insert_seaward_rows.py with `--shift-source date --n-rule measured
--domains measured`. Where its audit CSV is on disk this report checks itself
against it domain by domain and says so. So the report is a scoping document
for a decision that has a built candidate, not a proposal in the abstract - and
the thing still open is the FILL, not the footprint.

WHICH VERSION IS THE BASE
-------------------------
Resolved through hat_topo_version.topo_dirs, whose order is env override, then
the extractor's own VERSION literal, then `dune-topo/CURRENT`. That order bites:
the extractor literal currently says v3 while CURRENT says v5, so CURRENT is
never consulted and v3 wins. Whatever it resolves to is printed and written into
the report header, and `--base` overrides it. N itself does not depend on the
choice - interior row 0 cancels out of it - but `rows_now` and the figure do.

OUTPUTS (1-barrier3d-domains/1984-start/row-insert-scope/)
    HAT_row_insert_scope.txt        the report
    row_insert_scope_by_domain.csv  the same table, machine-readable
    figures/HAT_row_insert_grid.png the Barrier3D grid with the inserted cells
                                    drawn in their own colour
    figures/HAT_row_insert_plan.png the same thing in PLAN VIEW, on the DEM,
                                    in the layout of the island offset map

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
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch


def _find_project_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit(f"cannot find data/hatteras_init above {start}")


REPO = _find_project_root(Path(__file__).resolve())
INIT = REPO / "data" / "hatteras_init"
sys.path.insert(0, str(REPO / "scripts"))
from hat_topo_version import duneline_shift_dir, topo_dirs  # noqa: E402


PRODUCT = "1984-start"
START_DIR = INIT / "1-barrier3d-domains" / PRODUCT
OUT_DIR = START_DIR / "row-insert-scope"
FIG_DIR = OUT_DIR / "figures"

# THE authoritative measurement - the same file HAT_insert_seaward_rows.py reads
# under `--shift-source date`. Using anything else here would create a second
# standard for a number that already has one.
RETREAT_CSV = duneline_shift_dir(PRODUCT) / "duneline_retreat_1984_1997.csv"
# The independent easting-frame measurement, for the cross-check column only.
EASTING_CSV = (INIT / "0-elevation" / "2009-2014-1996-duneline"
               / "duneline_offset_by_domain.csv")
# The version that already applied this insert, if it is on disk. Its audit is
# the strongest check available: same rule, run by the script that does the
# work, months before this report existed.
BUILT_VERSION = "v5"
BUILT_AUDIT = (START_DIR / "dune-topo" / BUILT_VERSION
               / "HAT_seaward_row_insert_audit.csv")

CELL_M = 10.0            # the Barrier3D cross-shore cell
DUNE_ROWS = 2            # DuneWidth; row 1 is a copy of row 0
SENTINEL_DAM = -0.30     # the extractor's water sentinel, in decametres

# Cross-shore rows drawn below the dune. The interiors run to 178 rows; showing
# all of them makes a 4500 x 180 image in which the inserted block - at most 7
# rows - is invisible. The cap is a drawing choice and is stated on the figure.
ROWS_SHOWN = 46
DOMAINS_PER_STRIP = 30

# Cell classes in the grid figure.
C_BG, C_DUNE, C_INSERT, C_LAND, C_WATER = 0, 1, 2, 3, 4
CLASS_COLOURS = ["#FFFFFF",   # off the array
                 "#c8a165",   # dune rows
                 "#e6194b",   # THE INSERTED ROWS - the one thing to look at
                 "#f0e6c8",   # existing interior, land
                 "#8fb8d8"]   # existing interior, at or below the sentinel
CLASS_LABELS = ["off the array", f"dune rows ({DUNE_ROWS})",
                "INSERTED rows", "existing interior, land",
                "existing interior, water"]


def read_retreat(path):
    if not path.exists():
        raise SystemExit(
            f"retreat measurement not found:\n    {path}\n"
            f"  Write it with HAT_measure_duneline_shift.py "
            f"(the 1984 vs 1997 pass).")
    out = {}
    for r in csv.DictReader(open(path)):
        v = r["shift_m_median"]
        if v not in ("", "nan"):
            out[int(r["domain"])] = float(v)
    return out


def read_easting(path):
    if not path.exists():
        print(f"  NOTE: {path.name} absent - cross-check column left blank")
        return {}
    out = {}
    for r in csv.DictReader(open(path)):
        v = r["offset_med_m"]
        if v not in ("", "nan"):
            out[int(r["domain"])] = float(v)
    return out


def read_built_audit(path):
    """n_rows_inserted per domain from the version that already did this."""
    if not path.exists():
        return {}
    return {int(r["domain"]): int(r["n_rows_inserted"])
            for r in csv.DictReader(open(path))}


def n_rows_from(metres):
    """round(N/10), floored at 0. The rule HAT_insert_seaward_rows.py applies."""
    return max(int(round(metres / CELL_M)), 0)


# =============================================================================
# THE FIGURE
# =============================================================================

def build_strip(domains, topo_dir, n_by_dom):
    """
    One alongshore strip of the Barrier3D grid, classes not elevations.

    Column order is the model's: 50 alongshore cells per domain, domains in
    increasing order. Rows run DOWN the page from the dune, which is the
    convention HAT_plot_b3d_grid.py established - the dune does not move when
    rows are inserted, the ground behind it grows and everything landward is
    pushed down the page.
    """
    cols = []
    for d in domains:
        topo = np.load(topo_dir / f"domain_{d}_topography.npy")
        n = n_by_dom[d]
        col = np.full((ROWS_SHOWN, topo.shape[1]), C_BG, np.uint8)
        col[:DUNE_ROWS, :] = C_DUNE
        r0 = DUNE_ROWS + n
        col[DUNE_ROWS:r0, :] = C_INSERT
        take = min(topo.shape[0], ROWS_SHOWN - r0)
        if take > 0:
            body = topo[:take, :]
            col[r0:r0 + take, :] = np.where(body > SENTINEL_DAM,
                                            C_LAND, C_WATER)
        cols.append(col)
    return np.hstack(cols)


def fig_grid(rows, topo_dir, topo_name):
    n_by_dom = {r["domain"]: r["n_rows"] for r in rows}
    doms = sorted(n_by_dom)
    groups = [doms[i:i + DOMAINS_PER_STRIP]
              for i in range(0, len(doms), DOMAINS_PER_STRIP)]

    cmap = ListedColormap(CLASS_COLOURS)
    norm = BoundaryNorm(np.arange(-0.5, len(CLASS_COLOURS) + 0.5), cmap.N)

    fig, axes = plt.subplots(
        len(groups) + 1, 1, figsize=(16.5, 3.1 * len(groups) + 3.0),
        gridspec_kw=dict(height_ratios=[1.0] * len(groups) + [1.25]),
        constrained_layout=True)

    for ax, g in zip(axes[:len(groups)], groups):
        img = build_strip(g, topo_dir, n_by_dom)
        ax.imshow(img, cmap=cmap, norm=norm, aspect="auto",
                  interpolation="nearest", origin="upper",
                  extent=[g[0] - 0.5, g[-1] + 0.5, ROWS_SHOWN, 0])
        for d in g:
            ax.axvline(d + 0.5, color="0.45", linewidth=0.4)
            if n_by_dom[d] > 0:
                ax.text(d, -1.2, str(n_by_dom[d]), fontsize=7, ha="center",
                        va="bottom", color=CLASS_COLOURS[C_INSERT],
                        fontweight="bold")
        ax.set_xticks([d for d in g if d % 5 == 0])
        ax.set_xticklabels([str(d) for d in g if d % 5 == 0], fontsize=8)
        ax.set_ylabel("cross-shore cell\n(0 = behind the dune)", fontsize=8)
        ax.tick_params(labelsize=8)
        ax.set_xlim(g[0] - 0.5, g[-1] + 0.5)

    bx = axes[-1]
    nn = np.array([n_by_dom[d] for d in doms])
    bx.bar(doms, nn, width=0.8, color=CLASS_COLOURS[C_INSERT],
           edgecolor="none")
    bx.set_xlabel("domain (1 = south)")
    bx.set_ylabel("rows inserted")
    bx.set_xlim(doms[0] - 0.5, doms[-1] + 0.5)
    bx.set_yticks(range(0, int(nn.max()) + 1))
    bx.grid(axis="y", color="0.9", linewidth=0.6)
    bx.set_axisbelow(True)
    bx.set_title(f"{int((nn > 0).sum())} of {len(doms)} domains modified, "
                 f"{int(nn.sum())} rows in total, largest {int(nn.max())} "
                 f"(domain {doms[int(np.argmax(nn))]})", fontsize=10)

    axes[0].legend(handles=[Patch(facecolor=c, edgecolor="0.5", label=l)
                            for c, l in zip(CLASS_COLOURS, CLASS_LABELS)],
                   loc="upper center", bbox_to_anchor=(0.5, 1.42), ncol=5,
                   fontsize=8, framealpha=0.95)

    fig.suptitle(
        f"Where the 1984 seaward-row insert would touch the Barrier3D grid   "
        f"|   topography {topo_name}, scope only - no fill", fontsize=13)
    fig.text(0.5, -0.005,
             f"Cells, not elevations. The dune stays put and the ground behind "
             f"it grows, so the existing interior is pushed down the page by N "
             f"rows. Only the first {ROWS_SHOWN} cross-shore cells are drawn "
             f"(interiors run to 178) and cells are not square. N per domain is "
             f"round(1984 line - 1997 line / 10 m), floored at 0.",
             ha="center", va="top", fontsize=8, color="#444444", wrap=True)

    p = FIG_DIR / "HAT_row_insert_grid.png"
    fig.savefig(p, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return p


# =============================================================================
# THE SAME THING IN PLAN VIEW
# =============================================================================

def _plan_deps():
    """
    The 0-elevation figure module, imported rather than reimplemented.

    It already owns the island mosaic loader, the km axis formatter, the
    elevation panel, the dune-line styling and the three-panel split. Copying
    any of that here would give the two figures two ways to draw one island.
    Returns None if it is not importable, and the plan figure is skipped rather
    than the whole report failing.
    """
    d = (REPO / "scripts" / "input_prep" / "0-elevation" / "3-figures")
    if str(d) not in sys.path:
        sys.path.insert(0, str(d))
    try:
        import HAT_plot_duneline_offset as off  # noqa: E402
        return off
    except Exception as exc:                     # noqa: BLE001
        print(f"  NOTE: plan-view figure skipped - {type(exc).__name__}: {exc}")
        return None


def fig_plan(rows, off):
    """
    Where the inserts land, on the map, in the island-offset figure's layout.

    TWO ENCODINGS, BECAUSE ONE CANNOT WORK ALONE. A one-row insert is 10 m. On a
    panel 3 km wide that is a third of one percent - sub-pixel. So:

      * each DOMAIN BOX is shaded by N. Readable at island scale, and it is what
        answers "which domains would be modified" at a glance.
      * the true-scale BAND is drawn on top: the 1997 dune line offset seaward
        by exactly N x 10 m. Honest about size, visible only on a zoom.

    WHY THE BAND IS ANCHORED ON THE 1997 LINE. The rows go in at the seaward end
    of the interior, ahead of current row 0, and row 0 was picked from a surface
    whose dune is the 1996/1997 one. So 1997 is the map-space proxy for where
    the existing array starts, and the band runs from there toward the 1984
    position. It stops at N x 10 m rather than at the 1984 line itself because N
    is rounded to whole cells - the gap between the band's outer edge and the
    red line is the rounding, drawn rather than hidden.
    """
    print("  plan view: loading the island mosaic and re-sampling the lines ...")
    elev, _surv, extent, _n = off.m.load_mosaic()
    # off.load_domains(), not a bare read: domains.geojson lives on D: and
    # that drive can vanish mid-session. The loader falls back to rebuilding
    # the boxes from the resampled rasters, which reproduce it exactly.
    gdf = off.load_domains()
    lines = off.load_lines(gdf.crs)
    geoms = {yr: g.union_all() for yr, g in lines.items()}
    _rows, samples = off.measure(gdf, geoms)
    drawn = off.clip_for_drawing(lines, gdf.union_all())

    n_by = {r["domain"]: r["n_rows"] for r in rows}
    nmax = max(n_by.values())
    # Discrete, one step per row. A continuous ramp would imply N is continuous
    # and would make 1 and 2 rows indistinguishable, which is the distinction
    # the figure exists to carry.
    shades = plt.get_cmap("YlOrRd")(np.linspace(0.25, 0.95, nmax))

    vmin, vmax = off.m.elev_limits(elev)
    ids = gdf["domain_id"].astype(int).to_numpy()
    groups = np.array_split(np.sort(ids), off.N_PANELS)

    fig, axes = plt.subplots(1, off.N_PANELS,
                             figsize=(4.6 * off.N_PANELS, 17.0),
                             constrained_layout=True)
    im = None
    for ax, group in zip(np.atleast_1d(axes), groups):
        sub = gdf[gdf["domain_id"].astype(int).isin(group)]
        bx = sub.total_bounds
        im = off.m.panel_elev(ax, elev, extent, vmin, vmax,
                              f"domains {group.min()}-{group.max()}")

        for _, r in sub.iterrows():
            d = int(r["domain_id"])
            n = n_by.get(d, 0)
            b = r.geometry.bounds
            if n > 0:
                ax.add_patch(plt.Rectangle(
                    (b[0], b[1]), b[2] - b[0], b[3] - b[1],
                    facecolor=shades[n - 1], edgecolor="none", alpha=0.55,
                    zorder=3))
            ax.add_patch(plt.Rectangle(
                (b[0], b[1]), b[2] - b[0], b[3] - b[1], facecolor="none",
                edgecolor="0.25", linewidth=0.4, zorder=4))
            if n > 0:
                ax.text(b[0] - 90, (b[1] + b[3]) / 2, f"{d}  +{n}", fontsize=6.5,
                        ha="right", va="center", color="#7f0000",
                        fontweight="bold", zorder=6)
            elif d % off.LABEL_EVERY == 0 or d in (1, ids.max()):
                ax.text(b[0] - 90, (b[1] + b[3]) / 2, str(d), fontsize=6,
                        ha="right", va="center", color="0.45", zorder=6)

            # The true-scale band, from the 1997 line seaward by N * 10 m.
            if n > 0:
                m_ = samples["domain"] == d
                y = samples["y"][m_]
                x97 = samples["x1997"][m_]
                ok = np.isfinite(x97)
                if ok.sum() > 2:
                    poly = np.concatenate([
                        np.column_stack([x97[ok], y[ok]]),
                        np.column_stack([(x97[ok] + n * CELL_M)[::-1],
                                         y[ok][::-1]])])
                    ax.add_patch(plt.Polygon(
                        poly, closed=True, facecolor=CLASS_COLOURS[C_INSERT],
                        edgecolor=CLASS_COLOURS[C_INSERT], linewidth=0.5,
                        alpha=0.95, zorder=8))

        off.draw_lines(ax, drawn, scale=0.5)
        ax.set_xlim(bx[0] - off.PANEL_PAD_M, bx[2] + off.PANEL_PAD_M)
        ax.set_ylim(bx[1] - off.PANEL_PAD_M, bx[3] + off.PANEL_PAD_M)
        ax.set_aspect("equal")
        off.m.km_axes(ax, nx=3, ny=8)
        ax.set_xlabel("Easting (km)", fontsize=9)
    np.atleast_1d(axes)[0].set_ylabel("Northing (km)", fontsize=9)

    cb = fig.colorbar(im, ax=list(np.atleast_1d(axes)),
                      orientation="horizontal", fraction=0.02, pad=0.01,
                      aspect=55)
    cb.set_label("Elevation (m NAVD88)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    handles = [Patch(facecolor=shades[v - 1], alpha=0.55, edgecolor="0.4",
                     label=f"{v} row{'s' if v > 1 else ''}"
                           f" ({v * int(CELL_M)} m)")
               for v in range(1, nmax + 1)]
    handles += [Patch(facecolor=CLASS_COLOURS[C_INSERT], edgecolor="none",
                      label=f"inserted ground, true scale"),
                Patch(facecolor="none", edgecolor="0.25",
                      label="domain, no insert")]
    # Upper left, not lower: the southern domains 1-6 all carry inserts and
    # their labels sit at the bottom-left of the first panel, where a legend
    # lands on top of them.
    np.atleast_1d(axes)[0].legend(
        handles=handles + off.line_legend(), loc="upper left", fontsize=7.5,
        ncol=2, framealpha=0.94, title="rows inserted", title_fontsize=8)

    n_arr = np.array([n_by[d] for d in sorted(n_by)])
    fig.suptitle("Where the 1984 seaward-row insert would land, in plan view",
                 fontsize=13)
    fig.text(0.5, -0.008,
             f"{int((n_arr > 0).sum())} of {len(n_arr)} domains, "
             f"{int(n_arr.sum())} rows, largest {int(n_arr.max())}. Box shading "
             f"is N; the solid band is the SAME N at true scale - the 1997 dune "
             f"line offset seaward by N x 10 m - so one row is a 10 m sliver "
             f"and needs a zoom. The band stops short of the 1984 line by the "
             f"rounding in N = round((1984 line - 1997 line) / 10). Domains with "
             f"N = 0 are unshaded and untouched.",
             ha="center", va="top", fontsize=8, color="#444444", wrap=True)

    p = FIG_DIR / "HAT_row_insert_plan.png"
    fig.savefig(p, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return p


# =============================================================================
# THE REPORT
# =============================================================================

def version_provenance(resolved):
    """Why this version, in one line - the resolution order is a known trap."""
    cur = START_DIR / "dune-topo" / "CURRENT"
    said = cur.read_text(encoding="utf-8").strip() if cur.is_file() else "(none)"
    if said and said != resolved:
        return (f"resolved to {resolved}; dune-topo/CURRENT says {said}. "
                f"topo_dirs reads the extractor's VERSION literal BEFORE "
                f"CURRENT, so CURRENT was not consulted. N does not depend on "
                f"this; rows_now and the figure do.")
    return f"resolved to {resolved}; dune-topo/CURRENT agrees."


def write_report(rows, topo_name, path, fig_path):
    n = np.array([r["n_rows"] for r in rows])
    sh = np.array([r["shift_m"] for r in rows])
    mod = [r for r in rows if r["n_rows"] > 0]
    neg = [r for r in rows if r["shift_m"] < 0]
    agree = [r for r in rows if r["n_rows_easting"] == r["n_rows"]]
    differ = [r for r in rows if r["n_rows_easting"] != r["n_rows"]]

    L = []
    w = L.append
    w("=" * 78)
    w("1984 SEAWARD ROW INSERT - SCOPE REPORT")
    w("=" * 78)
    w(f"written        {datetime.now():%Y-%m-%d %H:%M}")
    w(f"topography     {PRODUCT} / {topo_name}   (read, never modified)")
    w(f"               {version_provenance(topo_name)}")
    w(f"N from         {RETREAT_CSV.relative_to(REPO)}")
    w(f"cross-check    {EASTING_CSV.relative_to(REPO)}")
    w(f"figure         {fig_path.relative_to(REPO)}")
    w("")
    w("SCOPE ONLY. No array is written, no elevation is fabricated. This says")
    w("WHERE the insert would land and HOW MANY cells, so the fill rule can be")
    w("argued separately against a known footprint.")
    w("")
    w("N_metres = shift_1984 - shift_1997 = line_1997 - line_1984")
    w("N_rows   = round(N_metres / 10), floored at 0")
    w("")
    w("  N > 0  the 1984 line is SEAWARD of the 1997 line; the island has")
    w("         retreated, so rows are added behind the dune.")
    w("  N < 0  the 1984 line is LANDWARD; floored to 0, domain left alone.")
    w("         NOT converted into a removal - see NEGATIVES below.")
    w("")
    w("-" * 78)
    w("SUMMARY")
    w("-" * 78)
    w(f"  domains examined                {len(rows)}")
    w(f"  domains modified (N >= 1)       {len(mod)}")
    w(f"  domains untouched (N == 0)      {len(rows) - len(mod)}")
    w(f"  rows inserted, total            {int(n.sum())}")
    w(f"  largest insert                  {int(n.max())} rows "
      f"(domain {rows[int(np.argmax(n))]['domain']}, "
      f"{rows[int(np.argmax(n))]['shift_m']:+.1f} m)")
    w(f"  mean over modified domains      {n[n > 0].mean():.2f} rows")
    w("")
    w("  N        domains")
    for v in range(int(n.max()) + 1):
        ids = [r["domain"] for r in rows if r["n_rows"] == v]
        w(f"  {v:>2}  {len(ids):>4}   {ids if v else ''}")
    w("")
    w("-" * 78)
    w("NEGATIVES - measured, floored, and left alone")
    w("-" * 78)
    w(f"  {len(neg)} of {len(rows)} domains have the 1984 line LANDWARD of the")
    w("  1997 line. Every one is floored to N = 0 and its topography is")
    w("  untouched. Had the rule been symmetric they would have implied")
    w(f"  {sum(abs(int(round(r['shift_m'] / CELL_M))) for r in neg)} rows of "
      f"REMOVAL, i.e. deleting surveyed interior cells.")
    w("")
    w("  That is not done. Progradation is not what a floored measurement")
    w("  licenses, and the line pair's own comparability is unresolved - see")
    w("  CAVEATS. The asymmetry is a choice and this is where it is recorded.")
    w("")
    w("  worst negatives:")
    for r in sorted(neg, key=lambda x: x["shift_m"])[:8]:
        w(f"    domain {r['domain']:>3}   {r['shift_m']:+7.1f} m   "
          f"({r['shift_m'] / CELL_M:+5.2f} cells)")
    w("")
    w("-" * 78)
    w("CROSS-CHECK against the independent easting-frame measurement")
    w("-" * 78)
    w("  Same two geojsons, completely different frame: raw easting inside the")
    w("  axis-aligned domain box, 1 m alongshore sampling, no extractor c0, no")
    w("  shear, no row 0. It bounds the FRAME, not the lines.")
    w("")
    d = np.array([r["shift_m"] - r["shift_m_easting"] for r in rows
                  if r["shift_m_easting"] == r["shift_m_easting"]])
    w(f"  identical N in                  {len(agree)} of {len(rows)} domains")
    w(f"  differ by one row in            {len(differ)}")
    w(f"  median difference in metres     {np.median(d):+.2f}")
    w(f"  mean absolute difference        {np.abs(d).mean():.2f} m "
      f"({np.abs(d).mean() / CELL_M:.2f} cells)")
    w(f"  largest difference              {np.abs(d).max():.1f} m")
    if differ:
        w("")
        w("  the domains where N differs by one:")
        for r in differ:
            w(f"    domain {r['domain']:>3}   extractor frame "
              f"{r['shift_m']:+7.1f} m -> {r['n_rows']}   "
              f"easting frame {r['shift_m_easting']:+7.1f} m -> "
              f"{r['n_rows_easting']}")
    w("")
    w("-" * 78)
    w(f"AGAINST {BUILT_VERSION}, WHICH ALREADY APPLIED THIS INSERT")
    w("-" * 78)
    built = read_built_audit(BUILT_AUDIT)
    if not built:
        w(f"  {BUILT_AUDIT.name} not on disk - check skipped.")
    else:
        same = [r for r in rows if built.get(r["domain"], 0) == r["n_rows"]]
        off = [(r["domain"], built.get(r["domain"], 0), r["n_rows"])
               for r in rows if built.get(r["domain"], 0) != r["n_rows"]]
        w(f"  dune-topo/{BUILT_VERSION} is the base plus exactly this insert,")
        w("  built by HAT_insert_seaward_rows.py --shift-source date")
        w("  --n-rule measured --domains measured. Its audit says:")
        w("")
        w(f"    domains with N >= 1   {sum(1 for v in built.values() if v > 0)}"
          f"   (this report: {len(mod)})")
        w(f"    rows inserted         {sum(built.values())}"
          f"   (this report: {int(n.sum())})")
        w(f"    per-domain agreement  {len(same)} of {len(rows)}")
        if off:
            w("")
            w("    domains that DISAGREE (domain, built, this report):")
            for t in off:
                w(f"      {t}")
            w("    A disagreement means the measurement underneath has changed")
            w(f"    since {BUILT_VERSION} was built. Trust neither number until")
            w("    that is explained.")
        else:
            w("")
            w("    Exact, every domain. The footprint in this report IS the")
            w(f"    footprint {BUILT_VERSION} already has on disk - so what is")
            w("    open is the fill and whether to adopt it, not the scope.")
    w("")
    w("-" * 78)
    w("AGAINST THE PUBLISHED BLOCK-SCOPE TABLE")
    w("-" * 78)
    w("  dune-topo/v4 is the same insert at the ten relocation-block domains")
    w("  only, and its published N values were:")
    w("")
    pub = {9: 0, 10: 1, 11: 2, 12: 2, 13: 1, 14: 0,
           84: 4, 85: 6, 86: 3, 87: 1}
    by = {r["domain"]: r for r in rows}
    ok = 0
    for k in sorted(pub):
        if k not in by:
            continue
        got = by[k]["n_rows"]
        ok += got == pub[k]
        w(f"    domain {k:>3}   published {pub[k]}   recomputed {got}   "
          f"{'ok' if got == pub[k] else 'MISMATCH'}")
    w("")
    w(f"  {ok} of {len(pub)} reproduce. v4 inserts at 8 of these 10 (9 and 14")
    w("  round to zero). This report's scope is wider only because it applies")
    w("  the same rule to all 90 domains rather than to the ten with a")
    w("  documented NC-12 relocation.")
    w("")
    w("-" * 78)
    w("CAVEATS - read before using the footprint")
    w("-" * 78)
    w("  1. THE ROWS ARE FABRICATED. No survey covers land that had eroded away")
    w("     by 1996. This report sizes the block; it does not fill it, and the")
    w("     fill rule is a separate decision with its own defence.")
    w("")
    w("  2. THE LINE IS 1997, THE DEM's BEACH IS 1996. One year of retreat is")
    w("     unaccounted for, under half a cell at the fastest measured rate.")
    w("     Recorded, not corrected - scaling by 12/13 would assume retreat was")
    w("     steady across 13 storm-driven years.")
    w("")
    w("  3. THE TWO LINES MAY NOT BE THE SAME FEATURE. duneline_1997 records")
    w("     its method as 'digitized from light/dark elevation break';")
    w("     duneline_1984 carries no metadata at all. In the true-scale zooms")
    w("     at 0-elevation/2009-2014-1996-duneline/figures/, the two excursion")
    w("     reaches look like both lines straddling ONE ridge from opposite")
    w("     sides and swapping which side they take - 1984 landward in 62-68,")
    w("     seaward in 78-85 - which is what a definitional difference looks")
    w("     like and not what one-directional retreat looks like. That")
    w("     observation is unresolved and it bears directly on N.")
    w("")
    w("  4. N IS NOT A RATE. Row 0 is a 1996 feature only where ALACE reached")
    w("     the dune and 2009 elsewhere, so the interval behind a per-domain N")
    w("     is not the same everywhere.")
    w("")
    w("-" * 78)
    w("PER DOMAIN")
    w("-" * 78)
    w(f"  {'dom':>4} {'N':>3} {'shift_m':>9} {'cells':>7} "
      f"{'easting_m':>10} {'N_e':>4} {'rows_now':>9} {'rows_after':>11}  note")
    for r in rows:
        note = ""
        if r["n_rows"] > 0:
            note = "INSERT"
        elif r["shift_m"] < 0:
            note = "negative, floored"
        e = (f"{r['shift_m_easting']:+10.1f}"
             if r["shift_m_easting"] == r["shift_m_easting"] else " " * 10)
        ne = (f"{r['n_rows_easting']:>4}"
              if r["shift_m_easting"] == r["shift_m_easting"] else "   -")
        w(f"  {r['domain']:>4} {r['n_rows']:>3} {r['shift_m']:>+9.1f} "
          f"{r['shift_m'] / CELL_M:>+7.2f} {e} {ne} "
          f"{r['rows_now']:>9} {r['rows_after']:>11}  {note}")
    w("")
    w("=" * 78)

    path.write_text("\n".join(L) + "\n", encoding="utf-8")
    return "\n".join(L)


# =============================================================================
# MAIN
# =============================================================================

def main(base=None):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    topo_dir, _dune_dir, topo_name = topo_dirs(PRODUCT, override=base)
    print(f"topography : {PRODUCT} / {topo_name}")
    print(f"             {version_provenance(topo_name)}")
    print(f"N from     : {RETREAT_CSV.relative_to(REPO)}")

    retreat = read_retreat(RETREAT_CSV)
    easting = read_easting(EASTING_CSV)

    rows = []
    for d in sorted(retreat):
        f = topo_dir / f"domain_{d}_topography.npy"
        if not f.exists():
            print(f"  WARNING: no topography array for domain {d} - skipped")
            continue
        rows_now = int(np.load(f).shape[0])
        n = n_rows_from(retreat[d])
        e = easting.get(d, float("nan"))
        rows.append({
            "domain": d,
            "shift_m": retreat[d],
            "n_rows": n,
            "shift_m_easting": e,
            "n_rows_easting": n_rows_from(e) if e == e else -1,
            "rows_now": rows_now,
            "rows_after": rows_now + n,
            "modified": int(n > 0),
        })

    cp = OUT_DIR / "row_insert_scope_by_domain.csv"
    with open(cp, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)

    fp = fig_grid(rows, topo_dir, topo_name)
    off = _plan_deps()
    pp = fig_plan(rows, off) if off is not None else None

    tp = OUT_DIR / "HAT_row_insert_scope.txt"
    text = write_report(rows, topo_name, tp, fp)
    print("\n" + "\n".join(text.splitlines()[:40]))

    print(f"\n  report : {tp}")
    print(f"  table  : {cp}")
    print(f"  figure : {fp}")
    if pp:
        print(f"  figure : {pp}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Scope the 1984 seaward row insert. Writes nothing into "
                    "a topography version.")
    ap.add_argument("--base", default=None,
                    help="topography version to read as the existing grid. "
                         "Default: whatever hat_topo_version resolves, which "
                         "is printed. N does not depend on it.")
    main(ap.parse_args().base)
