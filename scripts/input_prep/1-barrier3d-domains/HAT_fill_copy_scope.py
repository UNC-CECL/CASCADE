#!/usr/bin/env python3
r"""
HAT_fill_copy_scope.py
==============================================================================
What the added rows would CONTAIN under the copy fill, for the behind-the-road
placement of the 1984 footprint - scoped, drawn and audited, no version written.

THE FILL (Hannah's advisor; decided with Hannah 2026-09-07)
    The block of N rows goes in directly behind the model's two roadway rows
    AS PLACED under the 1984 setback (insert_row_behind_road in
    footprint_1984_by_domain.csv = int(setback_new/10) + 2; 2026-09-08). It is filled
    with a DIRECT COPY of the N interior rows immediately landward of the
    insert point - rows r..r+N-1 of the existing interior, in order, cell by
    cell across the 50 alongshore columns - so the block fabricates no value
    and reads as the backbarrier it stands beside. The window follows N.

    * seams: the seaward junction is continuous by construction (the block's
      first row is the row that used to follow the road). The only seam is at
      the LANDWARD end, where the block's last row (a copy of r+N-1) meets the
      original row r; `seam_jump_m` is that step, alongshore mean of |dz|.
    * no-road domains (GIS 2-5): the footprint puts the block BEHIND THE
      CREST ROW (anchor "crest", insert_row = crest + 1, 2026-09-08), so the
      rule is the same as behind the road - the window is the N rows that
      follow the insert point, the crest stays at the front, nothing is
      duplicated. (Until 2026-09-08 the block went in at row 0 and the window
      skipped the crest, which left the crest row stranded between block and
      source.)
    * outliers: GIS 82/83 windows hold 10-11 m cells (Rodanthe structures in
      the DEM). Copied as measured - they are already in the interior - and
      flagged.
    * water: no window holds a cell at or below MHW, so no rule is needed and
      none is applied. If a future footprint changes that, `window_water_frac`
      is the column to watch.
    * removals need no fill; they are listed for completeness with no window.

WHAT IS ASSUMED
    * the lost 1984 ground is booked behind the road, so the backbarrier next
      to the road is the analogue for it (what was lost was ocean-side beach
      and dune);
    * the 1996/2009 backbarrier behind the road stands for 1984 backbarrier
      (no accretion or subsidence in between);
    * alongshore structure in the window appears twice cross-shore;
    * the dune array, row 0, the road rows and the model's setback are as in
      the behind-road placement - untouched.

OUTPUTS  row-insert-scope/
    fill_copy_by_domain.csv          per domain: N, insert row, source rows,
                                     window stats, seam jump, flags
    HAT_fill_copy_scope.txt          the report
    figures/3-fill/HAT_fill_copy_grid_GIS<...>.png
                                     for the example domains: the near-road
                                     interior before and after, in elevation classes
    figures/3-fill/HAT_fill_copy_method_GIS<...>.png
                                     the method in three stages, model frame: the
                                     domain as extracted (dune rows + interior),
                                     the N rows inserted blank behind NC-12, the
                                     copy fill with the source window and the copy
                                     drawn as an arrow

USAGE
    python HAT_fill_copy_scope.py                    # examples 80, 85, 5, 49
    python HAT_fill_copy_scope.py --domains 80,73
==============================================================================
"""
from __future__ import annotations

import argparse
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


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
INIT = REPO / "data" / "hatteras_init"
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "0-elevation" / "3-figures"))
from hat_topo_version import array_name, dune_topo_root, insert_figures_dir, topo_dirs  # noqa: E402
from hat_figure_style import elevation_cmap  # noqa: E402
import HAT_plot_duneline_offset as off  # noqa: E402  the house style

PRODUCT = "1984-start"
CELL_M = 10.0
ROAD_ROWS = 2
BERM_EL_M = 1.7               # BermEl, Hatteras-CASCADE-parameters.yaml; the dune file is height above it
CREST_SEARCH_ROWS = 10        # the crest is looked for in the first rows of a no-road domain
OUTLIER_M = 8.0               # a copied cell above this is flagged (structures, not ground)
SCOPE_DIR = INIT / "1-barrier3d-domains" / PRODUCT / "row-insert-scope"
FOOTPRINT_CSV = SCOPE_DIR / "footprint_1984_by_domain.csv"
FIG_DIR = insert_figures_dir(PRODUCT, "3-fill")
BUILT_VERSION = "v3"          # the version HAT_build_footprint_version.py wrote; after-panels read it if present
C_ROAD_OLD = "0.35"
INK = off.INK
C_ADD, C_ADD_FILL = off.C_1984, off.C_1984_FILL
C_ROAD = "#1a1a1a"


# =============================================================================
# THE FILL, PER DOMAIN
# =============================================================================

def plan_fill(tab: pd.DataFrame, topo_dir: Path) -> tuple[pd.DataFrame, dict]:
    """The audit table, and per domain the arrays (before, after, source rows)."""
    rows, arrays = [], {}
    for d, r in tab.iterrows():
        n = int(r["n_cells"])
        if n == 0:
            continue
        z = np.load(topo_dir / array_name("topography", d)) * CELL_M     # m MHW
        rec = {"domain": int(d), "n_cells": n, "insert_anchor": r["insert_anchor"],
               "insert_row": int(r["insert_row_behind_road"])}
        if n < 0:
            rec.update({"fill": "none (rows removed)", "source_rows": "", "crest_row": np.nan,
                        "window_mean_m": np.nan, "window_min_m": np.nan, "window_max_m": np.nan,
                        "window_water_frac": np.nan, "seam_jump_m": np.nan, "flags": ""})
            rows.append(rec)
            continue
        ins = rec["insert_row"]
        flags = []
        # one rule for every anchor: the window is the N rows that follow the
        # insert point. For a no-road domain the footprint already put that
        # point behind the crest row.
        src0 = ins
        if r["insert_anchor"] == "crest":
            rec["crest_row"] = int(r["crest_row"])
            flags.append(f"NO_ROAD_behind_crest(row {int(r['crest_row'])})")
        else:
            rec["crest_row"] = np.nan
        src = z[src0:src0 + n]
        if src.shape[0] < n:
            flags.append("WINDOW_SHORT")
        block = src.copy()
        after = np.concatenate([z[:ins], block, z[ins:]], axis=0)
        # If the version has been built, the after-panel IS that version's
        # array - read it and check it is what the rule says.
        built = dune_topo_root(PRODUCT) / BUILT_VERSION / "topography" / array_name("topography", d)
        if built.is_file():
            zb = np.load(built) * CELL_M
            assert np.array_equal(zb, after), (d, "built version differs from the rule")
            after = zb
        seam = float(np.mean(np.abs(block[-1] - z[ins]))) if z.shape[0] > ins else np.nan
        rec.update({"fill": "copy of the next N rows", "source_rows": f"{src0}..{src0 + n - 1}",
                    "window_mean_m": round(float(block.mean()), 2),
                    "window_min_m": round(float(block.min()), 2),
                    "window_max_m": round(float(block.max()), 2),
                    "window_water_frac": round(float((block <= 0).mean()), 3),
                    "seam_jump_m": round(seam, 2)})
        if block.max() > OUTLIER_M:
            flags.append(f"OUTLIER_COPIED({block.max():.1f} m)")
        if (block <= 0).any():
            flags.append("WATER_IN_WINDOW")
        rec["flags"] = ",".join(flags)
        rows.append(rec)
        arrays[int(d)] = dict(before=z, after=after, insert=ins, n=n, src0=src0,
                              road_row=(int(r["setback_model_now_m"] // CELL_M)
                                        if np.isfinite(r["setback_model_now_m"]) else None),
                              # where v3's 1984 setback puts the model's road (on the block)
                              road_row_new=(int(r["setback_new_m"] // CELL_M)
                                            if np.isfinite(r["setback_new_m"]) else None))
    return pd.DataFrame(rows).set_index("domain"), arrays


# =============================================================================
# THE FIGURE: before / after / profile, for one domain
# =============================================================================

def _draw_road(ax, y: float, ncol: int, label: str, colour: str) -> None:
    ax.add_patch(Rectangle((-0.5, y - 0.5), ncol, ROAD_ROWS, facecolor=colour,
                           edgecolor=colour, lw=1.2, alpha=0.45, zorder=5))
    ax.text(ncol - 1.0, y + ROAD_ROWS / 2 - 0.5, label, ha="right", va="center",
            fontsize=7.5, fontweight="bold", color="white", zorder=6)


def fig_domain(d: int, a: dict, rows_shown: int = 40) -> Path:
    """Before and after, side by side, the road rows marked in both."""
    off.apply_style()
    cmap, norm, bounds = elevation_cmap()
    z0, z1, ins, n, src0, road = a["before"], a["after"], a["insert"], a["n"], a["src0"], a["road_row"]
    road_new = a["road_row_new"]
    R = min(z1.shape[0], max(rows_shown, ins + n + 14))     # the block and the road always in view
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 5.6), constrained_layout=True)
    for k, (ax, z, title) in enumerate(zip(axes, (z0, z1),
                                           ("as extracted (v2)", f"{BUILT_VERSION}: the {n}-row block, copy fill, 1984 setback"))):
        ax.imshow(z[:R], cmap=cmap, norm=norm, aspect="auto", interpolation="nearest",
                  origin="upper", extent=[-0.5, z.shape[1] - 0.5, R - 0.5, -0.5])
        if road is not None:
            # (a) NC-12 as the model holds it today, two rows at int(setback/10);
            # (b) the same rows outlined as the measured pavement, and the model
            #     road where v3's 1984 setback puts it - N rows inland, on the block
            if k == 0:
                _draw_road(ax, road, z.shape[1], "NC-12", C_ROAD)
            else:
                ax.add_patch(Rectangle((-0.5, road - 0.5), z.shape[1], ROAD_ROWS, facecolor="none",
                                       edgecolor=C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)), zorder=5))
                if road_new is not None:
                    _draw_road(ax, road_new, z.shape[1], f"NC-12 ({BUILT_VERSION})", C_ROAD)
        if k == 0:
            ax.add_patch(Rectangle((-0.5, src0 - 0.5), z.shape[1], n, facecolor="none",
                                   edgecolor=C_ADD, lw=1.4, ls=(0, (3, 2)), zorder=5))
            ax.axhline(ins - 0.5, color=C_ADD, lw=1.2, zorder=5)
        else:
            ax.add_patch(Rectangle((-0.5, ins - 0.5), z.shape[1], n, facecolor="none",
                                   edgecolor=C_ADD, lw=1.6, zorder=5))
        ax.set_xlabel("alongshore cell")
        if k == 0:
            ax.set_ylabel("cross-shore cell (0 = interior row 0)")
        ax.set_yticks(range(0, R, 5))
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(True)
        off._title(ax, k, f"GIS {d}: {title}")

    labels = ["below 0 (water)"] + [f"{lo:g}–{hi:g}" for lo, hi in zip(bounds[1:-2], bounds[2:-1])]         + [f"above {bounds[-2]:g}"]
    handles = [Patch(facecolor=cmap(i), edgecolor="0.4", lw=0.4, label=lab) for i, lab in enumerate(labels)]
    handles += [Line2D([0], [0], color=C_ADD, lw=1.4, ls=(0, (3, 2)), label="source window (copied)"),
                Line2D([0], [0], color=C_ADD, lw=1.6, label="the block"),
                Patch(facecolor=C_ROAD, alpha=0.45, edgecolor=C_ROAD, label="NC-12 rows as the model holds them"),
                Patch(facecolor="none", edgecolor=C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)),
                      label="measured pavement rows (today's setback)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=5, fontsize=8,
               title="elevation classes (m MHW)", title_fontsize=8)
    p = FIG_DIR / f"HAT_fill_copy_grid_GIS{d}.png"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


def fig_method(d: int, a: dict, dune_dir: Path, rows_shown: int = 40) -> Path:
    """
    The method in three stages, one domain, all in the MODEL's frame (dune rows
    on top, then the interior): (a) the Barrier3D domain as extracted, (b) the
    N rows identified by the footprint inserted BLANK behind the roadway rows,
    (c) the same rows filled by copying the N interior rows that follow them.
    """
    off.apply_style()
    cmap, norm, bounds = elevation_cmap()
    z0, z1, ins, n, src0, road = a["before"], a["after"], a["insert"], a["n"], a["src0"], a["road_row"]
    road_new = a["road_row_new"]
    dune = np.load(dune_dir / array_name("dune", d)) * CELL_M + BERM_EL_M     # m MHW
    dune_rows = np.tile(dune[None, :], (ROAD_ROWS, 1))                       # DuneWidth = 2 rows
    ncol = z0.shape[1]
    R = min(z1.shape[0], max(rows_shown, ins + n + 14)) + ROAD_ROWS

    def stack(z):
        return np.concatenate([dune_rows, z], axis=0)[:R]

    blank = np.full((n, ncol), np.nan)
    z_blank = np.concatenate([z0[:ins], blank, z0[ins:]], axis=0)
    panels = ((stack(z0), "as extracted (v2): dune rows + interior"),
              (stack(z_blank), f"{n} rows inserted blank behind {'NC-12 as placed' if road is not None else 'the crest row'}"),
              (stack(z1), f"{BUILT_VERSION}: the {n} rows filled, road at its 1984 setback"))
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 5.8), constrained_layout=True)
    y_road = None if road is None else road + ROAD_ROWS          # the dune rows shift everything by 2
    y_ins = ins + ROAD_ROWS
    for k, (ax, (img, title)) in enumerate(zip(axes, panels)):
        rgba = cmap(norm(img))
        rgba[np.isnan(img)] = (1, 1, 1, 1)
        ax.imshow(rgba, aspect="auto", interpolation="nearest", origin="upper",
                  extent=[-0.5, ncol - 0.5, R - 0.5, -0.5])
        # the dune rows, marked
        ax.axhline(ROAD_ROWS - 0.5, color=INK, lw=1.0, zorder=5)
        ax.text(ncol - 1.0, ROAD_ROWS / 2 - 0.5, "dune", ha="right", va="center", fontsize=7.5,
                fontweight="bold", color="white", zorder=6)
        if road is not None:
            if k == 0:
                _draw_road(ax, y_road, ncol, "NC-12 today", C_ROAD)
            else:
                # (b), (c): the model road at its 1984 setback - the block sits
                # directly behind it; the measured pavement rows outlined
                ax.add_patch(Rectangle((-0.5, y_road - 0.5), ncol, ROAD_ROWS, facecolor="none",
                                       edgecolor=C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)), zorder=5))
                if road_new is not None:
                    _draw_road(ax, road_new + ROAD_ROWS, ncol, "NC-12 at its 1984 setback", C_ROAD)
        if k == 0:
            ax.axhline(y_ins - 0.5, color=C_ADD, lw=1.4, zorder=5)
            ax.text(0.5, y_ins + 0.1, f"insert point: interior row {ins}", ha="left", va="top",
                    fontsize=7.5, color=C_ADD, fontweight="bold", zorder=6,
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.15"))
        else:
            ax.add_patch(Rectangle((-0.5, y_ins - 0.5), ncol, n, facecolor="none", edgecolor=C_ADD,
                                   lw=1.6, hatch="////" if k == 1 else None, zorder=5))
            if k == 1:
                ax.text(ncol / 2, y_ins + n / 2 - 0.5, f"{n} blank rows", ha="center", va="center",
                        fontsize=8, fontweight="bold", color=C_ADD, zorder=6,
                        bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.2"))
        if k == 2:
            # the source window, in the after frame, and the copy drawn as an arrow.
            # The source starts at src0, which is the insert point for a road
            # domain but one row behind the crest for a no-road domain (GIS 2-5)
            # - so it is src0 + n after the insert, not ins + n.
            y_src = src0 + n + ROAD_ROWS
            ax.add_patch(Rectangle((-0.5, y_src - 0.5), ncol, n, facecolor="none", edgecolor=C_ADD,
                                   lw=1.4, ls=(0, (3, 2)), zorder=5))
            ax.annotate("", xy=(ncol * 0.5, y_ins + n / 2 - 0.5), xytext=(ncol * 0.5, y_src + n / 2 - 0.5),
                        arrowprops=dict(arrowstyle="-|>", color=C_ADD, lw=1.6, mutation_scale=14,
                                        connectionstyle="arc3,rad=-0.4"), zorder=7)
            ax.text(ncol * 0.5 + 3, y_src + n / 2 - 0.5, "copied", ha="left", va="center", fontsize=7.5,
                    color=C_ADD, fontweight="bold", zorder=7)
        ax.set_xlabel("alongshore cell")
        if k == 0:
            ax.set_ylabel("cross-shore cell (0 = the dune)")
        ax.set_yticks(range(0, R, 5))
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(True)
        off._title(ax, k, f"GIS {d}: {title}")

    labels = ["below 0 (water)"] + [f"{lo:g}\u2013{hi:g}" for lo, hi in zip(bounds[1:-2], bounds[2:-1])] \
        + [f"above {bounds[-2]:g}"]
    handles = [Patch(facecolor=cmap(i), edgecolor="0.4", lw=0.4, label=lab) for i, lab in enumerate(labels)]
    handles += [Patch(facecolor="white", edgecolor=C_ADD, hatch="////", label="rows inserted, blank"),
                Line2D([0], [0], color=C_ADD, lw=1.6, label="the block, filled"),
                Line2D([0], [0], color=C_ADD, lw=1.4, ls=(0, (3, 2)), label="source window (the next N rows)"),
                Patch(facecolor=C_ROAD, alpha=0.45, edgecolor=C_ROAD, label="NC-12 rows as the model holds them"),
                Patch(facecolor="none", edgecolor=C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)),
                      label="measured pavement rows (today's setback)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=5, fontsize=8,
               title="elevation classes (m MHW); dune rows drawn at berm + dune height", title_fontsize=8)
    p = FIG_DIR / f"HAT_fill_copy_method_GIS{d}.png"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


# =============================================================================

def write_report(aud: pd.DataFrame, topo_name: str, figs: list[Path]) -> Path:
    add = aud[aud.n_cells > 0]
    L = []
    w = L.append
    w("=" * 78)
    w("COPY FILL for the rows behind the road - scope (no version written)")
    w("=" * 78)
    w(f"written        {datetime.now():%Y-%m-%d %H:%M}")
    w(f"topography     {PRODUCT} / {topo_name}   (read, never modified)")
    w(f"footprint      {FOOTPRINT_CSV.relative_to(REPO)}  (insert_row_behind_road, n_cells)")
    for f in figs:
        w(f"figure         {f.relative_to(REPO)}")
    w("")
    w("RULE (decided 2026-09-07)")
    w("  block = a direct copy of interior rows r..r+N-1, the N rows immediately landward")
    w("  of the insert point r = int(setback/10) + 2, cell by cell. No-road domains: r is")
    w("  the row behind the crest row. Outliers copied as measured and flagged. No")
    w("  window holds water, so no water rule is applied.")
    w("")
    w("SUMMARY")
    w(f"  domains filled                 {len(add)}   ({int(add.n_cells.sum())} rows copied)")
    w(f"  window elevation, m MHW        median of domain means {add.window_mean_m.median():.2f}; "
      f"range of means {add.window_mean_m.min():.2f}..{add.window_mean_m.max():.2f}")
    w(f"  seam at the landward end       mean {add.seam_jump_m.mean():.2f} m, max {add.seam_jump_m.max():.2f} m "
      f"(GIS {int(add.seam_jump_m.idxmax())})")
    w(f"  windows with water             {int((add.window_water_frac > 0).sum())}")
    w(f"  outliers copied                {add.index[add['flags'].str.contains('OUTLIER')].tolist()}")
    w(f"  no-road, block behind crest    {add.index[add.insert_anchor == 'crest'].tolist()}")
    w("")
    w("PER DOMAIN")
    w("  domain  N  anchor  insert  source rows  crest  mean   min    max   seam  flags")
    for d, r in aud.iterrows():
        if r.n_cells < 0:
            w(f"  {d:6d} {r.n_cells:+3d}  {r.insert_anchor:5s}  {r.insert_row:6d}  {'-':>11s}  {'-':>5s}  "
              f"{'-':>5s} {'-':>5s} {'-':>6s} {'-':>5s}  rows removed, no fill")
            continue
        crest = "-" if not np.isfinite(r.crest_row) else f"{int(r.crest_row)}"
        w(f"  {d:6d} {r.n_cells:+3d}  {r.insert_anchor:5s}  {r.insert_row:6d}  {r.source_rows:>11s}  {crest:>5s}  "
          f"{r.window_mean_m:5.2f} {r.window_min_m:5.2f} {r.window_max_m:6.2f} {r.seam_jump_m:5.2f}  {r['flags']}")
    w("")
    w("ASSUMED: the backbarrier beside the road stands for the lost 1984 ground (which was")
    w("ocean-side beach and dune); it has not accreted or subsided since 1984; alongshore")
    w("structure appears twice; dune array, row 0, road rows and setback untouched.")
    p = SCOPE_DIR / "HAT_fill_copy_scope.txt"
    p.write_text("\n".join(L) + "\n", encoding="utf-8")
    return p


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--domains", default="80,85,5,49",
                    help="example domains to draw (comma list)")
    args = ap.parse_args()
    if not FOOTPRINT_CSV.is_file():
        raise SystemExit(f"\n{FOOTPRINT_CSV} not found - run HAT_footprint_1984.py first\n")
    tab = pd.read_csv(FOOTPRINT_CSV).set_index("domain")
    topo_dir, dune_dir, topo_name = topo_dirs(PRODUCT)
    print(f"topography : {PRODUCT} / {topo_name}")
    aud, arrays = plan_fill(tab, topo_dir)
    aud.to_csv(SCOPE_DIR / "fill_copy_by_domain.csv")
    print(f"wrote {SCOPE_DIR / 'fill_copy_by_domain.csv'}")
    figs = []
    for d in [int(x) for x in args.domains.split(",") if x.strip()]:
        if d not in arrays:
            print(f"  GIS {d}: no rows added - nothing to fill")
            continue
        figs.append(fig_domain(d, arrays[d]))
        print(f"wrote {figs[-1]}")
        figs.append(fig_method(d, arrays[d], dune_dir))
        print(f"wrote {figs[-1]}")
    rep = write_report(aud, topo_name, figs)
    print(f"wrote {rep}")
    print(open(rep, encoding="utf-8").read().split("PER DOMAIN")[0])


if __name__ == "__main__":
    main()
