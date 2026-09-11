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
    * removals need no fill. The |N| rows directly seaward of today's roadway
      rows are deleted (insert_row = int(setback_v2/10) - |N|; Hannah,
      2026-09-08: the rows come out of the interior in front of the road);
      the after-panel is read from the built version and checked.

WHAT IS ASSUMED
    * the lost 1984 ground is booked behind the road, so the backbarrier next
      to the road is the analogue for it (what was lost was ocean-side beach
      and dune);
    * the 1996/2009 backbarrier behind the road stands for 1984 backbarrier
      (no accretion or subsidence in between);
    * alongshore structure in the window appears twice cross-shore;
    * the dune array, row 0, the road rows and the model's setback are as in
      the behind-road placement - untouched.

OUTPUTS  2-domain-reconstruction-1984/
    fill_copy_by_domain.csv          per domain: N, insert row, source rows,
                                     window stats, seam jump, flags
    HAT_fill_copy_scope.txt          the report
    figures/4-fill/rows-{added,removed}/HAT_fill_copy_grid_GIS<...>.png
                                     for the example domains: the near-road
                                     interior before and after, in elevation classes
    figures/4-fill/rows-{added,removed}/HAT_fill_copy_method_GIS<...>.png
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
from hat_topo_version import array_name, dune_topo_root, insert_figures_dir_for_domain, topo_dirs, insert_scope_step# noqa: E402
from hat_figure_style import C, caption, elevation_cmap, figsize, save  # noqa: E402
import HAT_plot_duneline_offset as off  # noqa: E402  the house style

PRODUCT = "1984-start"
CELL_M = 10.0
ROAD_ROWS = 2
BERM_EL_M = 1.7               # BermEl, Hatteras-CASCADE-parameters.yaml; the dune file is height above it
CREST_SEARCH_ROWS = 10        # the crest is looked for in the first rows of a no-road domain
OUTLIER_M = 8.0               # a copied cell above this is flagged (structures, not ground)
SCOPE_DIR = INIT / "1-barrier3d-domains" / PRODUCT / "2-domain-reconstruction-1984"
FOOTPRINT_CSV = insert_scope_step(PRODUCT, "2-extent") / "footprint_1984_by_domain.csv"
BUILT_VERSION = "v3"          # the version HAT_build_footprint_version.py wrote; after-panels read it if present
C_ROAD_OLD = C["BASE"]
INK = off.INK
C_ADD, C_ADD_FILL = off.C_1984, off.C_1984_FILL
C_REM = off.C_1997
C_ROAD = C["ROAD"]


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
            ins = rec["insert_row"]
            after = np.concatenate([z[:ins], z[ins - n:]], axis=0)
            built = dune_topo_root(PRODUCT) / BUILT_VERSION / "topography" / array_name("topography", d)
            if built.is_file():
                zb = np.load(built) * CELL_M
                assert np.array_equal(zb, after), (d, "built version differs from the rule")
                after = zb
            seam = float(np.mean(np.abs(z[ins - 1] - z[ins - n]))) if ins > 0 else np.nan
            rec.update({"fill": "none (rows removed in front of NC-12)", "source_rows": f"{ins}..{ins - n - 1} removed",
                        "crest_row": np.nan, "window_mean_m": np.nan, "window_min_m": np.nan,
                        "window_max_m": np.nan, "window_water_frac": np.nan,
                        "seam_jump_m": round(seam, 2), "flags": ""})
            rows.append(rec)
            arrays[int(d)] = dict(before=z, after=after, insert=ins, n=n, src0=None,
                                  road_row=(int(r["setback_model_now_m"] // CELL_M)
                                            if np.isfinite(r["setback_model_now_m"]) else None),
                                  road_row_new=(int(r["setback_new_m"] // CELL_M)
                                                if np.isfinite(r["setback_new_m"]) else None))
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
            fontsize=7, fontweight="bold", color="white", zorder=6)


def fig_domain(d: int, a: dict, rows_shown: int = 40) -> Path:
    """Before and after, side by side, the road rows marked in both."""
    off.apply_style()
    cmap, norm, bounds = elevation_cmap()
    z0, z1, ins, n, src0, road = a["before"], a["after"], a["insert"], a["n"], a["src0"], a["road_row"]
    road_new = a["road_row_new"]
    if n < 0:
        return _fig_domain_removal(d, a, rows_shown)
    R = min(z1.shape[0], max(rows_shown, ins + n + 14))     # the block and the road always in view
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", aspect=0.50), constrained_layout=True)
    for k, (ax, z, title) in enumerate(zip(axes, (z0, z1),
                                           ("original domain", f"+{n} rows, filled"))):
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
                    _draw_road(ax, road_new, z.shape[1], "NC-12, 1984 setback", C_ROAD)
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
        ax.set_yticks(range(0, R, 10))
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(True)
        off._title(ax, k, f"GIS {d}: {title}")

    labels = ["< 0 m (water)"] + [f"{lo:g}–{hi:g} m" for lo, hi in zip(bounds[1:-2], bounds[2:-1])]         + [f"> {bounds[-2]:g} m"]
    handles = [Patch(facecolor=cmap(i), edgecolor="0.4", lw=0.4, label=lab) for i, lab in enumerate(labels)]
    handles += [Line2D([0], [0], color=C_ADD, lw=1.4, ls=(0, (3, 2)), label="rows copied into them"),
                Line2D([0], [0], color=C_ADD, lw=1.6, label="rows inserted, filled"),
                Patch(facecolor=C_ROAD, alpha=0.45, edgecolor=C_ROAD, label="NC-12 at the 1984 setback"),
                Patch(facecolor="none", edgecolor=C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)),
                      label="NC-12 at the setback measured on the 1996 surface")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, fontsize=7, frameon=False,
               title="interior elevation (m above MHW)", title_fontsize=7)
    caption(fig, f"GIS {d}, the near-road interior in elevation classes, alongshore across and "
                 f"cross-shore down from interior row 0. (a) The domain as extracted from the 1996 "
                 f"surface, with NC-12 at the setback measured there and the {n} rows the fill copies "
                 f"outlined dashed; the solid line is the insertion point. (b) The same domain with "
                 f"{n} rows inserted directly behind NC-12 as the model places it under the 1984 setback "
                 f"and filled with a cell-by-cell copy of the {n} rows immediately landward of them, so "
                 f"the block fabricates no elevation. The road rows at the setback measured on the 1996 "
                 f"surface are outlined for comparison.")
    p = insert_figures_dir_for_domain(PRODUCT, "4-fill", d) / f"HAT_fill_copy_grid_GIS{d}.png"
    save(fig, p, vector=False, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


def _fig_domain_removal(d: int, a: dict, rows_shown: int = 40) -> Path:
    """Before and after for a removal domain: the |N| rows directly seaward of
    today's roadway rows go; the road's cells and everything behind are kept."""
    off.apply_style()
    cmap, norm, bounds = elevation_cmap()
    z0, z1, ins, n, road, road_new = a["before"], a["after"], a["insert"], a["n"], a["road_row"], a["road_row_new"]
    m = -n
    R = min(z0.shape[0], max(rows_shown, (road or 0) + 16))
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", aspect=0.50), constrained_layout=True)
    for k, (ax, z, title) in enumerate(zip(axes, (z0, z1),
                                           ("original domain", f"−{m} rows"))):
        ax.imshow(z[:R], cmap=cmap, norm=norm, aspect="auto", interpolation="nearest",
                  origin="upper", extent=[-0.5, z.shape[1] - 0.5, R - 0.5, -0.5])
        if k == 0:
            if road is not None:
                _draw_road(ax, road, z.shape[1], "NC-12, measured setback", C_ROAD)
            ax.add_patch(Rectangle((-0.5, ins - 0.5), z.shape[1], m, facecolor="none", edgecolor=C_REM,
                                   lw=1.6, hatch="////", zorder=5))
            ax.text(z.shape[1] / 2, ins + m / 2 - 0.5, f"rows {ins}\u2013{ins + m - 1} removed",
                    ha="center", va="center", fontsize=7, fontweight="bold", color=C_REM, zorder=6,
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.2"))
        else:
            if road is not None:
                # the old pavement rows, where the surviving rows put them
                ax.add_patch(Rectangle((-0.5, ins - 0.5), z.shape[1], ROAD_ROWS, facecolor="none",
                                       edgecolor=C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)), zorder=5))
            if road_new is not None:
                _draw_road(ax, road_new, z.shape[1], "NC-12, 1984 setback", C_ROAD)
            ax.axhline(ins - 0.5, xmax=0.60, color=C_REM, lw=2.0, ls=(0, (3, 1.5)), zorder=8)
            ax.annotate(f"seam: {m} rows removed",
                        xy=(z.shape[1] * 0.35, ins - 0.5), xytext=(z.shape[1] * 0.35, (road_new or ins) + ROAD_ROWS + 3.5),
                        ha="center", va="top", fontsize=7, color=C_REM, fontweight="bold", zorder=9,
                        arrowprops=dict(arrowstyle="-|>", color=C_REM, lw=1.2, mutation_scale=12, shrinkB=0),
                        bbox=dict(facecolor="white", alpha=0.9, edgecolor=C_REM, lw=0.8, boxstyle="square,pad=0.25"))
        ax.set_xlabel("alongshore cell")
        if k == 0:
            ax.set_ylabel("cross-shore cell (0 = interior row 0)")
        ax.set_yticks(range(0, R, 10))
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(True)
        off._title(ax, k, f"GIS {d}: {title}")
    labels = ["< 0 m (water)"] + [f"{lo:g}\u2013{hi:g} m" for lo, hi in zip(bounds[1:-2], bounds[2:-1])] \
        + [f"> {bounds[-2]:g} m"]
    handles = [Patch(facecolor=cmap(i), edgecolor="0.4", lw=0.4, label=lab) for i, lab in enumerate(labels)]
    handles += [Patch(facecolor="none", edgecolor=C_REM, hatch="////", label="rows removed seaward of NC-12"),
                Line2D([0], [0], color=C_REM, lw=2.0, ls=(0, (3, 1.5)), label="seam left by the removal"),
                Patch(facecolor=C_ROAD, alpha=0.45, edgecolor=C_ROAD, label="NC-12 at the 1984 setback"),
                Patch(facecolor="none", edgecolor=C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)),
                      label="NC-12 at the measured setback, after the removal")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, fontsize=7, frameon=False,
               title="interior elevation (m above MHW)", title_fontsize=7)
    caption(fig, f"GIS {d}, the near-road interior in elevation classes, alongshore across and "
                 f"cross-shore down from interior row 0. (a) The domain as extracted from the 1996 "
                 f"surface, with NC-12 at the setback measured there and the {m} rows directly seaward of "
                 f"the roadway rows hatched: those are the rows the footprint removes, the island having "
                 f"prograded here since 1984. (b) The same domain with those rows gone. Nothing is "
                 f"filled: the surviving rows close up, leaving one seam at the removal, and the road, "
                 f"its own cells and everything landward of them are kept as measured.")
    p = insert_figures_dir_for_domain(PRODUCT, "4-fill", d) / f"HAT_fill_copy_grid_GIS{d}.png"
    save(fig, p, vector=False, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


def _fig_method_removal(d: int, a: dict, dune_dir: Path, rows_shown: int = 40) -> Path:
    """The three stages for a removal domain, model frame: as extracted, the
    |N| rows identified (in front of NC-12), the rows removed and the road at
    its 1984 setback."""
    off.apply_style()
    cmap, norm, bounds = elevation_cmap()
    z0, z1, ins, n, road, road_new = a["before"], a["after"], a["insert"], a["n"], a["road_row"], a["road_row_new"]
    m = -n
    dune = np.load(dune_dir / array_name("dune", d)) * CELL_M + BERM_EL_M
    dune_rows = np.tile(dune[None, :], (ROAD_ROWS, 1))
    ncol = z0.shape[1]
    R = min(z0.shape[0], max(rows_shown, (road or 0) + 16)) + ROAD_ROWS

    def stack(z):
        return np.concatenate([dune_rows, z], axis=0)[:R]

    panels = ((stack(z0), "original domain"),
              (stack(z0), f"{m} rows identified"),
              (stack(z1), f"the {m} rows removed"))
    fig, axes = plt.subplots(1, 3, figsize=figsize("double", aspect=0.44), constrained_layout=True)
    y_road = None if road is None else road + ROAD_ROWS
    y_ins = ins + ROAD_ROWS
    for k, (ax, (img, title)) in enumerate(zip(axes, panels)):
        rgba = cmap(norm(img))
        rgba[np.isnan(img)] = (1, 1, 1, 1)
        ax.imshow(rgba, aspect="auto", interpolation="nearest", origin="upper",
                  extent=[-0.5, ncol - 0.5, R - 0.5, -0.5])
        ax.axhline(ROAD_ROWS - 0.5, color=INK, lw=1.0, zorder=5)
        ax.text(ncol - 1.0, ROAD_ROWS / 2 - 0.5, "dune", ha="right", va="center", fontsize=7.5,
                fontweight="bold", color="white", zorder=6)
        if k < 2 and road is not None:
            _draw_road(ax, y_road, ncol, "NC-12, measured", C_ROAD)
        if k == 1:
            ax.add_patch(Rectangle((-0.5, y_ins - 0.5), ncol, m, facecolor="none", edgecolor=C_REM,
                                   lw=1.6, hatch="////", zorder=5))
            ax.text(1.0, y_ins + m / 2 - 0.5, f"rows {ins}\u2013{ins + m - 1}",
                    ha="left", va="center", fontsize=7, fontweight="bold", color=C_REM, zorder=6,
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.2"))
        if k == 2:
            if road is not None:
                ax.add_patch(Rectangle((-0.5, y_ins - 0.5), ncol, ROAD_ROWS, facecolor="none",
                                       edgecolor=C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)), zorder=5))
            if road_new is not None:
                _draw_road(ax, road_new + ROAD_ROWS, ncol, "NC-12, 1984 setback", C_ROAD)
            ax.axhline(y_ins - 0.5, xmax=0.60, color=C_REM, lw=2.0, ls=(0, (3, 1.5)), zorder=8)
            ax.annotate(f"seam: {m} rows removed",
                        xy=(ncol * 0.35, y_ins - 0.5), xytext=(ncol * 0.35, (road_new or ins) + 2 * ROAD_ROWS + 3.5),
                        ha="center", va="top", fontsize=7, color=C_REM, fontweight="bold", zorder=9,
                        arrowprops=dict(arrowstyle="-|>", color=C_REM, lw=1.2, mutation_scale=12, shrinkB=0),
                        bbox=dict(facecolor="white", alpha=0.9, edgecolor=C_REM, lw=0.8, boxstyle="square,pad=0.25"))
        ax.set_xlabel("alongshore cell")
        if k == 0:
            ax.set_ylabel("cross-shore cell (0 = the dune)")
        ax.set_yticks(range(0, R, 10))
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(True)
        off._title(ax, k, title)
    labels = ["< 0 m (water)"] + [f"{lo:g}\u2013{hi:g} m" for lo, hi in zip(bounds[1:-2], bounds[2:-1])] \
        + [f"> {bounds[-2]:g} m"]
    handles = [Patch(facecolor=cmap(i), edgecolor="0.4", lw=0.4, label=lab) for i, lab in enumerate(labels)]
    handles += [Patch(facecolor="none", edgecolor=C_REM, hatch="////", label="rows removed seaward of NC-12"),
                Line2D([0], [0], color=C_REM, lw=2.0, ls=(0, (3, 1.5)), label="seam left by the removal"),
                Patch(facecolor=C_ROAD, alpha=0.45, edgecolor=C_ROAD, label="NC-12 at the 1984 setback"),
                Patch(facecolor="none", edgecolor=C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)),
                      label="NC-12 at the measured setback, after the removal")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, fontsize=7, frameon=False,
               title="interior elevation (m above MHW); dune rows drawn at berm + dune height",
               title_fontsize=7)
    caption(fig, f"GIS {d} in the model's own frame: the two dune rows on top, drawn at berm + dune "
                 f"height, then every interior row, in elevation classes. (a) The domain as extracted "
                 f"from the 1996 surface, NC-12 at the setback measured there. (b) The {m} rows the "
                 f"footprint removes, the |N| rows directly seaward of the roadway rows. (c) The domain "
                 f"with those rows gone and NC-12 at its 1984 setback, which lands on the old pavement's "
                 f"first row or the row seaward of it; the seam is where the surviving rows close up. "
                 f"Nothing is filled and no elevation is fabricated.")
    p = insert_figures_dir_for_domain(PRODUCT, "4-fill", d) / f"HAT_fill_copy_method_GIS{d}.png"
    save(fig, p, vector=False, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


def fig_method(d: int, a: dict, dune_dir: Path, rows_shown: int = 40) -> Path:
    """
    The method in three stages, one domain, all in the MODEL's frame (dune rows
    on top, then the interior): (a) the Barrier3D domain as extracted, (b) the
    N rows identified by the footprint inserted BLANK behind the roadway rows,
    (c) the same rows filled by copying the N interior rows that follow them.
    """
    if a["n"] < 0:
        return _fig_method_removal(d, a, dune_dir, rows_shown)
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
    panels = ((stack(z0), "original domain"),
              (stack(z_blank), f"{n} rows inserted (unfilled)"),
              (stack(z1), f"the {n} rows filled by copying"))
    fig, axes = plt.subplots(1, 3, figsize=figsize("double", aspect=0.44), constrained_layout=True)
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
                _draw_road(ax, y_road, ncol, "NC-12, measured", C_ROAD)
            else:
                # (b), (c): the model road at its 1984 setback - the block sits
                # directly behind it; the measured pavement rows outlined
                ax.add_patch(Rectangle((-0.5, y_road - 0.5), ncol, ROAD_ROWS, facecolor="none",
                                       edgecolor=C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)), zorder=5))
                if road_new is not None:
                    _draw_road(ax, road_new + ROAD_ROWS, ncol, "NC-12, 1984 setback", C_ROAD)
        if k == 0:
            ax.axhline(y_ins - 0.5, color=C_ADD, lw=1.4, zorder=5)
            ax.text(0.5, y_ins + 0.1, f"insertion point: row {ins}", ha="left", va="top",
                    fontsize=7, color=C_ADD, fontweight="bold", zorder=6,
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.15"))
        else:
            ax.add_patch(Rectangle((-0.5, y_ins - 0.5), ncol, n, facecolor="none", edgecolor=C_ADD,
                                   lw=1.6, hatch="////" if k == 1 else None, zorder=5))
            if k == 1:
                ax.text(1.0, y_ins + n / 2 - 0.5, f"+{n} rows", ha="left", va="center",
                        fontsize=7, fontweight="bold", color=C_ADD, zorder=6,
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
            ax.text(ncol * 0.5 + 3, y_src + n / 2 - 0.5, "copied", ha="left", va="center", fontsize=7,
                    color=C_ADD, fontweight="bold", zorder=7)
        ax.set_xlabel("alongshore cell")
        if k == 0:
            ax.set_ylabel("cross-shore cell (0 = the dune)")
        ax.set_yticks(range(0, R, 10))
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(True)
        off._title(ax, k, title)

    labels = ["< 0 m (water)"] + [f"{lo:g}\u2013{hi:g} m" for lo, hi in zip(bounds[1:-2], bounds[2:-1])] \
        + [f"> {bounds[-2]:g} m"]
    handles = [Patch(facecolor=cmap(i), edgecolor="0.4", lw=0.4, label=lab) for i, lab in enumerate(labels)]
    handles += [Patch(facecolor="white", edgecolor=C_ADD, hatch="////", label="rows inserted (unfilled)"),
                Line2D([0], [0], color=C_ADD, lw=1.6, label="rows inserted, filled"),
                Line2D([0], [0], color=C_ADD, lw=1.4, ls=(0, (3, 2)), label="rows copied into them"),
                Patch(facecolor=C_ROAD, alpha=0.45, edgecolor=C_ROAD, label="NC-12 at the 1984 setback"),
                Patch(facecolor="none", edgecolor=C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)),
                      label="NC-12 at the setback measured on the 1996 surface")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, fontsize=7, frameon=False,
               title="interior elevation (m above MHW); dune rows drawn at berm + dune height",
               title_fontsize=7)
    anchor = ("directly behind NC-12 as the model places it under the 1984 setback"
              if road is not None else "directly behind the dune crest row, this domain having no model road")
    caption(fig, f"GIS {d} in the model's own frame: the two dune rows on top, drawn at berm + dune "
                 f"height, then every interior row, in elevation classes. (a) The domain as extracted "
                 f"from the 1996 surface, NC-12 at the setback measured there; the line is the insertion "
                 f"point, {anchor}. (b) The {n} rows the footprint adds, inserted there and left blank. "
                 f"(c) The same rows filled with a cell-by-cell copy of the {n} interior rows immediately "
                 f"landward of them (dashed, arrow), so the block fabricates no elevation and reads as "
                 f"the backbarrier it stands beside. The only seam is at the landward end of the block.")
    p = insert_figures_dir_for_domain(PRODUCT, "4-fill", d) / f"HAT_fill_copy_method_GIS{d}.png"
    save(fig, p, vector=False, bbox_inches="tight", facecolor="white")
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
            w(f"  {d:6d} {r.n_cells:+3d}  {r.insert_anchor:5s}  {r.insert_row:6d}  {r.source_rows:>11s}  {'-':>5s}  "
              f"{'-':>5s} {'-':>5s} {'-':>6s} {r.seam_jump_m:5.2f}  rows removed in front of NC-12, no fill")
            continue
        crest = "-" if not np.isfinite(r.crest_row) else f"{int(r.crest_row)}"
        w(f"  {d:6d} {r.n_cells:+3d}  {r.insert_anchor:5s}  {r.insert_row:6d}  {r.source_rows:>11s}  {crest:>5s}  "
          f"{r.window_mean_m:5.2f} {r.window_min_m:5.2f} {r.window_max_m:6.2f} {r.seam_jump_m:5.2f}  {r['flags']}")
    w("")
    w("ASSUMED: the backbarrier beside the road stands for the lost 1984 ground (which was")
    w("ocean-side beach and dune); it has not accreted or subsided since 1984; alongshore")
    w("structure appears twice; dune array, row 0, road rows and setback untouched.")
    p = insert_scope_step(PRODUCT, "4-fill") / "HAT_fill_copy_scope.txt"
    p.write_text("\n".join(L) + "\n", encoding="utf-8")
    return p


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--domains", default="80,85,5,49,63",
                    help="example domains to draw (comma list)")
    args = ap.parse_args()
    if not FOOTPRINT_CSV.is_file():
        raise SystemExit(f"\n{FOOTPRINT_CSV} not found - run HAT_footprint_1984.py first\n")
    tab = pd.read_csv(FOOTPRINT_CSV).set_index("domain")
    topo_dir, dune_dir, topo_name = topo_dirs(PRODUCT)
    print(f"topography : {PRODUCT} / {topo_name}")
    aud, arrays = plan_fill(tab, topo_dir)
    aud.to_csv(insert_scope_step(PRODUCT, "4-fill") / "fill_copy_by_domain.csv")
    print(f"wrote {SCOPE_DIR / 'fill_copy_by_domain.csv'}")
    figs = []
    for d in [int(x) for x in args.domains.split(",") if x.strip()]:
        if d not in arrays:
            print(f"  GIS {d}: no rows added or removed - nothing to draw")
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
