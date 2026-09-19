#!/usr/bin/env python3
r"""
HAT_plot_version_figures.py
==============================================================================
The figures a BUILT dune-topo version gets, so that v3 sits beside v1 and v2
with a figure set of its own (Hannah, 2026-09-09: "v3 should have figures
here as well").

WHY NOT THE EXTRACTOR'S FIGURES. v1 and v2 carry `figures/qc/` and
`figures/gis_vs_processed/`: the raw DEM profile against the extracted one,
per domain, and the island summary of windows and dune heights. Those are
figures OF AN EXTRACTION - they need the picks, the raw profiles and the
straightening frame. A built version is not extracted: it is its source
version with rows inserted or removed and the road setback re-set, so the
questions its figures answer are different - what does each domain look like
now, beside what it looked like before; where did the rows go; what changed
island-wide. Nothing here re-measures anything: every number is read from the
version's own arrays, its setback CSV and its footprint audit.

WHAT IS WRITTEN, into dune-topo/<version>/
    figures/grid/domain_NNN_grid_<version>.png    one per domain: the source
        version and this version side by side as the model holds them - the
        two dune rows on top (drawn at berm + dune height), every interior row
        down the page, elevation classes (m MHW), NC-12's two rows at each
        version's setback, and the footprint: the inserted block outlined
        (add) or the removed rows hatched in the source and the seam marked
        in the version (remove). Unchanged domains are drawn too, so the set
        is complete; their two panels are identical.
    HAT_dune_topo_summary_<version>.png           every domain on one page:
        interior rows, the road setback, mean interior elevation and mean
        dune height, source against version, with the communities banded.
        The counterpart of the extractor's summary page.
    HAT_dune_topo_island_planview_<version>_<year>_{trimmed,padded}.png
        the island in plan view at the period's dune offsets, in the
        extractor's poster style, with NC-12 drawn where the MODEL places it
        (the version's setback), not from the GIS mask - a built version's
        interior frame is no longer the mask's frame. The counterpart of the
        extractor's plan views.
    figures/README.md                             what these are and are not

USAGE
    python HAT_plot_version_figures.py                      # v3, source from its manifest
    python HAT_plot_version_figures.py --version v3 --source v2
    python HAT_plot_version_figures.py --domains 85,63      # only the grid panels of these
    python HAT_plot_version_figures.py --no-grid            # summary and plan views only
==============================================================================
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer.hat_topo_version import array_name, dune_topo_root, require_version, year_for_product  # noqa: E402
from site_layer.hat_figure_style import (apply_style, C, C_1984, C_1997, elevation_cmap, spines_for_image,  # noqa: E402
                              figsize, save, DOMAIN_AXIS_LABEL, town_bands, open_frame, _title)

PRODUCT = "1984-start"
CELL_M = 10.0
DUNE_ROWS = 2
ROAD_ROWS = 2
BERM_EL_M = 1.7                 # BermEl, m MHW: the dune rows are drawn at berm + dune height
ROAD_OFFSET_SCRIPT = (REPO / "scripts" / "input_prep" / "4-mgmt-forcings" / "road_offset"
                      / "1-produce" / "HAT_road_offset_from_dune_start.py")
C_SRC_ROAD, C_ROAD = C["BASE"], C["ROAD"]
C_ADD, C_REM = C_1984, C_1997               # the RdBu pair of the reconstruction figures: red added, blue removed
# What the two versions are called on the figures (no working vocabulary; the
# version names themselves are in the file names and the README).
SRC_LABEL = "as extracted (1996 surface)"
VER_LABEL = "1984 reconstruction"


# =============================================================================
# READING A VERSION
# =============================================================================

def read_setback_csv(p: Path) -> dict[int, float]:
    """The two-row model-facing CSV: domain ids, then metres landward of row 0."""
    if not p.is_file():
        return {}
    rows = list(csv.reader(open(p, newline="")))
    ids = [int(float(x)) for x in rows[0] if x.strip()]
    vals = [float(x) for x in rows[1] if x.strip()]
    return dict(zip(ids, vals))


def read_audit(p: Path) -> dict[int, dict]:
    if not p.is_file():
        return {}
    return {int(r["domain"]): r for r in csv.DictReader(open(p, newline=""))}


def source_from_manifest(vdir: Path) -> str | None:
    m = vdir / "RUN_MANIFEST.txt"
    if not m.is_file():
        return None
    for line in m.read_text(encoding="utf-8").splitlines():
        hit = re.match(r"\s*source\s*:\s*\S+/(v\d+)", line)
        if hit:
            return hit.group(1)
    return None


class Version:
    def __init__(self, name: str):
        self.name = name
        self.dir = require_version(PRODUCT, name, "the version to draw")
        self.setback = read_setback_csv(self.dir / "RoadSetback_1984_dunestart.csv")
        self.audit = read_audit(self.dir / "HAT_footprint_audit.csv")

    def arrays(self, d: int) -> tuple[np.ndarray, np.ndarray]:
        """(interior m MHW, dune height m) for one domain."""
        topo = np.load(self.dir / "topography" / array_name("topography", d)) * CELL_M
        dune = np.load(self.dir / "dunes" / array_name("dune", d)) * CELL_M
        return topo, dune

    def road_row(self, d: int) -> int | None:
        sb = self.setback.get(d)
        return None if sb is None or not np.isfinite(sb) else int(sb // CELL_M)


# =============================================================================
# THE GRID, ONE DOMAIN, SOURCE BESIDE VERSION
# =============================================================================

def draw_grid(ax, topo, dune, road_row, nrows, k, title, ylabel=True):
    cmap, norm, _ = elevation_cmap()
    n_along = topo.shape[1]
    strip = np.tile(BERM_EL_M + dune[None, :n_along], (DUNE_ROWS, 1))
    grid = np.full((nrows + DUNE_ROWS, n_along), np.nan)
    grid[:DUNE_ROWS] = strip
    grid[DUNE_ROWS:DUNE_ROWS + topo.shape[0]] = topo[:nrows]
    ax.imshow(np.ma.masked_invalid(grid), cmap=cmap, norm=norm, aspect="auto",
              interpolation="nearest", origin="upper", zorder=1)
    ax.axhline(DUNE_ROWS - 0.5, color=C["INK"], lw=0.8, zorder=4)
    ax.text(n_along - 0.8, DUNE_ROWS / 2 - 0.5, "dune", fontsize=7, va="center", ha="right", zorder=6,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0))
    if road_row is not None:
        rs = road_row + DUNE_ROWS
        ax.add_patch(Rectangle((-0.5, rs - 0.5), n_along, ROAD_ROWS, fill=False, ec=C_ROAD, lw=1.2, zorder=6))
        ax.text(0.8, rs + ROAD_ROWS / 2 - 0.5, "NC-12", fontsize=7, ha="left", va="center", color=C_ROAD,
                zorder=7, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0))
    ax.set_xlim(-0.5, n_along - 0.5)
    ax.set_ylim(nrows + DUNE_ROWS - 0.5, -0.5)
    ax.set_xticks([0, 25, 49])
    ax.set_xlabel("alongshore cell")
    if ylabel:
        ax.set_ylabel("cross-shore row (0 = interior row 0)")
    _title(ax, k, title)
    spines_for_image(ax)


def fig_grid(d: int, src: Version, ver: Version, out_dir: Path) -> Path:
    t0, d0 = src.arrays(d)
    t1, d1 = ver.arrays(d)
    a = ver.audit.get(d)
    n = int(a["n_cells"]) if a else 0
    ins = int(a["insert_row"]) if a and n else None
    nrows = max(t0.shape[0], t1.shape[0])
    # a single-column figure per domain: the two panels side by side, the
    # height following the row count so a deep domain is not squashed
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=figsize("single", height=3.3 + nrows * 0.016), sharey=True,
                                   constrained_layout=True)
    draw_grid(ax0, t0, d0, src.road_row(d), nrows, 0, f"{t0.shape[0]} rows")
    draw_grid(ax1, t1, d1, ver.road_row(d), nrows, 1,
              f"{t1.shape[0]} rows ({'+' if n > 0 else '−'}{abs(n)})" if n else f"{t1.shape[0]} rows", ylabel=False)
    n_along = t0.shape[1]
    if n > 0:
        ax1.add_patch(Rectangle((-0.5, ins + DUNE_ROWS - 0.5), n_along, n, fill=False, ec=C_ADD, lw=1.4, zorder=5))
        ax1.text(n_along * 0.5, ins + DUNE_ROWS + n / 2 - 0.5, f"+{n} row{'s' if n != 1 else ''}",
                 fontsize=7, ha="center", va="center", color=C_ADD, zorder=6,
                 bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0))
        # the same rows in the source: the window the copy was taken from
        ax0.add_patch(Rectangle((-0.5, ins + DUNE_ROWS - 0.5), n_along, n, fill=False, ec=C_ADD, lw=1.0,
                                ls=(0, (3, 2)), zorder=5))
    elif n < 0:
        ax0.add_patch(Rectangle((-0.5, ins + DUNE_ROWS - 0.5), n_along, -n, facecolor=C_REM, alpha=0.3,
                                ec=C_REM, hatch="////", lw=1.0, zorder=5))
        ax0.text(n_along * 0.5, ins + DUNE_ROWS - n / 2 - 0.5, f"−{-n} row{'s' if n != -1 else ''}",
                 fontsize=7, ha="center", va="center", color=C_REM, zorder=6,
                 bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0))
        # the seam on the right half of the panel, clear of the NC-12 label at the left
        ax1.axhline(ins + DUNE_ROWS - 0.5, xmin=0.55, color=C_REM, lw=1.6, ls=(0, (3, 1.5)), zorder=8)
        ax1.text(n_along * 0.53, ins + DUNE_ROWS - 0.5, "seam", fontsize=7, ha="right",
                 va="center", color=C_REM, zorder=8, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0))
    cmap, norm, bounds = elevation_cmap()
    # one legend under the figure, two columns: the elevation classes (m MHW)
    # and the footprint marks; a single column is too narrow for two legends
    labels = ["< 0 m (water)"] + [f"{lo:g}–{hi:g} m" for lo, hi in zip(bounds[1:-2], bounds[2:-1])] + [f"> {bounds[-2]:g} m"]
    handles = [Patch(facecolor=cmap(i), edgecolor="0.4", lw=0.4, label=lab) for i, lab in enumerate(labels)]
    handles += [Patch(facecolor="none", edgecolor=C_ROAD, lw=1.2, label="NC-12 (two rows)")]
    if n > 0:
        handles += [Patch(facecolor="none", edgecolor=C_ADD, lw=1.4, label="rows inserted behind NC-12"),
                    Patch(facecolor="none", edgecolor=C_ADD, lw=1.0, ls=(0, (3, 2)), label="rows copied into them")]
    elif n < 0:
        handles += [Patch(facecolor=C_REM, alpha=0.3, edgecolor=C_REM, hatch="////", label="rows removed before NC-12"),
                    Line2D([0], [0], color=C_REM, lw=1.6, ls=(0, (3, 1.5)), label="seam left by the removal")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False, fontsize=7.5,
               columnspacing=1.0, handlelength=1.4)
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"domain_{d:03d}_grid_{ver.name}.png"
    save(fig, p, vector=False, facecolor="white")
    plt.close(fig)
    return p


# =============================================================================
# THE SUMMARY PAGE
# =============================================================================

def summary_counts(src: Version, ver: Version) -> dict:
    """The footprint's totals, for the README (they were in a panel title)."""
    ids = sorted(int(p.stem.split("_")[1]) for p in (ver.dir / "topography").glob("domain_*_topography.npy"))
    dn = np.array([ver.arrays(d)[0].shape[0] - src.arrays(d)[0].shape[0] for d in ids])
    return dict(n_add=int((dn > 0).sum()), rows_add=int(dn[dn > 0].sum()),
                n_rem=int((dn < 0).sum()), rows_rem=int(-dn[dn < 0].sum()), n_same=int((dn == 0).sum()))


def fig_summary(src: Version, ver: Version, sections, out: Path) -> Path:
    ids = sorted(int(p.stem.split("_")[1]) for p in (ver.dir / "topography").glob("domain_*_topography.npy"))
    rows0, rows1, sb0, sb1, z0, z1, h0, h1 = ([] for _ in range(8))
    for d in ids:
        t0, d0 = src.arrays(d)
        t1, d1 = ver.arrays(d)
        rows0.append(t0.shape[0]); rows1.append(t1.shape[0])
        sb0.append(src.setback.get(d, np.nan)); sb1.append(ver.setback.get(d, np.nan))
        land0, land1 = t0[t0 > -2.99], t1[t1 > -2.99]
        z0.append(float(land0.mean()) if land0.size else np.nan); z1.append(float(land1.mean()) if land1.size else np.nan)
        h0.append(float(np.mean(d0))); h1.append(float(np.mean(d1)))
    ids = np.array(ids)
    rows0, rows1, sb0, sb1, z0, z1, h0, h1 = map(np.array, (rows0, rows1, sb0, sb1, z0, z1, h0, h1))
    dn = rows1 - rows0

    apply_style()
    fig, axes = plt.subplots(4, 1, figsize=figsize("double", height=8.2), sharex=True, constrained_layout=True)
    a0, a1, a2, a3 = axes
    for k, ax in enumerate(axes):
        town_bands(ax, label=(k == 0))
        open_frame(ax)
        ax.grid(axis="y")
        ax.set_axisbelow(True)
    a0.bar(ids, dn, width=0.8, color=np.where(dn > 0, C_ADD, C_REM), zorder=3)
    a0.axhline(0, color=C["INK"], lw=0.6)
    a0.set_ylabel("rows added (+)\n/ removed (−)")
    _title(a0, 0, "interior rows added or removed")
    a1.plot(ids, rows0, "o-", ms=2.5, lw=0.8, color=C["BASE"], zorder=3)
    a1.plot(ids, rows1, "o-", ms=2.5, lw=0.8, color=C["ACCENT"], zorder=4)
    a1.set_ylabel("interior rows")
    _title(a1, 1, "interior extent (rows of 10 m)")
    a2.plot(ids, sb0, "o-", ms=2.5, lw=0.8, color=C["BASE"], zorder=3)
    a2.plot(ids, sb1, "o-", ms=2.5, lw=0.8, color=C["ACCENT"], zorder=4)
    a2.set_ylabel("NC-12 setback\n(m from row 0)")
    _title(a2, 2, "road setback the model receives")
    a3.plot(ids, z0, "o-", ms=2.5, lw=0.8, color=C["BASE"], zorder=3)
    a3.plot(ids, z1, "o-", ms=2.5, lw=0.8, color=C["ACCENT"], zorder=4)
    a3.plot(ids, h0 + BERM_EL_M, "s-", ms=2.5, lw=0.8, color=C["REF"], zorder=3)
    a3.set_ylabel("m MHW")
    _title(a3, 3, "mean interior elevation (land cells) and dune crest")
    a3.set_xlabel(DOMAIN_AXIS_LABEL)
    a3.set_xlim(0, ids.max() + 1)
    a3.set_xticks([1] + list(range(10, int(ids.max()) + 1, 10)))
    a3.set_xticks(ids, minor=True)
    handles = [Patch(facecolor=C_ADD, label="rows added"), Patch(facecolor=C_REM, label="rows removed"),
               Line2D([0], [0], marker="o", ms=2.5, lw=0.8, color=C["BASE"], label=SRC_LABEL),
               Line2D([0], [0], marker="o", ms=2.5, lw=0.8, color=C["ACCENT"], label=VER_LABEL),
               Line2D([0], [0], marker="s", ms=2.5, lw=0.8, color=C["REF"],
                      label="dune crest (berm + dune height), unchanged")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    save(fig, out, facecolor="white")
    plt.close(fig)
    return out


# =============================================================================
# THE ISLAND PLAN VIEW, IN THE EXTRACTOR'S STYLE
# =============================================================================

def load_extractor():
    spec = importlib.util.spec_from_file_location("hat_road_offset", ROAD_OFFSET_SCRIPT)
    ro = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ro)
    return ro.load_extractor(PRODUCT)


def fig_planview(ver: Version, ext, offsets: dict, year: int, mode: str, out: Path) -> Path | None:
    dom, off = offsets[year]
    omap = {int(a): float(b) for a, b in zip(dom, off)}
    ids = sorted(int(p.stem.split("_")[1]) for p in (ver.dir / "topography").glob("domain_*_topography.npy"))
    use = [d for d in ids if d in omap]
    if not use:
        return None
    grids, dunes, roads, off_cells = [], [], [], []
    pad = int(ext.ISLAND_PAD_ROWS)
    sentinel = float(ext.SENTINEL_WATER_M)
    for d in use:
        g, dh = ver.arrays(d)
        r = np.zeros(g.shape, dtype=bool)
        rr = ver.road_row(d)
        if rr is not None and rr < g.shape[0]:
            r[rr:rr + ROAD_ROWS] = True
        if mode == "padded":
            n = g.shape[0]
            if n < pad:
                g = np.vstack([g, np.full((pad - n, g.shape[1]), sentinel)])
                r = np.vstack([r, np.zeros((pad - n, r.shape[1]), dtype=bool)])
            else:
                g, r = g[:pad], r[:pad]
        grids.append(g)
        roads.append(r)
        dunes.append(dh + (float(ext.BERM_ELEV_NAVD_M) - float(ext.MHW_M)))
        off_cells.append(int(round(omap[d] / CELL_M)))
    max_rows = max(g.shape[0] for g in grids)
    n_rows = max(off_cells) + max_rows + 5
    n_cols = sum(g.shape[1] for g in grids)
    canvas = np.full((n_rows, n_cols), np.nan)
    road_canvas = np.zeros((n_rows, n_cols), dtype=bool)
    col, starts = 0, []
    for g, dh, r, oc in zip(grids, dunes, roads, off_cells):
        h, w = g.shape
        canvas[oc:oc + h, col:col + w] = g
        road_canvas[oc:oc + h, col:col + w] = r
        if ext.ISLAND_INCLUDE_DUNE and oc >= 1:
            canvas[oc - 1, col:col + w] = dh[:w]
        starts.append(col)
        col += w
    cmap, norm = ext._island_norm()
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.34), constrained_layout=True)
    ax.set_facecolor(ext.ISLAND_OCEAN_COLOR)
    # the villages as light bands, in the canvas's column frame: a translucent
    # white so they read on the ocean colour and vanish under the island
    col_of = {d: (starts[k], starts[k] + grids[k].shape[1]) for k, d in enumerate(use)}
    spans = {}
    try:
        from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS
        for name, (lo, hi) in HATTERAS_ANNOTATIONS.town_spans.items():
            inside = [d for d in use if lo <= d <= hi]
            if inside:
                spans[name] = (col_of[min(inside)][0] + 0.5, col_of[max(inside)][1] - 0.5)
    except ImportError:
        pass
    town_bands(ax, spans=spans, shade=(1.0, 1.0, 1.0, 0.35))
    # the extractor's cmap paints masked cells in the ocean colour, which would
    # cover the bands; here the axes background is the ocean and the mask is clear
    cmap = cmap.copy()
    cmap.set_bad((0.0, 0.0, 0.0, 0.0))
    im = ax.pcolormesh(np.ma.masked_invalid(canvas), cmap=cmap, norm=norm, shading="auto", rasterized=True)
    ax.pcolormesh(np.ma.masked_where(~road_canvas, road_canvas.astype(float)), cmap=ListedColormap([C_ROAD]),
                  vmin=0.0, vmax=1.0, shading="auto", rasterized=True, zorder=4)
    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.015)
    cbar.set_label("elevation (m MHW)")
    cbar.set_ticks([-1, 0, 1, 2, 3, 4])
    ticks, labels = [], []
    for k, d in enumerate(use):
        if d % 10 == 0 or d == 1:
            ticks.append(starts[k] + grids[k].shape[1] // 2)
            labels.append(str(d))
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("cross-shore cell")
    for k, d in enumerate(use):
        if d % 10 == 0:
            ax.axvline(starts[k] - 0.5, color="#aaaaaa", lw=0.4, alpha=0.5, zorder=2)
    open_frame(ax)
    fig.legend(handles=[Line2D([0], [0], color=C_ROAD, lw=3, label="NC-12 as the model places it (1984 setback)")],
               loc="outside lower center", frameon=False)
    save(fig, out, vector=False, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


# =============================================================================

def write_readme(ver: Version, src: Version, n_grid: int, figs: list[Path], counts: dict) -> Path:
    p = ver.dir / "figures" / "README.md"
    p.parent.mkdir(parents=True, exist_ok=True)
    rel = [str(f.relative_to(ver.dir)).replace("\\", "/") for f in figs]
    p.write_text(f'''# `{ver.name}/figures` — the figures of a built version

Written {datetime.now():%Y-%m-%d %H:%M} by `HAT_plot_version_figures.py`
(`scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/5-build/`).

`{ver.name}` was **built** from `{src.name}`, not extracted, so it does not get the
extractor's `qc/` and `gis_vs_processed/` sets (those compare the raw DEM
profile with the extraction and need the picks and the straightening frame).
It gets the counterparts a built version can answer, every number read from
its own arrays, its `RoadSetback_1984_dunestart.csv` and
`HAT_footprint_audit.csv`; nothing is re-measured.

On the figures `{src.name}` is called "{SRC_LABEL}" and `{ver.name}`
"{VER_LABEL}"; the images carry no titles or statistics, so the captions
below do. Every PNG has a PDF beside it except the raster-only grid panels.

| figure | what |
|---|---|
| `grid/domain_NNN_grid_{ver.name}.png` ({n_grid}) | one per domain, a single-column figure: (a) `{src.name}` beside (b) `{ver.name}` as the model holds them — the two dune rows on top, drawn at berm + dune height, every interior row down the page, elevation classes (m MHW), NC-12's two rows at each version's setback, and the footprint: the inserted block outlined in red (dashed in the source: the rows it copies) or the removed rows hatched blue in the source and the seam marked in the version. The panel titles give the row count and the change. Unchanged domains have two identical panels. |
| `../HAT_dune_topo_summary_{ver.name}.png` | every domain on one page, {src.name} (grey) against {ver.name} (purple), the villages banded: (a) interior rows added (red) or removed (blue) — {counts["n_add"]} domains gain {counts["rows_add"]} rows, {counts["n_rem"]} lose {counts["rows_rem"]}, {counts["n_same"]} are unchanged; (b) interior rows per domain; (c) the NC-12 setback, as measured on the 1996 surface and as the model receives it for 1984; (d) mean interior elevation over land cells, with the dune crest (berm + dune height, green), which the build does not change. Domain 1 is at Cape Point, 90 at Pea Island. The counterpart of the extractor's summary page. |
| `../HAT_dune_topo_island_planview_{ver.name}_<year>_{{trimmed,padded}}.png` | the island in plan view at the period's dune offsets (dune row plus interior; `padded` pads every domain to the extractor's cross-shore length, `trimmed` keeps each domain's own), NC-12 drawn where the **model** places it (the version's setback) rather than from the GIS mask, whose frame a built interior no longer shares. Villages as light bands over the water. |

Files:

{chr(10).join("- `" + r + "`" for r in rel[:3])}
- `figures/grid/` — {n_grid} domain panels
''', encoding="utf-8")
    return p


def main() -> None:
    ap = argparse.ArgumentParser(description="figures for a built dune-topo version")
    ap.add_argument("--version", default="v3")
    ap.add_argument("--source", default="", help="the version it was built from (default: from RUN_MANIFEST)")
    ap.add_argument("--domains", default="", help="comma-separated GIS ids for the grid panels (default: all)")
    ap.add_argument("--no-grid", action="store_true")
    ap.add_argument("--no-planview", action="store_true")
    a = ap.parse_args()
    apply_style()
    ver = Version(a.version)
    src_name = a.source or source_from_manifest(ver.dir)
    if not src_name:
        raise SystemExit(f"{ver.name}: no source version in RUN_MANIFEST.txt; pass --source")
    src = Version(src_name)
    print(f"{ver.name} (built from {src.name}): {len(ver.audit)} domains in the footprint audit, "
          f"{len(ver.setback)} setbacks")
    figs = []
    ext = load_extractor()
    year = year_for_product(PRODUCT, strict=False)

    figs.append(fig_summary(src, ver, ext.SECTIONS, ver.dir / f"HAT_dune_topo_summary_{ver.name}.png"))
    print(f"  wrote {figs[-1].name}")
    if not a.no_planview:
        offsets = ext.load_offsets()
        if year in offsets:
            for mode in ext.ISLAND_CROSS_SHORE_MODES:
                p = fig_planview(ver, ext, offsets, year, mode,
                                 ver.dir / f"HAT_dune_topo_island_planview_{ver.name}_{year}_{mode}.png")
                if p:
                    figs.append(p)
                    print(f"  wrote {p.name}")
        else:
            print(f"  no {year} offsets loaded ({sorted(offsets)}); plan views skipped")
    n_grid = 0
    if not a.no_grid:
        ids = ([int(x) for x in a.domains.split(",")] if a.domains else
               sorted(int(p.stem.split("_")[1]) for p in (ver.dir / "topography").glob("domain_*_topography.npy")))
        out_dir = ver.dir / "figures" / "grid"
        for k, d in enumerate(ids, 1):
            fig_grid(d, src, ver, out_dir)
            n_grid += 1
            if k % 15 == 0 or k == len(ids):
                print(f"  grid {k}/{len(ids)}", flush=True)
    write_readme(ver, src, n_grid if not a.domains else len(list((ver.dir / "figures" / "grid").glob("*.png"))), figs,
                 summary_counts(src, ver))
    print(f"done: {ver.dir}")


if __name__ == "__main__":
    main()
