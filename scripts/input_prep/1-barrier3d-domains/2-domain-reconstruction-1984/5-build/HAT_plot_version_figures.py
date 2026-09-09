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
from hat_topo_version import array_name, dune_topo_root, require_version, year_for_product  # noqa: E402
from hat_figure_style import apply_style, C, elevation_cmap, spines_for_image  # noqa: E402

PRODUCT = "1984-start"
CELL_M = 10.0
DUNE_ROWS = 2
ROAD_ROWS = 2
BERM_EL_M = 1.7                 # BermEl, m MHW: the dune rows are drawn at berm + dune height
ROAD_OFFSET_SCRIPT = (REPO / "scripts" / "input_prep" / "4-mgmt-forcings" / "road_offset"
                      / "1-produce" / "HAT_road_offset_from_dune_start.py")
C_SRC_ROAD, C_ROAD = "0.35", C["ROAD"]
C_ADD, C_REM = "#d62728", "#1f77b4"        # the RdBu pair of the reconstruction figures


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

def draw_grid(ax, topo, dune, road_row, nrows, title, ylabel=True):
    cmap, norm, _ = elevation_cmap()
    n_along = topo.shape[1]
    strip = np.tile(BERM_EL_M + dune[None, :n_along], (DUNE_ROWS, 1))
    grid = np.full((nrows + DUNE_ROWS, n_along), np.nan)
    grid[:DUNE_ROWS] = strip
    grid[DUNE_ROWS:DUNE_ROWS + topo.shape[0]] = topo[:nrows]
    ax.imshow(np.ma.masked_invalid(grid), cmap=cmap, norm=norm, aspect="auto",
              interpolation="nearest", origin="upper", zorder=1)
    ax.axhline(DUNE_ROWS - 0.5, color="#333333", lw=1.0, zorder=4)
    ax.text(n_along - 0.8, DUNE_ROWS / 2 - 0.5, "dune rows", fontsize=7, va="center", ha="right", zorder=6,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0))
    if road_row is not None:
        rs = road_row + DUNE_ROWS
        ax.add_patch(Rectangle((-0.5, rs - 0.5), n_along, ROAD_ROWS, fill=False, ec=C_ROAD, lw=1.5, zorder=6))
        ax.text(0.8, rs + ROAD_ROWS / 2 - 0.5, "NC-12", fontsize=7, ha="left", va="center", color=C_ROAD,
                zorder=7, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0))
    ax.set_xlim(-0.5, n_along - 0.5)
    ax.set_ylim(nrows + DUNE_ROWS - 0.5, -0.5)
    ax.set_xticks([0, 10, 20, 30, 40, 49])
    ax.set_xlabel("alongshore cell")
    if ylabel:
        ax.set_ylabel("cross-shore row (0 = interior row 0, behind the dune)")
    ax.set_title(title, loc="left")
    spines_for_image(ax)


def fig_grid(d: int, src: Version, ver: Version, out_dir: Path) -> Path:
    t0, d0 = src.arrays(d)
    t1, d1 = ver.arrays(d)
    a = ver.audit.get(d)
    n = int(a["n_cells"]) if a else 0
    ins = int(a["insert_row"]) if a and n else None
    nrows = max(t0.shape[0], t1.shape[0])
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11.0, 4.2 + nrows * 0.028), sharey=True,
                                   constrained_layout=True)
    draw_grid(ax0, t0, d0, src.road_row(d), nrows, f"(a) GIS {d}: {src.name}, the source ({t0.shape[0]} rows)")
    draw_grid(ax1, t1, d1, ver.road_row(d), nrows,
              f"(b) GIS {d}: {ver.name} ({t1.shape[0]} rows, {n:+d})" if n else f"(b) GIS {d}: {ver.name} (unchanged)",
              ylabel=False)
    n_along = t0.shape[1]
    if n > 0:
        ax1.add_patch(Rectangle((-0.5, ins + DUNE_ROWS - 0.5), n_along, n, fill=False, ec=C_ADD, lw=1.8, zorder=5))
        ax1.text(n_along * 0.5, ins + DUNE_ROWS + n / 2 - 0.5, f"+{n} rows inserted at row {ins}: {a['operation']}",
                 fontsize=7, ha="center", va="center", color=C_ADD, zorder=6,
                 bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0))
        # the same rows in the source: the window the copy was taken from
        ax0.add_patch(Rectangle((-0.5, ins + DUNE_ROWS - 0.5), n_along, n, fill=False, ec=C_ADD, lw=1.2,
                                ls=(0, (3, 2)), zorder=5))
    elif n < 0:
        ax0.add_patch(Rectangle((-0.5, ins + DUNE_ROWS - 0.5), n_along, -n, facecolor=C_REM, alpha=0.3,
                                ec=C_REM, hatch="////", lw=1.2, zorder=5))
        ax0.text(n_along * 0.5, ins + DUNE_ROWS - n / 2 - 0.5, f"{a['operation']}", fontsize=7, ha="center",
                 va="center", color=C_REM, zorder=6, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0))
        ax1.axhline(ins + DUNE_ROWS - 0.5, xmax=0.6, color=C_REM, lw=2.0, ls=(0, (3, 1.5)), zorder=8)
        ax1.text(n_along * 0.62, ins + DUNE_ROWS - 0.5, f"seam: {-n} rows removed", fontsize=7, ha="left",
                 va="center", color=C_REM, zorder=8, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0))
    cmap, norm, bounds = elevation_cmap()
    labels = ["below 0 (water)"] + [f"{lo:g}–{hi:g}" for lo, hi in zip(bounds[1:-2], bounds[2:-1])] + [f"above {bounds[-2]:g}"]
    handles = [Patch(facecolor=cmap(i), edgecolor="0.4", lw=0.4, label=lab) for i, lab in enumerate(labels)]
    handles += [Patch(facecolor="none", edgecolor=C_ROAD, lw=1.5, label="NC-12 rows at the version's setback")]
    if n > 0:
        handles += [Patch(facecolor="none", edgecolor=C_ADD, lw=1.8, label="inserted block (a copy of the rows below it)"),
                    Patch(facecolor="none", edgecolor=C_ADD, lw=1.2, ls=(0, (3, 2)), label="the source rows copied")]
    elif n < 0:
        handles += [Patch(facecolor=C_REM, alpha=0.3, edgecolor=C_REM, hatch="////", label="rows removed"),
                    Line2D([0], [0], color=C_REM, lw=2.0, ls=(0, (3, 1.5)), label="seam")]
    fig.legend(handles=handles, loc="outside lower center", ncol=6, fontsize=7.5,
               title="elevation classes (m MHW); dune rows drawn at berm + dune height", title_fontsize=7.5)
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"domain_{d:03d}_grid_{ver.name}.png"
    fig.savefig(p, dpi=170, facecolor="white")
    plt.close(fig)
    return p


# =============================================================================
# THE SUMMARY PAGE
# =============================================================================

def _bands(ax, sections):
    for k, ((lo, hi), label) in enumerate(sections):
        if k % 2:
            ax.axvspan(lo - 0.5, hi + 0.5, color="0.93", lw=0, zorder=0)
        ax.text((lo + hi) / 2, 0.98, label, transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=7, color="0.4")


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
    fig, axes = plt.subplots(4, 1, figsize=(13.5, 10.5), sharex=True, constrained_layout=True)
    a0, a1, a2, a3 = axes
    for ax in axes:
        _bands(ax, sections)
    a0.bar(ids, dn, width=0.8, color=np.where(dn > 0, C_ADD, C_REM), zorder=3)
    a0.axhline(0, color="0.2", lw=0.8)
    a0.set_ylabel("rows added / removed")
    a0.set_title(f"(a) the footprint: interior rows in {ver.name} minus {src.name}  "
                 f"({int((dn > 0).sum())} domains +{int(dn[dn > 0].sum())}, {int((dn < 0).sum())} domains {int(dn[dn < 0].sum())})", loc="left")
    a1.plot(ids, rows0, "o-", ms=3, lw=0.8, color="0.5", label=src.name, zorder=3)
    a1.plot(ids, rows1, "o-", ms=3, lw=0.8, color=C["ACCENT"], label=ver.name, zorder=4)
    a1.set_ylabel("interior rows")
    a1.set_title("(b) interior extent per domain (rows of 10 m)", loc="left")
    a1.legend(loc="upper left", ncol=2)
    a2.plot(ids, sb0, "o-", ms=3, lw=0.8, color="0.5", label=f"{src.name} (as measured on the surface)", zorder=3)
    a2.plot(ids, sb1, "o-", ms=3, lw=0.8, color=C["ACCENT"], label=f"{ver.name} (the model input)", zorder=4)
    a2.set_ylabel("NC-12 setback (m from row 0)")
    a2.set_title("(c) the road setback the model receives", loc="left")
    a2.legend(loc="upper left", ncol=2)
    a3.plot(ids, z0, "o-", ms=3, lw=0.8, color="0.5", label=f"mean interior elevation, {src.name}", zorder=3)
    a3.plot(ids, z1, "o-", ms=3, lw=0.8, color=C["ACCENT"], label=f"mean interior elevation, {ver.name}", zorder=4)
    a3.plot(ids, h0 + BERM_EL_M, "s-", ms=3, lw=0.8, color="0.3", label="dune crest (berm + height), both versions", zorder=3)
    a3.set_ylabel("m MHW")
    a3.set_title("(d) mean interior elevation (land cells) and dune crest; the dune array is unchanged", loc="left")
    a3.legend(loc="upper left", ncol=3)
    a3.set_xlabel("GIS domain (south at left)")
    a3.set_xlim(0, ids.max() + 1)
    a3.set_xticks(range(5, int(ids.max()) + 1, 5))
    fig.savefig(out, dpi=170, facecolor="white")
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
    fig_w = 20.0
    fig_h = min(max(4.5, fig_w * (n_rows / n_cols) * 1.8), 7.5)
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")
    ax = fig.add_axes([0.06, 0.18, 0.88, 0.68])
    ax.set_facecolor(ext.ISLAND_OCEAN_COLOR)
    im = ax.pcolormesh(np.ma.masked_invalid(canvas), cmap=cmap, norm=norm, shading="auto", rasterized=True)
    ax.pcolormesh(np.ma.masked_where(~road_canvas, road_canvas.astype(float)), cmap=ListedColormap([C_ROAD]),
                  vmin=0.0, vmax=1.0, shading="auto", rasterized=True, zorder=4)
    ax.plot([], [], color=C_ROAD, lw=3, label=f"NC-12 as the model places it ({ver.name} setback)")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    cax = fig.add_axes([0.955, 0.18, 0.013, 0.68])
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label("Elevation (m MHW)", fontsize=12, color="#1a1a2e", labelpad=10, rotation=270)
    cbar.set_ticks([-1, 0, 1, 2, 3, 4])
    ticks, labels = [], []
    for k, d in enumerate(use):
        if d % 5 == 0 or d == 1:
            ticks.append(starts[k] + grids[k].shape[1] // 2)
            labels.append(str(d))
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_xlabel("Domain (S → N,  Cape Hatteras to Rodanthe)", fontsize=12, labelpad=8)
    ax.set_ylabel("Cross-shore cell (raw_offset frame)", fontsize=12)
    for k, d in enumerate(use):
        if d % 10 == 0:
            ax.axvline(starts[k] - 0.5, color="#aaaaaa", lw=0.4, alpha=0.5, zorder=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    what = "dune + interior" if ext.ISLAND_INCLUDE_DUNE else "interior"
    note = (f"{what}, every domain padded to {pad} cells / {pad * CELL_M:.0f} m cross-shore" if mode == "padded"
            else f"{what}, each domain trimmed to its own island")
    ax.set_title(f"Hatteras Island — CASCADE Initialization  |  {year} offsets  |  {PRODUCT} {ver.name} "
                 f"(built from {source_from_manifest(ver.dir) or '?'}) {note}  ({len(use)} domains)",
                 fontsize=13, fontweight="bold", color="#1a1a2e", pad=12)
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


# =============================================================================

def write_readme(ver: Version, src: Version, n_grid: int, figs: list[Path]) -> Path:
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

| figure | what |
|---|---|
| `grid/domain_NNN_grid_{ver.name}.png` ({n_grid}) | one per domain: `{src.name}` beside `{ver.name}` as the model holds them — the two dune rows on top (berm + dune height), every interior row down the page, elevation classes (m MHW), NC-12's two rows at each version's setback, and the footprint: the inserted block outlined in red (add; dashed in the source: the rows it copies) or the removed rows hatched blue in the source and the seam marked in the version (remove). Unchanged domains have two identical panels. |
| `../HAT_dune_topo_summary_{ver.name}.png` | every domain on one page: rows added or removed, interior rows, the road setback the model receives, mean interior elevation and the dune crest, `{src.name}` against `{ver.name}`, communities banded. The counterpart of the extractor's summary page. |
| `../HAT_dune_topo_island_planview_{ver.name}_<year>_{{trimmed,padded}}.png` | the island in plan view at the period's dune offsets, in the extractor's poster style, NC-12 drawn where the **model** places it (the version's setback) rather than from the GIS mask, whose frame a built interior no longer shares. |

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
    write_readme(ver, src, n_grid if not a.domains else len(list((ver.dir / "figures" / "grid").glob("*.png"))), figs)
    print(f"done: {ver.dir}")


if __name__ == "__main__":
    main()
