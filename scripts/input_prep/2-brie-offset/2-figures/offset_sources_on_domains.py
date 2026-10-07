"""
The start-year island as the model builds it from each offset source: dune line above, shoreline below.

    python scripts/input_prep/2-brie-offset/2-figures/offset_sources_on_domains.py --year 1996

Places every Barrier3D domain grid of the start topography at its island offset, once
with the dune-line build and once with the shoreline build. One figure per 15 domains,
so the six together cover the island; each has the two sources stacked.
Written to domains/ inside the dune line vs shoreline comparison of that shoreline build.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-06
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from site_layer import hat_topo_version as _tv  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, ELEV_WATER, _title, apply_style, caption, elevation_cmap, figsize, save, town_bands)

# --- CONFIG ------------------------------------------------------------------
# One figure per 15 domains, the same sixths as the line overlay; the six cover GIS 1-90 (Hannah, 2026-10-06)
SECTIONS = tuple((lo, lo + 14) for lo in range(1, 91, 15))
CELL_M = 10.0                # Barrier3D cell, m
DAM = 10.0                   # decametre -> metre
SEAWARD_PAD_M = 250          # open water seaward of the most seaward shoreline, room for the village strip
LANDWARD_PAD_M = 60          # drawn landward of the last land row
PARAMETER_FILE = PROJECT_ROOT / "data" / "hatteras_init" / "Hatteras-CASCADE-parameters.yaml"
ROWS = {
    "duneline": {"label": "Dune-line offset", "colour": C["ACCENT"]},
    "shoreline": {"label": "Shoreline offset", "colour": C["LATE"]},
}
# -----------------------------------------------------------------------------


# Real-domain offsets (GIS 1-90) of one source's build, model frame (zero on its own minimum)
def _offsets(year, source, version=None):
    path = _tv.offset_file(year, "unpadded", source=source, version=version)
    return pd.read_csv(path).iloc[:, 1].to_numpy(float)


# Per-domain distance from the shared offshore datum (GIS 1-90), the frame where the beach width is real
def _datum(year, source, version):
    raw = next(_tv.offset_build_dir(year, version, source).glob("*_raw.csv"))
    return pd.read_csv(raw).groupby("domain_id")["ORIG_LEN"].mean().loc[1:90].to_numpy(float)


# One domain's grid in m MHW, cut to its last land row: the two dune rows, then the interior, ocean first
def _grid(elev_path, dune_path, berm_m):
    interior = np.load(elev_path) * DAM
    dune = np.load(dune_path) * DAM + berm_m
    grid = np.vstack([dune, dune, interior])
    land = np.flatnonzero((grid > 0).any(axis=1))
    return grid[:int(land[-1]) + 1 + int(LANDWARD_PAD_M / CELL_M)]


def _draw(ax, grids, x_s, other_x_s, beach_m, gis_lo, gis_hi, ylim):
    cmap, norm, _ = elevation_cmap()
    beach_col = cmap(norm(np.array([1.2])))[0]
    for g in range(gis_lo, gis_hi + 1):
        grid = np.ma.masked_less_equal(grids[g - 1], 0)
        y0 = x_s[g - 1] + beach_m
        ax.fill_between([g - 0.5, g + 0.5], x_s[g - 1], y0, color=beach_col, lw=0, zorder=1)
        ax.imshow(grid, cmap=cmap, norm=norm, origin="lower", interpolation="nearest",
                  aspect="auto", extent=(g - 0.5, g + 0.5, y0, y0 + len(grid) * CELL_M), zorder=2)
    gis = np.arange(gis_lo, gis_hi + 1)
    # Domain borders, faint and dashed, so each 500 m domain reads as its own column
    for edge in np.arange(gis_lo - 0.5, gis_hi + 1.0):
        ax.axvline(edge, color="white", lw=0.5, ls=(0, (2, 3)), alpha=0.55, zorder=3)
    # Where the shoreline WOULD be with the other source, one dash per domain: white over the water
    # (the other source is further seaward), black over the beach (further landward)
    for g in gis:
        seaward = other_x_s[g - 1] < x_s[g - 1]
        ax.hlines(other_x_s[g - 1], g - 0.45, g + 0.45, color="white" if seaward else "black",
                  lw=1.5, ls=(0, (3, 2)), zorder=5)
    ax.set_facecolor(ELEV_WATER)
    ax.set_xlim(gis_lo - 0.5, gis_hi + 0.5)
    ax.set_ylim(*ylim)
    _domain_ticks(ax, gis)


# What a bar in the shift panel is, top right of the panel
def _shift_note(ax, gap):
    ax.text(0.995, 0.97, f"Shift = {gap:.0f} m − beach width at the domain",
            transform=ax.transAxes, ha="right", va="top", fontsize=7, color=C["INK"],
            bbox=dict(facecolor="white", edgecolor="0.7", lw=0.5, pad=2))


# Red where the shoreline offset moves a domain seaward, blue where landward
def _shift_colours(shift):
    return [C["LATE"] if v > 0 else C["EARLY"] for v in shift]


# GIS numbers every 5 on every panel, a small tick for each domain between
def _domain_ticks(ax, gis):
    ax.set_xticks(gis[(gis % 5) == 0])
    ax.set_xticks(gis, minor=True)
    ax.tick_params(axis="x", which="both", labelbottom=True)
    ax.tick_params(axis="x", which="minor", length=2)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--year", type=int, default=1996)
    ap.add_argument("--duneline-version", default=None, help="default: the dune CURRENT")
    ap.add_argument("--shoreline-version", default=None, help="default: the shoreline CURRENT")
    args = ap.parse_args(argv)
    year = args.year
    apply_style()

    ver = {"duneline": args.duneline_version or _tv.offset_version(year, "duneline"),
           "shoreline": args.shoreline_version or _tv.offset_version(year, "shoreline")}
    off = {s: _offsets(year, s, ver[s]) for s in ROWS}
    diff = off["shoreline"] - off["duneline"]
    pos = off
    # Each build is zeroed on its own most seaward domain, so the shift is a constant minus the beach width:
    # shift = gap - beach, gap = the distance between the two zero points (dune line minus shoreline, datum frame)
    datum = {s: _datum(year, s, ver[s]) for s in ROWS}
    beach = datum["duneline"] - datum["shoreline"]
    gap = datum["duneline"].min() - datum["shoreline"].min()
    zero = {s: int(np.argmin(off[s])) + 1 for s in ROWS}
    if not np.allclose(diff, gap - beach, atol=1e-6):
        sys.exit("shift is not gap - beach width; the offsets and raw files disagree")

    product = _tv.product_for_year(year)
    dune_version = _tv.topo_dirs(product)[2]
    elev, dunes = _tv.domain_arrays(product, first_gis=1, last_gis=90)
    with open(PARAMETER_FILE) as f:
        p = yaml.safe_load(f)
    berm_m = float(p["BermEl"]) - float(p["MHW"])
    # Barrier3D's initial beach width: int(BermEl / beta) in dam, the same for every domain
    beach_m = int(berm_m / DAM / float(p["beta"])) * 10
    grids = [_grid(e, d, berm_m) for e, d in zip(elev, dunes)]
    # A subfolder of the comparison, so the eight figures do not crowd the overlay (Hannah, 2026-10-06)
    out_dir = _tv.offset_source_comparison_dir(year, "duneline_vs_shoreline", ver["shoreline"]) / "domains"
    cmap, norm, bounds = elevation_cmap()
    # Panel titles name the approach, not the build: colleagues need the source, not the version
    w0, w1 = (pd.Timestamp(t) for t in _tv.shoreline_window_for_year(year))
    approach = {
        "duneline": f"Dune line digitised from {_tv.dune_line_for_year(year)} aerial imagery",
        "shoreline": f"Mean satellite shoreline, {w0:%b %Y} – {w1:%b %Y} (CoastSat)",
    }

    for lo, hi in SECTIONS:
        sl = slice(lo - 1, hi)
        top = max(max(pos[s][g - 1] for s in ROWS) + beach_m + len(grids[g - 1]) * CELL_M
                  for g in range(lo, hi + 1))
        ylim = (min(pos[s][sl].min() for s in ROWS) - SEAWARD_PAD_M, top)

        fig, axes = plt.subplots(3, 1, figsize=figsize("double", height=8.4), sharex=True,
                                 gridspec_kw={"height_ratios": [1, 1, 0.38]}, constrained_layout=True)
        axes[1].sharey(axes[0])
        for i, src in enumerate(ROWS):
            other = "shoreline" if src == "duneline" else "duneline"
            ax = axes[i]
            _draw(ax, grids, pos[src], pos[other], beach_m, lo, hi, ylim)
            _title(ax, i, approach[src])
            ax.set_ylabel("Cross-shore position (m)")
            town_bands(ax, where="bottom", strip=0.05)
        d = diff[sl]
        # (c) the same gap as a number per domain, in the house red/blue pair: red seaward, blue landward
        gis = np.arange(lo, hi + 1)
        axb = axes[2]
        axb.bar(gis, d, width=0.8, color=_shift_colours(d), lw=0)
        axb.axhline(0, color=C["INK"], lw=0.7)
        lim = max(10.0, np.ceil(np.abs(d).max() / 10) * 10)
        axb.set_ylim(-lim * 1.25, lim * 1.25)
        axb.set_ylabel("Shift (m)")
        axb.set_xlabel("Model domain, 500 m each (south → north)")
        axb.grid(axis="y", lw=0.4, color="0.85")
        for edge in np.arange(lo - 0.5, hi + 1.0):
            axb.axvline(edge, color="0.8", lw=0.5, ls=(0, (2, 3)), zorder=0)
        _domain_ticks(axb, gis)
        for y, va, text in ((0.97, "top", "↑ Landward with the shoreline offset"),
                            (0.03, "bottom", "↓ Seaward with the shoreline offset")):
            axb.text(0.005, y, text, transform=axb.transAxes, ha="left", va=va, fontsize=7,
                     color=C["LATE"] if va == "top" else C["EARLY"])
        _shift_note(axb, gap)
        _title(axb, 2, "Cross-shore shift of each domain from (a) to (b)")
        for ax, letter in zip(axes[:2], ("b", "a")):
            # A mid-grey box, so both the white and the black dash read in the key
            ax.legend(handles=[plt.Line2D([], [], color="white", lw=1.5, ls=(0, (3, 2)),
                                          label=f"Shoreline as placed in ({letter}): seaward of this panel"),
                               plt.Line2D([], [], color="black", lw=1.5, ls=(0, (3, 2)),
                                          label=f"Shoreline as placed in ({letter}): landward of this panel")],
                      loc="upper right", fontsize=7, facecolor="0.62", edgecolor="0.3", framealpha=1.0)

        sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
        cb = fig.colorbar(sm, ax=axes[:2], ticks=bounds[1:-1], fraction=0.025, pad=0.01)
        cb.set_label("Elevation (m above MHW)")
        cb.outline.set_linewidth(0.5)

        nseaward = int((d < 0).sum())
        frame_clause = (
            f"Offsets are in the model frame, each build zeroed on its own most seaward domain (GIS "
            f"{zero['duneline']} for the dune line, GIS {zero['shoreline']} for the shoreline), which is what "
            f"Cascade receives up to one shared constant (checked against a run's year-0 positions to within "
            f"1 m). Because of that zeroing the shift in (c) is not the distance from shoreline to dune line, "
            f"which is always seaward (beach width {beach.min():.0f}–{beach.max():.0f} m): it is "
            f"{gap:.1f} m, the gap between the two zero points, minus the beach width at the domain. A domain "
            f"with a beach wider than {gap:.0f} m moves seaward under the shoreline offset, a narrower one "
            f"landward. ")
        caption(fig, (
            f"GIS {lo}–{hi} of the {year} start island as the model builds it from each offset source, "
            f"one of {len(SECTIONS)} figures covering GIS 1–90. Every Barrier3D domain of the {year} start "
            f"topography ({product}, dune-topo {dune_version}) is drawn at its island offset: (a) with the "
            f"dune-line build (duneline/{ver['duneline']}), (b) with the shoreline build "
            f"(shoreline/{ver['shoreline']}), on the same axes. The grids are identical in both panels; "
            f"only where each domain sits cross-shore changes. Each column of cells is one 500 m domain: "
            f"its seaward edge is the BRIE shoreline, then the model's {beach_m} m initial beach, the two "
            f"dune rows and the interior (water masked); landward is up. The dash across each domain in (a) "
            f"and (b) is where its shoreline WOULD be with the other source: white, over the ocean, where the "
            f"other source puts the domain further seaward; black, over the beach, where it puts it further "
            f"landward. (c) "
            f"gives that move in metres per domain, from (a) to (b): blue up is landward, red down is "
            f"seaward; the shoreline offset places {nseaward} of {hi - lo + 1} domains here further seaward. "
            f"Grey strips along the bottom of (a) and (b) mark the villages. {frame_clause}Vertical "
            f"exaggeration is large."))
        stem = f"offset_{year}_duneline_vs_shoreline_domains_GIS{lo:02d}-{hi:02d}"
        for path in save(fig, out_dir / stem, close=True):
            print(f"  wrote {path}")

    _overview(year, off, zero, approach, ver, out_dir)


# The whole island on one panel: both planforms and where the six section figures fall (Hannah, 2026-10-06)
def _overview(year, off, zero, approach, ver, out_dir):
    gis = np.arange(1, 91)
    fig, ax = plt.subplots(figsize=figsize("double", height=3.6), constrained_layout=True)
    for src in ROWS:
        ax.step(gis, off[src], where="mid", color=ROWS[src]["colour"], lw=1.2, label=approach[src])
    ax.set_ylabel("Cross-shore position (m)")
    ax.set_xlabel("Model domain, 500 m each (south → north)")
    ax.legend(loc="upper right", fontsize=7, frameon=False)
    ax.set_title("Island planform from each offset, and the six section figures")

    # The six sections: alternate shading, their names along the bottom, a border between each
    for k, (lo, hi) in enumerate(SECTIONS):
        if k % 2:
            ax.axvspan(lo - 0.5, hi + 0.5, color="0.95", lw=0, zorder=0)
        if k:
            ax.axvline(lo - 0.5, color="0.55", lw=0.6, ls=(0, (3, 2)), zorder=1)
        ax.text((lo + hi) / 2, 0.02, f"GIS {lo}–{hi}", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=7, color=C["INK_MUTED"])
    ax.set_xlim(0.5, 90.5)
    ax.grid(axis="y", lw=0.4, color="0.88")
    _domain_ticks(ax, gis)
    ax.tick_params(axis="x", which="minor", length=0)

    caption(fig, (
        f"The {year} island offset from each source over the whole island, the overview for the "
        f"{len(SECTIONS)} section figures beside it (shaded bands, named along the bottom). The cross-shore "
        f"position each offset gives every 500 m domain: {approach['duneline'].lower()} "
        f"(duneline/{ver['duneline']}) and {approach['shoreline'][0].lower() + approach['shoreline'][1:]} "
        f"(shoreline/{ver['shoreline']}), each in the model frame, zeroed on its own most seaward domain "
        f"(GIS {zero['duneline']} and GIS {zero['shoreline']}); landward is up. At this scale the two lines "
        f"overlap: they differ by tens of metres on a planform that spans about 6 km, which is why the "
        f"section figures zoom to 15 domains and give the shift per domain."))
    stem = f"offset_{year}_duneline_vs_shoreline_domains_overview"
    for path in save(fig, out_dir / stem, close=True):
        print(f"  wrote {path}")


if __name__ == "__main__":
    main()
