"""
The padded offset profile at true 1:1 scale
============================================
The buffer diagnostic that island_offset_hybrid.py draws stretches 7.5 km of
cross-shore offset across a panel that spans 72 km alongshore, so the island
reads as a deep V. This draws the same padded profile with equal axes: 1 m
alongshore is 1 m cross-shore, the way a map draws it (Hannah, 2026-09-16).

    python HAT_plot_offset_profile_1to1.py --year 1996 --geometry n115
    python HAT_plot_offset_profile_1to1.py --year 1996              # the CURRENT surveyed build
    python HAT_plot_offset_profile_1to1.py --file 1996/v1/Island_Dune_Offsets_1996_PADDED_120.csv
    python HAT_plot_offset_profile_1to1.py --all                    # every padded build on disk

Writes <build dir>/<the build's own stem>_buffer_diagnostic_1to1.png beside
the original diagnostic (exaggerated axes, kept), with the PDF and caption
under supporting/. Named `_v2` until 2026-09-23, which read as a build number
inside a v<n>/ folder.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from site_layer.hat_extension_domains import (BASE_GEOMETRY, GEOMETRIES,  # noqa: E402
                                   BUFFER_DOMAINS_PER_SIDE, DOMAIN_SPACING_M,
                                   SURVEYED_GIS, gis_bounds)
from site_layer.hat_topo_version import (BRIE_ROOT,  # noqa: E402
                              DEFAULT_OFFSET_SOURCE as SOURCE, OFFSET_SOURCES,
                              offset_basename, offset_file, offset_start_dir)


def padded_file(year, geometry):
    first, last = gis_bounds(geometry)
    n = (last - first + 1) + 2 * BUFFER_DOMAINS_PER_SIDE
    # Through hat_topo_version since 2026-09-22, when every build moved under
    # <year>/<source>/. This read CURRENT out of <year>/ itself and would now
    # find no CURRENT there, silently falling back to the year folder.
    if geometry == BASE_GEOMETRY:
        path = offset_file(year, "padded", n, source=SOURCE)
    else:
        path = (offset_start_dir(year, SOURCE) / "ext" / geometry
                / f"{offset_basename(year, SOURCE)}_PADDED_{n}.csv")
    if not path.is_file():
        sys.exit(f"no padded file at {path}")
    return path


def geometry_of(path):
    """(geometry, first, last) from a padded file's length; None if no
    geometry pads to that many domains."""
    n = len(pd.read_csv(path))
    for name, (first, last) in GEOMETRIES.items():
        if (last - first + 1) + 2 * BUFFER_DOMAINS_PER_SIDE == n:
            return name, first, last
    return None


def every_padded_file():
    """Every padded build under 2-brie-offset, superseded ones included.

    One glob per SOURCE (2026-09-22). This matched only the dune stem, so
    --all quietly skipped every shoreline build; the stems come from
    hat_topo_version so a new source cannot be forgotten here again.
    """
    out = []
    for src in OFFSET_SOURCES:
        stem = offset_basename(0, src).rsplit("_", 1)[0]   # drop the year
        out.extend(BRIE_ROOT.rglob(f"{stem}_*_PADDED_*.csv"))
    return sorted(set(out))


def draw(path, geometry, first, last):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import C, C_1997, INK_MUTED, apply_style, caption, save
    apply_style()
    year = int(re.search(r"Offsets_([0-9]{4})_", path.name).group(1))
    build = path.parent.relative_to(BRIE_ROOT).as_posix()

    offset = pd.read_csv(path).iloc[:, 0].to_numpy(dtype=float)
    n = offset.size
    along_km = np.arange(n) * DOMAIN_SPACING_M / 1000.0
    off_km = offset / 1000.0
    buf = BUFFER_DOMAINS_PER_SIDE
    real = slice(buf, n - buf)
    gis = np.arange(first, last + 1)
    lo_s, hi_s = SURVEYED_GIS
    surveyed = (gis >= lo_s) & (gis <= hi_s)

    width_in = 7.48                      # a double-column figure
    span_x = along_km[-1] - along_km[0] + 1.0
    span_y = off_km.max() - off_km.min() + 1.0
    fig, ax = plt.subplots(figsize=(width_in, width_in * span_y / span_x + 0.9),
                           constrained_layout=True)
    # buffer domains: the invented coast that closes BRIE's ring
    ax.plot(along_km[:buf + 1], off_km[:buf + 1], color=INK_MUTED, lw=0.8, ls=":",
            label="buffer (extrapolated, then bridged)")
    ax.plot(along_km[n - buf - 1:], off_km[n - buf - 1:], color=INK_MUTED, lw=0.8, ls=":")
    # the surveyed reach and the extension
    x_real, y_real = along_km[real], off_km[real]
    ax.plot(x_real[surveyed], y_real[surveyed], color=C["BASE"], lw=1.4,
            label=f"surveyed reach, GIS {lo_s}-{hi_s}")
    if (~surveyed).any():
        k = np.where(~surveyed)[0]
        k = np.r_[k[0] - 1, k] if k[0] > 0 else k      # joined to the reach
        ax.plot(x_real[k], y_real[k], color=C_1997, lw=1.4,
                label=f"extension, GIS {hi_s + 1}-{last}")
    for g in (1, 30, 60, 90, last):
        if first <= g <= last:
            i = buf + (g - first)
            ax.annotate(f"GIS {g}", (along_km[i], off_km[i]), xytext=(0, 6),
                        textcoords="offset points", ha="center", fontsize=6.5,
                        color=INK_MUTED)
    ax.set_aspect("equal")
    ax.set_xlabel("alongshore, km (padded index x 0.5 km)")
    ax.set_ylabel("offset, km")
    ax.set_xlim(along_km[0] - 0.5, along_km[-1] + 0.5)
    ax.set_ylim(off_km.min() - 0.5, off_km.max() + 0.5)
    ax.legend(loc="lower center", ncol=3, fontsize=6.5, frameon=False,
              bbox_to_anchor=(0.5, 1.0))
    theta = np.degrees(np.arctan2(np.diff(offset), DOMAIN_SPACING_M))
    mean_bearing = np.degrees(np.arctan2(off_km[real].max() - off_km[real].min(),
                                         x_real[-1] - x_real[0]))
    caption(fig, f"The padded {year} dune-line offset the model reads for geometry "
                 f"{geometry}, build {build} (GIS {first} to {last} plus {buf} buffer "
                 f"domains each side), at 1:1 scale: one metre alongshore is one metre "
                 f"cross-shore, as a map draws it. The offset is the distance from a "
                 f"north-south datum line east of the island to the dune line, so the "
                 f"slope is the coast's bearing relative to north, not its curvature; the "
                 f"mean bearing is {mean_bearing:.0f} degrees over the reach and the "
                 f"steepest domain-to-domain angle is {np.abs(theta[real][:-1]).max():.0f} "
                 f"degrees. The compressed planform every calibrated run uses divides "
                 f"these offsets by ten.")
    # Named from the file it was drawn FROM, not from the dune stem
    # (2026-09-22): this was hardcoded, so the shoreline build's diagnostic
    # landed in 1996/shoreline/v1/ calling itself Island_Dune_Offsets.
    out = path.parent / f"{path.stem.rsplit('_PADDED_', 1)[0]}_buffer_diagnostic_1to1.png"
    save(fig, out, close=True)
    print(f"wrote {out.relative_to(BRIE_ROOT).as_posix()}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--year", type=int, default=1996)
    ap.add_argument("--geometry", default=BASE_GEOMETRY, choices=sorted(GEOMETRIES))
    ap.add_argument("--file", default=None,
                    help="a padded CSV, relative to 2-brie-offset or a path")
    ap.add_argument("--all", action="store_true", help="every padded build on disk")
    a = ap.parse_args(argv)

    if a.all:
        files = every_padded_file()
    elif a.file:
        f = Path(a.file)
        files = [f if f.is_file() else BRIE_ROOT / a.file]
    else:
        files = [padded_file(a.year, a.geometry)]
    for path in files:
        found = geometry_of(path)
        if found is None:
            print(f"skipped {path}: no geometry pads to its length")
            continue
        draw(path, *found)
    return 0


if __name__ == "__main__":
    sys.exit(main())
