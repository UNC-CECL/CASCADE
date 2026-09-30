"""
The padded island offset profile at true 1:1 scale, so the buffer's shape reads as it is.

    python scripts/input_prep/2-brie-offset/2-figures/HAT_plot_offset_profile_1to1.py --year 1996
    python scripts/input_prep/2-brie-offset/2-figures/HAT_plot_offset_profile_1to1.py --all

One figure per padded build, beside the build. Details: scripts/input_prep/2-brie-offset/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-23
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


# The padded build for a year and geometry
def padded_file(year, geometry):
    first, last = gis_bounds(geometry)
    n = (last - first + 1) + 2 * BUFFER_DOMAINS_PER_SIDE
    # Through hat_topo_version since 2026-09-22, when every build moved under <year>/<source>/
    if geometry == BASE_GEOMETRY:
        path = offset_file(year, "padded", n, source=SOURCE)
    else:
        path = (offset_start_dir(year, SOURCE) / "ext" / geometry
                / f"{offset_basename(year, SOURCE)}_PADDED_{n}.csv")
    if not path.is_file():
        sys.exit(f"no padded file at {path}")
    return path


# (geometry, first, last) from a padded file's length
def geometry_of(path):
    n = len(pd.read_csv(path))
    for name, (first, last) in GEOMETRIES.items():
        if (last - first + 1) + 2 * BUFFER_DOMAINS_PER_SIDE == n:
            return name, first, last
    return None


# Every padded build under 2-brie-offset, superseded ones included
def every_padded_file():
    out = []
    for src in OFFSET_SOURCES:
        stem = offset_basename(0, src).rsplit("_", 1)[0]   # drop the year
        out.extend(BRIE_ROOT.rglob(f"{stem}_*_PADDED_*.csv"))
    return sorted(set(out))


# The profile at equal axes, 1 m alongshore to 1 m cross-shore
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
    # Named from the file it was drawn FROM, not from the dune stem (2026-09-22)
    out = path.parent / f"{path.stem.rsplit('_PADDED_', 1)[0]}_buffer_diagnostic_1to1.png"
    save(fig, out, close=True)
    print(f"wrote {out.relative_to(BRIE_ROOT).as_posix()}")


# Run: one build, or every padded build on disk
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
