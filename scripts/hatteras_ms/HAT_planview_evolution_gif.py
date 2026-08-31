#!/usr/bin/env python3
"""Animated plan view of the island through a run: elevation, not shoreline.

WHAT WAS MISSING, AND WHY THIS EXISTS
    Every matrix run already writes four GIFs and all four animate the
    SHORELINE -- a one-dimensional cross-shore position per domain, plotted
    against domain number. None of them shows elevation. The plan-view canvas
    that `init_planview` builds, the one the initialization figures use, was
    only ever drawn at t = 0.

    So there was no way to watch the barrier itself evolve: where the interior
    lowers, where overwash reaches, where the island narrows. This draws that
    canvas once per model year.

THE HISTORY IS ALREADY ON DISK -- no re-run is needed
    Barrier3D keeps `DomainTS`, one interior grid per year, and `x_s_TS`, the
    shoreline position per year, for all 120 padded domains. Both survive in
    the run's `.npz`, which is why those files are ~300 MB rather than the
    ~20 KB the shoreline matrix costs. This reads them and plots; it does not
    re-run anything.

TWO THINGS THAT HAVE TO BE RIGHT
    RAGGED GRIDS. `DomainTS[t]` is NOT a fixed shape -- Barrier3D trims water
    rows, so one domain-year is (174, 50) and another (200, 50), and the count
    changes as the island evolves. Each grid is put back on the full frame with
    `pad_cross_shore`, exactly as `load_domain_grids` does for the static
    figure, so the two are directly comparable.

    A FIXED FRAME. `build_canvas` sizes the canvas from the offsets it is
    given, so a per-year canvas is a per-year height and the animation would
    breathe. Offsets are taken against ONE reference for the whole run, and the
    axes get one y limit for every frame, so what moves in the GIF is the
    island and not the camera.

UNITS
    `DomainTS` is in decameters and `x_s_TS` likewise; the plan-view config
    carries `dam_to_m = 10.0` and a 10 m cell, so a decameter of shoreline
    movement is exactly one canvas row. Elevations are converted to meters to
    match the shared colorbar.

Usage:
    python HAT_planview_evolution_gif.py <run_directory> [--fps 3] [--out PATH]

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[2]
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no pyproject.toml.")
for _path in (PROJECT_BASE_DIR / "scripts",):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from hatteras_site_config import HATTERAS_DOMAINS as GEOMETRY  # noqa: E402
from cascade_pipeline.plotting.init_planview import (  # noqa: E402
    DEFAULT_PLAN_VIEW, build_canvas, pad_cross_shore, plot_canvas)

HOLD_FRAMES = 5      # frames held on the final year, so it can be read


def load_history(run_dir):
    """Per-year interior grids and shoreline positions from a run's .npz.

    Args:
        run_dir: A matrix run directory holding exactly one .npz.

    Returns:
        (grids_by_year, offsets_by_year, start_year) where grids_by_year[t] is
        a padded-order list of metre-valued arrays on the full cross-shore
        frame, and offsets_by_year[t] is the matching per-domain row origin.

    Raises:
        SystemExit: If no .npz is present, or it holds no domain history.
    """
    run_dir = Path(run_dir)
    matches = sorted(run_dir.glob("*.npz"))
    if not matches:
        raise SystemExit(
            f"no .npz model state in {run_dir}.\n"
            "The run was made with --no-model-state, and the per-year grids "
            "this needs were never written. Re-run that cell without it.")

    cascade = np.load(matches[0], allow_pickle=True)["cascade"].item()
    barrier3d = cascade._barrier3d
    inner = [getattr(model, "_model", model) for model in barrier3d]

    n_years = len(inner[0].DomainTS)
    config = DEFAULT_PLAN_VIEW

    # ONE reference for the whole run, so the canvas does not breathe.
    reference = min(float(model.x_s_TS[0]) for model in inner)

    grids_by_year, offsets_by_year = [], []
    for year in range(n_years):
        grids_by_year.append([
            pad_cross_shore(np.asarray(model.DomainTS[year], dtype=float)
                            * config.dam_to_m, config)
            for model in inner])
        offsets_by_year.append(np.array(
            [float(model.x_s_TS[year]) - reference for model in inner]))

    start_year = int(getattr(cascade, "_start_year", 0)) or None
    return grids_by_year, offsets_by_year, start_year


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run_dir", help="matrix run directory")
    parser.add_argument("--fps", type=float, default=3.0)
    parser.add_argument("--out", default=None)
    parser.add_argument("--start-year", type=int, default=None,
                        help="calendar year of frame 0, for the title")
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    grids_by_year, offsets_by_year, detected = load_history(run_dir)
    n_years = len(grids_by_year)

    # The run name carries the period, so the calendar year is recoverable
    # without trusting an attribute that may not exist.
    start_year = args.start_year or detected
    if start_year is None:
        for token in run_dir.name.split("_"):
            if token.isdigit() and len(token) == 4:
                start_year = int(token)
                break

    canvases = [build_canvas(grids, offsets, GEOMETRY)
                for grids, offsets in zip(grids_by_year, offsets_by_year)]
    # ONE y limit for every frame: what moves should be the island.
    frame_rows = max(canvas[0].shape[0] for canvas in canvases)

    figure, axis = plt.subplots(figsize=(16, 5.4))
    out = Path(args.out) if args.out else (
        run_dir / f"{run_dir.name}_planview_evolution.gif")

    def frame(index):
        year_index = min(index, n_years - 1)
        axis.clear()
        canvas, starts, per_domain, first_real = canvases[year_index]
        label = (f"{start_year + year_index}   (year {year_index})"
                 if start_year else f"year {year_index}")
        plot_canvas(canvas, starts, per_domain, first_real, GEOMETRY,
                    title=f"{run_dir.name}          {label}",
                    ax=axis, colorbar=False,
                    xlabel="Domain (S -> N,  Cape Hatteras to Rodanthe)")
        axis.set_ylim(0, frame_rows)
        axis.set_ylabel("Cross-shore cell  (10 m)")

    # The colorbar is drawn once, outside the frame loop: plot_canvas would
    # otherwise add a new one on every frame and shrink the axes each time.
    frame(0)
    mesh = axis.collections[0]
    bar = figure.colorbar(mesh, ax=axis, pad=0.01, fraction=0.02)
    bar.set_label("Elevation (m MHW)")
    bar.set_ticks([-1, 0, 1, 2, 3, 4])

    animation = FuncAnimation(figure, frame,
                              frames=n_years + HOLD_FRAMES, blit=False)
    animation.save(out, writer=PillowWriter(fps=args.fps))
    plt.close(figure)
    print(f"wrote {out}")
    print(f"  {n_years} model years + {HOLD_FRAMES} held, "
          f"{GEOMETRY.num_real_domains} real domains, frame {frame_rows} rows")


if __name__ == "__main__":
    main()
