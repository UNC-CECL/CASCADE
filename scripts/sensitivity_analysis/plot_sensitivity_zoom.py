#!/usr/bin/env python3
"""A zoomed view of one sweep axis over a value range, beside the full set.

plot_sensitivity.py draws every cell of an axis on one panel. A fine sweep
added later (e.g. high-angle fraction 0.51-0.54, 2026-09-29) would crowd that
panel, so this draws only the cells inside [--lo, --hi] with the plotter's own
functions and writes them to a sub-folder. The standard 01-06 set is untouched.

Usage:
    python plot_sensitivity_zoom.py --start-year 1996 \
        --sweep wave_angle_high_fraction --lo 0.5 --hi 0.55

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse

import plot_sensitivity as ps


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--start-year", type=int, required=True,
                        choices=sorted(ps.HATTERAS_PERIODS))
    parser.add_argument("--preset", default="edgeBE")
    parser.add_argument("--sweep", default="wave_angle_high_fraction",
                        choices=sorted(ps.SWEEPS))
    parser.add_argument("--lo", type=float, required=True)
    parser.add_argument("--hi", type=float, required=True)
    args = parser.parse_args()

    end_year = ps.HATTERAS_PERIODS[args.start_year].get(
        "end_year", args.start_year + 20)
    tag = f"{args.sweep}_{args.lo:g}_{args.hi:g}".replace(".", "p")
    out_dir = (ps.FIGURES_ROOT / f"{args.start_year}_{end_year}_{args.preset}"
               / f"zoom_{tag}")
    out_dir.mkdir(parents=True, exist_ok=True)

    index = ps.load_index()
    cells = ps.load_cells(args.start_year, args.preset)
    cells = cells[(cells.sweep == args.sweep)
                  & (cells.sort_key >= args.lo - 1e-9)
                  & (cells.sort_key <= args.hi + 1e-9)]
    if cells.empty:
        print("no cells in range")
        return 1
    ps.check_target_matches(index, list(cells.key))
    cs_series, target = ps.coastsat_layers(args.start_year)

    written = [
        ps.plot_skill_overview(cells, index, args.start_year, args.preset,
                               out_dir),
        ps.plot_alongshore(cells, args.sweep, args.start_year, args.preset,
                           cs_series, target, out_dir),
        ps.write_summary(cells, index, args.start_year, args.preset,
                         out_dir)[0],
    ]
    for path in written:
        if path is not None:
            print(f"  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
