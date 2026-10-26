#!/usr/bin/env python3
"""
A zoomed view of one sweep axis over a value range, beside the full set.

    python scripts/sensitivity_analysis/plot_sensitivity_zoom.py --start-year 1996 --sweep wave_angle_high_fraction --lo 0.5 --hi 0.55

Draws only the cells inside [--lo, --hi] with plot_sensitivity.py's own
functions, into a sub-folder of its figures; the standard set is untouched. Details: scripts/sensitivity_analysis/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-29
"""

from __future__ import annotations

import argparse

import plot_sensitivity as ps


# Run: the cells inside [--lo, --hi], drawn with plot_sensitivity's functions
def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
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
