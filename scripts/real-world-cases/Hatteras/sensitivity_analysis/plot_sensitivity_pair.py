#!/usr/bin/env python3
"""
Two settings of one sweep axis, stacked, each against the observed target.

    python scripts/sensitivity_analysis/plot_sensitivity_pair.py --start-year 1996 --sweep wave_angle_high_fraction --top 0.5 --bottom 0.51

For neighbouring values whose curves sit on top of each other in
plot_sensitivity.py. Writes pair_<...>.png beside that script's figures. Details: scripts/sensitivity_analysis/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import argparse

import numpy as np
import matplotlib.pyplot as plt

import plot_sensitivity as ps


# (run_dir, run_name, index row) for one value of the axis
def run_for(value, cells, index, sweep, start_year, preset, base_name):
    default = ps.normalise(ps.sweep_base_value(ps.SWEEPS[sweep]["setting"]))
    if abs(value - default) < 1e-9:
        run_dir = ps.find_run_dir(ps.RAW_RUNS, base_name,
                                  ps.period_component(start_year), preset,
                                  kind=ps.MATRIX_KIND)
        return run_dir, base_name, index.loc[(base_name, ps.MATRIX_KIND, "")]
    hit = cells[(cells.sweep == sweep) & (np.abs(cells.sort_key - value) < 1e-9)]
    if hit.empty:
        raise SystemExit(f"no completed cell at {sweep} = {value}")
    cell = hit.iloc[0]
    return cell.run_dir, cell.run_name, index.loc[cell.key]


# Run: the two cells, stacked on shared axes over the target
def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--start-year", type=int, required=True,
                        choices=sorted(ps.HATTERAS_PERIODS))
    parser.add_argument("--preset", default="edgeBE")
    parser.add_argument("--sweep", default="wave_angle_high_fraction",
                        choices=sorted(ps.SWEEPS))
    parser.add_argument("--top", type=float, required=True)
    parser.add_argument("--bottom", type=float, required=True)
    args = parser.parse_args()

    end_year = ps.HATTERAS_PERIODS[args.start_year].get(
        "end_year", args.start_year + 20)
    out_dir = ps.FIGURES_ROOT / f"{args.start_year}_{end_year}_{args.preset}"
    index = ps.load_index()
    cells = ps.load_cells(args.start_year, args.preset)
    block = cells[cells.sweep == args.sweep]
    base_name = block.iloc[0].base_name
    cs_series, target = ps.coastsat_layers(args.start_year)

    skip = ps.LOWESS_CONFIG.skip_southern_domains
    active = next(cs for cs in cs_series if cs["active"])
    south = np.asarray(active["transect_domains"]) <= skip
    dots_x = (np.asarray(active["transect_along_coast"])[south]
              / ps.HATTERAS_DOMAINS.domain_spacing_m
              + ps.HATTERAS_DOMAINS.first_gis_id)
    dots_y = np.asarray(active["transect_rates"])[south]

    label = ps.SWEEPS[args.sweep]["label"]
    default = ps.normalise(ps.sweep_base_value(ps.SWEEPS[args.sweep]["setting"]))
    fig, axes = plt.subplots(2, 1, figsize=(9.6, 6.6), sharex=True,
                             sharey=True, constrained_layout=True)
    colors = (ps.CURRENT_COLOR, ps.MODEL_RAMP(ps.RAMP_HI))
    for ax, value, color, letter in zip(axes, (args.top, args.bottom), colors,
                                        "ab"):
        run_dir, run_name, row = run_for(value, cells, index, args.sweep,
                                         args.start_year, args.preset,
                                         base_name)
        rates = ps.model_rates(run_dir, run_name)
        ax.scatter(dots_x, dots_y, color=ps.DEFAULT_RATE_COMPARISON.raw_color,
                   s=7, alpha=0.6, linewidths=0, zorder=2,
                   label=f"CoastSat transects, D1–{skip}")
        ax.plot(target.gis_domain, target.target_lrr_m_yr, color="#08306B",
                lw=1.8, zorder=5,
                label=f"CoastSat LRR, {ps.TARGET_WINDOW}-domain LOWESS "
                      f"(domain means D1–{skip})")
        ax.plot(rates.gis_domain, rates.lrr_m_yr, color=color, lw=2.0,
                zorder=6, label=f"Model, {label.lower()} {value:g}")
        ax.axhline(0.0, color=ps.INK_MUTED, lw=0.7, ls="--", zorder=1)
        tag = " (current)" if abs(value - default) < 1e-9 else ""
        ax.set_title(f"{label} {value:g}{tag}   ·   interior RMSE "
                     f"{row.rmse_interior_m_yr:.2f} m/yr, bias "
                     f"{row.mean_bias_interior_m_yr:+.2f} m/yr",
                     loc="left", pad=6)
        ax.set_ylabel("Shoreline change rate, LRR (m/yr)")
        ps.tidy(ax)
        ps.panel_label(ax, letter)
    axes[-1].set_xlim(ps.HATTERAS_DOMAINS.first_gis_id - 0.5,
                      ps.HATTERAS_DOMAINS.last_gis_id + 0.5)
    axes[-1].set_xlabel("Alongshore position (GIS domain, south → north)")
    handles, labels = axes[0].get_legend_handles_labels()
    bottom_h, bottom_l = axes[1].get_legend_handles_labels()
    handles = [handles[1], handles[2], bottom_h[2], handles[0]]
    labels = [labels[1], labels[2], bottom_l[2], labels[0]]
    fig.legend(handles, labels, loc="lower center",
               bbox_to_anchor=(0.5, -0.07), ncol=2, frameon=False)
    fig.suptitle(f"{label}, {args.start_year}–{end_year}, {args.preset}")

    stem = f"{args.sweep}_{args.top:g}_vs_{args.bottom:g}".replace(".", "p")
    path = out_dir / f"pair_{stem}.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
