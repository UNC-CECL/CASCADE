#!/usr/bin/env python3
"""Why the largest residual in the hindcast is deliberately left uncorrected.

At convergence the source/sink calibration leaves its biggest misfit at D6 --
2.00 m/yr in period 1 and 2.59 m/yr in period 2, roughly twice the next worst
domain in either. That looks like a calibration failure and is not one. D5-D7
are the Buxton groin's footprint, held in GROIN_RESERVED_DOMAINS, and the
residual there is the GROIN's shortfall rather than a background-erosion term.

WHY IT WOULD BE WRONG TO CORRECT IT
    The groin's trapping rate M and deterioration floor f were fitted against
    the observed shoreline, and the source/sink field is then derived from what
    the modules could NOT explain -- which is why the calibration runs against a
    groin-ON base run in the first place (GROIN_AWARE_BASE_RUN). Letting BE
    absorb the residual at D5-D7 would close the same gap twice: the groin would
    score as well-calibrated because a source term was quietly doing its work,
    and the M/f fit could never be falsified by the hindcast.

    So the number stays visible. It is the honest statement of what the groin
    module cannot do.

WHAT THE TWO SIGNS MEAN, AND WHY THEY ARE OPPOSITE
    period 1   residual POSITIVE -- observed is more seaward than modelled.
               The model does not build enough fillet. M = 60 m/yr is the most
               the sediment budget will support (719,000 m3/yr against a
               5-7e5 littoral drift), so this is a bound, not a missed fit.

    period 2   residual NEGATIVE -- modelled is more seaward than observed.
               The real fillet RELEASED after the 2003 storm damage; the module
               cannot, because trapping is bounded at >= 0, so it can stop
               adding sand but never remove it. This is outside the
               parameterisation at any (M, f), not a badly chosen one.

    The opposite signs are the point. A source/sink term fitted to close both
    would have to change sign between periods at the same domain, which is a
    fitted constant standing in for a structure that was built, damaged and
    left -- exactly the kind of thing the zone rules exist to keep out.

Usage:
    python HAT_groin_reserved_residual_figure.py

Reads  the converged calibBE full_management runs, groin on and off, plus the
       live GROIN_RESERVED_DOMAINS.
Writes output/fig_groin_reserved_residual.png

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import importlib.util
import pathlib
import sys

import numpy as np
import pandas as pd

_HERE = pathlib.Path(__file__).resolve()
PROJECT_BASE_DIR = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
OUTPUT_DIR = _HERE.parent / "output"
RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"

sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))

SHOW = list(range(2, 13))
PERIODS = {
    "1984_2004": dict(label="1984–2004", start=1984, scenario="road_bdm",
                      colour="#1565C0"),
    "2004_2024": dict(label="2004–2024", start=2004, scenario="road_bdm_nourish",
                      colour="#B71C1C"),
}
RESERVED_COLOUR = "#FF8C00"


def analysis_module():
    spec = importlib.util.spec_from_file_location(
        "_loess_analysis", _HERE.parent / "HAT_be_zone_LOESS_analysis.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_rates(period, scenario, groin):
    """One run's per-domain LRR, or None when the run is absent."""
    name = f"HAT_{period}_calibBE_{scenario}_{'groin' if groin else 'nogroin'}"
    path = RAW_RUNS / period / "calibBE" / name / f"{name}_shoreline_change_rate.csv"
    if not path.exists():
        return None
    return pd.read_csv(path).set_index("gis_domain")["lrr_m_yr"]


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    module = analysis_module()
    reserved = list(module.GROIN_RESERVED_DOMAINS)

    data = {}
    for key, meta in PERIODS.items():
        csv = module.P1_COASTSAT_CSV if meta["start"] == 1984 else module.P2_COASTSAT_CSV
        target = module.load_observed(meta["start"], csv)[1]
        on = run_rates(key, meta["scenario"], True)
        off = run_rates(key, meta["scenario"], False)
        if on is None or off is None:
            raise FileNotFoundError(
                f"{key}: need both the groin-on and groin-off calibBE "
                f"{meta['scenario']} runs; one is missing. Run "
                f"HAT_run_all.py --stages 2,6 --presets calibBE first.")
        data[key] = dict(target=target, on=on, off=off, **meta)

    figure = plt.figure(figsize=(15.5, 8.2))
    grid = figure.add_gridspec(2, 2, width_ratios=[1.25, 1], hspace=0.34,
                               wspace=0.22)
    axes = [figure.add_subplot(grid[0, 0]), figure.add_subplot(grid[1, 0])]
    bars = figure.add_subplot(grid[:, 1])

    # ---- LEFT: observed against the model, groin on and off ---------------
    for axis, (key, d) in zip(axes, data.items()):
        x = np.array(SHOW, dtype=float)
        tgt = np.array([d["target"].get(g, np.nan) for g in SHOW])
        on = np.array([d["on"].get(g, np.nan) for g in SHOW])
        off = np.array([d["off"].get(g, np.nan) for g in SHOW])

        axis.axvspan(min(reserved) - 0.5, max(reserved) + 0.5,
                     color=RESERVED_COLOUR, alpha=0.16, zorder=0)
        axis.text((min(reserved) + max(reserved)) / 2, 0.97,
                  "RESERVED\nfor the groin", ha="center", va="top", fontsize=8,
                  color="#B36200", weight="bold", zorder=6,
                  transform=axis.get_xaxis_transform())

        axis.plot(x, tgt, marker="s", markersize=7, linestyle="--",
                  color="#1A1A1A", linewidth=2.4, zorder=6, label="observed")
        axis.plot(x, on, marker="o", markersize=6, color=d["colour"],
                  linewidth=2.2, zorder=5, label="model, groin ON")
        axis.plot(x, off, marker="^", markersize=6, color="#888888",
                  linewidth=1.8, linestyle=":", zorder=4, label="model, groin OFF")

        worst = int(np.nanargmax(np.abs(tgt - on)))
        # Point at the MIDDLE of the gap, and park the text well clear of both
        # curves -- anchored on the model line it overlapped it.
        midpoint = (tgt[worst] + on[worst]) / 2.0
        axis.annotate(
            f"D{SHOW[worst]}  residual {tgt[worst] - on[worst]:+.2f} m/yr",
            xy=(x[worst], midpoint), xytext=(0.62, 0.12 if d["start"] == 1984 else 0.88),
            textcoords=axis.transAxes, fontsize=9.5, color=d["colour"],
            weight="bold", ha="left",
            va="bottom" if d["start"] == 1984 else "top",
            arrowprops=dict(arrowstyle="->", color=d["colour"], linewidth=1.4,
                            connectionstyle="arc3,rad=0.15"))
        axis.annotate("", xy=(x[worst], tgt[worst]), xytext=(x[worst], on[worst]),
                      arrowprops=dict(arrowstyle="<->", color=d["colour"],
                                      linewidth=1.6, alpha=0.75))

        axis.set_xticks(SHOW)
        axis.set_ylabel("shoreline rate (m/yr)\n[+ = seaward]", fontsize=9)
        axis.set_title(f"{d['label']}", fontsize=11, loc="left")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8, loc="best")
    axes[1].set_xlabel("GIS domain (south → north)")

    # ---- RIGHT: the residual at the reserved domains ----------------------
    width = 0.36
    idx = np.arange(len(reserved), dtype=float)
    for offset, (key, d) in zip((-width / 2, width / 2), data.items()):
        resid_on = [d["target"].get(g, np.nan) - d["on"].get(g, np.nan)
                    for g in reserved]
        resid_off = [d["target"].get(g, np.nan) - d["off"].get(g, np.nan)
                     for g in reserved]
        bars.bar(idx + offset, resid_on, width, color=d["colour"], alpha=0.9,
                 zorder=4, label=f"{d['label']}  groin ON")
        # groin-off as an outline behind: the gap between the two IS the groin
        # One legend entry only: the two periods' outlines are visually
        # identical, so labelling both just doubles the legend.
        bars.bar(idx + offset, resid_off, width, facecolor="none",
                 edgecolor="#333333", linewidth=1.3, linestyle="--", zorder=5,
                 label="same run, groin OFF" if offset < 0 else None)
        for i, (a, b) in enumerate(zip(resid_on, resid_off)):
            bars.annotate(f"{a:+.2f}", xy=(idx[i] + offset, a),
                          xytext=(0, 5 if a >= 0 else -13),
                          textcoords="offset points", ha="center", fontsize=8,
                          color=d["colour"], weight="bold")

    bars.axhline(0.0, color="#333333", linewidth=1.0, zorder=3)
    bars.set_xticks(idx)
    bars.set_xticklabels([f"D{g}" for g in reserved], fontsize=11)
    bars.set_ylabel("residual, observed − modelled (m/yr)")
    bars.set_title("WHAT IS LEFT AT THE RESERVED DOMAINS\n"
                   "the gap between filled and dashed is the groin's own "
                   "contribution", fontsize=11)
    bars.grid(alpha=0.25, axis="y")
    bars.legend(fontsize=8, loc="best")

    # Axes coords, hard against the left edge: in data coords these collided
    # with the D5 bars and their value labels.
    bars.text(0.015, 0.955,
              "observed more seaward\nmodel builds too LITTLE fillet",
              transform=bars.transAxes, ha="left", va="top", fontsize=8.5,
              color="#1565C0", style="italic")
    bars.text(0.015, 0.045,
              "model more seaward\nmodel cannot RELEASE the fillet",
              transform=bars.transAxes, ha="left", va="bottom", fontsize=8.5,
              color="#B71C1C", style="italic")

    figure.suptitle("The largest residual in the hindcast is the groin's, and is "
                    "left uncorrected on purpose", fontsize=13, y=0.985)
    figure.tight_layout(rect=(0, 0.115, 1, 0.965))
    figure.text(
        0.01, 0.012,
        "WHY IT IS NOT CORRECTED. M and f were fitted against the observed shoreline, and the source/sink field is then derived from what the "
        "modules could NOT explain -- which is why the calibration runs against a groin-ON base. Letting background erosion absorb D5-D7 would "
        "close the same gap twice: the groin would score as well-calibrated because a source term was doing its work, and the M/f fit could never "
        "be falsified by the hindcast.\n"
        "WHY THE SIGNS ARE OPPOSITE. Period 1 is positive -- observed is more seaward, so the model builds too little fillet, and M = 60 m/yr is "
        "already ~719,000 m3/yr against a 5-7e5 littoral drift, so it is a budget bound rather than a missed fit. Period 2 is negative -- the real "
        "fillet released after the 2003 storm damage and the module cannot, since trapping is bounded at >= 0. A single BE term closing both would "
        "have to change sign between periods at the same domain: a fitted constant standing in for a structure that was built, damaged and left.",
        fontsize=7.4, color="#333333", wrap=True)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / "fig_groin_reserved_residual.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)

    print(f"wrote {path}")
    for key, d in data.items():
        line = "  ".join(
            f"D{g} {d['target'].get(g, np.nan) - d['on'].get(g, np.nan):+.2f}"
            for g in reserved)
        print(f"  {d['label']}  residual at reserved domains:  {line}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
