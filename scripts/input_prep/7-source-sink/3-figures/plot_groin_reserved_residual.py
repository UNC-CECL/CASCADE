#!/usr/bin/env python3
"""
Why the largest residual in the hindcast, at D6, is deliberately left uncorrected.

    python scripts/input_prep/7-source-sink/3-figures/plot_groin_reserved_residual.py

The residual around the Buxton groin in both periods, from the calibBE runs. Details: scripts/input_prep/7-source-sink/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import importlib.util
import pathlib
import sys

import numpy as np
import pandas as pd

_HERE = pathlib.Path(__file__).resolve()
PROJECT_BASE_DIR = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
# --- CONFIG ------------------------------------------------------------------
RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
# The figure lives with the rest of the section 7 figures, in the data tree.
sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))
from site_layer import hat_source_sink as _be  # noqa: E402
FIG_DIR = _be.figures_dir()   # the default pair's (2026-09-18)

from cascade_pipeline.run_layout import resolve as resolve_run_file  # noqa: E402
from site_layer.hat_figure_style import (                                   # noqa: E402
    apply_style, figsize, save, caption, town_bands, open_frame,
    DOMAIN_AXIS_LABEL, C, C_1984, C_1997, INK, INK_MUTED, halo, _title)

SHOW = list(range(2, 13))
# The earlier period takes the red of the house vintage pair, the later the blue.
PERIODS = {
    "1984_2004": dict(label="1984–2004", start=1984, scenario="road_bdm",
                      colour=C_1984),
    "2004_2024": dict(label="2004–2024", start=2004, scenario="road_bdm_nourish",
                      colour=C_1997),
}
# -----------------------------------------------------------------------------


# be_zone_residual_fit.py, loaded as a module for its LOWESS helpers
def analysis_module():
    spec = importlib.util.spec_from_file_location(
        "_lowess_analysis",
        _HERE.parent.parent / "2-calibrate" / "be_zone_residual_fit.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# One run's per-domain LRR, or None when the run is absent
def run_rates(period, scenario, groin):
    name = f"HAT_{period}_calibBE_{scenario}_{'groin' if groin else 'nogroin'}"
    # Resolved, not joined: the rate CSV's path depends on the run layout
    path = resolve_run_file(RAW_RUNS / period / "calibBE" / name,
                            "rate_csv", name)
    if not path.exists():
        return None
    return pd.read_csv(path).set_index("gis_domain")["lrr_m_yr"]


# Run: the figure
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

    apply_style()
    figure = plt.figure(figsize=figsize("double", aspect=0.55),
                        constrained_layout=True)
    grid = figure.add_gridspec(2, 2, width_ratios=[1.3, 1])
    axes = [figure.add_subplot(grid[0, 0]), figure.add_subplot(grid[1, 0])]
    bars = figure.add_subplot(grid[:, 1])

    # Left: observed against the model, groin on and off
    for i, (axis, (key, d)) in enumerate(zip(axes, data.items())):
        x = np.array(SHOW, dtype=float)
        tgt = np.array([d["target"].get(g, np.nan) for g in SHOW])
        on = np.array([d["on"].get(g, np.nan) for g in SHOW])
        off = np.array([d["off"].get(g, np.nan) for g in SHOW])

        axis.plot(x, tgt, marker="s", markersize=3.5, linestyle=(0, (4, 2)),
                  color=C["REF"], linewidth=1.6, zorder=6, label="observed")
        axis.plot(x, on, marker="o", markersize=3.5, color=d["colour"],
                  linewidth=1.4, zorder=5, label="modelled, groin present")
        axis.plot(x, off, marker="^", markersize=3.5, color=C["BASE"],
                  linewidth=1.1, linestyle=":", zorder=4,
                  label="modelled, groin absent")

        worst = int(np.nanargmax(np.abs(tgt - on)))
        # The gap is drawn; its size goes in the caption and on the bars
        axis.annotate("", xy=(x[worst], tgt[worst]),
                      xytext=(x[worst], on[worst]),
                      arrowprops=dict(arrowstyle="<->", color=d["colour"],
                                      linewidth=1.0))

        axis.set_xticks(SHOW)
        axis.set_xlim(min(SHOW) - 0.6, max(SHOW) + 0.6)
        axis.set_ylabel("shoreline rate (m/yr)\npositive is seaward")
        _title(axis, i, d["label"])
        axis.axhline(0.0, color=INK, lw=0.7, zorder=2)
        axis.grid(axis="y")
        open_frame(axis)
        # The reserved reach is the groin module's own footprint, so it takes the accent tint rather than a grey
        axis.axvspan(min(reserved) - 0.5, max(reserved) + 0.5,
                     color=C["ACCENT_FILL"], alpha=0.40, lw=0, zorder=0)
        if i == 0:
            axis.text((min(reserved) + max(reserved)) / 2, 0.03,
                      "reserved for the groin", ha="center", va="bottom",
                      fontsize=7, color=INK_MUTED, zorder=6,
                      transform=axis.get_xaxis_transform())
        town_bands(axis, where="top", strip=0.10, fontsize=7)
        axis.legend(loc="upper right", frameon=False, fontsize=7)
    axes[1].set_xlabel(DOMAIN_AXIS_LABEL)

    # Right: the residual at the reserved domains
    width = 0.36
    idx = np.arange(len(reserved), dtype=float)
    for offset, (key, d) in zip((-width / 2, width / 2), data.items()):
        resid_on = [d["target"].get(g, np.nan) - d["on"].get(g, np.nan)
                    for g in reserved]
        resid_off = [d["target"].get(g, np.nan) - d["off"].get(g, np.nan)
                     for g in reserved]
        bars.bar(idx + offset, resid_on, width, color=d["colour"], zorder=4,
                 label=f"{d['label']}, groin present")
        # Groin-off as an outline behind
        bars.bar(idx + offset, resid_off, width, facecolor="none",
                 edgecolor=C["BASE"], linewidth=0.9, linestyle=(0, (3, 2)),
                 zorder=5,
                 label="the same run with the groin absent" if offset < 0 else None)
        for i, a in enumerate(resid_on):
            bars.annotate(f"{a:+.2f}", xy=(idx[i] + offset, a),
                          xytext=(0, 4 if a >= 0 else -12),
                          textcoords="offset points", ha="center", fontsize=7,
                          color=d["colour"], zorder=8,
                          path_effects=halo(2.0))

    bars.axhline(0.0, color=INK, linewidth=0.7, zorder=3)
    # Room under the bars for their labels, and a band above them for the key
    lo, hi = bars.get_ylim()
    bars.set_ylim(lo - 0.10 * (hi - lo), hi + 0.38 * (hi - lo))
    bars.set_xticks(idx)
    bars.set_xticklabels([f"D{g}" for g in reserved])
    bars.set_ylabel("residual, observed \u2212 modelled (m/yr)")
    bars.grid(axis="y")
    open_frame(bars)
    _title(bars, 2, "what is left at the reserved domains")
    bars.legend(loc="upper center", frameon=False, fontsize=7)

    caption(figure, (
        "The largest residual the source/sink calibration leaves anywhere in the "
        "hindcast sits at domain 6, and it is left uncorrected on purpose. "
        "Domains 5 to 7 are the Buxton groin's footprint, shaded here and held "
        "in the reserved set. (a, b) the observed shoreline rate against the "
        "calibrated model at the southern domains, with the groin module on and "
        "off; the double-headed arrow marks the domain where observed and "
        "modelled are furthest apart. (c) the residual at each reserved domain "
        "for both periods, filled with the groin present and outlined with it "
        "absent -- the gap between the two is the groin's own contribution. "
        "The groin's trapping rate and deterioration floor were fitted against "
        "the observed shoreline, and the source/sink field is then derived from "
        "what the modules could not explain, which is why the calibration runs "
        "against a base run with the groin on. Letting background erosion absorb "
        "domains 5 to 7 would close the same gap twice: the groin would score as "
        "well calibrated because a source term was quietly doing its work, and "
        "the trapping fit could never be falsified by the hindcast. The two "
        "periods' residuals have opposite signs, and that is the point. In "
        "1984-2004 the residual is positive -- observed is further seaward than "
        "modelled, so the model builds too little fillet -- and a trapping rate "
        "of 60 m/yr already moves about 719,000 m3/yr against a littoral drift "
        "of 5-7 x 10^5 m3/yr, so it is a budget bound rather than a missed fit. "
        "In 2004-2024 it is negative -- modelled is further seaward than "
        "observed -- because the real fillet released after the 2003 storm "
        "damage and the module cannot, since its trapping is bounded at zero "
        "from below. One background-erosion term closing both would have to "
        "change sign between the periods at the same domain: a fitted constant "
        "standing in for a structure that was built, damaged and left. Domain 1 "
        "is at Cape Point and domain 90 at Pea Island."))

    path = FIG_DIR / "3-limits" / "fig_groin_reserved_residual.png"
    save(figure, path)
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
