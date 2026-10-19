"""
Summary statistics of the dune crests at the two hindcast starts (1996, 2009), over every 10 m cell.

    python scripts/figure_making/pipeline/1-barrier3d-domains/dune_crest_statistics.py

Reads each start's dune arrays as dune_crest_figures.py does and writes, in m
NAVD88 and as height above the berm: a table of the cell statistics (minimum,
percentiles, median, mean, maximum, SD) and of the 500 m domain medians, minima
and maxima; the share of cells under the design height and the Dmaxel ceiling;
a markdown summary; and a figure of the two distributions. Written beside the
arrays in data/hatteras_init/1-barrier3d-domains/figures/dune_crest_statistics/
and the figure also in output/figures/3-model-inputs/1-domains/.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-19
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parent))
import dune_crest_figures as DC  # noqa: E402
from site_layer import hat_topo_version as tv  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, INK, INK_MUTED, GRID_C, figsize, figure_dir, save, record_caption, support_dir,
    _letter_inside, open_frame,
)
OUT = figure_dir("inputs", "1-domains", "dune_crest_statistics")
DATA_OUT = tv.DOMAIN_FIGURES_DIR / "dune_crest_statistics"
NAME = "dune_crest_statistics_1996_2009"
PERCENTILES = (10, 25, 50, 75, 90, 99)
BIN_M = 0.25            # histogram bin
COUNT_ABOVE_M = (10.0, 11.0)   # tail counts quoted in the summary
# -----------------------------------------------------------------------------


def stats_row(x):
    x = pd.Series(np.ravel(x))
    row = dict(n=len(x), min=x.min())
    row.update({f"p{q}": x.quantile(q / 100) for q in PERCENTILES if q != 50})
    row.update(median=x.median(), mean=x.mean(), max=x.max(), sd=x.std())
    return {k: row[k] for k in ("n", "min", "p10", "p25", "median", "mean", "p75", "p90", "p99", "max", "sd")}


# One start's table: every cell, then the per-domain medians, minima and maxima (m NAVD88)
def table(year, c, d):
    rows = [dict(start=year, statistic_of="10 m cells", **stats_row(c))]
    for label, f in (("500 m domain medians", np.median), ("500 m domain minima", np.min),
                     ("500 m domain maxima", np.max)):
        rows.append(dict(start=year, statistic_of=label, **stats_row(f(c, axis=1))))
    t = pd.DataFrame(rows)
    t["height_above_berm_offset_m"] = d["berm"]
    return t


def shares(year, c, d):
    design = DC.DESIGN_MHW + d["mhw"]
    return dict(start=year, below_design=float((c < design).mean()), below_ceiling=float((c < d["ceil_navd"]).mean()),
                **{f"cells_at_or_above_{v:g}_m": int((c >= v).sum()) for v in COUNT_ABOVE_M},
                gis_of_min=int(np.unravel_index(c.argmin(), c.shape)[0] + 1),
                gis_of_max=int(np.unravel_index(c.argmax(), c.shape)[0] + 1),
                lowest_domains=", ".join(f"GIS {g} ({v:.2f})" for g, v in
                                         sorted(enumerate(np.median(c, 1), 1), key=lambda gv: gv[1])[:5]),
                highest_domains=", ".join(f"GIS {g} ({v:.2f})" for g, v in
                                          sorted(enumerate(np.median(c, 1), 1), key=lambda gv: -gv[1])[:5]))


def markdown(stats, sh, d, data):
    cells = stats[stats.statistic_of == "10 m cells"].set_index("start")
    names = [("min", "Minimum"), ("p10", "10th percentile"), ("p25", "25th percentile"), ("median", "**Median**"),
             ("mean", "**Mean**"), ("p75", "75th percentile"), ("p90", "90th percentile"),
             ("p99", "99th percentile"), ("max", "**Maximum**"), ("sd", "Standard deviation")]
    years = list(DC.STARTS)
    lines = ["# Dune crest statistics at the two hindcast starts", "",
             "The dune crest (highest point across the dune) at every 10 m cell along Hatteras Island, GIS 1-90 "
             "(about 45 km), from the dune arrays CASCADE starts from; model year 0 equals these. "
             f"Height above berm = elevation - {d['berm']:.1f} m (the model's fixed BermEl, not a measured berm). "
             "Written by scripts/figure_making/pipeline/1-barrier3d-domains/dune_crest_statistics.py; "
             f"the full table is {NAME}.csv beside this file.", "",
             "| | " + " | ".join(f"{y}: m NAVD88 | {y}: above berm (m)" for y in years) + " |",
             "|---|" + "---|---|" * len(years)]
    for k, label in names:
        vals = []
        for y in years:
            v = cells.loc[y, k]
            vals += [f"{v:.2f}", f"{v:.2f}" if k == "sd" else f"{v - d['berm']:.2f}"]
        if label.startswith("**"):
            vals = [f"**{v}**" for v in vals]
        lines.append(f"| {label} | " + " | ".join(vals) + " |")
    lines.append("| Cells | " + " | ".join(f"{int(cells.loc[y, 'n']):,} | " for y in years) + "|")
    s = {r["start"]: r for r in sh}
    design = DC.DESIGN_MHW + d["mhw"]
    lines += ["", "| Share of cells | " + " | ".join(str(y) for y in years) + " |", "|---|" + "---|" * len(years),
              f"| Below the design height, {design:.1f} m | " + " | ".join(f"{s[y]['below_design']:.0%}" for y in years) + " |",
              f"| Below the Dmaxel ceiling, {d['ceil_navd']:.1f} m | " + " | ".join(f"{s[y]['below_ceiling']:.0%}" for y in years) + " |"]
    lines += [f"| Cells at or above {v:g} m | " + " | ".join(str(s[y][f'cells_at_or_above_{v:g}_m']) for y in years) + " |"
              for v in COUNT_ABOVE_M]
    lines += ["", "By domain (500 m, median crest, m NAVD88):", ""]
    for y in years:
        lines.append(f"- {y}: lowest {s[y]['lowest_domains']}; highest {s[y]['highest_domains']}; "
                     f"lowest cell in GIS {s[y]['gis_of_min']}, highest in GIS {s[y]['gis_of_max']}.")
    lines += ["", "The maxima are a few isolated peaks near the 12 m clip of the DEM product; the 99th percentile "
              "is the more robust high value.", ""]
    return "\n".join(lines)


# The two distributions (a) and their summary statistics as range bars (b)
def figure(data, d):
    berm = d["berm"]
    fig, (ax, axs) = plt.subplots(2, 1, figsize=figsize("double", height=5.8), sharex=True,
                                  gridspec_kw=dict(height_ratios=[2.2, 1.4], hspace=0.12))
    bins = np.arange(1.5, 12.75, BIN_M)
    for year, (c, *_) in data.items():
        col = DC.STARTS[year][0]
        ax.hist(c.ravel(), bins=bins, weights=np.full(c.size, 100 / c.size), histtype="step",
                color=col, lw=1.4, label=f"{year} start ({c.size:,} cells of 10 m)")
    for ref, txt, style, ha in ((DC.DESIGN_MHW + d["mhw"], "NC-12 design height", ":", "right"),
                                (d["ceil_navd"], "Dune ceiling (Dmaxel)", "--", "left")):
        for a in (ax, axs):
            a.axvline(ref, color=INK_MUTED, lw=0.8, ls=style, zorder=0)
        ax.text(ref + (0.08 if ha == "left" else -0.08), ax.get_ylim()[1] * 0.99, txt,
                color=INK_MUTED, fontsize=7, va="top", ha=ha)
    ax.set_ylabel(f"Share of cells (%) per {BIN_M} m")
    ax.grid(axis="y", color=GRID_C, lw=0.5)
    ax.legend(loc="upper right", frameon=False, fontsize=7.5)
    _letter_inside(ax, 0)
    top = ax.secondary_xaxis("top", functions=(lambda v: v - berm, lambda v: v + berm))
    top.set_xlabel("Height above berm (m)")

    for k, (year, (c, *_)) in enumerate(data.items()):
        x = c.ravel()
        q = dict(zip(("p10", "p25", "p50", "p75", "p90", "p99"), np.percentile(x, [10, 25, 50, 75, 90, 99])))
        y, col = len(data) - 1 - k, DC.STARTS[year][0]
        axs.plot([x.min(), x.max()], [y, y], color=col, lw=0.8, solid_capstyle="butt")
        axs.plot([q["p10"], q["p90"]], [y, y], color=col, lw=4, alpha=0.45, solid_capstyle="butt")
        axs.add_patch(plt.Rectangle((q["p25"], y - 0.18), q["p75"] - q["p25"], 0.36, fc="white", ec=col, lw=1.2))
        axs.plot([q["p50"]] * 2, [y - 0.18, y + 0.18], color=col, lw=2)
        axs.plot(x.mean(), y, marker="D", ms=4.5, color=INK, ls="none")
        for v in (x.min(), x.max()):
            axs.plot(v, y, marker="|", ms=8, mew=1.2, color=col)
        axs.plot(q["p99"], y, marker="|", ms=6, mew=1.2, color=INK)
        for v, txt, dy, ha in ((x.min(), f"min {x.min():.2f}", -0.32, "left"),
                               (q["p50"], f"median {q['p50']:.2f}", 0.30, "center"),
                               (q["p99"], f"P99 {q['p99']:.2f}", 0.30, "center"),
                               (x.max(), f"max {x.max():.2f}", -0.32, "right")):
            axs.text(v, y + dy, txt, ha=ha, va="center", fontsize=6.5, color=INK)
    axs.set_yticks(range(len(data)), [f"{y} start" for y in data][::-1])
    axs.set_ylim(-0.6, len(data) - 0.4)
    axs.set_xlim(bins[0], bins[-1])
    axs.set_xlabel("Dune crest elevation (m NAVD88)")
    _letter_inside(axs, 1)
    axs.legend(handles=[Patch(fc="white", ec="0.3", label="25th-75th percentile, median line"),
                        Line2D([], [], color="0.3", lw=4, alpha=0.45, label="10th-90th percentile"),
                        Line2D([], [], color="0.3", lw=0.8, label="Minimum to maximum"),
                        Line2D([], [], color=INK, marker="|", ls="none", ms=6, label="99th percentile"),
                        Line2D([], [], color=INK, marker="D", ls="none", ms=4.5, label="Mean")],
               loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=3, frameon=False, fontsize=6.5)
    fig.subplots_adjust(bottom=0.2)
    for a in (ax, axs):
        open_frame(a)
    return fig


def main():
    apply_style()
    d = DC.datum()
    data = {y: (c + d["mhw"], p, v) for y, (c, p, v) in ((y, DC.crests(y, d["berm_mhw"])) for y in DC.STARTS)}
    stats = pd.concat([table(y, c, d) for y, (c, *_) in data.items()], ignore_index=True)
    sh = [shares(y, c, d) for y, (c, *_) in data.items()]
    sup = support_dir(DATA_OUT)
    sup.mkdir(parents=True, exist_ok=True)
    stats.round(3).to_csv(sup / f"{NAME}.csv", index=False)
    pd.DataFrame(sh).round(3).to_csv(sup / f"{NAME}_shares.csv", index=False)
    (sup / f"{NAME}.md").write_text(markdown(stats, sh, d, data), encoding="utf-8")
    fig = figure(data, d)
    paths = [OUT / f"{NAME}.png", DATA_OUT / f"{NAME}.png"]
    for p in paths:
        save(fig, p)
    plt.close(fig)
    (p96, v96), (p09, v09) = data[1996][1:], data[2009][1:]
    for p in paths:
        record_caption(p, (
            f"Dune crest elevation at the two hindcast starts, every 10 m cell of GIS 1-90 (4,500 cells each): "
            f"1996 start ({p96} {v96}, the 1996 ALACE survey; red) and 2009 start ({p09} {v09}, the 2009-2014 "
            f"lidar; blue), the dune arrays CASCADE starts from, in m NAVD88. (a) Share of cells per {BIN_M} m bin; "
            f"the top axis is height above the berm ({d['berm']:.1f} m NAVD88). (b) Minimum to maximum (thin line), "
            f"10th-90th percentile (thick band), 25th-75th percentile (box) with the median (line), 99th "
            f"percentile (black tick) and mean (diamond). Dotted: NC-12 design height "
            f"{DC.DESIGN_MHW + d['mhw']:.1f} m; dashed: Dmaxel ceiling {d['ceil_navd']:.1f} m NAVD88. "
            f"The numbers are in supporting/{NAME}.csv and {NAME}.md beside the arrays."))
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(stats.round(2).to_string(index=False))
    print(*paths, sup / f"{NAME}.md", sep="\n")


if __name__ == "__main__":
    main()
