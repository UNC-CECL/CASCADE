"""
Does the survey-fit dune growth rate follow beach width? Step 2 of the alongshore growth-rate question.

    python scripts/hatteras_ms/experiments/HAT_dune_growth_beach_width_check.py

Why (2026-10-21, Hannah): a wider beach gives the wind more fetch, so dunes
behind it should grow faster. If the per-domain rates of
HAT_dune_growth_alongshore_survey_fit.py (--widened --intervals) track beach
width, the width could set the alongshore pattern with one scale factor.

Beach width per domain: dune line to CoastSat shoreline, mean of the domain's
transects (data/hatteras_init/5-scr/4-comparisons/duneline_positions/supporting/
beach_width_domains.csv): 1997 for 1996-2009, 2009 for 2009-2025. Also against
the starting dune height (a logistic rate should not depend on it) and the
CoastSat net change of the window. Domains a fill reached inside the window are
left out; Rodanthe and the S-curves (GIS 76-88, rebuilt by the road managers)
drawn hollow and scored both ways. Spearman rank correlation. Read-only.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-21
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path[:0] = [str(PROJECT_ROOT / "scripts"), str(_HERE.parent)]
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hat_observed_rates import net_change_dir  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
FIT_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / "dunes" / "2026-10-21-alongshore-dune-growth-survey-fit"
FIT_SUFFIX = "_widened_at_edge_intervals"
BEACH_WIDTH = (PROJECT_ROOT / "data" / "hatteras_init" / "5-scr" / "4-comparisons" / "duneline_positions"
               / "supporting" / "beach_width_domains.csv")
WINDOWS = {1996: (2009, 1997), 2009: (2025, 2009)}      # start -> (end year, beach-width vintage at the start)
REBUILT = range(76, 89)                                 # Rodanthe / S-curves
AMBER = "#b35806"                                       # the survey-fit figure's domain-median colour
OUT = FIT_DIR / "figures" / "dune_growth_rate_vs_beach_width.png"
# -----------------------------------------------------------------------------


def spearman(x, y):
    return round(float(x.rank().corr(y.rank())), 2)


def data(start):
    end, vintage = WINDOWS[start]
    t = pd.read_csv(FIT_DIR / "tables" / f"domain_rates_{start}_start{FIT_SUFFIX}.csv").set_index("gis")
    bw = pd.read_csv(BEACH_WIDTH)
    w = f"{start}_{end}"
    nc = pd.read_csv(net_change_dir(start, end) / f"transect_net_change_{w}.csv").groupby("domain_number").net_change_m.mean()
    t["beach_width_m"] = bw[bw.dune_vintage == vintage].set_index("gis").beach_width_m
    t["coastsat_net_change_m"] = nc
    t["rebuilt"] = t.index.isin(REBUILT)
    return t[t.r_median.notna() & ~t.fill_in_window]


def main():
    fs.apply_style()
    proxies = {"beach_width_m": "Beach width at the start (m)",
               "h0_median_m": "Starting dune height above the berm (m)"}
    fig, axes = plt.subplots(len(WINDOWS), len(proxies), figsize=fs.figsize("double", height=5.0),
                             layout="constrained", sharey="row")
    rows = []
    for row, start in zip(axes, WINDOWS):
        t = data(start)
        for col in ("beach_width_m", "h0_median_m", "coastsat_net_change_m"):
            rows.append(dict(window=f"{start}-{WINDOWS[start][0]}", against=col, domains=len(t),
                             spearman_all=spearman(t[col], t.r_median),
                             spearman_without_76_88=spearman(t[col][~t.rebuilt], t.r_median[~t.rebuilt])))
        for ax, (col, xlabel) in zip(row, proxies.items()):
            for rebuilt, mfc, lab in ((False, AMBER, "Domain median"), (True, "white", "Rodanthe / S-curves (GIS 76-88)")):
                m = t.rebuilt == rebuilt
                ax.plot(t[col][m], t.r_median[m], "o", ms=4, color=AMBER, mfc=mfc, mew=1.0, ls="none", label=lab)
            keep = ~t.rebuilt
            ax.text(0.98, 0.95, f"Spearman {spearman(t[col], t.r_median):+.2f}\n"
                                f"without GIS 76-88 {spearman(t[col][keep], t.r_median[keep]):+.2f}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=7)
            ax.set_xlabel(xlabel)
            ax.grid(color=fs.GRID_C, lw=0.5)
            fs.open_frame(ax)
        row[0].set_ylabel(f"Survey-fit r, {start}–{WINDOWS[start][0] if start == 1996 else 2019}")
    for i, ax in enumerate(axes.flat):
        fs._letter_inside(ax, i)
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="outside lower center", ncol=2, frameon=False, fontsize=7)
    fs.save(fig, OUT, close=True)
    s = pd.DataFrame(rows)
    s.to_csv(FIT_DIR / "tables" / "growth_rate_vs_beach_width.csv", index=False)
    print(s.to_string(index=False))
    print(OUT)


if __name__ == "__main__":
    main()
