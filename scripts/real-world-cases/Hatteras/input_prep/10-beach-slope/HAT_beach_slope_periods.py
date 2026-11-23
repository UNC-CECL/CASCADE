"""
The long-term beach slope of Hatteras over three periods, one folder each (2026-11-23).

    python scripts/real-world-cases/Hatteras/input_prep/10-beach-slope/HAT_beach_slope_periods.py

Reads (D:, hatteras_init/10-beach-slope/1-foreshore-mhw-to-berm/): foreshore_slope_domains.csv
(HAT_beach_slope_mhw_berm.py, every survey). For each period: the surveys in it; a survey counts as well covered when
at least MIN_DOMAINS of the 90 domains gave a slope (the others are listed, drawn hollow and left out of the
"well-covered" values). The long-term value is taken two ways, the mean and the median, over the surveys' island
medians; per domain, the mean over the surveys of that domain's median.

Writes, per period (2-hindcast-1996-2009/, 3-test-2011-2020/, 4-all-1996-2020/):
    slope_by_domain.csv     per domain: mean and median over the period's well-covered surveys, number of surveys
    slope_summary.csv       the long-term island value, mean and median, all surveys and well-covered only
    beach_slope_<period>.png (caption in CAPTIONS.md)

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
"""
import sys
from pathlib import Path

_HAT_DIR = next(p for p in Path(__file__).resolve().parents if (p / "hatteras_paths.py").exists())
sys.path.insert(0, str(_HAT_DIR))

import numpy as np
import pandas as pd
from hatteras_paths import WORKSPACE

ROOT = WORKSPACE / "data" / "hatteras_init" / "10-beach-slope"
MEAS = ROOT / "1-foreshore-mhw-to-berm"
STORM_SLOPE = 0.06          # the storm series' slope now (NC state team)
CHOSEN = 0.10               # Hannah, 2026-11-23: the 1996 island value, rounded
MIN_DOMAINS = 75            # a survey with fewer domains measured is patchy at the water's edge
PERIODS = {"2-hindcast-1996-2009": ("hindcast", 1996, 2009),
           "3-test-2011-2020": ("test", 2011, 2020),
           "4-all-1996-2020": ("all", 1996, 2020)}


def summary(isl, ok):
    rows = []
    for label, sel in (("all surveys", isl), ("well-covered surveys", isl[ok])):
        rows.append(dict(surveys=label, n=len(sel), years=" ".join(str(y) for y in sel.index),
                         mean=round(sel.mean(), 3), median=round(sel.median(), 3),
                         min=round(sel.min(), 3), max=round(sel.max(), 3)))
    return pd.DataFrame(rows)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as S
    S.apply_style()
    d = pd.read_csv(MEAS / "foreshore_slope_domains.csv")
    for folder, (name, y0, y1) in PERIODS.items():
        out_dir = ROOT / folder
        out_dir.mkdir(parents=True, exist_ok=True)
        p = d[d.survey.between(y0, y1)]
        cover = p.groupby("survey").size()
        isl = p.groupby("survey")["median"].median()
        ok = cover.reindex(isl.index) >= MIN_DOMAINS
        good = isl.index[ok]
        pd_good = p[p.survey.isin(good)]
        by_dom = (pd_good.groupby("gis")["median"].agg(mean="mean", median="median", surveys="size")
                  .round(4).reset_index())
        by_dom.to_csv(out_dir / "slope_by_domain.csv", index=False)
        sm = summary(isl, ok)
        sm.to_csv(out_dir / "slope_summary.csv", index=False)
        print(f"\n{folder}\n{sm.to_string(index=False)}")
        mean_good, med_good = isl[ok].mean(), isl[ok].median()

        fig, (a, b) = plt.subplots(2, 1, figsize=(S.FIG_W_DOUBLE, 5.2), layout="constrained",
                                   gridspec_kw=dict(height_ratios=(1.6, 1)))
        a.plot(by_dom.gis, by_dom["mean"], color=S.C_1997, lw=1.8,
               label=f"Domain mean over {len(good)} surveys")
        a.axhline(mean_good, color=S.C_1997, lw=1.0, ls=(0, (1, 1)), label=f"Island mean, {mean_good:.3f}")
        a.axhline(CHOSEN, color="#1A9850", lw=1.2, ls=(0, (6, 2)), label=f"Chosen, {CHOSEN:.2f}")
        a.axhline(STORM_SLOPE, color="#D7301F", lw=1.2, ls=(0, (4, 2)), label=f"Storm series now, {STORM_SLOPE}")
        a.set_xlim(1, 90)
        a.set_ylim(0, 0.16)
        a.set_xlabel("Alongshore position (GIS domain, south → north)")
        a.set_ylabel("Foreshore slope (–)")
        a.legend(loc="lower left", frameon=False, fontsize=7.5, ncol=2)
        which = "every survey" if name == "all" else f"the {name} surveys"
        a.set_title(f"(a) Beach slope along the island, mean of {which} ({y0}–{y1})", loc="left")

        x = np.arange(len(isl))
        q = p.groupby("survey")["median"]
        lo, hi = q.quantile(.1).reindex(isl.index), q.quantile(.9).reindex(isl.index)
        b.vlines(x, lo, hi, colors=S.C_1997, lw=1.2)
        b.scatter(x[ok.to_numpy()], isl[ok], color=S.C_1997, s=22, zorder=3, label="Well covered")
        if (~ok).any():
            b.scatter(x[~ok.to_numpy()], isl[~ok], facecolor="white", edgecolor=S.C_1997, s=22, zorder=3,
                      label=f"Fewer than {MIN_DOMAINS} domains")
        b.axhline(mean_good, color=S.C_1997, lw=1.0, ls=(0, (1, 1)))
        b.axhline(CHOSEN, color="#1A9850", lw=1.2, ls=(0, (6, 2)))
        b.axhline(STORM_SLOPE, color="#D7301F", lw=1.2, ls=(0, (4, 2)))
        b.set_xticks(x, [f"{y}\n({cover[y]})" for y in isl.index], fontsize=7)
        b.set_ylim(0, 0.16)
        b.set_ylabel("Foreshore slope (–)")
        b.legend(loc="lower left", frameon=False, fontsize=7.5, ncol=2)
        b.set_title("(b) Island median by survey (domains measured), with the 10–90% range of domains", loc="left")

        out = out_dir / f"beach_slope_{name}.png"
        S.save(fig, out, close=True)
        dropped = ", ".join(str(y) for y in isl.index[~ok]) or "none"
        S.record_caption(out, (
            f"Foreshore slope, MHW (0.36 m NAVD88) to the berm (1.7 m NAVD88), on the 100 m transects, the {name} "
            f"surveys {y0}-{y1} (HAT_beach_slope_mhw_berm.py). (a) Each domain's mean over the well-covered surveys "
            f"(at least {MIN_DOMAINS} of 90 domains measured; left out: {dropped}) of its transect median (dark "
            f"blue); dotted, the island mean over those surveys ({mean_good:.3f}; median {med_good:.3f}); green "
            f"dashed, the value chosen for the storm series ({CHOSEN:.2f}, the 1996 island median rounded); red "
            f"dashed, the slope the storm series uses now ({STORM_SLOPE}, NC state team). (b) Each survey's island "
            f"median (dot; hollow, fewer than {MIN_DOMAINS} domains) and the 10-90% range of its domains (line); "
            "the number of domains measured under each year."))
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
