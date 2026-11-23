"""
The beach slope along the island in 1996 and the island value through the surveys (2026-11-23).

    python scripts/real-world-cases/Hatteras/input_prep/10-beach-slope/HAT_beach_slope_figure.py

Reads (D:, hatteras_init/10-beach-slope/1-foreshore-mhw-to-berm/): foreshore_slope_domains.csv
(HAT_beach_slope_mhw_berm.py). The island value is the median of the 1996 domain medians: one slope for the storm
series, from the survey the model starts from.
Writes: island_slope.csv (the island value per survey), beach_slope_1996_and_surveys.png (caption in CAPTIONS.md).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
"""
import sys
from pathlib import Path

_HAT_DIR = next(p for p in Path(__file__).resolve().parents if (p / "hatteras_paths.py").exists())
sys.path.insert(0, str(_HAT_DIR))

import pandas as pd
from hatteras_paths import WORKSPACE

DIR = WORKSPACE / "data" / "hatteras_init" / "10-beach-slope" / "1-foreshore-mhw-to-berm"
STORM_SLOPE = 0.06          # the storm series' slope (NC state team)
START_SURVEY = 1996
SEASON = {1996: "Aug-Oct", 1997: "Sep-Oct", 2000: "Feb-Mar", 2001: "Jan-Mar", 2004: "Apr-Sep", 2005: "Aug-Nov",
          2009: "Aug"}      # the flights (hat_lidar_strips.MONTHS)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as S
    S.apply_style()
    d = pd.read_csv(DIR / "foreshore_slope_domains.csv")
    isl = (d.groupby("survey")["median"]
           .agg(island_median="median", p10=lambda x: x.quantile(.1), p90=lambda x: x.quantile(.9), domains="size")
           .round(3).reset_index())
    isl["flown"] = isl.survey.map(SEASON)
    isl.to_csv(DIR / "island_slope.csv", index=False)
    value = float(isl.loc[isl.survey == START_SURVEY, "island_median"].iloc[0])
    print(isl.to_string(index=False))

    fig, (a, b) = plt.subplots(2, 1, figsize=(S.FIG_W_DOUBLE, 5.2), layout="constrained",
                               gridspec_kw=dict(height_ratios=(1.6, 1)))
    s = d[d.survey == START_SURVEY].sort_values("gis")
    a.plot(s.gis, s["median"], color=S.C_1997, lw=1.8, label=f"{START_SURVEY} lidar, domain median")
    a.axhline(value, color=S.C_1997, lw=1.0, ls=(0, (1, 1)), label=f"Island value, {value:.3f}")
    a.axhline(STORM_SLOPE, color="#D7301F", lw=1.2, ls=(0, (4, 2)), label=f"Storm series now, {STORM_SLOPE}")
    a.set_xlim(1, 90)
    a.set_ylim(0, 0.16)
    a.set_xlabel("Alongshore position (GIS domain, south → north)")
    a.set_ylabel("Foreshore slope (–)")
    a.legend(loc="upper right", frameon=False, fontsize=7.5, ncol=3)
    a.set_title("(a) Beach slope along the island, 1996", loc="left")

    x = range(len(isl))
    cols = [S.C_1997 if y == START_SURVEY else "0.45" for y in isl.survey]
    b.vlines(x, isl.p10, isl.p90, colors=cols, lw=1.2)
    b.scatter(x, isl.island_median, color=cols, s=22, zorder=3)
    b.axhline(STORM_SLOPE, color="#D7301F", lw=1.2, ls=(0, (4, 2)))
    b.set_xticks(list(x), [f"{y}\n{SEASON.get(y, '')}" for y in isl.survey], fontsize=7.5)
    b.set_ylim(0, 0.16)
    b.set_ylabel("Foreshore slope (–)")
    b.set_title("(b) Island median by survey, with the 10–90% range of domains", loc="left")

    out = DIR / "beach_slope_1996_and_surveys.png"
    S.save(fig, out, close=True)
    S.record_caption(out, (
        "Foreshore slope, MHW (0.36 m NAVD88) to the berm (1.7 m NAVD88), measured on the 100 m transects "
        "(HAT_beach_slope_mhw_berm.py). (a) The 1996 lidar, the median of each domain's transects (dark blue); dotted, "
        f"the island value, the median of the 1996 domain medians ({value:.3f}); red dashed, the slope the storm series "
        f"uses now ({STORM_SLOPE}, from the NC state team). (b) The island median for every survey (dot) and the 10-90% "
        "range of the domain medians (line), with the months flown; 1996 in dark blue. The winter surveys (2000, 2001) "
        "are flatter, as winter beaches are; the 2009 grid is a different product (USACE DEM) with fewer usable "
        "domains."))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
