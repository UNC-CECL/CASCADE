"""
The village dune cap per village, measured from the lidar: the beach/dune manager's ceiling, one per village.

    python scripts/real-world-cases/Hatteras/input_prep/4-mgmt-forcings/village_dune_cap/HAT_village_dune_cap.py

The beach/dune manager clips each village dune cell to its cap every year, and the overwash it bulldozes
onto the dunes builds up to it. One island-wide cap of 8.0 m NAVD88 (hindcast config
village_dune_cap_navd88_m) was measured as the 95th percentile of the village dune cells. Per village the
lidar differs: Waves and Avon reach about 8 m, Salvo about 6.7 m, and Rodanthe's dunes are half as tall as
its neighbours'. Per village: the 95th percentile of its spike-screened lidar crest cells over the surveys
(2001 left out, as for 8.0 m), the cap; the 90th percentile and the median beside it. Rodanthe is also given
from GIS 78, since GIS 77 still carries Waves-like dunes. A measurement only: the runner still reads the one
config value.

Reads data/hatteras_init/1-barrier3d-domains/<start>/survey-crests/survey_dune_crests.csv. Writes
data/hatteras_init/4-mgmt-forcing/village_dune_cap/village_dune_caps.csv and draws
output/figures/3-model-inputs/4-management/village_dune_cap/village_dune_caps.png.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-29
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in sys.path:
    sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
import site_layer.hat_figure_style as fs  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
START = "1996-start"
SURVEY_CRESTS = (PROJECT_ROOT / "data" / "hatteras_init" / "1-barrier3d-domains" / START / "survey-crests"
                 / "survey_dune_crests.csv")
OUT_DIR = PROJECT_ROOT / "data" / "hatteras_init" / "4-mgmt-forcing" / "village_dune_cap"
SURVEYS_EXCLUDED = (2001,)           # as for the island-wide 8.0 m
CAP_QUANTILE = 0.95                  # the cap: a ceiling the tallest real dunes reach
ISLAND_CAP_M = 8.0                   # the current island-wide cap, m NAVD88
# The villages (HATTERAS_COMMUNITY_ZONES); the Tri-Village split at the midpoints between village centres
VILLAGES = {"Buxton": (7, 8), "Avon": (21, 31), "Salvo": (68, 71), "Waves": (72, 76), "Rodanthe": (77, 83)}
RODANTHE_ALT = (78, 83)              # without GIS 77's Waves-like dunes
# -----------------------------------------------------------------------------


# One row of statistics for the crest cells of GIS lo-hi
def village_row(sv, name, lo, hi):
    x = sv[sv.gis.between(lo, hi)].crest_m_navd88
    return dict(village=name, first_gis=lo, last_gis=hi, cells=len(x),
                p95_m_navd88=x.quantile(0.95), p90_m_navd88=x.quantile(0.90), median_m_navd88=x.median(),
                cap_m_navd88=round(float(x.quantile(CAP_QUANTILE)), 1))


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fs.apply_style()

    sv = pd.read_csv(SURVEY_CRESTS)
    sv = sv[(sv.status == "ok") & ~sv.survey.isin(SURVEYS_EXCLUDED)]

    # Measure each village, and Rodanthe from GIS 78
    rows = [village_row(sv, n, lo, hi) for n, (lo, hi) in VILLAGES.items()]
    rows.append(village_row(sv, "Rodanthe (from 78)", *RODANTHE_ALT))
    t = pd.DataFrame(rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t.round(2).to_csv(OUT_DIR / "village_dune_caps.csv", index=False)
    print(t.round(2).to_string(index=False))
    caps = dict(zip(t.village, t.cap_m_navd88))

    # Draw the crests along the island with each village's cap, and each village's spread (house paper style)
    per = sv.groupby("gis").crest_m_navd88.agg(median="median", p95=lambda a: a.quantile(0.95),
                                               p05=lambda a: a.quantile(0.05))
    fig, axs = plt.subplots(2, 1, figsize=fs.figsize("double", aspect=0.72), layout="constrained",
                            gridspec_kw={"height_ratios": [1.25, 1]})
    ylim = (0, 10.4)
    ax = axs[0]
    ax.set_xlim(0.5, 90.5)
    ax.set_ylim(*ylim)
    band = ax.fill_between(per.index, per.p05, per.p95, color=fs.C["BASE_FILL"], lw=0, zorder=1)
    (line,) = ax.plot(per.index, per["median"], color=fs.INK, lw=1.0, zorder=3)
    island = ax.axhline(ISLAND_CAP_M, color=fs.INK_MUTED, lw=0.9, ls="--", zorder=2)
    for name, (lo, hi) in VILLAGES.items():
        ax.hlines(caps[name], lo - 0.5, hi + 0.5, color=fs.C["ACCENT"], lw=2.2, zorder=4)
        ax.text(hi + 0.9, caps[name], f"{caps[name]:.1f}", va="center", ha="left", fontsize=7,
                color=fs.C["ACCENT"], zorder=5, bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
    # Villages as a strip along the top, drawn over the data: a full-height wash would merge with the band
    before = set(ax.patches) | set(ax.texts)
    fs.town_bands(ax, spans=VILLAGES, fontsize=7, strip=0.07, shade="0.82")
    for artist in (set(ax.patches) | set(ax.texts)) - before:
        artist.set_zorder(6)
    ax.set_xlabel(fs.DOMAIN_AXIS_LABEL)
    ax.set_ylabel("Dune crest elevation (m NAVD88)")
    fs._title(ax, 0, "Lidar dune crests along the island")
    ax = axs[1]
    names = list(VILLAGES)
    data = [sv[sv.gis.between(*VILLAGES[n])].crest_m_navd88.dropna().to_numpy() for n in names]
    ax.boxplot(data, whis=(5, 95), showfliers=False, widths=0.5, patch_artist=True,
               boxprops=dict(facecolor=fs.C["BASE_FILL"], edgecolor=fs.INK, lw=0.8),
               medianprops=dict(color=fs.INK, lw=1.2), whiskerprops=dict(color=fs.INK, lw=0.8),
               capprops=dict(color=fs.INK, lw=0.8))
    ax.axhline(ISLAND_CAP_M, color=fs.INK_MUTED, lw=0.9, ls="--", zorder=1)
    for i, n in enumerate(names, start=1):
        ax.hlines(caps[n], i - 0.36, i + 0.36, color=fs.C["ACCENT"], lw=2.2, zorder=4)
        ax.text(i + 0.4, caps[n], f"{caps[n]:.1f}", va="center", ha="left", fontsize=7, color=fs.C["ACCENT"],
                zorder=5, bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
    ax.set_xticks(range(1, len(names) + 1))
    ax.set_xticklabels([f"{n}\n(GIS {VILLAGES[n][0]}–{VILLAGES[n][1]})" for n in names])
    ax.set_xlim(0.4, len(names) + 0.75)
    ax.set_ylim(0, 9.5)
    ax.set_ylabel("Dune crest elevation (m NAVD88)")
    ax.grid(axis="y", color=fs.C["GRID"], lw=0.5)
    ax.set_axisbelow(True)
    fs._title(ax, 1, "Crest distribution by village")

    # One frameless legend under both panels
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [line, band, Patch(facecolor=fs.C["BASE_FILL"], edgecolor=fs.INK, lw=0.8),
               Line2D([], [], color=fs.C["ACCENT"], lw=2.2), island]
    labels = ["Domain median", "5th–95th percentile of cells",
              "Village cells: median, quartiles, 5th–95th percentile",
              f"Village cap: {CAP_QUANTILE * 100:.0f}th percentile", f"Island-wide cap ({ISLAND_CAP_M:.1f} m)"]
    fig.legend(handles, labels, loc="outside lower center", ncol=3, frameon=False, fontsize=8)

    # Write the figure (PNG, and a PDF under supporting/) and its caption
    out = fs.figure_dir("inputs", "4-management", "village_dune_cap") / "village_dune_caps.png"
    fs.save(fig, out, dpi=300, facecolor="white")
    fs.record_caption(out, (
        "Village dune caps measured from the lidar. (a) Dune crest elevation along the island (GIS 1 at Cape "
        "Point, 90 at the Pea Island boundary): the median of each domain's spike-screened lidar crest cells over "
        f"the {START} surveys ({', '.join(map(str, SURVEYS_EXCLUDED))} excluded), with the 5th-95th percentile of "
        f"the cells shaded. Purple bars: each village's cap, the {CAP_QUANTILE * 100:.0f}th percentile of its cells, "
        f"measured as the island-wide {ISLAND_CAP_M:.1f} m cap was (dashed). Villages named in the strip along the top. (b) The "
        "crest cells of each village: box, median and quartiles; whiskers, 5th-95th percentile; purple, the "
        "village cap. The Tri-Village is split at the midpoints between village centres (Salvo 68-71, Waves 72-76, "
        "Rodanthe 77-83)."))
    print(f"wrote {out} and {OUT_DIR / 'village_dune_caps.csv'}")


if __name__ == "__main__":
    main()
