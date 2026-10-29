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
COLOURS = {"Buxton": "#1B9E77", "Avon": "#D95F02", "Salvo": "#7570B3", "Waves": "#E7298A", "Rodanthe": "#66A61E"}
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

    # Draw the crests along the island with each village's cap, and each village's spread
    per = sv.groupby("gis").crest_m_navd88.agg(median="median", p95=lambda a: a.quantile(0.95),
                                               p05=lambda a: a.quantile(0.05))
    fig, axs = plt.subplots(2, 1, figsize=(12, 8.6), gridspec_kw={"height_ratios": [1.15, 1]})
    ax = axs[0]
    ax.fill_between(per.index, per.p05, per.p95, color="0.85", lw=0, label="lidar cells, 5th-95th percentile")
    ax.plot(per.index, per["median"], color="0.35", lw=1.4, label="lidar, domain median crest")
    ax.axhline(ISLAND_CAP_M, color="#B2182B", lw=1.4, ls="--", label=f"current island-wide cap {ISLAND_CAP_M:g} m")
    for name, (lo, hi) in VILLAGES.items():
        ax.axvspan(lo - 0.5, hi + 0.5, color=COLOURS[name], alpha=0.12, lw=0)
        ax.hlines(caps[name], lo - 0.5, hi + 0.5, color=COLOURS[name], lw=3)
        ax.text((lo + hi) / 2, 9.35, name, ha="center", fontsize=8, color=COLOURS[name], fontweight="bold")
        ax.text((lo + hi) / 2, caps[name] + 0.15, f"{caps[name]:g}", ha="center", fontsize=8, color=COLOURS[name])
    ax.plot([], [], color="0.3", lw=3, label=f"per-village cap: lidar {CAP_QUANTILE * 100:.0f}th percentile")
    ax.set_xlim(0.5, 90.5)
    ax.set_ylim(0, 9.8)
    ax.set_xlabel("GIS domain (south → north)")
    ax.set_ylabel("Dune crest (m NAVD88)")
    ax.legend(fontsize=7.5, loc="lower center", ncol=2)
    ax.set_title("(a)", loc="left", fontsize=9)
    ax = axs[1]
    names = list(VILLAGES)
    data = [sv[sv.gis.between(*VILLAGES[n])].crest_m_navd88.dropna().to_numpy() for n in names]
    parts = ax.violinplot(data, showextrema=False, widths=0.8)
    for body, n in zip(parts["bodies"], names):
        body.set_facecolor(COLOURS[n])
        body.set_alpha(0.35)
    for i, n in enumerate(names, start=1):
        ax.hlines(np.median(data[i - 1]), i - 0.25, i + 0.25, color="0.25", lw=1.5)
        ax.hlines(caps[n], i - 0.4, i + 0.4, color=COLOURS[n], lw=3)
        ax.text(i + 0.42, caps[n], f"cap {caps[n]:g}", va="center", fontsize=8, color=COLOURS[n])
    ax.axhline(ISLAND_CAP_M, color="#B2182B", lw=1.4, ls="--")
    ax.set_xticks(range(1, len(names) + 1))
    ax.set_xticklabels([f"{n}\nGIS {VILLAGES[n][0]}-{VILLAGES[n][1]}" for n in names])
    ax.set_ylabel("Lidar dune crest, all cells (m NAVD88)")
    ax.set_ylim(0, 9.8)
    ax.set_title("(b)", loc="left", fontsize=9)
    fig.tight_layout()

    # Write the figure and its caption
    out = fs.figure_dir("inputs", "4-management", "village_dune_cap") / "village_dune_caps.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, facecolor="white")
    fs.record_caption(out, (
        "The village dune cap per village, from the lidar. (a) Dune crest along the island: the median of each "
        f"domain's spike-screened lidar cells over the {START} surveys ({', '.join(map(str, SURVEYS_EXCLUDED))} "
        "left out) and their 5th-95th percentile band; coloured bars, each village's cap, the "
        f"{CAP_QUANTILE * 100:.0f}th percentile of its cells, measured as the island-wide {ISLAND_CAP_M:g} m was "
        "(dashed). (b) Each village's lidar crests (violins, median line) with its cap. Tri-Village split at the "
        "midpoints between village centres: Salvo 68-71, Waves 72-76, Rodanthe 77-83."))
    print(f"wrote {out} and {OUT_DIR / 'village_dune_caps.csv'}")


if __name__ == "__main__":
    main()
