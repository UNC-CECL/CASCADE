"""
The road's dune rebuild rule, three versions: where NC 12's dunes get rebuilt, against where the lidar shows they were.

    python scripts/real-world-cases/Hatteras/input_prep/4-mgmt-forcings/road_dune_rebuild/HAT_road_dune_rebuild.py

The roadway manager rebuilds a domain's dunes to the design height when its front dune row drops below a trigger,
but only where the road sits no farther behind the dune than a gate (obx-cascade's dune_mgmt_road_setback_trigger).
Three gates: the hindcast's 20 m; 150 m, just past the farthest road at a lidar rebuild; and NCDOT's own criterion,
dune maintenance where NC 12 is within 650 ft (198 m) of the shoreline (Bonner Bridge Replacement FEIS, NCDOT/FHWA
2008, p. 2-128, from the landward reach of Hurricane Isabel's washover fans), which on the manager's axis, the road
behind the dune, is 198 m less the domain's beach width. The same FEIS: minor dune maintenance 2-3 times a year,
most often where the dune is within 25 ft (7.6 m) of the road (p. 2-144, 2-128); CAMA Major Permit 86-99 (renewed
2021) authorises it on Hatteras in mapped areas.

A rebuild in the lidar: a domain's median spike-screened crest up >= 1 m at >= 0.4 m/yr between surveys (natural
growth here is under ~0.1 m/yr), outside the village zones (the beach/dune manager's) and where there is a road.
Surveys: the 1996- and 2009-start survey crests (1998 and 2001 left out, as for the toes and caps) and the 2011
post-Irene lidar read the same way (dune-storm-response/2026-10-29-collision-test-2009-2025, tables/
post_irene_2011_crests.csv). A measurement only: the runner still reads the one road_dune_rebuild_trigger_m.

Reads the road setbacks (4-mgmt-forcing/road_offset, RoadSetback_<year>_dunestart.csv), the beach widths
(1-barrier3d-domains/<start>/dune-toes/dune_toe_cells.csv) and the survey crests. Writes
data/hatteras_init/4-mgmt-forcing/road_dune_rebuild/road_dune_rebuild_gates.csv (per domain and start: setback,
beach width, the three gates, rebuild allowed under each) and observed_dune_rebuilds.csv, and draws
output/figures/3-model-inputs/4-management/road_dune_rebuild/road_dune_rebuild_rules.png.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-30
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
from site_layer.hat_topo_version import road_setback_file  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
INIT = PROJECT_ROOT / "data" / "hatteras_init"
DOMAINS = INIT / "1-barrier3d-domains"
POST_IRENE = (PROJECT_ROOT / "output" / "raw_runs" / "dune-storm-response" / "2026-10-29-collision-test-2009-2025"
              / "tables" / "post_irene_2011_crests.csv")
OUT_DIR = INIT / "4-mgmt-forcing" / "road_dune_rebuild"
STARTS = (1996, 2009)
SURVEYS_EXCLUDED = (1998, 2001)
VILLAGES = {"Buxton": (7, 8), "Avon": (21, 31), "Salvo": (68, 71), "Waves": (72, 76), "Rodanthe": (77, 83)}
JUMP_M, RATE_M_YR = 1.0, 0.4          # a rebuild in the lidar
GATE_NOW_M, GATE_FIT_M = 20.0, 150.0  # the hindcast's gate; just past the farthest road at a lidar rebuild (140 m)
SHORE_GATE_M = 198.0                  # NCDOT: NC 12 within 650 ft of the shoreline (Bonner Bridge FEIS 2008, p. 2-128)
MHW_M = 0.36                          # MHW, m NAVD88
BERM_MHW_M, MIN_DUNE_M = 1.34, 0.3    # Barrier3D berm (m MHW) and the roadway manager's absolute minimum dune
DESIGN_M_NAVD88 = 3.94 + MHW_M        # the hindcast's rebuild height
# -----------------------------------------------------------------------------


# Road setback behind the dune per GIS (m) from the 2-row RoadSetback file; nan where no road
def setbacks(year):
    a = np.loadtxt(road_setback_file(year), delimiter=",")
    s = pd.Series(a[1], index=a[0].astype(int))
    return s.where(s > 0)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fs.apply_style()

    # Per start: setback, beach width, the gates, and whether each rule lets the road manager rebuild there
    rows = []
    for y in STARTS:
        sb = setbacks(y)
        w = pd.read_csv(DOMAINS / f"{y}-start" / "dune-toes" / "dune_toe_cells.csv").groupby("gis").beach_width_m.median()
        for g in range(1, 91):
            s, bw = sb.get(g, np.nan), w.get(g, np.nan)
            shore = max(SHORE_GATE_M - bw, 0.0) if np.isfinite(bw) else np.nan
            rows.append(dict(start=y, gis=g, road_setback_m=s, beach_width_m=bw, gate_now_m=GATE_NOW_M,
                             gate_fit_m=GATE_FIT_M, gate_ncdot_m=shore,
                             rebuild_now=bool(s <= GATE_NOW_M) if np.isfinite(s) else False,
                             rebuild_fit=bool(s <= GATE_FIT_M) if np.isfinite(s) else False,
                             rebuild_ncdot=bool(s <= shore) if np.isfinite(s) and np.isfinite(shore) else False))
    gates = pd.DataFrame(rows)

    # Rebuilds in the lidar
    frames = []
    for y in STARTS:
        t = pd.read_csv(DOMAINS / f"{y}-start" / "survey-crests" / "survey_dune_crests.csv")
        t = t[(t.status == "ok") & ~t.survey.astype(int).isin(SURVEYS_EXCLUDED)]
        frames.append(t[["survey", "gis", "crest_m_navd88"]].assign(survey=t.survey.astype(int)))
    if POST_IRENE.exists():
        p = pd.read_csv(POST_IRENE)
        frames.append(p[p.status == "ok"].assign(survey=2011)[["survey", "gis", "crest_m_navd88"]])
    cr = pd.concat(frames).groupby(["gis", "survey"]).crest_m_navd88.median().unstack()
    years = list(cr.columns)
    reb = []
    for g in cr.index:
        if any(lo <= g <= hi for lo, hi in VILLAGES.values()):
            continue
        for a, b in zip(years[:-1], years[1:]):
            x0, x1 = cr.loc[g, a], cr.loc[g, b]
            if np.isfinite(x0) and np.isfinite(x1) and x1 - x0 >= JUMP_M and (x1 - x0) / (b - a) >= RATE_M_YR:
                st = 1996 if b <= 2009 else 2009
                row = gates[(gates.start == st) & (gates.gis == g)].iloc[0]
                if np.isfinite(row.road_setback_m):
                    reb.append(dict(gis=g, from_survey=a, to_survey=b, crest_before_m_navd88=x0, crest_after_m_navd88=x1,
                                    road_setback_m=row.road_setback_m, gate_ncdot_m=row.gate_ncdot_m))
    reb = pd.DataFrame(reb)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    gates.round(1).to_csv(OUT_DIR / "road_dune_rebuild_gates.csv", index=False)
    reb.round(2).to_csv(OUT_DIR / "observed_dune_rebuilds.csv", index=False)
    print(reb.round(2).to_string(index=False))
    for y in STARTS:
        g = gates[(gates.start == y) & gates.road_setback_m.notna()]
        print(f"{y} start: road domains {len(g)}; rebuild allowed now {g.rebuild_now.sum()}, 150 m {g.rebuild_fit.sum()}, "
              f"NCDOT 198 m {g.rebuild_ncdot.sum()}")
    print(f"lidar rebuilds inside the gate: now {(reb.road_setback_m <= GATE_NOW_M).sum()}, 150 m "
          f"{(reb.road_setback_m <= GATE_FIT_M).sum()}, NCDOT {(reb.road_setback_m <= reb.gate_ncdot_m).sum()} of {len(reb)}")

    # Figure: (a) the road behind the dune and the three gates along the island, the lidar rebuilds on it;
    # (b) the crest before and after each rebuild against the triggers and the design height
    fig, axs = plt.subplots(2, 1, figsize=fs.figsize("double", aspect=0.78), layout="constrained",
                            gridspec_kw={"height_ratios": [1.5, 1]})
    ax = axs[0]
    g09, g96 = (gates[gates.start == y].set_index("gis") for y in (2009, 1996))
    ax.fill_between(g09.index, 0, g09.road_setback_m, step="mid", color=fs.C["BASE_FILL"], lw=0, zorder=1)
    (road09,) = ax.step(g09.index, g09.road_setback_m, where="mid", color=fs.INK, lw=0.9, zorder=3)
    (road96,) = ax.step(g96.index, g96.road_setback_m, where="mid", color=fs.INK_MUTED, lw=0.8, ls=":", zorder=3)
    now = ax.axhline(GATE_NOW_M, color=fs.C["BASE"], lw=1.2, ls="--", zorder=4)
    fit = ax.axhline(GATE_FIT_M, color=fs.C["ACCENT"], lw=1.4, zorder=4)
    (shore,) = ax.step(g09.index, g09.gate_ncdot_m, where="mid", color="#e66101", lw=1.6, zorder=4)
    obs = ax.scatter(reb.gis, reb.road_setback_m, s=46, marker="*", color="#2166ac", edgecolor="white", lw=0.5, zorder=6)
    before = set(ax.patches) | set(ax.texts)
    fs.town_bands(ax, spans=VILLAGES, fontsize=7, strip=0.07, shade="0.82")
    for artist in (set(ax.patches) | set(ax.texts)) - before:
        artist.set_zorder(7)
    ax.set_xlim(0.5, 90.5)
    ax.set_ylim(0, 460)
    ax.set_xlabel(fs.DOMAIN_AXIS_LABEL)
    ax.set_ylabel("NC 12 behind the dune (m)")
    fs._title(ax, 0, "Where the road manager may rebuild the dunes")
    ax.legend([road09, road96, now, fit, shore, obs],
              ["Road behind the dune, 2009 start", "Road behind the dune, 1996 start", f"Gate now ({GATE_NOW_M:g} m)",
               f"Gate fitted to the lidar ({GATE_FIT_M:g} m)",
               f"NCDOT: road within {SHORE_GATE_M:g} m (650 ft) of the shoreline", "Rebuild in the lidar"],
              loc="upper left", bbox_to_anchor=(0.0, 0.92), fontsize=7, frameon=False, ncol=2)

    ax = axs[1]
    reb = reb.sort_values(["gis", "to_survey"]).reset_index(drop=True)
    x = np.arange(len(reb))
    for i, r in reb.iterrows():
        ax.annotate("", xy=(i, r.crest_after_m_navd88), xytext=(i, r.crest_before_m_navd88),
                    arrowprops=dict(arrowstyle="-|>", color="#2166ac", lw=1.3, shrinkA=0, shrinkB=0))
    ax.scatter(x, reb.crest_before_m_navd88, s=18, color="#2166ac", zorder=5)
    trig_now = BERM_MHW_M + MIN_DUNE_M + MHW_M
    ax.axhline(trig_now, color=fs.C["BASE"], lw=1.2, ls="--")
    ax.axhline(reb.crest_before_m_navd88.median(), color="#e66101", lw=1.4)
    ax.axhline(DESIGN_M_NAVD88, color=fs.INK, lw=1.0, ls="-.")
    ax.set_xticks(x)
    ax.set_xticklabels([f"GIS {r.gis}\n{r.from_survey}–{str(r.to_survey)[2:]}" for r in reb.itertuples()], fontsize=6.5)
    ax.set_xlim(-0.6, len(reb) - 0.4)
    ax.set_ylim(1, 6.5)
    ax.set_ylabel("Dune crest (m NAVD88)")
    ax.grid(axis="y", color=fs.C["GRID"], lw=0.5)
    ax.set_axisbelow(True)
    fs._title(ax, 1, "Crest before and after each rebuild in the lidar")
    ax.legend([Line2D([], [], color="#2166ac", marker="o", ms=4, lw=1.3), Line2D([], [], color=fs.C["BASE"], ls="--"),
               Line2D([], [], color="#e66101", lw=1.4), Line2D([], [], color=fs.INK, ls="-.")],
              ["Domain median crest, before → after", f"Trigger now (berm floor, {trig_now:.2f} m)",
               f"Trigger measured (median before, {reb.crest_before_m_navd88.median():.2f} m)",
               f"Design height ({DESIGN_M_NAVD88:.2f} m)"], loc="upper right", fontsize=7, frameon=False, ncol=2)

    out = fs.figure_dir("inputs", "4-management", "road_dune_rebuild") / "road_dune_rebuild_rules.png"
    fs.save(fig, out, dpi=300, facecolor="white")
    fs.record_caption(out, (
        "The road's dune rebuild rule. (a) NC 12's distance behind the dune per domain (GIS 1 at Cape Point, 90 at "
        "the Pea Island boundary; the model's road setback inputs, 2009 start shaded, 1996 dotted) against three "
        f"gates within which the roadway manager may rebuild: the hindcast's {GATE_NOW_M:g} m; {GATE_FIT_M:g} m, "
        "just past the farthest road at a lidar rebuild; and NCDOT's criterion, NC 12 within 650 ft (198 m) of the "
        "shoreline (Bonner Bridge Replacement FEIS 2008, p. 2-128), here 198 m less each domain's beach width. "
        "Stars: dune rebuilds in the lidar (domain-median crest up at least 1 m at 0.4 m/yr or more between "
        "surveys, village zones excluded). (b) The crest before and after each rebuild against the hindcast's "
        "trigger (the berm floor), the median crest before a rebuild, and the design height."))
    print(f"wrote {out} and {OUT_DIR}")


if __name__ == "__main__":
    main()
