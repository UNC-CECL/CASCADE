"""
Does the model build the dunes up in its first years?

    python scripts/hatteras_ms/experiments/HAT_dune_crest_dem_vs_model.py

Reads the dune crest alongshore from the two starting DEMs (the 1996 and 2009
dune files Barrier3D starts from) and from the current matrix runs at model
year 0 and model year 3, and draws them on one alongshore axis with the dune
ceiling and the NC-12 rebuild height. Read-only: no
run is made. Also writes, per window and model year, the share of dune cells
lower than that year's highest run-up, on the model dunes and on the DEM dunes
held fixed.

Heights: the dune files hold height above the berm (dam). The berm is
BermEl = 1.7 m NAVD88 everywhere, so crest elevation (m NAVD88) = height
above berm + 1.7.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-19
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(_p for _p in _HERE.parents  if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
sys.path.insert(0, str(_HAT_DIR))
sys.path.insert(0, str(_HAT_DIR / "input_prep" / "8-overwash-analysis" / "4-vs-model"))
import storms_vs_overwash as SV  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hat_topo_version import product_for_year, topo_dirs  # noqa: E402
from site_layer.hatteras_site_config import first_model_year  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-19-dune-crest-dem-vs-model"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
ARM = "managed"          # the production runs; "natural" is in the table too
LATE_YEAR = 3            # model year compared with the DEM
MHW_NAVD88 = 0.36        # m, as the runner
DESIGN_NAVD88 = 4.3      # m, the NC-12 rebuild height (3.94 m MHW)
CELLS = 50               # 10 m cells per 500 m domain
COLOURS = {1996: fs.C_1984, 2009: fs.C_1997}   # calibration red, test blue
# -----------------------------------------------------------------------------


# Crest height above berm (m) per 10 m cell, GIS 1-90 south to north: (years, 4500)
def model_crests(c):
    return np.stack([np.asarray(c.barrier3d[SV.DOM.gis_to_pad(g)].DuneDomain, float).max(axis=2)
                     for g in range(1, 91)], axis=1).reshape(-1, 90 * CELLS) * 10


def dem_crests(start):
    _, dune_dir, _ = topo_dirs(product_for_year(start))
    return np.concatenate([np.load(dune_dir / f"domain_{g}_dune.npy") for g in range(1, 91)]) * 10


def load_window(window, storms):
    start = window[0]
    runs = {}
    for arm in ("managed", "natural"):
        name = SV.ovm.WINDOWS[window][arm]
        c = np.load(SV.ovm.MATRIX / f"{window[0]}_{window[1]}" / "edgeBE" / name / f"{name}.npz",
                    allow_pickle=True)["cascade"][0]
        runs[arm] = model_crests(c)
    b = c.barrier3d[SV.DOM.gis_to_pad(1)]
    berm_navd = b._BermEl * 10 + MHW_NAVD88
    ceiling_navd = b._Dmaxel * 10 + MHW_NAVD88
    s = storms[storms.pi == list(SV.ovm.WINDOWS).index(window)]
    years = [first_model_year(start) + k for k in range(LATE_YEAR)]
    yearly = s[s.calendar_year.isin(years)].groupby("calendar_year").rhigh_m.max() + MHW_NAVD88
    top = s[s.calendar_year.isin(years)].sort_values("rhigh_m").iloc[-1]
    return dict(start=start, window=window, dem=dem_crests(start), runs=runs, berm=berm_navd,
                ceiling=ceiling_navd, yearly_runup=yearly, top_storm=top)


# Share of cells whose crest at the start of each model year sits below that year's highest run-up
def exposure_table(w):
    rows = []
    for k, (year, runup) in enumerate(w["yearly_runup"].items()):
        above_berm = runup - w["berm"]
        row = dict(window=f"{w['window'][0]}_{w['window'][1]}", model_year=k + 1, calendar_year=year,
                   max_runup_m_navd88=round(runup, 2),
                   dem_share_below=round(float((w["dem"] < above_berm).mean()), 3))
        for arm, cr in w["runs"].items():
            row[f"{arm}_share_below"] = round(float((cr[k] < above_berm).mean()), 3)
        rows.append(row)
    return pd.DataFrame(rows)


def summary_table(w):
    rows = []
    ceiling_h = w["ceiling"] - w["berm"]
    states = [("DEM", w["dem"])] + [(f"{arm} year {t}", w["runs"][arm][t])
                                     for arm in ("managed", "natural") for t in range(LATE_YEAR + 1)]
    for label, h in states:
        rows.append(dict(window=f"{w['window'][0]}_{w['window'][1]}", state=label,
                         p10_m_above_berm=np.percentile(h, 10), median_m_above_berm=np.median(h),
                         p90_m_above_berm=np.percentile(h, 90),
                         share_cells_below_2m_above_berm=float((h < 2.0).mean()),
                         share_cells_at_ceiling=float((h >= ceiling_h - 0.05).mean())))
    return pd.DataFrame(rows).round(3)


# Median of each 500 m domain, drawn as a step across the domain
def domain_median(h):
    return np.repeat(np.median(h.reshape(90, CELLS), axis=1), CELLS)


def draw(window, path):
    fs.apply_style()
    x = np.arange(90 * CELLS) / CELLS + 0.5           # GIS domain units, cell centres
    fig, axes = plt.subplots(2, 1, figsize=fs.figsize("double", height=4.9), sharex=True,
                             gridspec_kw=dict(height_ratios=[3, 1.2], hspace=0.12))
    for i, w in enumerate([window]):
        ax, axd = axes[2 * i], axes[2 * i + 1]
        col, start, b = COLOURS[w["start"]], w["start"], w["berm"]
        late = w["runs"][ARM][LATE_YEAR]
        y_late = first_model_year(start) + LATE_YEAR
        ax.plot(x, w["dem"] + b, color="0.55", lw=0.35, alpha=0.6, zorder=2)
        ax.plot(x, late + b, color=col, lw=0.35, alpha=0.45, zorder=3)
        ax.plot(x, domain_median(w["dem"]) + b, color=fs.INK, lw=1.4,
                label=f"{start} DEM = model year 0 (500 m median)", zorder=4)
        ax.plot(x, domain_median(late) + b, color=col, lw=1.4,
                label=f"Model year {LATE_YEAR}, 1 Jan {y_late} (500 m median)", zorder=5)
        ax.axhline(w["ceiling"], color=fs.INK, lw=0.8, ls="--", label="Dune ceiling (Dmaxel)", zorder=1)
        ax.axhline(DESIGN_NAVD88, color=fs.C["REF"], lw=0.8, ls=":", label="NC-12 rebuild height", zorder=1)
        ax.set_ylim(b, 10)
        ax.set_ylabel("Dune crest (m NAVD88)")
        sec = ax.secondary_yaxis("right", functions=(lambda v, b=b: v - b, lambda v, b=b: v + b))
        sec.set_ylabel("Height above berm (m)")
        ax.grid(axis="y", color=fs.GRID_C, lw=0.5)
        ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=3, frameon=False, fontsize=6.5,
                  title=f"({'abcd'[2 * i]}) {start} start, {start}-{w['window'][1]} window; faint lines are 10 m cells",
                  title_fontsize=7.5, alignment="left", handlelength=2.2, columnspacing=1.2)

        d = late - w["dem"]
        axd.fill_between(x, 0, d, where=d > 0, color="#4393c3", lw=0, step="mid", label="Gain")
        axd.fill_between(x, 0, d, where=d < 0, color="#d6604d", lw=0, step="mid", label="Loss")
        axd.axhline(0, color=fs.INK, lw=0.5)
        axd.set_ylabel(f"Year {LATE_YEAR} − DEM (m)")
        axd.grid(axis="y", color=fs.GRID_C, lw=0.5)
        axd.legend(loc="lower right", ncol=2, frameon=False, fontsize=6.5)
        fs._letter_inside(axd, 2 * i + 1)
    for ax in axes:
        fs.open_frame(ax)
        ax.set_xlim(0.5, 90.5)
        fs.town_bands(ax, label=ax is axes[0])
        fs.structures(ax, label=ax is axes[0])
    axes[-1].set_xlabel(fs.DOMAIN_AXIS_LABEL + "; 500 m domains")
    fs.save(fig, path, close=True)
    fs.record_caption(path, (
        "Dune crest along Hatteras Island from the starting DEM (black/grey) and the "
        f"{ARM} matrix run at model year {LATE_YEAR} (colour); faint lines are the 10 m cells, bold lines "
        "the 500 m domain medians. Model year 0 equals the DEM to within 6 mm, so it is the black line. "
        "Left axis is crest elevation; right axis is height above the berm (1.7 m NAVD88), the "
        "quantity Barrier3D carries. Dashed: the dune ceiling (Dmaxel 5.5 m NAVD88, 3.8 m above "
        "the berm); dotted: NC-12 rebuild height (4.3 m NAVD88). "
        f"b: crest change per 10 m cell from the DEM to model year "
        f"{LATE_YEAR} (blue gain, red loss)."))


def main():
    forcing = SV.sf.load_forcing()
    storms = SV.sf.classify(SV.sf.load_record(SV.env.DEFAULT_STORM_VARIANT, forcing), forcing)
    windows = [load_window(w, storms) for w in SV.ovm.WINDOWS]
    for w in windows:
        gap = max(np.abs(w["runs"][arm][0] - w["dem"]).max() for arm in w["runs"])
        print(f"{w['start']}: model year 0 vs DEM, max |difference| {gap:.3f} m")
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    summ = pd.concat([summary_table(w) for w in windows], ignore_index=True)
    expo = pd.concat([exposure_table(w) for w in windows], ignore_index=True)
    summ.to_csv(EXP_DIR / "tables" / "dune_crest_summary.csv", index=False)
    expo.to_csv(EXP_DIR / "tables" / "cells_below_runup.csv", index=False)
    pd.concat([pd.DataFrame({"gis": np.repeat(np.arange(1, 91), CELLS), "cell": np.tile(np.arange(CELLS), 90),
                             f"dem_{w['start']}_m_navd88": (w["dem"] + w["berm"]).round(3),
                             **{f"{arm}_{w['start']}_year{t}_m_navd88": (w["runs"][arm][t] + w["berm"]).round(3)
                                for arm in ("managed", "natural") for t in (0, LATE_YEAR)}})
               .set_index(["gis", "cell"]) for w in windows], axis=1) \
        .to_csv(EXP_DIR / "tables" / "dune_crest_by_cell.csv")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(summ.to_string(index=False))
        print(expo.to_string(index=False))
    for w in windows:
        path = EXP_DIR / "figures" / f"dune_crest_dem_vs_model_year0_year3_{w['start']}_start.png"
        draw(w, path)
        print(f"figure: {path}")


if __name__ == "__main__":
    main()
