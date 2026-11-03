"""
Where along the island, and when, the model without the collision module misses the lidar dunes.

    python scripts/real-world-cases/Hatteras/dune-collision/1-observations/HAT_collision_motivation_alongshore.py

Five figures beside collision_motivation.png, each with a calibration (1996-2009) and a test (2009-2025) part,
both from the module-off controls:

    collision_offset_heatmap.png     model minus lidar crest per GIS domain and survey; villages marked under the axis
    collision_hotspot_crests.png     crest through time at each hotspot, model band against lidar dots
    collision_alongshore_crests.png  crest along the island in three surveys per period
    collision_exposure_offset.png    storms that put each domain's dunes in collision, above its crest offset
    collision_beach_offset.png       each domain's beach width, above its crest offset

Model: cal_control_g16 (2026-10-26-dune-collision; Hs 3.0, StormStart 1, groin 0.16, module off) and the
2009-2025 test's control (2026-10-29-collision-test-2009-2025, module off, the 2014, 2017 and 2022 fills in).
Lidar: survey_dune_crests.csv of each start (HAT_survey_dune_crests.py), read and set against the model state
its model_state names as cascade_pipeline/plotting/dune_crest_comparison.py does (a cell whose dune window is
now beach counts as the berm). The 2009-start surveys are 2011 (post-Irene), 2012 (post-Sandy), 2014,
2016-2019 and 2020 (Feb; the last public lidar of the island, searched 2026-11-03).

Hotspots: the runs of neighbouring domains whose crest offset, averaged over the surveys after each start, exceeds
HOT_M in either period (gaps of one domain bridged), the HOT_N largest by summed offset.
Only too-high runs: the module lowers dunes, it cannot raise them.

Collision exposure: per cell and storm of the run's own storm series, the storm's peak runup (R_high, m NAVD88)
above the cell's dune toe (the module's toe, dune-toes/dune_toe_cells.csv) and below the model's crest at the
start of that storm's year (Sallenger's 2000 collision regime), counted over the run and averaged over the
domain's cells: the storms the model saw reach its dunes and, with the module off, did nothing about.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-03
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
for _p in (_HAT_DIR, _HERE.parent):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import FixedLocator, FormatStrFormatter, MultipleLocator  # noqa: E402

from cascade_pipeline.plotting import dune_crest_comparison as dc  # noqa: E402
from cascade_pipeline.saved_runs import load_saved_cascade  # noqa: E402
from site_layer import hat_topo_version as tv  # noqa: E402
from site_layer.hat_figure_style import (DOMAIN_AXIS_LABEL, INK, INK_MUTED, apply_style, record_caption,  # noqa: E402
                                         save, town_bands)
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS  # noqa: E402
import HAT_collision_motivation as M  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
# The topic folders under 1-observations/ (3-mechanism is HAT_collision_motivation.py's)
WHERE_DIR = M.ROOT_DIR / "1-mismatch-along-island"     # the heatmap and the crests along the island
HOT_DIR = M.ROOT_DIR / "2-hotspots"
DRIVER_DIR = M.ROOT_DIR / "4-drivers"                   # collision exposure and beach width
RUNS = M.RUNS
PERIODS = {
    "calibration": dict(label="Calibration, 1996–2009", start="1996-start",
                        run=RUNS / "2026-10-26-dune-collision" / "runs" / "cal_control_g16" / "1996_2009" / "edgeBE"
                        / "HAT_1996-2009_managed_groin_edgeBE"),
    "test": dict(label="Test, 2009–2025", start="2009-start",
                 run=RUNS / "2026-10-29-collision-test-2009-2025" / "runs" / "control" / "2009_2025" / "edgeBE"
                 / "HAT_2009-2025_managed_groin_edgeBE_nocollision"),
}
MIN_CELLS = 10                   # a domain-survey counts with this many surveyed cells
HOT_M, HOT_N = 1.0, 4            # hotspots: mean offset above this, the largest this many
# The three surveys per period drawn along the island
PROFILE_SURVEYS = {"calibration": (1999, 2005, 2009), "test": (2012, 2017, 2020)}
UNCLASSIFIED = {2001}
SURVEY_NOTES = {1999: "post-Dennis", 2011: "post-Irene", 2012: "post-Sandy", 2004: "post-Isabel"}
# Offset bands as the run overview's (dune_crest_comparison), RdBu with a neutral middle
OFFSET_BINS_M = dc.OFFSET_BINS_M
NO_SURVEY = dc.NO_SURVEY_COLOR
MODEL_C, LIDAR_C, HIGH_C, LOW_C = M.MODEL_C, M.LIDAR_C, "#B2182B", "#2166AC"
EXPOSURE_C = M.WATER_C
BEACH_C = "#B08D57"              # sand
HOT_SHADE = "#FBE3D6"
OFFSET_CLIP_M = 4.0              # the driver figures' offset axis; a bar beyond it is cut and labelled
HOT_DEF = ("the stretches of neighbouring GIS domains where the model's crests stand most above the lidar's: "
           "the domain's crest offset, averaged over a period's surveys after its start, above {m:g} m in either "
           "period (a one-domain gap bridged), the {n} largest by summed offset")
# -----------------------------------------------------------------------------


def load(key):
    p = PERIODS[key]
    run = p["run"]
    meta = json.loads((run / f"{run.name}_run_metadata.json").read_text(encoding="utf-8"))
    c = load_saved_cascade(run / f"{run.name}.npz")
    crests = dc.model_crests(c)                                   # (states, 90, cells), m NAVD88
    dom = HATTERAS_DOMAINS.for_total(len(c.barrier3d))
    b0 = c.barrier3d[dom.gis_to_pad(1)]
    mhw = b0._MHW * 10
    berm = (b0._BermEl + b0._MHW) * 10
    fmy = int(meta["period"]["first_model_year"])
    start_year = int(meta["period"]["start_year"])
    surveys = dc.survey_table(start_year, berm)
    table = dc.compare(crests, surveys, fmy)
    table = table[table.cells >= MIN_CELLS].copy()
    table["when"] = table.survey.map(surveys.groupby("survey").when.first())
    storms = dc.big_storms(PROJECT_ROOT / "data" / "hatteras_init" / meta["identity"]["storm_file"])
    toes = pd.read_csv(tv.product_dir(p["start"]) / "dune-toes" / "dune_toe_cells.csv").set_index(["gis", "cell"])
    ss = np.asarray(b0._StormSeries, float)
    return dict(key=key, label=p["label"], crests=crests, years=fmy + np.arange(crests.shape[0]), fmy=fmy,
                start_year=start_year, surveys=surveys, table=table, storms=storms, toes=toes, mhw=mhw,
                storm_series=ss)


def mean_offset(d):
    """Per GIS, the crest offset averaged over the surveys after the start."""
    t = d["table"]
    return t[t.survey > d["start_year"]].groupby("gis").model_minus_survey_m.mean().reindex(range(1, 91))


def span(h):
    return f"GIS {h['lo']}" if h["lo"] == h["hi"] else f"GIS {h['lo']}–{h['hi']}"


def hotspots(data):
    score = pd.concat([mean_offset(d) for d in data.values()], axis=1).max(axis=1)
    hot = (score > HOT_M).to_numpy()
    # bridge single-domain gaps
    for i in range(1, len(hot) - 1):
        if not hot[i] and hot[i - 1] and hot[i + 1]:
            hot[i] = True
    runs, i = [], 0
    while i < len(hot):
        if hot[i]:
            j = i
            while j + 1 < len(hot) and hot[j + 1]:
                j += 1
            g = list(range(i + 1, j + 2))
            runs.append((g[0], g[-1], float(score.loc[g].clip(lower=0).sum()), float(score.loc[g].mean())))
            i = j + 1
        else:
            i += 1
    runs = sorted(sorted(runs, key=lambda r: -r[2])[:HOT_N])
    return [dict(name="ABCDEFG"[k], lo=r[0], hi=r[1], summed_m=r[2], mean_m=r[3]) for k, r in enumerate(runs)], score


def exposure(d):
    """Per GIS: storms per cell in collision (toe < R_high < model crest at the start of the storm's year)."""
    ss, crests, toes = d["storm_series"], d["crests"], d["toes"]
    rows = []
    for g in range(1, 91):
        toe = np.array([toes.loc[(g, c), "toe_m_navd88"] if (g, c) in toes.index else np.nan
                        for c in range(crests.shape[2])])
        n_col = np.zeros(crests.shape[2])
        n_over = np.zeros(crests.shape[2])
        for t, r_high in zip(ss[:, 0].astype(int), ss[:, 1]):
            if not 1 <= t <= crests.shape[0] - 1:
                continue
            r = r_high * 10 + d["mhw"]
            crest = crests[t - 1, g - 1]
            n_col += (toe < r) & (r < crest)
            n_over += r >= crest
        ok = np.isfinite(toe)
        rows.append(dict(gis=g, collision_storms_per_cell=float(n_col[ok].mean()) if ok.any() else np.nan,
                         overtopping_storms_per_cell=float(n_over.mean())))
    return pd.DataFrame(rows).set_index("gis")


def beach_width(d):
    return d["toes"].groupby("gis").beach_width_m.median().reindex(range(1, 91))


def gis_axis(ax, label=True):
    ax.set_xlim(0.5, 90.5)
    ax.xaxis.set_major_locator(FixedLocator([1, 10, 20, 30, 40, 50, 60, 70, 80, 90]))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    if label:
        ax.set_xlabel(DOMAIN_AXIS_LABEL)


def shade_hotspots(ax, hots, label=False):
    for h in hots:
        ax.axvspan(h["lo"] - 0.5, h["hi"] + 0.5, color=HOT_SHADE, lw=0, zorder=0)
        if label:
            ax.text((h["lo"] + h["hi"]) / 2, 1.0, span(h).replace("GIS ", ""), transform=ax.get_xaxis_transform(),
                    ha="center", va="bottom", fontsize=6.5, color=INK)


def title(ax, letter, text):
    ax.set_title(f"({letter})", loc="left", fontweight="bold")
    ax.set_title(text, loc="center")


# ---- 1: heatmap ----------------------------------------------------------------------------------------------
def fig_heatmap(data, hots):
    colors = plt.get_cmap("RdBu_r", len(OFFSET_BINS_M) + 1)(np.arange(len(OFFSET_BINS_M) + 1))
    colors[len(colors) // 2] = (0.97, 0.97, 0.97, 1)
    cmap = ListedColormap(colors[1:-1])
    cmap.set_under(colors[0])
    cmap.set_over(colors[-1])
    cmap.set_bad(NO_SURVEY)
    norm = BoundaryNorm(OFFSET_BINS_M, cmap.N)
    n_rows = [d["table"].survey.nunique() for d in data.values()]
    fig, axs = plt.subplots(2, 1, figsize=(7.48, 6.6), sharex=True,
                            gridspec_kw=dict(height_ratios=n_rows, hspace=0.32))
    out = []
    for ax, (letter, d) in zip(axs, zip("ab", data.values())):
        piv = d["table"].pivot(index="survey", columns="gis", values="model_minus_survey_m").reindex(
            columns=range(1, 91))
        piv = piv.sort_index()
        im = ax.pcolormesh(np.arange(0.5, 91), np.arange(len(piv) + 1), np.ma.masked_invalid(piv.to_numpy()),
                           cmap=cmap, norm=norm, edgecolors="white", linewidth=0.3)
        ax.set_ylim(len(piv), 0)                                   # time runs down
        ax.set_yticks(np.arange(len(piv)) + 0.5)
        ax.set_yticklabels([f"{s}{'*' if s in UNCLASSIFIED else ''}"
                            + (f" {SURVEY_NOTES[s]}" if s in SURVEY_NOTES else "") for s in piv.index], fontsize=7)
        gis_axis(ax, label=letter == "b")
        ax.set_ylabel("Lidar survey")
        title(ax, letter, d["label"])
        out.append(piv.stack(dropna=False).rename("model_minus_lidar_m").reset_index().assign(period=d["key"]))
    # villages under the bottom panel's GIS ticks
    for name, (lo, hi) in HATTERAS_ANNOTATIONS.town_spans.items():
        for ax in axs:
            ax.plot([lo - 0.5, hi + 0.5], [-0.02, -0.02], transform=ax.get_xaxis_transform(), color=INK_MUTED,
                    lw=2.2, clip_on=False, solid_capstyle="butt")
        axs[1].text((lo + hi) / 2, -0.30, name, transform=axs[1].get_xaxis_transform(), ha="center", va="top",
                    fontsize=6.5, color=INK_MUTED)
    cb = fig.colorbar(im, ax=axs, orientation="horizontal", fraction=0.035, pad=0.13, extend="both",
                      ticks=OFFSET_BINS_M)
    cb.set_label("Model minus lidar crest, domain median (m): red, model dunes too high")
    cb.ax.tick_params(labelsize=7)
    out_png = WHERE_DIR / "collision_offset_heatmap.png"
    save(fig, out_png, dpi=300, facecolor="white", bbox_inches="tight", close=True)
    pd.concat(out).round(3).to_csv(WHERE_DIR / "supporting" / "collision_offset_heatmap.csv", index=False)
    record_caption(out_png, (
        "Where and when the modelled dune crests leave the lidar's, without the collision module. Each cell: the "
        "median modelled crest minus the median lidar crest over a GIS domain's surveyed 10 m cells (at least "
        f"{MIN_CELLS}), against the model state of the 1 January the survey is compared with; red, the model's "
        "dunes too high; white, within 0.5 m; grey, fewer surveyed cells. (a) Calibration, 1996–2009, the "
        "calibration control (Hs 3.0, StormStart 1, groin 0.16); (b) test, 2009–2025, the test control (the "
        "2014, 2017 and 2022 fills in). Both with the module off. A cell whose dune window is now beach counts as "
        "the berm. *2001, unclassified lidar. 2011 and 2012 were gridded from the ground returns at 1 m; 2020 "
        "(flown Feb 2020 at Buxton) is the last public lidar of the island. Bars under the axis: Buxton, Avon and the Tri-Village "
        "(Salvo, Waves, Rodanthe)."))
    print(f"wrote {out_png}")


# ---- 2: hotspot time series -------------------------------------------------------------------------------------
def fig_hotspots(data, hots):
    fig, axs = plt.subplots(len(hots), 2, figsize=(7.48, 1.95 * len(hots) + 0.6), sharey=True,
                            gridspec_kw=dict(width_ratios=[13, 16]))
    axs = np.atleast_2d(axs)
    rows = []
    for i, h in enumerate(hots):
        gs = range(h["lo"], h["hi"] + 1)
        for j, d in enumerate(data.values()):
            ax = axs[i, j]
            mod = np.concatenate([d["crests"][:, g - 1, :] for g in gs], axis=1)
            q = np.percentile(mod, [25, 50, 75], axis=1)
            ax.fill_between(d["years"], q[0], q[2], color=MODEL_C, alpha=0.22, lw=0)
            ax.plot(d["years"], q[1], color=MODEL_C, lw=1.8)
            s = d["surveys"][d["surveys"].gis.isin(gs)]
            for sv, grp in s.groupby("survey"):
                qq = grp.crest.quantile([0.25, 0.5, 0.75]).to_numpy()
                ax.errorbar(grp.when.iloc[0], qq[1], yerr=[[qq[1] - qq[0]], [qq[2] - qq[1]]], fmt="o", ms=4,
                            color=LIDAR_C, mfc="white" if sv in UNCLASSIFIED else LIDAR_C, elinewidth=1,
                            capsize=0, zorder=5)
                rows.append(dict(hotspot=span(h), gis=f"{h['lo']}-{h['hi']}", period=d["key"], survey=sv,
                                 lidar_p25=qq[0], lidar_median=qq[1], lidar_p75=qq[2]))
            for t, lab in d["storms"]:
                ax.axvline(t, color=INK_MUTED, lw=0.6, ls=":")
                if i == 0:
                    ax.text(t - 0.15, 0.97, lab, rotation=90, ha="right", va="top", fontsize=6, color=INK_MUTED,
                            transform=ax.get_xaxis_transform())
            lo_y = d["years"][0] - 1
            ax.set_xlim(lo_y - 0.3, d["years"][-1] + 0.3)
            ax.xaxis.set_major_locator(FixedLocator(range(int(np.ceil(lo_y / 2) * 2), int(d["years"][-1]) + 1, 2)))
            ax.xaxis.set_minor_locator(MultipleLocator(1))
            ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
            if i == 0:
                ax.set_title(d["label"], loc="center")
            if j == 0:
                ax.set_ylabel(f"{span(h)}\ncrest (m NAVD88)")
            ax.set_title(f"({'abcdefghij'[2 * i + j]})", loc="left", fontweight="bold")
    fig.legend(handles=[Line2D([], [], color=MODEL_C, lw=2, label="Model without the module (median; band 25th–75th)"),
                        Line2D([], [], color=LIDAR_C, marker="o", ls="none", label="Lidar (median; bar 25th–75th)")],
               loc="lower center", ncol=2, frameon=False, fontsize=7.5)
    fig.tight_layout(rect=(0, 0.035, 1, 1), h_pad=0.8)
    out = HOT_DIR / "collision_hotspot_crests.png"
    save(fig, out, dpi=300, facecolor="white", close=True)
    pd.DataFrame(rows).round(3).to_csv(HOT_DIR / "supporting" / "collision_hotspot_crests.csv", index=False)
    record_caption(out, (
        f"Dune crests at the hotspots, {HOT_DEF.format(m=HOT_M, n=HOT_N)}; calibration (left) and test (right). Line and "
        "band: the model without the collision module, median and 25th–75th percentile over the hotspot's 10 m "
        "cells, each 1 January; dots and bars: each lidar survey on the same cells, at the middle of the months "
        "flown (hollow, 2001 unclassified lidar; a cell whose dune window is now beach counts as the berm). Dotted: "
        "each period's three largest storms by R_high in the model's storm series. The test period has no lidar "
        "after Feb 2020."))
    print(f"wrote {out}")


# ---- 3: alongshore crests ---------------------------------------------------------------------------------------
def fig_alongshore(data, hots):
    n = max(len(v) for v in PROFILE_SURVEYS.values())
    fig, axs = plt.subplots(n, 2, figsize=(7.48, 1.9 * n + 0.6), sharex=True, sharey=True)
    rows = []
    for j, d in enumerate(data.values()):
        for i, sv in enumerate(PROFILE_SURVEYS[d["key"]]):
            ax = axs[i, j]
            s = d["surveys"][d["surveys"].survey == sv]
            if s.empty:
                ax.text(0.5, 0.5, f"no {sv} survey", transform=ax.transAxes, ha="center", color=INK_MUTED)
                continue
            state = int(s.state_year.iloc[0]) - d["fmy"]
            g_rows = []
            for g, grp in s.groupby("gis"):
                if len(grp) < MIN_CELLS:
                    continue
                m = d["crests"][state, g - 1, grp.cell.to_numpy()]
                g_rows.append(dict(gis=g, model_p25=np.percentile(m, 25), model_median=np.median(m),
                                   model_p75=np.percentile(m, 75), lidar_median=grp.crest.median()))
            t = pd.DataFrame(g_rows).set_index("gis").reindex(range(1, 91))
            ax.fill_between(t.index, t.model_p25, t.model_p75, color=MODEL_C, alpha=0.22, lw=0, step="mid")
            ax.plot(t.index, t.model_median, color=MODEL_C, lw=1.5, drawstyle="steps-mid")
            ax.plot(t.index, t.lidar_median, "o", ms=2.6, color=LIDAR_C, zorder=5)
            gis_axis(ax, label=i == n - 1)
            shade_hotspots(ax, hots, label=i == 0)
            note = f", {SURVEY_NOTES[sv]}" if sv in SURVEY_NOTES else ""
            ax.text(0.01, 0.96, f"{sv} lidar{note}", transform=ax.transAxes, ha="left", va="top", fontsize=7.5,
                    color=INK)
            ax.set_title(f"({'abcdefgh'[2 * i + j]})", loc="left", fontweight="bold")
            if i == 0:
                ax.set_title(d["label"], loc="center", pad=12)
            if j == 0:
                ax.set_ylabel("Crest (m NAVD88)")
            rows += [dict(period=d["key"], survey=sv, gis=g, **r) for g, r in t.dropna(how="all").iterrows()]
    axs[0, 0].set_ylim(1, 12)
    fig.legend(handles=[Line2D([], [], color=MODEL_C, lw=2, label="Model without the module (median; band 25th–75th)"),
                        Line2D([], [], color=LIDAR_C, marker="o", ls="none", ms=4, label="Lidar (domain median)"),
                        Patch(color=HOT_SHADE, label="Hotspot")],
               loc="lower center", ncol=3, frameon=False, fontsize=7.5)
    fig.tight_layout(rect=(0, 0.035, 1, 1), h_pad=0.6)
    out = WHERE_DIR / "collision_alongshore_crests.png"
    save(fig, out, dpi=300, facecolor="white", close=True)
    pd.DataFrame(rows).round(3).to_csv(WHERE_DIR / "supporting" / "collision_alongshore_crests.csv", index=False)
    record_caption(out, (
        "Dune crests along the island in three lidar surveys of each period: the model without the collision "
        "module (line, median over each GIS domain's surveyed 10 m cells; band, 25th–75th percentile) at the "
        "1 January the survey is compared with, and the lidar (dot, the domain's median), in domains with at least "
        f"{MIN_CELLS} surveyed cells. Calibration (left): {', '.join(map(str, PROFILE_SURVEYS['calibration']))}; "
        f"test (right): {', '.join(map(str, PROFILE_SURVEYS['test']))}. Shaded, the hotspots: "
        f"{HOT_DEF.format(m=HOT_M, n=HOT_N)}."))
    print(f"wrote {out}")


# ---- 4 and 5: a driver above the offset --------------------------------------------------------------------------
def fig_driver(data, hots, kind):
    fig, axs = plt.subplots(2, 2, figsize=(7.48, 4.6), sharex=True, sharey="row",
                            gridspec_kw=dict(height_ratios=[1, 1.15]))
    rows = []
    stats = []
    for j, d in enumerate(data.values()):
        off = mean_offset(d)
        if kind == "exposure":
            e = exposure(d)
            top = e.collision_storms_per_cell
            colour, ylab = EXPOSURE_C, "Storms in collision\nper cell"
        else:
            top = beach_width(d)
            colour, ylab = BEACH_C, "Dry beach width (m)"
        ax = axs[0, j]
        shade_hotspots(ax, hots, label=True)
        ax.bar(top.index, top.to_numpy(), width=0.8, color=colour, lw=0, zorder=2)
        if kind == "beach":
            ax.axhline(M.W0_M, color=INK, lw=0.8, ls="--", zorder=3)
            ax.text(90.3, M.W0_M, f"{M.W0_M:.0f} m", ha="right", va="bottom", fontsize=6.5, color=INK)
        ax.set_title(d["label"], loc="center", pad=18)
        ax.set_title(f"({'abcd'[j]})", loc="left", fontweight="bold", pad=18)
        if j == 0:
            ax.set_ylabel(ylab)
        ax = axs[1, j]
        shade_hotspots(ax, hots)
        v = off.to_numpy()
        ax.bar(off.index, np.clip(v, -OFFSET_CLIP_M, OFFSET_CLIP_M), width=0.8,
               color=np.where(v > 0, HIGH_C, LOW_C), lw=0, zorder=2)
        for g, x in off[off.abs() > OFFSET_CLIP_M].items():      # clipped bars carry their value
            ax.text(g - 0.8, np.sign(x) * (OFFSET_CLIP_M - 0.15), f"{x:.1f}", ha="right",
                    va="bottom" if x < 0 else "top", fontsize=6.5, color=INK)
        ax.set_ylim(-OFFSET_CLIP_M, OFFSET_CLIP_M)
        ax.axhline(0, color=INK, lw=0.7, zorder=3)
        ax.set_title(f"({'abcd'[2 + j]})", loc="left", fontweight="bold")
        gis_axis(ax)
        if j == 0:
            ax.set_ylabel("Model minus lidar\ncrest (m)")
        ok = top.notna() & off.notna()
        rho = pd.Series(top[ok]).rank().corr(pd.Series(off[ok]).rank())
        stats.append(f"{d['key']} Spearman r {rho:.2f} ({int(ok.sum())} domains)")
        rows.append(pd.DataFrame({"period": d["key"], "gis": top.index, kind: top.to_numpy(),
                                  "mean_offset_m": off.to_numpy()}))
    for ax in axs[0]:
        town_bands(ax, where="top", strip=0.09, label=True, fontsize=6)
    fig.tight_layout(h_pad=0.5)
    name = "collision_exposure_offset" if kind == "exposure" else "collision_beach_offset"
    out = DRIVER_DIR / f"{name}.png"
    save(fig, out, dpi=300, facecolor="white", close=True)
    pd.concat(rows).round(3).to_csv(DRIVER_DIR / "supporting" / f"{name}.csv", index=False)
    print(f"{name}: " + "; ".join(stats))
    off_txt = ("(c, d) Model minus lidar crest per GIS domain, averaged over the surveys after each period's start "
               "(the domain median of the surveyed cells, at least "
               f"{MIN_CELLS}); red, the model's dunes too high; a bar beyond ±{OFFSET_CLIP_M:g} m is cut and its value written. "
               f"Model: the module-off controls. Shaded, the "
               f"hotspots: {HOT_DEF.format(m=HOT_M, n=HOT_N)}. Rank correlation of the top row with the offset: "
               f"{'; '.join(stats)}.")
    if kind == "exposure":
        record_caption(out, (
            "How often the storms reached each domain's dunes without overtopping them. (a, b) Storms per 10 m cell "
            "in the collision regime (Sallenger 2000): the storm's peak runup in the run's own storm series "
            "(R_high, m NAVD88) above the cell's dune toe (the module's toe) and below the model's crest at the "
            "start of that storm's year; counted over the run, averaged over the domain's cells. Barrier3D lowers "
            "none of these dunes. " + off_txt))
    else:
        record_caption(out, (
            "The beach in front of each domain's dunes. (a, b) The dry beach width the collision module reads "
            "(beach start, 0.5 m above MHW, to the dune toe), the domain's median over its cells; dashed, "
            f"the {M.W0_M:.0f} m beyond which the module's beach factor is zero. " + off_txt))
    print(f"wrote {out}")


def main():
    apply_style()
    for d in (WHERE_DIR, HOT_DIR, DRIVER_DIR):
        (d / "supporting").mkdir(parents=True, exist_ok=True)
    data = {k: load(k) for k in PERIODS}
    hots, score = hotspots(data)
    print("hotspots:", [(h["name"], h["lo"], h["hi"], round(h["mean_m"], 2)) for h in hots])
    pd.DataFrame(hots).round(3).to_csv(HOT_DIR / "supporting" / "collision_hotspots.csv", index=False)
    fig_heatmap(data, hots)
    fig_hotspots(data, hots)
    fig_alongshore(data, hots)
    fig_driver(data, hots, "exposure")
    fig_driver(data, hots, "beach")


if __name__ == "__main__":
    main()
