"""
How the collision hotspots were found, what lowered their dunes, and what the lowering looked like.

    python scripts/real-world-cases/Hatteras/dune-collision/1-observations/HAT_collision_hotspot_diagnostics.py

Three figures in the workspace's dune-collision/1-observations/2-hotspots/, beside collision_hotspot_crests.png,
from the same module-off controls, lidar and hotspot rule as HAT_collision_motivation_alongshore.py (imported, A):

    collision_hotspot_classification.png  every domain's mean overprediction against the crest lowering the model
                                          missed, per period, the rule's thresholds drawn and each domain classed
    collision_hotspot_storms.png          per hotspot, the lidar's and the model's crest change between successive
                                          surveys, the largest storm of each interval named
    collision_hotspot_profiles.png        per hotspot, one representative 10 m cell: the lidar profile the period
                                          started from and its last survey, the model's crest at both

Profiles: the start's own 1 m DEM clip (what the model started from) and the period's last lidar on the same
grid, as HAT_survey_dune_crests.py reads them (its vertical shifts applied), oriented ocean first; the median of
the cell's ten 1 m transects. The representative cell is the hotspot cell whose lidar crest loss, start to last
survey, is nearest the hotspot's median.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-03
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
if str(_HERE.parent) not in sys.path:
    sys.path.insert(0, str(_HERE.parent))
import HAT_collision_motivation_alongshore as A  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import FixedLocator, FormatStrFormatter, MultipleLocator  # noqa: E402

from site_layer.hat_figure_style import INK, INK_MUTED, apply_style, record_caption, save  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT_DIR = A.HOT_DIR
SUPPORT_DIR = OUT_DIR / "supporting"
# Storms named in the figures, by start (year, month) in the storm series; any other is named by its date
NAMED_STORMS = {(1998, 8): "Bonnie", (1999, 8): "Dennis", (1999, 9): "Floyd", (2003, 9): "Isabel",
                (2005, 9): "Ophelia", (2006, 9): "Ernesto", (2009, 11): "Nor'Ida", (2010, 9): "Earl",
                (2011, 8): "Irene", (2012, 10): "Sandy", (2014, 7): "Arthur", (2015, 10): "Joaquin",
                (2016, 10): "Matthew", (2017, 9): "Jose", (2018, 9): "Florence", (2019, 9): "Dorian",
                (2020, 8): "Isaias"}
LABEL_LOSS_M = 0.3              # a lidar crest loss at least this large gets its interval's storm named
PROFILE_HALF_M = (-70, 45)      # profile window, m landward (-) and seaward (+) of the start's crest
SPIKE_M = 20.0                  # 1 m returns above this are lidar spikes
CLASS_X = (-4.2, 4.2)           # the classification figure's x-range; a domain beyond is named at the edge
CLASS_C = {"hotspot": A.HIGH_C, "single": A.HIGH_C, "other": "0.35", "none": "0.80"}
START_C, END_C = "0.55", A.LIDAR_C
# -----------------------------------------------------------------------------


def survey_module():
    sys.path.insert(0, str(A._HAT_DIR / "input_prep" / "1-barrier3d-domains" / "2-survey-crests"))
    import HAT_survey_dune_crests as S
    return S


def storm_summary(d):
    run = A.PERIODS[d["key"]]["run"]
    meta = json.loads((run / f"{run.name}_run_metadata.json").read_text(encoding="utf-8"))
    sf = A.PROJECT_ROOT / "data" / "hatteras_init" / meta["identity"]["storm_file"]
    s = pd.read_csv(str(sf).replace(".npy", "_summary.csv"), parse_dates=["StartTime"])
    s["when"] = s.StartTime.dt.year + (s.StartTime.dt.dayofyear - 0.5) / 365.25
    s["rhigh_m_navd88"] = s.Rhigh * 10 + d["mhw"]
    s["name"] = [NAMED_STORMS.get((t.year, t.month), t.strftime("%b %Y")) for t in s.StartTime]
    return s


def classes(d, hots):
    """Per domain: the rule's numbers and its class (hotspot; single, a lone domain meeting the rule; other, too
    high without missed lowering; none)."""
    m, _, _ = A.collision_hotspots(d)
    in_hot = np.zeros(90, bool)
    for h in hots:
        in_hot[h["lo"] - 1:h["hi"]] = True
    m["class"] = np.select([in_hot, m["class"].eq("collision").to_numpy(), m["class"].eq("too_high_other").to_numpy()],
                           ["hotspot", "single", "other"], "none")
    return m


def period_grid(data, hots, height=1.95):
    keys = list(data)
    n = max(len(hots[k]) for k in keys)
    fig, axs = plt.subplots(n, 2, figsize=(7.48, height * n + 0.6), squeeze=False,
                            gridspec_kw=dict(width_ratios=[13, 16]))
    return fig, axs, keys, n


def letters(hots, keys, n):
    """Panel letters row by row (a, b across the first row ...), as collision_hotspot_crests.png."""
    out, it = {}, iter("abcdefghij")
    for i in range(n):
        for key in keys:
            if i < len(hots[key]):
                out[(i, key)] = next(it)
    return out


def column_titles(fig, axs, data, keys):
    for j, key in enumerate(keys):
        p = axs[0, j].get_position()
        fig.text((p.x0 + p.x1) / 2, 0.985, data[key]["label"], ha="center", va="top", fontsize=10.5, color=INK)


def year_axis(ax, d):
    lo_y = d["years"][0] - 1
    ax.set_xlim(lo_y - 0.3, d["years"][-1] + 0.3)
    ax.xaxis.set_major_locator(FixedLocator(range(int(np.ceil(lo_y / 2) * 2), int(d["years"][-1]) + 1, 2)))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))


# ---- 1: classification -------------------------------------------------------------------------------------------
def fig_classification(data, hots):
    fig, axs = plt.subplots(1, 2, figsize=(7.48, 3.6), sharex=True, sharey=True)
    rows = []
    for ax, letter, (key, d) in zip(axs, "ab", data.items()):
        m = classes(d, hots[key])
        ax.axvspan(A.HOT_OFF_M, 99, ymin=0, ymax=1, color="none")
        ax.fill_between([A.HOT_OFF_M, 99], A.HOT_DROP_M, 99, color=A.HOT_SHADE, lw=0, zorder=0)
        ax.axvline(A.HOT_OFF_M, color=INK_MUTED, lw=0.8, ls="--", zorder=1)
        ax.axhline(A.HOT_DROP_M, color=INK_MUTED, lw=0.8, ls="--", zorder=1)
        ax.axhline(0, color="0.85", lw=0.6, zorder=0)
        ax.axvline(0, color="0.85", lw=0.6, zorder=0)
        for cls, kw in (("none", dict(s=14)), ("other", dict(s=22)), ("single", dict(s=26, facecolor="white")),
                        ("hotspot", dict(s=26))):
            sel = m[m["class"] == cls]
            ax.scatter(sel.mean_m, sel.missed_loss_m, color=CLASS_C[cls], edgecolor=CLASS_C[cls], lw=0.8,
                       zorder=3 if cls != "none" else 2, **kw)
            if cls in ("other", "single"):
                for g, r in sel.iterrows():
                    ax.annotate(str(g), (r.mean_m, r.missed_loss_m), xytext=(3, 2), textcoords="offset points",
                                fontsize=6, color=INK if cls in ("hotspot", "single") else INK_MUTED,
                                path_effects=A.STRIP_HALO, zorder=4)
        for h in hots[key]:                                        # each hotspot named once, beside its points
            pts = m.loc[h["lo"]:h["hi"]]
            ax.annotate(A.span(h).replace("GIS ", ""), (pts.mean_m.max(), pts.missed_loss_m.max()), xytext=(5, 0),
                        textcoords="offset points", fontsize=6.5, fontweight="bold", color=A.HIGH_C, va="center",
                        path_effects=A.STRIP_HALO, zorder=5)
        off = m[m.mean_m < CLASS_X[0]]
        for g, r in off.iterrows():                                # off scale to the left: named at the edge
            ax.annotate(f"GIS {g} ({r.mean_m:.1f} m)", (CLASS_X[0], r.missed_loss_m), xytext=(4, 0),
                        textcoords="offset points", fontsize=6, color=INK_MUTED, va="center")
        ax.set_title(f"({letter})", loc="left", fontweight="bold")
        ax.set_title(d["label"], loc="center")
        ax.set_xlabel("Mean overprediction of the crest (m)")
        if letter == "a":
            ax.set_ylabel("Observed crest lowering not reproduced (m)")
        rows.append(m.assign(period=key))
    lim_y = max(np.nanmax(r.missed_loss_m) for r in rows) + 0.5
    axs[0].set_xlim(*CLASS_X)
    axs[0].set_ylim(min(np.nanmin(r.missed_loss_m) for r in rows) - 0.3, lim_y)
    fig.legend(handles=[
        Line2D([], [], marker="o", ls="none", color=CLASS_C["hotspot"], ms=5, label="In a collision hotspot"),
        Line2D([], [], marker="o", ls="none", mfc="white", mec=CLASS_C["single"], ms=5,
               label="Meets the rule, single domain"),
        Line2D([], [], marker="o", ls="none", color=CLASS_C["other"], ms=5, label="Overpredicted, no missed lowering"),
        Line2D([], [], marker="o", ls="none", color=CLASS_C["none"], ms=4, label="Neither")],
        loc="lower center", ncol=4, frameon=False, fontsize=7)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    out = OUT_DIR / "collision_hotspot_classification.png"
    save(fig, out, dpi=300, facecolor="white", close=True)
    pd.concat(rows).round(3).to_csv(SUPPORT_DIR / "collision_hotspot_classification.csv")
    record_caption(out, (
        "How the collision hotspots were found. Each point is a GIS domain: its mean overprediction of the dune "
        "crest (model minus lidar, averaged over the surveys after the period's start) against the crest lowering "
        "the model did not reproduce (the lidar's crest loss from the period's first survey to its last, less the "
        "model's over the same interval). Dashed, the rule's thresholds "
        f"({A.HOT_OFF_M:g} m each); shaded, where both are met. A collision hotspot also needs the overprediction "
        f"above {A.HOT_OFF_M:g} m in at least half the surveys, a lidar crest loss of at least {A.HOT_DROP_M:g} m, "
        "and two or more neighbouring domains meeting the rule (a one-domain gap bridged): red, in a hotspot; "
        "hollow, meets the rule alone; dark grey, persistently overpredicted without missed lowering (the model's "
        "dunes grew, or were rebuilt, faster than the lidar's; the collision module cannot correct these); light "
        "grey, neither. Numbers: GIS domains. Model: the module-off controls."))
    print(f"wrote {out}")


# ---- 2: storms between surveys --------------------------------------------------------------------------------------
def fig_storms(data, hots):
    fig, axs, keys, n = period_grid(data, hots)
    rows, let = [], letters(hots, keys, n)
    for j, key in enumerate(keys):
        d = data[key]
        storms = storm_summary(d)
        for i in range(n):
            ax = axs[i, j]
            if i >= len(hots[key]):
                ax.set_visible(False)
                continue
            h = hots[key][i]
            gs = range(h["lo"], h["hi"] + 1)
            s = d["surveys"][d["surveys"].gis.isin(gs)]
            per = []
            for sv, grp in s.groupby("survey"):
                state = int(grp.state_year.iloc[0]) - d["fmy"]
                per.append(dict(survey=sv, when=grp.when.iloc[0], lidar=grp.crest.median(),
                                model=np.median(d["crests"][state, grp.gis.to_numpy() - 1, grp.cell.to_numpy()])))
            # the start: the model's own first state is the start DEM's crest
            start = dict(survey=d["start_year"], when=float(d["start_year"]) + 0.6,
                         lidar=np.median(d["crests"][0, np.array(list(gs)) - 1, :]),
                         model=np.median(d["crests"][0, np.array(list(gs)) - 1, :]))
            p = pd.DataFrame([start] + [r for r in per if r["survey"] > d["start_year"]]).sort_values("when")
            for (_, a), (_, b) in zip(p.iloc[:-1].iterrows(), p.iloc[1:].iterrows()):
                width = b.when - a.when
                dl, dm = b.lidar - a.lidar, b.model - a.model
                ax.bar(a.when + width * 0.30, dl, width=width * 0.42, color=A.LIDAR_C, lw=0, zorder=3)
                ax.bar(a.when + width * 0.70, dm, width=width * 0.42, facecolor="white", edgecolor=A.MODEL_C,
                       lw=0.9, zorder=3)
                inside = storms[(storms.when > a.when) & (storms.when <= b.when)]
                big = inside.loc[inside.rhigh_m_navd88.idxmax()] if len(inside) else None
                if big is not None and dl <= -LABEL_LOSS_M:
                    ax.text(a.when + width * 0.30, dl - 0.12, big["name"], ha="center", va="top", fontsize=6,
                            color=INK, rotation=90 if width < 1.6 else 0, path_effects=A.STRIP_HALO, zorder=4)
                rows.append(dict(period=key, hotspot=A.span(h), from_survey=int(a.survey), to_survey=int(b.survey),
                                 lidar_change_m=dl, model_change_m=dm,
                                 largest_storm=None if big is None else big["name"],
                                 largest_storm_start=None if big is None else big.StartTime,
                                 largest_storm_rhigh_m_navd88=None if big is None else big.rhigh_m_navd88))
            ax.axhline(0, color=INK, lw=0.6, zorder=2)
            year_axis(ax, d)
            mine = [r["lidar_change_m"] for r in rows if r["hotspot"] == A.span(h) and r["period"] == key]
            lo = min(-1.0, np.nanmin(mine) - 0.9)                 # room under the deepest bar for its storm
            ax.set_ylim(lo, 1.2)
            ax.set_title(f"({let[(i, key)]}) {A.span(h)}", loc="left", fontweight="bold", fontsize=8.5)
            ax.set_title(f"missed lowering {h['missed_loss_m']:.1f} m", loc="right", fontsize=7, color=INK_MUTED)
            ax.set_xlim(right=p.when.max() + 0.8)                  # to the last survey (Feb 2020 in the test)
            if j == 0 or not axs[i, 0].get_visible():
                ax.set_ylabel("Crest change (m)")
    fig.legend(handles=[Patch(color=A.LIDAR_C, label="Lidar crest change (median over the hotspot's cells)"),
                        Patch(facecolor="white", edgecolor=A.MODEL_C, label="Model without the module, same interval")],
               loc="lower center", ncol=2, frameon=False, fontsize=7.5)
    fig.tight_layout(rect=(0, 0.035, 1, 0.965), h_pad=0.8)
    column_titles(fig, axs, data, keys)
    out = OUT_DIR / "collision_hotspot_storms.png"
    save(fig, out, dpi=300, facecolor="white", close=True)
    pd.DataFrame(rows).round(3).to_csv(SUPPORT_DIR / "collision_hotspot_storms.csv", index=False)
    record_caption(out, (
        "What lowered the dunes at each collision hotspot. Per interval between successive surveys (the first from "
        "the start DEM the model began with), the change in the median dune crest over the hotspot's 10 m cells: "
        "the lidar (filled) and the model without the collision module (open) at the 1 January each survey is "
        f"compared with. Where the lidar lost at least {LABEL_LOSS_M:g} m, the interval's largest storm by peak "
        "runup (R_high) in the model's storm series is named (or dated). Missed lowering, top right: the lidar's "
        "crest loss over the period less the model's, averaged over the hotspot's domains."))
    print(f"wrote {out}")


# ---- 4: profiles ----------------------------------------------------------------------------------------------------
def representative_cell(d, h):
    """The hotspot cell whose lidar crest loss (start crest to the last survey's) is nearest the hotspot's median."""
    last = max(d["surveys"].survey)
    s = d["surveys"][(d["surveys"].survey == last) & d["surveys"].gis.between(h["lo"], h["hi"])].copy()
    s["start"] = [d["crests"][0, g - 1, c] for g, c in zip(s.gis, s.cell)]
    s["loss"] = s.start - s.crest                          # crest: the berm where the dune became beach
    s = s[s.loss.notna()]
    r = s.iloc[(s.loss - s.loss.median()).abs().argmin()]
    return int(r.gis), int(r.cell), int(last), float(r.crest), str(r.status)


def cell_profile(start, g, cell, year, S):
    """The median of a cell's ten 1 m transects, ocean first: the start's DEM clip (year None) or a survey."""
    import rasterio
    prod, surveys = S.STARTS[start]
    clip = S.product(prod).gapfill_1m / f"clip_domain_{g}_filled.tif"
    grid = S.clip_grid(clip)
    if year is None:
        with rasterio.open(clip) as r:
            a = S.clean(r.read(1), r.nodata)
    else:
        kind, path, *_ = surveys[year]
        files = sorted(path.parent.glob(path.name)) if "*" in path.name else [path]
        srcs = [rasterio.open(f) for f in files]
        try:
            if kind == "mosaic":
                a = S.mosaic_on_clip(srcs[0], grid)
            else:
                unit = S.unit_of(srcs, start, {k: S.clip_grid(S.product(prod).gapfill_1m /
                                                                f"clip_domain_{k}_filled.tif") for k in (10, 30, 50, 70)})[0]
                a = S.raster_on_clip(srcs, grid) * (S.US_FT if unit == "ft" else 1.0)
            a = a + S.SURVEY_VERTICAL_OFFSET_M.get(start, {}).get(year, 0.0)
        finally:
            for f in srcs:
                f.close()
    a = S.X.orient_ocean_right(a, S.X.OCEAN_LOC)[:, ::-1]
    rows = a[10 * cell:10 * cell + 10]
    return np.nanmedian(np.where(rows < SPIKE_M, rows, np.nan), axis=0)


def fig_profiles(data, hots):
    S = survey_module()
    apply_style()                                                  # the survey script's import resets rcParams
    warnings.filterwarnings("ignore")
    fig, axs, keys, n = period_grid(data, hots, height=1.85)
    rows, let = [], letters(hots, keys, n)
    for j, key in enumerate(keys):
        d = data[key]
        start = A.PERIODS[key]["start"]
        for i in range(n):
            ax = axs[i, j]
            if i >= len(hots[key]):
                ax.set_visible(False)
                continue
            h = hots[key][i]
            g, cell, last, lid1, status = representative_cell(d, h)
            z0 = cell_profile(start, g, cell, None, S)
            z1 = cell_profile(start, g, cell, last, S)
            state = int(d["surveys"].query("survey == @last").state_year.iloc[0]) - d["fmy"]
            m0, m1 = d["crests"][0, g - 1, cell], d["crests"][state, g - 1, cell]
            x = np.arange(len(z0), dtype=float)
            # the window around the start's crest: the highest point within the first 150 m of land
            land = np.where(np.isfinite(z0) & (z0 > 1.5))[0]
            seek = slice(land[0], land[0] + 150) if land.size else slice(None)
            xc = int(np.nanargmax(np.where(np.isfinite(z0[seek]), z0[seek], -9)) + (seek.start or 0))
            xr = x - xc                                            # m from the start crest, seaward positive
            xr = -xr                                               # ocean first, so flip: seaward to the right
            keep = (xr >= PROFILE_HALF_M[0]) & (xr <= PROFILE_HALF_M[1])
            ax.plot(xr[keep], z0[keep], color=START_C, lw=1.3, label="start")
            ax.plot(xr[keep], z1[keep], color=END_C, lw=1.5, label="last")
            for mz, c in ((m0, START_C), (m1, A.MODEL_C)):
                ax.plot([-6, 6], [mz, mz], color=c, lw=3.2, solid_capstyle="butt", alpha=0.9, zorder=4)
            # the lidar's crest at the last survey, in the start's dune window (the berm where it became beach)
            ax.plot([-6, 6], [lid1, lid1], color=END_C, lw=1.2, ls=(0, (2, 1)), zorder=4)
            ax.annotate("", xy=(10, m1), xytext=(10, lid1),
                        arrowprops=dict(arrowstyle="<->", color=A.HIGH_C, lw=1.0, shrinkA=0, shrinkB=0), zorder=5)
            ax.text(12, (m1 + lid1) / 2, f"{m1 - lid1:+.1f} m", fontsize=6.5, color=A.HIGH_C, va="center",
                    path_effects=A.STRIP_HALO)
            if status != "ok":
                ax.text(0.02, 0.05, "dune now beach", transform=ax.transAxes, fontsize=6.5, color=INK_MUTED)
            ax.set_xlim(*PROFILE_HALF_M)
            ax.set_title(f"({let[(i, key)]}) {A.span(h)}: GIS {g}, cell {cell}", loc="left", fontweight="bold",
                         fontsize=8.5)
            ax.set_title(f"{d['start_year']} → {last}", loc="right", fontsize=7, color=INK_MUTED)
            if j == 0 or not axs[i, 0].get_visible():
                ax.set_ylabel("Elevation (m NAVD88)")
            if i == len(hots[key]) - 1:
                ax.set_xlabel("Distance from the start crest (m; ocean right)")
            rows.append(pd.DataFrame({"period": key, "hotspot": A.span(h), "gis": g, "cell": cell,
                                      "x_from_start_crest_m": xr[keep], "z_start_m_navd88": z0[keep],
                                      f"z_last_m_navd88": z1[keep], "last_survey": last,
                                      "model_crest_start_m": m0, "model_crest_last_m": m1}))
    fig.legend(handles=[Line2D([], [], color=START_C, lw=1.5, label="Lidar the period started from (the start DEM)"),
                        Line2D([], [], color=END_C, lw=1.5, label="Lidar, the period's last survey"),
                        Line2D([], [], color=START_C, lw=3.2, label="Model crest at the start"),
                        Line2D([], [], color=A.MODEL_C, lw=3.2, label="Model crest at the last survey"),
                        Line2D([], [], color=END_C, lw=1.2, ls=(0, (2, 1)), label="Lidar crest at the last survey")],
               loc="lower center", ncol=2, frameon=False, fontsize=7.5)
    fig.tight_layout(rect=(0, 0.06, 1, 0.965), h_pad=0.9)
    column_titles(fig, axs, data, keys)
    out = OUT_DIR / "collision_hotspot_profiles.png"
    save(fig, out, dpi=300, facecolor="white", close=True)
    pd.concat(rows).round(3).to_csv(SUPPORT_DIR / "collision_hotspot_profiles.csv", index=False)
    record_caption(out, (
        "What the missed crest lowering looked like at each collision hotspot. One representative 10 m cell per "
        "hotspot (the cell whose lidar crest loss is nearest the hotspot's median): the median of its ten 1 m "
        "cross-shore transects in the DEM the model started from (grey) and in the period's last lidar survey "
        "(blue), against distance from the starting crest, the ocean to the right. Bars: the model's crest "
        "without the collision module at the start (grey) and at the last survey (dark grey); dashed, the lidar's "
        "crest at the last survey, the highest point in the start's dune window (the berm where the dune has "
        "become beach); the arrow, the model's crest minus the lidar's."))
    print(f"wrote {out}")


def main():
    apply_style()
    SUPPORT_DIR.mkdir(parents=True, exist_ok=True)
    data = {k: A.load(k) for k in A.PERIODS}
    hots, _ = A.hotspots(data)
    fig_classification(data, hots)
    fig_storms(data, hots)
    fig_profiles(data, hots)


if __name__ == "__main__":
    main()
