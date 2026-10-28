"""
Where and when did the model overwash, and when was washover seen?

    python scripts/input_prep/8-overwash-analysis/4-vs-model/overwash_timeline.py

Years across and GIS domains up: the managed runs' overwash per model year as
shading, each image as a narrow column on its date (washover seen / none seen /
not assessed). A strip above counts domains per image and per model year; a strip
to the right gives each domain's share of images and of model years with overwash.
No scoring; overwash_vs_model.py does that. Details: scripts/input_prep/8-overwash-analysis/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-16
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in Path(__file__).resolve().parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as REPO  # noqa: E402
sys.path.insert(0, str(_HAT_DIR))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import overwash_vs_model as ovm  # noqa: E402
from site_layer import hat_overwash as ow  # noqa: E402
from site_layer import hat_env_forcings as env  # noqa: E402
from site_layer.hat_figure_style import (C_1984, C_1997, C_1984_FILL, C_1997_FILL, INK,  # noqa: E402
                                         INK_MUTED, apply_style, figsize, open_frame, record_caption,
                                         save, _title)
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
ARM = "managed"
BINS = [0.0, 1.0, 10.0, 50.0, 1e6]                    # m3/m; > 0 is overwash
BIN_LABELS = ["trace\n(< 1)", "1–10", "10–50", "> 50"]
QOW_COLOURS = ["#c7eae5", "#80cdc1", "#35978f", "#01665e"]   # one hue, light to dark
OBS_SEEN, OBS_NONE = INK, "0.80"                       # image columns; not assessed = hatched white
MOD_C = QOW_COLOURS[2]
COL_W = 0.22                                           # one column of an image's pair, years
BACK_ALPHA = 0.35                                      # the per-year shading behind the pairs
# Storms named along the top: (label, date of the peak). Hurricanes by name, nor'easters by month
STORMS = [("Bonnie", "1998-08-27", 0), ("Dennis", "1999-09-05", 1), ("Isabel", "2003-09-18", 2),
          ("Ophelia", "2005-09-14", 0), ("Nor'Ida", "2009-11-13", 1), ("Earl", "2010-09-03", 2),
          ("Irene", "2011-08-27", 0), ("Sandy", "2012-10-29", 1), ("Matthew", "2016-10-09", 0),
          ("Mar 2018 NE", "2018-03-04", 2), ("Florence", "2018-09-14", 1), ("Dorian", "2019-09-06", 0),
          ("Dec 2023 NE", "2023-12-18", 0)]
STORM_LEFT = {"Florence"}                              # labels set right of their line to clear a neighbour
# Places on the y axis: (GIS row it sits at, label)
PLACES = [(1, "Cape Point  1"), (7.5, "Buxton  7–8"), (26, "Avon  21–31"), (69, "Salvo  69"),
          (74, "Waves  74"), (80, "Rodanthe  80"), (90, "Pea Island  90")]
OUT = ow.VS_MODEL / "overwash_timeline_1997_2025.png"
# -----------------------------------------------------------------------------


def frac_year(t):
    t = pd.Timestamp(t)
    return t.year + (t.dayofyear - 1) / 365.25


# Each run's overwash per domain and model year: {first model year: (90, years) array}
def model_overwash():
    out, names = {}, []
    for window in ovm.WINDOWS:
        q, name = ovm.load_qow(window, ARM)
        pads = [DOM.for_total(q.shape[1]).gis_to_pad(g) for g in range(1, 91)]   # the run's own buffers, not the config's
        y0, y1 = ovm.model_years(window)
        out[y0] = q[1:y1 - y0 + 2][:, pads].T
        names.append(name)
    return out, names


# Per image and domain: observed, and the model's overwash since the previous image (overwash_vs_model's
# matching: model years dated by their largest storm, 7-day grace), from the run that covers the image
def pairs(obs):
    d = pd.concat([ovm.compare(w, ARM, obs, 0.0) for w in ovm.WINDOWS], ignore_index=True)
    d["date"] = pd.to_datetime(d.image)
    return d


# The header: one band per run, naming it
def panel_runs(ax, first, last):
    for window, fill, edge in zip(ovm.WINDOWS, (C_1984_FILL, C_1997_FILL), (C_1984, C_1997)):
        y0, y1 = ovm.model_years(window)
        ax.add_patch(Rectangle((y0, 0), y1 + 1 - y0, 1, fc=fill, ec="white", lw=1.5, alpha=0.55))
        ax.add_patch(Rectangle((y0, 0), y1 + 1 - y0, 0.12, fc=edge, lw=0))
        kind = "Calibration" if window == list(ovm.WINDOWS)[0] else "Test"
        ax.text((y0 + y1 + 1) / 2, 0.56, f"{kind} run · {window[0]} DEM · storms {y0}–{y1}",
                ha="center", va="center", fontsize=7, color=INK)
    ax.set_xlim(first, last + 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


# Storm names at their dates, over three rows (set in STORMS) so neighbours do not collide
def panel_storms(ax, first, last):
    for name, date, row in STORMS:
        x = frac_year(date)
        ax.plot([x, x], [0, 0.22 + 0.3 * row], color=INK_MUTED, lw=0.5)
        left = name in STORM_LEFT
        ax.text(x - 0.05 if left else x, 0.25 + 0.3 * row, name, ha="left" if left else "center",
                va="bottom", fontsize=6.5, color=INK)
    ax.set_xlim(first, last + 1)
    ax.set_ylim(0, 1.15)
    ax.axis("off")


# Domains overwashed per image: the model since the previous image (teal), the image (black)
def panel_counts(ax, pr, first, last):
    for date, g in pr.groupby("date"):
        x = frac_year(date)
        a = g.dropna(subset=["observed"])                    # only the domains the image assessed
        ax.bar(x - COL_W / 2, a.model.sum(), width=COL_W, color=MOD_C, lw=0, zorder=2)
        ax.bar(x + COL_W / 2, a.observed.sum(), width=COL_W, color=OBS_SEEN, lw=0, zorder=2)
    ax.set_xlim(first, last + 1)
    ax.set_ylim(0, 90)
    ax.set_yticks([0, 45, 90])
    ax.set_ylabel("Domains")
    ax.tick_params(labelbottom=False)
    ax.grid(axis="y", color="0.9", lw=0.5)
    open_frame(ax)


# The main panel: faded per-year shading, then each image as a pair (model since previous image | image)
def panel_main(ax, qow, pr, first, last, cmap, norm):
    for y0, v in qow.items():
        cls = np.digitize(v, BINS[1:-1]).astype(float)
        cls[v <= 0] = np.nan
        ax.pcolormesh(np.arange(y0, y0 + v.shape[1] + 1), np.arange(0.5, 91.5),
                      np.ma.masked_invalid(cls), cmap=cmap, norm=norm, alpha=BACK_ALPHA, zorder=1)
    for date, g in pr.groupby("date"):
        x = frac_year(date)
        g = g.set_index("gis")
        # left: the model since the previous image, overwashed or not (as the scoring reads it)
        ax.add_patch(Rectangle((x - COL_W, 0.5), COL_W, 90, fc="white", lw=0, zorder=2))
        for gis, m in g.model.items():
            if m:
                ax.add_patch(Rectangle((x - COL_W, gis - 0.5), COL_W, 1, fc=MOD_C, lw=0, zorder=3))
        # right: the image
        ax.add_patch(Rectangle((x, 0.5), COL_W, 90, fc="white", ec="0.6", hatch="/////", lw=0, zorder=2))
        for gis, val in g.observed.dropna().items():    # NaN = not assessed: the hatching shows through
            ax.add_patch(Rectangle((x, gis - 0.5), COL_W, 1, lw=0, zorder=3,
                                   fc=OBS_SEEN if val == 1 else OBS_NONE))
        ax.add_patch(Rectangle((x - COL_W, 0.5), 2 * COL_W, 90, fc="none", ec="0.55", lw=0.4, zorder=4))
    split = ovm.model_years(list(ovm.WINDOWS)[1])[0]
    ax.axvline(split, color=INK, lw=0.8, ls=":", zorder=4)
    ax.set_xlim(first, last + 1)
    ax.set_ylim(0.5, 90.5)
    ax.set_yticks([g for g, _ in PLACES], [lab for _, lab in PLACES], fontsize=7)
    ax.yaxis.set_minor_locator(MultipleLocator(5))
    ax.set_ylabel("Domain (GIS, south → north)")
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.set_xlabel("Year")
    # image dates as ticks along the top edge
    top = ax.secondary_xaxis("top")
    top.set_xticks(sorted({frac_year(d) for d in pr.date}))
    top.tick_params(length=3, labeltop=False, color=INK_MUTED)


# Per domain: share of images showing washover, share of model years with overwash
def panel_domains(ax, qow, imgs):
    gis = np.arange(1, 91)
    obs = imgs.groupby("domain").overwash.mean().reindex(gis) * 100
    v = np.concatenate(list(qow.values()), axis=1)
    mod = (v > 0).mean(axis=1) * 100
    ax.step(mod, gis, where="mid", color=MOD_C, lw=1.4, zorder=2)
    ax.step(obs.values, gis, where="mid", color=OBS_SEEN, lw=1.1, zorder=3)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel("Share (%)")
    ax.set_ylim(0.5, 90.5)
    ax.tick_params(labelleft=False)
    ax.grid(axis="x", color="0.9", lw=0.5)
    open_frame(ax)


def main():
    apply_style()
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    obs["date"] = pd.to_datetime(obs.Imagery_Date)
    qow, names = model_overwash()
    first = min(qow)
    last = max(y0 + v.shape[1] - 1 for y0, v in qow.items())
    imgs = obs[(obs.date > f"{first}-01-01") & (obs.date < f"{last + 1}-01-01")]
    pr = pairs(obs)
    cmap = ListedColormap(QOW_COLOURS)
    norm = BoundaryNorm(np.arange(len(BINS)), cmap.N)

    fig = plt.figure(figsize=figsize("double", height=7.6), constrained_layout=True)
    gs = fig.add_gridspec(4, 2, width_ratios=[5.2, 1], height_ratios=[0.32, 0.62, 1.2, 5.2])
    ax_r = fig.add_subplot(gs[0, 0])
    ax_s = fig.add_subplot(gs[1, 0])
    ax_main = fig.add_subplot(gs[3, 0])
    ax_c = fig.add_subplot(gs[2, 0], sharex=ax_main)
    ax_d = fig.add_subplot(gs[3, 1], sharey=ax_main)

    panel_runs(ax_r, first, last)
    panel_storms(ax_s, first, last)
    panel_counts(ax_c, pr, first, last)
    _title(ax_c, 0, "When: domains overwashed, per image")
    panel_main(ax_main, qow, pr, first, last, cmap, norm)
    _title(ax_main, 1, "Where and when")
    panel_domains(ax_d, qow, imgs)
    _title(ax_d, 2, "Where")

    # key: the model's stepped colour scale, then the observation and line keys
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cax = fig.add_subplot(gs[2, 1])
    cax.set_axis_off()
    inset = cax.inset_axes([0.05, 0.62, 0.9, 0.16])
    cb = fig.colorbar(sm, cax=inset, orientation="horizontal", ticks=np.arange(len(BINS) - 1) + 0.5)
    cb.ax.set_xticklabels(BIN_LABELS, fontsize=6)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(length=0)
    inset.set_title("Modelled overwash\n(m$^3$ m$^{-1}$ per model year)", fontsize=6.8, pad=3)
    handles = [Patch(fc=MOD_C, lw=0, label="Model overwashed since previous image (left of pair)"),
               Patch(fc=OBS_SEEN, lw=0, label="Image: washover seen"),
               Patch(fc=OBS_NONE, lw=0, label="Image: none seen"),
               Patch(fc="white", ec="0.6", hatch="/////", lw=0.4, label="Image: not assessed"),
               Patch(fc=MOD_C, alpha=BACK_ALPHA, lw=0, label="Model, per year (background)"),
               Line2D([], [], color=OBS_SEEN, lw=1.1, label="Images with washover (c)"),
               Line2D([], [], color=MOD_C, lw=1.4, label="Model years with overwash (c)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, frameon=False, fontsize=6.8,
               handlelength=1.4, columnspacing=1.0)

    out = save(fig, OUT, close=True)
    n_img = imgs.date.nunique()
    partial = " and ".join(f"{d:%b %Y}" for d in sorted(pr[pr.partial].date.unique()))
    record_caption(out[0],
        f"Where and when the model overwashed, against the {n_img} images of the overwash record, storm years "
        f"{first}–{last}. Managed runs, groin on ({', '.join(names)}; storm series {env.DEFAULT_STORM_VARIANT}). "
        "(b) Each image is a framed pair on its date. Left half: teal where the model overwashed the domain "
        "(any overwash > 0) since the previous image, matched as overwash_vs_model scores it: a model year counts toward the first image after its largest storm, with "
        "the 7-day grace of the observed storm table. Right half: the image, black where it shows washover, grey "
        "where it was assessed and shows none, hatched where the domain was not assessed. So a storm and the "
        "image that recorded it sit side by side however late the imagery came (Isabel, Sept 2003, pairs with "
        f"May 2004). For the first image of each run ({partial}) the interval starts before the run, so the model "
        "covers only part of it. Faded behind "
        "the pairs: the model's overwash per model year, filling its calendar year because Barrier3D records "
        "overwash per year, not per storm. Ticks on the top edge mark the image dates; the dotted line is where "
        "the calibration run (1996 DEM) gives way to the test run (2009 DEM), named in the header band. Storms "
        "named above: the hurricanes and nor'easters (NE) with the highest water levels in the model's storm "
        "series, at their peak dates. (a) For each image, the domains the model overwashed since the previous "
        "image (teal) and the domains the image shows overwashed (black), both over the domains the image "
        "assessed. (c) Each domain's share of images showing washover (black) and of model years with any "
        "overwash (teal). GIS 1 is Cape Point, GIS 90 Pea Island.")
    print(out[0])


if __name__ == "__main__":
    main()
