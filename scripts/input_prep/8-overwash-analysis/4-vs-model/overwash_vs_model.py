"""
Does the model overwash where and when the imagery shows washover?

    python scripts/input_prep/8-overwash-analysis/4-vs-model/overwash_vs_model.py

The observed record against Barrier3D's overwash in the 1996-2010 and 2010-2024
hindcast runs, each image against the storms since the previous one. Details: scripts/input_prep/8-overwash-analysis/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-29
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import ListedColormap, BoundaryNorm  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))

from site_layer import hat_overwash as ow  # noqa: E402
from site_layer import hat_env_forcings as env  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, C, DOMAIN_AXIS_LABEL, figsize, save, record_caption,
    _title, open_frame, town_bands,
)
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
MATRIX = REPO / "output" / "raw_runs" / "matrix"
WINDOWS = {
    (1996, 2010): {"managed": "HAT_1996_2010_edgeBE_offsetmetres_road_bdm_nogroin",
                   "natural": "HAT_1996_2010_edgeBE_offsetmetres_noroad_nobdm_nogroin"},
    (2010, 2024): {"managed": "HAT_2010_2024_edgeBE_offsetmetres_road_bdm_nourish_nogroin",
                   "natural": "HAT_2010_2024_edgeBE_offsetmetres_noroad_nobdm_nogroin"},
}
THRESHOLDS = (0.0, 1.0, 5.0)        # m3/m
# -----------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "1-observations"))
from overwash_data import CAPTURE_GRACE_DAYS  # noqa: E402
GRACE = pd.Timedelta(days=CAPTURE_GRACE_DAYS)   # the observed storm table's own rule
HEADLINE_THR = 0.0

# cell classes for the matrix figure
BOTH, OBS_ONLY, MOD_ONLY, NEITHER, UNASSESSED = 0, 1, 2, 3, 4
CLASS_COLOURS = ["0.2", C["ACCENT"], C["ADDED"], "0.93", "white"]
CLASS_LABELS = ["both", "observed only", "model only", "neither", "not assessed"]


# A run's overwash per model year and domain (QowTS, m3/m), and its name
def load_qow(window, arm):
    name = WINDOWS[window][arm]
    path = MATRIX / f"{window[0]}_{window[1]}" / "edgeBE" / name / f"{name}.npz"
    c = np.load(path, allow_pickle=True)["cascade"][0]
    q = np.array([np.asarray(b.QowTS) for b in c.barrier3d]).T      # (nt, pads) m3/m, row t = model year t
    return q, name


# Per model year: the end date of its largest storm, and of every storm
def storm_dates(window):
    f = env.storm_summary_file(*window)          # the series the matrix runs on
    s = pd.read_csv(f, parse_dates=["StartTime", "EndTime"])
    s["dated"] = s.EndTime - GRACE
    largest = s.loc[s.groupby("time").Rhigh.idxmax()].set_index("time").dated
    every = s.groupby("time").dated.apply(list)
    return largest, every, f


# Images whose storms the run covers, each with its comparison window
def images_in(window, obs):
    start = pd.Timestamp(f"{window[0]}-01-01")
    end = pd.Timestamp(f"{window[1]}-01-01")        # the last modelled storm year ends here
    imgs = (obs[["Obs_ID", "Imagery_Date"]].drop_duplicates()
            .assign(date=lambda d: pd.to_datetime(d.Imagery_Date)).sort_values("date")
            .reset_index(drop=True))
    imgs["prev"] = imgs.date.shift(1)
    keep = imgs[(imgs.date > start) & (imgs.date <= end)].copy()
    keep["from"] = keep.prev.where(keep.prev >= start, start)
    keep["partial"] = keep.prev < start
    return keep.reset_index(drop=True)


# Observed against modelled overwash per image and domain, for one window, arm and threshold
def compare(window, arm, obs, thr):
    q, name = load_qow(window, arm)
    largest, every, _ = storm_dates(window)
    imgs = images_in(window, obs)
    rows = []
    for _, im in imgs.iterrows():
        o = obs[obs.Obs_ID == im.Obs_ID].set_index("domain").overwash
        for gis in range(1, 91):
            p = DOM.gis_to_pad(gis)
            yrs = [t for t in range(1, q.shape[0]) if q[t, p] > thr]
            m_big = any(im["from"] < largest[t] <= im.date for t in yrs if t in largest)
            m_any = any(im["from"] < d <= im.date for t in yrs if t in every for d in every[t])
            ob = o.get(gis, np.nan)
            rows.append(dict(window=f"{window[0]}-{window[1]}", arm=arm, run=name, threshold=thr,
                             obs_id=im.Obs_ID, image=im.date.date(), since=im["from"].date(),
                             partial=bool(im.partial), gis=gis, observed=ob,
                             model=int(m_big), model_any_storm=int(m_any),
                             qow_in_window=float(sum(q[t, p] for t in range(1, q.shape[0])
                                                     if t in largest and im["from"] < largest[t] <= im.date))))
    return pd.DataFrame(rows)


# Hits, misses, false alarms and correct negatives, per cell and per domain
def scores(df, model_col="model"):
    a = df.dropna(subset=["observed"])
    o, m = a.observed.astype(int), a[model_col].astype(int)
    hits, miss = int(((o == 1) & (m == 1)).sum()), int(((o == 1) & (m == 0)).sum())
    fa, cn = int(((o == 0) & (m == 1)).sum()), int(((o == 0) & (m == 0)).sum())
    ever = a.groupby("gis").agg(o=("observed", "max"), m=(model_col, "max"))
    return dict(cells=len(a), both=hits, observed_only=miss, model_only=fa, neither=cn,
                hit_rate=hits / max(hits + miss, 1), model_hits_confirmed=hits / max(hits + fa, 1),
                domains_obs_ever=int((ever.o == 1).sum()), domains_model_ever=int((ever.m == 1).sum()),
                domains_both_ever=int(((ever.o == 1) & (ever.m == 1)).sum()))


# A cell's agreement class for the matrix figure
def cell_class(o, m):
    if np.isnan(o):
        return UNASSESSED
    return {(1, 1): BOTH, (1, 0): OBS_ONLY, (0, 1): MOD_ONLY, (0, 0): NEITHER}[(int(o), int(m))]


# The image x domain agreement matrix for one window (managed, headline threshold)
def figure(df, window):
    d = df[(df.window == f"{window[0]}-{window[1]}") & (df.arm == "managed") & (df.threshold == HEADLINE_THR)]
    imgs = d[["obs_id", "image", "since", "partial"]].drop_duplicates().reset_index(drop=True)
    grid = np.full((len(imgs), 90), UNASSESSED)
    for i, im in imgs.iterrows():
        sub = d[d.obs_id == im.obs_id].set_index("gis")
        for g in range(1, 91):
            grid[i, g - 1] = cell_class(sub.observed[g], sub.model[g])
    cmap = ListedColormap(CLASS_COLOURS)
    norm = BoundaryNorm(np.arange(-0.5, 5.5), cmap.N)

    fig = plt.figure(figsize=figsize("double", height=2.2 + 0.24 * len(imgs) + 1.8), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 0.18], height_ratios=[0.24 * len(imgs) + 0.6, 1.5])
    ax = fig.add_subplot(gs[0, 0])
    ax.imshow(grid, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest",
              extent=(0.5, 90.5, len(imgs) - 0.5, -0.5))
    ax.set_yticks(range(len(imgs)))
    ax.set_yticklabels([f"{im.image:%Y-%m-%d}" + (" *" if im.partial else "") for _, im in imgs.iterrows()],
                       fontsize=7)
    ax.set_xlim(0.5, 90.5)
    ax.tick_params(labelbottom=False)
    for y in np.arange(len(imgs)) + 0.5:
        ax.axhline(y, color="white", lw=0.6)
    _title(ax, 0, f"Each image against the model since the previous image, {window[0]}-{window[1]}")

    ax_n = fig.add_subplot(gs[0, 1], sharey=ax)
    n_o = [np.nansum(d[d.obs_id == im.obs_id].observed) for _, im in imgs.iterrows()]
    n_m = [d[(d.obs_id == im.obs_id) & d.observed.notna()].model.sum() for _, im in imgs.iterrows()]
    y = np.arange(len(imgs))
    ax_n.barh(y - 0.2, n_o, height=0.38, color=C["ACCENT"], label="observed")
    ax_n.barh(y + 0.2, n_m, height=0.38, color=C["ADDED"], label="model")
    ax_n.tick_params(labelleft=False)
    ax_n.set_xlabel("domains")
    open_frame(ax_n)
    ax_n.legend(frameon=False, fontsize=7, loc="lower right")
    _title(ax_n, 1, "count")

    ax_f = fig.add_subplot(gs[1, 0], sharex=ax)
    a = d.dropna(subset=["observed"])
    by = a.groupby("gis").agg(o=("observed", "mean"), m=("model", "mean"))
    ax_f.plot(by.index, by.o, color=C["ACCENT"], lw=1.3, label="observed")
    ax_f.plot(by.index, by.m, color=C["ADDED"], lw=1.3, label="model (largest storm)")
    anyst = a.groupby("gis").model_any_storm.mean()
    ax_f.plot(anyst.index, anyst, color=C["ADDED"], lw=0.8, ls=(0, (2, 1.5)), label="model (any storm, upper bound)")
    ax_f.set_ylabel("share of images\nwith overwash")
    ax_f.set_xlabel(DOMAIN_AXIS_LABEL)
    ax_f.set_ylim(0, 1)
    town_bands(ax_f)
    open_frame(ax_f)
    ax_f.legend(frameon=False, fontsize=7, ncol=3, loc="upper left")
    _title(ax_f, 2, "How often each domain overwashes")
    fig.legend(handles=[Patch(fc=c_, ec="0.6", lw=0.4, label=l_) for c_, l_ in zip(CLASS_COLOURS, CLASS_LABELS)],
               loc="outside lower center", ncol=5, frameon=False, fontsize=7.5)

    ow.VS_MODEL.mkdir(parents=True, exist_ok=True)
    out = save(fig, ow.VS_MODEL / f"overwash_vs_model_{window[0]}_{window[1]}.png")
    plt.close(fig)
    s = scores(d)
    record_caption(out[0],
        f"Observed washover against modelled overwash, {window[0]}-{window[1]}, the managed run "
        f"({d.run.iloc[0]}). (a) One row per image (* = its window starts before the run, so the model covers "
        "only part of it); one column per GIS domain. Each image is compared with the model's overwash "
        "since the previous image; a model year's overwash is dated by its largest storm, with the 7-day grace of the observed storm table. Dark: both show "
        "overwash; purple: observed only; amber: model only; light grey: neither; white: domain not "
        "assessed in that image. (b) Domains with overwash per image. (c) The share of images with overwash "
        "per domain; the dashed line dates each model year's overwash by every storm in it instead of the "
        f"largest, an upper bound. Of {s['cells']} assessed image-domain cells: {s['both']} both, "
        f"{s['observed_only']} observed only, {s['model_only']} model only; the model catches "
        f"{s['hit_rate']:.0%} of observed overwash. Washover fades and images are months to years apart, so "
        "'model only' is not by itself a model error. GIS 1 is Cape Point, GIS 90 Pea Island.")
    return out


# Run: every window, arm and threshold, the tables and the figures
def main():
    apply_style()
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    frames = []
    for window in WINDOWS:
        for arm in ("managed", "natural"):
            for thr in THRESHOLDS:
                frames.append(compare(window, arm, obs, thr))
    df = pd.concat(frames, ignore_index=True)
    ow.VS_MODEL_TABLES.mkdir(parents=True, exist_ok=True)
    df.to_csv(ow.VS_MODEL_TABLES / "overwash_vs_model_cells.csv", index=False)
    summ = []
    for (w, arm, thr), g in df.groupby(["window", "arm", "threshold"]):
        for col in ("model", "model_any_storm"):
            summ.append(dict(window=w, arm=arm, threshold_m3_per_m=thr, dating=col, **scores(g, col)))
    summ = pd.DataFrame(summ)
    summ.to_csv(ow.VS_MODEL_TABLES / "overwash_vs_model_summary.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(summ.round(2).to_string(index=False))
    for window in WINDOWS:
        out = figure(df, window)
        print(" ->", out[0].relative_to(REPO))


if __name__ == "__main__":
    main()
