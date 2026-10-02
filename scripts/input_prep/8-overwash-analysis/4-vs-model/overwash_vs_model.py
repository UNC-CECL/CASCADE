"""
Does the model overwash where and when the imagery shows washover?

    python scripts/input_prep/8-overwash-analysis/4-vs-model/overwash_vs_model.py

The observed record against Barrier3D's overwash in the 1996-2010 and 2010-2024
hindcast runs, each image against the storms since the previous one. Details: scripts/input_prep/8-overwash-analysis/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
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
    apply_style, C, INK_MUTED, figsize, save, record_caption,
    _title, open_frame, spines_for_image, _scalebar, _north_arrow,
)
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS  # noqa: E402

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
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "2-record"))
from overwash_map_periods import draw_island, load_geometry  # noqa: E402
GRACE = pd.Timedelta(days=CAPTURE_GRACE_DAYS)   # the observed storm table's own rule
HEADLINE_THR = 0.0

# cell classes: teal = the model overwashed (dark when the image confirms it), amber = it missed one
BOTH, OBS_ONLY, MOD_ONLY, NEITHER, UNASSESSED = 0, 1, 2, 3, 4
CLASS_COLOURS = ["#01665e", "#e6ab02", "#80cdc1", "#ececec", "white"]
CLASS_LABELS = ["Observed and modelled", "Observed, not modelled", "Modelled, not observed",
                "Neither", "Not assessed"]
C_OBS, C_MOD = "0.15", "#01665e"       # observed and modelled lines and bars


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


# A domain's class over the whole window: agreed at least once, else missed, else model only
def domain_class(g):
    a = g.dropna(subset=["observed"])
    if a.empty:
        return UNASSESSED
    o, m = a.observed.astype(int), a.model.astype(int)
    if ((o == 1) & (m == 1)).any():
        return BOTH
    if (o == 1).any():
        return OBS_ONLY
    if (m == 1).any():
        return MOD_ONLY
    return NEITHER


# Village spans as brackets beside a GIS-up axis, so they do not tint the cells
def town_rows(ax):
    tr = ax.get_yaxis_transform()
    for name, (lo, hi) in HATTERAS_ANNOTATIONS.town_spans.items():
        ax.plot([1.015, 1.015], [lo - 0.4, hi + 0.4], transform=tr, color=INK_MUTED, lw=0.8,
                clip_on=False, solid_capstyle="butt")
        ax.text(1.03, (lo + hi) / 2, name, transform=tr, rotation=90, ha="left", va="center",
                fontsize=6.5, color=INK_MUTED)


# Map + per-image counts + image-by-domain matrix, for one window (managed, headline threshold)
def figure(df, window, geom):
    d = df[(df.window == f"{window[0]}-{window[1]}") & (df.arm == "managed") & (df.threshold == HEADLINE_THR)]
    imgs = d[["obs_id", "image", "since", "partial"]].drop_duplicates().reset_index(drop=True)
    n = len(imgs)
    grid = np.full((90, n), UNASSESSED)
    for i, im in imgs.iterrows():
        sub = d[d.obs_id == im.obs_id].set_index("gis")
        for g in range(1, 91):
            grid[g - 1, i] = cell_class(sub.observed[g], sub.model[g])
    cmap = ListedColormap(CLASS_COLOURS[:4])
    norm = BoundaryNorm(np.arange(-0.5, 4.5), cmap.N)

    fig = plt.figure(figsize=figsize("double", height=8.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.9], height_ratios=[1, 6.2])

    # (a) the island, each domain by how it scored over the window
    dom, land, road, bounds = geom
    per = {g: domain_class(sub) for g, sub in d.groupby("gis")}
    ax_m = fig.add_subplot(gs[:, 0])
    draw_island(ax_m, dom, land, road, bounds, {g: CLASS_COLOURS[c] for g, c in per.items()}, "a",
                "Whole window", reach_labels=False, pad_w=2500, pad_e=2600, pad_s=5000)
    _scalebar(ax_m, 5000)
    _north_arrow(ax_m, x=0.85, y=0.06)
    cnt = {c: sum(v == c for v in per.values()) for c in range(5)}

    # (c) every image against the model since the previous image; GIS runs up the page like the map
    ax = fig.add_subplot(gs[1, 1])
    ax.add_patch(plt.Rectangle((-0.5, 0.5), n, 90, fc="white", ec="0.75", hatch="////", lw=0, zorder=0))
    ax.imshow(np.ma.masked_equal(grid, UNASSESSED), cmap=cmap, norm=norm, aspect="auto",
              interpolation="nearest", origin="lower", extent=(-0.5, n - 0.5, 0.5, 90.5), zorder=1)
    for x in np.arange(n - 1) + 0.5:
        ax.axvline(x, color="white", lw=1.2, zorder=2)
    town_rows(ax)
    ax.set_xticks(range(n))
    ax.set_xticklabels([f"{im.image:%b %Y}" + (" *" if im.partial else "") for _, im in imgs.iterrows()],
                       rotation=55, ha="right", rotation_mode="anchor", fontsize=7)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0.5, 90.5)
    ax.set_yticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90])
    ax.set_ylabel("GIS domain (south → north)")
    ax.set_xlabel("Image date")
    spines_for_image(ax)
    _title(ax, 2, "Each image, domain by domain")

    # (b) how many domains each image shows overwashed, and the model
    ax_n = fig.add_subplot(gs[0, 1], sharex=ax)
    n_o = [np.nansum(d[d.obs_id == im.obs_id].observed) for _, im in imgs.iterrows()]
    n_m = [d[(d.obs_id == im.obs_id) & d.observed.notna()].model.sum() for _, im in imgs.iterrows()]
    x = np.arange(n)
    ax_n.bar(x - 0.19, n_o, width=0.36, color=C_OBS, label="Observed")
    ax_n.bar(x + 0.19, n_m, width=0.36, color=C_MOD, label="Modelled")
    ax_n.set_ylim(0, 90)
    ax_n.set_yticks([0, 45, 90])
    ax_n.set_ylabel("domains")
    ax_n.tick_params(labelbottom=False)
    ax_n.grid(axis="y", color="0.9", lw=0.5)
    open_frame(ax_n)
    ax_n.legend(loc="upper right", frameon=False, fontsize=6.8, ncol=2)
    _title(ax_n, 1, "Domains with overwash")

    fig.legend(handles=[Patch(fc=c_, ec="0.6", lw=0.4, label=l_)
                        for c_, l_ in zip(CLASS_COLOURS[:4], CLASS_LABELS[:4])]
               + [Patch(fc="white", ec="0.6", hatch="////", lw=0.4, label=CLASS_LABELS[4]),
                  plt.Line2D([], [], color=C["ROAD"], lw=0.8, label="NC-12")],
               loc="outside lower center", ncol=6, frameon=False, fontsize=7, handlelength=1.4,
               columnspacing=1.2)

    ow.VS_MODEL.mkdir(parents=True, exist_ok=True)
    out = save(fig, ow.VS_MODEL / f"overwash_vs_model_{window[0]}_{window[1]}.png")
    plt.close(fig)
    s = scores(d)
    record_caption(out[0],
        f"Washover seen in the imagery against overwash in the model, {window[0]}-{window[1]}, the managed "
        f"run ({d.run.iloc[0]}, storm series {env.DEFAULT_STORM_VARIANT}). Teal means the model overwashed a "
        "domain: dark when the image also shows washover there, light when it does not. Amber is washover "
        "the model missed. Light grey: neither; hatched: domain not assessed in that image. "
        "(a) Each domain over the whole window: dark teal if model and image agreed in at least one image, "
        "else amber if washover was ever seen, else light teal if only the model overwashed "
        f"({cnt[BOTH]}, {cnt[OBS_ONLY]} and {cnt[MOD_ONLY]} domains; {cnt[NEITHER]} neither). A domain "
        "agreeing in one image can still be missed in others: (c) shows each image. NC-12 in black. (b) Domains with overwash in each image, observed and modelled. "
        "(c) One column per image, one row per domain; each image is scored against the model's overwash "
        "since the previous image (* = that interval starts before the run, so the model covers only part "
        "of it). A model year's overwash is dated by its largest storm, with the 7-day grace of the observed "
        "storm table. Brackets on the right mark Buxton, Avon and the Tri-Village. "
        f"Of {s['cells']} assessed image-domain cells: {s['both']} both, "
        f"{s['observed_only']} observed only, {s['model_only']} modelled only; the model catches "
        f"{s['hit_rate']:.0%} of observed washover. Washover fades and images are months to years apart, so "
        "'modelled, not observed' is not by itself a model error. GIS 1 is Cape Point, GIS 90 Pea Island.")
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
    geom = load_geometry()
    for window in WINDOWS:
        out = figure(df, window, geom)
        print(" ->", out[0].relative_to(REPO))


if __name__ == "__main__":
    main()
