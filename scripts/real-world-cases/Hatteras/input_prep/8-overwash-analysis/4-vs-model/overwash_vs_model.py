"""
Does the model overwash where and when the imagery shows washover?

    python scripts/input_prep/8-overwash-analysis/4-vs-model/overwash_vs_model.py
    python scripts/input_prep/8-overwash-analysis/4-vs-model/overwash_vs_model.py --arm natural --window 1996
    python scripts/input_prep/8-overwash-analysis/4-vs-model/overwash_vs_model.py --run-dir output/raw_runs/.../HAT_...

The tables always hold both runs and both windows; the figures default to the managed
run in both windows (--arm, --window pick others; the natural figure gets a _natural suffix).
--run-dir draws one run's figure into its own figures/ folder instead (run_figure), as the
runner does at the end of every run; it can take several run folders.

The observed record against Barrier3D's overwash in the calendar-year calibration
(storms 1997-2009) and test (storms 2010-2025) runs, each image against the storms
since the previous one. Details: scripts/input_prep/8-overwash-analysis/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-22
"""

from __future__ import annotations

import argparse
import json

import sys
from pathlib import Path

import matplotlib

if __name__ == "__main__":
    matplotlib.use("Agg")   # the runner imports run_figure and keeps its own backend
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import ListedColormap, BoundaryNorm  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

# The scripts are in Hatteras/ (the folder holding hatteras_paths.py); data and output in the workspace
HATTERAS = next(p for p in Path(__file__).resolve().parents if (p / "hatteras_paths.py").exists())
sys.path.insert(0, str(HATTERAS))
from hatteras_paths import WORKSPACE as REPO  # noqa: E402

from site_layer import hat_overwash as ow  # noqa: E402
from site_layer import hat_env_forcings as env  # noqa: E402
from cascade_pipeline.run_layout import write_path  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, C, INK_MUTED, figsize, save, record_caption,
    _title, open_frame, spines_for_image, _scalebar, _north_arrow,
)
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS  # noqa: E402
from site_layer.hatteras_site_config import first_model_year, last_model_year  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
MATRIX = REPO / "output" / "raw_runs" / "matrix"
# Keyed by the matrix folder (DEM year, end year); the storms run first_model_year..last_model_year
WINDOWS = {
    (1996, 2009): {"managed": "HAT_1996_2009_edgeBE_offsetmetres_road_bdm_groin",
                   "natural": "HAT_1996_2009_edgeBE_offsetmetres_noroad_nobdm_nogroin"},
    (2009, 2025): {"managed": "HAT_2009_2025_edgeBE_offsetmetres_road_bdm_nourish_groin",
                   "natural": "HAT_2009_2025_edgeBE_offsetmetres_noroad_nobdm_nogroin"},
}
THRESHOLDS = (0.0, 1.0, 5.0)        # m3/m
# How each run is named in the captions
ARM_TEXT = {"managed": "managed run (groin on", "natural": "natural run (no road, no management, no groin"}
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


# The calendar years a window's run steps through (its storm file's window)
def model_years(window):
    return first_model_year(window[0]), last_model_year(window[0])


# How a window is written in captions and table rows: the storm years it runs
def label(window):
    a, b = model_years(window)
    return f"{a}-{b}"


# A cascade's overwash per model year and domain (QowTS, m3/m): (nt, pads), row t = model year t
def qow(cascade):
    return np.array([np.asarray(b.QowTS) for b in cascade.barrier3d]).T


# A matrix run's overwash and its name
def load_qow(window, arm):
    name = WINDOWS[window][arm]
    path = MATRIX / f"{window[0]}_{window[1]}" / "edgeBE" / name / f"{name}.npz"
    return qow(np.load(path, allow_pickle=True)["cascade"][0]), name


# Per model year: the end date of its largest storm, and of every storm
def storm_dates(window, f=None):
    f = f or env.storm_summary_file(*model_years(window))   # default: the series the matrix runs on
    s = pd.read_csv(f, parse_dates=["StartTime", "EndTime"])
    s["dated"] = s.EndTime - GRACE
    largest = s.loc[s.groupby("time").Rhigh.idxmax()].set_index("time").dated
    every = s.groupby("time").dated.apply(list)
    return largest, every, f


# Images whose storms the run covers, each with its comparison window
def images_in(window, obs):
    first, last = model_years(window)
    start = pd.Timestamp(f"{first}-01-01")
    end = pd.Timestamp(f"{last + 1}-01-01")         # the last modelled storm year ends here
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
    return compare_q(q, name, window, arm, obs, thr, storm_dates(window))


# The same for any run's overwash q, dated by that run's own storms
def compare_q(q, name, window, arm, obs, thr, dates):
    largest, every, _ = dates
    imgs = images_in(window, obs)
    nbuf = (q.shape[1] - 90) // 2       # from the run itself: runs before 2026-10-21 have 15 a side, not 30
    rows = []
    for _, im in imgs.iterrows():
        o = obs[obs.Obs_ID == im.Obs_ID].set_index("domain").overwash
        for gis in range(1, 91):
            p = nbuf + gis - DOM.first_gis_id
            yrs = [t for t in range(1, q.shape[0]) if q[t, p] > thr]
            m_big = any(im["from"] < largest[t] <= im.date for t in yrs if t in largest)
            m_any = any(im["from"] < d <= im.date for t in yrs if t in every for d in every[t])
            ob = o.get(gis, np.nan)
            rows.append(dict(window=label(window), arm=arm, run=name, threshold=thr,
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


# One matrix arm's figure in 4-vs-model (headline threshold)
def figure(df, window, geom, arm="managed"):
    d = df[(df.window == label(window)) & (df.arm == arm) & (df.threshold == HEADLINE_THR)]
    suffix = "" if arm == "managed" else f"_{arm}"
    out = ow.VS_MODEL / f"overwash_vs_model_{label(window).replace('-', '_')}{suffix}.png"
    return draw(d, window, geom, out,
                f"the {ARM_TEXT[arm]}, {d.run.iloc[0]}, storm series {env.DEFAULT_STORM_VARIANT})")


# Map + per-image counts + image-by-domain matrix, for one window and run; `run_desc` names the run
def draw(d, window, geom, out_path, run_desc):
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

    out = save(fig, out_path)
    plt.close(fig)
    s = scores(d)
    record_caption(out[0],
        f"Washover seen in the imagery against overwash in the model, storm years {label(window)}, "
        f"{run_desc}. Teal means the model overwashed a "
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


# One run's figure in its own figures/ folder, scored against the storms it ran; None when no image falls in its window
def run_figure(q, run_dir, run_name, start_year, storm_file):
    window = (start_year, None)
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    if images_in(window, obs).empty:
        print(f"  overwash vs observed  no images in storm years {label(window)}, skipped")
        return None
    storm_file = Path(storm_file)
    dates = storm_dates(window, storm_file.with_name(storm_file.stem + "_summary.csv"))
    d = compare_q(q, run_name, window, "run", obs, HEADLINE_THR, dates)
    variant = storm_file.stem.split("_storms_")[-1]
    with plt.rc_context():          # the house style for this figure only, not the caller's later ones
        apply_style()
        out = draw(d, window, load_geometry(), write_path(Path(run_dir), "figure_overwash", run_name),
                   f"run {run_name} (storm series {variant})")
    s = scores(d)
    fa = s["model_only"] / max(s["model_only"] + s["neither"], 1)
    print(f"  overwash vs observed  hit {s['hit_rate']:.2f}  false alarm {fa:.2f}  "
          f"skill {s['hit_rate'] - fa:.2f}  -> {out[0].name}")
    return out


# A finished run's figure from its folder: the archive, plus the period and storm file its metadata records
def run_figure_from_dir(run_dir):
    run_dir = Path(run_dir)
    name = run_dir.name
    meta = json.loads((run_dir / f"{name}_run_metadata.json").read_text(encoding="utf-8"))
    c = np.load(run_dir / f"{name}.npz", allow_pickle=True)["cascade"][0]
    storm_file = env.INIT_ROOT / Path(meta["identity"]["storm_file"].replace("\\", "/"))
    return run_figure(qow(c), run_dir, name, meta["period"]["start_year"], storm_file)


# Run: every window, arm and threshold, the tables and the figures
def main():
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--arm", choices=("managed", "natural"), default="managed", help="run to draw")
    ap.add_argument("--window", type=int, choices=[w[0] for w in WINDOWS],
                    help="DEM year of the window to draw (default: both)")
    ap.add_argument("--run-dir", nargs="+", type=Path,
                    help="run folder(s): draw each run's own figure into its figures/ and stop")
    args = ap.parse_args()
    if args.run_dir:
        for rd in args.run_dir:
            print(rd.name)
            run_figure_from_dir(rd)
        return
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
        if args.window and window[0] != args.window:
            continue
        out = figure(df, window, geom, args.arm)
        print(" ->", out[0].relative_to(REPO))


if __name__ == "__main__":
    main()
