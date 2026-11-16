"""
pybeach (Beuzen 2019, "pybeach: A Python package for extracting the location of dune toes on beach profile
transects", Journal of Open Source Software 4(44), 1890; the machine-learning toe from Beuzen et al. 2019, "Ensemble
models from machine learning: an example of wave runup and coastal dune erosion", Natural Hazards and Earth System
Sciences 19, 2295-2309) on the 1996 transects. Its four dune toe methods:
- perpendicular distance (pd): the point farthest below the line from the dune crest to the shoreline (our toe rule);
- maximum curvature (mc): the largest second derivative of the smoothed profile between the crest and the shoreline
  (Stockdon et al. 2007);
- relative relief (rr): the most seaward point between the crest and the shoreline where the profile's relief
  relative to its neighbourhood (21-point window) is below 0.2, else the lowest relief there (Wernette et al. 2016),
  with mean high water as the water level;
- machine learning (ml): a random forest trained on profiles with surveyed toes gives each point a probability of
  being the toe, from the profile's slope in a window around it; the toe is the most probable point seaward of the
  crest. Three pre-trained classifiers ship with pybeach: barrier island (used), mixed and wave-embayed.
pd, mc and rr are searched between our shoreline and crest (HAT_toe_method_line), ml seaward of our crest. pybeach puts
the profile on a 0.5 m grid with a 5-point moving average. Each toe is then raised to the berm as ours is (the raw
toe is kept as *_raw). Per domain, the median of the kept transects when at least 3.

Two steps, in two environments:
1. The toes, in pybeach's own environment (its pre-trained classifiers need scikit-learn 0.24, so Python 3.9):
       conda create -p <env> python=3.9 scikit-learn=0.24.2 "pandas<2" "joblib<2" scipy "numpy<1.24"
       <env>/python -m pip install --no-deps pybeach==0.1.3
       (Windows: put <env>, <env>/Library/bin and <env>/Scripts first on PATH, else numpy's DLLs fail to load)
       <env>/python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_method_pybeach.py
2. The method figure, in the project environment:
       python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_method_pybeach.py --figure

Reads 1-toes/viewer_data_1996.json (HAT_toe_pickers.py --viewer-data): each transect's profile from 40 m seaward to
150 m landward of its shoreline (1 m, m NAVD88) and its crest. Writes, in
data/hatteras_init/9-dune-toes/1-transects-1996/3-methods/4-pybeach/ (D:): pybeach_transects_1996.csv,
pybeach_domains_1996.csv, pybeach_examples_1996.npz (the working of each method on the example transects) and
pybeach_method_1996.png (the method figure).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-05
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# --- CONFIG ------------------------------------------------------------------
# As HAT_toe_common (not imported here: pybeach's environment has none of the project's packages)
OUT = Path(r"D:\research\UNC\hatteras_ch1\model_runs\data\hatteras_init\9-dune-toes\1-transects-1996")
VIEWER = OUT / "1-toes" / "viewer_data_1996.json"
DEST = OUT / "3-methods" / "4-pybeach"
BERM_M = 1.7                      # as HAT_toe_common.BERM_M
MIN_KEPT = 3
EXAMPLES = (119, 326)             # as HAT_toe_common.EXAMPLES
BEHIND_CREST_M = 5.0              # the profile handed to pybeach ends this far landward of the crest
CLASSIFIERS = {"ml": "barrier_island_clf", "ml_mixed": "mixed_clf", "ml_embayed": "wave_embayed_clf"}
RR_WINDOW, RR_THRESHOLD = 21, 0.2  # pybeach's relative relief defaults
METHODS = {"pd": "perpendicular distance", "mc": "maximum curvature", "rr": "relative relief",
           "ml": "machine learning (barrier island)", "ml_mixed": "machine learning (mixed)",
           "ml_embayed": "machine learning (wave-embayed)"}
# -----------------------------------------------------------------------------


def toes():
    warnings.filterwarnings("ignore")
    from pybeach.beach import Profile
    from pybeach.support import data_support as ds
    D = json.loads(VIEWER.read_text(encoding="utf-8"))
    mhw, x0 = D["mhw"], D["x0"]
    rows, ex = [], {}
    for t in D["transects"]:
        row = dict(transect=t["id"], gis=t["gis"], ours=t["ft"][1] if t["kept"] and t.get("ft") else np.nan)
        if not (t["kept"] and t.get("z") and t.get("fc")):
            rows.append(row)
            continue
        x = np.arange(x0, x0 + len(t["z"]), 1.0)                 # m landward of the shoreline
        z = np.array([np.nan if v is None else v for v in t["z"]])
        ok = np.isfinite(z) & (x <= t["fc"][0] + BEHIND_CREST_M)  # pybeach wants the profile to fall seaward:
        x, z = x[ok], z[ok]                                      # beach and foredune only, not the land behind
        g = np.arange(x.min(), x.max() + 0.25, 0.5)              # pybeach's 0.5 m grid, built here so indices match
        zg = np.interp(g, x, z)
        xl, zl = (g.max() - g)[::-1], zg[::-1]                   # land first (left), sea last
        dist = lambda k: float(g.max() - xl[k])                  # index -> m landward of the shoreline
        idx = lambda xm: int(np.argmin(np.abs(xl - (g.max() - xm))))
        crest, shore = idx(t["fc"][0]), idx(0.0)
        try:
            p = Profile(xl, zl)
        except Warning:                                          # pybeach's own orientation check refuses it
            row["refused"] = True
            rows.append(row)
            continue
        floor = lambda k: next((j for j in range(k, crest - 1, -1) if zl[j] >= BERM_M), crest)  # landward = left
        probs = {}
        with contextlib.redirect_stdout(io.StringIO()):          # relative relief prints a line per fallback
            picks = [("pd", p.predict_dunetoe_pd(dune_crest=crest, shoreline=shore)[0]),
                     ("mc", p.predict_dunetoe_mc(dune_crest=crest, shoreline=shore)[0]),
                     ("rr", p.predict_dunetoe_rr(water_level=mhw, dune_crest=crest, shoreline=shore)[0])]
            for col, clf in CLASSIFIERS.items():
                k, pr = p.predict_dunetoe_ml(clf, dune_crest=int(crest))
                picks.append((col, np.atleast_1d(k)[0]))
                probs[col] = pr[0]
        for name, k in picks:
            k = int(k)
            row[f"{name}_raw"] = float(zl[k])
            row[f"{name}_x"] = dist(k)
            row[name] = float(zl[k]) if zl[k] >= BERM_M else float(zl[floor(k)])
        rows.append(row)
        if t["id"] in EXAMPLES:                                  # the working of each method, for the figure
            zs = p.z_interp[0]
            d2 = np.full(len(zs), np.nan)
            d2[crest:shore + 1] = np.gradient(np.gradient(zs[crest:shore + 1], 0.5), 0.5)
            ex[t["id"]] = dict(x=np.array([dist(k) for k in range(len(xl))]), z=zl, zs=zs, d2=d2,
                               rr=ds.relative_relief(zs, np.atleast_1d(RR_WINDOW), mhw),
                               crest=crest, shore=shore, **{f"p_{c}": v for c, v in probs.items()},
                               **{f"k_{n}": int(k) for n, k in picks})
    c = pd.DataFrame(rows)
    print(f"pybeach refused {int(c.get('refused', pd.Series(dtype=bool)).fillna(False).sum())} profiles "
          f"(its land-left orientation check); {int(c.pd.notna().sum())} run")
    DEST.mkdir(parents=True, exist_ok=True)
    c.round(3).to_csv(DEST / "pybeach_transects_1996.csv", index=False)
    med = lambda col: c.groupby("gis")[col].median().where(c.groupby("gis")[col].count() >= MIN_KEPT)
    d = pd.DataFrame({k: med(k) for k in ["ours"] + list(METHODS)})
    d.index.name = "gis"
    d.round(3).to_csv(DEST / "pybeach_domains_1996.csv")
    np.savez(DEST / "pybeach_examples_1996.npz",
             **{f"{tid}_{k}": np.asarray(v) for tid, e in ex.items() for k, v in e.items()})
    for k, name in METHODS.items():
        dt = (c[k] - c.ours).dropna()
        dd = (d[k] - d.ours).dropna()
        print(f"pybeach {name:34s} transects within 0.25 m of ours {100 * (dt.abs() <= 0.25).mean():3.0f}%, "
              f"within 0.5 m {100 * (dt.abs() <= 0.5).mean():3.0f}% (n {len(dt)}); domain medians: mean "
              f"{dd.mean():+.2f} m, RMSE {np.sqrt((dd ** 2).mean()):.2f} m, > 0.5 m apart in "
              f"{int((dd.abs() > 0.5).sum())} of {len(dd)}")
    print(f"wrote {DEST}")


def figure():
    """The four methods on the two example transects (columns): (row 1) the profile with each toe and the
    perpendicular distance construction; (2) the curvature; (3) the relative relief; (4) the toe probability."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from HAT_toe_common import CREST, FILL, MHW, PROFILE, WATER   # first: it puts site_layer on the path
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    E = np.load(DEST / "pybeach_examples_1996.npz")
    col = {"pd": "#6E6E6E", "mc": "#C9A227", "rr": "#2E8B57", "ml": "#C2185B"}
    fig, axes = plt.subplots(4, 2, figsize=(fs.figsize("double")[0], 8.2), layout="constrained", sharex="col",
                             gridspec_kw=dict(height_ratios=(1.6, 1, 1, 1)))
    for j, tid in enumerate(EXAMPLES):
        g = lambda k: E[f"{tid}_{k}"]
        x, z = g("x"), g("z")
        cr, sh = int(g("crest")), int(g("shore"))
        k = {m: int(g(f"k_{m}")) for m in ("pd", "mc", "rr", "ml")}
        ax = axes[0, j]
        ax.fill_between(x, z, -1, color=FILL, lw=0)
        ax.plot(x, z, color=PROFILE, lw=1.0)
        ax.axhline(MHW, color=WATER, lw=0.6)
        ax.axhline(BERM_M, color=CREST, lw=0.6, ls=(0, (2, 2)))
        ax.plot(x[cr], z[cr], "^", ms=7, color=CREST, zorder=6, label="Our crest")
        ax.plot(x[sh], z[sh], "o", ms=5, color=WATER, zorder=6, label="Shoreline")
        ax.plot([x[sh], x[cr]], [z[sh], z[cr]], color=col["pd"], lw=0.8, ls=(0, (4, 2)))
        P0, P1, Q = np.array([x[sh], z[sh]]), np.array([x[cr], z[cr]]), np.array([x[k["pd"]], z[k["pd"]]])
        u = (P1 - P0) / np.hypot(*(P1 - P0))
        F = P0 + u * np.dot(Q - P0, u)
        ax.plot([Q[0], F[0]], [Q[1], F[1]], color=col["pd"], lw=0.8)
        for m, mk, lab in (("pd", "v", "Perpendicular distance"), ("mc", "^", "Maximum curvature"),
                           ("rr", "P", "Relative relief"), ("ml", "*", "Machine learning")):
            ax.plot(x[k[m]], z[k[m]], mk, ms=8 if m == "ml" else 6.5, color=col[m], mec="white", mew=0.5,
                    zorder=7, label=lab)
        ax.set_ylim(-1, np.nanmax(z) + 0.8)
        ax.text(0.02, 0.97, f"Transect {tid}", transform=ax.transAxes, fontsize=8, va="top")
        # the working of the other three
        a2 = axes[1, j]
        a2.plot(x, g("d2"), color=col["mc"], lw=1.0)
        a2.axvline(x[k["mc"]], color=col["mc"], lw=0.7, ls=(0, (3, 2)))
        a2.axhline(0, color="0.7", lw=0.5)
        a3 = axes[2, j]
        a3.plot(x, g("rr"), color=col["rr"], lw=1.0)
        a3.axhline(RR_THRESHOLD, color=col["rr"], lw=0.6, ls=(0, (2, 2)))
        a3.axvline(x[k["rr"]], color=col["rr"], lw=0.7, ls=(0, (3, 2)))
        a3.set_ylim(-0.05, 1.05)
        a4 = axes[3, j]
        for c_, ls, lab in (("ml", "-", "Barrier island (used)"), ("ml_mixed", (0, (3, 2)), "Mixed"),
                            ("ml_embayed", (0, (1, 1.5)), "Wave-embayed")):
            a4.plot(x, g(f"p_{c_}"), color=col["ml"], lw=1.0, ls=ls, label=lab)
        a4.axvline(x[k["ml"]], color=col["ml"], lw=0.7, ls=(0, (3, 2)))
        a4.set_xlabel("Distance landward of the shoreline (m)")
        for a in axes[:, j]:
            a.set_xlim(x.min(), x.max())
            fs.open_frame(a)
            a.axvspan(x[cr], x.max(), color="0.92", lw=0, zorder=0)   # landward of the crest: not searched
    axes[0, 0].set_ylabel("Elevation (m NAVD88)")
    axes[1, 0].set_ylabel("d²z/dx² (1/m)")
    axes[2, 0].set_ylabel("Relative relief")
    axes[3, 0].set_ylabel("Toe probability")
    for r, txt in enumerate(("Perpendicular distance: farthest below the shoreline-to-crest line",
                             "Maximum curvature: largest d²z/dx² between shoreline and crest",
                             f"Relative relief: most seaward point below {RR_THRESHOLD:g}",
                             "Machine learning: most probable point seaward of the crest")):
        axes[r, 0].text(0, 1.06, txt, transform=axes[r, 0].transAxes, fontsize=8, fontweight="bold")
    for r, a in enumerate(axes.flat):
        a.text(-0.02, 1.02, "abcdefgh"[r], transform=a.transAxes, fontsize=10, fontweight="bold", ha="right",
               va="bottom")
    h0, l0 = axes[0, 0].get_legend_handles_labels()
    h1, l1 = axes[3, 0].get_legend_handles_labels()
    fig.legend(h0 + h1, l0 + l1, loc="outside lower center", ncol=5, fontsize=7, frameon=False)
    png = DEST / "pybeach_method_1996.png"
    fs.save(fig, png, close=True)
    fs.record_caption(png, (
        f"pybeach's four dune toe methods (Beuzen 2019) on two 1996 lidar transects, {EXAMPLES[0]} (left, a single "
        f"foredune) and {EXAMPLES[1]} (right, a dune field), on pybeach's 0.5 m grid with its 5-point moving "
        f"average, given our shoreline and crest; grey, landward of the crest, not searched. (a, b) The profile with "
        f"each method's toe, before the berm floor (dotted, {BERM_M:g} m NAVD88), and the perpendicular-distance "
        f"construction: the toe is the point farthest below the shoreline-to-crest line. (c, d) Maximum curvature "
        f"(Stockdon et al. 2007): the second derivative of the profile; the toe is its largest value. (e, f) "
        f"Relative relief (Wernette et al. 2016; {RR_WINDOW}-point window): the toe is the most seaward point below "
        f"{RR_THRESHOLD:g} (dashed). (g, h) Machine learning (Beuzen et al. 2019): each point's probability of being the toe "
        f"from a random forest trained on surveyed profiles, for pybeach's three pre-trained classifiers; the toe "
        f"is the barrier island classifier's most probable point seaward of the crest."))
    print(f"wrote {png}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure", action="store_true", help="the method figure (project environment)")
    figure() if ap.parse_args().figure else toes()
