"""
Version 2 of the crest picker view (2026-11-05), for a paper: one figure per domain showing only the final crests.

    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_crest_picks_paper.py [--gis 22]

Per transect, the 1996 profile and the dune crest used: your picked or accepted crest, else the automatic one (the
front dune). A rejected transect is drawn in grey and marked as excluded. No toes, dune window, automatic-vs-picked
marks or status: the picker views (crest_picker_views/) keep those.

Reads viewer_data_1996.json and crest_picks_1996.json; writes 5-crest-pick-figures/crest_picks_gis<NN>.png
(PDF under supporting/, captions in supporting/CAPTIONS.md). Without --gis, every domain.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-05
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import HAT_toe_common as DT  # noqa: E402

OUT_DIR = DT.OUT / "5-crest-pick-figures"
PANEL_H = 1.1             # in, per transect
X_LIM = (-30, 150)        # m landward of the shoreline; the viewer profiles end at 150 m


def final_crests(rs, picks):
    """Per transect: (record, crest (x, z) or None, excluded)."""
    out = []
    for r in rs:
        e = picks.get(str(r["id"]), {})
        if e.get("status") == "rejected" or not r.get("z"):
            out.append((r, None, True))
            continue
        crest = (e["crest_x"], e["crest_z"]) if e.get("crest_x") is not None else (tuple(r["fc"]) if r.get("fc") else None)
        out.append((r, crest, crest is None))
    return out


def figure(gis, D, picks, quiet=False, wide=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from site_layer import hat_figure_style as fs
    fs.apply_style()

    x0 = D["x0"]
    rs = [r for r in D["transects"] if r["gis"] == gis]
    if not rs:
        return
    rows = final_crests(rs, picks)
    xlim = (X_LIM[0], wide or X_LIM[1])
    if wide:                                     # the viewer profiles end at 150 m: read the lidar out to `wide`
        import pandas as pd
        tt = pd.read_csv(DT.TOES / f"transect_toes_{DT.YEAR}.csv").set_index("transect")
        (tr, _), (src, scale, offset) = DT.transects(), DT.survey()
        for r in rs:
            s_, xy, _, _ = DT.line_points(tr[tr.ObjectID == r["id"]].geometry.iloc[0])
            sh = tt.shoreline_s_m[r["id"]]
            keep = (s_ - sh >= xlim[0]) & (s_ - sh <= xlim[1])
            zz = DT.sample(src, scale, offset, xy[keep], tr.crs)
            r["z"], r["xs"] = [None if not np.isfinite(v) else float(v) for v in zz], s_[keep] - sh
    off = [r["id"] for r, c, ex in rows if c and not xlim[0] <= c[0] <= xlim[1]]
    zs = [v for r in rs for v in (r.get("z") or []) if v is not None]
    ylo, yhi = np.floor(min([-0.5] + zs)), np.ceil(max([4.0] + zs) + 0.5)
    w, _ = fs.figsize("double")
    fig, axes = plt.subplots(len(rows), 1, figsize=(w, PANEL_H * len(rows) + 0.75), sharex=True, sharey=True,
                             squeeze=False, layout="constrained")
    axes = axes[:, 0]
    for k, (ax, (r, crest, ex)) in enumerate(zip(axes, rows)):
        if r.get("z"):
            z = np.array([np.nan if v is None else v for v in r["z"]])
            xs = r["xs"] if "xs" in r else np.arange(x0, x0 + len(z))
            ax.fill_between(xs, z, ylo, color=DT.FILL, lw=0, zorder=1)
            ax.plot(xs, z, color=DT.DROPPED if ex else DT.PROFILE, lw=1.0, zorder=2)
        if wide and r.get("w"):
            ax.axvspan(*r["w"], color=DT.WINDOW, alpha=0.13, lw=0, zorder=0.5)
        if crest and not ex and xlim[0] <= crest[0] <= xlim[1]:
            ax.plot(*crest, "^", ms=7.5, color=DT.CREST, mec="white", mew=0.6, zorder=5)
        ax.axhline(D["mhw"], color=DT.WATER, lw=0.6, zorder=1.5)
        ax.set_ylim(ylo, yhi)
        ax.yaxis.set_major_locator(MaxNLocator(3, integer=True))
        ax.set_xlim(*xlim)
        fs.open_frame(ax)
        ax.text(0.008, 0.96, f"({chr(97 + k)}) Transect {r['id']}", transform=ax.transAxes, ha="left", va="top",
                fontsize=8, fontweight="bold")
        info = "excluded" if ex else (f"crest {crest[1]:.2f} m" + (f" at {crest[0]:.0f} m, beyond the view"
                                                                    if r["id"] in off else ""))
        ax.text(0.008, 0.70, info, transform=ax.transAxes, ha="left", va="top", fontsize=7.5,
                color=fs.INK_MUTED if ex else fs.INK, style="italic" if ex else "normal")
    axes[-1].text(xlim[0] + 1.5, D["mhw"] + 0.15, "MHW", fontsize=6.5, color=fs.INK_MUTED, ha="left", va="bottom")
    axes[-1].set_xlabel("Distance landward of the shoreline (m)")
    fig.supylabel("Elevation (m NAVD88)", fontsize=9)
    hd = [plt.Line2D([], [], color=DT.PROFILE, lw=1.0, label=f"{DT.YEAR} lidar profile"),
          plt.Line2D([], [], marker="^", ms=7.5, lw=0, color=DT.CREST, label="Dune crest")]
    if wide:
        from matplotlib.patches import Patch
        hd.append(Patch(color=DT.WINDOW, alpha=0.25, lw=0, label="Model's dune search window"))
    fig.suptitle(f"GIS {gis}", x=0.01, ha="left", fontsize=10, fontweight="bold")
    fig.legend(handles=hd, loc="outside upper center", ncol=len(hd), frameon=False, fontsize=7.5,
               handlelength=1.8, columnspacing=1.6)

    png = OUT_DIR / f"crest_picks_gis{gis:02d}{'_wide' if wide else ''}.png"
    fs.save(fig, png, close=True)
    n_ex = sum(ex for *_, ex in rows)
    fs.record_caption(png, (
        f"Dune crests on the {len(rows)} shore-normal transects of GIS domain {gis}, {DT.YEAR} lidar (elevations m "
        f"NAVD88; distance from the shoreline, the most seaward point at mean high water, thin line). Triangle: the "
        f"crest of the front dune on each profile."
        + (" Shaded: the model's dune search window, the crest is its highest point." if wide else "")
        + (f" {n_ex} transect{'s' if n_ex != 1 else ''} excluded (grey)." if n_ex else "")))
    if not quiet:
        print(f"wrote {png}" + (f"   crest beyond {X_LIM[1]} m on transects {off}" if off else ""))
    return off


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gis", type=int, help="one domain; every domain without it")
    ap.add_argument("--wide", type=float, help="with --gis: profiles out to this many m from the shoreline, read "
                    "from the lidar, with the model's dune search window (figure *_wide.png)")
    a = ap.parse_args()
    D = json.loads((DT.TOES / f"viewer_data_{DT.YEAR}.json").read_text(encoding="utf-8"))
    picks = DT.crest_picks()
    if a.gis is not None:
        return figure(a.gis, D, picks, wide=a.wide)
    beyond = {}
    for g in range(1, 91):
        off = figure(g, D, picks, quiet=True)
        if off:
            beyond[g] = off
    print(f"wrote the domain figures in {OUT_DIR}")
    if beyond:
        print(f"crests beyond {X_LIM[1]} m of the shoreline (past the viewer profiles): {beyond}")


if __name__ == "__main__":
    main()
