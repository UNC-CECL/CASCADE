"""
When and where the collision module misses the measured dune crest drops (2026-11-11). For every domain and every
pair of successive lidar surveys, the change in the domain-median crest in the lidar against the change in the model
over the same interval (the model state at 1 January after each survey), for the current setup: the toe at the berm,
the road fix, Cs 2e-6 (cal_toeberm_cs2); the module-off run is listed beside it.

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_missed_drops.py

A drop is a fall of 0.5 m or more in the domain median between two surveys. Pairs where either survey covers fewer
than 10 of a domain's 50 cells are left out (grey), as are domains 1-5, where the module does not run.

Writes, in the workspace's dune-collision/4-hindcast/1-module-on-off/ (D:):
- missed_crest_drops.png: (a) lidar change, (b) model change, (c) each pair classed: lidar drop missed (model under
  half of it), lidar drop caught (half or more), model drop the lidar does not show, neither
- supporting/missed_crest_drops.csv: the pairs, with the module-off and Cs 1e-6 changes too

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-11
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_collision_on_off_position as P  # noqa: E402
import HAT_collision_row_mode_sweep as R  # noqa: E402

H = R.H

# --- CONFIG ------------------------------------------------------------------
OUT = P.OUT / "1-module-on-off"
ON, OFF, ON_LOW = "cal_toeberm_cs2", "cal_step2_growth_noscrape", "cal_toeberm_cs1"
SURVEYS = (1996, 1997, 1998, 1999, 2000, 2001, 2004, 2005, 2009)
STORMS = {(1998, 1999): "Bonnie, Dennis", (2001, 2004): "Isabel"}   # the two the testing log reads (sections 4, 11)
NO_STEP = {(1999, 2000)}   # both surveys map to the 1 Jan 2000 model state: the model cannot change between them
GIS = range(6, 91)                                   # where the module runs
DROP_M, CAUGHT, MIN_CELLS = 0.5, 0.5, 10
VILLAGES = {"Buxton": (7, 8), "Avon": (21, 31), "Tri-Village": (68, 83)}
CLASSES = (("Lidar drop missed", "#D7301F"), ("Lidar drop caught", "#1A9850"),
           ("Model drop, not in lidar", "#2171B5"), ("No drop in either", "#F2F2F2"))
# -----------------------------------------------------------------------------


def pairs(v):
    from cascade_pipeline.plotting.dune_crest_comparison import TABLE, TABLES
    t = pd.read_csv(H.run_dir(v) / TABLES / TABLE)
    t = t[t.gis.isin(GIS)].set_index(["gis", "survey"])
    rows = []
    for g in GIS:
        for a, b in zip(SURVEYS[:-1], SURVEYS[1:]):
            if (g, a) not in t.index or (g, b) not in t.index:
                continue
            A, B = t.loc[(g, a)], t.loc[(g, b)]
            rows.append(dict(gis=g, a=a, b=b, cells=int(min(A.cells, B.cells)),
                             lidar=B.survey_median_m - A.survey_median_m, model=B.model_median_m - A.model_median_m))
    return pd.DataFrame(rows)


def classify(r):
    if (r.a, r.b) in NO_STEP:
        return 4
    if r.cells < MIN_CELLS:
        return np.nan
    if r.lidar <= -DROP_M:
        return 1 if r.model <= CAUGHT * r.lidar else 0
    if r.model <= -DROP_M:
        return 2
    return 3


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    from site_layer import hat_figure_style as S
    S.apply_style()
    (OUT / "supporting").mkdir(parents=True, exist_ok=True)
    d = pairs(ON)
    d = d.merge(pairs(OFF)[["gis", "a", "b", "model"]].rename(columns={"model": "model_off"}), on=["gis", "a", "b"])
    d = d.merge(pairs(ON_LOW)[["gis", "a", "b", "model"]].rename(columns={"model": "model_cs1"}), on=["gis", "a", "b"])
    d["cls"] = d.apply(classify, axis=1)
    d["class"] = d.cls.map(lambda c: "no model step" if c == 4 else CLASSES[int(c)][0] if np.isfinite(c)
                           else "too few cells")
    d.round(3).to_csv(OUT / "supporting" / "missed_crest_drops.csv", index=False)

    ivs = list(zip(SURVEYS[:-1], SURVEYS[1:]))
    gis = list(GIS)

    def grid(col, mask_cells=True):
        z = np.full((len(ivs), len(gis)), np.nan)
        for r in d.itertuples():
            if (r.a, r.b) in ivs and (not mask_cells or r.cells >= MIN_CELLS):
                z[ivs.index((r.a, r.b)), gis.index(r.gis)] = getattr(r, col)
        return z

    ext = (gis[0] - 0.5, gis[-1] + 0.5, len(ivs) - 0.5, -0.5)
    fig, axs = plt.subplots(3, 1, figsize=(S.FIG_W_DOUBLE, 7.6), sharex=True, layout="constrained")
    cmap = plt.get_cmap("RdBu").copy()
    cmap.set_bad("0.82")
    for ax, col, lab in ((axs[0], "lidar", "(a) Lidar"), (axs[1], "model", "(b) Model, Cs 2e-6")):
        im = ax.imshow(np.ma.masked_invalid(grid(col)), extent=ext, aspect="auto", cmap=cmap, vmin=-2, vmax=2,
                       interpolation="nearest")
        ax.set_title(lab, loc="left", fontsize=9, fontweight="bold", pad=14 if col == "lidar" else 6)
    cb = fig.colorbar(im, ax=axs[:2], shrink=0.8, aspect=25, pad=0.01, extend="both")
    cb.set_label("Crest change between surveys,\ndomain median (m)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ccmap = ListedColormap([c for _, c in CLASSES] + ["white"])
    ccmap.set_bad("0.82")
    axs[2].imshow(np.ma.masked_invalid(grid("cls")), extent=ext, aspect="auto", cmap=ccmap,
                  norm=BoundaryNorm(np.arange(-0.5, 5.5), 5), interpolation="nearest")
    for a, b in NO_STEP:
        k = ivs.index((a, b))
        for ax in axs[1:]:
            ax.axhspan(k - 0.5, k + 0.5, facecolor="white", edgecolor="0.6", hatch="////", lw=0, zorder=3)
            ax.text(np.mean(gis), k, "no model step between these surveys (both read the 1 Jan 2000 state)",
                    ha="center", va="center", fontsize=6.5, color="0.3", zorder=4,
                    bbox=dict(facecolor="white", edgecolor="none", pad=1))
    axs[2].set_title(f"(c) Drops of {DROP_M:g} m or more: caught = model at least half the lidar drop", loc="left",
                     fontsize=9, fontweight="bold")
    labels = [f"{a}–{str(b)[2:]}" + (f"  {STORMS[(a, b)]}" if (a, b) in STORMS else "") for a, b in ivs]
    for ax in axs:
        ax.set_yticks(range(len(ivs)), labels=labels, fontsize=6.5)
        for x in np.arange(gis[0] - 0.5, gis[-1] + 1, 1):
            ax.axvline(x, color="white", lw=0.25)
        for name, (lo, hi) in VILLAGES.items():
            ax.axvline(lo - 0.5, color="0.2", lw=0.6, ls=(0, (3, 2)))
            ax.axvline(hi + 0.5, color="0.2", lw=0.6, ls=(0, (3, 2)))
        ax.tick_params(labelsize=7)
    for name, (lo, hi) in VILLAGES.items():
        axs[0].text((lo + hi) / 2, 1.02, name, transform=axs[0].get_xaxis_transform(), ha="center", fontsize=7,
                    color="0.25")
    axs[-1].set_xticks(range(10, 91, 10))
    axs[-1].set_xlabel("GIS domain (south to north)")
    counts = d[d.cls < 4].cls.astype(int).value_counts()
    h = [Patch(color=c, label=f"{n} ({counts.get(k, 0)})") for k, (n, c) in enumerate(CLASSES)]
    h.append(Patch(color="0.82", label=f"Survey covers < {MIN_CELLS} cells"))
    fig.legend(handles=h, loc="outside lower center", ncol=5, fontsize=7, frameon=False)
    out = OUT / "missed_crest_drops.png"
    S.save(fig, out, bbox_inches="tight", close=True)
    S.record_caption(out, (
        "Change in the domain-median dune crest between successive lidar surveys (rows; Bonnie and Dennis "
        "fall in 1998-99, Isabel in 2001-04), GIS 6-90, (a) in the lidar and (b) in the model (toe at the berm, road fix, Cs 2e-6; model "
        "state at 1 January after each survey); red lower, blue higher, clipped at +-2 m. (c) Each domain and interval "
        f"classed: a lidar drop of {DROP_M:g} m or more that the model misses (falls by less than half of it) or "
        "catches, a model drop of that size the lidar does not show, or neither. Grey: either survey covers fewer "
        f"than {MIN_CELLS} of the domain's 50 cells. Hatched: 1999-2000, where both surveys read the same model state "
        "(1 January 2000), so the model cannot change; left out of the counts. Dashed: the villages (Buxton GIS 7-8, Avon 21-31, Tri-Village "
        "68-83), where grading, construction or buildings in the crest pick can also lower a measured crest."))
    print(f"wrote {out}")
    ok = d[d.cls < 4]
    print(ok.groupby("class").size())
    drop = ok[ok.lidar <= -DROP_M]
    print("lidar drops", len(drop), "caught", int((drop.cls == 1).sum()),
          "model off caught", int((drop.model_off <= CAUGHT * drop.lidar).sum()),
          "Cs 1e-6 caught", int((drop.model_cs1 <= CAUGHT * drop.lidar).sum()))


if __name__ == "__main__":
    main()
