"""
How the beach slope is measured, on two 1996 transects: a steep beach north of Buxton and a gentle one (2026-11-23).

    python scripts/real-world-cases/Hatteras/input_prep/10-beach-slope/HAT_beach_slope_method_figure.py

For each example domain the transect whose slope is nearest that domain's 1996 median is drawn: the lidar sampled at
1 m (grey), the 5 m median-filtered profile the slope is read from (blue), MHW and the berm, the two crossings and
the slope between them. Reads foreshore_slope_transects.csv and the 1996 lidar (as HAT_beach_slope_mhw_berm.py).
Writes (D:, 10-beach-slope/1-foreshore-mhw-to-berm/): beach_slope_method.png (caption in CAPTIONS.md).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
"""
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
_HAT_DIR = next(p for p in _HERE.parents if (p / "hatteras_paths.py").exists())
sys.path.insert(0, str(_HAT_DIR))
sys.path.insert(0, str(_HERE.parent))

import numpy as np
import pandas as pd
import HAT_beach_slope_mhw_berm as M
from site_layer import hat_lidar_strips as L

SURVEY = 1996
EXAMPLES = ((10, "steep: GIS 10, north of Buxton"), (25, "gentle: GIS 25, south of Avon"))


def main():
    import geopandas as gpd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.ndimage import map_coordinates
    from site_layer import hat_figure_style as S
    S.apply_style()
    t_all = pd.read_csv(M.OUT / "foreshore_slope_transects.csv")
    t_all = t_all[t_all.survey == SURVEY]
    entries, tabs = {}, []
    for start, (_, ent, *_) in L.S.STARTS.items():
        entries.update(ent)
        tabs.append(pd.read_csv(L.S.tv.product_dir(start) / "survey-crests" / "survey_sources.csv"))
    src_tab = pd.concat(tabs).drop_duplicates("survey").set_index("survey")

    fig, axs = plt.subplots(1, 2, figsize=(S.FIG_W_DOUBLE, 3.2), sharey=True, layout="constrained")
    for ax, (gis, label) in zip(axs, EXAMPLES):
        tt = t_all[t_all.gis == gis]
        med = tt.slope.median()
        line_id = tt.iloc[(tt.slope - med).abs().argsort().iloc[0]].transect
        g = L.grid(L.Region("one", "one", range(gis, gis + 1), (gis,)))
        tf, _, crs = g
        z = M.load(SURVEY, g, src_tab, entries)
        tr = gpd.read_file(L.TRANSECTS).to_crs(crs)
        geom = tr[tr["Transects_100m.LineID"] == line_id].geometry.iloc[0]
        c = np.asarray(geom.coords, float)
        u = (c[-1] - c[0]) / np.hypot(*(c[-1] - c[0]))
        s = np.arange(0, geom.length, 1.0)
        xy = c[0] + s[:, None] * u
        col, rr = (xy[:, 0] - tf.c) / tf.a - 0.5, (xy[:, 1] - tf.f) / tf.e - 0.5
        raw = map_coordinates(z, [rr, col], order=1, mode="constant", cval=np.nan)
        p, i_m, i_b = M.foreshore_points(raw)
        x = s - s[i_m]                                   # m landward of the MHW crossing
        w = x[i_b]
        slope = (M.BERM_NAVD88_M - M.MHW_NAVD88_M) / w

        ax.plot(x, raw, color="0.6", lw=0.8, label="Lidar, 1 m samples")
        ax.plot(x, p, color=S.C_1997, lw=1.6, label="Filtered (5 m median)")
        ax.axhline(M.MHW_NAVD88_M, color="#2171B5", lw=0.8, ls=(0, (4, 2)))
        ax.axhline(M.BERM_NAVD88_M, color="0.35", lw=0.8, ls=(0, (1, 1)))
        ax.text(-28, M.MHW_NAVD88_M + 0.12, "MHW 0.36 m", fontsize=7, color="#2171B5")
        ax.text(-28, M.BERM_NAVD88_M + 0.12, "Berm 1.7 m", fontsize=7, color="0.35")
        ax.plot([0, w], [M.MHW_NAVD88_M, M.BERM_NAVD88_M], color="#D7301F", lw=2.2, zorder=4,
                label="Foreshore slope")
        ax.scatter([0, w], [M.MHW_NAVD88_M, M.BERM_NAVD88_M], color="#D7301F", s=24, zorder=5)
        ax.annotate("", xy=(w, -0.6), xytext=(0, -0.6),
                    arrowprops=dict(arrowstyle="<->", color="0.2", lw=0.8))
        ax.text(w / 2, -1.05, f"{w:.0f} m", ha="center", fontsize=7.5)
        ax.text(w + 3, 1.0, f"slope = 1.34 m / {w:.0f} m\n= {slope:.3f}", fontsize=7.5, color="#D7301F",
                va="center")
        ax.set_xlim(-30, 110)
        ax.set_ylim(-1.5, 8)
        ax.set_xlabel("Distance landward of the MHW crossing (m)")
        ax.set_title(f"({'ab'[EXAMPLES.index((gis, label))]}) {label}\ntransect {int(line_id)}, domain median "
                     f"{med:.3f}", loc="left", fontsize=8.5)
    axs[0].set_ylabel("Elevation (m NAVD88)")
    axs[0].legend(loc="upper left", frameon=False, fontsize=7)

    out = M.OUT / "beach_slope_method.png"
    S.save(fig, out, close=True)
    S.record_caption(out, (
        f"How the beach slope is measured, on the {SURVEY} lidar. Each 100 m transect is sampled every 1 m from the "
        "ocean end (grey) and median-filtered over 5 m to remove single spikes (blue). From the ocean, the first point "
        "above MHW (0.36 m NAVD88) with real data just seaward of it (so a gap in the lidar is not taken for the "
        "water's edge), then the first point above the berm (1.7 m NAVD88) landward of it (red dots). The foreshore "
        "slope is the 1.34 m rise between them over their distance (red line); slopes outside 0.01-0.3 are dropped. "
        "Each domain's slope is the median of its transects. Shown: the transect nearest its domain's median, (a) a "
        "steep beach north of Buxton, (b) a gentle beach south of Avon. HAT_beach_slope_mhw_berm.py."))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
