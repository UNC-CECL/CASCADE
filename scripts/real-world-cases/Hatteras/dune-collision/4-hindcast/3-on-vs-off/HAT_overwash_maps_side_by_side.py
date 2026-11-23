"""
Modelled and observed overwash on two island maps of the same size (2026-11-19): left, a run's net overwash per
GIS domain over 1996-2009; right, the number of aerial images showing washover per domain.

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_overwash_maps_side_by_side.py
    python .../HAT_overwash_maps_side_by_side.py grid_lowerfront_equal_cs200_vdefault 1996_2025

Arguments: the run (default Lower front row, equal rows, Cs 2e-4) and the observation window (default 1997_2009,
the run's own period; 1996_2025 is the full record of overwash_where_record_1996_2025.png). Both maps come from
the same geometry (overwash_map_periods.load_geometry / draw_island), so they share extent and scale. Writes
3-collision-design/0-all-designs/overwash_model_vs_observed_<run>_<window>.png.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-19
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_HAT = next(p for p in _HERE.parents if (p / "hatteras_paths.py").exists())
for p in (_HAT, _HERE.parent, _HAT / "input_prep" / "8-overwash-analysis" / "2-record",
          _HAT / "input_prep" / "8-overwash-analysis" / "1-observations"):
    sys.path.insert(0, str(p))
import collision_designs as CD  # noqa: E402
import HAT_collision_row_mode_sweep as R  # noqa: E402

H = R.H

# --- CONFIG ------------------------------------------------------------------
DEFAULT_RUN = "grid_lowerfront_equal_cs200_vdefault"
DEFAULT_WINDOW = "1997_2009"
MODEL_MAX = 1000.0                        # m3/m, top of the model colour scale (square-root)
C_NONE = "#ececec"
# -----------------------------------------------------------------------------


def main(run=DEFAULT_RUN, window=DEFAULT_WINDOW):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import BoundaryNorm, ListedColormap, PowerNorm
    from overwash_map_periods import draw_island, load_geometry, scalebar_and_north
    from overwash_data import load_observations, period_images, period_label
    from overwash_record_heatmap import count_cmap, per_domain
    from site_layer import hat_figure_style as S
    S.apply_style()

    # model: the run's own per-domain total (written by overwash_total.py with the run)
    tot = pd.read_csv(H.run_dir(run) / "tables" / "overwash_total.csv").set_index("gis").net_overwash_m3_per_m
    cmap_m = ListedColormap(plt.get_cmap("YlOrRd")(np.linspace(0.10, 0.80, 256)))   # yellow -> red, no brown
    norm_m = PowerNorm(0.5, vmin=0.0, vmax=MODEL_MAX, clip=True)
    fill_m = {int(g): (cmap_m(norm_m(v)) if v > 0 else C_NONE) for g, v in tot.items()}

    # observed: images with washover per domain, in the window
    obs, domains, matrix = load_observations()
    keep = period_images(obs, window)
    obs, matrix = obs.loc[keep].reset_index(drop=True), matrix[keep]
    rows = per_domain(obs, domains, matrix)
    vmax = max(r[1] for r in rows)
    cols = count_cmap(vmax)
    fill_o = {d: cols[h - 1] for d, h, _ in rows if h > 0}

    geom = load_geometry()
    dom, land, road, bounds = geom
    fig, (ax_m, ax_o) = plt.subplots(1, 2, figsize=(S.FIG_W_DOUBLE * 0.8, 8.4), layout="constrained")
    from HAT_collision_design_cs_sweep import cs_label
    label = next((d.name for d in CD.DESIGNS if d.run == run), run)
    for (key, c), v in H.DESIGN_SWEEP.items():
        if v == run:
            e, s = key.split("_")
            label = f"{CD.EROSION[e][0]}, {CD.START[s][0].lower()}, Cs {cs_label(c)}"
    draw_island(ax_m, dom, land, road, bounds, fill_m, "a", "Model, 1996–2009", reach_labels=False, pad_w=1500.0)
    draw_island(ax_o, dom, land, road, bounds, fill_o, "b", f"Observed, {period_label(window)}",
                reach_labels=False, pad_w=1500.0)
    # the scale bar east of the island, clear of Cape Point (as overwash_record_heatmap.py places it)
    x0, x1 = ax_m.get_xlim()
    ax_m.set_xlim(x0 + 0.44 * (x1 - x0), x1 + 0.44 * (x1 - x0))
    scalebar_and_north(ax_m)
    ax_m.set_xlim(x0, x1)
    cb = fig.colorbar(ScalarMappable(norm=norm_m, cmap=cmap_m), ax=ax_m, orientation="horizontal",
                      location="bottom", shrink=0.85, pad=0.01, aspect=22)
    cb.set_ticks([0, 50, 100, 250, 500, 1000])
    cb.set_label("m³ per m of shore (square-root scale); grey: none", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    bn = BoundaryNorm(np.arange(0.5, vmax + 1.5), vmax)
    cb2 = fig.colorbar(ScalarMappable(norm=bn, cmap=ListedColormap(cols)), ax=ax_o, orientation="horizontal",
                       location="bottom", shrink=0.85, pad=0.01, aspect=22, ticks=range(1, vmax + 1))
    cb2.set_label(f"images with washover (of {len(obs)}); bare land: none", fontsize=7.5)
    cb2.ax.tick_params(labelsize=7)
    fig.suptitle(label, fontsize=10, fontweight="bold")
    o = CD.ALL_DIR / f"overwash_model_vs_observed_{run}_{window}.png"
    S.save(fig, o, vector=False, close=True)
    S.record_caption(o, (
        f"Overwash along Hatteras Island on two maps at one extent and scale: (a) the run's net overwash per GIS "
        f"domain over 1996-2009 (run {run}; its tables/overwash_total.csv), yellow to red on a square-root scale, "
        f"grey where it never overwashed; (b) the number of the {len(obs)} aerial images of {period_label(window)} "
        "showing washover in each domain (Hatteras_Overwash_Data.xlsx), purples, bare land where none did."))
    print(f"wrote {o}")


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0] if a else DEFAULT_RUN, a[1] if len(a) > 1 else DEFAULT_WINDOW)
