"""
How much overwash one run made over its whole model period, where and when.

    python scripts/input_prep/8-overwash-analysis/4-vs-model/overwash_total.py --run-dir output/raw_runs/.../HAT_...

Draws figures/overwash_total.png and tables/overwash_total.csv in each run folder, as the
runner does at the end of every run (run_figure, from run_figures.py). The overwash is the
run's QowTS: net of what the beach-dune manager cleared off the road and filtered, in a run
that manages. Details: scripts/input_prep/8-overwash-analysis/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-27
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
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.colors import ListedColormap, PowerNorm  # noqa: E402

import overwash_vs_model as ovm  # noqa: E402  (beside this file; puts Hatteras/ on the path)
from site_layer import hat_overwash as ow  # noqa: E402
from cascade_pipeline.run_layout import write_path  # noqa: E402
from cascade_pipeline.saved_runs import load_saved_cascade  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, figsize, save, record_caption, _title, open_frame, _scalebar, _north_arrow,
)
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM  # noqa: E402
from overwash_map_periods import draw_island, load_geometry  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
# One scale for every run of a window (keyed by DEM year), so its runs compare by eye: (m3/m total
# per domain, map colour and (c) axis; 10^3 m3 island total per model year, (b) axis). The largest
# of the 2026-10-26 runs: 1996 910 and 3365, 2009 2165 and 4797 (no road, no management)
SCALES = {1996: (1000.0, 3500.0), 2009: (2500.0, 5000.0)}
# -----------------------------------------------------------------------------

# sand on a sqrt scale, as the volume bars of overwash_vs_observed; grey where the run never overwashed
CMAP = ListedColormap(plt.get_cmap("YlOrBr")(np.linspace(0.12, 1.0, 256)))
C_NONE = "#ececec"
C_YEAR = "0.1"                      # the island total per year, as overwash_vs_observed's volume panel
C_SEEN = "0.15"                     # washover seen in the imagery


# Per GIS domain (1-90): the run's total; per model year: the island total in 10^3 m3
def totals(q):
    nbuf = (q.shape[1] - 90) // 2               # from the run itself: 15 a side before 2026-10-21, 30 since
    real = q[1:, nbuf + 1 - DOM.first_gis_id: nbuf + 91 - DOM.first_gis_id]   # row 0 is the initial state
    return real.sum(axis=0), real.sum(axis=1) * DOM.domain_spacing_m / 1e3


# A window's scales; another window gets its own run's, rounded up
def scales(start_year, by_dom, by_year):
    if start_year in SCALES:
        return SCALES[start_year]
    up = lambda v, step: max(step, np.ceil(v / step) * step)
    return up(by_dom.max(), 500.0), up(by_year.max(), 500.0)


# GIS domains with washover in any image the run's storms cover
def seen_domains(start_year):
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    ids = ovm.images_in((start_year, None), obs).Obs_ID
    return set(obs[obs.Obs_ID.isin(ids) & (obs.overwash == 1)].domain.astype(int))


# Map + island total per year + total per domain, for one run
def draw(by_dom, by_year, years, seen, geom, out_path, run_desc, period, scale):
    dmax, ymax = scale
    norm = PowerNorm(0.5, vmin=0.0, vmax=dmax, clip=True)
    fig = plt.figure(figsize=figsize("double", height=8.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.9], height_ratios=[1, 6.2])
    gis = np.arange(1, 91)
    fill = [CMAP(norm(v)) if v > 0 else C_NONE for v in by_dom]

    # (a) the island, each domain by its total
    dom, land, road, bounds = geom
    ax_m = fig.add_subplot(gs[:, 0])
    draw_island(ax_m, dom, land, road, bounds, dict(zip(gis, fill)), "a",
                "Total over the run", reach_labels=False, pad_w=2500, pad_e=2600, pad_s=5000)
    _scalebar(ax_m, 5000)
    _north_arrow(ax_m, x=0.85, y=0.06)

    # (c) the same totals as bars, GIS up the page like the map
    ax = fig.add_subplot(gs[1, 1])
    ax.barh(gis, np.clip(by_dom, 0, None), height=0.8, color=fill, ec="0.6", lw=0.2)
    s = sorted(seen)
    ax.plot(np.full(len(s), dmax * 0.985), s, ls="none", marker="s", ms=2.6, color=C_SEEN,
            label="Washover seen in an image of the run's period")
    ovm.town_rows(ax)
    ax.set_xlim(0, dmax)
    ax.set_ylim(0.5, 90.5)
    ax.set_yticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90])
    ax.set_ylabel("GIS domain (south → north)")
    ax.set_xlabel("Net overwash over the run (m³/m)")
    ax.grid(axis="x", color="0.9", lw=0.5)
    open_frame(ax)
    ax.legend(loc="upper right", bbox_to_anchor=(0.95, 1.0), frameon=False, fontsize=6.8,
              handletextpad=0.3)
    _title(ax, 2, "Each domain")
    sm = ScalarMappable(norm=norm, cmap=CMAP)
    cb = fig.colorbar(sm, ax=ax_m, orientation="horizontal", location="bottom", shrink=0.8,
                      pad=0.01, aspect=25)
    cb.set_ticks([t for t in (0, 100, 250, 500, 1000, 2500) if t <= dmax])
    cb.set_label("m³/m (square-root scale)", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)

    # (b) the whole island, each model year
    ax_y = fig.add_subplot(gs[0, 1])
    ax_y.bar(years, by_year, width=0.6, color=C_YEAR)
    ax_y.set_xlim(years[0] - 0.6, years[-1] + 0.6)
    ax_y.set_xticks(years[::2] if len(years) > 14 else years)
    ax_y.tick_params(axis="x", labelsize=6.5)
    ax_y.set_ylim(0, ymax)
    ax_y.set_yticks([0, ymax / 2, ymax])
    ax_y.set_ylabel("10³ m³")
    ax_y.grid(axis="y", color="0.9", lw=0.5)
    open_frame(ax_y)
    _title(ax_y, 1, f"Whole island, each year: {by_year.sum() / 1e3:.1f} × 10⁶ m³ in all")

    out = save(fig, out_path)
    plt.close(fig)
    top = np.argsort(by_dom)[::-1][:3] + 1
    record_caption(out[0],
        f"Net overwash in the model over storm years {period}, {run_desc}. Net: Barrier3D's overwash "
        "less the sand the beach-dune manager cleared off the road and filtered, in a run that manages "
        "(the run's QowTS); in a run that does not, all of it. "
        "(a) Each domain's total over the run, in m3 per metre of shore, on a square-root scale; grey: "
        "the run never overwashed it. NC-12 in black. "
        f"(b) The whole island's overwash in each model year (each domain's m3/m times its "
        f"{DOM.domain_spacing_m:.0f} m length, summed over GIS 1-90), {by_year.sum() / 1e3:.2f} x 10^6 m3 "
        f"in all. (c) The totals of (a) as bars; black squares mark domains with washover in at least "
        f"one image the run's storms cover ({len(seen)} domains). Brackets mark Buxton, Avon and the "
        f"Tri-Village. The most: GIS {', '.join(map(str, top))}. Every run of this window is drawn on the "
        f"same scales ({dmax:.0f} m3/m, {ymax:.0f} x 10^3 m3 a year). GIS 1 is Cape Point, GIS 90 Pea Island.")
    return out


# One run's figure and table in its own folder; q is the run's QowTS (model years x padded domains)
def run_figure(q, run_dir, run_name, start_year):
    by_dom, by_year = totals(np.asarray(q, float))
    first = ovm.first_model_year(start_year)
    years = np.arange(first, first + len(by_year))
    run_dir = Path(run_dir)
    csv = write_path(run_dir, "overwash_total_csv", run_name)
    csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(dict(gis=np.arange(1, 91), net_overwash_m3_per_m=by_dom)).to_csv(csv, index=False)
    with plt.rc_context():          # the house style for this figure only, not the caller's later ones
        apply_style()
        out = draw(by_dom, by_year, years, seen_domains(start_year), load_geometry(),
                   write_path(run_dir, "figure_overwash_total", run_name),
                   f"run {run_name}", f"{years[0]}-{years[-1]}",
                   scales(int(start_year), by_dom, by_year))
    print(f"  overwash total        {by_year.sum() / 1e3:.2f} x 10^6 m3, most in GIS "
          f"{int(np.argmax(by_dom)) + 1} ({by_dom.max():.0f} m3/m)  -> {out[0].name}")
    return out


# A finished run's figure from its folder: the archive, or for runs saved without one its _overwash.npy
def run_figure_from_dir(run_dir):
    run_dir = Path(run_dir)
    name = run_dir.name
    meta = json.loads((run_dir / f"{name}_run_metadata.json").read_text(encoding="utf-8"))
    archive = run_dir / f"{name}.npz"
    q = ovm.qow(load_saved_cascade(archive)) if archive.exists() else np.load(run_dir / f"{name}_overwash.npy")
    return run_figure(q, run_dir, name, meta["period"]["start_year"])


# Run: each run folder named on the command line
def main():
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--run-dir", nargs="+", type=Path, required=True, help="run folder(s)")
    for rd in ap.parse_args().run_dir:
        print(rd.name)
        run_figure_from_dir(rd)


if __name__ == "__main__":
    main()
