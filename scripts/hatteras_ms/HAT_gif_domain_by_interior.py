#!/usr/bin/env python3
r"""
HAT_gif_domain_by_interior.py
==============================================================================
One domain, every interior of the row-insert set, one frame per model year.

WHAT A FRAME SHOWS
    One panel per interior: the near-dune part of the
                domain (the first CROSS_ROWS interior rows, 10 m cells), with
                Barrier3D's two dune rows in front of it, NC-12 outlined at its
                current setback and width, and the outline flashed in the year
                the module relocates it.

CROSS-SHORE FRAME
    Landward-positive metres from the YEAR-0 dune line of each interior. The
    dune line moves landward through the run (x_s_TS), so the island slides
    right across frames while a road that has not been relocated stays put --
    the mechanism the relocation module acts on. Different interiors start
    with row 0 at different absolute positions (the insert moves it N cells
    seaward), so 0 is each interior's OWN year-0 dune line, not a shared
    coordinate; what is comparable across panels is the road-to-dune distance
    and the elevations.

ELEVATIONS
    DomainTS[t] is the interior in dam above MHW; the dune rows are BermEl +
    DuneDomain[t]. The colour scale is the project's shared discrete one
    (hat_figure_style.elevation_cmap), water at or below MHW in blue.

READS
    the relocation-OFF run of each arm under
    output/raw_runs/row-insert/<arm>/1984_2004/calibBE/ (the model .npz).
    No re-run.

USAGE
    python HAT_gif_domain_by_interior.py                 # GIS 85, all arms
    python HAT_gif_domain_by_interior.py --domain 86 --arms none,median,platform
    python HAT_gif_domain_by_interior.py --stills 1984,1989,1995,2004
==============================================================================
"""
from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(HERE))
from hat_figure_style import elevation_cmap  # noqa: E402
from hatteras_site_config import HATTERAS_DOMAINS  # noqa: E402
from cascade_pipeline.run_registry import preset_dir_for  # noqa: E402
from HAT_run_row_insert_set import ARMS  # noqa: E402

# Only the two unmodified arms remain in ARMS (the layers and the set's runs
# were deleted 2026-09-07; see HAT_run_row_insert_set.py).
LABEL = {
    "original": "v1 original picks",
    "none": "v2 no insert",
}
COLOUR = {"original": "0.6", "none": "0.25"}
SET = "row-insert"
RAW = REPO / "output" / "raw_runs"
OUT = REPO / "output" / "experiments" / "row_insert_set" / "gifs"
START_YEAR = 1984
CELL = 10.0
CROSS_ROWS = 40          # interior rows drawn (400 m landward of the dune)
FPS = 1


def load_arm(arm):
    d = preset_dir_for(RAW, "1984_2004", "calibBE", arm="{}/{}".format(SET, arm))
    runs = [p for p in d.glob("*") if p.is_dir() and "_reloc_" not in p.name
            and list(p.glob("*.npz"))]
    if len(runs) != 1:
        raise SystemExit("expected one relocation-off run under {}, found {}".format(
            d, [r.name for r in runs]))
    npz = next(runs[0].glob("*.npz"))
    return np.load(npz, allow_pickle=True)["cascade"][0]


def year_state(c, i, t):
    """(raster m MHW with dune rows first, x0 of raster in m, road (x, ele, w, flash))."""
    b3d = c.barrier3d[i]
    interior = np.asarray(b3d.DomainTS[t], dtype=float)[:CROSS_ROWS] * CELL
    dune = (float(b3d.BermEl) + np.asarray(b3d.DuneDomain[t], dtype=float)) * CELL  # (along, 2)
    raster = np.vstack([dune.T[::-1], interior])          # seaward dune row first
    xs = np.asarray(b3d.x_s_TS, dtype=float)
    shift = (xs[t] - xs[0]) * CELL                          # landward-positive m
    x0 = shift - 2 * CELL                                  # first raster row = seaward dune row
    road = None
    if bool(np.asarray(c.roadway_management_module)[i]):
        m = c.roadways[i]
        sb = np.asarray(m._road_setback_TS, dtype=float)
        ele = np.asarray(m._road_ele_TS, dtype=float)
        wid = np.asarray(m._road_width_TS, dtype=float)
        rel = np.asarray(m._road_relocated_TS, dtype=float)
        if t < sb.size and ele[t] > 0:
            road = (shift + sb[t], float(ele[t]), float(wid[t]), bool(rel[t] > 0))
    return raster, x0, road


def render(states, arms, D, year, x_lim, cmap, norm, bounds, fig_w=16.0):
    n = len(arms)
    fig = plt.figure(figsize=(fig_w, 5.6))
    gs = fig.add_gridspec(1, n, wspace=0.12, top=0.84, bottom=0.18)
    fig.text(0.01, 0.975, "GIS {}  by interior".format(D), fontsize=15, fontweight="bold", va="top")
    fig.text(0.99, 0.975, "{}".format(year), fontsize=22, fontweight="bold", va="top", ha="right")
    fig.text(0.01, 0.935, "cross-shore: landward-positive metres from each interior's own "
             "year-0 dune line; each panel = first {} interior rows + the dune; "
             "box = NC-12 (red = relocated this year)".format(CROSS_ROWS),
             fontsize=8.5, va="top", color="#333333")
    for k, arm in enumerate(arms):
        raster, x0, road = states[arm]
        ax = fig.add_subplot(gs[0, k])
        nrow, nal = raster.shape
        ax.imshow(raster, cmap=cmap, norm=norm, aspect="auto", origin="upper",
                  extent=(-0.5, nal - 0.5, x0 + nrow * CELL, x0))
        ax.set_ylim(x_lim[1], x_lim[0])
        ax.set_title(LABEL[arm], fontsize=9.5, color=COLOUR[arm], loc="left")
        ax.set_xticks([])
        if k == 0:
            ax.set_ylabel("cross-shore (m)")
        else:
            ax.set_yticklabels([])
        ax.axhline(x0 + 2 * CELL, color="k", lw=.6, alpha=.6)    # interior row 0
        if road is not None:
            rx, rele, rw, flash = road
            ax.add_patch(Rectangle((-0.5, rx), nal, rw, fill=False,
                                   ec="#d62728" if flash else "k",
                                   lw=2.2 if flash else 1.2))
    # colour key, in the space the profile panel used to take
    cax = fig.add_axes([0.35, 0.085, 0.3, 0.02])
    inner = [b for b in bounds if abs(b) < 50]     # drop the sentinel ends
    fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax,
                 orientation="horizontal", ticks=inner)
    cax.set_xlabel("elevation (m MHW); leftmost class is water", fontsize=8)
    cax.tick_params(labelsize=7)
    return fig


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=int, default=85)
    ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--stills", default="",
                    help="comma list of years to also save as PNG")
    args = ap.parse_args()
    arms = [a.strip() for a in args.arms.split(",") if a.strip()]
    D = args.domain
    i = HATTERAS_DOMAINS.gis_to_pad(D)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cmap, norm, bounds = elevation_cmap()

    cas = {arm: load_arm(arm) for arm in arms}
    nt = len(np.asarray(cas[arms[0]].barrier3d[i].x_s_TS))
    # a fixed cross-shore window for the whole run: from just seaward of the
    # year-0 dune to CROSS_ROWS landward of the final dune line
    max_shift = max((np.asarray(cas[a].barrier3d[i].x_s_TS)[nt - 1]
                     - np.asarray(cas[a].barrier3d[i].x_s_TS)[0]) * CELL for a in arms)
    x_lim = (-3 * CELL, max_shift + CROSS_ROWS * CELL)

    stills = {int(y) for y in args.stills.split(",") if y.strip()}
    frames = []
    for t in range(nt):
        year = START_YEAR + t
        states = {arm: year_state(cas[arm], i, t) for arm in arms}
        fig = render(states, arms, D, year, x_lim, cmap, norm, bounds)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=80)
        if year in stills:
            fig.savefig(out / "GIS{}_by_interior_{}.png".format(D, year), dpi=130)
        plt.close(fig)
        buf.seek(0)
        frames.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE))
    path = out / "GIS{}_by_interior.gif".format(D)
    frames[0].save(path, save_all=True, append_images=frames[1:],
                   duration=int(1000 / FPS), loop=0)
    print("wrote {}  ({} frames, {} fps)".format(path, len(frames), FPS))


if __name__ == "__main__":
    main()
