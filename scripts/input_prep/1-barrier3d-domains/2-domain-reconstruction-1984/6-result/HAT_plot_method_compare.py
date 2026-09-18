#!/usr/bin/env python3
r"""
HAT_plot_method_compare.py
==============================================================================
One domain, two methodologies, side by side in the model's frame:

    v2   (default base) the re-pick extraction v3 is built on. The 1984 road
         was measured against interior row 0, came out NEGATIVE (seaward of
         row 0, in the dune), and was FLOORED to 0 - so the model placed NC-12
         on rows 0-1, at the dune, and relocated it in year 1. No rows added.
         Same pick set as v3, so the panels differ ONLY by the rows and the
         setback. (--base v1 shows the original extraction instead, which
         also differs by the 2026-09-02 re-pick.)
    v3   the re-pick base (v2) + the symmetric 1984 footprint placed directly
         behind the road AS PLACED under its 1984 setback (road against the
         1984 dune line, row-0 convention; no floor) and filled by copying the
         rows that follow. The road sits on measured cells; the block is behind it.

Each panel: the two dune rows (berm + dune height) on top, the interior below
in elevation classes, a metres axis on the right; NC-12 as the model places it
(dark band); the measured 1984 road position (outlined); and the inserted rows
outlined in the accent colour. The interior depth, the retreat the block stands
for and what each version's hindcast did with the road (v1: arm
pea1989basenoreloc; v3: arm behindroad-copy; both calibBE, full management,
prescribed relocations off) go to the CAPTIONS.md beside the figure, not onto
the canvas; the figure is drawn double-column in the house style of
hat_figure_style.

USAGE
    python HAT_plot_method_compare.py                 # GIS 85
    python HAT_plot_method_compare.py --domain 86 --rows 40
==============================================================================
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
INIT = REPO / "data" / "hatteras_init"
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "0-elevation" / "3-figures"))
from site_layer.hat_topo_version import array_name, dune_topo_root, insert_figures_dir_for_domain, topo_dirs, insert_scope_step# noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C as STYLE_C, INK, apply_style, elevation_cmap, figsize, record_caption, save, spines_for_image, _title,
)

PRODUCT = "1984-start"
CELL_M = 10.0
ROAD_ROWS = 2
DUNE_ROWS = 2
BERM_EL_M = 1.7
BUFFER = 15
START = 1984
RECORDED = {84: 1989, 85: 1989, 86: 1989, 87: 1989, 9: 1999, 10: 1999, 11: 1999, 12: 1999, 13: 1999, 14: 1999}
SCOPE = INIT / "1-barrier3d-domains" / PRODUCT / "2-domain-reconstruction-1984"
RUN = "HAT_1984_2004_calibBE_road_bdm_groin"
RUNS = {"v1": REPO / "output/raw_runs/pea1989basenoreloc/1984_2004/calibBE" / RUN,
        "v2": REPO / "output/raw_runs/1984_2004/calibBE" / RUN,
        "v3": REPO / "output/raw_runs/behindroad-copy/1984_2004/calibBE" / RUN}
from site_layer import hat_topo_version as _tv  # noqa: E402
MEASURED = {  # the unfloored measurement each version's CSV was floored from
    # v1 named dunestart_offset_ARCHIVE_1984start_v1/, which became a dated
    # superseded folder; resolved since 2026-09-18.
    "v1": _tv.SETBACK_1984_V1_DIR / "RoadOffset_1984_domains.csv",
    "v2": _tv.road_setback_dir(1984) / "RoadOffset_1984_domains.csv",
    "v3": _tv.road_setback_dir(1984) / "RoadOffset_1984_domains.csv",
}
C_ADD = STYLE_C["ACCENT"]        # the modification under test: the inserted rows
C_ROAD = STYLE_C["ROAD"]         # NC-12 as the model places it
C_OLD = STYLE_C["REF"]           # the measured 1984 road position, an observation


def model_setback(version: str, d: int) -> float:
    rows = list(csv.reader(open(dune_topo_root(PRODUCT) / version / "RoadSetback_1984_dunestart.csv", newline="")))
    ids = [int(float(x)) for x in rows[0] if x.strip()]
    vals = [float(x) for x in rows[1] if x.strip()]
    return dict(zip(ids, vals))[d]


def measured_setback(version: str, d: int) -> float:
    return float(pd.read_csv(MEASURED[version]).set_index("domain").loc[d, "setback_dunestart_m"])


def stack(version: str, d: int):
    topo, dune, _ = topo_dirs(PRODUCT, override=version)
    z = np.load(topo / array_name("topography", d)) * CELL_M
    dn = np.load(dune / array_name("dune", d)) * CELL_M + BERM_EL_M
    return np.concatenate([np.tile(dn[None, :], (DUNE_ROWS, 1)), z], axis=0), z.shape[0]


def run_years(version: str, d: int):
    """Emergent relocation years the hindcast on this version produced, or None."""
    p = RUNS.get(version)
    if p is None or not p.is_dir():
        return None
    c = np.load(next(p.glob("*.npz")), allow_pickle=True)["cascade"][0]
    mgr = c.roadways[d + BUFFER - 1]
    if mgr is None:
        return []
    rel = np.asarray(getattr(mgr, "_road_relocated_TS", []), float).ravel()
    return [START + int(k) for k in np.flatnonzero(np.nan_to_num(rel) > 0)]


def draw_road(ax, y: float, ncol: int, label: str) -> None:
    ax.add_patch(Rectangle((-0.5, y - 0.5), ncol, ROAD_ROWS, facecolor=C_ROAD, edgecolor=C_ROAD,
                           lw=1.2, alpha=0.5, zorder=5))
    ax.text(ncol - 1.0, y + ROAD_ROWS / 2 - 0.5, label, ha="right", va="center", fontsize=8,
            fontweight="bold", color="white", zorder=6)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--domain", type=int, default=85)
    ap.add_argument("--rows", type=int, default=40)
    ap.add_argument("--base", choices=("v2", "v1"), default="v2",
                    help="v2 (default): same pick set as v3, so only the rows and the setback differ. "
                         "v1: the original extraction, which also differs by the 2026-09-02 re-pick")
    args = ap.parse_args()
    d = args.domain
    apply_style()
    cmap, norm, bounds = elevation_cmap()
    fp = pd.read_csv(insert_scope_step(PRODUCT, "2-extent") / "footprint_1984_by_domain.csv").set_index("domain").loc[d]
    n = int(fp["n_cells"])
    ins = int(fp["insert_row_behind_road"]) if n != 0 else -1
    shift = float(fp["shift_m_median"])

    base_title = {"v1": "as first extracted (1996 surface)",
                  "v2": "as extracted (1996 surface)"}[args.base]
    panels = [(args.base, base_title),
              ("v3", "1984 reconstruction")]
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=5.2), constrained_layout=True)
    numbers = {}
    for k, (ax, (ver, title)) in enumerate(zip(axes, panels)):
        img, rows_total = stack(ver, d)
        sb_meas, sb_model = measured_setback(ver, d), model_setback(ver, d)
        R = min(img.shape[0], args.rows + DUNE_ROWS)
        ncol = img.shape[1]
        ax.imshow(cmap(norm(img[:R])), aspect="auto", interpolation="nearest", origin="upper",
                  extent=[-0.5, ncol - 0.5, R - 0.5, -0.5])
        ax.axhline(DUNE_ROWS - 0.5, color=INK, lw=1.0, zorder=5)
        ax.text(ncol - 1.0, DUNE_ROWS / 2 - 0.5, "dune", ha="right", va="center", fontsize=8,
                fontweight="bold", color="white", zorder=6)

        # the road as the model places it: int(setback/10) rows behind row 0
        y_model = DUNE_ROWS + sb_model // CELL_M
        draw_road(ax, y_model, ncol, f"NC-12 as placed: {sb_model:.0f} m")
        # the measured 1984 position where it differs: a short box at the left
        y_meas = DUNE_ROWS + sb_meas / CELL_M
        if abs(y_meas - y_model) > 0.05:
            wbox = 12
            ax.add_patch(Rectangle((-0.5, y_meas - 0.5), wbox, ROAD_ROWS, facecolor="none", edgecolor="white",
                                   lw=2.6, zorder=6))
            ax.add_patch(Rectangle((-0.5, y_meas - 0.5), wbox, ROAD_ROWS, facecolor="none", edgecolor=C_OLD,
                                   lw=1.2, ls=(0, (2, 2)), zorder=7))
            ax.text(wbox + 0.5, y_meas + ROAD_ROWS / 2 - 0.5,
                    f"measured: {sb_meas:+.0f} m" + (", floored to 0" if sb_model == 0 and sb_meas < 0 else ""),
                    ha="left", va="center", fontsize=8, color=C_OLD,
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.15"), zorder=7)

        # the inserted rows, with a bracket for the retreat they stand for
        if ver == "v3" and n > 0:
            y0, y1 = DUNE_ROWS + ins - 0.5, DUNE_ROWS + ins + n - 0.5
            ax.add_patch(Rectangle((-0.5, y0), ncol, n, facecolor="none", edgecolor=C_ADD, lw=1.4, zorder=6))
            ax.text(0.5, y1 + 0.3, f"+{n} rows inserted (copy fill)",
                    ha="left", va="top", fontsize=8, color=C_ADD, fontweight="bold",
                    bbox=dict(facecolor="white", alpha=0.88, edgecolor="none", boxstyle="square,pad=0.2"), zorder=7)
        numbers[ver] = dict(rows_total=rows_total, sb_model=sb_model, sb_meas=sb_meas,
                            years=run_years(ver, d))

        ax.set_xlim(-0.5, ncol - 0.5)
        ax.set_xlabel("alongshore cell")
        if k == 0:
            ax.set_ylabel("cross-shore cell (0 = the dune)")
        ax.set_yticks(range(0, R, 5))
        sec = ax.secondary_yaxis("right", functions=(lambda c: (c - DUNE_ROWS) * CELL_M,
                                                     lambda m: m / CELL_M + DUNE_ROWS))
        if k == len(panels) - 1:            # one metres label, at the right edge
            sec.set_ylabel("m landward of interior row 0")
        sec.set_yticks(range(0, int((R - DUNE_ROWS) * CELL_M), 100))
        spines_for_image(ax)
        _title(ax, k, title)

    labels = ["below 0 (water)"] + [f"{lo:g}–{hi:g}" for lo, hi in zip(bounds[1:-2], bounds[2:-1])] \
        + [f"above {bounds[-2]:g}"]
    handles = [Patch(facecolor=cmap(i), edgecolor="0.4", lw=0.4, label=lab) for i, lab in enumerate(labels)]
    handles += [Patch(facecolor=C_ROAD, alpha=0.5, edgecolor=C_ROAD, label="NC-12 as the model places it"),
                Patch(facecolor="none", edgecolor=C_OLD, lw=1.2, ls=(0, (2, 2)), label="measured 1984 road position"),
                Line2D([0], [0], color=C_ADD, lw=1.4, label="rows inserted (copy fill)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=6, frameon=False,
               title="elevation classes (m MHW); the dune rows are drawn at berm + dune height")
    p = insert_figures_dir_for_domain(PRODUCT, "6-result", d) / f"HAT_method_compare_{args.base}_v3_GIS{d}.png"
    save(fig, p, vector=False, bbox_inches="tight")
    plt.close(fig)

    def _road(v):
        b = numbers[v]
        floored = " (the measurement floored to 0)" if b["sb_model"] == 0 and b["sb_meas"] < 0 else ""
        did = ("relocated " + ", ".join(str(y) for y in b["years"]) if b["years"] else "no relocation") \
            if b["years"] is not None else "no hindcast on this version"
        return (f"{b['rows_total']} interior rows, the road at {b['sb_model']:.0f} m{floored} on rows "
                f"{int(b['sb_model'] // CELL_M)}–{int(b['sb_model'] // CELL_M) + 1}, measured "
                f"{b['sb_meas']:+.0f} m; the hindcast with prescribed relocations off {did}")

    base_name = {"v1": "the original extraction", "v2": "the extraction the reconstruction is built on"}[args.base]
    rec = RECORDED.get(d)
    record_caption(p, f"GIS {d} under the two methodologies, in the model's frame: the two dune rows on top (berm "
                      "plus dune height), the interior below in elevation classes (m above MHW), cells on the left "
                      "axis and metres landward of interior row 0 on the right; NC-12 as the model places it is the "
                      "dark band, the measured 1984 position the outlined one. (a) The domains as extracted from the "
                      f"1996 surface ({base_name}, dune-topo {args.base}): {_road(args.base)}. (b) The 1984 "
                      f"reconstruction (v3): the same extraction with {n} rows inserted directly behind the road at "
                      f"interior row {ins} and filled by copying the {n} rows that follow them, standing for the "
                      f"{shift:.0f} m the dune line retreated between 1984 and 1997 in whole cells, and the road at "
                      f"its 1984 setback; {_road('v3')}."
                      + (f" The recorded NC-12 relocation here is {rec}." if rec else "")
                      + (" Both panels share the 2026-09-02 pick set, so they differ only by the block and the "
                         "setback." if args.base == "v2" else " The base is the original extraction, so the "
                         "2026-09-02 re-pick is part of the difference as well."))
    print(f"wrote {p}")


if __name__ == "__main__":
    main()
