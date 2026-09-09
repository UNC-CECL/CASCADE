#!/usr/bin/env python3
r"""
HAT_plot_footprint_result.py
==============================================================================
The first RESULT of the 1984 footprint: the hindcast on v3 (rows behind the
road, copy fill, 1984 setbacks) against the same run on v2, at the road.

Reads the two runs' saved model state (the roadway objects in the .npz) and
draws, for every road domain, the setback the model started with and every
year it relocated NC-12, v2 above v3. The island-wide skill of both runs is
in the caption, not on the figure.

    v2   output/raw_runs/1984_2004/calibBE/<run>              the calibration tree
    v3   output/raw_runs/behindroad-copy/1984_2004/calibBE/<run>   arm behindroad-copy

USAGE
    python HAT_plot_footprint_result.py
==============================================================================
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "0-elevation" / "3-figures"))
from hat_topo_version import insert_figures_dir  # noqa: E402
import HAT_plot_duneline_offset as off  # noqa: E402

NAME = "HAT_1984_2004_calibBE_road_bdm_groin"
RUNS = {"v2": REPO / "output/raw_runs/1984_2004/calibBE" / NAME,
        "v3": REPO / "output/raw_runs/behindroad-copy/1984_2004/calibBE" / NAME}
BUFFER = 15
START = 1984
FIG_DIR = insert_figures_dir("1984-start", "6-result", "island")
C = {"v2": "0.45", "v3": off.C_1984}
INK = off.INK
EVENTS = {1989: (84, 87), 1999: (9, 14)}       # the recorded NC-12 relocations


def load(d: Path):
    meta = json.load(open(d / f"{NAME}_run_metadata.json"))
    c = np.load(next(d.glob("*.npz")), allow_pickle=True)["cascade"][0]
    out = {}
    for g in range(1, 91):
        mgr = c.roadways[g + BUFFER - 1]
        if mgr is None:
            continue
        sb = np.asarray(getattr(mgr, "_road_setback_TS", []), float).ravel()
        rel = np.asarray(getattr(mgr, "_road_relocated_TS", []), float).ravel()
        out[g] = dict(sb0=float(sb[0]) if sb.size else np.nan,
                      years=[START + int(k) for k in np.flatnonzero(np.nan_to_num(rel) > 0)])
    return meta, out


def main() -> None:
    off.apply_style()
    data = {k: load(d) for k, d in RUNS.items()}
    doms = sorted(set(data["v2"][1]) | set(data["v3"][1]))
    fig, (a, b) = plt.subplots(2, 1, figsize=(13.0, 7.4), sharex=True, constrained_layout=True,
                               gridspec_kw=dict(height_ratios=[0.9, 1.1]))
    ann = off.HATTERAS_ANNOTATIONS
    for ax in (a, b):
        for name, (lo, hi) in ann.town_spans.items():
            ax.axvspan(lo - .5, hi + .5, color="0.93", zorder=0)
        for yr, (lo, hi) in EVENTS.items():
            ax.axvspan(lo - .5, hi + .5, facecolor="none", edgecolor="#2c6e49", lw=0.9, ls=(0, (3, 2)), zorder=1)
    for name, (lo, hi) in ann.town_spans.items():
        a.text((lo + hi) / 2, 0.985, name, transform=a.get_xaxis_transform(), ha="center", va="top",
               fontsize=7.5, color=off.INK_MUTED)

    # (a) the setback the model starts with
    for k, dx in (("v2", -0.2), ("v3", 0.2)):
        r = data[k][1]
        a.bar([g + dx for g in doms], [r.get(g, {}).get("sb0", np.nan) for g in doms], 0.4, color=C[k], zorder=3,
              label=f"{k}: setback at 1984")
    a.axhline(0, color=INK, lw=0.6)
    a.set_yscale("symlog", linthresh=50, linscale=1.0)
    a.set_yticks([0, 10, 20, 50, 100, 200, 500])
    a.set_yticklabels(["0", "10", "20", "50", "100", "200", "500"])
    a.set_ylabel("NC-12 setback at 1984\n(m landward of row 0)")
    a.grid(axis="y", color="0.92", lw=0.5)
    a.set_axisbelow(True)
    a.legend(loc="upper center", ncol=2)
    off._title(a, 0, "the setback the model starts with: today's (v2, floored at 0) and the 1984 one (v3)")

    # (b) every relocation the model made, by year
    for k, dx, mk in (("v2", -0.18, "o"), ("v3", 0.18, "s")):
        r = data[k][1]
        xs, ys = [], []
        for g in doms:
            for y in r.get(g, {}).get("years", []):
                xs.append(g + dx); ys.append(y)
        b.plot(xs, ys, mk, ms=4.2, color=C[k], mec="white", mew=0.5, zorder=4, linestyle="none",
               label=f"{k}: relocation ({len(xs)} in all)")
    for yr, (lo, hi) in EVENTS.items():
        b.hlines(yr, lo - .5, hi + .5, color="#2c6e49", lw=1.6, zorder=3)
        b.text(hi + 0.7, yr, f"NC-12 relocated {yr}", va="center", fontsize=7, color="#2c6e49")
    b.set_ylim(START - 0.5, START + 20.5)
    b.set_yticks(range(START, START + 21, 4))
    b.set_ylabel("model year of an emergent relocation")
    b.set_xlabel("domain (1 = south, Cape Hatteras)")
    b.set_xlim(0.2, 90.8)
    b.set_xticks([1] + list(range(10, 91, 10)))
    b.set_xticks(doms, minor=True)
    b.grid(axis="y", color="0.92", lw=0.5)
    b.set_axisbelow(True)
    h, l = b.get_legend_handles_labels()
    h += [Line2D([0], [0], color="#2c6e49", lw=1.6, label="recorded relocation, its block outlined")]
    b.legend(handles=h, loc="upper left", ncol=3, fontsize=7.5)
    off._title(b, 1, "when the model relocated NC-12, v2 against v3")

    p = FIG_DIR / "HAT_footprint_v3_vs_v2_road.png"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    for k in RUNS:
        m = data[k][0]
        print(f"{k}: topo {m['identity'].get('topo_dune_version')}  interior RMSE {m['skill'].get('rmse_interior_m_yr')}  "
              f"bias {m['skill'].get('mean_bias_interior_m_yr')}  drowned {m['verification'].get('roads_drowned')}  "
              f"relocations {sum(len(v['years']) for v in data[k][1].values())}")
    print(f"wrote {p}")


if __name__ == "__main__":
    main()
