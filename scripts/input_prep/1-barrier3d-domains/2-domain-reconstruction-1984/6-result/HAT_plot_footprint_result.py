#!/usr/bin/env python3
r"""
HAT_plot_footprint_result.py
==============================================================================
The first RESULT of the 1984 footprint: the hindcast on v3 (rows behind the
road, copy fill, 1984 setbacks) against the same run on v2, at the road.

Reads the two runs' saved model state (the roadway objects in the .npz) and
draws, for every road domain, the setback the model started with and every
year it relocated NC-12, the two versions together. The island-wide skill of
both runs is in the caption beside the figure, not on it: the figure is drawn
double-column in the house style of hat_figure_style (v2 grey C["BASE"], v3
purple C["ACCENT"], the recorded events C["REF"]).

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
from cascade_pipeline.run_layout import resolve as resolve_run_file  # noqa: E402
from site_layer.hat_topo_version import insert_figures_dir  # noqa: E402
from cascade_pipeline.run_registry import find_run_dir  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C as STYLE_C, INK, DOMAIN_AXIS_LABEL, apply_style, figsize, open_frame, record_caption, save, town_bands, _title,
)

NAME = "HAT_1984_2004_calibBE_road_bdm_groin"
# Asked for rather than spelled: the arm moved under arms/ on 2026-09-10 and
# a hand-built path missed it. find_run_dir reads either layout.
_RAW_RUNS = REPO / "output" / "raw_runs"
RUNS = {"v2": find_run_dir(_RAW_RUNS, NAME, "1984_2004", "calibBE"),
        "v3": find_run_dir(_RAW_RUNS, NAME, "1984_2004", "calibBE",
                           "behindroad-copy")}
BUFFER = 15
START = 1984
FIG_DIR = insert_figures_dir("1984-start", "6-result", "island")
VS = ("v2", "v3")
C = {"v2": STYLE_C["BASE"], "v3": STYLE_C["ACCENT"]}       # the unmodified input, the modification under test
C_REF = STYLE_C["REF"]                                     # the recorded events
MK = {"v2": "o", "v3": "s"}
LABEL = {"v2": "as extracted (1996 surface)", "v3": "1984 reconstruction"}
EVENTS = {1989: (84, 87), 1999: (9, 14)}       # the recorded NC-12 relocations


def load(d: Path):
    # RESOLVED, NOT JOINED: run_layout knows where a run folder keeps each of
    # its files, in the new layout and the old flat one alike.
    meta = json.load(open(resolve_run_file(d, "metadata_json", NAME)))
    c = np.load(resolve_run_file(d, "archive", NAME, must_exist=True),
                allow_pickle=True)["cascade"][0]
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
    apply_style()
    data = {k: load(d) for k, d in RUNS.items()}
    doms = sorted(set(data["v2"][1]) | set(data["v3"][1]))
    fig, (a, b) = plt.subplots(2, 1, figsize=figsize("double", height=5.4), sharex=True,
                               constrained_layout=True, gridspec_kw=dict(height_ratios=[0.9, 1.1]))
    town_bands(a)
    town_bands(b, label=False)
    for ax in (a, b):
        for yr, (lo, hi) in EVENTS.items():
            ax.axvspan(lo - .5, hi + .5, facecolor="none", edgecolor=C_REF, lw=0.8, ls=(0, (3, 2)), zorder=1)

    # (a) the setback the model starts with: one thin line with markers per version
    for k in VS:
        r = data[k][1]
        a.plot(doms, [r.get(g, {}).get("sb0", np.nan) for g in doms], ls="-", lw=0.9, marker=MK[k], ms=2.6,
               mew=0, color=C[k], zorder=3, label=LABEL[k])
    a.axhline(0, color=INK, lw=0.5)
    a.set_yscale("symlog", linthresh=20, linscale=0.8)
    a.set_yticks([0, 10, 20, 50, 100, 200, 500])
    a.set_yticklabels(["0", "10", "20", "50", "100", "200", "500"])
    a.set_ylim(-2, 1500)
    a.set_ylabel("NC-12 setback at 1984\n(m landward of row 0)")
    a.grid(axis="y")
    a.set_axisbelow(True)
    open_frame(a)
    a.legend(loc="lower center", ncol=2, frameon=False, handlelength=2.2)
    _title(a, 0, "the setback the model starts with")

    # (b) every relocation the model made, by year
    handles = []
    for k, dx in (("v2", -0.18), ("v3", 0.18)):
        r = data[k][1]
        xs, ys = [], []
        for g in doms:
            for y in r.get(g, {}).get("years", []):
                xs.append(g + dx); ys.append(y)
        b.plot(xs, ys, MK[k], ms=4.0, color=C[k], mec="white", mew=0.5, zorder=4, linestyle="none")
        handles.append(Line2D([0], [0], marker=MK[k], ms=4.0, color=C[k], mec="white", mew=0.5, ls="none",
                              label=f"{LABEL[k]}: {len(xs)} relocations"))
    for yr, (lo, hi) in EVENTS.items():
        b.hlines(yr, lo - .5, hi + .5, color=C_REF, lw=1.6, zorder=3)
        right = hi < 60
        b.text(hi + 1.6 if right else lo - 1.6, yr, f"recorded {yr}", va="center",
               ha="left" if right else "right", fontsize=7, color=C_REF, zorder=6)
    handles.append(Line2D([0], [0], color=C_REF, lw=1.6, label="recorded relocation, its block outlined"))
    b.set_ylim(START - 0.5, START + 20.5)
    b.set_yticks(range(START, START + 21, 4))
    b.set_ylabel("year of a relocation")
    b.set_xlabel(DOMAIN_AXIS_LABEL)
    b.set_xlim(0.2, 90.8)
    b.set_xticks([1] + list(range(10, 91, 10)))
    b.grid(axis="y")
    b.set_axisbelow(True)
    open_frame(b)
    b.legend(handles=handles, loc="upper center", ncol=2, frameon=False)
    _title(b, 1, "when the model relocated NC-12")

    p = FIG_DIR / "HAT_footprint_v3_vs_v2_road.png"
    save(fig, p, bbox_inches="tight")
    m = {k: data[k][0] for k in RUNS}
    n = {k: sum(len(v["years"]) for v in data[k][1].values()) for k in RUNS}
    record_caption(p, "The 1984-2004 calibBE full-management hindcast on the 1984 reconstruction (purple; "
                      "dune-topo v3, the footprint behind NC-12, copy fill, 1984 setbacks) against the same run on "
                      "the domains as extracted from the 1996 surface (grey; v2, setbacks floored at 0), at the "
                      "road. (a) The setback each run starts with, per road domain, on a symmetric-log axis (linear "
                      "below 20 m); the reconstruction's are the 1984 measurements and none is floored. (b) Every "
                      "relocation the model made on its own, by year and domain; the recorded 1989 Pea Island and "
                      "1999 inter-village relocations are green bars across their outlined blocks. The run relocates "
                      f"{n['v2']} times as extracted and {n['v3']} times on the reconstruction. Interior shoreline "
                      f"RMSE {m['v2']['skill'].get('rmse_interior_m_yr')} (as extracted) and "
                      f"{m['v3']['skill'].get('rmse_interior_m_yr')} (reconstruction) m/yr; roads drowned "
                      f"{m['v2']['verification'].get('roads_drowned')} and "
                      f"{m['v3']['verification'].get('roads_drowned')}. Domain 1 is Cape Point, 90 Pea Island; "
                      "villages banded. Run `" + NAME + "`.")
    plt.close(fig)
    for k in RUNS:
        mm = data[k][0]
        print(f"{k}: topo {mm['identity'].get('topo_dune_version')}  interior RMSE {mm['skill'].get('rmse_interior_m_yr')}  "
              f"bias {mm['skill'].get('mean_bias_interior_m_yr')}  drowned {mm['verification'].get('roads_drowned')}  "
              f"relocations {n[k]}")
    print(f"wrote {p}")


if __name__ == "__main__":
    main()
