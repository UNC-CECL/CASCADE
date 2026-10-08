#!/usr/bin/env python3
"""
Animate the three groin modules year by year on the 1996-2009 calibration run.

    python compare_groin_modules.py [--conserving b0.60_f0.6]   ->  figures/groin_module_comparison_1996_2009.gif + .png

The dipole (a fixed source/sink pair), the pinned blocking groin (traps a fraction of the
transport, each flank with its own diffusion number) and the conserving blocking groin
(one face coefficient, so GIS 6 gains exactly what GIS 5 loses). Each frame: shoreline
change from 1996 along GIS 1-12 for each module against no groin, the change each module
applied at GIS 5 and 6 that year with the running net, and the GIS 5|6 gap against the
photos. Every run: full management, the solved edgeBE ends, relocations off, failure
instant from the 2004 step.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-08
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path[:0] = [str(REPO / "hard-structures" / "groin" / "2-module-tests" / "3-real-planform"),
                str(REPO / "scripts" / "hatteras_ms" / "groin-sweep"),
                str(REPO / "scripts" / "hatteras_ms"), str(REPO / "scripts")]
from score_instant_grid import observed_series  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _title, apply_style, figsize, open_frame, record_caption, save)
from site_layer.hat_observed_rates import net_change_domain_csv  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RAW = REPO / "output" / "raw_runs" / "experiments" / "groin"
PIN = RAW / "2026-10-05-blocking-fit-dem-to-dem"
CONS = RAW / "2026-10-08-conserving-groin"
NOGROIN = (REPO / "output" / "raw_runs" / "matrix" / "1996_2009" / "edgeBE"
           / "HAT_1996_2009_edgeBE_offsetmetres_road_bdm_nogroin")
WINDOW = "1996_2009"
START = 1996
UP, DOWN = 15 + 6 - 1, 15 + 5 - 1          # GIS 6 / GIS 5 in the padded array
PAD = 15                                   # padded index of GIS 1
ZOOM = (1, 12)
FPS = 1.5
OUT = HERE / "figures" / "groin_module_comparison_1996_2009"
# -----------------------------------------------------------------------------


def run_dir(root):
    hits = sorted((root / WINDOW / "edgeBE").glob("*/"))
    if not hits:
        raise FileNotFoundError(f"no run under {root / WINDOW}")
    return hits[-1]


def matrix(rd):
    return np.load(next(rd.glob("*_shoreline_matrix.npy")))


# Seaward change the module applied each model year at GIS 6 and GIS 5
def applied(rd):
    d = pd.read_csv(rd / "tables" / "groin_diagnostics.csv")
    return -d["applied_dx_updrift_m"].to_numpy(), -d["applied_dx_downdrift_m"].to_numpy()


def versions(conserving):
    b, f = conserving[1:].split("_f")
    return (
        dict(name="Source/sink (dipole)", detail="M 12, f 0.3: fixed +M / -M",
             rd=run_dir(PIN / "dipole_M12_f0.3"), col=C["ADDED"]),
        dict(name="Trapping (pinned)", detail="b 0.6, f 0.6: own r each side",
             rd=run_dir(PIN / "b0.60_f0.6"), col=C["ACCENT"]),
        dict(name="Trapping, conserving", detail=f"b {float(b):.2f}, f {f}: one face r",
             rd=run_dir(CONS / conserving), col=C["REF"]),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--conserving", default="b0.60_f0.6", help="member folder of the conserving run")
    a = ap.parse_args()

    apply_style()
    vs = versions(a.conserving)
    base = matrix(NOGROIN)
    for v in vs:
        v["x"] = matrix(v["rd"])
        v["up"], v["down"] = applied(v["rd"])
    n = base.shape[0]                                   # 1 Jan 1996 .. 1 Jan 2009
    gis = np.arange(ZOOM[0], ZOOM[1] + 1)
    cols = PAD + gis - 1

    def change(x, k):
        return -(x[k, cols] - x[0, cols])               # seaward positive

    def gap(x):
        g = x[:, DOWN] - x[:, UP]
        return g - g[0]

    net_obs = pd.read_csv(net_change_domain_csv(START, START + n - 1),
                          index_col=0)["net_change_m"].loc[gis].to_numpy()
    lim_pos = max(np.abs(change(v["x"], k)).max() for v in [*vs, dict(x=base)] for k in range(n))
    lim_pos = np.ceil(max(lim_pos, np.abs(net_obs).max()) / 20) * 20
    lim_app = np.ceil(max(np.abs(np.r_[v["up"], v["down"]]).max() for v in vs) / 10) * 10
    obs = observed_series()
    g0 = np.interp(START, obs.index, obs.values)
    obs_yrs = [y for y in obs.index if START < y < START + n]
    gaps = [gap(base)] + [gap(v["x"]) for v in vs]
    lim_gap = np.ceil(max(np.abs(np.r_[np.concatenate(gaps),
                                        [obs[y] - g0 for y in obs_yrs]]).max(), 20) / 20) * 20

    fig = plt.figure(figsize=figsize("double", height=7.4))
    gs = fig.add_gridspec(3, 3, height_ratios=(1.25, 0.8, 1.0), hspace=0.55, wspace=0.28,
                          left=0.11, right=0.98, top=0.88, bottom=0.08)
    ax_pos = [fig.add_subplot(gs[0, i]) for i in range(3)]
    ax_app = [fig.add_subplot(gs[1, i]) for i in range(3)]
    ax_gap = fig.add_subplot(gs[2, :])

    def draw(k):
        year = START + k
        for i, v in enumerate(vs):
            ax = ax_pos[i]
            ax.clear()
            ax.axhline(0, color=INK_MUTED, lw=0.6)
            ax.axvline(5.5, color=C["GROIN"], lw=0.9)
            ax.plot(gis, change(base, k), color=C["BASE"], lw=1.3)
            ax.plot(gis, net_obs, color=INK, lw=0.9, ls=(0, (3, 2)), marker="o", ms=3.5,
                    mfc=INK, mec="white", mew=0.5, zorder=5)
            ax.plot(gis, change(v["x"], k), color=v["col"], lw=1.8)
            ax.set_xlim(*ZOOM)
            ax.set_ylim(-lim_pos, lim_pos)
            ax.set_xticks(gis[::1])
            ax.tick_params(labelsize=6.5)
            if i == 0:
                ax.set_ylabel("Shoreline change\nfrom 1996 (m, seaward +)")
            ax.set_xlabel("GIS domain (S to N)", fontsize=7)
            ax.grid(axis="y")
            open_frame(ax)
            _title(ax, i, v["name"])
            ax.text(0.03, 0.97, v["detail"], transform=ax.transAxes, ha="left", va="top",
                    fontsize=6.5, color=INK_MUTED)

            ax = ax_app[i]
            ax.clear()
            j = k - 1                                       # the model year just finished
            up = v["up"][j] if j >= 0 else 0.0
            dn = v["down"][j] if j >= 0 else 0.0
            net_cum = float(np.sum(v["up"][:k] + v["down"][:k])) if k else 0.0
            ax.bar([0, 1], [dn, up], color=[v["col"], v["col"]], width=0.6,
                   alpha=0.9, edgecolor="none")
            ax.axhline(0, color=INK_MUTED, lw=0.6)
            ax.set_xticks([0, 1], ["GIS 5\n(south)", "GIS 6\n(north)"], fontsize=6.5)
            ax.set_xlim(-0.6, 1.6)
            ax.set_ylim(-lim_app, lim_app)
            ax.tick_params(labelsize=6.5)
            if i == 0:
                ax.set_ylabel("Applied this year\n(m, seaward +)")
            ax.text(0.5, 0.96, f"net this year {up + dn:+.1f} m\nnet so far {net_cum:+.0f} m",
                    transform=ax.transAxes, ha="center", va="top", fontsize=6.5,
                    color=INK if abs(net_cum) < 0.5 else C["GROIN"])
            ax.grid(axis="y")
            open_frame(ax)

        ax = ax_gap
        ax.clear()
        t = START + np.arange(k + 1)
        ax.plot(t, gaps[0][:k + 1], color=C["BASE"], lw=1.3)
        for v, g in zip(vs, gaps[1:]):
            ax.plot(t, g[:k + 1], color=v["col"], lw=1.8)
        shown = [y for y in obs_yrs if y <= year]
        ax.plot(shown, [obs[y] - g0 for y in shown], "o", ms=5, mfc=INK, mec="white", mew=0.6,
                zorder=6)
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.axvline(2004, color=INK_MUTED, lw=0.6, ls=(0, (2, 2)))
        ax.text(2004.1, -lim_gap * 0.9, "2004: failure", fontsize=6.5, color=INK_MUTED, va="bottom")
        ax.set_xlim(START, START + n - 1)
        ax.set_ylim(-lim_gap, lim_gap)
        ax.set_ylabel("GIS 5|6 gap change\nfrom 1996 (m)")
        ax.set_xlabel("Year (1 January)")
        ax.grid(axis="y")
        open_frame(ax)
        _title(ax, 3, "The step at the groin; dots are the wet/dry photos")
        fig.suptitle(f"1 January {year}", fontsize=10, y=0.985)

    handles = [Line2D([], [], color=C["BASE"], lw=1.3, label="no groin")]
    handles += [Line2D([], [], color=v["col"], lw=1.8, label=v["name"]) for v in vs]
    handles += [Line2D([], [], color=C["GROIN"], lw=0.9, label="groin (GIS 5|6)"),
                Line2D([], [], color=INK, lw=0.9, ls=(0, (3, 2)), marker="o", ms=3.5, mfc=INK,
                       mec="white", label="CoastSat 1996-2009"),
                Line2D([], [], ls="", marker="o", ms=5, mfc=INK, mec="white", label="photos")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=7,
               fontsize=6.8, frameon=False)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    writer = PillowWriter(fps=FPS)
    with writer.saving(fig, str(OUT.with_suffix(".gif")), dpi=130):
        for k in range(n):
            draw(k)
            writer.grab_frame()
        for _ in range(2):                              # hold the last frame
            writer.grab_frame()
    png = OUT.with_suffix(".png")
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        "The three groin modules on the 1996-2009 calibration run (final frame of "
        f"{OUT.name}.gif). Full management, the solved edgeBE ends, relocations off, failure "
        "instant from the 2004 step. Columns: the source/sink dipole (M 12, f 0.3, the best "
        "option-A dipole), the pinned blocking groin (b 0.6, f 0.6), and the conserving blocking "
        f"groin ({a.conserving}). Top: shoreline change from 1 Jan 1996 along GIS 1-12 against no "
        "groin (grey), and the observed CoastSat net change 1996-2009 (dashed, domain means, shown "
        "on every frame as the end-of-run target); red line, the groin. Middle: the shoreline change the module applied at "
        "GIS 5 and GIS 6 in the model year just finished, with the net for that year and so far; "
        "a non-zero net means the module created or removed sand. Bottom: the GIS 5|6 gap change "
        "against the wet/dry photo gaps from the same start."))
    print(f"wrote {OUT.with_suffix('.gif').relative_to(REPO)} and {png.name}")


if __name__ == "__main__":
    main()
