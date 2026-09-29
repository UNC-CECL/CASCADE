"""
be_method_figures.py
==============================================================================
The source/sink (background erosion, BE) method figures for the current
1996-2010 / 2010-2024 pair, which existed only for 1984/2004.

    python scripts/figure_making/pipeline/7-source-sink/be_method_figures.py

Writes output/figures/3-model-inputs/7-source-sink/:
    be_end_solve.png   how the two end values the current (edgeBE) runs carry
                       were solved: residual per Newton step at GIS 1 and 90,
                       and the direct probes where 2010 GIS 1 stopped
                       converging
    be_zone_field.png  the zone-by-zone calibration (calibBE) for the pair:
                       residuals, which domains were eligible, and the rates.
                       Made 2026-09-18, BEFORE the metres offset and wave
                       option A; the current runs do not use it

WHAT EXISTS AND WHAT DOES NOT
    The 1984/2004 convergence figure reads convergence_history.json, the
    per-pass RMSE of the iterated zone calibration. The 1996/2010 pair has no
    such file: its zone field (2-calibrate/1996_2010__2010_2024/) is a single
    pass. So no zone-iteration convergence figure can be drawn for it, and
    none is invented. What the current runs DO carry is the edge-only preset,
    solved by Newton steps on 2026-09-27; its step-by-step record is the solve
    log of output/raw_runs/experiments/end-domain-boundaries/
    2026-09-27-ends-resolved-metres-offset/tables/, which be_end_solve.png draws.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))

from site_layer.hatteras_site_config import HATTERAS_BE_EDGE_ONLY  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, C, C_1984, C_1997, INK, INK_MUTED, DOMAIN_AXIS_LABEL, figsize, figure_dir, save,
    record_caption, _title, open_frame,
)

OUT = figure_dir("inputs", "7-source-sink")
SOLVE = (REPO / "output" / "raw_runs" / "experiments" / "end-domain-boundaries"
         / "2026-09-27-ends-resolved-metres-offset" / "tables")
# The ends the runs carry since 2026-09-28: re-solved against the LOESS-7 target,
# seeded at the 09-27 values drawn here (only GIS 90, and 2010 GIS 1 slightly, moved).
CURRENT = (REPO / "output" / "raw_runs" / "experiments" / "end-domain-boundaries"
           / "2026-09-28-ends-resolved-loess7" / "tables")
WAVES = "Hs2.0_period7.5_asymmetry0.6_highangle0.5"
ZONES = REPO / "data" / "hatteras_init" / "7-source-sink" / "2-calibrate" / "1996_2010__2010_2024"
TOL = 0.05


def fig_end_solve():
    ends = json.loads((SOLVE / "ends.json").read_text())["ends_m_yr"]
    now = json.loads((CURRENT / "ends.json").read_text())["ends_m_yr"]
    logs = {1996: pd.read_csv(SOLVE / f"solve_log_1996_{WAVES}.csv"),
            2010: pd.concat([pd.read_csv(SOLVE / f"solve_log_2010_{WAVES}.csv"),
                             pd.read_csv(SOLVE / f"solve_log_2010_probes_{WAVES}.csv")], ignore_index=True)}
    for y in logs:        # the values the runs carry must be the current solve's
        cfg = HATTERAS_BE_EDGE_ONLY[y]
        if abs(cfg[0] - now[str(y)]["1"]) > 1e-3 or abs(cfg[1] - now[str(y)]["90"]) > 1e-3:
            raise RuntimeError(f"{y}: config ends {cfg} differ from the solve record {now[str(y)]}")

    fig = plt.figure(figsize=figsize("double", height=5.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 3)
    for i, (y, lg) in enumerate(logs.items()):
        n_solve = len(pd.read_csv(SOLVE / f"solve_log_{y}_{WAVES}.csv"))
        for j, (col, gis) in enumerate([("gis1", 1), ("gis90", 90)]):
            ax = fig.add_subplot(gs[i, j])
            ax.axhspan(-TOL, TOL, color="0.92", lw=0)
            ax.axhline(0, color=INK_MUTED, lw=0.5)
            solve = lg.iloc[:n_solve]
            probe = lg.iloc[n_solve:]
            colr = C_1984 if y == 1996 else C_1997
            ax.plot(solve.step, solve[f"{col}_residual"], "o-", color=colr, ms=4, lw=1.1, label="Newton step")
            if len(probe):
                ax.plot(probe.step, probe[f"{col}_residual"], "s", color=C["ACCENT"], ms=4.5, label="direct probe")
            ax.set_ylabel("model − target LRR (m/yr)")
            if i == 1:
                ax.set_xlabel("step")
            open_frame(ax)
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
            _title(ax, 2 * i + j, f"{y} · GIS {gis}")
            if i == 1 and j == 0:
                ax.legend(frameon=False, fontsize=7)
    ax = fig.add_subplot(gs[:, 2])
    lg = logs[2010]
    ax.axhspan(-TOL, TOL, color="0.92", lw=0)
    ax.axhline(0, color=INK_MUTED, lw=0.5)
    ax.plot(lg.gis1_imposed, lg.gis1_residual, "o", color=C_1997, ms=4)
    for _, r in lg.iterrows():
        ax.annotate(f"{int(r.step)}", (r.gis1_imposed, r.gis1_residual), xytext=(3, 3),
                    textcoords="offset points", fontsize=6.5, color=INK_MUTED)
    ax.axvline(ends["2010"]["1"], color=INK, lw=0.8, ls=(0, (3, 2)))
    ax.set_xlabel("GIS 1 source/sink imposed (m/yr)")
    ax.set_ylabel("GIS 1 residual (m/yr)")
    open_frame(ax)
    _title(ax, 4, "2010 · GIS 1 response")
    out = save(fig, OUT / "be_end_solve.png")
    plt.close(fig)
    record_caption(out[0],
        "How the two end-domain source/sink values that the current (edgeBE) runs carry were solved, "
        "2026-09-27, at wave option A (Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle fraction 0.5) on the "
        "metres island offset, full-management run, against each window's CoastSat LRR (GIS 1 the raw domain "
        "mean, GIS 90 the LOESS-10 value, the target then; see observed_target). Only GIS 1 and GIS 90 carry a rate: they are "
        "boundary-artefact absorbers at the open ends of the reach, not a sediment budget. (a, b) 1996-2010, "
        "(c, d) 2010-2024: the residual, model minus target LRR, after each safeguarded Newton step; the band "
        f"is ±{TOL} m/yr. 1996 converged in four steps to GIS 1 {ends['1996']['1']:+.4f} and GIS 90 "
        f"{ends['1996']['90']:+.3f} m/yr. In 2010 GIS 90 converged ({ends['2010']['90']:+.3f}) but GIS 1 stalled "
        "at +0.21 m/yr, and was set by direct probes (purple squares). (e) Why: the 2010 GIS 1 residual against "
        "the rate imposed there, every step and probe labelled with its number; the response is not monotonic "
        "within ~0.1 m/yr, so the closest probe was taken, "
        f"{ends['2010']['1']:+.1f} m/yr (dashed; residual -0.045). An imposed edge rate is roughly ten times the "
        "misfit it closes, because BRIE diffuses most of it alongshore (be_edge_domain_solve.py). Record: "
        "output/raw_runs/experiments/end-domain-boundaries/2026-09-27-ends-resolved-metres-offset/. "
        "On 2026-09-28 the target moved to LOESS 7 and the ends were re-solved from these values "
        f"(end-domain-boundaries/2026-09-28-ends-resolved-loess7/): 1996 GIS 90 {now['1996']['90']:+.4f}, "
        f"2010 GIS 1 {now['2010']['1']:+.4f} and GIS 90 {now['2010']['90']:+.4f} m/yr, residuals within "
        "0.02 m/yr; those are the values the runs now carry.")
    return out


def fig_zone_field():
    m = pd.read_csv(ZONES / "be_zone_metrics.csv")
    g = m.domain.to_numpy()
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.4), sharex=True, constrained_layout=True)
    ax = axes[0]
    zones = m.groupby((m.physical_zone != m.physical_zone.shift()).cumsum())
    for k, (_, z) in enumerate(zones):
        ax.axvspan(z.domain.min() - 0.5, z.domain.max() + 0.5, color="0.95" if k % 2 else "white", lw=0,
                   zorder=0)
        axes[1].axvspan(z.domain.min() - 0.5, z.domain.max() + 0.5, color="0.95" if k % 2 else "white",
                        lw=0, zorder=0)
        name = z.physical_zone.iloc[0].replace("–", "-").replace("�", "-").replace(" / ", "/\n")
        name = name.replace(" Influence", "\ninfluence").replace(" Transition", "\ntransition")
        ax.text(z.domain.mean(), 0.98, name, transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=6.3, color=INK_MUTED, linespacing=1.0)
    for p, col, lab in [("p1", C_1984, "1996-2010"), ("p2", C_1997, "2010-2024")]:
        ax.plot(g, m[f"smooth_residual_{p}"], color=col, lw=1.2, label=lab)
        w = m[f"correction_warranted_{p}"].astype(bool)
        ax.plot(g[w], m.loc[w, f"smooth_residual_{p}"], "o", color=col, ms=3)
    ax.axhline(0, color=INK_MUTED, lw=0.5)
    ax.set_ylim(None, ax.get_ylim()[1] + 1.6)
    ax.set_ylabel("observed − base model LRR (m/yr)")
    open_frame(ax)
    ax.legend(frameon=False, fontsize=7.5, loc="center", bbox_to_anchor=(0.5, 0.74), ncol=2)
    _title(ax, 0, "")
    ax = axes[1]
    wdt = 0.4
    ax.bar(g - wdt / 2, m.be_hindcast_p1, width=wdt, color=C_1984, label="1996-2010")
    ax.bar(g + wdt / 2, m.be_hindcast_p2, width=wdt, color=C_1997, label="2010-2024")
    for strat, mk in [("groin-reserved", "x"), ("locked", "D")]:
        s = m.strategy == strat
        ax.plot(g[s], np.zeros(s.sum()), mk, color=C["ACCENT"] if strat == "locked" else C["ADDED"], ms=5,
                label=strat)
    ax.axhline(0, color=INK_MUTED, lw=0.5)
    ax.set_ylabel("source/sink rate (m/yr)")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_xlim(0.5, 90.5)
    open_frame(ax)
    fig.legend(*ax.get_legend_handles_labels(), loc="outside lower center", ncol=4, frameon=False, fontsize=7.5)
    _title(ax, 1, "")
    out = save(fig, OUT / "be_zone_field.png")
    plt.close(fig)
    record_caption(out[0],
        "The zone-by-zone source/sink calibration (calibBE) for the 1996-2010 / 2010-2024 pair, as made on "
        "2026-09-18 (7-source-sink/2-calibrate/1996_2010__2010_2024/): a single pass, BEFORE the island offset "
        "went to metres and before wave option A. The current matrix runs use the edge-only preset instead "
        "(be_end_solve), so this field is the method, not an input in use. (a) The residual, observed CoastSat "
        "LRR minus the base model's, smoothed, per domain and window, over the physical zones the calibration "
        "names (bands); points mark domains where a correction was judged warranted. (b) The rate each domain "
        "was given per window: shifting zones take the window's own residual, stable ones a shared value; "
        "GIS 1 and 90 are locked (purple diamonds; their values come from the separate end solve) and "
        "GIS 5-7 are reserved for the Buxton groin (amber crosses). No convergence figure exists for this "
        "pair: it has no convergence_history.json, unlike 1984/2004. GIS 1 is Cape Point, GIS 90 Pea Island.")
    return out


def main():
    apply_style()
    print(fig_end_solve()[0].relative_to(REPO))
    print(fig_zone_field()[0].relative_to(REPO))


if __name__ == "__main__":
    main()
