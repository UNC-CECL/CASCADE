"""
How the source/sink (BE) end rates are solved and what zone field the runs carry.

    python scripts/figure_making/pipeline/7-source-sink/be_method_figures.py

Reads the end-domain solve experiments and the site config's BE field; writes
to output/figures/3-model-inputs/7-source-sink/. Details: scripts/figure_making/pipeline/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-29
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


# --- CONFIG ------------------------------------------------------------------
REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer.hatteras_site_config import HATTERAS_BE_EDGE_ONLY  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, C, C_1984, C_1997, INK, INK_MUTED, DOMAIN_AXIS_LABEL, figsize, figure_dir, save,
    record_caption, _title, open_frame,
)
OUT = figure_dir("inputs", "7-source-sink")
EDB = REPO / "output" / "raw_runs" / "experiments" / "end-domain-boundaries"
# The end solves the runs carry, oldest to newest (history in README)
ADOPTED = EDB / "2026-09-28-ends-resolved-adopted"
DUNECAP = EDB / "2026-09-28-ends-resolved-dunecap"
SPLIT12 = EDB / "2026-09-29-ends-resolved-split12"
WAVES = "Hs2.0_period7.5_asymmetry0.6_highangle0.5"
ZONES = REPO / "data" / "hatteras_init" / "7-source-sink" / "2-calibrate" / "1996_2010__2010_2024"
TOL = 0.02
# -----------------------------------------------------------------------------


# The 2010-2024 target at GIS 1: the raw CoastSat domain mean
def _gis1_target_2010():
    from cascade_pipeline.coastsat_lowess import CoastSatDataset, LowessConfig, build_coastsat_series
    from cascade_pipeline.hindcast import build_target_table
    from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS
    ds = CoastSatDataset(label="CoastSat LRR (2010-2024)", period_start=2010,
                         csv_path=str(COASTSAT_LRR_ROOT / "2010_2024" / "transect_lrr_full.csv"))
    cfg = LowessConfig(window_domains=(7,), skip_southern_domains=10)
    cs = build_coastsat_series([ds], active_period_start=2010, lowess_config=cfg,
                               domains=HATTERAS_DOMAINS)[0]
    return float(build_target_table(cs, cfg, HATTERAS_DOMAINS, 7)
                 .set_index("gis_domain")["target_lrr_m_yr"].loc[1])


# The adopted solve's direct GIS 1 probes: imposed rate and residual
def _direct_probes():
    target = _gis1_target_2010()
    rows = []
    for d in sorted((ADOPTED / "runs").glob("direct_gis1_*")):
        rate = next(d.rglob("tables/shoreline_change_rate.csv"), None)
        if rate is None:
            continue
        lrr = pd.read_csv(rate).set_index("gis_domain")["lrr_m_yr"].loc[1]
        rows.append(dict(imposed=float(d.name.split("_")[-1]), residual=float(lrr) - target))
    return pd.DataFrame(rows).sort_values("imposed").reset_index(drop=True)


# The end-rate solve, step by step, for both windows
def fig_end_solve():
    adopted = json.loads((ADOPTED / "tables" / "ends.json").read_text())["ends_m_yr"]
    dunecap = json.loads((DUNECAP / "tables" / "ends.json").read_text())["ends_m_yr"]
    carried = json.loads((SPLIT12 / "tables" / "ends.json").read_text())["ends_m_yr"]
    carried = {int(y): e for y, e in carried.items()}
    for y, e in carried.items():     # the values the runs carry must be these records'
        cfg = HATTERAS_BE_EDGE_ONLY[y]
        if abs(cfg[0] - e["1"]) > 1e-3 or abs(cfg[1] - e["90"]) > 1e-3:
            raise RuntimeError(f"{y}: config ends {cfg} differ from the solve record {e}")
    log = pd.read_csv(ADOPTED / "tables" / f"solve_log_1996_2010_{WAVES}.csv")
    cap = pd.read_csv(DUNECAP / "tables" / f"solve_log_2010_{WAVES}.csv")
    probes = _direct_probes()
    n_adopt = int(log[log.period_start == 2010].step.max())
    cap = cap.assign(step=cap.step + n_adopt)       # drawn after the adopted steps
    # split12: the probe log (step 1 = the trim24 ends) plus the accepted step, from the ends table
    s12 = pd.read_csv(SPLIT12 / "tables" / f"solve_log_1996_2010_{WAVES}.csv")
    s12_end = pd.read_csv(SPLIT12 / "tables" / f"ends_1996_2010_{WAVES}.csv")
    s12 = pd.concat([s12, pd.DataFrame(dict(
        period_start=s12_end.period.str[:4].astype(int), step=s12_end.steps,
        gis1_imposed=s12_end.gis1, gis90_imposed=s12_end.gis90,
        gis1_residual=s12_end.gis1_residual, gis90_residual=s12_end.gis90_residual))])
    row_last = {1996: int(log[log.period_start == 1996].step.max()),
                2010: int(max(cap.step.max(), n_adopt + len(probes)))}

    fig = plt.figure(figsize=figsize("double", height=5.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 3)
    for i, y in enumerate((1996, 2010)):
        lg = log[log.period_start == y]
        colr = C_1984 if y == 1996 else C_1997
        for j, (col, gis) in enumerate([("gis1", 1), ("gis90", 90)]):
            ax = fig.add_subplot(gs[i, j])
            ax.axhspan(-TOL, TOL, color="0.92", lw=0)
            ax.axhline(0, color=INK_MUTED, lw=0.5)
            ax.plot(lg.step, lg[f"{col}_residual"], "o-", color=colr, ms=4, lw=1.1,
                    label="secant step, adopted model")
            if y == 2010 and gis == 1 and len(probes):
                xs = range(n_adopt + 1, n_adopt + 1 + len(probes))
                ax.plot(list(xs), probes.residual, "s", color=C["ACCENT"], ms=4.5, label="direct probe")
            if y == 2010 and gis == 90:
                ax.plot(cap.step, cap.gis90_residual, "^-", color=INK, ms=4.5, lw=1.0,
                        label="after the dune-cap fix")
            sg = s12[s12.period_start == y]
            ax.plot(sg.step + row_last[y], sg[f"{col}_residual"], "D-", color=INK_MUTED, ms=4, lw=1.0,
                    mfc="white", label="split12 storms")
            ax.set_ylabel("model − target LRR (m/yr)")
            if i == 1:
                ax.set_xlabel("step")
            open_frame(ax)
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
            _title(ax, 2 * i + j, f"{y} · GIS {gis}")
            if i == 1:
                ax.legend(frameon=False, fontsize=7, loc="center left" if j == 0 else "best")
    ax = fig.add_subplot(gs[:, 2])
    lg = log[log.period_start == 2010]
    ax.axhspan(-TOL, TOL, color="0.92", lw=0)
    ax.axhline(0, color=INK_MUTED, lw=0.5)
    ax.plot(lg.gis1_imposed, lg.gis1_residual, "o", color=C_1997, ms=4, label="secant step")
    ax.plot(probes.imposed, probes.residual, "s", color=C["ACCENT"], ms=4.5, label="direct probe")
    sg = s12[s12.period_start == 2010]
    ax.plot(sg.gis1_imposed, sg.gis1_residual, "D", color=INK_MUTED, ms=4, mfc="white",
            label="split12 storms")
    ax.axvline(carried[2010]["1"], color=INK, lw=0.8, ls=(0, (3, 2)))
    ax.set_xlabel("GIS 1 source/sink imposed (m/yr)")
    ax.set_ylabel("GIS 1 residual (m/yr)")
    ax.legend(frameon=False, fontsize=7)
    open_frame(ax)
    _title(ax, 4, "2010 · GIS 1 response")
    out = save(fig, OUT / "be_end_solve.png")
    plt.close(fig)
    best = probes.iloc[(probes.imposed - carried[2010]["1"]).abs().argmin()]
    record_caption(out[0],
        "How the two end-domain source/sink values that the current (edgeBE) runs carry were solved, on the "
        "adopted model (Barrier3D hatteras/adopted, storms v3_trim24, then re-solved on v3_split12_trim24) at wave option A (Hs 2.0 m, Tp 7.5 s, "
        "asymmetry 0.6, high-angle fraction 0.5), metres island offset, full-management run, against each "
        "window's CoastSat LRR (GIS 1 the raw domain mean, GIS 90 the LOWESS-7 value; see observed_target). "
        "Only GIS 1 and GIS 90 carry a rate: they are boundary-artefact absorbers at the open ends of the reach, "
        "not a sediment budget. (a, b) 1996-2010, (c, d) 2010-2024: the residual, model minus target LRR, "
        f"after each secant step from the previous (LOWESS-7, pre-adoption) ends; the band is ±{TOL} m/yr. "
        f"1996 converged in four steps to GIS 1 {adopted['1996']['1']:+.4f} and GIS 90 {adopted['1996']['90']:+.4f} m/yr. "
        "In 2010 the secant stalled at GIS 1, whose residual stays near +0.8 m/yr for any rate from about +12 "
        "upward, so GIS 1 was set by direct probes (purple squares) with GIS 90 held at its secant value; "
        f"(e) shows why: the residual against the rate imposed, crossing zero near +8.0 m/yr "
        f"(probe residual {best.residual:+.3f}). The dune-cap fix of 2026-09-28 changed the managed "
        "interior, not the ends, and GIS 90 of 2010 was re-solved on it with GIS 1 held (black triangles in d), "
        f"to {dunecap['2010']['90']:+.4f} m/yr (it was {adopted['2010']['90']:+.4f}). "
        "On the split12 storms of 2026-09-29 (open diamonds) the secant, seeded at those ends, moved only GIS 1, "
        f"to {carried[1996]['1']:+.4f} (1996) and {carried[2010]['1']:+.4f} (2010; dashed in e) m/yr; GIS 90 held "
        f"at {carried[1996]['90']:+.4f} and {carried[2010]['90']:+.4f}. These are the values the runs carry. An imposed edge rate is "
        "roughly ten times the misfit it closes, because BRIE diffuses most of it alongshore "
        "(be_edge_domain_solve.py). Records: output/raw_runs/experiments/end-domain-boundaries/"
        "2026-09-28-ends-resolved-adopted/, .../2026-09-28-ends-resolved-dunecap/ and "
        ".../2026-09-29-ends-resolved-split12/.")
    return out


# The source/sink field the runs carry, by zone
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


# Run: draw both figures
def main():
    apply_style()
    print(fig_end_solve()[0].relative_to(REPO))
    print(fig_zone_field()[0].relative_to(REPO))


if __name__ == "__main__":
    main()
