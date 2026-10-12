#!/usr/bin/env python3
"""
The runs from the ±1-yr start (matrix) and the ±6-month start (6mo), overlaid with their targets.

    python model_1yr_vs_6mo.py  ->  figures/model_1yr_vs_6mo.png
                                    figures/model_1yr_vs_6mo_edgeBE.png
                                    figures/model_1yr_vs_6mo_edgeBE_fit.png
                                    figures/model_1yr_vs_6mo_edgeBE_fit_overlay.png

Reads domain_net_change.csv (written by six_month_starts.py score); no runs.
Columns: calibration 1996-2009, test 2009-2025. Rows: zeroBE, edgeBE, domainBE + groin.
Red ±1 yr, blue ±6 months; solid the observed target from that start, dashed the run.
The _edgeBE figure is edgeBE alone on a tighter axis: per period, targets then runs.
The _edgeBE_fit figure is edgeBE run against its own target, one panel per period and start;
_fit_overlay puts both starts in one panel per period (thick pale target, thin run).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-10
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer.hat_figure_style import (  # noqa: E402
    C_1984, C_1997, DOMAIN_AXIS_LABEL, INK_MUTED, _title, apply_style, figsize,
    mark_offaxis, offaxis_clause, record_caption, save, town_bands)

# --- CONFIG ------------------------------------------------------------------
POS_HALF = 100.0                                 # m, the net-change figures' fixed axis
PERIODS = {"1996": "Calibration 1996-2009", "2009": "Test 2009-2025"}
# The end window every target in a period shares (the config's); only the start differs
END = {"1996": "1 Jul 2007 to 30 Jun 2009", "2009": "17 Feb 2025 to 17 Feb 2026"}
PRESETS = {"zeroBE": "zeroBE, no groin", "edgeBE": "edgeBE, no groin",
           "domainBE": "domainBE + groin"}
# arm -> (legend name, colour, start window per period)
ARMS = {"matrix": ("±1 yr start", C_1984, {"1996": "12 Oct 1995 to 12 Oct 1997",
                                           "2009": "17 Aug 2008 to 17 Aug 2010"}),
        "6mo": ("±6 month start", C_1997, {"1996": "12 Apr 1996 to 12 Apr 1997",
                                           "2009": "17 Feb 2009 to 17 Feb 2010"})}
PNG = HERE / "figures" / "model_1yr_vs_6mo.png"
EDGE_PNG = HERE / "figures" / "model_1yr_vs_6mo_edgeBE.png"
FIT_PNG = HERE / "figures" / "model_1yr_vs_6mo_edgeBE_fit.png"
OVERLAY_PNG = HERE / "figures" / "model_1yr_vs_6mo_edgeBE_fit_overlay.png"
EDGE_HALF = 60.0                                 # m, the edgeBE figure's tighter axis
# -----------------------------------------------------------------------------


def main():
    nc = pd.read_csv(HERE / "domain_net_change.csv", index_col="gis_domain")
    sc = pd.read_csv(HERE / "scores.csv")
    x = nc.index.to_numpy()
    apply_style()
    fig, axes = plt.subplots(len(PRESETS), len(PERIODS), sharex=True, sharey=True,
                             figsize=figsize("double", height=7.0), constrained_layout=True)
    off = []
    for i, (preset, run) in enumerate(PRESETS.items()):
        for j, (yr, title) in enumerate(PERIODS.items()):
            a = axes[i, j]
            a.axhline(0, color=INK_MUTED, lw=0.6)
            a.set_ylim(-POS_HALF, POS_HALF)
            a.set_xlim(x.min(), x.max())
            a.grid(axis="y")
            score_txt = []
            for arm, (name, c, _) in ARMS.items():
                tgt, mod = nc[f"{yr}_target_{arm}_lowess7_m"], nc[f"{yr}_{preset}_{arm}_lowess7_m"]
                a.plot(x, tgt, color=c, lw=1.2, alpha=0.8)
                a.plot(x, mod, color=c, lw=1.2, ls="--")
                if i == 0:                       # the targets are the same in every row
                    off.append((f"{title} {name} target", mark_offaxis(a, x, tgt.values, POS_HALF, color=c)))
                off.append((f"{title} {run} {name} run", mark_offaxis(a, x, mod.values, POS_HALF, color=c)))
                s = sc[(sc.period.str.startswith(yr)) & (sc.run == preset) & (sc.arm == arm)].iloc[0]
                score_txt.append(f"{name}: bias {s.bias_m:+.1f}, RMSE {s.rmse_m:.1f}, r {s.r:.2f}")
            a.text(0.01, 0.03, "\n".join(score_txt), transform=a.transAxes, ha="left",
                   va="bottom", fontsize=6.5, color=INK_MUTED)
            _title(a, i * len(PERIODS) + j, f"{title}, {run}")
            town_bands(a, where="top", strip=0.06, label=(i == 0))
            if j == 0:
                a.set_ylabel("Net change (m)")
            if i == len(PRESETS) - 1:
                a.set_xlabel(DOMAIN_AXIS_LABEL)
    handles = [Line2D([], [], color=c, lw=1.2, label=f"{name} ({w['1996']}; {w['2009']})")
               for name, c, w in ARMS.values()]
    handles += [Line2D([], [], color="0.3", lw=1.2, label="solid: observed target from that start"),
                Line2D([], [], color="0.3", lw=1.2, ls="--", label="dashed: CASCADE run from that start")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, fontsize=7, frameon=False)
    save(fig, PNG, dpi=300, close=True)
    record_caption(PNG, (
        "Net shoreline change per GIS domain from the two candidate start lines. Red: the "
        "±1-yr CoastSat mean around the DEM date (offset v2, the production matrix runs); blue: "
        "the ±6-month mean (offset v3, the 6mo runs). Solid: the observed target, the end window "
        "(calibration 1 Jul 2007 to 30 Jun 2009, test 17 Feb 2025 to 17 Feb 2026) minus the "
        "mean over that start window; dashed: the run started from that line. Scores (m, "
        "interior GIS 2-89) in each panel. All at 7-domain LOWESS. Seaward positive."
        + offaxis_clause(off, POS_HALF, "m")))
    print(f"wrote {PNG}")
    edge_only(nc, sc)
    edge_fit(nc, sc)
    edge_fit_overlay(nc, sc)


# edgeBE alone: per period, the two targets in one panel and the two runs in the next
def edge_only(nc, sc):
    x = nc.index.to_numpy()
    apply_style()
    fig, axes = plt.subplots(4, 1, sharex=True, sharey=True,
                             figsize=figsize("double", height=8.5), constrained_layout=True)
    off = []
    for i, (yr, title) in enumerate(PERIODS.items()):
        for k, (what, col) in enumerate([("observed target", "target"),
                                         ("CASCADE edgeBE run", "edgeBE")]):
            a = axes[2 * i + k]
            a.axhline(0, color=INK_MUTED, lw=0.6)
            a.set_ylim(-EDGE_HALF, EDGE_HALF)
            a.set_xlim(x.min(), x.max())
            a.set_xticks(range(5, 91, 5), minor=True)
            a.grid(axis="y")
            a.grid(axis="x", which="both", color="0.93", lw=0.5)
            for arm, (name, c, w) in ARMS.items():
                y = nc[f"{yr}_{col}_{arm}_lowess7_m"]
                a.plot(x, y, color=c, lw=1.4, label=f"from {name} mean ({w[yr]})")
                off.append((f"{title} {name} {what}", mark_offaxis(a, x, y.values, EDGE_HALF, color=c)))
            if col == "edgeBE":                  # scores belong with the run
                txt = []
                for arm, (name, _, _) in ARMS.items():
                    s = sc[(sc.period.str.startswith(yr)) & (sc.run == "edgeBE") & (sc.arm == arm)].iloc[0]
                    txt.append(f"{name}: bias {s.bias_m:+.1f} m, RMSE {s.rmse_m:.1f} m, r {s.r:.2f}")
                a.text(0.99, 0.03, "\n".join(txt), transform=a.transAxes, ha="right",
                       va="bottom", fontsize=7, color=INK_MUTED)
            a.legend(loc="lower center", fontsize=7, frameon=False)   # clear of the GIS 6 dip
            _title(a, 2 * i + k, f"{title}: {what}"
                   + (f", change to end mean {END[yr]}" if col == "target" else ""))
            town_bands(a, where="top", strip=0.07, label=True)
            a.set_ylabel("Net change (m)")
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    save(fig, EDGE_PNG, dpi=300, close=True)
    record_caption(EDGE_PNG, (
        "Net shoreline change per GIS domain from the two candidate start lines, edgeBE with no "
        "groin. (a) Calibration 1996-2009, observed target; (b) calibration, model run; (c) test "
        "2009-2025, observed target; (d) test, model run. Red: from the ±1-yr CoastSat mean "
        "around the DEM date (offset v2, the production matrix run); blue: from the ±6-month mean "
        "(offset v3, the 6mo run); start dates in each legend. Each target is the config's end "
        "window minus the mean over its own start window. Scores (m, interior GIS 2-89) under "
        "each run, each run against the target from its own start. All at 7-domain LOWESS. "
        f"Seaward positive. Axis ±{EDGE_HALF:.0f} m." + offaxis_clause(off, EDGE_HALF, "m")))
    print(f"wrote {EDGE_PNG}")


# edgeBE against its own target: one panel per period and start, two lines each
def edge_fit(nc, sc):
    x = nc.index.to_numpy()
    apply_style()
    fig, axes = plt.subplots(4, 1, sharex=True, sharey=True,
                             figsize=figsize("double", height=8.5), constrained_layout=True)
    off = []
    panels = [(yr, title, arm) for yr, title in PERIODS.items() for arm in ARMS]
    for n, (yr, title, arm) in enumerate(panels):
        name, c, w = ARMS[arm]
        a = axes[n]
        a.axhline(0, color=INK_MUTED, lw=0.6)
        a.set_ylim(-EDGE_HALF, EDGE_HALF)
        a.set_xlim(x.min(), x.max())
        a.set_xticks(range(5, 91, 5), minor=True)
        a.grid(axis="y")
        a.grid(axis="x", which="both", color="0.93", lw=0.5)
        tgt, mod = nc[f"{yr}_target_{arm}_lowess7_m"], nc[f"{yr}_edgeBE_{arm}_lowess7_m"]
        a.plot(x, tgt, color="0.15", lw=1.4, label=f"Observed change to end mean ({END[yr]})")
        a.plot(x, mod, color=c, lw=1.4, label="CASCADE edgeBE run")
        off.append((f"{title} {name} target", mark_offaxis(a, x, tgt.values, EDGE_HALF, color="0.15")))
        off.append((f"{title} {name} run", mark_offaxis(a, x, mod.values, EDGE_HALF, color=c)))
        s = sc[(sc.period.str.startswith(yr)) & (sc.run == "edgeBE") & (sc.arm == arm)].iloc[0]
        leg = a.legend(loc="lower center", fontsize=7, frameon=False, ncol=2, title_fontsize=7,
                       title=f"bias {s.bias_m:+.1f} m, RMSE {s.rmse_m:.1f} m, r {s.r:.2f}")
        leg.get_title().set_color(INK_MUTED)
        _title(a, n, f"{title}, from the {name} mean ({w[yr]})")
        town_bands(a, where="top", strip=0.07, label=True)
        a.set_ylabel("Net change (m)")
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    save(fig, FIT_PNG, dpi=300, close=True)
    record_caption(FIT_PNG, (
        "edgeBE (no groin) against the observed target from the same start line. (a) Calibration "
        "1996-2009 from the ±1-yr start, (b) from the ±6-month start; (c) test 2009-2025 from the "
        "±1-yr start, (d) from the ±6-month start. Black: observed net change, the config's end "
        "window minus the CoastSat mean over the start window; red (±1 yr, offset v2, the "
        "production matrix run) or blue (±6 months, offset v3, the 6mo run): the run started from "
        "that line. Scores (m, interior GIS 2-89) in each panel. All at 7-domain LOWESS. Seaward "
        f"positive. Axis ±{EDGE_HALF:.0f} m." + offaxis_clause(off, EDGE_HALF, "m")))
    print(f"wrote {FIT_PNG}")


# Both starts in one panel per period: hue is the start, weight is target (thick, pale) or run
def edge_fit_overlay(nc, sc):
    x = nc.index.to_numpy()
    apply_style()
    fig, axes = plt.subplots(2, 1, sharex=True, sharey=True,
                             figsize=figsize("double", height=6.5), constrained_layout=True)
    off = []
    for i, (yr, title) in enumerate(PERIODS.items()):
        a = axes[i]
        a.axhline(0, color=INK_MUTED, lw=0.6)
        a.set_ylim(-EDGE_HALF, EDGE_HALF)
        a.set_xlim(x.min(), x.max())
        a.set_xticks(range(5, 91, 5), minor=True)
        a.grid(axis="y")
        a.grid(axis="x", which="both", color="0.93", lw=0.5)
        handles = []
        for arm, (name, c, w) in ARMS.items():
            tgt, mod = nc[f"{yr}_target_{arm}_lowess7_m"], nc[f"{yr}_edgeBE_{arm}_lowess7_m"]
            s = sc[(sc.period.str.startswith(yr)) & (sc.run == "edgeBE") & (sc.arm == arm)].iloc[0]
            handles.append(a.plot(x, tgt, color=c, lw=4.0, alpha=0.3, solid_capstyle="round",
                                  label=f"Observed change from {name} mean ({w[yr]})")[0])
            handles.append(a.plot(x, mod, color=c, lw=1.3,
                                  label=f"edgeBE run, {name}: bias {s.bias_m:+.1f} m, "
                                        f"RMSE {s.rmse_m:.1f} m, r {s.r:.2f}")[0])
            off.append((f"{title} {name} target", mark_offaxis(a, x, tgt.values, EDGE_HALF, color=c)))
            off.append((f"{title} {name} run", mark_offaxis(a, x, mod.values, EDGE_HALF, color=c)))
        # one column over GIS ~38-68, clear of the GIS 35 and Rodanthe troughs
        a.legend(handles=handles[0::2] + handles[1::2], loc="lower center",
                 bbox_to_anchor=(0.58, 0.0), fontsize=6.5, frameon=False)
        _title(a, i, f"{title}, edgeBE, no groin: observed change to end mean {END[yr]}")
        town_bands(a, where="top", strip=0.06, label=True)
        a.set_ylabel("Net change (m)")
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    save(fig, OVERLAY_PNG, dpi=300, close=True)
    record_caption(OVERLAY_PNG, (
        "edgeBE (no groin) against the observed target, both start lines in one panel. "
        "(a) Calibration 1996-2009, (b) test 2009-2025. Colour is the start: red the ±1-yr "
        "CoastSat mean around the DEM date (offset v2, the production matrix run), blue the "
        "±6-month mean (offset v3, the 6mo run). Thick pale lines: the observed target from that "
        "start; thin lines: the run started from that line. Scores (m, interior GIS 2-89), each run "
        "against its own target, in the legend. All at 7-domain LOWESS. Seaward positive. "
        f"Axis ±{EDGE_HALF:.0f} m." + offaxis_clause(off, EDGE_HALF, "m")))
    print(f"wrote {OVERLAY_PNG}")


if __name__ == "__main__":
    main()
