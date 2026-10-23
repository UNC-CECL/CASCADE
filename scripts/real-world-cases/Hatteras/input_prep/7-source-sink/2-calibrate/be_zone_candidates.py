#!/usr/bin/env python3
"""
Step 2 of the source/sink zones: candidate zones from the step 1 residuals, before any physics.

    python scripts/input_prep/7-source-sink/2-calibrate/be_zone_candidates.py

Reads 2-calibrate/1996_2009__2009_2025/zone_residuals_edgeBE.csv (be_zone_residuals.py) and writes
zone_candidates_edgeBE.csv beside it and 3-figures/1996_2009__2009_2025/2-method/
zone_candidates_edgeBE.png.

A candidate zone is a run of at least MIN_WIDTH adjacent domains whose residual is beyond
THRESHOLD_M_YR with one sign. Three sets:
  shared       beyond the threshold in BOTH periods with the same sign (scenario A)
  calibration  beyond it in calibration alone (scenario B, 1996-2009)
  test         beyond it in test alone (scenario B, 2009-2025)
GIS 1 and 90 (solved end rates) and 5-7 (the groin, still undecided) are never candidates.
Each zone also records whether it still qualifies at the stricter STRICT_M_YR. Nothing is
fitted and nothing is screened for physical sense yet (step 3).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-12
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
import be_zone_residuals as step1  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK_MUTED, _title, apply_style, figsize, open_frame, record_caption,
    save, structures, town_bands)
from site_layer.hat_source_sink import calibrate_dir, figures_dir  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_PERIODS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
THRESHOLD_M_YR = 0.5            # 2026-10-12: the earlier zone method's significance threshold
STRICT_M_YR = 1.0               # reported per zone, not used to select
MIN_WIDTH = 3                   # domains (1.5 km): no correction narrower than BRIE can carry
EXCLUDED = {1: "solved end rate", 90: "solved end rate",
            5: "groin", 6: "groin", 7: "groin"}
TABLE = calibrate_dir(step1.PAIR) / "zone_candidates_edgeBE.csv"
PNG = figures_dir(step1.PAIR) / "2-method" / "zone_candidates_edgeBE.png"
# -----------------------------------------------------------------------------


# Runs of consecutive allowed domains where `flag` holds with one sign: [(sign, [gis, ...]), ...]
def runs(flag, sign):
    out, cur = [], []
    for g in flag.index:
        ok = bool(flag[g]) and g not in EXCLUDED
        if ok and cur and sign[g] == sign[cur[-1]] and g == cur[-1] + 1:
            cur.append(g)
            continue
        if len(cur) >= MIN_WIDTH:
            out.append((int(sign[cur[0]]), cur))
        cur = [g] if ok else []
    if len(cur) >= MIN_WIDTH:
        out.append((int(sign[cur[0]]), cur))
    return out


def candidates(r):
    cal, test = r[1996], r[2009]
    sets = {
        "shared": ((cal.abs() > THRESHOLD_M_YR) & (test.abs() > THRESHOLD_M_YR)
                   & (np.sign(cal) == np.sign(test)), np.sign(cal)),
        "calibration": (cal.abs() > THRESHOLD_M_YR, np.sign(cal)),
        "test": (test.abs() > THRESHOLD_M_YR, np.sign(test)),
    }
    rows = []
    for name, (flag, sign) in sets.items():
        for k, (s, gis) in enumerate(runs(flag, sign), start=1):
            c, t = cal.loc[gis], test.loc[gis]
            # what decides the strict test: both periods for shared, the one period otherwise
            judge = (pd.concat([c, t], axis=1).abs().min(axis=1) if name == "shared"
                     else (c if name == "calibration" else t).abs())
            rows.append(dict(
                set=name, zone=f"{name[0].upper()}{k}", kind="source" if s > 0 else "sink",
                first_gis=gis[0], last_gis=gis[-1], width=len(gis),
                cal_mean_m_yr=round(c.mean(), 3), test_mean_m_yr=round(t.mean(), 3),
                cal_peak_m_yr=round(c.loc[c.abs().idxmax()], 3),
                test_peak_m_yr=round(t.loc[t.abs().idxmax()], 3),
                domains_beyond_strict=int((judge > STRICT_M_YR).sum()),
                passes_strict=bool(len(runs(judge > STRICT_M_YR, pd.Series(s, index=judge.index))))))
    return pd.DataFrame(rows)


def figure(r, z):
    apply_style()
    panels = [("shared", "A. Shared: same sign beyond ±0.5 m/yr in both periods", (1996, 2009)),
              ("calibration", "B. Calibration only, 1996–2009", (1996,)),
              ("test", "B. Test only, 2009–2025", (2009,))]
    fig, axes = plt.subplots(3, 1, sharex=True, constrained_layout=True,
                             figsize=figsize("double", height=8.4))
    for k, (ax, (name, title, periods)) in enumerate(zip(axes, panels)):
        for row in z[z.set == name].itertuples():
            c = step1.C_SOURCE if row.kind == "source" else step1.C_SINK
            ax.axvspan(row.first_gis - 0.5, row.last_gis + 0.5, color=c, alpha=0.16, lw=0, zorder=1)
            ax.text((row.first_gis + row.last_gis) / 2, step1.YLIM[1] * 0.88, row.zone, ha="center", va="top",
                    fontsize=7, fontweight="bold", color=c, zorder=7)
        for g in EXCLUDED:
            ax.axvspan(g - 0.5, g + 0.5, facecolor="none", hatch="xxx", edgecolor="0.6", lw=0, zorder=1)
        for s in (-1, 1):
            ax.axhline(s * THRESHOLD_M_YR, color=INK_MUTED, lw=0.9, ls=":", zorder=2)
            ax.axhline(s * STRICT_M_YR, color=INK_MUTED, lw=0.8, ls="-.", zorder=2)
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        for p in periods:
            ax.plot(r.index, r[p], zorder=4, **step1.PERIOD_STYLE[p])
        ax.set_xlim(0.5, 90.5)
        ax.set_ylim(*step1.YLIM)
        ax.grid(axis="y")
        open_frame(ax)
        town_bands(ax, label=(k == 0))
        structures(ax, label=(k == 2))
        ax.set_ylabel("Residual (m/yr)")
        _title(ax, k, title)
        n = z[z.set == name]
        ax.set_title(f"{len(n)} zones, {int(n.width.sum())} domains; {int(n.passes_strict.sum())} "
                     f"also at ±{STRICT_M_YR:g}", loc="right", fontsize=7, color=INK_MUTED)
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    fig.legend(handles=[Line2D([], [], **step1.PERIOD_STYLE[1996], label="Calibration 1996–2009"),
                        Line2D([], [], **step1.PERIOD_STYLE[2009], label="Test 2009–2025"),
                        Patch(color=step1.C_SOURCE, alpha=0.3, label="Candidate source zone"),
                        Patch(color=step1.C_SINK, alpha=0.3, label="Candidate sink zone"),
                        Patch(facecolor="none", hatch="xxx", edgecolor="0.7",
                              label="Excluded (end domains, groin)"),
                        Line2D([], [], color=INK_MUTED, ls=":", lw=0.9, label=f"±{THRESHOLD_M_YR} m/yr"),
                        Line2D([], [], color=INK_MUTED, ls="-.", lw=0.8, label=f"±{STRICT_M_YR:g} m/yr")],
               loc="outside lower center", ncol=4, fontsize=7, frameon=False)
    save(fig, PNG, dpi=300, close=True)
    record_caption(PNG, (
        "Step 2 of the source/sink zones: candidate zones from the step 1 residuals (edgeBE "
        "full-management run, no groin; observed minus modelled net change over the run length, "
        "both at 7-domain LOWESS). A candidate is a run of at least three adjacent domains beyond "
        f"±{THRESHOLD_M_YR} m/yr (dotted) with one sign; blue a source, red a sink, labelled as in "
        "zone_candidates_edgeBE.csv. (a) Scenario A, shared: beyond the threshold in both periods "
        "with the same sign. (b, c) Scenario B, each period on its own. Hatched: GIS 1 and 90 "
        "(solved end rates) and GIS 5-7 (the groin), never candidates. ±1 m/yr (dash-dot) is "
        "reported per zone but does not select. No physical screen has been applied yet."))
    print(f"wrote {PNG}")


def main():
    t = pd.read_csv(step1.TABLE)
    r = t.pivot(index="gis_domain", columns="period", values="residual_m_yr")
    z = candidates(r)
    TABLE.parent.mkdir(parents=True, exist_ok=True)
    z.to_csv(TABLE, index=False)
    pd.set_option("display.width", 200)
    print(z.to_string(index=False))
    print(f"wrote {TABLE}")
    figure(r, z)


if __name__ == "__main__":
    main()
