#!/usr/bin/env python3
"""
matrix_management_ladder.py
==============================================================================
The matrix runs of one window and preset in order of increasing management,
so the effect of each setting can be seen building up (Hannah, 2026-09-29:
"all of the runs in sequential order of increasing management so we can see
the effect of each setting"; then "remove the step panels ... show all of the
curves in each progressive step, also do a lrr version and position change
version").

THE LADDER
    natural -> road management only -> + beach and dune management
            -> + nourishment fills
    Each rung is the matrix run that adds one layer to the rung above. A rung
    the window does not have is skipped rather than drawn twice: 1996-2010 has
    no fills (the driver skips full_no_fill there as identical to
    full_management). What each rung adds is read from the runs themselves
    (scenario and the fills the run applied), not typed. Left off the ladder:
    beach/dune-only, a side branch, and the historical-relocation runs
    (Hannah, 2026-09-29: "remove the historical relocations panel"; in
    1996-2010 it matched full management to 0.0002 m/yr, since relocation
    moves the road, not the shoreline; 2010-2024 has no relocation event).

ONE PANEL PER RUNG, CUMULATIVE
    Panel k draws rungs 1..k, each in its own shade of the house blue ramp
    (lighter = less managed), the newest rung heaviest, over the observation.
    A rung keeps its colour in every panel, so a curve can be followed down.

TWO VERSIONS
    rate      the OLS rate (lrr_m_yr) against the CoastSat LRR scoring target
              (7-domain LOWESS, raw means GIS 1-10), as the runner scores it;
              +/-7.5 m/yr, as the runner's own figure
    position  the position change over the window (endpoint rate x 14 yr)
              against the observed CoastSat change (mean position in the last
              calendar year minus the first, smoothed at 7 domains);
              +/-130 m, as output/comparisons/matrix_vs_observed/
    Interior (GIS 2-89) RMSE and bias of the newest rung in each row title,
    every rung's in supporting/ladder_scores.csv.

WRITES
    output/raw_runs/matrix/figures/<window>/management_ladder_{rate,position}_<preset>_<window>.png
    output/raw_runs/matrix/figures/<window>/supporting/  PDF, CAPTIONS.md
    output/raw_runs/matrix/figures/supporting/ladder_scores.csv

USAGE
    python scripts/analyze_output/compare_runs/matrix_management_ladder.py
==============================================================================

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

_REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(_REPO / "scripts"))
sys.path.insert(0, str(_REPO / "scripts" / "hatteras_ms"))
sys.path.insert(0, str(_REPO / "scripts" / "hatteras_ms" / "experiments"))

import HAT_metres_1_offset_units as common  # noqa: E402
from matrix_vs_observed import (YEARS, matrix_runs, observed_change, rates,  # noqa: E402
                                score, window)
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK, INK_MUTED, SMOOTH_RAMP, _title, apply_style, figsize,
    open_frame, record_caption, save, structures, support_dir, town_bands)

OUT = _REPO / "output" / "raw_runs" / "matrix" / "figures"
PRESETS = ("edgeBE", "zeroBE")
OBSERVED = dict(color=INK, lw=1.6)
# One colour per layer, the same in every window (Hannah, 2026-09-29: "keep
# colours consistent"): a rung takes the colour of the last layer it adds,
# light -> dark with management, so 1996-2010 (no fills) simply never uses the
# darkest.
LAYER_COLOUR = dict(zip(("none", "road", "bdm", "fills"), SMOOTH_RAMP))
LW_NEWEST, LW_EARLIER = 1.7, 1.0

PRESET_TEXT = {"edgeBE": "end rates solved on CoastSat (edgeBE)",
               "zeroBE": "no imposed end rates (zeroBE)"}
LAYER_TEXT = {"road": "road management", "bdm": "beach and dune management",
              "fills": "nourishment fills", "reloc": "historical relocations"}

VERSIONS = {
    "rate": dict(
        ylim=(-7.5, 7.5), ylabel="Shoreline change rate,\nLRR (m/yr)", unit="m/yr", fmt="{:+.2f}",
        rmse_fmt="{:.2f}", obs_label="CoastSat LRR target",
        model=lambda rt: rt["lrr_m_yr"],
        obs_text="the CoastSat LRR scoring target (black; 7-domain LOWESS, raw domain means GIS 1–10)",
        what="modelled OLS shoreline-change rate"),
    "position": dict(
        ylim=(-130, 130), ylabel=None, unit="m", fmt="{:+.1f}", rmse_fmt="{:.1f}",
        obs_label="CoastSat observed change",
        model=lambda rt: rt["change_rate_m_yr"] * YEARS,
        obs_text=("the observed CoastSat change (black; mean position over the last calendar year "
                  "minus the first, smoothed at 7 domains)"),
        what="modelled position change over the window (endpoint rate × 14 yr)"),
}


def layers(r):
    """The management layers a matrix run has, from the run itself."""
    md = json.loads(next(Path(r.run_dir).glob("*_run_metadata.json")).read_text(encoding="utf-8"))
    fills = re.match(r"\s*(\d+)", str(md["scenario"].get("nourishment fills", "0")))
    return {"road": r.scenario in ("roadway_only", "full_no_fill", "full_management"),
            "bdm": r.scenario in ("beachdune_only", "full_no_fill", "full_management"),
            "fills": r.scenario == "full_management" and fills is not None and int(fills.group(1)) > 0,
            "reloc": bool(r.relocations_enabled)}


def ladder(runs):
    """The rungs for one window and preset, least managed first, each adding
    one or more layers to the one above; a run that adds nothing is dropped."""
    order = [("natural", False), ("roadway_only", False), ("full_no_fill", False),
             ("full_management", False)]
    rungs, have = [], None
    for scen, reloc in order:
        hit = runs[(runs.scenario == scen) & (runs.relocations_enabled == reloc)]
        if hit.empty:
            continue
        r = hit.iloc[0]
        lay = layers(r)
        if have is not None and lay == have:
            continue
        added = [k for k in LAYER_TEXT if lay[k] and not (have or {}).get(k)]
        rungs.append((r, added))
        have = lay
    return rungs


def rung_label(i, added):
    return "Natural" if i == 0 else "+ " + " and ".join(LAYER_TEXT[k] for k in added)


def draw(version, period, preset, rungs, observed, rows_out):
    v = VERSIONS[version]
    n = len(rungs)
    fig, axes = plt.subplots(n, 1, sharex=True, constrained_layout=True,
                             figsize=figsize("double", height=min(1.75 * n + 0.8, 9.4)))
    series = [v["model"](rates(r.run_dir)) for r, _ in rungs]
    labels = [rung_label(i, added) for i, (_, added) in enumerate(rungs)]
    colours = [LAYER_COLOUR[added[-1] if added else "none"] for _, added in rungs]
    for k, ax in enumerate(axes):
        ax.plot(observed.index, observed.values, zorder=2, **OBSERVED)
        for i in range(k + 1):
            newest = i == k
            ax.plot(series[i].index, series[i].values, color=colours[i],
                    lw=LW_NEWEST if newest else LW_EARLIER, zorder=4 if newest else 3)
        bias, rmse = score(series[k], observed)
        _title(ax, k, labels[k] if k == 0 else f"{labels[k]}")
        ax.set_title(f"RMSE {v['rmse_fmt'].format(rmse)}, bias {v['fmt'].format(bias)} {v['unit']}",
                     loc="right", color=INK_MUTED, fontsize=7)
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_xlim(1, 90)
        ax.set_ylim(*v["ylim"])
        ax.grid(axis="y")
        open_frame(ax)
        town_bands(ax, label=(k == 0))
        structures(ax, label=(k == n - 1))
        ax.set_ylabel(v["ylabel"] or f"Position change,\n{period + YEARS} minus {period} (m)")
        r = rungs[k][0]
        row = next((x for x in rows_out if x["run_name"] == r.run_name
                    and x["preset"] == preset), None)
        if row is None:
            row = dict(window=window(period), preset=preset, rung=k + 1, run_name=r.run_name,
                       adds=", ".join(rungs[k][1]) or "none")
            rows_out.append(row)
        row[f"{version}_bias"] = round(bias, 4)
        row[f"{version}_rmse"] = round(rmse, 4)
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    fig.legend(handles=[Line2D([], [], label=v["obs_label"], **OBSERVED)]
               + [Line2D([], [], color=c, lw=LW_NEWEST, label=lab) for c, lab in zip(colours, labels)],
               loc="outside lower center", ncol=min(n + 1, 3), frameon=False)
    w = window(period)
    png = OUT / w / f"management_ladder_{version}_{preset}_{w}.png"
    save(fig, png, dpi=300, close=True)
    runs = "; ".join(f"({chr(97 + i)}) {r.run_name}" for i, (r, _) in enumerate(rungs))
    record_caption(png, (
        f"The {w.replace('_', '–')} option A matrix runs, {PRESET_TEXT[preset]}, in order of "
        f"increasing management. Each panel adds one layer to the panel above and draws every run "
        f"up to that point: the {v['what']} of each rung in its own blue (lighter = less managed; "
        f"the newest rung heaviest, and each rung keeps its colour down the figure), over "
        f"{v['obs_text']}. RMSE and bias of the newest rung over the interior GIS 2–89 (GIS 1 at "
        f"Cape Point, 90 at Pea Island); seaward positive. Rungs the window lacks are skipped: "
        f"1996–2010 has no nourishment fills. The historical-relocation runs are not on the ladder "
        f"(relocation moves the road, not the shoreline). Runs: {runs}."))
    return png


def main():
    apply_style()
    runs = matrix_runs()
    for old in OUT.glob("*/management_ladder_*.png"):   # the step-panel version, 09-29
        if not old.stem.startswith(("management_ladder_rate_", "management_ladder_position_")):
            old.unlink()
            pdf = old.parent / "supporting" / (old.stem + ".pdf")
            if pdf.exists():
                pdf.unlink()
    rows = []
    for period in sorted(int(p) for p in runs.period.unique()):
        obs = {"rate": common.coastsat_target(period), "position": observed_change(period)}
        for preset in PRESETS:
            sub = runs[(runs.period == period) & (runs.source_sink_preset == preset)]
            if sub.empty:
                continue
            rungs = ladder(sub)
            for version in VERSIONS:
                png = draw(version, period, preset, rungs, obs[version], rows)
                print(f"wrote {png.relative_to(_REPO)}  ({len(rungs)} rungs)")
    table = support_dir(OUT) / "ladder_scores.csv"
    table.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(table, index=False)
    print(f"wrote {table.relative_to(_REPO)}")


if __name__ == "__main__":
    main()
