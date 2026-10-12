#!/usr/bin/env python3
"""
The matrix runs of each window and preset in order of increasing management, one panel per rung.

    python scripts/analyze_output/compare_runs/matrix_vs_observed/matrix_management_ladder.py

Two versions per window and preset: rate (LRR vs the CoastSat target) and net
change (final minus initial state vs the DEM-to-DEM net-change target, since
2026-10-12; before it, position change vs the calendar-year total change).
A preset with only one rung is skipped. Writes output/raw_runs/matrix/figures/.
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-12
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
from matplotlib.transforms import blended_transform_factory  # noqa: E402

_REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(_REPO / "scripts"))
sys.path.insert(0, str(_REPO / "scripts" / "hatteras_ms"))
sys.path.insert(0, str(_REPO / "scripts" / "hatteras_ms" / "experiments"))

import HAT_metres_1_offset_units as common  # noqa: E402
from matrix_vs_observed import (end_year, matrix_runs, model_end_change, rates,  # noqa: E402
                                score, window)
from site_layer.hat_observed_rates import NET_CHANGE_WINDOWS, net_change_domain_csv  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK, INK_MUTED, _title, apply_style, figsize, mark_offaxis,
    open_frame, record_caption, save, structures, support_dir, town_bands)

# --- CONFIG ------------------------------------------------------------------
OUT = _REPO / "output" / "raw_runs" / "matrix" / "figures"
PRESETS = ("edgeBE", "zeroBE")
OBSERVED = dict(color="black", lw=2.0)
# One colour per layer, the same in every window; a rung takes its last layer's
# Okabe-Ito, colour-blind safe and well clear of the black observation (a blue ramp, until
# 2026-10-12, ran its darkest rung into the black)
LAYER_COLOUR = {"none": "#E69F00", "road": "#56B4E9", "bdm": "#009E73", "fills": "#D55E00"}
LW_NEWEST, LW_EARLIER = 1.7, 1.0

PRESET_TEXT = {"edgeBE": "end rates solved on CoastSat (edgeBE)",
               "zeroBE": "no imposed end rates (zeroBE)"}
LAYER_TEXT = {"road": "road management", "bdm": "beach and dune management",
              "fills": "nourishment fills", "reloc": "historical relocations"}

VERSIONS = {
    "rate": dict(
        ylim=(-7.5, 7.5), ylabel="Shoreline change rate,\nLRR (m/yr)", unit="m/yr", fmt="{:+.2f}",
        rmse_fmt="{:.2f}", obs_label="CoastSat LRR target",
        model=lambda r: rates(r.run_dir)["lrr_m_yr"],
        obs_text="the CoastSat LRR scoring target (black; 7-domain LOWESS, raw domain means GIS 1–10)",
        what="modelled OLS shoreline-change rate"),
    "net_change": dict(
        ylim=(-100, 100), ylabel="Net shoreline change (m)", unit="m", fmt="{:+.1f}", rmse_fmt="{:.1f}",
        obs_label="CoastSat net-change target",
        model=lambda r: model_end_change(r.run_dir),
        obs_text=None,           # set per window from NET_CHANGE_WINDOWS
        what="modelled net change (final minus initial state)"),
}
# -----------------------------------------------------------------------------


# The shoal zones (site config), hatched gold as in add_geographic_annotations, named at the foot;
# drawn after town_bands, over the grey and under the data
def shoals(ax, label=True):
    from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS as A
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    for name, (lo, hi) in A.shoal_zones.items():
        ax.axvspan(lo - 0.5, hi + 0.5, color=A.color_shoal, alpha=0.10, zorder=1, hatch="///",
                   edgecolor=A.color_shoal, linewidth=0)
        if label:
            ax.text((lo + hi) / 2.0, 0.04, name, transform=trans, ha="center", va="bottom",
                    fontsize=7, color="#7A5800", style="italic", zorder=7,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))


# The management layers a matrix run has, read from the run itself
def layers(r):
    md = json.loads(next(Path(r.run_dir).glob("*_run_metadata.json")).read_text(encoding="utf-8"))
    fills = re.match(r"\s*(\d+)", str(md["scenario"].get("nourishment fills", "0")))
    return {"road": r.scenario in ("roadway_only", "full_no_fill", "full_management"),
            "bdm": r.scenario in ("beachdune_only", "full_no_fill", "full_management"),
            "fills": r.scenario == "full_management" and fills is not None and int(fills.group(1)) > 0,
            "reloc": bool(r.relocations_enabled)}


# The rungs for one window and preset, least managed first; a run adding nothing is dropped
def ladder(runs):
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


# Panel title: what this rung adds
def rung_label(i, added):
    return "Natural" if i == 0 else "+ " + " and ".join(LAYER_TEXT[k] for k in added)


# One figure: panel k draws rungs 1..k over the observation; records each rung's scores
def draw(version, period, preset, rungs, observed, rows_out):
    v = VERSIONS[version]
    n = len(rungs)
    fig, axes = plt.subplots(n, 1, sharex=True, constrained_layout=True,
                             figsize=figsize("double", height=min(1.75 * n + 0.8, 9.4)))
    # Each rung smoothed as the observation is
    series = [common.smooth_like_target(v["model"](r)) for r, _ in rungs]
    labels = [rung_label(i, added) for i, (_, added) in enumerate(rungs)]
    colours = [LAYER_COLOUR[added[-1] if added else "none"] for _, added in rungs]
    for k, ax in enumerate(axes):
        ax.plot(observed.index, observed.values, zorder=10, **OBSERVED)
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
        mark_offaxis(ax, observed.index.to_numpy(float), observed.to_numpy(float),
                     v["ylim"][1], color=OBSERVED["color"])
        for i in range(k + 1):     # a held end domain beyond the axis, marked at the edge
            mark_offaxis(ax, series[i].index.to_numpy(float), series[i].to_numpy(float),
                         v["ylim"][1], color=colours[i])
        ax.grid(axis="y")
        open_frame(ax)
        town_bands(ax, label=(k == 0))
        shoals(ax, label=(k == 0))
        structures(ax, label=(k == n - 1))
        ax.set_ylabel(v["ylabel"])
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
    obs_text = v["obs_text"] or net_change_text(period)
    record_caption(png, (
        f"The {w.replace('_', '–')} matrix runs, {PRESET_TEXT[preset]}, in order of "
        f"increasing management. Each panel adds one layer to the panel above and draws every run "
        f"up to that point: the {v['what']} of each rung in its own colour (Okabe-Ito; "
        f"the newest rung heaviest, and each rung keeps its colour down the figure), under "
        f"{obs_text}. RMSE and bias of the newest rung over the interior GIS 2–89 (GIS 1 at "
        f"Cape Point, 90 at Pea Island); seaward positive. Every model line is smoothed as the "
        f"observation is (7-domain LOWESS, GIS 1–10 raw). A rung that adds nothing in this window "
        f"is skipped. Hatched gold: the Avon (GIS 24–39) and Wimble (GIS 60–74) shoal zones. A value beyond the axis (GIS 1, observed and modelled) is marked with a triangle at the edge. The historical-relocation runs are not on the ladder "
        f"(relocation moves the road, not the shoreline). Runs: {runs}."))
    return png


# The net-change target of a window, LOWESS-7 (GIS 1-10 raw), seaward +
def net_change(period):
    d = pd.read_csv(net_change_domain_csv(period, end_year(period)))
    return d.set_index(d.columns[0])["net_change_lowess7_m"]


def net_change_text(period):
    (s0, s1), (e0, e1) = NET_CHANGE_WINDOWS[(period, end_year(period))]
    return (f"the CoastSat net-change target (black; the mean shoreline over {e0} to {e1} minus "
            f"the mean over {s0} to {s1}, the island offset's window; 7-domain LOWESS, GIS 1–10 raw)")


# Every rung overlaid in one panel per window, net change only; the windows stacked
def draw_overlay(preset, ladders):
    v = VERSIONS["net_change"]
    n = len(ladders)
    fig, axes = plt.subplots(n, 1, sharex=True, constrained_layout=True,
                             figsize=figsize("double", height=3.6 * n + 0.6), squeeze=False)
    parts = []
    for k, (ax, (period, (rungs, observed))) in enumerate(zip(axes[:, 0], ladders.items())):
        ax.plot(observed.index, observed.values, zorder=10, **OBSERVED)
        handles = [Line2D([], [], label=v["obs_label"], **OBSERVED)]
        for i, (r, added) in enumerate(rungs):
            s = common.smooth_like_target(v["model"](r))
            c = LAYER_COLOUR[added[-1] if added else "none"]
            last = i == len(rungs) - 1
            ax.plot(s.index, s.values, color=c, lw=LW_NEWEST if last else 1.3, zorder=3 + i)
            mark_offaxis(ax, s.index.to_numpy(float), s.to_numpy(float), v["ylim"][1], color=c)
            bias, rmse = score(s, observed)
            handles.append(Line2D([], [], color=c, lw=LW_NEWEST if last else 1.3,
                                  label=f"{rung_label(i, added)}: RMSE {rmse:.1f}, bias {bias:+.1f} m"))
        mark_offaxis(ax, observed.index.to_numpy(float), observed.to_numpy(float), v["ylim"][1],
                     color=OBSERVED["color"])
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_xlim(1, 90)
        ax.set_ylim(*v["ylim"])
        ax.grid(axis="y")
        open_frame(ax)
        town_bands(ax, label=True)
        shoals(ax)
        structures(ax, label=(k == n - 1))
        ax.set_ylabel(v["ylabel"])
        w = window(period)
        _title(ax, k, f"{'Calibration' if period == 1996 else 'Test'} {w.replace('_', '–')}")
        ax.legend(handles=handles, loc="upper center", ncol=2, fontsize=6.5, frameon=False,
                  bbox_to_anchor=(0.5, -0.04 if k < n - 1 else -0.13))
        parts.append(f"({chr(97 + k)}) {w.replace('_', '–')}: {net_change_text(period)}")
    axes[-1, 0].set_xlabel(DOMAIN_AXIS_LABEL)
    png = OUT / f"management_ladder_overlay_net_change_{preset}.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        f"The management ladder overlaid, {PRESET_TEXT[preset]}: every run of a window in one "
        f"panel, natural to full management, in the ladder figures' colours (orange natural, sky blue "
        f"+ roads, green + beach and dune management, vermillion + fills; the last rung heaviest) "
        f"under the black target. The {v['what']}, smoothed as the target is. RMSE and bias "
        f"over the interior GIS 2–89 in the legend; seaward positive. Each rung uses the window's "
        f"full-management end rates. Targets: " + "; ".join(parts) + ". Hatched gold: the Avon (GIS 24–39) and Wimble (GIS 60–74) "
        f"shoal zones. A value beyond the axis "
        f"is marked with a triangle at the edge."))
    return png


# Run: clear the retired step-panel figures, draw every ladder, write the scores
def main():
    apply_style()
    runs = matrix_runs()
    for old in OUT.glob("*/management_ladder_*.png"):   # the step-panel version, 09-29
        if not old.stem.startswith(("management_ladder_rate_", "management_ladder_position_",
                                       "management_ladder_net_change_",
                                       "management_ladder_overlay_")):
            old.unlink()
            pdf = old.parent / "supporting" / (old.stem + ".pdf")
            if pdf.exists():
                pdf.unlink()
    # Both versions for every window and preset
    rows, overlay = [], {}
    for period in sorted(int(p) for p in runs.period.unique()):
        obs = {"rate": common.coastsat_target(period), "net_change": net_change(period)}
        for preset in PRESETS:
            sub = runs[(runs.period == period) & (runs.source_sink_preset == preset)]
            if sub.empty:
                continue
            rungs = ladder(sub)
            if len(rungs) < 2:             # nothing to build up
                continue
            overlay.setdefault(preset, {})[period] = (rungs, obs["net_change"])
            for version in VERSIONS:
                png = draw(version, period, preset, rungs, obs[version], rows)
                print(f"wrote {png.relative_to(_REPO)}  ({len(rungs)} rungs)")
    for preset, ladders in overlay.items():
        print(f"wrote {draw_overlay(preset, ladders).relative_to(_REPO)}")
    # Every rung's scores in one table
    table = support_dir(OUT) / "ladder_scores.csv"
    table.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(table, index=False)
    print(f"wrote {table.relative_to(_REPO)}")


if __name__ == "__main__":
    main()
