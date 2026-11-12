"""
Shoreline position change alongshore with the dune collision module on and off, calibration and test windows; and
the 2009-2025 wave-height sweep with the module on.

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_on_off_position.py [onoff|hs|sand|upkeep|village|toes|toesnow]

Drawn as the run figures draw position change (end minus start, seaward positive, the model smoothed with the
7-domain LOWESS, CoastSat net change as the observation), from runs already made (nothing rerun here).

onoff  one figure per window: the module-off control and the module on as it now stands (dry beach, W0 25 m,
       Cs 1.02e-4, first-dune toes, off at Cape Point GIS 1-5, the 198 m road-dune rule):
           1996-2009  cal_control_g16_rd198 and cal_g16_nocape_rd198_dry (1-calibration-1996-2009/HAT_dune_collision.py)
           2009-2025  control_rd198 and nocape_rd198_dry (2-test-2009-2025/HAT_collision_test_2009.py)
       Both module-off runs have the 198 m road rule too, so the difference is the module alone.
hs     2009-2025, the module on at Hs 1.5-3.5 (nocape, nocape_hs*; HAT_collision_test_2009.py). These ran on the
       settings before 2026-11-03 (width to the crest, W0 45 m, Cs 1.24e-4, deepest-dip toes, no road rule), and
       only Hs 3.0 has solved ends: the others are on ENDS_BY_HS, interpolated or extrapolated.

Writes, under the workspace's dune-collision/4-hindcast/ (each figure's topic folder): collision_on_off_position_<start>_<end>.png and
collision_hs_sweep_position_2009_2025.png, with the captions under supporting/.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-03
"""
from __future__ import annotations

import argparse
import dataclasses
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
for _p in (_HAT_DIR, _HERE.parents[1] / "1-calibration-1996-2009", _HERE.parents[1] / "2-test-2009-2025"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
from hatteras_paths import WORKSPACE  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = WORKSPACE / "dune-collision" / "4-hindcast"
# (start, end, script module, off variant, on variant)
WINDOWS = ((1996, 2009, "HAT_dune_collision", "cal_control_g16_rd198", "cal_g16_nocape_rd198_dry"),
           (2009, 2025, "HAT_collision_test_2009", "control_rd198", "nocape_rd198_dry"))
# The off runs that have the 198 m road rule too, so on minus off is the module alone (control_rd198, 2026-11-03)
OFF_HAS_ROAD_RULE = ("control_rd198", "cal_control_g16_rd198")
OFF_STYLE = ("Module off", "#7f7f7f", 1.6)
ON_STYLE = ("Module on (dry beach, 25 m, Cs 1.02e-4; 198 m road rule)", "#7b3294", 2.0)
# The Hs sweep, module on (variant, label, colour, width); the first is the one the run figure's call draws
HS_SWEEP = (("nocape", "Hs 3.0 m (solved ends)", "#54278F", 2.2),
            ("nocape_hs1p5", "Hs 1.5 m", "#E7D4E8", 1.5),
            ("nocape_hs2", "Hs 2.0 m", "#C2A5CF", 1.5),
            ("nocape_hs2p5", "Hs 2.5 m", "#9970AB", 1.6),
            ("nocape_hs3p5", "Hs 3.5 m", "#E08214", 1.8))
Y_HALF_M = 100.0
# -----------------------------------------------------------------------------


# One position-change figure: runs = ((variant, label, colour, width), ...), the first drawn by the run figure's call
def figure(start, end, module, runs, title, out_name, caption, height_scale=1.0, legend_runs=None):
    import importlib
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from cascade_pipeline.coastsat_lowess import CoastSatDataset, build_coastsat_series
    from cascade_pipeline.plotting.comparison_style import restyle
    from cascade_pipeline.plotting.net_change_comparison import LOWESS_CONFIG, smooth_model
    from cascade_pipeline.plotting.rate_comparison import DEFAULT_RATE_COMPARISON, plot_annotated_rate_comparison
    from cascade_pipeline.run_info import RunInfo
    from cascade_pipeline.run_layout import resolve
    from cascade_pipeline.shoreline import compute_change_rate
    from site_layer.hat_figure_style import mark_offaxis, offaxis_clause, record_caption
    from site_layer.hat_observed_rates import NET_CHANGE_WINDOWS, WINDOW_ROLE, net_change_dir
    from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS

    m = importlib.import_module(module)
    changes = {}
    for v, *_ in runs:
        rd = m.run_dir(v)
        if rd is None:
            raise SystemExit(f"missing run {v} ({module})")
        changes[v] = (rd, smooth_model(compute_change_rate(np.load(resolve(rd, "matrix", rd.name)), span_years=1,
                                                           flip_sign=True)))
    tag = f"{start}_{end}"
    cs = build_coastsat_series(
        [CoastSatDataset(label=f"CoastSat net change ({start}–{end})", period_start=start,
                         csv_path=str(net_change_dir(start, end) / f"transect_net_change_{tag}.csv"),
                         rate_col="net_change_m")],
        active_period_start=start, lowess_config=LOWESS_CONFIG, domains=HATTERAS_DOMAINS)
    # each run: (variant, label, colour, line width[, line style])
    runs = [tuple(r) + (("-",) if len(r) == 4 else ()) for r in runs]
    base, base_label, base_color, base_lw, base_ls = runs[0]
    base_dir, base_change = changes[base]
    gis = np.arange(HATTERAS_DOMAINS.first_gis_id, HATTERAS_DOMAINS.last_gis_id + 1)
    g = HATTERAS_DOMAINS.for_total(len(base_change))
    real = slice(g.start_real_index, g.end_real_index)
    subtitle = " · ".join((WINDOW_ROLE.get((start, end), ""), "whole island", "full management + groin, edgeBE"))
    config = dataclasses.replace(
        DEFAULT_RATE_COMPARISON, quantity="position", publication_text=True, show_features=True,
        plot_domain_means=False, raw_lrr_southern_only=True, smoothing_marker=False,
        observed_label="CoastSat (7-domain LOWESS)", model_label=base_label, title=title,
        ylim=(-Y_HALF_M, Y_HALF_M), subtitle=subtitle)
    run = RunInfo(run_name=base_dir.name, run_dir=str(base_dir), start_year=start, end_year=end, Hs=3.0,
                  flip_sign_model=True)
    fig, ax = plot_annotated_rate_comparison(
        base_change, cs, run, estimator=None, domains=HATTERAS_DOMAINS,
        annotations=dataclasses.replace(HATTERAS_ANNOTATIONS, model_color=base_color),
        lowess_config=LOWESS_CONFIG, config=config)
    for line in ax.get_lines():
        if line.get_label() == base_label:
            line.set_linewidth(base_lw)
            line.set_linestyle(base_ls)
            line.set_zorder(5)
    named = [(base_label, mark_offaxis(ax, gis, base_change[real], Y_HALF_M, color=base_color))]
    for v, label, color, lw, ls in runs[1:]:
        ch = changes[v][1]
        ax.plot(gis, ch[real], color=color, lw=lw, ls=ls, zorder=6)
        named.append((label, mark_offaxis(ax, gis, ch[real], Y_HALF_M, color=color)))
    old = fig.legends.pop()
    entries = [(h, t.get_text()) for h, t in zip(old.legend_handles, old.get_texts()) if t.get_text() != base_label]
    # legend_runs: every run of a slide build-up; the ones not drawn yet keep their space but are invisible, so every
    # level of the build-up has the same legend and figure size
    shown = {r[0] for r in runs}
    slots = [tuple(r) + (("-",) if len(r) == 4 else ()) for r in (legend_runs or runs)]
    handles = [Line2D([0], [0], color=c if v in shown else (0, 0, 0, 0), lw=w + 0.4, ls=ls, label=lab)
               for v, lab, c, w, ls in slots]
    fig.legend(handles + [h for h, _ in entries], [h.get_label() for h in handles] + [t for _, t in entries],
               loc="outside lower center")
    restyle(fig, ax, subtitle)
    hidden = {lab for v, lab, *_ in slots if v not in shown}
    for leg in fig.legends:                       # after restyle, which redraws the legend text
        for text in leg.get_texts():
            if text.get_text() in hidden:
                text.set_alpha(0)
    if height_scale != 1.0:                       # taller plot, same axis range: closer lines pull apart
        w, h = fig.get_size_inches()
        fig.set_size_inches(w, h * height_scale)
    out = OUT / out_name                          # out_name may start with a topic subfolder
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close("all")
    (s0, s1), (e0, e1) = NET_CHANGE_WINDOWS[(start, end)]
    record_caption(out, (
        f"Shoreline position change, {start}–{end}, whole island; positive is seaward. Modelled (CASCADE): the "
        f"shoreline at the end of each run minus the start, smoothed with a 7-domain LOWESS north of domain 10, "
        f"for the full-management edgeBE run with the groin, {caption} Observed (CoastSat, blue): per transect, the "
        f"mean position over {e0} to {e1} minus the mean over {s0} to {s1}, smoothed with a 7-domain LOWESS north "
        f"of domain 10; dots south of it. The y axis is fixed at ±{Y_HALF_M:g} m."
        + offaxis_clause([(k, n) for k, n in named if n], Y_HALF_M)))
    return out, gis, {v: c[1][real] for v, c in changes.items()}


def onoff():
    for start, end, module, off, on in WINDOWS:
        out, gis, ch = figure(
            start, end, module, ((off, OFF_STYLE[0] + (" (198 m road rule)" if off in OFF_HAS_ROAD_RULE else ""), *OFF_STYLE[1:]),
                                  (on, *ON_STYLE)),
            f"Shoreline position change, dune collision module on and off, {start}–{end}",
            f"1-module-on-off/collision_on_off_position_{start}_{end}.png",
            f"the dune collision module off ({off}, grey) and on ({on}, purple: Larson et al. 2004 per storm times "
            f"max(0, 1 - W/25 m) on the dry beach, Cs 1.02e-4, first-dune toes, off at Cape Point GIS 1-5, with "
            f"the 198 m road-dune rule, which the module-off run "
            + ("has too)." if off in OFF_HAS_ROAD_RULE else "does not have)."))
        d = ch[on] - ch[off]
        print(f"wrote {out}\n  on minus off: mean {d.mean():+.1f} m, range {d.min():+.1f} (GIS {gis[d.argmin()]}) to "
              f"{d.max():+.1f} m (GIS {gis[d.argmax()]})")


def hs():
    out, gis, ch = figure(
        2009, 2025, "HAT_collision_test_2009", HS_SWEEP,
        "Shoreline position change with the collision module on, wave height sweep, 2009–2025",
        "8-other-forcings/collision_hs_sweep_position_2009_2025.png",
        "the dune collision module on (off at Cape Point GIS 1-5) at deep-water wave heights (alongshore "
        "diffusivity) of 1.5-3.5 m. These runs are on the module's settings before 2026-11-03 (beach width to the "
        "crest, W0 45 m, Cs 1.24e-4, deepest-dip toes, no 198 m road rule). Only Hs 3.0 runs on end conditions "
        "solved for it; the others are interpolated (2.5) or extrapolated (1.5, 3.5) from the Hs 2.0 and 3.0 solves.")
    print(f"wrote {out}")
    for v, label, *_ in HS_SWEEP:
        print(f"  {label:24s} GIS 24 {ch[v][23]:+6.1f}  GIS 30 {ch[v][29]:+6.1f}  GIS 80 {ch[v][79]:+6.1f} m")


def sand():
    runs = (("control_rd198", "Module off (198 m road rule)", *OFF_STYLE[1:]),
            ("nocape_rd198_dry", "Module on, eroded sand deleted (as now)", *ON_STYLE[1:]),
            ("nocape_rd198_dry_sand", "Module on, eroded sand onto the beach", "#E08214", 2.0))
    out, gis, ch = figure(
        2009, 2025, "HAT_collision_test_2009", runs,
        "Shoreline position change, collision-eroded sand deleted or put on the beach, 2009–2025",
        "8-other-forcings/collision_sand_to_beach_position_2009_2025.png",
        "all with the 198 m road-dune rule: the dune collision module off (grey), on as it is, the eroded dune "
        "sand deleted (purple), and on with that sand placed on the beach and shoreface after each storm as CASCADE "
        "places a nourishment (orange; shoreline seaward by 2V/(2h_b + D_shoreface), Ashton and Lorenzo-Trueba 2018).")
    print(f"wrote {out}")
    for v in ("nocape_rd198_dry", "nocape_rd198_dry_sand"):
        d = ch[v] - ch["control_rd198"]
        print(f"  {v:24s} minus off: mean {d.mean():+.1f} m, range {d.min():+.1f} to {d.max():+.1f} m")


# Rodanthe village dune upkeep at GIS 77-83 (2026-11-03), the module off, all with the 198 m road rule
def upkeep():
    runs = (("control_rd198", "Module off (198 m road rule)", *OFF_STYLE[1:]),
            ("control_rd198_ow30", "Rodanthe: 30% of overwash pushed onto the dune", "#74C476", 1.8),
            ("control_rd198_ow59", "Rodanthe: 59% of overwash pushed onto the dune", "#006D2C", 2.0),
            ("control_rd198_upkeep", "Rodanthe: dunes kept to 2009 height, sand off the beach", "#E08214", 2.0))
    out, gis, ch = figure(
        2009, 2025, "HAT_collision_test_2009", runs,
        "Shoreline position change with Rodanthe village dune upkeep, 2009–2025",
        "8-other-forcings/rodanthe_dune_upkeep_position_2009_2025.png",
        "the dune collision module off and the 198 m road-dune rule in all; only Rodanthe village (GIS 77-83) "
        "differs. Grey: the village default (9% of each year's overwash bulldozed onto the dune, 40% returned to "
        "the shoreface). Greens: 30% and 59% bulldozed onto the dune. Orange: a dune cell more than 0.5 m below its "
        "2009 starting height rebuilt to it after each year's storms, the sand taken off the beach (shoreline "
        "landward by 2V/(2h_b + D_shoreface)).")
    print(f"wrote {out}")
    for v in ("control_rd198_ow30", "control_rd198_ow59", "control_rd198_upkeep"):
        d = np.asarray(ch[v] - ch["control_rd198"])          # GIS 1-90 in order
        rod = (gis >= 77) & (gis <= 83)
        print(f"  {v:22s} minus off at GIS 77-83: mean {d[rod].mean():+.1f} m; elsewhere max |d| "
              f"{np.abs(d[~rod]).max():.1f} m")


# The village overwash share returned to the shoreface, 40% (default) and 90%, module off and on (2026-11-04)
def village():
    runs = (("control_rd198", "Module off, villages 40% (default)", *OFF_STYLE[1:]),
            ("control_rd198_vf90", "Module off, villages 90%", "#BDBDBD", 1.8),
            ("nocape_rd198_dry", "Module on, villages 40% (default)", *ON_STYLE[1:]),
            ("nocape_rd198_dry_vf90", "Module on, villages 90%", "#C2A5CF", 1.8))
    out, gis, ch = figure(
        2009, 2025, "HAT_collision_test_2009", runs,
        "Shoreline position change with more village overwash returned to the beach, 2009–2025",
        "8-other-forcings/village_overwash_position_2009_2025.png",
        "the 198 m road-dune rule in all. In the managed villages each year's overwash is cleared and a share "
        "returned to the shoreface: 40% (the default) or 90% (Rogers et al. 2015 give 40-90%, residential to "
        "commercial). Greys: the dune collision module off; purples: on (dry beach, W0 25 m, Cs 1.02e-4, "
        "first-dune toes, off at Cape Point GIS 1-5).")
    print(f"wrote {out}")
    for a, b in (("control_rd198_vf90", "control_rd198"), ("nocape_rd198_dry_vf90", "nocape_rd198_dry"),
                 ("nocape_rd198_dry", "control_rd198")):
        d = np.asarray(ch[a] - ch[b])
        print(f"  {a:22s} minus {b:18s}: mean {d.mean():+.1f} m, range {d.min():+.1f} (GIS {gis[d.argmin()]}) "
              f"to {d.max():+.1f} m (GIS {gis[d.argmax()]})")


# Toes from the 2009 start profile alone vs the default median of the 2014-2019 surveys, the module on at the
# settings before 2026-11-03 (both runs on them; HAT_collision_test_2009.py nocape, nocape_toestart)
def toes():
    runs = (("control", "Module off", *OFF_STYLE[1:]),
            ("nocape", "Module on, toes: median of the 2014-2019 surveys", *ON_STYLE[1:]),
            ("nocape_toestart", "Module on, toes: 2009 start profile only", "#E08214", 2.0))
    out, gis, ch = figure(
        2009, 2025, "HAT_collision_test_2009", runs,
        "Shoreline position change, dune toes from the start vs later surveys, 2009–2025",
        "8-other-forcings/toes_start_vs_median_position_2009_2025.png",
        "the dune collision module off (grey) and on with two sets of dune toes, fixed for the whole run: the "
        "per-cell median of the 2014-2019 lidar surveys (purple, the default) and the 2009 start profile alone "
        "(orange). Both on the module's settings before 2026-11-03 (beach width to the crest, W0 45 m, Cs 1.24e-4, "
        "no 198 m road rule), off at Cape Point GIS 1-5.")
    print(f"wrote {out}")
    for a, b in (("nocape", "control"), ("nocape_toestart", "control"), ("nocape_toestart", "nocape")):
        d = np.asarray(ch[a] - ch[b])
        print(f"  {a:16s} minus {b:8s}: mean {d.mean():+.1f} m, range {d.min():+.1f} (GIS {gis[d.argmin()]}) "
              f"to {d.max():+.1f} m (GIS {gis[d.argmax()]})")


# The same on the current settings (dry beach, W0 25 m, Cs 1.02e-4, 198 m road rule, Rodanthe 90%), 2026-11-04
def toesnow():
    runs = (("control_rd198", "Module off", *OFF_STYLE[1:]),
            ("nocape_rd198_dry", "Module on, toes: median of the 2014-2019 surveys", *ON_STYLE[1:]),
            ("nocape_rd198_dry_ts", "Module on, toes: 2009 start profile only", "#E08214", 2.0))
    out, gis, ch = figure(
        2009, 2025, "HAT_collision_test_2009", runs,
        "Shoreline position change, dune toes from the start vs later surveys, 2009–2025",
        "8-other-forcings/toes_start_vs_median_position_2009_2025_current.png",
        "the dune collision module off (grey) and on with two sets of dune toes and beach widths, fixed for the "
        "whole run: the per-cell median of the 2014-2019 lidar surveys (purple, the default) and the 2009 start "
        "profile alone (orange). Current settings in all: dry beach (beach start to toe), W0 25 m, Cs 1.02e-4, "
        "off at Cape Point GIS 1-5, the 198 m road-dune rule, villages 40% of overwash to the shoreface and "
        "Rodanthe (GIS 77-83) 90%.")
    print(f"wrote {out}")
    for a, b in (("nocape_rd198_dry", "control_rd198"), ("nocape_rd198_dry_ts", "control_rd198"),
                 ("nocape_rd198_dry_ts", "nocape_rd198_dry")):
        d = np.asarray(ch[a] - ch[b])
        print(f"  {a:20s} minus {b:16s}: mean {d.mean():+.1f} m, range {d.min():+.1f} (GIS {gis[d.argmin()]}) "
              f"to {d.max():+.1f} m (GIS {gis[d.argmax()]})")


# The beach input at the domain level (2026-11-04): domain toes with the per-cell beach widths against the domain
# beach factor (the mean over the cells of max(0, 1 - W / 25 m)); module + 198 m road rule + Rodanthe 90% in both
def widths():
    for start, end, module, base, cw, bf in ((1996, 2009, "HAT_dune_collision", "cal_step1_base", "cal_step6_toeavg_cw",
                                              "cal_step6_toeavg"),
                                             (2009, 2025, "HAT_collision_test_2009", "step1_base",
                                              "nocape_rd198_dry_toeavg_cw", "nocape_rd198_dry_toeavg")):
        runs = ((base, "Base (module off)", *OFF_STYLE[1:]),
                (cw, "Module on: per-cell beach widths", "#1B9E77", 2.0),
                (bf, "Module on: one beach factor per domain", "#7B3294", 2.2))
        out, gis, ch = figure(
            start, end, module, runs,
            f"Shoreline position change, per-cell beach widths vs one beach factor per domain, {start}–{end}",
            f"8-other-forcings/beach_widths_vs_domain_factor_position_{start}_{end}.png",
            "the base (dune collision module off, no 198 m road rule, villages 40% of overwash to the shoreface; grey) "
            "and the module with the 198 m road-dune rule and Rodanthe (GIS 77-83) 90% of overwash to the "
            "shoreface, one toe per domain from the lidar nearest the start, and the beach either per 10 m cell (the "
            "median dry beach width over the period's surveys; green) or one per domain (the mean over its cells in "
            "the start's surveys of the beach factor max(0, 1 - W / 25 m); purple).")
        d = np.asarray(ch[bf] - ch[cw])
        print(f"wrote {out}")
        print(f"  domain factor minus per-cell widths: mean {d.mean():+.1f} m, range {d.min():+.1f} "
              f"(GIS {gis[d.argmin()]}) to {d.max():+.1f} m (GIS {gis[d.argmax()]})")


# The starting dunes from the hand-picked windows (dune-topo v3) against the earlier fixed search, 8 cells landward of
# each profile's beach start (v3_win8, HAT_dune_topo_win8.py), module off and on at Cs 2e-6, toe at the berm (2026-11-11)
def win8():
    runs = (("cal_step2_growth_noscrape", "Module off, hand-picked dune windows (v3)", "#7f7f7f", 1.8),
            ("cal_step2_growth_noscrape_win8", "Module off, 8-cell dune search", "#7f7f7f", 1.6, (0, (4, 2))),
            ("cal_toeberm_cs2", "Module on (Cs 2e-6), hand-picked dune windows (v3)", "#7B3294", 2.0),
            ("cal_toeberm_cs2_win8", "Module on (Cs 2e-6), 8-cell dune search", "#E66101", 1.8, (0, (4, 2))))
    out, gis, ch = figure(
        1996, 2009, "HAT_dune_collision", runs,
        "Shoreline position change, starting dunes from hand-picked windows vs an 8-cell search, 1996–2009",
        "3-collision-design/dune-search-window/start_dunes_win8_position_1996_2009.png",
        "growth 0.02-0.04, groin 0.16, the road fix, no road rule, Rodanthe 40%, with the starting dunes and interior "
        "from the 1996-start DEM picked two ways: the highest cell in each domain's hand-picked window (dune-topo v3, "
        "the default; solid) or the highest cell within 8 cells (80 m) of each profile's beach start, the earlier "
        "Benton-parity search (v3_win8; dashed), with the road setback raised by how far seaward the 8-cell dune "
        "sits so NC-12 stays in place; the dune collision module off (grey) and on (toe at the berm, Cs 2e-6).")
    print(f"wrote {out}")
    for a, b in (("cal_step2_growth_noscrape_win8", "cal_step2_growth_noscrape"), ("cal_toeberm_cs2_win8", "cal_toeberm_cs2")):
        d = np.asarray(ch[a] - ch[b])
        print(f"  {a:32s} minus {b:26s}: mean {d.mean():+.1f} m, range {d.min():+.1f} (GIS {gis[d.argmin()]}) "
              f"to {d.max():+.1f} m (GIS {gis[d.argmax()]})")


# Two equal dune rows at the window maximum (v3) against a back row at the window median (v3_maxmed,
# HAT_dune_topo_maxmed.py), the module eroding both rows in proportion or the front row only (2026-11-11)
MAXMED_RUNS = (("cal_step2_growth_noscrape", "Module off, both rows at the window max (v3)", "#7f7f7f", 1.8),
               ("cal_step2_growth_noscrape_maxmed", "Module off, back row at the window median", "#7f7f7f", 1.6,
                (0, (4, 2))),
               ("cal_toeberm_cs2", "Module on, both rows at max, both rows eroded (current)", "#7B3294", 2.0),
               ("cal_toeberm_cs2_front", "Module on, both rows at max, front row eroded", "#2171B5", 1.8),
               ("cal_toeberm_cs2_maxmed", "Module on, back row at median, front row eroded", "#E66101", 1.8,
                (0, (4, 2))))


def maxmed():
    out, gis, ch = figure(
        1996, 2009, "HAT_dune_collision", MAXMED_RUNS,
        "Shoreline position change, a back dune row at the window median, 1996–2009",
        "3-collision-design/dune-rows/dune_rows_maxmed_position_1996_2009.png",
        "growth 0.02-0.04, groin 0.16, the road fix, no road rule, Rodanthe 40%, v3's starting dunes and interior. "
        "Barrier3D's two dune rows: both at the highest cell of each profile's hand-picked window (v3, as now), or "
        "the seaward row there and the landward row at the median elevation of the same window (dashed); the dune "
        "collision module off (grey) or on (toe at the berm, Cs 2e-6), its eroded volume taken from both rows in "
        "proportion to their height (purple, as now) or from the seaward row only (blue, orange).")
    print(f"wrote {out}")
    for a, b in (("cal_step2_growth_noscrape_maxmed", "cal_step2_growth_noscrape"), ("cal_toeberm_cs2_front", "cal_toeberm_cs2"),
                 ("cal_toeberm_cs2_maxmed", "cal_toeberm_cs2_front"), ("cal_toeberm_cs2_maxmed", "cal_toeberm_cs2")):
        d = np.asarray(ch[a] - ch[b])
        print(f"  {a:34s} minus {b:26s}: mean {d.mean():+.1f} m, range {d.min():+.1f} (GIS {gis[d.argmin()]}) "
              f"to {d.max():+.1f} m (GIS {gis[d.argmax()]})")


# The median back row with lateral erosion (the face retreats through the front row, then the back), 2026-11-11
MAXMED_LAT_RUNS = (("cal_step2_growth_noscrape", "Module off, both rows at the window max (v3)", "#7f7f7f", 1.8),
                   ("cal_toeberm_cs2", "Module on, both rows at max, proportional, Cs 2e-6 (current)", "#7B3294", 2.0),
                   ("cal_toeberm_lat_cs10", "Module on, both rows at max, lateral, Cs 1e-5", "#2171B5", 1.8),
                   ("cal_toeberm_lat_cs3_maxmed", "Back row at median, lateral, Cs 3e-6", "#FDAE61", 1.6, (0, (4, 2))),
                   ("cal_toeberm_lat_cs10_maxmed", "Back row at median, lateral, Cs 1e-5", "#F46D43", 1.8, (0, (4, 2))),
                   ("cal_toeberm_lat_cs20_maxmed", "Back row at median, lateral, Cs 2e-5", "#D73027", 1.8, (0, (4, 2))))


def maxmedlat():
    out, gis, ch = figure(
        1996, 2009, "HAT_dune_collision", MAXMED_LAT_RUNS,
        "Shoreline position change, a back dune row at the window median with lateral erosion, 1996–2009",
        "3-collision-design/dune-rows/dune_rows_maxmed_lateral_position_1996_2009.png",
        "growth 0.02-0.04, groin 0.16, the road fix, no road rule, Rodanthe 40%, v3's starting dunes and interior, "
        "toe at the berm. Grey: the dune collision module off; purple: on, both dune rows at the window maximum, the "
        "eroded volume from both rows in proportion (Cs 2e-6, the current setting); blue: both rows at the maximum, "
        "lateral erosion (the face retreats through the front row, which keeps its height until its sand above the "
        "toe is gone, then the back row; Cs 1e-5); oranges to red (dashed): the back row at the median of each "
        "profile's window, lateral erosion, Cs 3e-6, 1e-5, 2e-5.")
    print(f"wrote {out}")
    for v in ("cal_toeberm_lat_cs3_maxmed", "cal_toeberm_lat_cs10_maxmed", "cal_toeberm_lat_cs20_maxmed"):
        for b in ("cal_toeberm_cs2", "cal_toeberm_lat_cs10"):
            d = np.asarray(ch[v] - ch[b])
            print(f"  {v:28s} minus {b:22s}: mean {d.mean():+.1f} m, range {d.min():+.1f} (GIS {gis[d.argmin()]}) "
                  f"to {d.max():+.1f} m (GIS {gis[d.argmax()]})")


# Background erosion in the missed crest-drop hotspots on a source/sink sink (GIS 10-14, 21-24, 82-89), 25-100% of
# each design's own CoastSat residual (2026-11-12; HAT_be_hotspot_zones.py, whose DESIGNS give the plain names)
def behot():
    import HAT_be_hotspot_zones as Z
    for base, (name, tag, _) in Z.DESIGNS.items():
        if any(H_run_missing(v) for v, _ in Z.runs_for(base)):
            print(f"  skipping {name}: runs missing")
            continue
        runs = (("cal_step2_growth_noscrape", "No collision, equal rows", "#7f7f7f", 1.8),) * (
            base != "cal_step2_growth_noscrape") + (
            (base, name, "#1A1A1A", 1.8),
            *((f"{base}_be{p}", f"  + {p}% hotspot erosion", Z.COLOURS[p], 1.6) for p in Z.FRACS))
        out, gis, ch = figure(
            1996, 2009, "HAT_dune_collision", runs,
            f"Shoreline position change, hotspot erosion on {name.split(',')[0].lower()}, "
            f"{name.split(', ', 1)[1]}, 1996–2009",
            f"6-missed-drops/2-by-design/{Z.DESIGN_FOLDERS[base]}/hotspot_erosion_{tag}_position_1996_2009.png",
            f"{Z.SHARED}. Black: the design {name}"
            + ("; grey: no collision, equal rows" if base != "cal_step2_growth_noscrape" else "")
            + "; purples, light to dark: the same with 25, 50, 75 and 100% of the background erosion that closes "
            "its shoreline gap to CoastSat (smoothed CoastSat LRR minus its own) added in the hotspots of missed "
            "lidar crest drops that lie on a source/sink sink, GIS 10-14, 21-24 and 82-89.")
        print(f"wrote {out}")
        for v, p in Z.runs_for(base)[1:]:
            d = np.asarray(ch[v] - ch[base])
            print(f"  {name} +{p}% minus base: mean {d.mean():+.1f} m, range {d.min():+.1f} (GIS {gis[d.argmin()]}) "
                  f"to {d.max():+.1f} m (GIS {gis[d.argmax()]})")


def _simple_inputs():
    """For the plain position plots: the figure style, CoastSat's LOWESS 1996-2009 change per domain, a function
    giving a calibration run's smoothed change per domain, and the CoastSat start and end windows."""
    import matplotlib
    matplotlib.use("Agg")
    from cascade_pipeline.coastsat_lowess import CoastSatDataset, build_coastsat_series
    from cascade_pipeline.plotting.net_change_comparison import LOWESS_CONFIG, smooth_model
    from cascade_pipeline.run_layout import resolve
    from cascade_pipeline.shoreline import compute_change_rate
    from site_layer import hat_figure_style as S
    from site_layer.hat_observed_rates import NET_CHANGE_WINDOWS, net_change_dir
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS
    import HAT_be_hotspot_zones as Z
    S.apply_style()
    cs = build_coastsat_series(
        [CoastSatDataset(label="CoastSat", period_start=1996, rate_col="net_change_m",
                         csv_path=str(net_change_dir(1996, 2009) / "transect_net_change_1996_2009.csv"))],
        active_period_start=1996, lowess_config=LOWESS_CONFIG, domains=HATTERAS_DOMAINS)[0]["windows"][0]

    def change(v):
        rd = Z.H.run_dir(v)
        c = smooth_model(compute_change_rate(np.load(resolve(rd, "matrix", rd.name)), span_years=1, flip_sign=True))
        g = HATTERAS_DOMAINS.for_total(len(c))
        return c[g.start_real_index:g.end_real_index]

    return (S, cs, change) + tuple(NET_CHANGE_WINDOWS[(1996, 2009)])


def _shade_zones(Z, axes, label_ax=None, names=True, zones=None):
    """The three hotspot-erosion zones shaded on every axis, named at the top of the first (or label_ax)."""
    for zn, zr in (zones or Z.ZONES).items():
        for a in axes:
            a.axvspan(zr[0] - 0.5, zr[-1] + 0.5, color=Z.ZONE_COLOURS[zn], alpha=0.13, lw=0, zorder=0)
        if not names:
            continue
        (label_ax or axes[0]).text((zr[0] + zr[-1]) / 2, 0.985, zn.split(",")[0],
                                   transform=(label_ax or axes[0]).get_xaxis_transform(), ha="center", va="top",
                                   fontsize=8, fontweight="bold", color=Z.ZONE_COLOURS[zn])


OBS_C, OFF_C, ON_C = "#D7301F", "#7F7F7F", "#1A1A1A"     # the plain plots: observed dashed red, off grey, on black


def _position_axes(ax, S, Z, cs, gis, c_off, c_on, off, on, shade, zones=None):
    """Position change, CoastSat against the two runs, the zones shaded or not (the same top margin either way)."""
    if shade:
        _shade_zones(Z, (ax,), zones=zones)
    else:
        ax.text(0.5, 0.985, " ", transform=ax.transAxes, fontsize=8, fontweight="bold", va="top")
    ax.axhline(0, color="0.6", lw=0.6, zorder=1)
    ax.plot(cs["gis_x"], cs["smoothed"], color=OBS_C, lw=1.8, ls="--", zorder=4, label="CoastSat")
    ax.plot(gis, c_off, color=OFF_C, lw=1.8, zorder=5, label=Z.DESIGNS[off][0])
    ax.plot(gis, c_on, color=ON_C, lw=1.8, zorder=6, label=Z.DESIGNS[on][0])
    ax.set_ylabel("Shoreline position change (m)")
    ax.set_ylim(-60, 40)
    ax.set_xlim(0.5, 90.5)
    ax.set_xticks([1] + list(range(10, 91, 10)))
    ax.grid(axis="y", color="0.92", lw=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower left", ncol=3, frameon=False, fontsize=7.5)
    ax.text(1.005, 0.97, "seaward", transform=ax.transAxes, fontsize=7, color=S.INK_MUTED, va="top", rotation=90)
    ax.text(1.005, 0.03, "landward", transform=ax.transAxes, fontsize=7, color=S.INK_MUTED, va="bottom", rotation=90)


SHOAL_FC, SHOAL_EC, SHOAL_TXT = "#FBEFC4", "#D9B53A", "#8A6D0B"   # the run figures' shoal zone, yellow hatched


def _context_strip(ax, S):
    """Above an alongshore panel on the same GIS axis: the villages (grey bars) and the offshore shoals' alongshore
    footprints (yellow hatched), from the site config, named."""
    from matplotlib.patches import Rectangle
    from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS as A
    ax.set_ylim(0, 2.3)
    ax.axis("off")
    names = {"Tri-Village": "Salvo · Waves · Rodanthe"}
    for name, (lo, hi) in A.town_spans.items():
        ax.add_patch(Rectangle((lo - 0.5, 1.45), hi - lo + 1, 0.32, facecolor="0.62", edgecolor="none"))
        ax.text((lo + hi) / 2, 1.85, names.get(name, name), ha="center", va="bottom", fontsize=7.5, color=S.INK)
    for name, (lo, hi) in A.shoal_zones.items():
        ax.add_patch(Rectangle((lo - 0.5, 0.25), hi - lo + 1, 0.32, facecolor=SHOAL_FC, edgecolor=SHOAL_EC,
                               hatch="////", lw=0.6))
        ax.text((lo + hi) / 2, 0.65, name, ha="center", va="bottom", fontsize=7.5, color=SHOAL_TXT,
                fontstyle="italic")


# Buxton's rebuilt dune (2026-11-12): from 2004 the picked lidar crest reads about 0.8 m lower than the highest point
# between NC-12 and the beach, beside the dune built along the relocated road
# (1-observations/8-dune-change-by-survey/crest_offset_picked_vs_independent.png), so these cells are hatched
REBUILT_HATCH = dict(gis=(9, 14), surveys=(2004, 2005, 2009))
REBUILT_NOTE = ("Hatched, Buxton GIS 9-14 in 2004-2009: beside the dune rebuilt along the relocated NC-12 the picked "
                "lidar crest reads about 0.8 m low (crest_offset_picked_vs_independent.png), so the model looks too "
                "high there partly because of the measurement.")


def _hatch_rebuilt(ax, surveys):
    """Hatches REBUILT_HATCH's cells on a crest heatmap whose rows are `surveys` (one row per survey, from 0)."""
    from matplotlib.patches import Rectangle
    lo, hi = REBUILT_HATCH["gis"]
    rows = [i for i, y in enumerate(surveys) if y in REBUILT_HATCH["surveys"]]
    for i in rows:
        ax.add_patch(Rectangle((lo - 0.5, i - 0.5), hi - lo + 1, 1, facecolor="none", edgecolor="0.15",
                               hatch="////", lw=0.0, zorder=6))


def _rebuilt_legend(fig):
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(facecolor="white", edgecolor="0.15", hatch="////", lw=0.6,
                              label="Buxton 2004-09: rebuilt dune, lidar crest reads low")],
               loc="outside lower center", fontsize=7.5, frameon=False)


def crest_offset_style():
    """The crest-offset colour scale of HAT_collision_crest_compare.py: RdBu bands, white within +-0.5 m."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from cascade_pipeline.plotting.dune_crest_comparison import NO_SURVEY_COLOR, OFFSET_BINS_M
    half = (len(OFFSET_BINS_M) - 2) // 2
    blues = plt.get_cmap("Blues_r")(np.linspace(0.0, 0.75, half + 1))
    reds = plt.get_cmap("Reds")(np.linspace(0.25, 1.0, half + 1))
    cmap = ListedColormap(list(blues[1:]) + ["white"] + list(reds[:-1]))
    cmap.set_under(blues[0])
    cmap.set_over(reds[-1])
    cmap.set_bad(NO_SURVEY_COLOR)
    return cmap, BoundaryNorm(OFFSET_BINS_M, cmap.N), OFFSET_BINS_M


def zonescrest():
    """Both versions for a slide build: without the zones, then with them marked."""
    for both in (True, False):                     # (a)-(c), and (a)-(b) with the no-collision crests only
        _zonescrest(shade=False, both=both)
        _zonescrest(shade=True, both=both)


def _zonescrest(shade=True, both=True):
    """In 6-missed-drops/1-zones/: (a) shoreline position change 1996-2009 without the module and with the original
    design, against CoastSat; (b, c) each run's domain-median dune crest minus the lidar's at every survey, in the
    style of crest_offset_base_vs_module_1996_2009.png: too high (red) where the zones are."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    import pandas as pd
    import HAT_be_hotspot_zones as Z
    from cascade_pipeline.plotting.dune_crest_comparison import TABLE, TABLES
    S, cs, change, (s0, s1), (e0, e1) = _simple_inputs()
    gis = np.arange(1, 91)
    off, on = "cal_step2_growth_noscrape", "cal_toeberm_cs2"
    lo_g, hi_g = 6, 90                                        # where the module runs, as the crest-offset figure
    tabs = {v: (lambda t: t[t.gis.between(lo_g, hi_g)])(pd.read_csv(Z.H.run_dir(v) / TABLES / TABLE)) for v in (off, on)}
    surveys = sorted(tabs[off].survey.unique())
    cmap, norm, bins = crest_offset_style()
    crest_runs = (off, on) if both else (off,)
    n = len(crest_runs)
    fig = plt.figure(figsize=(S.FIG_W_DOUBLE, (7.4 if both else 5.4) + 0.45), layout="constrained")
    gs = fig.add_gridspec(n + 3, 1, height_ratios=(0.26, 1.55) + (1.0,) * n + (0.07,))
    ax = fig.add_subplot(gs[1])
    axs = fig.add_subplot(gs[0], sharex=ax)
    _context_strip(axs, S)
    _position_axes(ax, S, Z, cs, gis, change(off), change(on), off, on, shade, zones=Z.zones_of(on))
    ax.tick_params(axis="x", labelbottom=False)
    ax.set_title("(a) Shoreline position change", loc="left", fontsize=9, fontweight="bold")
    for k, v in enumerate(crest_runs):
        axh = fig.add_subplot(gs[k + 2], sharex=ax)
        grid = tabs[v].pivot_table(index="survey", columns="gis", values="model_minus_survey_m").reindex(
            index=surveys, columns=range(lo_g, hi_g + 1))
        im = axh.imshow(np.ma.masked_invalid(grid.values), aspect="auto", cmap=cmap, norm=norm, origin="lower",
                        extent=(lo_g - 0.5, hi_g + 0.5, -0.5, len(surveys) - 0.5), interpolation="nearest")
        _hatch_rebuilt(axh, surveys)
        axh.set_yticks(range(len(surveys)), [str(s) for s in surveys], fontsize=7)
        axh.set_ylabel("Lidar survey", fontsize=8)
        axh.tick_params(axis="x", labelbottom=(k == n - 1), labelsize=7.5)
        axh.set_title(f"({'bc'[k]}) Dune crest, {Z.DESIGNS[v][0][0].lower() + Z.DESIGNS[v][0][1:]}", loc="left", fontsize=9,
                      fontweight="bold")
        d = tabs[v].model_minus_survey_m
        axh.set_title(f"RMSE {np.sqrt((d ** 2).mean()):.2f} m, within ±0.5 m {100 * (d.abs() <= 0.5).mean():.0f}%",
                      loc="right", fontsize=8, color=S.INK_MUTED)
        if shade:                                            # outlined, so the colours underneath stay readable
            for zn, zr in Z.zones_of(on).items():
                axh.add_patch(Rectangle((zr[0] - 0.5, -0.5), len(zr), len(surveys), fill=False,
                                        edgecolor=Z.ZONE_COLOURS[zn], lw=1.8, zorder=5))
    axh.set_xlabel(S.DOMAIN_AXIS_LABEL)
    ax.set_xlim(0.5, 90.5)
    cax = fig.add_subplot(gs[n + 2])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend="both", ticks=bins)
    cb.set_label("Model − lidar, domain-median dune crest (m)   ← model lower   |   model higher →", fontsize=8)
    cb.ax.tick_params(labelsize=7.5)
    _rebuilt_legend(fig)
    stem = "hotspot_erosion_zones_position_and_crest" + ("" if both else "_nocollision") + "_1996_2009"
    out = Z.ZONES_DIR / (stem + ("" if shade else "_nozones") + ".png")
    S.save(fig, out, dpi=300, facecolor="white", close=True)
    S.record_caption(out, (
        "Shoreline position and dune crest, 1996–2009, without the collision module and with its original design, "
        f"before any background erosion ({Z.SHARED}). Above, the villages (grey; Buxton GIS 7-8, Avon 21-31, Salvo to Rodanthe 68-83) and the alongshore footprints of Avon Shoals (26-38) and Wimble Shoals (61-74), yellow hatched. (a) Shoreline position change, positive seaward. "
        "Dashed red, "
        f"CoastSat: per transect the mean position over {e0} to {e1} minus the mean over {s0} to {s1}, smoothed with "
        "a 7-domain LOWESS; grey, no collision, the two dune rows equal in height; black, the original design, lower "
        "both rows (each collision takes a volume from both dune rows, the rows equal in height, Cs 2e-6); model "
        "lines smoothed as CoastSat north of domain 10. " + ("(b) No collision and (c) the original design" if both
        else "(b) No collision") + ": the model's "
        "domain-median dune crest minus the lidar's at each survey, GIS 6-90 (where the module runs), as "
        "crest_offset_base_vs_module_1996_2009.png; red, the model too high; white, within ±0.5 m; grey, no "
        "survey; top right, the crest RMSE and the share within ±0.5 m over all domains and surveys. " + (
            f"Shaded in (a) and outlined in {'(b, c)' if both else '(b)'}, the three zones that then got background erosion (Buxton GIS "
            "9-14, Avon south 23-26, Rodanthe 82-89): missed lidar crest-drop stretches on a source/sink sink, Buxton from GIS 9 (its sink without the module) and Avon south moved two domains north of its sink (21-24) onto its missed drops."
            if shade else "The same figure without the zones marked, for a slide build.") + " " + REBUILT_NOTE))
    for zn, zr in Z.zones_of(on).items():
        r = list(zr)
        print(f"  {zn}: model - lidar, mean over surveys after 1996: " + ", ".join(
            f"{Z.DESIGNS[v][0].split(',')[0]} {tabs[v][tabs[v].gis.isin(r) & (tabs[v].survey > 1996)].model_minus_survey_m.mean():+.2f} m"
            for v in (off, on)))
    print(f"wrote {out}")


ZN2 = "cal_toeberm_cs2_zn2_be"                                 # the original design's runs on the second zone set
# The sensitivity series: (0% run, % run prefix, percentages, file tag, caption note)
SENS = {"zn2": ("cal_toeberm_cs2", "cal_toeberm_cs2_zn2_be", (25, 50, 75, 100), "zn2",
                "runs cal_toeberm_cs2_zn2_be25 to _be100"),
        "vmig": ("cal_toeberm_cs2_vmig", "cal_toeberm_cs2_zn2_vmig_be", (50, 100), "zn2_vmig",
                 "with dune migration on in the villages (their beach managers otherwise hold the dune line fixed "
                 "until the beach is gone); runs cal_toeberm_cs2_vmig, cal_toeberm_cs2_zn2_vmig_be50 and _be100")}


def zonessens():
    _zonessens("zn2")


def zonessensvmig():
    _zonessens("vmig")


def _zonessens(series):
    """In 6-missed-drops/2-by-design/<original design>/: the % sensitivity on the original design, zones 9-14, 23-26,
    82-89, in the style of the zones figures: (a) shoreline position change, CoastSat against 0-100%; (b-f) the crest
    offset from the lidar at each %."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    import pandas as pd
    import HAT_be_hotspot_zones as Z
    from cascade_pipeline.plotting.dune_crest_comparison import TABLE, TABLES
    S, cs, change, (s0, s1), (e0, e1) = _simple_inputs()
    gis = np.arange(1, 91)
    base = Z.ORIGINAL
    zones = Z.zones_of(base)
    run0, prefix, fracs, ftag, note = SENS[series]
    runs = [(run0, 0)] + [(f"{prefix}{p}", p) for p in fracs]
    missing = [v for v, _ in runs if H_run_missing(v)]
    if missing:
        raise SystemExit(f"runs missing: {missing}")
    lo_g, hi_g = 6, 90
    tabs = {v: (lambda t: t[t.gis.between(lo_g, hi_g)])(pd.read_csv(Z.H.run_dir(v) / TABLES / TABLE)) for v, _ in runs}
    surveys = sorted(tabs[run0].survey.unique())
    cmap, norm, bins = crest_offset_style()
    n = len(runs)
    fig = plt.figure(figsize=(S.FIG_W_DOUBLE, 3.0 + 1.25 * n), layout="constrained")
    gs = fig.add_gridspec(n + 3, 1, height_ratios=(0.26, 1.7) + (0.82,) * n + (0.07,))
    ax = fig.add_subplot(gs[1])
    _context_strip(fig.add_subplot(gs[0], sharex=ax), S)
    _shade_zones(Z, (ax,), zones=zones)
    ax.axhline(0, color="0.6", lw=0.6, zorder=1)
    ax.plot(cs["gis_x"], cs["smoothed"], color=OBS_C, lw=1.8, ls="--", zorder=4, label="CoastSat")
    rows = []
    for v, p in runs:
        c = change(v)
        ax.plot(gis, c, color=Z.COLOURS[p], lw=1.9 if p == 0 else 1.5, zorder=6 if p == 0 else 5,
                label=(Z.DESIGNS[base][0] + (", village dunes migrate" if series == "vmig" else "")) if p == 0
                else f"+{p}% erosion")
        for zn, zr in zones.items():
            i = np.array(list(zr)) - 1
            rows.append(dict(pct=p, zone=zn, model_m=float(c[i].mean()),
                             coastsat_m=float(np.interp(np.array(list(zr)), cs["gis_x"], cs["smoothed"]).mean())))
    ax.set_ylabel("Shoreline position change (m)")
    ax.set_ylim(-60, 40)
    ax.set_xlim(0.5, 90.5)
    ax.set_xticks([1] + list(range(10, 91, 10)))
    ax.tick_params(axis="x", labelbottom=False)
    ax.grid(axis="y", color="0.92", lw=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower left", ncol=3, frameon=False, fontsize=7.5)
    ax.text(1.005, 0.97, "seaward", transform=ax.transAxes, fontsize=7, color=S.INK_MUTED, va="top", rotation=90)
    ax.text(1.005, 0.03, "landward", transform=ax.transAxes, fontsize=7, color=S.INK_MUTED, va="bottom", rotation=90)
    ax.set_title("(a) Shoreline position change", loc="left", fontsize=9, fontweight="bold")
    for k, (v, p) in enumerate(runs):
        axh = fig.add_subplot(gs[k + 2], sharex=ax)
        grid = tabs[v].pivot_table(index="survey", columns="gis", values="model_minus_survey_m").reindex(
            index=surveys, columns=range(lo_g, hi_g + 1))
        im = axh.imshow(np.ma.masked_invalid(grid.values), aspect="auto", cmap=cmap, norm=norm, origin="lower",
                        extent=(lo_g - 0.5, hi_g + 0.5, -0.5, len(surveys) - 0.5), interpolation="nearest")
        _hatch_rebuilt(axh, surveys)
        axh.set_yticks(range(len(surveys)), [str(s) if s in (1996, 1999, 2001, 2005, 2009) else "" for s in surveys],
                       fontsize=6.5)
        axh.set_ylabel("Survey", fontsize=7.5)
        axh.tick_params(axis="x", labelbottom=(k == n - 1), labelsize=7.5)
        axh.set_title(f"({chr(98 + k)}) Dune crest, " + ("no added erosion" if p == 0 else f"+{p}% erosion"),
                      loc="left", fontsize=8.5, fontweight="bold")
        d = tabs[v].model_minus_survey_m
        inz = tabs[v][tabs[v].gis.isin([g for r in zones.values() for g in r]) & (tabs[v].survey > 1996)]
        axh.set_title(f"RMSE {np.sqrt((d ** 2).mean()):.2f} m island, {np.sqrt((inz.model_minus_survey_m ** 2).mean()):.2f} m "
                      f"in the zones", loc="right", fontsize=7.5, color=S.INK_MUTED)
        for zn, zr in zones.items():
            axh.add_patch(Rectangle((zr[0] - 0.5, -0.5), len(zr), len(surveys), fill=False,
                                    edgecolor=Z.ZONE_COLOURS[zn], lw=1.6, zorder=5))
        for zn, zr in zones.items():
            r = tabs[v][tabs[v].gis.isin(list(zr)) & (tabs[v].survey > 1996)].model_minus_survey_m
            rows.append(dict(pct=p, zone=zn, crest_offset_m=float(r.mean()), crest_rmse_m=float(np.sqrt((r ** 2).mean()))))
    axh.set_xlabel(S.DOMAIN_AXIS_LABEL)
    cax = fig.add_subplot(gs[n + 2])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend="both", ticks=bins)
    cb.set_label("Model − lidar, domain-median dune crest (m)   ← model lower   |   model higher →", fontsize=8)
    cb.ax.tick_params(labelsize=7.5)
    _rebuilt_legend(fig)
    out = Z.design_dir(base) / f"hotspot_erosion_{ftag}_position_and_crest_sensitivity_1996_2009.png"
    S.save(fig, out, dpi=300, facecolor="white", close=True)
    t = pd.DataFrame(rows).groupby(["pct", "zone"]).first().reset_index()
    (out.parent / "supporting").mkdir(parents=True, exist_ok=True)
    t.round(3).to_csv(out.parent / "supporting" / out.name.replace(".png", ".csv"), index=False)
    S.record_caption(out, (
        "Background erosion added in the hotspot zones of the original design, lower both rows, equal rows, Cs 2e-6 "
        f"({Z.SHARED}): {', '.join(str(p) for _, p in runs)}% of the erosion that closes its shoreline gap to CoastSat (smoothed "
        "CoastSat LRR minus its own) in Buxton GIS 9-14, Avon south 23-26 (moved two domains north of its sink onto "
        f"its missed crest drops, with the 21-24 rates) and Rodanthe 82-89; {note}. "
        "Above, the villages (grey) and the alongshore footprints of Avon and Wimble Shoals (yellow hatched). (a) "
        f"Shoreline position change, 1996–2009, positive seaward; dashed red, CoastSat (the mean over {e0} to {e1} "
        f"minus the mean over {s0} to {s1}, 7-domain LOWESS); black, no added erosion; purples, light to dark, "
        f"{fracs[0]}-{fracs[-1]}%; zones shaded. (b-{chr(97 + len(runs))}) The model's domain-median dune crest minus the lidar's at each survey, GIS 6-90, "
        "at each %: red, the model too high; white, within ±0.5 m; grey, no survey; zones outlined; top right, the "
        "crest RMSE over the island and over the zones (surveys after 1996). " + REBUILT_NOTE))
    print(t.round(2).to_string(index=False))
    print(f"wrote {out}")


def villagemig():
    """In 6-missed-drops/2-by-design/<original design>/: the original design with the village dune lines held fixed
    (CASCADE's default) against migrating, no added erosion: (a) position change against CoastSat, (b) the position
    difference, (c) the domain-median dune crest difference at each survey; the villages shaded."""
    import matplotlib.pyplot as plt
    import pandas as pd
    import HAT_be_hotspot_zones as Z
    from cascade_pipeline.plotting.dune_crest_comparison import TABLE, TABLES
    from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS as A
    S, cs, change, (s0, s1), (e0, e1) = _simple_inputs()
    gis = np.arange(1, 91)
    fixed, mig = "cal_toeberm_cs2", "cal_toeberm_cs2_vmig"
    c_f, c_m = change(fixed), change(mig)
    tabs = {v: pd.read_csv(Z.H.run_dir(v) / TABLES / TABLE) for v in (fixed, mig)}
    surveys = sorted(tabs[fixed].survey.unique())
    grid = {v: tabs[v].pivot_table(index="survey", columns="gis", values="model_median_m").reindex(
        index=surveys, columns=range(1, 91)) for v in (fixed, mig)}
    dcrest = grid[mig] - grid[fixed]
    fig = plt.figure(figsize=(S.FIG_W_DOUBLE, 6.9), layout="constrained")
    gs = fig.add_gridspec(5, 1, height_ratios=(0.26, 1.5, 0.75, 1.0, 0.07))
    ax = fig.add_subplot(gs[1])
    _context_strip(fig.add_subplot(gs[0], sharex=ax), S)
    axd = fig.add_subplot(gs[2], sharex=ax)
    axh = fig.add_subplot(gs[3], sharex=ax)
    for a in (ax, axd, axh):
        for name, (lo, hi) in A.town_spans.items():
            a.axvspan(lo - 0.5, hi + 0.5, color="0.90", lw=0, zorder=0)
    ax.axhline(0, color="0.6", lw=0.6, zorder=1)
    ax.plot(cs["gis_x"], cs["smoothed"], color=OBS_C, lw=1.8, ls="--", zorder=4, label="CoastSat")
    ax.plot(gis, c_f, color=OFF_C, lw=1.8, zorder=5, label="Village dune lines fixed (CASCADE default)")
    ax.plot(gis, c_m, color=ON_C, lw=1.5, zorder=6, label="Village dune lines migrate")
    ax.set_ylabel("Shoreline position\nchange (m)")
    ax.set_ylim(-60, 40)
    ax.set_xlim(0.5, 90.5)
    ax.set_xticks([1] + list(range(10, 91, 10)))
    ax.legend(loc="lower left", ncol=3, frameon=False, fontsize=7.5)
    ax.set_title("(a) Shoreline position change, 1996–2009", loc="left", fontsize=9, fontweight="bold")
    d = c_m - c_f
    axd.axhline(0, color="0.6", lw=0.6)
    axd.bar(gis, d, width=0.8, color=np.where(d < 0, "#B2182B", "#2166AC"), lw=0, zorder=3)
    axd.set_ylabel("Migrate − fixed\n(m; − landward)")
    lim = max(2.0, float(np.ceil(np.nanmax(np.abs(d)))))
    axd.set_ylim(-lim, lim)
    axd.set_title("(b) Shoreline position difference", loc="left", fontsize=9, fontweight="bold")
    vmax = 1.0
    im = axh.imshow(np.ma.masked_invalid(dcrest.values), aspect="auto", cmap="RdBu", vmin=-vmax, vmax=vmax,
                    origin="lower", extent=(0.5, 90.5, -0.5, len(surveys) - 0.5), interpolation="nearest", zorder=2)
    axh.set_yticks(range(len(surveys)), [str(s) for s in surveys], fontsize=6.5)
    axh.set_ylabel("Survey", fontsize=8)
    axh.set_title("(c) Dune crest difference, migrate − fixed", loc="left", fontsize=9, fontweight="bold")
    axh.set_xlabel(S.DOMAIN_AXIS_LABEL)
    for a in (ax, axd):
        a.tick_params(axis="x", labelbottom=False)
        a.grid(axis="y", color="0.92", lw=0.5)
        a.spines[["top", "right"]].set_visible(False)
    cax = fig.add_subplot(gs[4])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend="both")
    cb.set_label("Model dune crest, migrate − fixed (m; red lower when the dune line migrates)", fontsize=8)
    cb.ax.tick_params(labelsize=7.5)
    out = Z.design_dir(Z.ORIGINAL) / "village_dune_migration_position_and_crest_1996_2009.png"
    S.save(fig, out, dpi=300, facecolor="white", close=True)
    S.record_caption(out, (
        "The original design (lower both rows, equal rows, Cs 2e-6; " + Z.SHARED + "), no added erosion, with the "
        "village dune lines held fixed, CASCADE's default (each village beach manager keeps dune migration off until "
        "the beach width reaches its threshold, 0 m), against migrating with the shoreline as everywhere else "
        "(cal_toeberm_cs2_vmig). Grey bands, the villages (Buxton GIS 7-8, Avon 21-31, Salvo to Rodanthe 68-83); "
        "above, the villages and the shoals. (a) Shoreline position change, positive seaward; dashed red, CoastSat "
        f"(the mean over {e0} to {e1} minus the mean over {s0} to {s1}, 7-domain LOWESS). (b) Migrate minus fixed, "
        "per domain (red, the migrating run's shoreline further landward). (c) The model's domain-median dune crest, "
        "migrate minus fixed, at each lidar survey (red, lower when the dune line migrates); grey, no survey."))
    vill = [g for lo, hi in A.town_spans.values() for g in range(lo, hi + 1)]
    out_v = [g for g in range(1, 91) if g not in vill]
    for nm, gg in (("villages", vill), ("outside", out_v)):
        i = np.array(gg) - 1
        print(f"  {nm}: position diff mean {d[i].mean():+.2f} m (min {d[i].min():+.2f}), "
              f"crest diff mean {np.nanmean(dcrest[gg].values):+.2f} m (min {np.nanmin(dcrest[gg].values):+.2f})")
    for g in (21, 23, 24, 25, 26, 30, 70, 76, 80, 82, 83):
        print(f"  GIS {g}: position {d[g - 1]:+.1f} m, crest 2009 {dcrest.loc[2009, g]:+.2f} m")
    print(f"wrote {out}")


VMIG_TS_GIS = (8, 24, 25, 76, 77, 78, 80, 82, 83)     # where the village dune line moved (the crest difference map)


def villagemigcrests():
    """Beside village_dune_migration_position_and_crest: the domain-median crest through 1996-2009, village dune lines
    fixed against migrating (no added erosion), against the lidar, in the domains where the dune line moved."""
    import matplotlib.pyplot as plt
    import pandas as pd
    import HAT_be_hotspot_zones as Z
    from site_layer import hat_figure_style as S
    S.apply_style()
    runs = (("cal_toeberm_cs2", "Village dune lines fixed (CASCADE default)", OFF_C, 2.0, "-"),
            ("cal_toeberm_cs2_vmig", "Village dune lines migrate", ON_C, 1.4, "-"))
    crest = {v: Z.per_gis(v)[1] for v, *_ in runs}
    sv = pd.read_csv(Z.H.SURVEY_CRESTS)
    sv = sv[(sv.status == "ok") & sv.gis.isin(VMIG_TS_GIS)].copy()
    tt = pd.to_datetime(sv.flown.astype(str).str[:7], errors="coerce")
    sv["yr"] = tt.dt.year + (tt.dt.month - 0.5) / 12
    lid = sv.groupby(["gis", "survey"]).agg(z=("crest_m_navd88", "median"), yr=("yr", "first")).reset_index()
    fig, axs = plt.subplots(3, 3, figsize=(S.FIG_W_DOUBLE, 6.4), sharex=True, sharey=True, layout="constrained")
    for ax, gg in zip(axs.ravel(), VMIG_TS_GIS):
        for v, label, colour, lw, ls in runs:
            c = crest[v][gg]
            yrs = Z.H.RUN_START + 1 + np.arange(c.shape[0])
            ax.plot(yrs, np.median(c, axis=1), color=colour, lw=lw, ls=ls, label=label, zorder=4)
        Ld = lid[lid.gis == gg].sort_values("yr")
        ax.plot(Ld.yr, Ld.z, "o", color=OBS_C, ms=4, label="Lidar (domain median)", zorder=5)
        ax.axhline(Z.H.TOE_AT_BERM_M, color="0.6", lw=0.7, label="Berm, 1.7 m (the toe)")
        for x, _ in Z.EVENTS:
            ax.axvline(x, color="0.75", lw=0.6, ls=":")
        where = next((n for n, (lo, hi) in (("Buxton", (7, 8)), ("Avon", (21, 31)), ("Tri-Village", (68, 83)))
                      if lo <= gg <= hi), "")
        ax.set_title(f"GIS {gg}" + (f", {where}" if where else ""), loc="left", fontsize=8.5, fontweight="bold")
        ax.grid(axis="y", color="0.92", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    for x, n in Z.EVENTS:
        axs[0, 0].text(x, 0.97, f" {n} " if n != "Bonnie" else f"{n} ", transform=axs[0, 0].get_xaxis_transform(),
                       fontsize=6.5, color="0.4", va="top", ha="right" if n == "Bonnie" else "left")
    for ax in axs[:, 0]:
        ax.set_ylabel("Crest (m NAVD88)")
    for ax in axs[-1]:
        ax.set_xlabel("Year")
    axs[0, 0].set_xlim(1996, 2010)
    axs[0, 0].set_xticks(range(1996, 2011, 2))
    axs[0, 0].set_ylim(0, 9)
    h, lab = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="outside lower center", ncol=4, fontsize=7.5, frameon=False)
    out = Z.design_dir(Z.ORIGINAL) / "village_dune_migration_crests_over_time.png"
    S.save(fig, out, dpi=300, bbox_inches="tight", close=True)
    S.record_caption(out, (
        "Domain-median dune crest (m NAVD88; the higher of Barrier3D's two dune rows) through 1996-2009, the original "
        f"design (lower both rows, equal rows, Cs 2e-6; {Z.SHARED}) with no added erosion: the village dune lines held "
        "fixed, CASCADE's default (grey), against migrating with the shoreline (black; cal_toeberm_cs2_vmig), in the "
        "domains where the dune line moved (village_dune_migration_position_and_crest_1996_2009.png); red, the lidar "
        "domain median at each survey. Model state at 1 January each year; dotted verticals: Bonnie, Dennis, Isabel."))
    for gg in VMIG_TS_GIS:
        m = {v: np.median(crest[v][gg], axis=1) for v, *_ in runs}
        dd = m[runs[1][0]] - m[runs[0][0]]
        k = np.flatnonzero(np.abs(dd) > 0.05)
        print(f"  GIS {gg}: first differs {Z.H.RUN_START + 1 + k[0] if len(k) else '-'}, 2009 diff {dd[-1]:+.2f} m")
    print(f"wrote {out}")


def zonesonoff():
    """Both versions for a slide build: without the zones, then with them shaded."""
    _zonesonoff(shade=False)
    _zonesonoff(shade=True)


def _zonesonoff(shade=True):
    """In 6-missed-drops/1-zones/: shoreline position change 1996-2009 without the collision module and with the
    original design (collision a volume loss from both dune rows, the rows equal in height, Cs 2e-6), against
    CoastSat, before any background erosion; the zones shaded."""
    import matplotlib.pyplot as plt
    import HAT_be_hotspot_zones as Z
    S, cs, change, (s0, s1), (e0, e1) = _simple_inputs()
    gis = np.arange(1, 91)
    off, on = "cal_step2_growth_noscrape", "cal_toeberm_cs2"          # no collision; the original design
    c_off, c_on = change(off), change(on)
    fig = plt.figure(figsize=(S.FIG_W_DOUBLE, 4.05), layout="constrained")
    gs = fig.add_gridspec(2, 1, height_ratios=(0.26, 1.55))
    ax = fig.add_subplot(gs[1])
    _context_strip(fig.add_subplot(gs[0], sharex=ax), S)
    _position_axes(ax, S, Z, cs, gis, c_off, c_on, off, on, shade, zones=Z.zones_of(on))
    ax.set_xlabel(S.DOMAIN_AXIS_LABEL)
    out = Z.ZONES_DIR / ("hotspot_erosion_zones_position_off_vs_on_1996_2009.png" if shade else
                         "hotspot_erosion_zones_position_off_vs_on_1996_2009_nozones.png")
    S.save(fig, out, dpi=300, facecolor="white", close=True)
    S.record_caption(out, (
        "Shoreline position change, 1996–2009, without the collision module and with its original design, before "
        f"any background erosion is added ({Z.SHARED}). Above, the villages (grey; Buxton GIS 7-8, Avon 21-31, Salvo to Rodanthe 68-83) and the alongshore footprints of Avon Shoals (26-38) and Wimble Shoals (61-74), yellow hatched. Positive is seaward. Dashed red, CoastSat: per "
        "transect the "
        f"mean position over {e0} to {e1} minus the mean over {s0} to {s1}, smoothed with a 7-domain LOWESS; grey, "
        "no collision, the two dune rows equal in height; black, the original design, lower both rows: each "
        "collision takes a volume from both dune rows, the rows equal in height, Cs 2e-6; model lines smoothed as "
        "CoastSat north of domain 10. " + (
            "Shaded, the three zones that then got background erosion (Buxton GIS 9-14, Avon south 23-26, Rodanthe "
            "82-89; hotspot_erosion_zones_lowerboth_equal_cs2e-6.png): missed lidar crest-drop stretches on a "
            "source/sink sink, Buxton from GIS 9 (its sink without the module) and Avon south moved two domains north of its sink (21-24) onto its missed drops." if shade else
            "The same figure without the zones shaded, for a slide build.")))
    for zn, zr in Z.zones_of(on).items():
        i = np.array(list(zr)) - 1
        print(f"  {zn}: CoastSat {np.mean(np.interp(np.array(list(zr)), cs['gis_x'], cs['smoothed'])):+.1f} m, "
              f"off {c_off[i].mean():+.1f}, on {c_on[i].mean():+.1f}")
    print(f"wrote {out}")


def behotsimple():
    """A plainer version of behot's figures: (a) position change, CoastSat against the design and its 25-100% hotspot
    erosion, no map features, the three zones where erosion was added shaded; (b) the rate added per domain at 100%."""
    import matplotlib.pyplot as plt
    import HAT_be_hotspot_zones as Z
    import HAT_dune_collision as DC
    S, cs, change, (s0, s1), (e0, e1) = _simple_inputs()
    gis = np.arange(1, 91)
    for base, (name, tag, _) in Z.DESIGNS.items():
        if any(H_run_missing(v) for v, _ in Z.runs_for(base)):
            print(f"  skipping {name}: runs missing")
            continue
        fig, (ax, axr) = plt.subplots(2, 1, figsize=(S.FIG_W_DOUBLE, 4.6), sharex=True, layout="constrained",
                                      gridspec_kw=dict(height_ratios=(3.2, 1.0)))
        _shade_zones(Z, (ax, axr))
        ax.axhline(0, color="0.6", lw=0.6, zorder=1)
        ax.plot(cs["gis_x"], cs["smoothed"], color="#D7301F", lw=1.8, ls="--", zorder=4, label="CoastSat")
        ax.plot(gis, change(base), color="#1A1A1A", lw=1.8, zorder=6, label=name)
        for v, p in Z.runs_for(base)[1:]:
            ax.plot(gis, change(v), color=Z.COLOURS[p], lw=1.4, zorder=5, label=f"+{p}% hotspot erosion")
        ax.set_ylabel("Shoreline position\nchange (m)")
        ax.set_ylim(-60, 40)
        ax.grid(axis="y", color="0.92", lw=0.5)
        ax.legend(loc="lower left", ncol=3, frameon=False, fontsize=7, bbox_to_anchor=(0, 0.0))
        ax.text(1.005, 0.97, "seaward", transform=ax.transAxes, fontsize=7, color=S.INK_MUTED, va="top", rotation=90)
        ax.text(1.005, 0.03, "landward", transform=ax.transAxes, fontsize=7, color=S.INK_MUTED, va="bottom",
                rotation=90)
        rates = DC.BE_HOTSPOT_100[base]
        cols = [next(Z.ZONE_COLOURS[zn] for zn, zr in Z.ZONES.items() if g in zr) for g in rates]
        axr.bar(list(rates), [-r for r in rates.values()], color=cols, width=0.8, lw=0)
        axr.set_ylabel("Erosion added\nat 100% (m/yr)")
        axr.set_ylim(0, max(-r for r in rates.values()) * 1.25)
        axr.set_xlim(0.5, 90.5)
        axr.set_xticks([1] + list(range(10, 91, 10)))
        axr.set_xlabel(S.DOMAIN_AXIS_LABEL)
        for a, L in zip((ax, axr), "ab"):
            a.spines[["top", "right"]].set_visible(False)
            a.set_title(f"({L})", loc="left", fontweight="bold")
        out = Z.design_dir(base) / f"hotspot_erosion_{tag}_position_simple_1996_2009.png"
        S.save(fig, out, dpi=300, facecolor="white", close=True)
        S.record_caption(out, (
            f"Shoreline position change, 1996–2009, with background erosion added in the hotspot zones, on the design "
            f"{name} ({Z.SHARED}). (a) Positive is seaward. Dashed red, CoastSat: per transect the mean position over "
            f"{e0} to {e1} minus the mean over {s0} to {s1}, smoothed with a 7-domain LOWESS; black, the design; "
            "purples, light to dark, the same with 25, 50, 75 and 100% of the added erosion; model lines smoothed as "
            "CoastSat north of domain 10. Shaded, the three zones where erosion was added: missed lidar crest-drop "
            "stretches on a source/sink sink (Buxton GIS 10-14, Avon south 21-24, Rodanthe 82-89; "
            "hotspot_erosion_zones.png). (b) The rate added per domain at 100%: the design's own residual, smoothed "
            "CoastSat LRR minus its LRR, where negative (the model eroding too little); 25-75% scale it down. A "
            "plainer version of hotspot_erosion_<design>_position_1996_2009.png, without the map features."))
        print(f"wrote {out}")


def H_run_missing(v):
    import HAT_be_hotspot_zones as Z
    return Z.H.run_dir(v) is None


def main():
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("action", nargs="?", default="onoff", choices=["onoff", "hs", "sand", "upkeep", "village", "toes", "toesnow", "widths", "win8", "maxmed", "maxmedlat", "behot", "behotsimple", "zonesonoff", "zonescrest", "zonessens", "zonessensvmig", "villagemig", "villagemigcrests"])
    {"onoff": onoff, "hs": hs, "sand": sand, "upkeep": upkeep, "village": village, "toes": toes, "toesnow": toesnow, "widths": widths, "win8": win8, "maxmed": maxmed, "maxmedlat": maxmedlat, "behot": behot, "behotsimple": behotsimple, "zonesonoff": zonesonoff, "zonescrest": zonescrest, "zonessens": zonessens, "zonessensvmig": zonessensvmig, "villagemig": villagemig, "villagemigcrests": villagemigcrests}[ap.parse_args().action]()


if __name__ == "__main__":
    main()
