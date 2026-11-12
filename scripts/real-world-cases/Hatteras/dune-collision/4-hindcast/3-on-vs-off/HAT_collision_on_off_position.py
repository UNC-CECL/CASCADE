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

Writes, in the workspace's dune-collision/4-hindcast/3-on-vs-off/: collision_on_off_position_<start>_<end>.png and
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
OUT = WORKSPACE / "dune-collision" / "4-hindcast" / "3-on-vs-off"
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
        "6-other-tests/collision_hs_sweep_position_2009_2025.png",
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
        "6-other-tests/collision_sand_to_beach_position_2009_2025.png",
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
        "6-other-tests/rodanthe_dune_upkeep_position_2009_2025.png",
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
        "6-other-tests/village_overwash_position_2009_2025.png",
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
        "6-other-tests/toes_start_vs_median_position_2009_2025.png",
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
        "6-other-tests/toes_start_vs_median_position_2009_2025_current.png",
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
            f"6-other-tests/beach_widths_vs_domain_factor_position_{start}_{end}.png",
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
        "6-other-tests/dune-search-win8/start_dunes_win8_position_1996_2009.png",
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
        "6-other-tests/dune-rows-maxmed/dune_rows_maxmed_position_1996_2009.png",
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
        "6-other-tests/dune-rows-maxmed/dune_rows_maxmed_lateral_position_1996_2009.png",
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
            f"6-other-tests/hotspot-erosion/hotspot_erosion_{tag}_position_1996_2009.png",
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


def H_run_missing(v):
    import HAT_be_hotspot_zones as Z
    return Z.H.run_dir(v) is None


def main():
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("action", nargs="?", default="onoff", choices=["onoff", "hs", "sand", "upkeep", "village", "toes", "toesnow", "widths", "win8", "maxmed", "maxmedlat", "behot"])
    {"onoff": onoff, "hs": hs, "sand": sand, "upkeep": upkeep, "village": village, "toes": toes, "toesnow": toesnow, "widths": widths, "win8": win8, "maxmed": maxmed, "maxmedlat": maxmedlat, "behot": behot}[ap.parse_args().action]()


if __name__ == "__main__":
    main()
