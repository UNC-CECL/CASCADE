#!/usr/bin/env python3
"""
Step 1 of the source/sink zones: what the full-management run leaves unexplained, per domain.

    python scripts/input_prep/7-source-sink/2-calibrate/be_zone_residuals.py

Writes data/hatteras_init/7-source-sink/2-calibrate/1996_2009__2009_2025/zone_residuals_edgeBE.csv
and, in 3-figures/1996_2009__2009_2025/1-field/, each figure twice: 1-rate/ (m/yr, what the zones
are picked on) and 2-net-change/ (metres from the starting shoreline, what the target is in).
Per view: both periods (zone_residuals_<view>_edgeBE.png), each period alone (..._calibration_
1996_2009, ..._test_2009_2025) and where they agree (zone_residuals_<view>_overlap_edgeBE.png).
Paths from hat_source_sink.

Residual = observed net change (CoastSat target) minus modelled net change (final minus initial
state), both at 7-domain LOWESS with GIS 1-10 raw, for the edgeBE full-management runs (no
groin) of calibration 1996-2009 and test 2009-2025. Positive: the observed shoreline ended more
seaward than the model's, so the model is missing a SOURCE there; negative, a SINK. Divided by
the run length it reads as the background rate (m/yr) a domain would need, before BRIE's
alongshore spreading, which is the unit the source/sink field is in. No new runs.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-13
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
from matplotlib.transforms import blended_transform_factory  # noqa: E402

_HERE = Path(__file__).resolve().parent
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in _HERE.parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as REPO  # noqa: E402
NC = _HAT_DIR / "input_prep" / "5-scr" / "3-rates" / "coastsat" / "net_change"
sys.path[:0] = [str(_HAT_DIR), str(NC)]
from coastsat_net_change import smooth_like_model  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK_MUTED, _title, apply_style, figsize, open_frame, record_caption,
    save, structures, town_bands)
from site_layer.hat_observed_rates import NET_CHANGE_WINDOWS, net_change_domain_csv  # noqa: E402
from site_layer.hat_source_sink import calibrate_dir, figures_dir, pair_tag  # noqa: E402
from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_ANNOTATIONS as A, HATTERAS_DOMAINS as D, HATTERAS_PERIODS, run_years)

# --- CONFIG ------------------------------------------------------------------
RAW = REPO / "output" / "raw_runs" / "matrix"
GIS = pd.RangeIndex(1, 91)
INTERIOR = slice(2, 89)
PERIODS = {
    1996: ("Calibration", RAW / "1996_2009/edgeBE/HAT_1996_2009_edgeBE_offsetmetres_road_bdm_nogroin"),
    2009: ("Test", RAW / "2009_2025/edgeBE/HAT_2009_2025_edgeBE_offsetmetres_road_bdm_nourish_nogroin"),
}
REFERENCE_M_YR = 0.5            # the old zone method's significance threshold, drawn for reference
WIDER_M_YR = 1.0                # a second reference line on the overlap figure (2026-10-12)
C_SOURCE, C_SINK = "#2166ac", "#b2182b"
YLIM = (-5, 5)                  # m/yr; test GIS 2 is +4.52
PAIR = pair_tag(1996, 2009, 2009, 2025)
TABLE = calibrate_dir(PAIR) / "zone_residuals_edgeBE.csv"
FIELD = figures_dir(PAIR) / "1-field"
# The two views of the same residual (2026-10-13): the rate the zones are picked on, and the
# net position change from the starting shoreline that the target and its scores are in
VIEWS = {
    "rate": dict(folder="1-rate", column="residual_m_yr", unit="m/yr", ylim=YLIM,
                 noun="observed minus modelled net change divided by the run length"),
    "net_change": dict(folder="2-net-change", column="residual_m", unit="m", ylim=(-80, 80),
                       noun="observed minus modelled net change from the starting shoreline, in metres"),
}
# The periods told apart by line style, not colour, so the only colours are the source/sink
# classes (blue and orange lines until 2026-10-12)
PERIOD_STYLE = {1996: dict(color="black", ls="-", lw=1.6),
                2009: dict(color="black", ls=(0, (4, 2)), lw=1.6)}
# Agreement classes for the strip: (label, colour)
AGREE = {"source": ("Source in both", C_SOURCE), "sink": ("Sink in both", C_SINK),
         "opposite": ("Opposite signs, both beyond the threshold", "0.35"),
         "weak": ("Within the threshold in at least one period", "0.90")}
# -----------------------------------------------------------------------------


# e.g. 1-rate/zone_residuals_rate_edgeBE.png (both periods), ..._edgeBE_test_2009_2025.png (one),
# zone_residuals_rate_overlap_edgeBE.png (overlap)
def png_path(view, start=None, overlap=False):
    name = f"zone_residuals_{view}{'_overlap' if overlap else ''}_edgeBE"
    if start is not None:
        name += f"_{PERIODS[start][0].lower()}_{start}_{HATTERAS_PERIODS[start]['end_year']}"
    return FIELD / VIEWS[view]["folder"] / f"{name}.png"


def observed(start):
    end = HATTERAS_PERIODS[start]["end_year"]
    d = pd.read_csv(net_change_domain_csv(start, end))
    return d.set_index("gis_domain")["net_change_lowess7_m"].reindex(GIS)


def modelled(run_dir, start):
    m = np.load(next(Path(run_dir).glob("*_shoreline_matrix.npy")))
    assert m.shape[0] == run_years(start) + 1, f"{run_dir}: unexpected run length"
    d = D.for_total(m.shape[1])   # the run's own buffers, not the config's
    return smooth_like_model(pd.Series(-(m[-1] - m[0])[d.start_real_index:d.end_real_index], index=GIS))


def shoals(ax):
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    for name, (lo, hi) in A.shoal_zones.items():
        ax.axvspan(lo - 0.5, hi + 0.5, color=A.color_shoal, alpha=0.10, zorder=1, hatch="///",
                   edgecolor=A.color_shoal, linewidth=0)
        ax.text((lo + hi) / 2.0, 0.04, name, transform=trans, ha="center", va="bottom", fontsize=7,
                color="#7A5800", style="italic", zorder=7,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))


# Both periods stacked, or one period on its own, as rates or as net change in metres. In metres
# the ±0.5 m/yr reference is that rate times each period's run length.
def figure(table, view, starts=tuple(PERIODS)):
    v = VIEWS[view]
    png = png_path(view, starts[0] if len(starts) == 1 else None)
    apply_style()
    fig, axes = plt.subplots(len(starts), 1, sharex=True, constrained_layout=True, squeeze=False,
                             figsize=figsize("double", height=3.6 + 2.8 * (len(starts) - 1)))
    axes = axes[:, 0]
    for k, (ax, start) in enumerate(zip(axes, starts)):
        role = PERIODS[start][0]
        t = table[table.period == start].set_index("gis_domain")
        r = t[v["column"]]
        ref = REFERENCE_M_YR * (run_years(start) if view == "net_change" else 1)
        ax.bar(r.index, r.values, width=0.85, zorder=3,
               color=[C_SOURCE if x >= 0 else C_SINK for x in r.values])
        for s in (-1, 1):
            ax.axhline(s * ref, color=INK_MUTED, lw=0.8, ls="--", zorder=2)
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_xlim(0.5, 90.5)
        ax.set_ylim(*v["ylim"])
        ax.grid(axis="y")
        open_frame(ax)
        town_bands(ax, label=True)
        shoals(ax)
        structures(ax, label=(k == len(starts) - 1))
        ax.set_ylabel(f"Residual ({v['unit']})\nobserved − model")
        i = r.loc[2:89]
        title = f"{role} {start}–{HATTERAS_PERIODS[start]['end_year']}"
        if len(starts) > 1:
            _title(ax, k, title)
        else:
            ax.set_title(title, loc="center")
        n_big = int((t["residual_m_yr"].loc[2:89].abs() > REFERENCE_M_YR).sum())
        stats = (f"interior mean {i.mean():+.2f} m/yr" if view == "rate" else
                 f"mean {i.mean():+.1f} m, RMS {np.sqrt((i ** 2).mean()):.1f} m")
        ax.set_title(f"{stats}, {n_big}/88 beyond ±{REFERENCE_M_YR} m/yr",
                     loc="right", fontsize=7, color=INK_MUTED)
        top = v["ylim"][1]
        for gis, x in r.items():     # beyond the axis: name the value at the edge
            if abs(x) > top:
                ax.annotate(f"{x:+.1f}", (gis, np.sign(x) * top * 0.93), ha="center",
                            fontsize=6.5, color=C_SOURCE if x > 0 else C_SINK)
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    ref_label = (f"±{REFERENCE_M_YR} m/yr" if view == "rate" else
                 f"±{REFERENCE_M_YR} m/yr × run length ("
                 + ", ".join(f"±{REFERENCE_M_YR * run_years(s):g} m {PERIODS[s][0].lower()}"
                             for s in starts) + ")")
    fig.legend(handles=[Patch(color=C_SOURCE, label="Model too landward: a missing source"),
                        Patch(color=C_SINK, label="Model too seaward: a missing sink"),
                        Line2D([], [], color=INK_MUTED, ls="--", lw=0.8, label=ref_label)],
               loc="outside lower center", ncol=3, fontsize=7, frameon=False)
    save(fig, png, dpi=300, close=True)
    parts = []
    for start in starts:
        role, run = PERIODS[start]
        (s0, s1), (e0, e1) = NET_CHANGE_WINDOWS[(start, HATTERAS_PERIODS[start]["end_year"])]
        parts.append(f"{role.lower()}, mean shoreline {e0} to {e1} minus {s0} to {s1}, run {run.name} "
                     f"({run_years(start)} yr)")
    record_caption(png, (
        "Step 1 of the source/sink zones: the net shoreline change the edgeBE full-management run "
        f"(no groin) does not explain, per GIS domain, as {v['noun']}. Blue: the observed "
        "shoreline ended more seaward than the model's, a missing source; red: more landward, a "
        "missing sink. Both sides at 7-domain LOWESS, GIS 1-10 raw. GIS 1 and 90 carry the solved "
        f"end rates. Dashed: {ref_label}, the significance threshold of the earlier zone method, "
        "for reference only. Hatched gold: the Avon and Wimble shoal zones. A bar beyond the axis "
        "is labelled with its value. Periods: " + "; ".join(parts) + "."))
    print(f"wrote {png}")


# Each domain's agreement class, from both periods' rates and the reference threshold
def agreement(cal, test, thr=REFERENCE_M_YR):
    big = (cal.abs() > thr) & (test.abs() > thr)
    out = pd.Series("weak", index=cal.index)
    out[big & (cal > 0) & (test > 0)] = "source"
    out[big & (cal < 0) & (test < 0)] = "sink"
    out[big & (np.sign(cal) != np.sign(test))] = "opposite"
    return out


# Both periods' residuals in one panel, and a strip of where they agree. The strip is always
# classed on the rates, so both views mark the same domains; in metres the reference lines are
# left out (each period's would sit at a different height).
def figure_overlap(table, view):
    v = VIEWS[view]
    png = png_path(view, overlap=True)
    apply_style()
    fig, (ax, sx) = plt.subplots(2, 1, sharex=True, constrained_layout=True,
                                 figsize=figsize("double", height=4.6),
                                 gridspec_kw=dict(height_ratios=[5, 0.6]))
    rates = table.pivot(index="gis_domain", columns="period", values="residual_m_yr")
    r = table.pivot(index="gis_domain", columns="period", values=v["column"])
    cls = agreement(rates[1996], rates[2009])
    if view == "rate":
        for s in (-1, 1):
            ax.axhline(s * REFERENCE_M_YR, color=INK_MUTED, lw=0.9, ls=":", zorder=2)
            ax.axhline(s * WIDER_M_YR, color=INK_MUTED, lw=0.8, ls="-.", zorder=2)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    for start, (role, _) in PERIODS.items():
        ax.plot(r.index, r[start], zorder=4, **PERIOD_STYLE[start],
                label=f"{role} {start}–{HATTERAS_PERIODS[start]['end_year']}")
    ax.set_xlim(0.5, 90.5)
    ax.set_ylim(*v["ylim"])
    ax.grid(axis="y")
    open_frame(ax)
    town_bands(ax, label=True)
    shoals(ax)
    structures(ax, label=True)
    ax.set_ylabel(f"Residual ({v['unit']})\nobserved − model")
    i = r.loc[2:89]
    n = cls.loc[2:89].value_counts()
    ax.set_title(f"r = {i[1996].corr(i[2009]):.2f};  source in both {n.get('source', 0)}, sink in both "
                 f"{n.get('sink', 0)}, opposite {n.get('opposite', 0)} of 88 domains",
                 loc="right", fontsize=7, color=INK_MUTED)
    for gis, c in cls.items():
        sx.axvspan(gis - 0.5, gis + 0.5, color=AGREE[c][1], lw=0)
    sx.set_ylim(0, 1)
    sx.set_yticks([])
    for side in ("left", "right", "top"):
        sx.spines[side].set_visible(False)
    sx.set_ylabel("Agree?", rotation=0, ha="right", va="center", fontsize=7)
    sx.set_xlabel(DOMAIN_AXIS_LABEL)
    # One legend in three titled columns (filled column by column): lines, references, strip
    def header(text):
        return Line2D([], [], alpha=0, label=text)
    blank = Line2D([], [], alpha=0, label=" ")
    columns = [
        [header("Residual")] + [Line2D([], [], **PERIOD_STYLE[s], label=f"{role} {s}–"
                                       f"{HATTERAS_PERIODS[s]['end_year']}")
                                for s, (role, _) in PERIODS.items()] + [blank],
        [header("Reference"), Line2D([], [], color=INK_MUTED, ls=":", lw=0.9, label=f"±{REFERENCE_M_YR} m/yr"),
         Line2D([], [], color=INK_MUTED, ls="-.", lw=0.8, label=f"±{WIDER_M_YR:g} m/yr"), blank],
        [header(f"Strip: beyond ±{REFERENCE_M_YR} m/yr in both periods")]
        + [Patch(color=AGREE[k][1], label=lab)
           for k, lab in (("source", "Source in both"), ("sink", "Sink in both"),
                          ("opposite", "Opposite signs"))],
    ]
    if view != "rate":
        del columns[1]
    leg = fig.legend(handles=[h for col in columns for h in col], loc="outside lower center",
                     ncol=len(columns), fontsize=7, frameon=False, columnspacing=3.0, handlelength=2.6)
    for txt in leg.get_texts():
        if txt.get_text() in {"Residual", "Reference"} or txt.get_text().startswith("Strip:"):
            txt.set_fontweight("bold")
            txt.set_x(-28)          # headers sit over their column, not after an empty handle
    save(fig, png, dpi=300, close=True)
    refs = ("±0.5 m/yr (dotted) is the earlier zone method's threshold, used for reference and "
            "for the strip; ±1 m/yr (dash-dot) is a second reference. " if view == "rate" else
            "No reference lines: ±0.5 m/yr is a different distance in each period (±"
            + " m and ±".join(f"{REFERENCE_M_YR * run_years(s):g}" for s in PERIODS) + " m). ")
    record_caption(png, (
        "Step 1 of the source/sink zones, the two periods together: the net-change residual the "
        f"edgeBE full-management run (no groin) leaves per GIS domain, as {v['noun']}, for "
        "calibration 1996-2009 (solid) and test 2009-2025 (dashed). Positive is a missing source, "
        "negative a missing sink. The strip classes each domain on the rates: a source in both "
        "periods, a sink in both, opposite signs with both beyond ±0.5 m/yr, or within ±0.5 m/yr in "
        "at least one period (light grey). " + refs + "r is the interior (GIS 2-89) correlation of "
        "the two patterns. Both sides at 7-domain LOWESS, GIS 1-10 raw; GIS 1 and 90 carry the "
        "solved end rates. Hatched gold: the Avon and Wimble shoal zones."))
    print(f"wrote {png}")
    return cls


def main():
    rows = []
    for start, (role, run) in PERIODS.items():
        obs, mod = observed(start), modelled(run, start)
        years = run_years(start)
        rows.append(pd.DataFrame({
            "period": start, "role": role, "gis_domain": GIS, "observed_m": obs.values,
            "model_m": mod.values, "residual_m": (obs - mod).values,
            "residual_m_yr": ((obs - mod) / years).values}))
    table = pd.concat(rows, ignore_index=True)
    TABLE.parent.mkdir(parents=True, exist_ok=True)
    table.round(3).to_csv(TABLE, index=False)
    print(f"wrote {TABLE}")
    for start, (role, _) in PERIODS.items():
        i = table[(table.period == start) & table.gis_domain.between(2, 89)]
        print(f"{role}: interior residual mean {i.residual_m.mean():+.1f} m "
              f"({i.residual_m_yr.mean():+.2f} m/yr), RMS {np.sqrt((i.residual_m ** 2).mean()):.1f} m, "
              f"|rate| > {REFERENCE_M_YR}: {(i.residual_m_yr.abs() > REFERENCE_M_YR).sum()} of 88")
    c = table.pivot(index="gis_domain", columns="period", values="residual_m_yr").loc[2:89]
    print(f"calibration vs test residual pattern, r = {c[1996].corr(c[2009]):.2f}")
    for view in VIEWS:
        figure(table, view)
        for start in PERIODS:
            figure(table, view, (start,))
        cls = figure_overlap(table, view)
    for name in ("source", "sink", "opposite"):
        print(f"{name:9s} GIS {', '.join(str(g) for g in cls.loc[2:89][cls.loc[2:89] == name].index)}")


if __name__ == "__main__":
    main()
