"""
Observed against modelled shoreline change rate on the candidate windows, both smoothed at the same widths.

    python scripts/hatteras_ms/experiments/HAT_candidate_windows_plot.py

Reads the runs HAT_candidate_windows.py made (zeroBE, and the last edgeBE solve
step) and writes figures and a score table to the experiment's figures/ folder.
Details: output/raw_runs/experiments/candidate-windows/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-02
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "input_prep" / "5-scr" / "lib"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "input_prep" / "5-scr" / "3-rates"
                       / "coastsat" / "lrr_smoothed"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

import coastsat_lrr_smoothed as ls  # noqa: E402
import coastsat_lrr_smoothing_windows as sw  # noqa: E402
import coastsat_lrr_windows as cw  # noqa: E402
from cascade_pipeline.coastsat_lowess import (  # noqa: E402
    DEFAULT_DOMAINS, DEFAULT_LOWESS, spliced_lowess_series,
)
from site_layer.hat_figure_style import (  # noqa: E402
    C, DOMAIN_AXIS_LABEL, INK, _title, apply_style, caption, figsize, save,
    structures, support_dir,
)

# --- CONFIG ------------------------------------------------------------------
EXP_DIR = (PROJECT_ROOT / "output" / "raw_runs" / "experiments" / "candidate-windows"
           / "2026-10-02-candidate-windows")
OUT = EXP_DIR / "figures"
# Observed window -> model period (start, end)
WINDOWS = {"1996_2015": (1996, 2015), "2010_2026": (2010, 2026)}
# Run-folder prefix per period start, as HAT_candidate_windows.RUN_KEY (2010-2024 runs superseded 2026-10-02)
RUN_KEY = {1996: "1996", 2010: "2010_2026"}
REFERENCE = "1996_2026"
WIDTHS = (5, 7)
# The interior the run index scores, GIS 2-89
INTERIOR = (2, 89)
OBS_C = C["LATE"]
# Dark green, the model colour of the target-comparison figures
MODEL_C = "#1b6b3a"
MODEL_LS = {"edgeBE": "-", "zeroBE": (0, (4, 2.2))}
MODEL_LABEL = {"edgeBE": "Model, ends solved (edgeBE)",
               "zeroBE": "Model, no source/sink (zeroBE)"}
WINDOW_C = {"1996_2015": INK, "2010_2026": C["ACCENT"]}
# -----------------------------------------------------------------------------


# The adopted 2010 GIS 1 end (HATTERAS_BE_EDGE_ONLY), the yardstick for a forced Cape Point
from site_layer.hatteras_site_config import HATTERAS_BE_EDGE_ONLY  # noqa: E402
ADOPTED_2010_GIS1 = HATTERAS_BE_EDGE_ONLY[2010][0]


# The run directory: zeroBE's single run, or edgeBE's last solve step
def run_dir(start, preset):
    runs = EXP_DIR / "runs"
    if preset == "zeroBE":
        cands = [runs / f"{RUN_KEY[start]}_zeroBE_zero"]
    else:
        cands = sorted(runs.glob(f"{RUN_KEY[start]}_edgeBE_step*"),
                       key=lambda p: int(p.name.rsplit("step", 1)[1]))
    for d in reversed(cands):
        hits = list(d.glob("*/*/HAT_*/tables/shoreline_change_rate.csv"))
        if hits:
            return hits[0].parent.parent
    return None


# The model's per-domain LRR, raw and LOWESS-smoothed at one width the way the observations are
def model_curve(start, preset, width):
    d = run_dir(start, preset)
    if d is None:
        return None, None
    t = pd.read_csv(d / "tables" / "shoreline_change_rate.csv")
    t = t[t["gis_domain"].between(1, 90)]
    ids = t["gis_domain"].to_numpy(int)
    centres = (ids - 0.5) * DEFAULT_DOMAINS.domain_spacing_m
    s, _ = spliced_lowess_series(ids, centres, t["lrr_m_yr"].to_numpy(float), width)
    return s, d.name


# The end rates (GIS 1, GIS 90) a run imposed, from its index row
def end_rates(run_path):
    idx = pd.read_csv(PROJECT_ROOT / "output" / "raw_runs" / "run_index.csv")
    tag = ("candidate-windows/2026-10-02-candidate-windows/runs/"
           + run_path.parents[2].name)
    row = idx[(idx["run_name"] == run_path.name) & (idx["tag"] == tag)].iloc[-1]
    return float(row["be_rate_gis1_m_yr"]), float(row["be_rate_gis90_m_yr"])


# Mean and RMS of model minus observed over the interior
def score(model, obs):
    lo, hi = INTERIOR
    d = (model - obs).loc[lo:hi].dropna()
    return float(d.mean()), float(np.sqrt((d ** 2).mean()))


# The shared panel furniture
def frame(ax, half, fills, label=True):
    f = pd.DataFrame({"domain_number": np.arange(1, cw.N_DOMAINS + 1),
                      "mean_lrr": np.nan, "std_lrr": 0.0})
    cw.draw_panel(ax, f, half, std=False, label=label)
    cw.draw_shoals(ax, label=label)
    if fills:
        cw.draw_fills(ax, fills, half)
    ax.yaxis.set_major_locator(MultipleLocator(cw.Y_TICK_M))




# Panels are clipped to +-CLIP m/yr (Hannah, 2026-10-02); a value beyond gets an edge marker and its number
CLIP = 5.0
MISFIT_CLIP = 4.0
PRESETS = ("edgeBE", "zeroBE")
# edgeBE solid dark green, zeroBE light green dotted on top: they coincide away from the ends
PRESET_LS = {"edgeBE": "-", "zeroBE": (0, (1, 1.6))}
PRESET_C = {"edgeBE": MODEL_C, "zeroBE": "#74c476"}
PRESET_LABEL = {"edgeBE": "Model, ends solved (edgeBE)",
                "zeroBE": "Model, no source/sink (zeroBE)"}
# Misfit lines: colour is the width, line style the preset
WIDTH_C = {5: "#74c476", 7: MODEL_C}
MISFIT_LS = {"edgeBE": "-", "zeroBE": (0, (1, 1.6))}


# Edge markers for the points of one curve beyond the clip, labelled at its most extreme point
def mark_clipped(ax, s, clip, colour):
    s = s.dropna()
    out = s[s.abs() > clip]
    if out.empty:
        return
    up, down = out[out > 0], out[out < 0]
    ax.scatter(up.index, np.full(len(up), clip * 0.97), marker="^", s=18,
               color=colour, zorder=20, clip_on=False)
    ax.scatter(down.index, np.full(len(down), -clip * 0.97), marker="v", s=18,
               color=colour, zorder=20, clip_on=False)
    g = out.abs().idxmax()
    # One label per spot: skip if another curve already labelled within 4 domains on this edge
    done = ax.__dict__.setdefault("_clip_labels", [])
    if any(abs(g - g0) < 4 and np.sign(out[g]) == s0 for g0, s0 in done):
        return
    done.append((g, np.sign(out[g])))
    ax.annotate(f"{out[g]:+.1f}", (g, np.sign(out[g]) * clip * 0.97),
                xytext=(5, -9 if out[g] > 0 else 3), textcoords="offset points",
                fontsize=7, color=colour, zorder=20)


# The panel's furniture on the clipped axis
def clipped_frame(ax, label):
    frame(ax, CLIP, None, label=label)
    ax.set_ylim(-CLIP, CLIP)
    ax.yaxis.set_major_locator(MultipleLocator(1))


# A quiet alongshore frame: village bands, faint shoal tints, structure hairlines; names only where asked
def quiet_frame(ax, clip, label):
    ax.set_xlim(0.5, cw.N_DOMAINS + 0.5)
    ax.set_ylim(-clip, clip)
    cw.town_bands(ax, label=label, shade="0.95")
    for name, (lo, hi) in cw.HATTERAS_ANNOTATIONS.shoal_zones.items():
        ax.axvspan(lo - 0.5, hi + 0.5, color=cw.SHOAL_C, alpha=0.10, lw=0, zorder=0.5)
        if label:
            ax.text((lo + hi) / 2, 0.03, name, transform=ax.get_xaxis_transform(),
                    ha="center", va="bottom", fontsize=7, color=cw.SHOAL_TEXT)
    ax.axhline(0, color=cw.INK_MUTED, lw=0.7, zorder=2)
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.grid(True, which="major", color="0.88", lw=0.6, zorder=0)
    ax.yaxis.grid(True, which="minor", color="0.94", lw=0.4, zorder=0)
    ax.set_axisbelow(True)
    cw.open_frame(ax)


# The role each candidate window plays (Hannah, 2026-10-02)
WINDOW_ROLE = {"1996_2015": "Calibration window", "2010_2026": "Test window"}


# Figure 1: (a) the calibration window above (b) the test window, both widths in each (Hannah, 2026-10-02)
# Observed blue, edgeBE green, zeroBE grey; the wider width dark, the narrower light
FIT_C = {("obs", 7): OBS_C, ("obs", 5): "#9ecae1",
         ("edgeBE", 7): MODEL_C, ("edgeBE", 5): "#74c476",
         ("zeroBE", 7): "0.35", ("zeroBE", 5): "0.70",
         ("ref", 7): "0.10", ("ref", 5): "0.55"}
REF_DASH = (0, (5, 2.5))


# The figure folders by model setup (Hannah, 2026-10-02)
BOTH = "edgeBE_and_zeroBE"
# Inside each setup folder: one subfolder per window, one for the figures showing both
SUBDIR = {"1996_2015": "1-calibration_1996_2015", "2010_2026": "2-test_2010_2026",
          None: "3-both_windows"}


# A preset set's folder and file token
def preset_tag(presets):
    return BOTH if len(presets) > 1 else presets[0]


# A model line's colour: grey for zeroBE beside edgeBE, otherwise the model greens
def model_c(preset, w, presets):
    if len(presets) > 1:
        return FIT_C[(preset, w)]
    return FIT_C[("edgeBE", w)]


# `only` draws one window on its own figure (Hannah, 2026-10-02)
def fit_panels(with_reference=False, presets=PRESETS, only=None):
    a, b = WIDTHS
    tag = preset_tag(presets)
    rows = [(stem, m) for stem, m in WINDOWS.items() if only in (None, stem)]
    fig, axes = plt.subplots(len(rows), 1, sharex=True, squeeze=False,
                             figsize=figsize("double", height=5.6 if len(rows) > 1 else 3.6),
                             constrained_layout=True)
    axes = axes[:, 0]
    for k, (ax, (stem, (m0, m1))) in enumerate(zip(axes, rows)):
        quiet_frame(ax, CLIP, label=(k == 0))
        for w in WIDTHS:
            obs, unc = ls.curve_and_band(stem, w)
            x = obs.index.to_numpy(float)
            if w == b:
                ax.fill_between(x, obs - unc, obs + unc, color=OBS_C, alpha=0.18,
                                lw=0, zorder=8)
            ax.plot(x, obs, color=FIT_C[("obs", w)], lw=1.6 if w == b else 1.3,
                    zorder=12 if w == b else 11, solid_capstyle="round")
            mark_clipped(ax, obs, CLIP, FIT_C[("obs", w)])
            if with_reference:
                ref, _ = ls.curve_and_band(REFERENCE, w)
                ax.plot(ref.index.to_numpy(float), ref, color=FIT_C[("ref", w)],
                        lw=1.2, ls=REF_DASH, zorder=11.5)
                mark_clipped(ax, ref, CLIP, FIT_C[("ref", w)])
            for preset, z in (("zeroBE", 9), ("edgeBE", 10)):
                if preset not in presets:
                    continue
                c = model_c(preset, w, presets)
                thin = preset == "zeroBE" and len(presets) > 1
                m, _ = model_curve(m0, preset, w)
                ax.plot(m.index.to_numpy(float), m, color=c,
                        lw=(1.0 if thin else 1.4) * (1 if w == b else 0.9),
                        zorder=z + (0.5 if w == b else 0), solid_capstyle="round")
                mark_clipped(ax, m, CLIP, c)
        model_txt = "" if stem == f"{m0}_{m1}" else f"  (model {m0}–{m1})"
        letter = f"({chr(97 + k)})  " if len(rows) > 1 else ""
        ax.set_title(f"{letter}{WINDOW_ROLE[stem]}, "
                     f"{stem.replace('_', '–')}{model_txt}", loc="left", fontsize=9)
        structures(ax, label=(k == 0))
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    if len(rows) > 1:
        fig.supylabel(cw.Y_LABEL, fontsize=9)
    else:
        axes[0].set_ylabel(cw.Y_LABEL, fontsize=9)
    fig.suptitle(f"Observed and modelled ({' and '.join(presets)}) shoreline change "
                 f"rate, LOWESS {sw.km_of(a):g} km and {sw.km_of(b):g} km", fontsize=10)
    h, labels = [], []
    keys = [("obs", "Observed (CoastSat)"),
            ("edgeBE", "Model, ends solved (edgeBE)"),
            ("zeroBE", "Model, no source/sink (zeroBE)")]
    if with_reference:
        # Four groups need short names to fit across the page
        keys = [("obs", "Observed, window"),
                ("ref", f"Observed {REFERENCE.replace('_', '–')}"),
                ("edgeBE", "Model edgeBE"), ("zeroBE", "Model zeroBE")]
    keys = [(key, name) for key, name in keys if key not in PRESETS or key in presets]
    for key, name in keys:
        for w in (b, a):
            c = model_c(key, w, presets) if key in PRESETS else FIT_C[(key, w)]
            thin = (key == "zeroBE" and len(presets) > 1) or key == "ref"
            h.append(Line2D([], [], color=c, lw=1.2 if thin else 1.6,
                            ls=REF_DASH if key == "ref" else "-"))
            labels.append(f"{name}, {w} domains")
    fig.legend(h, labels, loc="outside lower center", ncol=len(keys),
               frameon=False, fontsize=8 if len(keys) <= 3 else 7)
    e96, e10 = end_rates(run_dir(1996, "edgeBE")), end_rates(run_dir(2010, "edgeBE"))
    caption(fig, (
        "Observed and modelled shoreline change rate by GIS domain (1 at Cape "
        "Point, 90 at Pea Island) on "
        + ("the calibration window (a, 1996–2015) and the test window (b, 2010–2026)"
           if len(rows) > 1 else
           f"the {WINDOW_ROLE[rows[0][0]].lower()} ({rows[0][0].replace('_', '–')})")
        + ", each LOWESS-smoothed alongshore at "
        f"{sw.km_of(b):g} km ({b} domains, dark) and {sw.km_of(a):g} km ({a} "
        "domains, light), both sides the same way, with raw domain means kept over "
        f"GIS 1–{DEFAULT_LOWESS.skip_southern_domains}. Blue is the CoastSat linear "
        f"regression rate over the calendar window, with the 95% band on the {b}-domain "
        "curve (ordinary least-squares, too narrow for a serially correlated record). "
        + ends_txt(presets, e96, e10)
        + (f"The dashed line in both panels is the long-term observed "
           f"{REFERENCE.replace('_', '–')} rate, fitted and smoothed the same way. "
           if with_reference else "")
        + "Both model periods match their observed windows; the test-window run is "
        "forced through 2025 with the Duck water levels and WIS waves extended past "
        "the 1984–2024 records. "
        + (f"Its edgeBE GIS 1 end is {e10[0] / ADOPTED_2010_GIS1:.0f} times the "
           "adopted value, so that model's agreement at Cape Point is forced. "
           if "edgeBE" in presets and e10[0] > 3 * ADOPTED_2010_GIS1 else "")
        + 
        "All runs are full management (fills in the window included), no groin, no "
        "relocation. Grey bands are the villages, amber tints the offshore shoals; the "
        "solid hairline is the Buxton groin, the dotted ones the Avon and Rodanthe "
        f"piers. The y axis is clipped at ±{CLIP:g} m/yr; a value beyond it is marked "
        "by a triangle at the edge, with the largest labelled."))
    span = "_and_".join(st for st, _ in rows)
    stem = f"fit_{tag}_{span}_lowess5_and_lowess7"
    out = save(fig, OUT / tag / SUBDIR[only] / (stem + (f"_with_{REFERENCE}"
                                                        if with_reference else "")))
    plt.close(fig)
    return out


# What the model lines are, for a fit caption
def ends_txt(presets, e96, e10):
    ends = (f"(1996–2015: GIS 1 {e96[0]:+.2f}, GIS 90 {e96[1]:+.2f}; 2010–2026: "
            f"GIS 1 {e10[0]:+.2f}, GIS 90 {e10[1]:+.2f} m/yr)")
    if len(presets) > 1:
        return ("Green is edgeBE, which carries end rates solved against the observed "
                f"window {ends}; grey is zeroBE, no source or sink, drawn underneath "
                "and visible only where it leaves edgeBE near the ends. ")
    if presets[0] == "edgeBE":
        return ("Green is the model with end rates solved against the observed window "
                f"(edgeBE) {ends}. ")
    return "Green is the model with no source or sink anywhere (zeroBE). "


# Figure 2: model minus observed, one panel per window; colour the width, line style the preset
def misfit(presets=PRESETS):
    tag = preset_tag(presets)
    fig, axes = plt.subplots(2, 1, sharex=True, sharey=True,
                             figsize=figsize("double", height=4.6),
                             constrained_layout=True)
    for k, (ax, (stem, (m0, m1))) in enumerate(zip(axes, WINDOWS.items())):
        ax.set_xlim(0.5, cw.N_DOMAINS + 0.5)
        ax.set_ylim(-MISFIT_CLIP, MISFIT_CLIP)
        cw.town_bands(ax, label=(k == 0))
        cw.draw_shoals(ax, label=(k == 0))
        ax.axhline(0, color=cw.INK_MUTED, lw=0.8, zorder=2)
        for w in WIDTHS:
            obs, _ = ls.curve_and_band(stem, w)
            for preset in presets:
                m, _ = model_curve(m0, preset, w)
                d = (m - obs).dropna()
                ax.plot(d.index.to_numpy(float), d, color=WIDTH_C[w], lw=1.3,
                        ls=MISFIT_LS[preset], zorder=10 + w)
                mark_clipped(ax, d, MISFIT_CLIP, WIDTH_C[w])
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_minor_locator(MultipleLocator(5))
        ax.yaxis.grid(True, zorder=0)
        ax.set_axisbelow(True)
        cw.open_frame(ax)
        _title(ax, k, f"Observed {stem.replace('_', '–')}, model {m0}–{m1}")
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    fig.supylabel("Model minus observed (m/yr)", fontsize=9)
    fig.suptitle(f"Model ({' and '.join(presets)}) minus observed shoreline change "
                 "rate on the candidate windows", fontsize=10)
    h, labels = [], []
    for preset in presets:
        for w in WIDTHS:
            h.append(Line2D([], [], color=WIDTH_C[w], lw=1.3, ls=MISFIT_LS[preset]))
            labels.append(f"{preset}, {w} domains")
    fig.legend(h, labels, loc="outside lower center", ncol=len(h), frameon=False)
    what = ("line style is the model, solid edgeBE (ends solved against the window) "
            "and dotted zeroBE (no source or sink). Away from the ends the two models "
            "coincide and" if len(presets) > 1 else
            ("the model is edgeBE (ends solved against the window)." if presets[0]
             == "edgeBE" else "the model is zeroBE (no source or sink).")
            + " Away from the ends")
    caption(fig, (
        "Model minus observed shoreline change rate by GIS domain, both sides "
        "LOWESS-smoothed at the same width: below zero the model erodes more (or "
        "accretes less) than the shoreline did. Colour is the width, light green "
        f"2.5 km (5 domains) and dark green 3.5 km (7 domains); {what} the width "
        "moves a line by tenths of a metre a year. The y axis is clipped at "
        f"±{MISFIT_CLIP:g} m/yr; a value beyond it is marked by a triangle at the "
        "edge, with the largest labelled."))
    out = save(fig, OUT / tag / SUBDIR[None]
               / f"misfit_{tag}_1996_2015_and_2010_2026_lowess5_and_lowess7")
    plt.close(fig)
    return out


# The score table: bias and RMSE over GIS 2-89 for every window x preset x width
def score_table():
    rows = []
    for stem, (m0, m1) in WINDOWS.items():
        for w in WIDTHS:
            obs, _ = ls.curve_and_band(stem, w)
            for preset in PRESETS:
                m, name = model_curve(m0, preset, w)
                bias, rms = score(m, obs)
                rows.append(dict(window=stem, model_period=f"{m0}_{m1}",
                                 preset=preset, width_domains=w,
                                 width_km=sw.km_of(w), run=name,
                                 mean_bias_m_yr=round(bias, 3),
                                 rmse_m_yr=round(rms, 3)))
    return pd.DataFrame(rows)


# Figure 3: the scores as dots; rows window x width, filled edgeBE, open zeroBE
def score_dots(table):
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=figsize("double", aspect=0.32),
                             constrained_layout=True)
    keys = [(stem, w) for stem in WINDOWS for w in WIDTHS]
    ylab = [f"{stem.replace('_', '–')}, {w} domains" for stem, w in keys]
    y = np.arange(len(keys))[::-1]
    for k, (ax, col, title) in enumerate(zip(
            axes, ("mean_bias_m_yr", "rmse_m_yr"),
            ("Mean bias, model minus observed", "RMSE"))):
        for (stem, w), yy in zip(keys, y):
            for preset in PRESETS:
                r = table[(table.window == stem) & (table.width_domains == w)
                          & (table.preset == preset)]
                v = float(r[col].iloc[0])
                ax.scatter(v, yy, s=36, lw=1.2, zorder=5,
                           facecolors=WINDOW_C[stem] if preset == "edgeBE" else "white",
                           edgecolors=WINDOW_C[stem])
        if col == "mean_bias_m_yr":
            ax.axvline(0, color=cw.INK_MUTED, lw=0.8, zorder=1)
        ax.set_yticks(y)
        ax.set_yticklabels(ylab)
        ax.xaxis.grid(True, zorder=0)
        ax.set_axisbelow(True)
        ax.set_xlabel("m/yr")
        cw.open_frame(ax)
        _title(ax, k, title)
    h = [Line2D([], [], marker="o", lw=0, color=INK, markerfacecolor=INK),
         Line2D([], [], marker="o", lw=0, color=INK, markerfacecolor="white")]
    fig.legend(h, ["edgeBE", "zeroBE"], loc="outside lower center", ncol=2,
               frameon=False)
    caption(fig, (
        "Model skill on the candidate windows over GIS 2–89, after both sides are "
        "LOWESS-smoothed at the width in the row label: (a) mean of model minus "
        "observed, (b) root mean square of it, m/yr. Filled dots are edgeBE (ends "
        "solved against the window), open dots zeroBE (no source or sink); black is "
        "1996–2015, purple 2010–2026. Values are in "
        "supporting/observed_vs_model_scores.csv."))
    out = save(fig, OUT / BOTH / SUBDIR[None]
               / f"scores_{BOTH}_1996_2015_and_2010_2026_lowess5_and_lowess7")
    plt.close(fig)
    return out


# Run: the score table and the three figures
def main() -> int:
    apply_style()
    table = score_table()
    csv = support_dir(OUT / BOTH / SUBDIR[None]) / "observed_vs_model_scores.csv"
    table.to_csv(csv, index=False)
    written = []
    for presets in (PRESETS, ("edgeBE",), ("zeroBE",)):
        for only in (None, *WINDOWS):
            written += fit_panels(presets=presets, only=only)
            written += fit_panels(with_reference=True, presets=presets, only=only)
        written += misfit(presets)
    for p in written + score_dots(table):
        print("wrote", p.relative_to(PROJECT_ROOT))
    print(table.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
