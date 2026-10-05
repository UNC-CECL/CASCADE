"""
The 1996-2025 full-window run against CoastSat: the LRR rate, the net position change, and a yearly GIF.

    python scripts/analyze_output/compare_runs/hindcast_vs_observed/full_window.py

Rate: the modelled OLS rate against the CoastSat LRR 1996-2025. Net change: the
model's 1 Jan 2026 shoreline minus its start, against the CoastSat calendar-2025
mean minus the DEM-centred 1996 start mean. GIF: one frame per year, model change
since the start beside CoastSat's calendar-year mean change since the start.
Both sides 7-domain LOWESS (GIS 1-10 raw) in the smoothed reading. Writes
output/comparisons/full_window/1996_2025/.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-05
"""
from __future__ import annotations

import io
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from PIL import Image  # noqa: E402

_REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts" / "analyze_output" / "compare_runs" / "matrix_vs_observed"))
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import matrix_vs_observed as mvo  # noqa: E402  (rate table, smoothing, scoring helpers)
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)
from coastsat_vs_duneline import load_chainage  # noqa: E402

from site_layer.hat_figure_style import (  # noqa: E402
    COMPARISONS_ROOT, DOMAIN_AXIS_LABEL, INK_MUTED, _title, apply_style,
    compare_header, figsize, mark_offaxis, offaxis_clause, open_frame,
    record_caption, save, structures, town_bands)
from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT  # noqa: E402
from site_layer.hat_topo_version import INIT_ROOT  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS  # noqa: E402

common = mvo.common

# --- CONFIG ------------------------------------------------------------------
START, END = 1996, 2025
WINDOW = f"{START}_{END}"
RUN_ROOT = _REPO / "output" / "raw_runs" / "matrix" / WINDOW
SCENARIO_RUN = "HAT_1996_2025_{preset}_offsetmetres_road_bdm_nourish_nogroin"
PRESETS = ("zeroBE", "edgeBE")
# The DEM-centred start mean (it built the shoreline offset) and the end year (calendar 2025, Hannah 2026-10-05)
START_MEAN = "1995-10-12_1997-10-12"
OBS_END_YEAR = 2025
MEAN_ROOT = INIT_ROOT / "5-scr" / "1-observations" / "mean_shoreline"
OUT = COMPARISONS_ROOT / "full_window" / WINDOW
N_DOMAINS = 90
RATE_HALF = 7.5                     # m/yr, as matrix_vs_observed
POS_HALF = 100.0                    # m, as target_comparison
C_OBS = "#92c5de"                   # observed pale and thick, model dark and thin
C_MODEL = "#2166ac"
LW_OBS, LW_MODEL = 2.8, 1.2
GIF_MS = 450                        # per frame; the last frame holds three times as long
PRESET_TEXT = {"zeroBE": "no source/sink in any domain (zeroBE)",
               "edgeBE": "end rates at GIS 1 and 90 solved on the 1996-2025 CoastSat LRR (edgeBE)"}
# -----------------------------------------------------------------------------


# The run directory of a preset, or None until it has run
def run_dir(preset):
    d = RUN_ROOT / preset / SCENARIO_RUN.format(preset=preset)
    return d if (d / "tables" / "shoreline_change_rate.csv").is_file() else None


# Per-domain series indexed by GIS
def _gis(values):
    return pd.Series(values, index=pd.RangeIndex(1, N_DOMAINS + 1, name="gis_domain"), dtype=float)


# Every annual model shoreline as change from the start, seaward +; row k is 1 Jan of START + k
def model_changes(rd):
    m = np.load(next(rd.glob("*_shoreline_matrix.npy")))
    D = HATTERAS_DOMAINS
    return [_gis(-(m[k] - m[0])[D.start_real_index:D.end_real_index]) for k in range(len(m))]


# Per transect: the start mean and every calendar-year mean, as change from the start
def observed_yearly():
    t = pd.read_csv(MEAN_ROOT / START_MEAN / f"transect_means_{START_MEAN}.csv")
    t = t[t["domain_number"].between(1, N_DOMAINS)].set_index("transect_id")
    start = t["mean_chainage_m"].where(t["included"].astype(bool))
    rows = {}
    for tid in start.index:
        df = load_chainage(tid)
        if df is None or not np.isfinite(start[tid]):
            continue
        y = df.groupby(df["date"].dt.year)["chainage"].agg(["mean", "size"])
        rows[tid] = (y["mean"] - start[tid]).reindex(range(START, END + 1))
    change = pd.DataFrame(rows).T
    change.insert(0, "domain_number", t["domain_number"].astype(int).reindex(change.index))
    return change


# Domain means of one year's transect changes
def domain_year(change, year):
    return change.groupby("domain_number")[year].mean().reindex(range(1, N_DOMAINS + 1)).pipe(_gis)


# Interior bias, RMSE and r of model minus observation
def skill(model, obs):
    m, o = common.interior(model), common.interior(obs)
    ok = m.notna() & o.notna()
    d = (m - o)[ok]
    return dict(n=int(ok.sum()), bias=float(d.mean()), rmse=float(np.sqrt((d ** 2).mean())),
                r=float(np.corrcoef(m[ok], o[ok])[0, 1]))


# Alongshore axis: zero line, GIS 1-90, grid, village bands
def _axis(ax, half, label_towns):
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(1, N_DOMAINS)
    ax.set_ylim(-half, half)
    ax.grid(axis="y")
    open_frame(ax)
    town_bands(ax, label=label_towns)


# Observed and model on one panel; returns the off-axis clause for the caption
def _panel(ax, i, title, obs, model, half, unit, s, obs_label):
    _axis(ax, half, i == 0)
    fmt = ".2f" if unit == "m/yr" else ".1f"
    ax.plot(obs.index, obs.values, color=C_OBS, lw=LW_OBS, zorder=5)
    ax.plot(model.index, model.values, color=C_MODEL, lw=LW_MODEL, zorder=6)
    handles = [Line2D([], [], color=C_OBS, lw=LW_OBS, label=obs_label),
               Line2D([], [], color=C_MODEL, lw=LW_MODEL,
                      label=f"CASCADE: bias {s['bias']:+{fmt}} {unit}, RMSE {s['rmse']:{fmt}}, "
                            f"r {s['r']:.2f}")]
    off = [("observed", mark_offaxis(ax, obs.index, obs.values, half, color=C_OBS)),
           ("model", mark_offaxis(ax, model.index, model.values, half, color=C_MODEL))]
    _title(ax, i, title)
    ax.legend(handles=handles, loc="lower center", fontsize=6.5, frameon=False)
    return offaxis_clause(off, half, unit)


# Rate and net change for one preset and one reading
def figure(preset, rd, obs, model, smoothed):
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.6), sharex=True,
                             constrained_layout=True)
    reading = "7-domain LOWESS, GIS 1-10 raw" if smoothed else "domain means"
    clauses = [
        _panel(axes[0], 0, "Shoreline change rate, 1996-2025", obs["rate"], model["rate"],
               RATE_HALF, "m/yr", skill(model["rate"], obs["rate"]), "CoastSat LRR 1996-2025"),
        _panel(axes[1], 1, "Net position change, 1996 start to 2025", obs["net"], model["net"],
               POS_HALF, "m", skill(model["net"], obs["net"]),
               "CoastSat: calendar-2025 mean minus the 1996 start mean"),
    ]
    axes[0].set_ylabel("Change rate,\nLRR (m/yr)")
    axes[1].set_ylabel("Position change,\nfrom 1996 start (m)")
    structures(axes[0], label=False)
    structures(axes[1], label=True)
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    compare_header(fig, f"Full management, no groin, {preset}, 1996-2025: CASCADE against "
                        f"CoastSat ({reading}, both sides)")
    tag = "smoothed" if smoothed else "domain_means"
    png = OUT / tag / f"full_window_{preset}_full_management_{WINDOW}_{tag}.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        f"Full management (road, beach and dune manager, fills), no groin, no relocation; "
        f"{PRESET_TEXT[preset]}. Run {rd.name}: 30 model years, 1996 through 2025, final state "
        f"1 Jan 2026, starting from the shoreline offset built on the CoastSat mean over "
        f"{START_MEAN.replace('_', ' to ')} (the 1996 ALACE survey +/-1 yr). Pale line "
        "CoastSat, dark line CASCADE. (a) The modelled OLS rate against the CoastSat LRR "
        "1996-2025 (observations 1996-01-26 to 2025-12-28). (b) The modelled end minus start "
        "shoreline against the CoastSat calendar-2025 mean minus the start mean; the "
        "calendar-2025 mean is centred about half a year before the model end, and is a 1-yr "
        f"mean against a 2-yr start mean. Both sides are {reading}. Seaward positive; scores "
        "over the interior GIS 2-89." + "".join(clauses)))
    return png


# One GIF frame as an image: model and observed change since the start for one year
def _frame(year, model, obs, end_obs, preset):
    fig, ax = plt.subplots(figsize=figsize("double", height=3.5), layout="constrained")
    _axis(ax, POS_HALF, True)
    ax.plot(end_obs.index, end_obs.values, color=INK_MUTED, lw=0.9, ls=(0, (4, 2)), zorder=4)
    if obs is not None:
        ax.plot(obs.index, obs.values, color=C_OBS, lw=LW_OBS, zorder=5)
        mark_offaxis(ax, obs.index, obs.values, POS_HALF, color=C_OBS)
    ax.plot(model.index, model.values, color=C_MODEL, lw=LW_MODEL, zorder=6)
    mark_offaxis(ax, model.index, model.values, POS_HALF, color=C_MODEL)
    structures(ax, label=True)
    ax.set_ylabel("Position change,\nfrom 1996 start (m)")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    # Below the axes, so it never sits on the structure labels
    fig.legend(handles=[
        Line2D([], [], color=C_OBS, lw=LW_OBS, label=f"CoastSat, {year} mean"),
        Line2D([], [], color=C_MODEL, lw=LW_MODEL, label=f"CASCADE, 1 Jan {year + 1}"),
        Line2D([], [], color=INK_MUTED, lw=0.9, ls=(0, (4, 2)),
               label="CoastSat, 2025 mean (the end target)")],
        loc="outside lower center", fontsize=7, frameon=False, ncol=3)
    ax.set_title(f"{year}", loc="left", fontsize=10)
    compare_header(fig, f"Full management, no groin, {preset}: change since the 1996 start "
                        "(7-domain LOWESS, GIS 1-10 raw, both sides)")
    buf = io.BytesIO()
    fig.savefig(buf, dpi=150, format="png")
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGB")


# The yearly GIF: frame for calendar year Y pairs CoastSat's Y mean with the model on 1 Jan Y+1
def gif(preset, rd, models, change):
    end_obs = mvo.smoothed(domain_year(change, OBS_END_YEAR))
    frames = []
    for year in range(START, END + 1):
        o = domain_year(change, year)
        obs = mvo.smoothed(o) if o.notna().sum() > 10 else None
        frames.append(_frame(year, mvo.smoothed(models[year - START + 1]), obs, end_obs, preset))
    path = OUT / "gif" / f"full_window_{preset}_full_management_{WINDOW}_yearly_change.gif"
    path.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(path, save_all=True, append_images=frames[1:], loop=0,
                   duration=[GIF_MS] * (len(frames) - 1) + [GIF_MS * 3])
    record_caption(path.with_suffix(".png"), (
        f"GIF `{path.name}`. Full management, no groin, {PRESET_TEXT[preset]}; run {rd.name}. "
        "One frame per calendar year 1996-2025: pale, the CoastSat mean over that calendar "
        "year minus the 1996 start mean; dark, the modelled shoreline on 1 Jan of the next "
        "year minus the model start; dashed, the calendar-2025 CoastSat mean (the end "
        "target). Each CoastSat year is centred about half a year before its model state. "
        "Both sides 7-domain LOWESS, GIS 1-10 raw. Seaward positive; y held at +/-100 m, "
        "off-axis values marked at the edge."))
    return path


# Run: every preset that has run, both readings, the GIF, scores and the domain table
def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "tables").mkdir(exist_ok=True)
    change = observed_yearly()
    change.round(4).to_csv(OUT / "tables" / "transect_yearly_change_from_1996_start.csv")
    lrr = pd.read_csv(COASTSAT_LRR_ROOT / WINDOW / "domain_lrr_summary.csv")
    obs_raw = {"rate": _gis(lrr.set_index("domain_number")["mean_lrr"]
                            .reindex(range(1, N_DOMAINS + 1)).values),
               "net": domain_year(change, OBS_END_YEAR)}
    obs_smooth = {"rate": common.coastsat_target(START, END).reindex(obs_raw["rate"].index),
                  "net": mvo.smoothed(obs_raw["net"])}
    rows, dom, written = [], {}, []
    for k in obs_raw:
        dom[f"coastsat_{k}_raw"], dom[f"coastsat_{k}_smoothed"] = obs_raw[k], obs_smooth[k]
    for preset in PRESETS:
        rd = run_dir(preset)
        if rd is None:
            print(f"skip  {preset}: not run yet")
            continue
        models = model_changes(rd)
        m_raw = {"rate": _gis(mvo.rates(rd)["lrr_m_yr"].reindex(range(1, N_DOMAINS + 1)).values),
                 "net": models[-1]}
        m_smooth = {k: mvo.smoothed(v) for k, v in m_raw.items()}
        for smoothed, obs, mod in ((False, obs_raw, m_raw), (True, obs_smooth, m_smooth)):
            written.append(figure(preset, rd, obs, mod, smoothed))
            for k in obs:
                rows.append(dict(run_name=rd.name, preset=preset, quantity=k,
                                 smoothed=smoothed, **skill(mod[k], obs[k])))
        for k in m_raw:
            dom[f"{preset}_{k}_raw"], dom[f"{preset}_{k}_smoothed"] = m_raw[k], m_smooth[k]
        written.append(gif(preset, rd, models, change))
    pd.DataFrame(rows).round(4).to_csv(OUT / "tables" / "scores.csv", index=False)
    pd.DataFrame(dom).round(4).to_csv(OUT / "tables" / "domain_series.csv")
    for p in written:
        print(f"wrote {Path(p).relative_to(_REPO)}")


if __name__ == "__main__":
    main()
