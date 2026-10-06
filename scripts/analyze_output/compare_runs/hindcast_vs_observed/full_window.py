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

import dataclasses
import io
import json
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
    C_1984_FILL, C_1997_FILL, COMPARISONS_ROOT, DOMAIN_AXIS_LABEL, INK, INK_MUTED, _title, apply_style,
    compare_header, figsize, mark_offaxis, offaxis_clause, open_frame,
    record_caption, save, structures, town_bands)
from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT  # noqa: E402
from site_layer.hat_topo_version import INIT_ROOT  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS  # noqa: E402

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
# The six alongshore sections of the zoomed positions figure, GIS inclusive
# The 15-domain sections of the GIF-zoom position figures, GIS inclusive
SECTIONS = [(1, 15), (16, 30), (31, 45), (46, 60), (61, 75), (76, 90)]
# Fixed y axes for the runner-style figures: rate in m/yr, observed position change in m (None fits the data)
RATE_YLIM = (-5.0, 5.0)
POSITION_YLIM = (-100.0, 100.0)
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

# Start, observed end and modelled end in the model frame, landward up as in the runner's position GIFs
def positions_figure(preset, rd, obs_net_raw):
    m = np.load(next(rd.glob("*_shoreline_matrix.npy")))
    D = HATTERAS_DOMAINS
    # Landward positive, so a plain axis puts landward up and the ocean at the bottom
    start = _gis(m[0][D.start_real_index:D.end_real_index])
    model_end = _gis(m[-1][D.start_real_index:D.end_real_index])
    obs_end = start - obs_net_raw
    s = skill(start - model_end, obs_net_raw)
    fig, ax = plt.subplots(figsize=figsize("double", height=3.8), layout="constrained")
    ax.plot(start.index, start.values, color=INK_MUTED, lw=2.2, zorder=3)
    ax.plot(obs_end.index, obs_end.values, color=C_OBS, lw=1.3, zorder=5, marker="o", ms=2.2)
    ax.plot(model_end.index, model_end.values, color=C_MODEL, lw=1.2, zorder=6)
    ax.set_xlim(1, N_DOMAINS)
    # Headroom above the planform, so the village labels clear the lines
    lo = float(np.nanmin([start.min(), obs_end.min(), model_end.min()]))
    hi = float(np.nanmax([start.max(), obs_end.max(), model_end.max()]))
    ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.18 * (hi - lo))
    ax.grid(axis="y")
    open_frame(ax)
    town_bands(ax, label=True)
    ax.set_ylabel("Cross-shore position,\nmodel frame (m), landward ▲")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    fig.legend(handles=[
        Line2D([], [], color=INK_MUTED, lw=2.2, label="Start position, 1996 (model year 0)"),
        Line2D([], [], color=C_OBS, lw=1.3, marker="o", ms=2.2,
               label="CoastSat 2025 (calendar-year mean)"),
        Line2D([], [], color=C_MODEL, lw=1.2,
               label=f"Modelled 1 Jan 2026: bias {s['bias']:+.1f} m, RMSE {s['rmse']:.1f}, "
                     f"r {s['r']:.2f}")],
        loc="outside lower center", ncol=3, frameon=False, fontsize=7)
    compare_header(fig, f"Full management, no groin, {preset}: the 1996 start, CoastSat 2025 "
                        "and the modelled end shoreline")
    # Last, once the legend and header have fixed the layout
    structures(ax, label=True)
    caption = (
        f"Full management, no groin, no relocation; {PRESET_TEXT[preset]}. Run {rd.name}. "
        "Shoreline position in the model's cross-shore frame, landward up and the ocean at the "
        "bottom, as in the runner's position GIFs: grey the 1996 start (the shoreline offset "
        f"built on the CoastSat mean over {START_MEAN.replace('_', ' to ')}), blue CoastSat "
        "2025 (the start moved by the calendar-2025 mean minus the start mean, per domain), "
        "dark the modelled shoreline on 1 Jan 2026. Domain values, unsmoothed. The CoastSat "
        "mean is centred about half a year before the model end. Scores: modelled minus "
        "observed change, seaward positive, over the interior GIS 2-89.")
    pngs = [rd / "figures" / "start_and_end_positions_1996_2025.png",
            OUT / "positions" / f"full_window_{preset}_full_management_{WINDOW}_start_and_end_positions.png"]
    for png in pngs:
        save(fig, png, dpi=300, close=False)
        record_caption(png, caption)
    plt.close(fig)
    return pngs


# One GIF-style frame per 15-domain section: the start, the modelled end with its change shaded, CoastSat 2025 on top
def position_section_figures(preset, rd, obs_net_raw):
    m = np.load(next(rd.glob("*_shoreline_matrix.npy")))
    D = HATTERAS_DOMAINS
    # Landward positive, as the runner's position GIFs draw it (ocean at the bottom)
    start = _gis(m[0][D.start_real_index:D.end_real_index])
    model_end = _gis(m[-1][D.start_real_index:D.end_real_index])
    obs_end = start - obs_net_raw
    pngs = []
    for g0, g1 in SECTIONS:
        ref = start.loc[g0:g1].mean()
        s, mo, ob = (v.loc[g0:g1] - ref for v in (start, model_end, obs_end))
        fig, ax = plt.subplots(figsize=figsize("double", aspect=0.56), layout="constrained")
        x = s.index.to_numpy(float)
        ax.fill_between(x, mo, s, where=(mo <= s), interpolate=True, color=C_1997_FILL,
                        alpha=0.7, lw=0, zorder=1)
        ax.fill_between(x, mo, s, where=(mo > s), interpolate=True, color=C_1984_FILL,
                        alpha=0.7, lw=0, zorder=1)
        ax.plot(x, s.values, color=INK_MUTED, ls=(0, (4, 3)), lw=0.9, zorder=3)
        ax.plot(x, ob.values, color=C_OBS, lw=2.4, zorder=4, marker="o", ms=3.0)
        ax.plot(x, mo.values, color=INK, lw=1.6, zorder=5)
        vals = np.concatenate([s.values, mo.values, ob.values])
        pad = (np.nanmax(vals) - np.nanmin(vals)) * 0.10
        # Values are landward positive, so a plain axis puts landward up and the ocean at the bottom
        ax.set_ylim(np.nanmin(vals) - pad, np.nanmax(vals) + pad)
        ax.set_xlim(g0 - 0.5, g1 + 0.5)
        ax.grid(axis="y")
        open_frame(ax)
        town_bands(ax, label=True)
        ax.set_ylabel("Cross-shore position (m, rel. section start mean)\nlandward ▲")
        ax.set_xlabel(DOMAIN_AXIS_LABEL)
        ax.set_title(f"GIS {g0}-{g1}", loc="left", fontsize=10)
        fig.legend(handles=[
            Line2D([], [], color=INK_MUTED, ls=(0, (4, 3)), lw=0.9, label="1996 start"),
            Line2D([], [], color=INK, lw=1.6, label="Modelled 1 Jan 2026"),
            plt.Rectangle((0, 0), 1, 1, color=C_1997_FILL, alpha=0.7, lw=0,
                          label="model accretion, seaward of 1996"),
            plt.Rectangle((0, 0), 1, 1, color=C_1984_FILL, alpha=0.7, lw=0,
                          label="model erosion, landward of 1996"),
            Line2D([], [], color=C_OBS, lw=2.4, marker="o", ms=3.0, label="CoastSat 2025")],
            loc="outside lower center", ncol=3, frameon=False, fontsize=7)
        compare_header(fig, f"Full management, no groin, {preset}: 1996 start, CoastSat 2025 and "
                            "the modelled end shoreline")
        structures(ax, label=True)
        png = rd / "figures" / "position_sections" / f"start_and_end_positions_1996_2025_gis_{g0:02d}-{g1:02d}.png"
        save(fig, png, dpi=300, close=True)
        record_caption(png, (
            f"Full management, no groin, no relocation; {PRESET_TEXT[preset]}. Run {rd.name}. "
            f"GIS {g0}-{g1} at the zoom of the runner's position GIFs: cross-shore position "
            "relative to the section's mean 1996 position, landward up and the ocean at the "
            "bottom. Dashed the 1996 start (the shoreline offset built on the CoastSat mean over "
            f"{START_MEAN.replace('_', ' to ')}); black the modelled shoreline on 1 Jan 2026, "
            "shaded blue where it lies seaward of the start and red where landward; light blue "
            "CoastSat 2025 (the start moved by the calendar-2025 mean minus the start mean). "
            "Domain values, unsmoothed, true scale."))
        pngs.append(png)
    return pngs

# Rename a runner rate figure's rate wording for a position change in metres
def _relabel_as_position(fig, ax):
    ax.set_ylabel("Shoreline Position Change (m)")
    swaps = (("CoastSat LRR per 500 m domain",
              "CoastSat observed change per 500 m domain (calendar-2025 mean minus the 1996 start mean)"),
             ("transect LRR", "transect change"))
    texts = list(fig.texts) + list(ax.texts) + [ax.title]
    for leg in [ax.get_legend(), *fig.legends]:
        if leg is not None:
            texts += list(leg.get_texts())
    for t in texts:
        new = t.get_text()
        for a, b in swaps:
            new = new.replace(a, b)
        if new != t.get_text():
            t.set_text(new)

# The model line smoothed as the observation is (7-domain LOWESS, GIS 1-10 raw), buffers left as they are
def _smooth_model(values):
    D = HATTERAS_DOMAINS
    out = np.asarray(values, float).copy()
    out[D.start_real_index:D.end_real_index] = mvo.smoothed(
        _gis(out[D.start_real_index:D.end_real_index])).to_numpy(float)
    return out


# Draw one runner figure with the smoothed model, put the raw model faint behind it, rename, save
def _draw_smoothed_model(draw, raw, png, position):
    D = HATTERAS_DOMAINS
    fig, ax = draw(_smooth_model(raw))[:2]
    gis = np.arange(D.first_gis_id, D.last_gis_id + 1)
    ax.plot(gis, np.asarray(raw, float)[D.start_real_index:D.end_real_index],
            color=HATTERAS_ANNOTATIONS.model_color, lw=0.9, alpha=0.45, zorder=5.5)
    if position:
        _relabel_as_position(fig, ax)
    # The legend's model entry names both lines
    for leg in [ax.get_legend(), *fig.legends]:
        if leg is None:
            continue
        for t in leg.get_texts():
            if t.get_text().startswith("CASCADE"):
                t.set_text(t.get_text() + " (7-domain LOWESS; faint: unsmoothed)")
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# The runner's own rate figures redrawn at a fixed y axis, and the same style for the observed position change
def runner_style_figures(preset, rd, change):
    from cascade_pipeline.coastsat_lowess import (CoastSatDataset, LowessConfig,
                                                  build_coastsat_series)
    from cascade_pipeline.plotting.rate_comparison import (
        DEFAULT_RATE_COMPARISON, plot_annotated_rate_comparison, plot_rate_comparison)
    from cascade_pipeline.run_info import RunInfo
    from cascade_pipeline.run_layout import resolve
    from cascade_pipeline.shoreline import compute_change_rate, compute_lrr
    run_name = rd.name
    meta = json.loads((rd / f"{run_name}_run_metadata.json").read_text(encoding="utf-8"))
    wave, src = meta.get("wave climate", {}), meta.get("source/sink", {})
    run = RunInfo(run_name=run_name, run_dir=str(rd), start_year=START, end_year=END,
                  Hs=wave.get("wave_height_m"), flip_sign_model=True,
                  background_erosion_on=bool(src.get("background_erosion_on", True)),
                  wave_climate=meta.get("scenario", {}).get("wave climate"))
    m = np.load(next(rd.glob("*_shoreline_matrix.npy")))
    years = m.shape[0] - 1
    # As the runner: 7-domain LOWESS, GIS 1-10 raw
    lowess = LowessConfig(window_domains=(7,), skip_southern_domains=10)
    kw = dict(domains=HATTERAS_DOMAINS, annotations=HATTERAS_ANNOTATIONS, lowess_config=lowess)

    rate, _ = compute_lrr(m, span_years=years, flip_sign=True)
    cs_rate = build_coastsat_series(
        [CoastSatDataset(label=f"CoastSat LRR ({START}-{END})", period_start=START,
                         csv_path=str(COASTSAT_LRR_ROOT / WINDOW / "transect_lrr_full.csv"))],
        active_period_start=START, lowess_config=lowess, domains=HATTERAS_DOMAINS)
    cfg = dataclasses.replace(DEFAULT_RATE_COMPARISON, ylim=RATE_YLIM, ylim_real=RATE_YLIM)
    pngs = [resolve(rd, "figure_rate", run_name), resolve(rd, "figure_rate_buffers", run_name)]
    _draw_smoothed_model(lambda v: plot_rate_comparison(
        v, cs_rate, run, real_domains_only=True, estimator="lrr", show=False, config=cfg, **kw),
        rate, pngs[0], position=False)
    _draw_smoothed_model(lambda v: plot_annotated_rate_comparison(
        v, cs_rate, run, estimator="lrr", show=False, config=cfg, **kw),
        rate, pngs[1], position=False)

    # Observed change per transect (no rate): the calendar-2025 mean minus the DEM-centred start mean
    obs = change[["domain_number", OBS_END_YEAR]].rename(columns={OBS_END_YEAR: "change_m"})
    obs_csv = OUT / "tables" / "transect_observed_change_1996_start_to_2025.csv"
    obs.rename_axis("transect_id").reset_index().round(4).to_csv(obs_csv, index=False)
    cs_pos = build_coastsat_series(
        [CoastSatDataset(label="CoastSat observed change",
                         period_start=START, csv_path=str(obs_csv), rate_col="change_m")],
        active_period_start=START, lowess_config=lowess, domains=HATTERAS_DOMAINS)
    position = compute_change_rate(m, span_years=1, flip_sign=True)
    pcfg = dataclasses.replace(
        DEFAULT_RATE_COMPARISON, quantity="position",
        observed_label="CoastSat observed change (2025 mean minus 1996 start mean, 7-domain LOWESS)",
        observed_description=(f"observed change, the calendar-2025 CoastSat mean minus the mean over "
                              f"{START_MEAN.replace('_', ' to ')} (no rate)"),
        ylim=POSITION_YLIM, ylim_real=POSITION_YLIM)
    pdir = rd / "figures" / "position_change" / "observed"
    pdir.mkdir(parents=True, exist_ok=True)
    ppngs = [pdir / "shoreline_position_change.png",
             pdir / "shoreline_position_change_with_buffers.png"]
    _draw_smoothed_model(lambda v: plot_rate_comparison(
        v, cs_pos, run, real_domains_only=True, estimator=None, show=False, config=pcfg, **kw),
        position, ppngs[0], position=True)
    _draw_smoothed_model(lambda v: plot_annotated_rate_comparison(
        v, cs_pos, run, estimator=None, show=False, config=pcfg, **kw),
        position, ppngs[1], position=True)
    return pngs + ppngs


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
        written.extend(positions_figure(preset, rd, obs_raw["net"]))
        written.extend(runner_style_figures(preset, rd, change))
        written.extend(position_section_figures(preset, rd, obs_raw["net"]))
    pd.DataFrame(rows).round(4).to_csv(OUT / "tables" / "scores.csv", index=False)
    pd.DataFrame(dom).round(4).to_csv(OUT / "tables" / "domain_series.csv")
    for p in written:
        print(f"wrote {Path(p).relative_to(_REPO)}")


if __name__ == "__main__":
    main()
