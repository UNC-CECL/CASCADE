"""
Every input a hindcast period reads, on one page: the island, its offset, the road, the target, the ends, the storms, sea level and the wave climate.

    python scripts/figure_making/pipeline/inputs_overview_figure.py

Reads each input through the same resolvers the runner uses (HATTERAS_PERIODS, hat_topo_version,
HAT_hindcast_config defaults); writes output/figures/3-model-inputs/inputs_overview_<start>_<end>.png
for the 1996 and 2010 periods. Details: scripts/figure_making/pipeline/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""
from __future__ import annotations

import functools
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402


# --- CONFIG ------------------------------------------------------------------
REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "figure_making" / "island"))
sys.path.insert(0, str(REPO / "scripts" / "hatteras_ms"))
from site_layer import hat_topo_version as tv  # noqa: E402
from site_layer import hat_env_forcings as env  # noqa: E402
from site_layer import hat_observed_rates as obs  # noqa: E402
from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_PERIODS, HATTERAS_BE_EDGE_ONLY, HATTERAS_COMMUNITY_ZONES, HATTERAS_ROAD_EVENTS,
    HATTERAS_NOURISHMENT_PROJECTS, island_offset_version,
)
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, C, C_1984, C_1997, INK, INK_MUTED, GRID_C, DOMAIN_AXIS_LABEL, figsize, figure_dir, save,
    record_caption, _title, open_frame, spines_for_image, town_bands,
)
from cascade_pipeline.coastsat_lowess import (  # noqa: E402
    CoastSatDataset, load_transect_data, spliced_lowess_series, DEFAULT_LOWESS,
)
from HAT_hindcast_config import field_default  # noqa: E402
import initialization_figures as init  # noqa: E402
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "3-env-forcings" / "3-storms"))
import storm_figures as sf  # noqa: E402
OUT = figure_dir("inputs")
# Literal names, so figure_index.py finds this script as their producer
FIGURE_NAMES = {1996: "inputs_overview_1996_2010.png", 2010: "inputs_overview_2010_2024.png"}
# The earlier period red, the later blue (house vintage pair)
PERIOD_COLOUR = {1996: C_1984, 2010: C_1997}
# Where each period's setback file comes from (4-management/road_setback_inputs)
SETBACK_SOURCE = {
    1996: "the 1984 measurement (the 1978-digitised line on the 1984-start extraction) with the 1989 "
          "Pea Island relocation applied",
    2010: "the 2004 measurement (the 2008-digitised line on the 2004-start extraction), unchanged because "
          "no relocation falls between 2004 and 2010",
}
TARGET_WINDOW = 7          # LOWESS width in domains, as the runner scores
# Hannah, 2026-10-01: the terrain ramp of the initial-island figures, not the classes
ELEVATION_SCHEME = "terrain"
# Hannah, 2026-10-01: name the major hurricanes, not the highest events by month
HURRICANES_NAMED = 6
# Hannah, 2026-10-01: Isabel (5.2 m) stretched the 1996-2010 axis; above the cap it is drawn at the edge
STORM_AXIS_CAP_M = {1996: 4.0}
GIS = np.arange(1, 91)
XLIM = (0.5, 90.5)
# -----------------------------------------------------------------------------


# A two-row CSV (domain ids, values) as a Series
def two_row(path):
    a = np.loadtxt(path, delimiter=",")
    return pd.Series(a[1], index=a[0].astype(int))


# The scoring target: LOWESS-7 averaged to domains, raw means at the south end
def observed_target(start, end):
    tag = f"{start}_{end}"
    csv = obs.COASTSAT_LRR_ROOT / tag / "transect_lrr_full.csv"
    ds = CoastSatDataset(label=tag, period_start=start, csv_path=str(csv))
    dom_ids, lrr, along = load_transect_data(ds)
    target, _frac = spliced_lowess_series(dom_ids, along, lrr, TARGET_WINDOW)
    return target, DEFAULT_LOWESS.skip_southern_domains


# One period's storm events at their peak hour, typed against HURDAT2 as the storm record figure types them
@functools.lru_cache(maxsize=1)
def _storm_record():
    forcing = sf.load_forcing()
    return sf.classify(sf.load_record(env.DEFAULT_STORM_VARIANT, forcing), forcing)


def storm_events(start, end):
    df = _storm_record()
    if (start, end) not in sf.PERIODS:
        raise ValueError(f"{start}-{end} is not a storm-record period {sf.PERIODS}")
    return df[df.pi == sf.PERIODS.index((start, end))].copy()


# The hurricanes to name: hurricane strength at closest approach, the highest runup, once per storm
def hurricanes_named(storms):
    h = storms[(storms.type == "tropical") & (storms.tc_status == "HU")]
    h = h.sort_values("rhigh_m", ascending=False).drop_duplicates(["tc_name", "calendar_year"])
    return h.head(HURRICANES_NAMED).sort_values("frac_year")


# The Duck record over the window, and this window's fit
def sea_level(start, end):
    rec = env.RSLR_RECORD_FILE
    lines = rec.read_text(encoding="utf-8").splitlines()
    first = next(i for i, ln in enumerate(lines) if ln.strip().startswith("Year"))
    d = pd.read_csv(rec, skiprows=first + 1, header=None, usecols=[0, 1, 2],
                    names=["year", "month", "msl"]).dropna()
    d["t"] = d.year + (d.month - 0.5) / 12
    d = d[d.t.between(start, end)]
    fits = pd.read_csv(env.RSLR_RATES_CSV)
    fit = fits[(fits.start_year == start) & (fits.end_year == end)].iloc[0]
    return d, fit


# Where a storm name may sit, tried in order: (dx, dy) in points from the dot, ha, va
_NEAR_DOT = ((0, 5, "center", "bottom"), (6, 0, "left", "center"), (-6, 0, "right", "center"),
             (5, 5, "left", "bottom"), (-5, 5, "right", "bottom"), (0, -6, "center", "top"),
             (5, -5, "left", "top"), (-5, -5, "right", "top"))


# Each storm name beside its own dot: the first spot clear of the other names, the dots and the axes edge
def label_near(ax, named, events, fontsize=7):
    from matplotlib.text import Text
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    box = ax.get_window_extent(rend)
    dots = ax.transData.transform(np.c_[events.frac_year, events.rhigh_m])
    r_dot = 2.6 * fig.dpi / 72
    placed = []
    # Highest first, so the biggest storms get the spot right above them
    for _, r in named.sort_values("rhigh_m", ascending=False).iterrows():
        own = ax.transData.transform((r.frac_year, r.rhigh_m))
        best = None
        for dx, dy, ha, va in _NEAR_DOT:
            ann = ax.annotate(sf.label_text(r), (r.frac_year, r.rhigh_m), xytext=(dx, dy),
                              textcoords="offset points", ha=ha, va=va, fontsize=fontsize)
            ann.update_positions(rend)
            bb = Text.get_window_extent(ann, rend).expanded(1.04, 1.1)
            ann.remove()
            inside = box.x0 <= bb.x0 and bb.x1 <= box.x1 and box.y0 <= bb.y0 and bb.y1 <= box.y1
            hits = sum(bb.overlaps(o) for o in placed)
            near = ((dots[:, 0] > bb.x0 - r_dot) & (dots[:, 0] < bb.x1 + r_dot)
                    & (dots[:, 1] > bb.y0 - r_dot) & (dots[:, 1] < bb.y1 + r_dot)
                    & (np.hypot(*(dots - own).T) > 1.0))
            cost = (0 if inside else 1e6) + 1e4 * hits + 50 * int(near.sum())
            if best is None or cost < best[0]:
                best = (cost, dx, dy, ha, va, bb)
            if cost == 0:
                break
        _, dx, dy, ha, va, bb = best
        placed.append(bb)
        ax.annotate(sf.label_text(r), (r.frac_year, r.rhigh_m), xytext=(dx, dy), textcoords="offset points",
                    ha=ha, va=va, fontsize=fontsize, color=INK, zorder=6,
                    path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])


# A labelled management span on the road panel
def mark_span(ax, gis, y, text, **patch):
    ax.axvspan(min(gis) - 0.5, max(gis) + 0.5, lw=0, zorder=0, **patch)
    ax.text((min(gis) + max(gis)) / 2, y, text, transform=ax.get_xaxis_transform(), ha="center",
            va="center", fontsize=6.8, color=INK, zorder=4, linespacing=1.0,
            bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", boxstyle="square,pad=0.1"))


# The page for one period
def fig_overview(start):
    period = HATTERAS_PERIODS[start]
    end = period["end_year"]
    colour = PERIOD_COLOUR[start]
    # The run's own offset file, checked against what the resolvers hand the figure
    offset_path = tv.INIT_ROOT / period["island_offset_file"]
    padded = np.loadtxt(offset_path, skiprows=1, delimiter=",")
    if np.abs(padded - init.load_offsets_m(start)).max() > 1e-6:
        raise RuntimeError(f"offset file the run reads ({offset_path}) differs from the figure's")
    offset = padded[init.START_REAL_INDEX:init.END_REAL_INDEX]
    setback = two_row(tv.INIT_ROOT / period["road_setback_file"]).reindex(GIS)
    target, skip = observed_target(start, end)
    ends = HATTERAS_BE_EDGE_ONLY[start]
    storms = storm_events(start, end)
    sl, fit = sea_level(start, end)
    canvas = init.detrended_canvas(start, include_buffers=False)
    _, product, topo_version, _rows = init.product_for(start)
    cmap, norm, cb_ticks, water = init.scheme_colours(ELEVATION_SCHEME)
    # The end year is a boundary, not a simulated year
    relocs = [e for e in HATTERAS_ROAD_EVENTS if hasattr(e, "displacement_m") and start <= e.year < end]
    bridges = [e for e in HATTERAS_ROAD_EVENTS if hasattr(e, "gis_domains") and start <= e.year < end]
    fills = [p for p in HATTERAS_NOURISHMENT_PROJECTS if start <= p.year < end]
    if bool(fills) != period["enable_nourishment"]:
        raise RuntimeError(f"{start}: nourishment projects {fills} disagree with enable_nourishment")
    waves = {k: field_default(k) for k in ("hs", "wave_period_s", "wave_asymmetry",
                                           "wave_angle_high_fraction", "relocation_setback_m")}

    # The storms get a full-width row of their own, so the hurricane names have room
    fig = plt.figure(figsize=figsize("double", height=9.4), constrained_layout=True)
    gs = fig.add_gridspec(7, 3, height_ratios=[0.95, 0.62, 0.78, 0.85, 0.4, 1.05, 1.15],
                          width_ratios=[1, 1, 0.9])
    ax_e = fig.add_subplot(gs[0, :])
    ax_o = fig.add_subplot(gs[1, :], sharex=ax_e)
    ax_r = fig.add_subplot(gs[2, :], sharex=ax_e)
    ax_t = fig.add_subplot(gs[3, :], sharex=ax_e)
    ax_b = fig.add_subplot(gs[4, :], sharex=ax_e)
    ax_s = fig.add_subplot(gs[5, :])
    ax_l = fig.add_subplot(gs[6, :2])
    ax_p = fig.add_subplot(gs[6, 2])

    # (a) the island at t=0, detrended, seaward at the bottom
    rows_km = canvas.shape[0] * init.CELL_SIZE_M / 1000.0
    im = ax_e.imshow(np.ma.masked_less(np.ma.masked_invalid(canvas), 0.0), cmap=cmap, norm=norm,
                     origin="lower", extent=(XLIM[0], XLIM[1], 0.0, rows_km),
                     interpolation="nearest", aspect="auto")
    ax_e.set_facecolor(water)
    ax_e.set_ylim(0, rows_km)
    ax_e.set_yticks([0, 1, 2])
    ax_e.set_ylabel("cross-shore\n(km)")
    spines_for_image(ax_e)
    cb = fig.colorbar(im, ax=ax_e, ticks=cb_ticks, pad=0.01, fraction=0.025, aspect=8)
    cb.set_label("m MHW", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    cb.outline.set_linewidth(0.5)
    _title(ax_e, 0, f"Initial island, {start}")

    # (b) the BRIE shoreline offset each domain starts at
    ax_o.plot(GIS, offset, color=INK, lw=1.3)
    ax_o.set_ylabel("offset\n(m, landward +)")
    open_frame(ax_o)
    _title(ax_o, 1, "BRIE shoreline offset")

    # (c) NC-12: setback, managed zones, and the events the window fires
    ax_r.step(GIS, setback, where="mid", color=C["ROAD"], lw=1.3)
    handles = [Line2D([], [], color=C["ROAD"], lw=1.3, label="NC-12 setback")]
    for e in relocs:
        mark_span(ax_r, e.displacement_m, 0.42, f"relocated\n{e.year}", color=C["ACCENT_FILL"], alpha=0.6)
    if relocs:
        handles.append(Patch(facecolor=C["ACCENT_FILL"], label="relocation"))
    for p in fills:
        mark_span(ax_r, p.gis_domains, 0.58, f"{p.name.split()[0]}\n{p.year}", color=C["ADDED_FILL"])
    if fills:
        handles.append(Patch(facecolor=C["ADDED_FILL"], label="nourishment"))
    for e in bridges:
        g = e.gis_domains
        ax_r.axvspan(min(g) - 0.5, max(g) + 0.5, facecolor="none", edgecolor=INK_MUTED, hatch="////",
                     lw=0, zorder=0)
        ax_r.text((min(g) + max(g)) / 2, 0.30, f"bridge\n{e.year}", transform=ax_r.get_xaxis_transform(),
                  ha="center", va="center", fontsize=6.8, color=INK, zorder=4, linespacing=1.0,
                  bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", boxstyle="square,pad=0.1"))
    if bridges:
        handles.append(Patch(facecolor="none", edgecolor=INK_MUTED, hatch="////", label="bridge (road removed)"))
    for lo, hi in HATTERAS_COMMUNITY_ZONES:
        ax_r.axvspan(lo - 0.5, hi + 0.5, ymax=0.08, color=C["BASE"], lw=0, zorder=1)
    handles.append(Patch(facecolor=C["BASE"], label="managed beach/dune"))
    ax_r.set_ylabel("setback\n(m)")
    ax_r.set_ylim(0, np.nanmax(setback) * 1.4)
    ax_r.grid(axis="y", color=GRID_C, lw=0.5)
    open_frame(ax_r)
    ax_r.legend(handles=handles, loc="upper left", frameon=False, fontsize=7, ncol=len(handles),
                handlelength=1.4, columnspacing=1.0)
    _title(ax_r, 2, "Road and management")

    # (d) the observed target the run is scored against
    # One line through the whole target: raw means at the south end run into the LOWESS
    tn, ts_ = target[target.index > skip], target[target.index <= skip]
    ax_t.axhline(0, color=INK_MUTED, lw=0.5)
    ax_t.plot(tn.index, tn.values, "o-", color=colour, lw=1.3, ms=2.2, zorder=3,
              label=f"LOWESS ({TARGET_WINDOW} domains)")
    joined = target.loc[:tn.index.min()]       # the raw means and the first LOWESS domain
    ax_t.plot(joined.index, joined.values, "-", color=C["ACCENT"], lw=1.3, zorder=2)
    ax_t.plot(ts_.index, ts_.values, "s-", color=C["ACCENT"], lw=1.3, ms=2.2, zorder=4,
              label="raw domain mean")
    ax_t.set_ylabel("LRR\n(m/yr, + accretion)")
    ax_t.grid(axis="y", color=GRID_C, lw=0.5)
    open_frame(ax_t)
    ax_t.legend(frameon=False, fontsize=7, loc="upper right", ncol=2)
    _title(ax_t, 3, f"Observed target, CoastSat {start}–{end}")

    # (e) the source/sink the edgeBE arm carries: the two ends, nothing between
    ax_b.bar([1, 90], ends, width=1.6, color=C["REF"])
    for g, v in zip((1, 90), ends):
        ax_b.text(g + (1.6 if g == 1 else -1.6), v / 2, f"{v:+.1f} m/yr", ha="left" if g == 1 else "right",
                  va="center", fontsize=7, color=INK)
    ax_b.axhline(0, color=INK_MUTED, lw=0.5)
    ax_b.set_ylabel("source/sink\n(m/yr)")
    ax_b.set_ylim(min(0, min(ends)) * 1.15, max(0, max(ends)) * 1.15)
    ax_b.set_xlim(*XLIM)
    ax_b.set_xticks([1, 15, 30, 45, 60, 75, 90])
    ax_b.set_xlabel(DOMAIN_AXIS_LABEL)
    open_frame(ax_b)
    _title(ax_b, 4, "End source/sink (edgeBE runs)")
    for ax in (ax_e, ax_o, ax_r, ax_t):
        ax.tick_params(labelbottom=False)
    for ax in (ax_o, ax_r, ax_t, ax_b):
        town_bands(ax, label=ax is ax_o)

    # (f) every storm event at its peak hour, coloured by type; the biggest hurricanes named
    # A capped axis: anything above the cap is a triangle at the top edge, named with its height
    cap = STORM_AXIS_CAP_M.get(start)
    shown = storms if cap is None else storms[storms.rhigh_m <= cap]
    above = storms.iloc[0:0] if cap is None else storms[storms.rhigh_m > cap]
    for k, size in (("other", 11), ("tropical", 15)):
        d = shown[shown.type == k]
        ax_s.scatter(d.frac_year, d.rhigh_m, s=size, color=sf.TYPE_COLOURS[k], alpha=0.9,
                     edgecolors="white", linewidths=0.35, zorder=3)
    named = hurricanes_named(storms)
    named = named[named.index.isin(shown.index)]
    ax_s.set_xlim(start, end)
    ax_s.set_ylim(1.0, cap if cap is not None else storms.rhigh_m.max() + 0.7)
    for _, r in above.sort_values("rhigh_m", ascending=False).drop_duplicates(["tc_name", "calendar_year"]).iterrows():
        ax_s.scatter([r.frac_year], [cap], marker="^", s=30, color=sf.TYPE_COLOURS[r.type],
                     edgecolors="white", linewidths=0.35, zorder=5, clip_on=False)
        name = sf.label_text(r) if r.type == "tropical" else r.peak.strftime("%b %Y")
        ax_s.annotate(f"{name}, {r.rhigh_m:.1f} m", (r.frac_year, cap), xytext=(6, -2),
                      textcoords="offset points", ha="left", va="top", fontsize=7, color=INK, zorder=6,
                      path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])
    ax_s.set_ylabel("R$_{high}$ (m above MHW)")
    ax_s.set_xlabel("year")
    ax_s.grid(axis="y", color=GRID_C, lw=0.5)
    ax_s.set_axisbelow(True)
    open_frame(ax_s)
    ax_s.legend(handles=[Line2D([], [], marker="o", ls="", ms=4, color=sf.TYPE_COLOURS["tropical"],
                                label="tropical cyclone"),
                         Line2D([], [], marker="o", ls="", ms=3.5, color=sf.TYPE_COLOURS["other"],
                                label="other (mostly nor'easters)")],
                loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, frameon=False, fontsize=7,
                handletextpad=0.1, borderaxespad=0.2)
    _title(ax_s, 5, "Storms")
    label_near(ax_s, named, storms)

    # (g) sea level: the record and the rate the model is given
    ann = sl.groupby("year")["msl"].mean()
    ax_l.plot(sl.t, sl.msl, color=C["BASE_FILL"], lw=0.6)
    ax_l.plot(ann.index + 0.5, ann.values, color=INK, lw=1.0, label="Duck, annual mean")
    tt = np.array([start, end], float)
    rate = period["sea_level_rise_rate"]
    mid = fit.slope_m_yr * tt.mean() + fit.intercept_m
    ax_l.plot(tt, mid + rate * (tt - tt.mean()), color=colour, lw=1.6,
              label=f"model: {rate * 1000:.0f} mm/yr")
    ax_l.set_xlim(start, end)
    ax_l.set_ylabel("mean sea level (m)")
    ax_l.set_xlabel("year")
    ax_l.grid(axis="y", color=GRID_C, lw=0.5)
    open_frame(ax_l)
    ax_l.legend(frameon=False, fontsize=7, loc="upper left")
    _title(ax_l, 6, "Sea level")

    # (h) the scalar settings, as a table
    rowsx = [("wave height H$_s$", f"{waves['hs']:g} m"),
             ("wave period T$_p$", f"{waves['wave_period_s']:g} s"),
             ("wave asymmetry", f"{waves['wave_asymmetry']:g}"),
             ("high-angle fraction", f"{waves['wave_angle_high_fraction']:g}"),
             ("sea-level rise", f"{rate * 1000:.0f} mm/yr"),
             ("relocation setback", f"{waves['relocation_setback_m']:g} m"),
             ("nourishment", "none" if not fills else f"{len(fills)} projects"),
             ("domains", "90 + 2 × 15"),
             ("years simulated", f"{start}–{end - 1}")]
    ax_p.set_axis_off()
    for k, (name, val) in enumerate(rowsx):
        y = 0.95 - k * 0.105
        ax_p.text(0.0, y, name, transform=ax_p.transAxes, fontsize=7.3, color=INK_MUTED, va="top")
        ax_p.text(1.0, y, val, transform=ax_p.transAxes, fontsize=7.3, color=INK, va="top", ha="right")
        ax_p.plot([0, 1], [y - 0.085] * 2, transform=ax_p.transAxes, color=GRID_C, lw=0.5)
    _title(ax_p, 7, "Settings")

    out = save(fig, OUT / FIGURE_NAMES[start])
    plt.close(fig)
    zones = ", ".join(f"GIS {lo}–{hi}" for lo, hi in HATTERAS_COMMUNITY_ZONES)
    events = [f"purple, the road relocation of {e.year} at GIS {min(e.displacement_m)}–{max(e.displacement_m)}, "
              f"{min(e.displacement_m.values()):.0f}–{max(e.displacement_m.values()):.0f} m landward"
              for e in relocs]
    events += [f"orange, the {p.name} of {p.year} over GIS {min(p.gis_domains)}–{max(p.gis_domains)} "
               f"({p.volume_cubic_yards / 1e6:.1f} million cubic yards reported)" for p in fills]
    events += [f"hatched, the bridge of {e.year} that takes the road off GIS {min(e.gis_domains)}–"
               f"{max(e.gis_domains)}" for e in bridges]
    no_fill = "" if fills else f" No nourishment falls in {start}–{end - 1}."
    record_caption(out[0],
        f"Every input the {start}–{end} hindcast reads, on one page; each panel is drawn from the file the "
        "runner itself resolves (HATTERAS_PERIODS, hat_topo_version, the HAT_hindcast_config defaults), and the "
        "step figures in the numbered folders beside this one show how each is built. "
        f"(a) The island at t=0: every domain's elevation array from the {product} {topo_version} extraction, "
        f"{init.SCHEME_NOTE[ELEVATION_SCHEME]}; detrended, each domain on its own frame origin, "
        "seaward at the bottom, cross-shore exaggerated against alongshore, domains in GIS order south to north "
        "with each one's 50 cells in the same direction (1-domains/initial_island/). "
        f"(b) The BRIE shoreline offset each domain starts at, metres landward of the most seaward domain "
        f"(offset build {island_offset_version(start)}, from the {tv.dune_line_for_year(start)} dune line; the "
        "15 buffer domains each end are not drawn; 2-brie-offset/). "
        f"(c) The NC-12 setback, metres landward of interior row 0: {SETBACK_SOURCE[start]} "
        f"(4-management/road_setback_inputs). Against it, the management the window fires: "
        f"{'; '.join(events)}; grey along the axis, the community zones where the beach/dune manager runs "
        f"({zones}).{no_fill} "
        f"(d) The target the run is scored against: CoastSat linear regression rates {start}–{end}, a "
        f"{TARGET_WINDOW}-domain LOWESS averaged to domains for GIS {skip + 1}–90 (circles) and the raw domain "
        f"means for GIS 1–{skip} (purple squares and line), joined into the one series the runner scores "
        "(5-observed-target/). "
        f"(e) The source/sink rates the edgeBE arm carries, GIS 1 {ends[0]:+.4f} and GIS 90 {ends[1]:+.4f} m/yr, "
        "solved so the ends match the target; every other domain is 0, and the zeroBE arm carries none "
        "(7-source-sink/be_end_solve). "
        f"(f) Every storm event ({len(storms)}) in the series the run reads ({Path(period['storm_file']).name}), "
        "at the hour of its peak total water level, by peak runup R_high above MHW. Orange: a tropical or "
        f"subtropical cyclone in the NHC best tracks (HURDAT2) was within {sf.TC_NEAR_KM} km of Cape Hatteras "
        f"within {sf.TC_WINDOW_H} h of the peak ({int((storms.type == 'tropical').sum())} events); grey: every "
        "other event, mostly nor'easters, identified only by the absence of a nearby cyclone. Named: the "
        f"{len(named)} hurricanes (hurricane strength at their closest fix) that raised the highest water, "
        f"each once ({', '.join(sf.label_text(r) for _, r in named.iterrows())}); the storm record figure "
        "names every tropical cyclone (3-env-forcings/3-storms/figures/). "
        + ("" if above.empty else
           f"The axis stops at {cap:g} m so the other events are not squeezed; "
           + "; ".join(f"{sf.label_text(r) if r.type == 'tropical' else r.peak.strftime('%B %Y')} "
                       f"({r.rhigh_m:.2f} m)" for _, r in above.iterrows())
           + " lies above it, drawn as a triangle at the top edge. ")
        + "(g) Monthly (grey) and annual (black) mean sea level at Duck (NOAA 8651370) and the rate the model is "
        f"given, {rate:.3f} m/yr, the {fit.slope_mm_yr:.2f} ± {fit.ci95_mm_yr:.2f} mm/yr fit over the window "
        "rounded to 0.001, drawn through the fit's mid-window level (3-forcing/forcing_timeline_1996). "
        "(h) The scalar settings: the wave climate (option A, adopted 2026-09-27), the sea-level rate, where a "
        "relocated road is placed, and the domain layout. The Buxton groin is a separate on/off arm of every "
        "scenario and is not drawn. "
        "GIS 1 is Cape Point, GIS 90 Pea Island.")
    return out


# Run: one page per period
def main():
    apply_style()
    for start in FIGURE_NAMES:
        print(fig_overview(start)[0].relative_to(REPO))


if __name__ == "__main__":
    main()
