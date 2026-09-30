"""
model_mechanics_figures.py
==============================================================================
How the models work, drawn from a finished Hatteras run rather than from a
cartoon: the Barrier3D grid through one storm year, one domain's cross-shore
budget over a window, BRIE's alongshore diffusion, how CASCADE splits a
shoreline's change between the two, the whole island as the coupled model
holds it, the annual coupling loop with the unit handed across each join, and
what the two management modules do to the grid.

    python scripts/figure_making/model/model_mechanics_figures.py [--only NAME]

Writes to output/figures/4-model-mechanics/<model>/ (PNG at the top, PDF + CAPTIONS.md under
supporting/):

    barrier3d_storm_year.png     one domain, one storm year: the year's storms
                                 against the dune, the grid before and after,
                                 and what moved
    barrier3d_domain_budget.png  one domain over 1996-2010: the profile at the
                                 two ends, the shoreface toe / shoreline /
                                 back-barrier, and the annual fluxes
    brie_diffusion.png           BRIE's angle-dependent diffusivity and where
                                 the Hatteras shoreline sits on it
    brie_domain_order.png        BRIE alone on the 1996 offset with the
                                 domains fed south->north (as run) and
                                 reversed: the order matters through the
                                 wave asymmetry
    brie_domain_orientation.png  north-up map of the numbered domains, the
                                 wave asymmetry and net drift, beside BRIE's
                                 array (index <-> GIS) in the same orientation
    brie_asymmetry_explained.png what the asymmetry counts, why its waves are
                                 head-on to positive-θ links, what that is on
                                 Hatteras, and what it does to a cape (smoothing
                                 rate, not drift; the step does not conserve sand)
    cascade_shoreline_split.png  each domain's shoreline change split into the
                                 Barrier3D cross-shore part, the source/sink
                                 (BE) part and the BRIE alongshore part
    cascade_island_grids.png     all 90 Barrier3D grids placed on the BRIE
                                 shoreline, the island as CASCADE holds it
    cascade_coupling_loop.png    the annual loop: who runs, what is handed
                                 across, and the unit it is handed in
    management_modules.png       the roadway manager (overwash cleared, dunes
                                 rebuilt) and a nourishment spreading
                                 alongshore

THE RUNS
    Everything is read from the saved Cascade object in each run's .npz
    (output/raw_runs/matrix/...), so a figure and a run cannot disagree. The
    natural run (no road, no beach/dune manager, edgeBE, 1996-2010) carries
    the physics figures because nothing human touches its grids; the
    management figure pairs a managed run with the natural run of the same
    window. The unit contract behind the loop figure is UNITS.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-29
"""

from __future__ import annotations

import argparse
import functools
import sys
from pathlib import Path

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))

from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, C, C_1984, C_1997, INK, INK_MUTED, CELL_M, DOMAIN_AXIS_LABEL,
    figsize, figure_dir, save, record_caption, _title, open_frame,
    elevation_cmap, town_bands,
)
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent))
from storm_replay import replay  # noqa: E402

OUT = figure_dir("mechanics")   # one sub-folder per model: barrier3d/, brie/, cascade/, storm_routing/
MATRIX = REPO / "output" / "raw_runs" / "matrix"
DAM = 10.0                      # Barrier3D decametre -> metre

RUN_NATURAL = ("1996_2010", "edgeBE", "HAT_1996_2010_edgeBE_offsetmetres_noroad_nobdm_nogroin")
RUN_ROAD = ("1996_2010", "edgeBE", "HAT_1996_2010_edgeBE_offsetmetres_road_nobdm_nogroin")
# the fill pair differs ONLY in the fills: same road, same beach/dune manager
RUN_NOURISH = ("2010_2024", "edgeBE", "HAT_2010_2024_edgeBE_offsetmetres_road_bdm_nourish_nogroin")
RUN_NOURISH_OFF = ("2010_2024", "edgeBE", "HAT_2010_2024_edgeBE_offsetmetres_road_bdm_nonourish_nogroin")

STORM_GIS, STORM_T = 6, 11      # 2006 storms: an intact dune (mean crest 2.9 m) that 40% of storms overtop
BUDGET_GIS = 45                 # the project's example domain (site/domain_schematic)
ROAD_GIS = 13                   # most overwash cleared off NC-12 in the road run
NOURISH_GIS, NOURISH_YEAR = 86, 2014   # Rodanthe emergency fill, GIS 84-89


# =============================================================================
# LOADING
# =============================================================================

@functools.lru_cache(maxsize=4)
def load_run(spec):
    window, preset, name = spec
    path = MATRIX / window / preset / name / f"{name}.npz"
    return np.load(path, allow_pickle=True)["cascade"][0]


def start_year(spec):
    return int(spec[0].split("_")[0])


def pad(gis):
    return DOM.gis_to_pad(gis)


def real_pads():
    return np.arange(pad(1), pad(90) + 1)


def plan_grid(b3d, t):
    """(cross-shore rows, alongshore columns) in m MHW: the two dune rows on
    top of the berm, then the interior. Ocean first, as Barrier3D stores it."""
    dune = (b3d.DuneDomain[t].T + b3d.BermEl) * DAM            # (2, 50)
    return np.vstack([dune, np.asarray(b3d.DomainTS[t]) * DAM])


def land_rows(*grids, margin=8):
    last = max(int(np.nonzero((g > 0).any(axis=1))[0].max()) for g in grids)
    return last + 1 + margin


def beach_width_m(cascade):
    return float(cascade._initial_beach_width[0])


# =============================================================================
# 1. BARRIER3D: ONE STORM YEAR
# =============================================================================

def fig_barrier3d_storm_year():
    spec = RUN_NATURAL
    c = load_run(spec)
    b = c.barrier3d[pad(STORM_GIS)]
    t = STORM_T
    year = start_year(spec) + t - 1          # update t applies storms with time == t
    berm = b.BermEl * DAM
    storms = b.StormSeries[b.StormSeries[:, 0] == t]
    # The storms as the model met them: replayed through update() from the
    # saved grid, which reproduces the run exactly (storm_routing_check). The
    # crest they are tested against is the dune AFTER the year's growth,
    # computed once before the first storm (storm_replay.py).
    got, _ = replay(b, t, [(r[1] * DAM, r[2] * DAM, r[3], int(r[4])) for r in storms])
    crest = got[0]["crest_pre"]
    crest_mean, crest_min = crest.mean(), crest.min()
    crest_jan = (b.DuneDomain[t - 1].max(axis=1) + b.BermEl) * DAM + b._RSLR[t] * DAM
    over = np.array([bool(g["gaps"]) for g in got])

    pre, post = plan_grid(b, t - 1), plan_grid(b, t)
    n = min(len(pre), len(post))
    moved = np.nonzero((np.abs(post[:n] - pre[:n]) > 0.05).any(axis=1))[0]
    rows = min(n, max(40, int(moved.max()) + 15), land_rows(pre[:n], post[:n]))
    pre, post = pre[:rows], post[:rows]
    diff = post - pre

    fig = plt.figure(figsize=figsize("double", height=6.6), constrained_layout=True)
    gs = fig.add_gridspec(4, 2, width_ratios=[1, 0.025], height_ratios=[1.05, 1, 1, 1])
    ax_s = fig.add_subplot(gs[0, 0])
    x = np.arange(len(got))
    for i, g in enumerate(got):
        ax_s.vlines(i, g["rlow"], g["rhigh"], color=C["ACCENT"] if over[i] else C["BASE"], lw=3.2)
        if over[i]:
            ax_s.text(i, g["rhigh"] + 0.08, f"{g['owloss']:.0f} m$^3$/m", ha="center", va="bottom", fontsize=7)
    ax_s.axhline(berm, color=INK_MUTED, lw=0.8, ls=(0, (3, 2)))
    ax_s.axhspan(crest_min, crest_mean, color=C["ADDED_FILL"], lw=0, zorder=0)
    ax_s.axhline(crest_mean, color=C["ADDED"], lw=1.0)
    ax_s.axhline(crest_min, color=C["ADDED"], lw=0.8, ls=(0, (1, 1.5)))
    ax_s.axhline(crest_jan.min(), color=INK_MUTED, lw=0.6, ls=(0, (1, 2)))
    ax_s.set_xticks(x)
    ax_s.set_xticklabels([f"{int(d)} h" for d in storms[:, 4]], fontsize=7)
    ax_s.set_xlabel("the year's storms in order (hours above the berm)")
    ax_s.set_ylabel("water level (m MHW)")
    ax_s.set_xlim(-0.6, len(x) - 0.4)
    ax_s.set_ylim(berm - 0.3, max(max(g["rhigh"] for g in got), crest_mean) + 1.6)
    open_frame(ax_s)
    ax_s.legend(handles=[Line2D([], [], color=C["ACCENT"], lw=3.2, label="Rhigh over a dune gap: overwash"),
                         Line2D([], [], color=C["BASE"], lw=3.2, label="below every crest: collision"),
                         Line2D([], [], color=C["ADDED"], lw=1.0, label="mean crest, after growth"),
                         Line2D([], [], color=C["ADDED"], lw=0.8, ls=(0, (1, 1.5)), label="lowest crest, after growth"),
                         Line2D([], [], color=INK_MUTED, lw=0.6, ls=(0, (1, 2)), label="lowest crest, 1 January"),
                         Line2D([], [], color=INK_MUTED, lw=0.8, ls=(0, (3, 2)), label="berm")],
                loc="upper left", frameon=False, fontsize=7.5, ncol=3)
    _title(ax_s, 0, f"GIS {STORM_GIS}: the {year} storms, Rlow to Rhigh")

    cmap, norm, bounds = elevation_cmap()
    ext = (0, rows * CELL_M, 0, pre.shape[1] * CELL_M)
    axes = []
    for i, (grid, label) in enumerate([(pre, f"Before: the grid entering {year}"),
                                       (post, "After: storms, overwash and dune growth")]):
        ax = fig.add_subplot(gs[1 + i, 0])
        im = ax.imshow(grid.T, cmap=cmap, norm=norm, extent=ext, origin="lower",
                       interpolation="nearest", aspect="auto")
        ax.axvline(2 * CELL_M, color=INK, lw=0.5, ls=(0, (2, 1.5)))
        ax.set_ylabel("alongshore (m)")
        ax.set_yticks([0, 250, 500])
        ax.tick_params(labelbottom=False)
        ax.invert_xaxis()                        # ocean at the right
        open_frame(ax)
        _title(ax, 1 + i, label)
        axes.append(ax)
    cax = fig.add_subplot(gs[1:3, 1])
    cb = fig.colorbar(im, cax=cax, ticks=bounds[1:-1])
    cb.set_label("elevation (m MHW)")
    cb.outline.set_linewidth(0.5)

    ax_d = fig.add_subplot(gs[3, 0])
    lim = max(0.25, np.nanpercentile(np.abs(diff), 99.5))
    im_d = ax_d.imshow(diff.T, cmap="BrBG_r", norm=TwoSlopeNorm(0, -lim, lim), extent=ext,
                       origin="lower", interpolation="nearest", aspect="auto")
    ax_d.axvline(2 * CELL_M, color=INK, lw=0.5, ls=(0, (2, 1.5)))
    ax_d.set_ylabel("alongshore (m)")
    ax_d.set_yticks([0, 250, 500])
    ax_d.set_xlabel("landward of the first dune row (m)  ·  ocean at the right")
    ax_d.invert_xaxis()
    open_frame(ax_d)
    _title(ax_d, 3, "What moved: deposition brown, erosion green")
    cax_d = fig.add_subplot(gs[3, 1])
    cb = fig.colorbar(im_d, cax=cax_d)
    cb.set_label("change (m)")
    cb.outline.set_linewidth(0.5)

    q = np.asarray(b.QowTS)[t]
    out = save(fig, OUT / "barrier3d" / "barrier3d_storm_year.png")
    plt.close(fig)
    record_caption(out[0],
        f"Barrier3D through one storm year: GIS {STORM_GIS}, model year {year} (the natural 1996-2010 run, "
        "edgeBE, no road, no beach/dune management). (a) The year's storms in the order the model applies them: "
        "each bar runs from Rlow to Rhigh, the total water level range above MHW over the hours the water "
        f"stood above the berm ({berm:.2f} m MHW). Barrier3D grows the dune first and then tests EVERY storm of "
        "the year against that one crest (barrier3d.py:1358), although each storm also lowers the dune it "
        f"erodes. The shaded band spans that crest's lowest cell ({crest_min:.2f} m) to its mean "
        f"({crest_mean:.2f} m); the dotted grey line is the lowest crest on 1 January, before growth "
        f"({crest_jan.min():.2f} m). A storm whose Rhigh clears a gap overwashes through it (purple, labelled "
        f"with its overwash volume); {int(over.sum())} of {len(got)} did. Read from a replay of the year "
        "through the model's own update, which reproduces the saved grid exactly (storm_routing_check). "
        "(b, c) The domain's 10 m grid entering and leaving the year: the two dune rows (right of the dashed "
        "line) as crest elevation, then the interior, ocean at the right; elevation classes in m MHW, water "
        "below 0. (d) Their difference: overwash carries sand through the dune gaps into the interior (brown) "
        "and scours the throats it passes through (green); the dune rows lose height to storm erosion. "
        f"The year's overwash flux was {q:.0f} m3/m. The year was chosen so that the shoreline did not step a "
        "whole cell, so the grid did not shift between the two panels.")
    return out


# =============================================================================
# 2. BARRIER3D: ONE DOMAIN'S BUDGET OVER A WINDOW
# =============================================================================

def fig_barrier3d_domain_budget():
    spec = RUN_NATURAL
    c = load_run(spec)
    b = c.barrier3d[pad(BUDGET_GIS)]
    y0 = start_year(spec)
    nt = len(b.x_s_TS)
    years = y0 + np.arange(nt)
    xs, xt, xb = (np.asarray(v) * DAM for v in (b.x_s_TS, b.x_t_TS, b.x_b_TS))
    ref = xs[0]
    bw = beach_width_m(c)
    d_sf = b.DShoreface * DAM

    fig = plt.figure(figsize=figsize("double", height=6.4), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.25, 1, 1])
    ax_p = fig.add_subplot(gs[0, :])
    for t, col, lab in [(0, C_1984, str(years[0])), (nt - 1, C_1997, str(years[-1]))]:
        interior = np.asarray(b.DomainTS[t]).mean(axis=1) * DAM
        dune = (b.DuneDomain[t].mean() + b.BermEl) * DAM
        x_dune = xs[t] - ref + bw
        x_int = x_dune + 2 * CELL_M + (np.arange(len(interior)) + 0.5) * CELL_M
        px = np.r_[xt[t] - ref, xs[t] - ref, x_dune, x_dune, x_dune + 2 * CELL_M, x_int]
        pz = np.r_[-d_sf, 0, b.BermEl * DAM, dune, dune, interior]
        ax_p.plot(px, pz, color=col, lw=1.1, label=lab)
    ax_p.axhline(0, color=INK_MUTED, lw=0.6)
    ax_p.text(-390, 0.15, "MHW", fontsize=8, color=INK_MUTED, va="bottom")
    ax_p.annotate("shoreline x$_s$", xy=(0, 0), xytext=(-60, 3.0), fontsize=8, ha="right",
                  arrowprops=dict(arrowstyle="-", lw=0.5, color=INK_MUTED))
    ax_p.set_xlim(-400, xb[0] - ref + 100)
    ax_p.set_ylim(-5, None)
    ax_p.text(-395, -4.6, f"the shoreface continues to its toe x$_t$, {d_sf:.1f} m deep "
              f"and {xs[0] - xt[0]:.0f} m offshore", fontsize=8, color=INK_MUTED, va="bottom")
    ax_p.set_ylabel("elevation (m MHW)")
    ax_p.set_xlabel(f"cross-shore from the {years[0]} shoreline (m), ocean at the left")
    open_frame(ax_p)
    ax_p.legend(loc="lower right", frameon=False, title="alongshore-mean profile", title_fontsize=8)
    _title(ax_p, 0, f"GIS {BUDGET_GIS}: the profile at both ends of the window")

    ax_x = fig.add_subplot(gs[1, 0])
    for v, lab, col in [(xt, "shoreface toe x$_t$", C["BASE"]), (xs, "shoreline x$_s$", INK),
                        (xb, "back-barrier x$_b$", C["REF"])]:
        ax_x.plot(years, v - v[0], color=col, lw=1.3, label=lab)
    ax_x.axhline(0, color=INK_MUTED, lw=0.5)
    ax_x.set_ylabel("landward movement (m)")
    open_frame(ax_x)
    ax_x.legend(frameon=False, fontsize=7.5, loc="upper left")
    _title(ax_x, 1, "The three moving boundaries")

    ax_q = fig.add_subplot(gs[1, 1])
    qow, qsf = np.asarray(b.QowTS), np.asarray(b.QsfTS)
    w = 0.38
    ax_q.bar(years[1:] - w / 2, qow[1:], width=w, color=C["ACCENT"], label="overwash Q$_{ow}$")
    ax_q.bar(years[1:] + w / 2, qsf[1:], width=w, color=C["BASE"], label="shoreface Q$_{sf}$")
    ax_q.axhline(0, color=INK, lw=0.6)
    ax_q.set_ylabel("flux (m$^3$/m per year)")
    open_frame(ax_q)
    ax_q.legend(frameon=False, fontsize=7.5, loc="upper left")
    _title(ax_q, 2, "Annual sediment fluxes")

    ax_d = fig.add_subplot(gs[2, 0], sharex=ax_x)
    crest = [(d.max(axis=1).mean() + b.BermEl) * DAM for d in b.DuneDomain[:nt]]
    ax_d.plot(years, crest, color=C["ADDED"], lw=1.3)
    ax_d.set_ylabel("mean dune crest (m MHW)")
    ax_d.set_xlabel("year")
    open_frame(ax_d)
    _title(ax_d, 3, "Dune crest")

    ax_h = fig.add_subplot(gs[2, 1], sharex=ax_q)
    ax_h.plot(years, np.asarray(b.h_b_TS) * DAM, color=INK, lw=1.3)
    ax_h.set_ylabel("mean interior height h$_b$ (m MHW)")
    ax_h.set_xlabel("year")
    open_frame(ax_h)
    _title(ax_h, 4, "Barrier height")

    out = save(fig, OUT / "barrier3d" / "barrier3d_domain_budget.png")
    plt.close(fig)
    record_caption(out[0],
        f"One Barrier3D domain over a hindcast window: GIS {BUDGET_GIS}, {years[0]}-{years[-1]}, the natural "
        "run (edgeBE, no management). (a) The alongshore-mean cross-section the model carries at the start "
        f"(red) and end (blue): the shoreface from its toe x_t ({d_sf:.1f} m deep, set by BRIE from the wave "
        f"height) up to the shoreline x_s, a {bw:.0f} m beach at the berm, the two dune rows and the interior. "
        "(b) The shoreface toe, shoreline and back-barrier shoreline, as landward movement since the start. "
        "(c) The two fluxes that move them each year (Lorenzo-Trueba and Ashton 2014): overwash Q_ow carries "
        "sand from the front of the barrier to the interior, and the shoreface flux Q_sf moves sand onshore "
        "(positive) or offshore to keep the shoreface at its equilibrium slope. (d) The mean dune crest, which "
        "grows logistically between storms and is cut down by them. (e) The mean interior height above MHW. "
        "The shoreline here also carries BRIE's alongshore transport, which Barrier3D does not compute "
        "(see cascade_shoreline_split).")
    return out


# =============================================================================
# 3. BRIE: THE ALONGSHORE DIFFUSIVITY
# =============================================================================

def shoreline_angle_deg(x_s_m, dy):
    """BRIE's angle between each domain and the next (brie.py:820)."""
    return np.degrees(np.arctan2(np.roll(x_s_m, -1) - x_s_m, dy))


def fig_brie_diffusion():
    spec = RUN_NATURAL
    c = load_run(spec)
    br = c.brie
    dy = float(br._dy)
    cd = np.asarray(br._coast_diff)                   # m2/yr, index i <-> theta = 90 - i
    theta_axis = 90 - np.arange(len(cd))
    ang = np.deg2rad(np.linspace(-89.5, 89.5, 360))
    raw = np.cos(ang) ** 0.2 * (1.2 * np.sin(ang) ** 2 - np.cos(ang) ** 2)

    x_s0 = np.array([b.x_s_TS[0] for b in c.barrier3d]) * DAM
    theta = shoreline_angle_deg(x_s0, dy)
    idx = np.clip(np.round(90 - theta).astype(int), 1, br._wave_climl)
    k_dom = np.maximum(0, cd[idx])
    rp = real_pads()
    gis_axis = np.arange(1, 91)

    fig = plt.figure(figsize=figsize("double", height=6.2), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1, 1])
    ax_r = fig.add_subplot(gs[0, 0])
    ax_r.plot(np.degrees(ang), -raw, color=INK, lw=1.3)
    ax_r.axhline(0, color=INK_MUTED, lw=0.6)
    for s in (-1, 1):
        ax_r.axvline(s * 42.4, color=C["ACCENT"], lw=0.7, ls=(0, (3, 2)))
    ax_r.fill_between(np.degrees(ang), 0, -raw, where=-raw < 0, color=C["ACCENT_FILL"], lw=0)
    ax_r.text(0, 0.35, "diffusive:\nbumps smooth out", ha="center", fontsize=8)
    ax_r.text(66, -0.5, "anti-\ndiffusive", ha="center", fontsize=8, color=C["ACCENT"])
    ax_r.set_xlabel("wave angle relative to the shoreline (°)")
    ax_r.set_ylabel("relative diffusivity")
    ax_r.set_xlim(-90, 90)
    open_frame(ax_r)
    _title(ax_r, 0, "One wave direction")

    ax_c = fig.add_subplot(gs[0, 1])
    ax_c.plot(theta_axis, cd / 1e6, color=INK, lw=1.3, label="wave-climate averaged")
    ax_c.axhline(0, color=INK_MUTED, lw=0.6)
    th_r = theta[rp]
    ax_c.plot(th_r, np.full_like(th_r, -0.012), "|", color=C_1997, ms=7, mew=0.8,
              label="Hatteras domains (1996)", clip_on=False)
    ax_c.set_ylim(-0.02, None)
    ax_c.set_xlim(-60, 60)
    ax_c.set_xlabel("shoreline angle between neighbouring domains (°)")
    ax_c.set_ylabel("diffusivity (10$^6$ m$^2$/yr)")
    open_frame(ax_c)
    ax_c.legend(frameon=False, fontsize=7.5, loc="upper left")
    _title(ax_c, 1, "Averaged over the wave climate")

    ax_x = fig.add_subplot(gs[1, :])
    ax_x.plot(gis_axis, x_s0[rp] - x_s0[rp].min(), color=INK, lw=1.3)
    ax_x.set_ylabel("shoreline x$_s$ (m, landward up)")
    ax_x.set_xlim(0.5, 90.5)
    ax_x.tick_params(labelbottom=False)
    town_bands(ax_x)
    open_frame(ax_x)
    _title(ax_x, 2, "The BRIE shoreline: one point per domain, 500 m apart")

    ax_k = fig.add_subplot(gs[2, :], sharex=ax_x)
    ax_k.plot(gis_axis, theta[rp], color=C["BASE"], lw=1.0, label="angle to the next domain (°)")
    ax_k.set_ylabel("angle (°)")
    ax_k.set_xlabel(DOMAIN_AXIS_LABEL)
    open_frame(ax_k)
    r_ipl = k_dom * 1.0 / 2 / dy ** 2
    ax_k.scatter(gis_axis, theta[rp], c=np.where(r_ipl[rp] > 0, C_1997, C["ACCENT"]), s=9, zorder=3)
    ax_k.legend(handles=[Line2D([], [], color=C_1997, marker="o", ls="", ms=3.5, label="diffuses"),
                         Line2D([], [], color=C["ACCENT"], marker="o", ls="", ms=3.5,
                                label="negative diffusivity, clamped to 0")],
                frameon=False, fontsize=7.5, loc="upper left", ncol=2)
    _title(ax_k, 3, "Where each domain sits on that curve")

    out = save(fig, OUT / "brie" / "brie_diffusion.png")
    plt.close(fig)
    record_caption(out[0],
        "How BRIE moves sand alongshore in CASCADE. (a) The angle dependence of alongshore transport "
        "divergence for a single wave direction, cos^0.2(θ)(1.2 sin²θ - cos²θ) with the sign flipped so that "
        "positive is diffusive: waves within ~42° of shore-normal smooth a bump in the shoreline; at higher "
        "angles the shoreline is anti-diffusive and bumps grow. (b) The same term convolved with the wave "
        f"climate (Hs {c._wave_height} m, Tp {c._wave_period} s, asymmetry {c._wave_asymmetry}, high-angle "
        f"fraction {c._wave_angle_high_fraction}), as a diffusivity in m²/yr against the shoreline angle between "
        "neighbouring domains; ticks are the Hatteras domains at the start of the 1996 window. BRIE clamps a "
        "negative diffusivity to zero (brie.py:1293). (c) The BRIE shoreline, one position per 500 m domain, "
        "which carries the island's measured planform (the shoreline offset, in metres). (d) The angle between "
        "each domain and the next; the implicit solve each year is x_s(t+1) = x_s + r ∇²x_s + Δx_s(Barrier3D), "
        "with r = K Δt / 2Δy². GIS 1 is Cape Point, GIS 90 Pea Island.")
    return out


ORDER_START, ORDER_YEARS = 1996, 14


def brie_alone(offset_m, asymmetry, c):
    """BRIE's alongshore diffusion on its own (no Barrier3D, no source/sink):
    the padded offset added to BRIE's straight initial shoreline, run
    ORDER_YEARS annual steps. Returns the change in x_s, m, landward positive."""
    import warnings
    from brie import Brie
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        br = Brie(barrier_model=False, ast_model=True, inlet_model=False, b3d=True,
                  wave_height=c._wave_height, wave_period=c._wave_period,
                  wave_asymmetry=asymmetry,
                  wave_angle_high_fraction=c._wave_angle_high_fraction,
                  alongshore_section_length=c.brie._dy,
                  alongshore_section_count=offset_m.size,
                  time_step=1, time_step_count=ORDER_YEARS + 1)
    start = br.x_s + offset_m
    br.x_s[:] = start
    # the width check prints "Barrier Drowned" (x_b is not offset); it does not stop the solve
    import contextlib, io
    with contextlib.redirect_stdout(io.StringIO()):
        for _ in range(ORDER_YEARS):
            br.update()
    return br.x_s - start


def fig_brie_domain_order():
    from site_layer.hat_topo_version import offset_file
    c = load_run(RUN_NATURAL)
    a = float(c._wave_asymmetry)
    path = offset_file(ORDER_START)
    off = np.loadtxt(path, skiprows=1, delimiter=",")
    rp = real_pads()
    gis_axis = np.arange(1, 91)

    fwd = brie_alone(off, a, c)[rp]
    rev = brie_alone(off[::-1], a, c)[::-1][rp]
    rev_mirror = brie_alone(off[::-1], 1 - a, c)[::-1][rp]
    rms = lambda v: float(np.sqrt(np.mean(v ** 2)))  # noqa: E731

    fig, axes = plt.subplots(2, 1, sharex=True, constrained_layout=True,
                             figsize=figsize("double", height=4.6),
                             gridspec_kw={"height_ratios": [1.25, 1]})
    ax = axes[0]
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.plot(gis_axis, fwd, color=INK, lw=1.4,
            label=f"south → north (as run), asymmetry {a:g}")
    ax.plot(gis_axis, rev, color=C["ACCENT"], lw=1.2,
            label=f"north → south, asymmetry {a:g}")
    ax.plot(gis_axis, rev_mirror, color=C["BASE"], lw=1.1, ls=(0, (4, 2)),
            label=f"north → south, asymmetry {1 - a:g}")
    ax.set_ylabel(f"{ORDER_YEARS}-yr shoreline change (m,\nlandward up)")
    ax.set_xlim(0.5, 90.5)
    town_bands(ax)
    open_frame(ax)
    _title(ax, 0, "BRIE alongshore diffusion, domains fed in either order")

    ax = axes[1]
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.plot(gis_axis, rev - fwd, color=C["ACCENT"], lw=1.2)
    ax.plot(gis_axis, rev_mirror - fwd, color=C["BASE"], lw=1.1, ls=(0, (4, 2)))
    ax.set_ylabel("difference from\nas run (m)")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    open_frame(ax)
    _title(ax, 1, "Difference from the order as run")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3, frameon=False, fontsize=7.5)

    out = save(fig, OUT / "brie" / "brie_domain_order.png")
    plt.close(fig)
    record_caption(out[0],
        "Does the order of the domains matter to BRIE? BRIE's alongshore step alone (no Barrier3D, no "
        f"source/sink), run for {ORDER_YEARS} annual steps on the {ORDER_START} island offset "
        f"({path.relative_to(REPO).as_posix()}), with the wave climate of the adopted setup (Hs "
        f"{c._wave_height} m, Tp {c._wave_period} s, high-angle fraction {c._wave_angle_high_fraction}). "
        "(a) Shoreline change, landward positive, with the domains handed to BRIE south to north as the "
        f"model runs (black), reversed with the same asymmetry {a:g} (purple), and reversed with the "
        f"asymmetry swapped to {1 - a:g} (grey dashed); reversed runs are flipped back onto the GIS axis. "
        "(b) Each reversed run minus the order as run. BRIE takes the shoreline angle between each domain "
        "and the next, and the asymmetry decides which facing of the shoreline diffuses faster, so "
        f"reversing the order with asymmetry {a:g} changes the result by as much as the result itself "
        f"(RMS {rms(rev - fwd):.0f} m against {rms(fwd):.0f} m), while reversing it and swapping the "
        f"asymmetry recovers it to RMS {rms(rev_mirror - fwd):.0f} m; the rest is the one-sided angle "
        f"(domain to next domain). With the order as run, asymmetry {a:g} is the fraction of waves from "
        "the north, driving sand south toward Cape Point. GIS 1 is Cape Point, GIS 90 Pea Island; "
        "the 15 buffer domains at each end are run but not drawn.")
    return out


def fig_brie_domain_orientation():
    """North-up map of the 90 domains beside BRIE's array, same orientation."""
    import geopandas as gpd
    from site_layer import hat_map_layers as ml
    from site_layer.hat_observed_rates import DOMAIN_BOXES
    crs = "EPSG:26918"
    a = float(load_run(RUN_NATURAL)._wave_asymmetry)
    dom = gpd.read_file(DOMAIN_BOXES).to_crs(crs)
    dom = dom[(dom.domain_id >= 1) & (dom.domain_id <= 90)].sort_values("domain_id")
    outline = gpd.read_file(ml.ISLAND_OUTLINE).to_crs(crs)
    cen = np.c_[dom.geometry.centroid.x, dom.geometry.centroid.y]
    n_pad, first_pad = DOM.total_domains, pad(1)

    fig = plt.figure(figsize=figsize("double", height=7.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1])

    # ---- (a) the map, north up ------------------------------------------------
    ax = fig.add_subplot(gs[0])
    x0, x1, y0, y1 = 441_000, 470_000, 3_895_500, 3_947_500
    ax.set_facecolor("#e9eff4")
    outline.plot(ax=ax, color="#ede9df", edgecolor="0.55", lw=0.4, zorder=1)
    dom.boundary.plot(ax=ax, color=INK_MUTED, lw=0.35, zorder=2)
    for gis, (cx, cy), g in zip(dom.domain_id, cen, dom.geometry):
        if gis in (1, 90) or gis % 10 == 0:
            ax.add_patch(plt.Polygon(np.asarray(g.exterior.coords)[:, :2], closed=True,
                                     facecolor=C["ACCENT_FILL"], edgecolor=C["ACCENT"], lw=0.6,
                                     zorder=3))
            ax.text(g.bounds[2] + 250, cy, f"GIS {gis}", fontsize=7.5, va="center",
                    color=INK, zorder=5)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(True)
    ax.text(466_500, 3_925_000, "Atlantic\nOcean", ha="center", fontsize=8.5,
            color=INK_MUTED, style="italic")
    ax.text(444_000, 3_935_000, "Pamlico\nSound", ha="center", fontsize=8.5,
            color=INK_MUTED, style="italic")
    ax.text(cen[0, 0] - 1500, cen[0, 1] - 2200, "Cape Point", ha="center", fontsize=8)
    ax.text(cen[-1, 0] - 3500, cen[-1, 1] + 1500, "Pea Island", ha="center", fontsize=8)

    # waves offshore: the asymmetry's share from the north, the rest from the south
    wx, wy = 465_000, 3_909_000
    for frac, dy_sign, label in ((a, 1, f"{a:.0%} of waves\nfrom the north"),
                                 (1 - a, -1, f"{1 - a:.0%} from\nthe south")):
        ax.annotate("", xy=(wx - 2600, wy), xytext=(wx + 1400, wy + dy_sign * 4000),
                    arrowprops=dict(arrowstyle="-|>", color=C_1997, lw=5 * frac,
                                    mutation_scale=8 + 10 * frac, shrinkA=0, shrinkB=0))
        ax.text(wx + 1400, wy + dy_sign * 4700, label, ha="center",
                va="bottom" if dy_sign > 0 else "top", fontsize=7.5, color=C_1997)
    # the drift this wave climate implies physically; BRIE's diffusion has no net-flux term
    ax.annotate("", xy=(cen[8, 0] + 6200, cen[8, 1]), xytext=(cen[38, 0] + 6200, cen[38, 1]),
                arrowprops=dict(arrowstyle="-|>", color=C["ADDED"], lw=2.2, mutation_scale=16,
                                ls=(0, (4, 2))))
    ax.text(cen[38, 0] + 6200, cen[38, 1] + 700,
            "wave climate favours\nsouthward drift\n(not a BRIE flux)", fontsize=7,
            ha="center", color=C["ADDED"], va="bottom")
    # north arrow and a 5 km bar
    ax.annotate("", xy=(0.08, 0.97), xytext=(0.08, 0.90), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.0, mutation_scale=11))
    ax.text(0.08, 0.985, "N", transform=ax.transAxes, ha="center", va="bottom",
            fontweight="bold", fontsize=8.5)
    ax.plot([x1 - 6500, x1 - 1500], [y0 + 1500] * 2, color=INK, lw=2)
    ax.text(x1 - 4000, y0 + 2000, "5 km", ha="center", fontsize=7.5)
    _title(ax, 0, "The 90 domains, north up")

    # ---- (b) BRIE's array, same orientation ----------------------------------
    ax = fig.add_subplot(gs[1])
    for i in range(n_pad):
        real = first_pad <= i <= pad(90)
        marked = real and (i - first_pad + 1) in (1, 90) or (real and (i - first_pad + 1) % 10 == 0)
        fc = C["ACCENT_FILL"] if marked else ("white" if real else "0.85")
        ax.add_patch(Rectangle((0, i), 1, 1, facecolor=fc, edgecolor=INK_MUTED, lw=0.25))
    for i, txt in ((0, "index 0"),
                   (first_pad, f"index {first_pad} = GIS 1"),
                   (pad(90), f"index {pad(90)} = GIS 90"), (n_pad - 1, f"index {n_pad - 1}")):
        ax.text(1.25, i + 0.5, txt, va="center", fontsize=7.5)
    for gis in range(10, 90, 10):
        ax.text(1.25, pad(gis) + 0.5, f"index {pad(gis)} = GIS {gis}", va="center",
                fontsize=7, color=INK_MUTED)
    ax.text(-0.3, (first_pad - 1) / 2, f"{DOM.num_buffer_domains} buffer\ndomains",
            ha="right", va="center", fontsize=7.5, color=INK_MUTED)
    ax.text(-0.3, (pad(90) + n_pad) / 2, f"{DOM.num_buffer_domains} buffer\ndomains",
            ha="right", va="center", fontsize=7.5, color=INK_MUTED)
    # the ring: the last node's neighbour is the first
    ax.annotate("", xy=(-0.15, 0.5), xytext=(-0.15, n_pad - 0.5),
                arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=0.8, ls=(0, (3, 2)),
                                connectionstyle="bar,fraction=0.08", mutation_scale=9))
    ax.text(-2.0, n_pad / 2, "ring: index 119's\nneighbour is index 0", rotation=90,
            ha="center", va="center", fontsize=7, color=INK_MUTED)
    # the angle BRIE reads: node i to node i+1, i.e. looking north
    ax.annotate("", xy=(0.5, 62.5), xytext=(0.5, 58.5),
                arrowprops=dict(arrowstyle="-|>", color=C["ACCENT"], lw=1.2, mutation_scale=10))
    ax.set_xlim(-2.6, 5.2)
    ax.set_ylim(-1, n_pad + 1)
    ax.axis("off")
    _title(ax, 1, "BRIE's shoreline array")

    out = save(fig, OUT / "brie" / "brie_domain_orientation.png")
    plt.close(fig)
    record_caption(out[0],
        "Which way the domains run, and what BRIE's wave asymmetry means in that frame. (a) The 90 model "
        "domains (500 m alongshore, 2 km cross-shore, 5-scr/2-transect-frame/transect_domains/"
        "HAT_domains.json) on the island outline, north up; every tenth domain and the two ends are "
        "shaded and numbered. GIS 1 is at Cape Point and the numbers increase NORTHWARD to GIS 90 at "
        f"Pea Island. The blue arrows are the adopted wave asymmetry {a:g}: BRIE's asymmetry is the "
        "fraction of waves from the left looking offshore, which from this beach (looking east) is the "
        f"north, so {a:.0%} of waves come from the north and {1 - a:.0%} from the south, and physically "
        "that climate favours drift south toward Cape Point (amber, dashed). BRIE does NOT compute that "
        "drift: its alongshore step is a diffusion with no net-flux term, so a straight shoreline moves no "
        "sand under any asymmetry; the asymmetry only sets which shoreline orientations smooth fastest "
        "(brie_asymmetry_explained.png). (b) The array BRIE holds, one shoreline position per "
        f"domain, drawn in the same orientation: {DOM.num_buffer_domains} buffer domains at each end "
        f"(grey), GIS 1 at index {first_pad} and GIS 90 at index {pad(90)}, so index increases northward "
        "exactly as the GIS numbers do. The array is a ring (the last node's neighbour is the first), "
        "which is why the buffers exist. The purple arrow is the direction BRIE reads the shoreline "
        "angle, from each node to the next one (brie.py:820), i.e. northward. Reversing the order (GIS 90 "
        "at the bottom of the array) would make the same asymmetry mean waves from the SOUTH; see "
        "brie_domain_order.png for what that does to the result.")
    return out


def brie_diffusivity(asymmetry, high_fraction, theta_deg, c):
    """BRIE's wave-climate diffusivity (m2/yr) at shoreline angles theta, read
    from its own table exactly as the solve does (brie.py:1293, before the
    clamp at zero)."""
    import warnings
    from brie import Brie
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        br = Brie(barrier_model=False, ast_model=True, inlet_model=False, b3d=True,
                  wave_height=c._wave_height, wave_period=c._wave_period,
                  wave_asymmetry=asymmetry, wave_angle_high_fraction=high_fraction,
                  alongshore_section_length=500, alongshore_section_count=4,
                  time_step=1, time_step_count=3)
    cd, cl = np.asarray(br._coast_diff), br._wave_climl
    idx = np.clip(np.round(90 - np.asarray(theta_deg, float)).astype(int), 1, cl)
    return cd[idx]


def brie_cape(asymmetry, high_fraction, c, years=20, ny=40, amp=600.0, sigma=2.0,
              reverse=False):
    """A seaward cape (x_s negative = seaward) run through BRIE alone."""
    import contextlib, io, warnings
    from brie import Brie
    y = np.arange(ny)
    cape = -amp * np.exp(-0.5 * ((y - ny // 2) / sigma) ** 2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        br = Brie(barrier_model=False, ast_model=True, inlet_model=False, b3d=True,
                  wave_height=c._wave_height, wave_period=c._wave_period,
                  wave_asymmetry=asymmetry, wave_angle_high_fraction=high_fraction,
                  alongshore_section_length=500, alongshore_section_count=ny,
                  time_step=1, time_step_count=max(years, 1) + 1)
    base = br.x_s.copy()
    br.x_s[:] = base + (cape[::-1] if reverse else cape)
    with contextlib.redirect_stdout(io.StringIO()):
        for _ in range(years):
            br.update()
    out = br.x_s - base
    return cape, (out[::-1] if reverse else out)


def brie_cape_reversed(asymmetry, high_fraction, c, **kw):
    """The cape run with the array reversed and flipped back: its mirror image."""
    return brie_cape(asymmetry, high_fraction, c, reverse=True, **kw)[1]


def fig_brie_asymmetry_explained():
    c = load_run(RUN_NATURAL)
    a, h = float(c._wave_asymmetry), float(c._wave_angle_high_fraction)
    C_N, C_S = C_1997, C_1984          # waves from the north side / the south side
    fig = plt.figure(figsize=figsize("double", height=9.0), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.15, 1.25])

    # ---- (a) what the asymmetry counts --------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    edges = [-90, -45, 0, 45, 90]
    share = [a * h, a * (1 - h), (1 - a) * (1 - h), (1 - a) * h]
    for (lo, hi), s in zip(zip(edges[:-1], edges[1:]), share):
        ax.bar((lo + hi) / 2, s, width=44, color=C_N if hi <= 0 else C_S, alpha=0.85)
        ax.text((lo + hi) / 2, s + 0.01, f"{s:.2f}", ha="center", va="bottom", fontsize=7.5)
    ax.text(-45, max(share) + 0.09, f"a = {a:g} of all waves", ha="center", color=C_N, fontsize=8)
    ax.text(45, max(share) + 0.09, f"1 − a = {1 - a:g}", ha="center", color=C_S, fontsize=8)
    ax.axvline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(-90, 90)
    ax.set_ylim(0, max(share) + 0.17)
    ax.set_xticks([-90, -45, 0, 45, 90])
    ax.set_xlabel(r"wave angle $\varphi_0$ in BRIE (°)")
    ax.set_ylabel("share of waves")
    open_frame(ax)
    _title(ax, 0, "What the asymmetry a counts")

    # ---- (b) one wave: head-on smooths fastest -------------------------------
    ax = fig.add_subplot(gs[0, 1])
    rel = np.linspace(-89.5, 89.5, 359)
    r = np.deg2rad(rel)
    psi = -(np.cos(r) ** 0.2 * (1.2 * np.sin(r) ** 2 - np.cos(r) ** 2))
    ax.plot(rel, psi, color=INK, lw=1.3)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.annotate("head-on:\nsmooths fastest", xy=(0, 1), xytext=(28, 0.78), fontsize=7.5,
                arrowprops=dict(arrowstyle="-", color=INK_MUTED, lw=0.6))
    ax.set_xlim(-90, 90)
    ax.set_xticks([-90, -45, 0, 45, 90])
    ax.set_xlabel("angle between the wave and the shoreline's normal (°)")
    ax.set_ylabel("relative diffusivity")
    open_frame(ax)
    _title(ax, 1, "One wave: head-on smooths most")

    # ---- (c) the whole climate: which shoreline angle is head-on ------------
    ax = fig.add_subplot(gs[1, 0])
    th = np.arange(-60, 61)
    curves = ((1.0, C_N, r"a = 1 (all waves at $-\varphi_0$)"), (a, INK, f"a = {a:g} (adopted)"),
              (0.5, C["BASE"], "a = 0.5"), (0.0, C_S, r"a = 0 (all waves at $+\varphi_0$)"))
    peaks = {}
    for aa, col, lab in curves:
        d = brie_diffusivity(aa, 0.0, th, c) / 1e6
        ax.plot(th, d, color=col, lw=1.4 if aa == a else 1.1, label=lab,
                ls=(0, (4, 2)) if aa == 0.5 else "-")
        k = int(np.argmax(d))
        peaks[aa] = int(th[k])
        if aa != 0.5:
            ax.plot(th[k], d[k], "o", color=col, ms=4)
            ax.text(th[k] + (3 if aa == a else 0), d[k] + 0.03, f"{th[k]:+d}°",
                    ha="left" if aa == a else "center", va="bottom", fontsize=7.5, color=col)
    ax.axvline(0, color=INK_MUTED, lw=0.6)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlabel("shoreline angle θ from domain i to i+1 (°)")
    ax.set_ylabel("diffusivity (10$^6$ m$^2$/yr)")
    ax.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    ax.set_ylim(-0.45, 0.85)
    open_frame(ax)
    _title(ax, 2, "BRIE's diffusivity table")

    # ---- (d) the geometry on Hatteras, north up ------------------------------
    ax = fig.add_subplot(gs[1, 1])
    ax.set_aspect("equal")
    ax.axis("off")
    L = 1.0
    for x0, theta, col, side, share, faces in ((0.0, 23, C_N, "WEST", "north side\n(share a)", "ENE"),
                                               (3.1, -22, C_S, "EAST", "south side\n(share 1 − a)", "ESE")):
        t = np.deg2rad(theta)
        # node i (south) at the bottom; node i+1 500 m north; x_s landward = WEST
        p0 = np.array([x0, 0.0])
        p1 = p0 + np.array([-np.tan(t) * L, L])
        seg = p1 - p0
        nrm = np.array([seg[1], -seg[0]]) / np.hypot(*seg)          # the seaward normal
        land = np.array([p0 + [-1.2, 0], p1 + [-1.2, 0], p1, p0])
        ax.add_patch(plt.Polygon(land, closed=True, facecolor="#ede9df", edgecolor="none"))
        ax.plot([x0, x0], [0, L], color=INK_MUTED, lw=0.6, ls=(0, (2, 2)))
        ax.plot(*np.c_[p0, p1], color=INK, lw=2.0)
        for p, lab in ((p0, "i"), (p1, "i+1")):
            ax.plot(*p, "o", color=INK, ms=4)
            ax.text(p[0] - 0.08, p[1], lab, ha="right", va="center", fontsize=8)
        mid = (p0 + p1) / 2
        ax.annotate("", xy=mid + nrm * 0.05, xytext=mid + nrm * 1.0,
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=2.2, mutation_scale=14))
        ax.text(*(mid + nrm * 1.08), f"waves from\nthe {share}", color=col, fontsize=7,
                ha="left", va="center")
        ax.text(x0 - 0.1, -0.14, f"θ = {theta:+d}°", ha="center", va="top", fontsize=8.5,
                fontweight="bold")
        ax.text(x0 - 0.1, -0.42, f"steps {side} going north,\nso it faces {faces}",
                ha="center", va="top", fontsize=7, color=col)
    ax.annotate("", xy=(-1.25, 1.55), xytext=(-1.25, 1.1),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.0, mutation_scale=10))
    ax.text(-1.25, 1.6, "N", ha="center", va="bottom", fontweight="bold", fontsize=8.5)
    ax.text(-0.75, 1.12, "land\n(x$_s$ larger)", fontsize=7, color=INK_MUTED, va="bottom")
    ax.text(0.55, 0.02, "ocean", fontsize=7, color=INK_MUTED, va="bottom")
    ax.set_xlim(-1.5, 5.0)
    ax.set_ylim(-1.0, 1.9)
    _title(ax, 3, "θ on Hatteras, north up")

    # ---- (e) a cape: which flank BRIE smooths faster ------------------------
    ny = 40
    km = np.arange(ny) * 0.5
    cape = brie_cape(1.0, 0.0, c, years=0, ny=ny)[0]
    link_theta = np.degrees(np.arctan2(np.roll(cape, -1) - cape, 500.0))[:-1]
    ax = fig.add_subplot(gs[2, 0])
    d_link = brie_diffusivity(1.0, 0.0, link_theta, c) / 1e6
    norm = plt.Normalize(0, float(d_link.max()))
    cmap = plt.get_cmap("Greys")
    for i in range(ny - 1):
        ax.plot(-cape[i:i + 2], km[i:i + 2], color=cmap(0.25 + 0.75 * norm(max(d_link[i], 0))),
                lw=3.2, solid_capstyle="round")
    ax.annotate("", xy=(560, 12.2), xytext=(1150, 15.8),
                arrowprops=dict(arrowstyle="-|>", color=C_N, lw=2.0, mutation_scale=12))
    ax.text(1150, 16.1, "a = 1: all waves from\nthe north side", color=C_N, fontsize=7.2,
            ha="center", va="bottom")
    ax.text(660, 11.4, "north flank faces the waves:\nfast smoothing (dark)", fontsize=7.2,
            ha="left", va="center")
    ax.text(660, 8.6, "south flank faces away:\nslow smoothing (light)", fontsize=7.2,
            ha="left", va="center")
    ax.set_xlim(-100, 1500)
    ax.set_ylim(-0.5, 19.5)
    ax.set_xlabel("distance seaward (m), ocean to the right")
    ax.set_ylabel("alongshore (km), north up")
    open_frame(ax)
    _title(ax, 4, "What BRIE's table does to a cape")

    # ---- (f) what the asymmetry changes after 20 years ----------------------
    ax = fig.add_subplot(gs[2, 1], sharey=fig.axes[-1])
    runs = {aa: brie_cape(aa, 0.0, c, years=20, ny=ny)[1] for aa in (1.0, 0.5, 0.0)}
    diff1 = -(runs[1.0] - runs[0.5])      # seaward positive
    diff0 = -(runs[0.0] - runs[0.5])
    ax.axvline(0, color=INK_MUTED, lw=0.6)
    ax.axhspan(8.5, 11.5, color="0.94", lw=0, zorder=0)
    ax.text(0.02, 10.0, "cape", transform=ax.get_yaxis_transform(), ha="left", va="center",
            fontsize=7, color=INK_MUTED)
    ax.plot(diff1, km, color=C_N, lw=1.4, label="a = 1, all from the north side")
    ax.plot(diff0, km, color=C_S, lw=1.4, label="a = 0, all from the south side")
    ax.set_xlabel("20-yr shoreline minus the a = 0.5 run (m, seaward +)")
    ax.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1)
    ax.tick_params(labelleft=False)
    open_frame(ax)
    _title(ax, 5, "What the asymmetry changes in 20 years")

    out = save(fig, OUT / "brie" / "brie_asymmetry_explained.png")
    plt.close(fig)
    straight = float(np.abs(brie_cape(1.0, 0.0, c, years=20, ny=ny, amp=0.0)[1]).max())
    mirror = float(np.abs(runs[0.0] - brie_cape_reversed(1.0, 0.0, c, years=20, ny=ny)).max())
    record_caption(out[0],
        "How BRIE's wave asymmetry acts, and which way it points on Hatteras. (a) BRIE draws wave angles "
        f"φ₀ from four 45° bins (brie/waves.py); the asymmetry a is the share of waves at NEGATIVE angles, "
        f"here the adopted a = {a:g} with high-angle fraction {h:g}. (b) For one wave, the Ashton & Murray "
        "(2006) diffusivity term cos^0.2(1.2 sin² − cos²), sign flipped so positive smooths the shoreline: "
        "a wave hitting the shoreline head-on smooths it fastest; beyond ~42° it sharpens it. (c) BRIE "
        "averages (b) over the wave climate into a table of diffusivity against the shoreline angle θ "
        "between neighbouring domains (brie.py:380-400) and reads it at θ each year (brie.py:1293). The "
        "peak is the orientation the waves hit head-on: θ ≈ "
        f"{peaks[1.0]:+d}° with every wave at a negative angle (a = 1), {peaks[0.0]:+d}° with every wave at a "
        f"positive angle (a = 0), {peaks[0.5]:+d}° at a = 0.5 (the table's 1° rounding) and "
        f"{peaks[a]:+d}° at the adopted a = {a:g}. So the waves counted by a are head-on to links with "
        "POSITIVE θ. (d) What a positive θ is on Hatteras. θ = atan((x_s[i+1] − x_s[i]) / 500 m); x_s grows "
        "LANDWARD (west) and index i grows NORTHWARD (GIS 1 at Cape Point), both checked against the "
        "dune-line coordinates. A positive θ is a link that steps west going north, which faces "
        "east-north-east, so the waves counted by a come from the NORTH side of shore-normal and the rest "
        "from the south side. (e) A 600 m seaward cape coloured by BRIE's diffusivity with a = 1: the "
        "north flank (positive θ) faces the waves and smooths fast (dark), the south flank faces away "
        "(light). (f) What the asymmetry changes: the same cape after 20 years at a = 1 and at a = 0, "
        "each minus the a = 0.5 run. The difference is mostly a near-uniform shift of the whole shoreline "
        "(landward at a = 1, seaward at a = 0) with a small tilt between the flanks, and a = 0 is not an "
        f"exact mirror of a = 1 (they differ by up to {mirror:.0f} m once reflected). Both follow from how "
        "BRIE writes the step: x_s changes by D·∂²x_s/∂y² with D read on the link from each domain to "
        "the next one north, not as the divergence of a flux, so when D differs between flanks sand is "
        "not conserved and the scheme is not mirror-symmetric. It has no net-flux term at all: a straight "
        f"shoreline moves at most {straight:.1f} m in 20 years under a = 1. The asymmetry therefore does "
        "not drive an alongshore drift in BRIE; it sets which shoreline orientations smooth fastest. "
        "Panels (c), (e) and (f) set the high-angle fraction to 0 so that only the asymmetry differs; "
        f"Hs {c._wave_height} m and Tp {c._wave_period} s as adopted. Reversing the domain order (GIS 90 "
        "at the low index) flips the sign of every θ, close to swapping a for 1 − a "
        "(brie_domain_order.png).")
    return out


# =============================================================================
# 4. CASCADE: CROSS-SHORE VS ALONGSHORE
# =============================================================================

def split_shoreline_change(c):
    """Invert BRIE's implicit solve year by year to recover the Barrier3D
    shoreline change it was handed, so each domain's total change splits
    exactly into what Barrier3D did (cross-shore) and what BRIE did
    (alongshore). Only valid for a run with no management: the managers move
    x_s between the two solves. All in metres, + = landward."""
    br = c.brie
    dy, dt = float(br._dy), 1.0
    cd = np.asarray(br._coast_diff)
    xs = np.array([np.asarray(b.x_s_TS) for b in c.barrier3d]).T * DAM     # (nt, ny)
    hb = np.array([np.asarray(b.h_b_TS) for b in c.barrier3d]).T            # dam
    qat = np.array([b._Qat for b in c.barrier3d])                            # dam3/dam/yr
    d_sf = c.barrier3d[0].DShoreface
    nt, ny = xs.shape
    b3d = np.zeros((nt - 1, ny))
    be = np.zeros((nt - 1, ny))
    for t in range(1, nt):
        old, new = xs[t - 1], xs[t]
        theta = shoreline_angle_deg(old, dy)
        idx = np.maximum(1, np.minimum(br._wave_climl, np.round(90 - theta).astype(int)))
        r = np.maximum(0, cd[idx] * dt / 2 / dy ** 2)
        lap = np.roll(old, -1) - 2 * old + np.roll(old, 1)
        # brie.py:1310 builds A row by row from that row's own r (periodic):
        # A x = (1 + 2r) x - r (x[i-1] + x[i+1]) = x - r lap(x).
        a_new = new - r * (np.roll(new, -1) - 2 * new + np.roll(new, 1))
        b3d[t - 1] = a_new - old - r * lap
        be[t - 1] = 2 * qat / (2 * hb[t] + d_sf) * DAM
    total = xs[-1] - xs[0]
    return total, b3d.sum(0), be.sum(0), total - b3d.sum(0), xs


def fig_cascade_shoreline_split():
    spec = RUN_NATURAL
    c = load_run(spec)
    total, b3d, be, ast, xs = split_shoreline_change(c)
    rp = real_pads()
    g = np.arange(1, 91)
    y0 = start_year(spec)
    yrs = xs.shape[0] - 1
    flip = -1.0      # + = seaward, the convention every rate figure uses

    fig = plt.figure(figsize=figsize("double", height=5.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.2, 1])
    ax = fig.add_subplot(gs[0])
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.plot(g, flip * (b3d - be)[rp], color=C["BASE"], lw=1.1, label="Barrier3D cross-shore (storms, shoreface, dunes)")
    ax.plot(g, flip * be[rp], color=C["ADDED"], lw=1.1, label="source/sink (BE)")
    ax.plot(g, flip * ast[rp], color=C_1997, lw=1.1, label="BRIE alongshore diffusion")
    ax.plot(g, flip * total[rp], color=INK, lw=1.8, label="total")
    ax.set_ylabel(f"shoreline change {y0}-{y0 + yrs} (m)")
    clip = 80.0
    ax.set_ylim(-clip, clip)
    for comp, col in [((b3d - be), C["BASE"]), (be, C["ADDED"]), (ast, C_1997), (total, INK)]:
        v = flip * comp[rp]
        for i in np.nonzero(np.abs(v) > clip)[0]:
            ax.annotate(f"{v[i]:+.0f}", xy=(g[i], np.sign(v[i]) * clip), xytext=(-4, -8 * np.sign(v[i])),
                        textcoords="offset points", ha="right", va="center", fontsize=7, color=col)
    ax.set_xlim(0.5, 90.5)
    ax.tick_params(labelbottom=False)
    town_bands(ax)
    open_frame(ax)
    ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="lower left")
    _title(ax, 0, "Each domain's net change, split by the model that made it (+ seaward)")

    ax_h = fig.add_subplot(gs[1], sharex=ax)
    cum = flip * (xs[:, rp] - xs[0, rp])
    lim = np.nanpercentile(np.abs(cum), 99)
    im = ax_h.imshow(cum, aspect="auto", origin="lower", cmap="RdBu", norm=TwoSlopeNorm(0, -lim, lim),
                     extent=(0.5, 90.5, y0 - 0.5, y0 + yrs + 0.5), interpolation="nearest")
    ax_h.set_xlabel(DOMAIN_AXIS_LABEL)
    ax_h.set_ylabel("year")
    cb = fig.colorbar(im, ax=ax_h, fraction=0.025, pad=0.01)
    cb.set_label("change since start (m)")
    cb.outline.set_linewidth(0.5)
    _title(ax_h, 1, "Cumulative shoreline change, blue seaward, red landward")

    out = save(fig, OUT / "cascade" / "cascade_shoreline_split.png")
    plt.close(fig)
    resid = np.abs(total - (b3d + ast)).max()
    record_caption(out[0],
        f"How CASCADE divides a shoreline's change between its two models: the natural run, {y0}-{y0 + yrs} "
        "(edgeBE, no management), positive seaward. (a) Each domain's net shoreline change (black) split into "
        "the Barrier3D cross-shore change (grey: overwash, shoreface flux and dune effects), the prescribed "
        "source/sink (amber, the edgeBE rates converted through the Barrier3D shoreline equation, so a rate "
        "of R m/yr moves the shoreline by 2R·D/(2h_b+D) m/yr rather than R) and BRIE's alongshore diffusion "
        "(blue). The split is exact: each year's BRIE implicit solve is inverted to recover the Barrier3D change "
        f"it was handed (residual {resid:.1e} m). Alongshore diffusion only moves sand between domains, so its "
        "blue line sums to zero over the padded array; it carries the planform curvature of the island and "
        "the ends where the source/sink rates sit. Values beyond ±80 m (the two end domains, where the edge rates are fed in and BRIE carries most of it straight into the buffers) are printed at the axis edge. (b) The cumulative total change year by year. "
        "GIS 1 is Cape Point, GIS 90 Pea Island.")
    return out


# =============================================================================
# 5. CASCADE: THE ISLAND AS THE COUPLED MODEL HOLDS IT
# =============================================================================

def draw_island(ax, c, t, pads, alpha_pad=0.35, borders=False):
    cmap, norm, _ = elevation_cmap()
    bw = beach_width_m(c)
    gis0 = pad(1)
    beach_col = cmap(norm(np.array([1.2])))[0]
    for p in pads:
        b = c.barrier3d[p]
        grid = plan_grid(b, t)
        rows = land_rows(grid, margin=2)
        grid = np.ma.masked_less_equal(grid[:rows], 0)
        xs = b.x_s_TS[t] * DAM
        x0 = (p - gis0) * 0.5
        a = 1.0 if pad(1) <= p <= pad(90) else alpha_pad
        ax.add_patch(Rectangle((x0, xs), 0.5, bw, facecolor=beach_col, lw=0, alpha=a, zorder=1))
        ax.imshow(grid, cmap=cmap, norm=norm, origin="lower", interpolation="nearest", aspect="auto",
                  extent=(x0, x0 + 0.5, xs + bw, xs + bw + rows * CELL_M), alpha=a, zorder=2)
        if borders:
            ax.axvline(x0, color="white", lw=0.8, zorder=3)
    xs_line = np.array([c.barrier3d[p].x_s_TS[t] for p in pads]) * DAM
    ax.plot((np.asarray(pads) - gis0) * 0.5 + 0.25, xs_line, color=INK, lw=0.8, zorder=4)
    ax.set_facecolor(C["WATER"])


def fig_cascade_island_grids():
    spec = RUN_NATURAL
    c = load_run(spec)
    y0 = start_year(spec)
    pads_all = np.arange(len(c.barrier3d))
    fig = plt.figure(figsize=figsize("double", height=6.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1])
    ax = fig.add_subplot(gs[0])
    draw_island(ax, c, 0, pads_all)
    x_lo = (pads_all[0] - pad(1)) * 0.5
    x_hi = (pads_all[-1] + 1 - pad(1)) * 0.5
    ax.set_xlim(x_lo, x_hi)
    ys = np.array([b.x_s_TS[0] for b in c.barrier3d]) * DAM
    ax.set_ylim(ys.min() - 300, ys.max() + 1900)
    ax.axvspan(x_lo, 0, color="white", alpha=0.25, lw=0, zorder=5)
    ax.axvspan(45, x_hi, color="white", alpha=0.25, lw=0, zorder=5)
    ax.text((x_lo + 0) / 2, ys.max() + 1700, "buffer", ha="center", va="top", fontsize=8)
    ax.text((45 + x_hi) / 2, ys.max() + 1700, "buffer", ha="center", va="top", fontsize=8)
    ax.set_xlabel("alongshore from the south end of GIS 1 (km)")
    ax.set_ylabel("cross-shore (m, landward up)")
    zoom = (pad(43), pad(48))
    zx0, zx1 = (zoom[0] - pad(1)) * 0.5, (zoom[1] + 1 - pad(1)) * 0.5
    ax.add_patch(Rectangle((zx0, ys[zoom[0]:zoom[1] + 1].min() - 60), zx1 - zx0, 1900, fill=False,
                           ec=INK, lw=0.8, zorder=6))
    open_frame(ax)
    _title(ax, 0, f"All {len(pads_all)} Barrier3D grids on the BRIE shoreline, {y0}")

    ax_z = fig.add_subplot(gs[1])
    zp = np.arange(zoom[0], zoom[1] + 1)
    draw_island(ax_z, c, 0, zp, borders=True)
    zys = ys[zp]
    ax_z.set_xlim(zx0, zx1)
    ax_z.set_ylim(zys.min() - 100, zys.max() + 1700)
    for p in zp:
        ax_z.text((p - pad(1)) * 0.5 + 0.25, zys.max() + 1650, f"GIS {DOM.pad_to_gis(p)}",
                  ha="center", va="top", fontsize=8)
    ax_z.set_xlabel("alongshore (km)")
    ax_z.set_ylabel("cross-shore (m)")
    open_frame(ax_z)
    _title(ax_z, 1, "Six domains: each a 50-column grid, stepped to its own shoreline")
    cmap, norm, bounds = elevation_cmap()
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=[ax, ax_z], ticks=bounds[1:-1], fraction=0.02, pad=0.01)
    cb.set_label("elevation (m MHW)")
    cb.outline.set_linewidth(0.5)

    out = save(fig, OUT / "cascade" / "cascade_island_grids.png", vector=False)
    plt.close(fig)
    record_caption(out[0],
        f"The island as CASCADE holds it at the start of the {y0} window. (a) Every Barrier3D domain's 10 m "
        "grid (two dune rows, then the interior, water masked) placed at its BRIE shoreline position (black "
        f"line) plus a {beach_width_m(c):.0f} m beach at the berm (the model's initial beach width). Each "
        "domain is 500 m alongshore; the 90 Hatteras domains run from Cape Point (0 km) to Pea Island (45 km), "
        "with buffer domains, faded, at both ends that close BRIE's periodic line. The planform is the "
        "measured shoreline offset: the island bends ~6 km cross-shore over the reach. Vertical exaggeration "
        "is large; the cross-shore axis is in metres, the alongshore in km. (b) Six domains near Avon at full "
        "width: the grids do not overlap or share cells; they are coupled only through the one shoreline "
        "number each hands BRIE every year.")
    return out


# =============================================================================
# 6. CASCADE: THE ANNUAL COUPLING LOOP
# =============================================================================

def fig_cascade_coupling_loop():
    fig = plt.figure(figsize=figsize("double", height=5.2))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 70)
    ax.axis("off")

    def box(x, y, w, h, title, body, fc):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                    fc=fc, ec=INK, lw=0.7))
        ax.text(x + w / 2, y + h - 1.6, title, ha="center", va="top", fontsize=9, fontweight="bold")
        ax.text(x + w / 2, y + h - 5.2, body, ha="center", va="top", fontsize=7.6, linespacing=1.35)

    def arrow(p, q, label, lx, ly, ha="center", rad=0.0):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=10, lw=0.9, color=INK,
                                     connectionstyle=f"arc3,rad={rad}"))
        if label:
            ax.text(lx, ly, label, ha=ha, va="center", fontsize=7.3, color=INK,
                    bbox=dict(fc="white", ec="none", pad=0.6))

    b3d_fc, brie_fc, hum_fc, in_fc = "#f0e6bf", "#dbe8f3", "#ece2f0", "0.95"
    box(2, 48, 26, 17, "Inputs (once)",
        "topography, dunes: dam MHW\nberm, MHW: m NAVD88\nshoreline offset: m\nstorms: dam MHW, hours\n"
        "RSLR, BE: m/yr   Hs, Tp: m, s", in_fc)
    box(37, 48, 26, 17, "1  Barrier3D, every domain",
        "each domain on its own core\nstorms: dune erosion, overwash\nshoreface flux, dune growth\n"
        "RSLR lowers the grid\nsource/sink Q$_{at}$\nworks in dam", b3d_fc)
    box(72, 48, 26, 17, "2  BRIE",
        "one shoreline, 500 m spacing\nangle-dependent diffusivity\nimplicit solve, periodic\n"
        "x$_s$ ← x$_s$ + r$\\nabla^2$x$_s$ + Δx$_s$\nworks in m", brie_fc)
    box(72, 14, 26, 17, "3  Back to Barrier3D",
        "shoreline moved by BRIE\nwhole-cell steps migrate the\ndune line and interior\n"
        "(update_dune_domain)\ndrowning check", b3d_fc)
    box(37, 14, 26, 17, "4  Human modules",
        "roadway manager: clear overwash\noff NC-12, rebuild dunes, relocate\n"
        "beach/dune manager: nourish,\nfilter overwash, rebuild dunes\n(per domain, where switched on)", hum_fc)
    box(2, 14, 26, 17, "5  Sync BRIE",
        "managed x$_t$, x$_s$, x$_b$, h$_b$\nwritten back to BRIE\n(dam × 10 → m)\nfor next year's solve",
        brie_fc)

    arrow((28.6, 57), (36.4, 57), "init: /10, − MHW", 32.5, 60.5)
    arrow((63.6, 57), (71.4, 57), "Δx$_s$, Δx$_t$, Δh$_b$\ndam × 10 → m", 67.5, 62)
    arrow((85, 47.4), (85, 31.6), "new x$_s$\nm ÷ 10 → dam", 85, 39.5, ha="center")
    arrow((71.4, 22.5), (63.6, 22.5), "grid, dune\n(dam)", 67.5, 27.5)
    arrow((36.4, 22.5), (28.6, 22.5), "", 0, 0)
    arrow((15, 31.6), (40, 47.4), "", 0, 0, rad=-0.25)
    ax.text(20, 41.5, "next year", fontsize=7.3, ha="center")
    ax.text(50, 6, "Groin callback (when on): after step 1, adds a source and sink to Δx$_s$ in metres "
            "at the groin's two domains.   Units: UNITS.md.", ha="center", fontsize=7.3, color=INK_MUTED)

    out = save(fig, OUT / "cascade" / "cascade_coupling_loop.png")
    plt.close(fig)
    record_caption(out[0],
        "CASCADE's annual loop (cascade/cascade_groin.py, Cascade.update) and the unit each exchange is made "
        "in. (1) Every Barrier3D domain runs its year independently: the year's storms erode dunes and "
        "overwash the interior, the shoreface relaxes toward its equilibrium slope, dunes grow, sea-level rise "
        "lowers the grid, and the source/sink rate adds its flux. Each returns the change in shoreface toe, "
        "shoreline and barrier height, converted from decametres to metres. (2) BRIE adds those changes to "
        "its single shoreline and diffuses it alongshore with an implicit, periodic solve. (3) The new "
        "shoreline goes back to each Barrier3D domain in decametres; when it has moved a whole 10 m cell the "
        "dune line and interior shift with it. (4) Where switched on, the roadway and beach/dune managers act "
        "on the grid. (5) The managed geometry is written back to BRIE for the next year. Barrier3D's "
        "elevations are decametres relative to MHW throughout; BRIE works in metres. The full unit contract "
        "is in UNITS.md at the repository root.")
    return out


# =============================================================================
# 7. THE MANAGEMENT MODULES
# =============================================================================

def fig_management_modules():
    nat, road = load_run(RUN_NATURAL), load_run(RUN_ROAD)
    y0 = start_year(RUN_ROAD)
    p = pad(ROAD_GIS)
    bn, br_ = nat.barrier3d[p], road.barrier3d[p]
    rw = road.roadways[p]
    t_end = len(br_.x_s_TS) - 1
    g_nat, g_road = plan_grid(bn, t_end), plan_grid(br_, t_end)
    n = min(len(g_nat), len(g_road))
    rows = min(land_rows(g_nat[:n], g_road[:n]), n)
    d = (g_road[:n] - g_nat[:n])[:rows]

    non, noff = load_run(RUN_NOURISH), load_run(RUN_NOURISH_OFF)
    y1 = start_year(RUN_NOURISH)
    xs_on = np.array([np.asarray(b.x_s_TS) for b in non.barrier3d]).T * DAM
    xs_off = np.array([np.asarray(b.x_s_TS) for b in noff.barrier3d]).T * DAM
    diff = -(xs_on - xs_off)                  # + = seaward of the unmanaged run
    t_fill = NOURISH_YEAR - y1 + 1           # first state after the fill year
    years = y1 + np.arange(xs_on.shape[0])

    fig = plt.figure(figsize=figsize("double", height=6.4), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.1], width_ratios=[1.35, 1])
    ax_m = fig.add_subplot(gs[0, 0])
    lim = max(0.3, np.nanpercentile(np.abs(d), 99.5))
    ext = (0, rows * CELL_M, 0, d.shape[1] * CELL_M)
    im = ax_m.imshow(d.T, cmap="BrBG_r", norm=TwoSlopeNorm(0, -lim, lim), extent=ext, origin="lower",
                     interpolation="nearest", aspect="auto")
    sb = float(rw._road_setback) if np.ndim(rw._road_setback) == 0 else float(np.asarray(rw._road_setback)[0])
    x_road = 2 * CELL_M + sb
    ax_m.axvspan(x_road, x_road + float(rw._road_width), color=C["ROAD"], alpha=0.12, lw=0)
    ax_m.text(x_road + float(rw._road_width) / 2, d.shape[1] * CELL_M * 0.97, "NC-12", ha="center", va="top", fontsize=8)
    ax_m.axvline(2 * CELL_M, color=INK, lw=0.5, ls=(0, (2, 1.5)))
    ax_m.set_ylabel("alongshore (m)")
    ax_m.set_yticks([0, 250, 500])
    ax_m.set_xlabel("landward of the first dune row (m)  ·  ocean at the right")
    ax_m.invert_xaxis()
    open_frame(ax_m)
    cb = fig.colorbar(im, ax=ax_m, fraction=0.04, pad=0.01)
    cb.set_label("road run − natural (m)")
    cb.outline.set_linewidth(0.5)
    _title(ax_m, 0, f"GIS {ROAD_GIS}, {y0 + t_end}: road run − natural")

    ax_c = fig.add_subplot(gs[0, 1])
    # State t is the grid on 1 January of y0 + t; model year t (storms of
    # calendar year y0 + t - 1) runs between states t - 1 and t. The roadway
    # series are written at index t by update t (roadway_manager.py:780).
    yrs = y0 + np.arange(t_end + 1)
    cr_n = [(dd.max(axis=1).mean() + bn.BermEl) * DAM for dd in bn.DuneDomain[:t_end + 1]]
    cr_r = [(dd.max(axis=1).mean() + br_.BermEl) * DAM for dd in br_.DuneDomain[:t_end + 1]]
    ax_c.plot(yrs, cr_n, color=C["BASE"], lw=1.2, label="natural")
    ax_c.plot(yrs, cr_r, color=C["ACCENT"], lw=1.2, label="road managed")
    rb = np.nonzero(np.asarray(rw._dunes_rebuilt_TS)[:t_end + 1])[0]
    ax_c.plot(yrs[rb], np.asarray(cr_r)[rb], "o", color=C["ACCENT"], ms=4, label="dune rebuilt")
    ax_c.set_ylabel("mean dune crest (m MHW)")
    open_frame(ax_c)
    ax_c.legend(frameon=False, fontsize=7.5)
    _title(ax_c, 1, "Dune crest")

    ax_v = fig.add_subplot(gs[1, 0])
    ov = np.asarray(rw._road_overwash_volume)[:t_end + 1]
    ax_v.bar(yrs[1:] - 0.5, ov[1:], color=C["ROAD"], width=0.7)
    ax_v.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax_v.set_ylabel("overwash cleared (m$^3$)")
    ax_v.set_xlabel("year")
    open_frame(ax_v)
    _title(ax_v, 2, "Overwash cleared off NC-12")

    ax_s = fig.add_subplot(gs[1, 1])
    ax_s.plot(years, diff[:, pad(NOURISH_GIS)], color=C["ACCENT"], lw=1.3)
    ax_s.axvspan(NOURISH_YEAR, NOURISH_YEAR + 1, color=C["ACCENT_FILL"], alpha=0.4, lw=0)
    ax_s.axhline(0, color=INK_MUTED, lw=0.5)
    ax_s.set_ylabel("seaward of no-fill run (m)")
    ax_s.set_xlabel("year")
    open_frame(ax_s)
    _title(ax_s, 3, f"GIS {NOURISH_GIS}: the {NOURISH_YEAR} fill")

    ax_a = fig.add_subplot(gs[2, :])
    rp = real_pads()
    gis = np.arange(1, 91)
    cols = [C["BASE_FILL"], "#c2a5cf", "#9970ab", C["ACCENT"], "#40004b"]
    show = [t_fill - 1, t_fill, t_fill + 2, t_fill + 5, xs_on.shape[0] - 1]
    for col, t in zip(cols, show):
        if 0 <= t < xs_on.shape[0]:
            ax_a.plot(gis, diff[t, rp], color=col, lw=1.3, label=str(years[t] - 1))
    ax_a.axhline(0, color=INK_MUTED, lw=0.5)
    ax_a.set_xlim(55.5, 90.5)
    ax_a.set_xlabel(DOMAIN_AXIS_LABEL)
    ax_a.set_ylabel("seaward of no-fill run (m)")
    town_bands(ax_a)
    open_frame(ax_a)
    ax_a.legend(frameon=False, fontsize=7.5, ncol=5, loc="upper left", title="end of year", title_fontsize=7.5)
    _title(ax_a, 4, "BRIE spreads the fill alongshore")

    out = save(fig, OUT / "cascade" / "management_modules.png")
    plt.close(fig)
    record_caption(out[0],
        "What the two CASCADE management modules do to the grid. (a-c) The roadway manager, GIS "
        f"{ROAD_GIS}, {y0}-{y0 + t_end}: the road run minus the natural run of the same window on its final grid, ocean at the right "
        "(a), where overwash deposits have been bulldozed off NC-12 (grey band; green = lower than natural) "
        "and the dune rebuilt to its design height (brown along the dune rows); the mean dune crest in the two "
        "runs with the years a rebuild fired (b); and the volume cleared off the road each year (c), which "
        f"the module returns to the beach. (d, e) The beach/dune manager's nourishment: the {NOURISH_YEAR} "
        "Rodanthe emergency fill (GIS 84-89) as the shoreline position relative to a run identical except "
        f"that it has no fills (same road and beach/dune management), at GIS {NOURISH_GIS} through time with the fill year shaded (d) and along the northern reach at the end of "
        "selected years (e). The fill is placed as a volume per metre on the shoreface, which moves the "
        "shoreline seaward in one step; BRIE's alongshore diffusion then spreads it into the neighbouring "
        "domains, and the step rings at the grid scale for a few years (the Crank-Nicolson solve).")
    return out


# =============================================================================

FIGURES = {
    "barrier3d_storm_year": fig_barrier3d_storm_year,
    "barrier3d_domain_budget": fig_barrier3d_domain_budget,
    "brie_diffusion": fig_brie_diffusion,
    "brie_domain_order": fig_brie_domain_order,
    "brie_domain_orientation": fig_brie_domain_orientation,
    "brie_asymmetry_explained": fig_brie_asymmetry_explained,
    "cascade_shoreline_split": fig_cascade_shoreline_split,
    "cascade_island_grids": fig_cascade_island_grids,
    "cascade_coupling_loop": fig_cascade_coupling_loop,
    "management_modules": fig_management_modules,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", choices=sorted(FIGURES))
    args = ap.parse_args()
    apply_style()
    for name in args.only or FIGURES:
        out = FIGURES[name]()
        print(f"{name:28s} -> {out[0].relative_to(REPO)}")


if __name__ == "__main__":
    main()
