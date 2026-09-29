"""
overwash_routing_figures.py
==============================================================================
Where the water and the sand go during a storm, and whether Barrier3D's
overwash routing does what it should as a storm gets stronger.

    python scripts/figure_making/model/overwash_routing_figures.py [--only NAME]

Barrier3D routes each storm hour by hour over the domain grid (plus one dune
row in front and a strip of bay behind), tracking discharge and sediment flux
in every cell. None of that is saved. This script REPLAYS storms through the
model's own `Barrier3d.update()` from a saved run's grid, and reads the
routing arrays out of the running update at the point the storm finishes (a
line trace on that one frame; the model code is not copied or modified).

    storm_routing_check.png     the replay of the year's real storms against
                                the grid the run saved: if they differ, the
                                replay is not the model and nothing below
                                should be trusted
    storm_routing_ladder.png    the same grid hit by four storms of rising
                                strength, collision -> run-up through the
                                gaps -> run-up over the dune -> inundation:
                                the water that crossed each cell, the
                                elevation change, and the cross-shore
                                deposition profile
    storm_routing_hours.png     one run-up storm hour by hour: where the
                                water is and what the bed has done so far

Writes to output/figures/4-model-mechanics/<model>/, ocean at the RIGHT in every plan panel.
The domain and year are the storm-year example of model_mechanics_figures
(GIS 6, the 2006 storms, natural 1996-2010 run).
"""

from __future__ import annotations

import argparse
import functools
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LogNorm, TwoSlopeNorm  # noqa: E402

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, C, INK, INK_MUTED, CELL_M, figsize, figure_dir, save, record_caption,
    _title, open_frame,
)
import model_mechanics_figures as mm  # noqa: E402
from storm_replay import replay, regime  # noqa: E402

OUT = figure_dir("mechanics")   # one sub-folder per model: barrier3d/, brie/, cascade/, storm_routing/
DAM = 10.0
GIS, T = mm.STORM_GIS, mm.STORM_T
DURATION_H = 24
PERIOD_S = 10.0


# =============================================================================
# THE REPLAY
# =============================================================================

@functools.lru_cache(maxsize=1)
def saved_domain():
    c = mm.load_run(mm.RUN_NATURAL)
    return c.barrier3d[mm.pad(GIS)], mm.start_year(mm.RUN_NATURAL)


def real_storms():
    b, _ = saved_domain()
    s = b.StormSeries[b.StormSeries[:, 0] == T]
    return [(r[1] * DAM, r[2] * DAM, r[3], int(r[4])) for r in s]


def ladder():
    """Four storms against this grid's dune: below its lowest crest, between
    its lowest and mean crest, over its mean crest, and with Rlow above the
    gaps (inundation). Levels are set from the grid, so the ladder means the
    same thing on any domain."""
    b, _ = saved_domain()
    crest = (b.DuneDomain[T - 1].max(axis=1) + b.BermEl) * DAM
    lo, mean = crest.min(), crest.mean()
    berm = b.BermEl * DAM
    return [
        ("below every crest", lo - 0.25, berm + 0.05),
        ("through the gaps", 0.5 * (lo + mean), berm + 0.05),
        ("over the dune", mean + 0.5, berm + 0.2),
        ("Rlow above the gaps", mean + 1.0, mean + 0.3),
    ]


def interior_rows(s):
    """The routing domain is the dune row, the interior, then a bay strip."""
    return s["elevation"].shape[1]


# =============================================================================
# FIGURES
# =============================================================================

def extent(nrows, ncols):
    """Plan extent with the dune row at x = 0 and landward positive; the axis
    is then inverted so the ocean sits at the right."""
    return (-0.5 * CELL_M, (nrows - 0.5) * CELL_M, 0, ncols * CELL_M)


def fig_storm_routing_check():
    b, y0 = saved_domain()
    storms = real_storms()
    got, _ = replay(b, T, storms)
    last = got[-1]["elevation"][-1, 1:, :]
    saved = np.asarray(b.DomainTS[T]) * DAM
    n = min(len(saved), len(last))
    d = last[:n] - saved[:n]
    rows = mm.land_rows(saved[:n], margin=10)
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=2.6), constrained_layout=True,
                             gridspec_kw=dict(width_ratios=[1.6, 1]))
    ax = axes[0]
    lim = max(1e-3, np.abs(d[:rows]).max())
    im = ax.imshow(d[:rows].T, cmap="BrBG_r", norm=TwoSlopeNorm(0, -lim, lim), origin="lower",
                   extent=extent(rows, d.shape[1]), aspect="auto", interpolation="nearest")
    ax.invert_xaxis()
    ax.set_xlabel("landward of the dune (m)  ·  ocean at the right")
    ax.set_ylabel("alongshore (m)")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.01)
    cb.set_label("replay − saved (m)")
    _title(ax, 0, "Replay minus the saved grid")
    ax2 = axes[1]
    qs = [s["owloss"] for s in got]
    ax2.bar(np.arange(len(qs)), qs, color=[C["ACCENT"] if s["gaps"] else C["BASE"] for s in got])
    ax2.set_xticks(np.arange(len(qs)))
    ax2.set_xticklabels([f"{s['rhigh']:.1f}" for s in got], fontsize=7)
    ax2.set_xlabel("storm Rhigh (m MHW), in file order")
    ax2.set_ylabel("overwash (m$^3$/m)")
    open_frame(ax2)
    _title(ax2, 1, "Overwash per storm")
    out = save(fig, OUT / "storm_routing" / "storm_routing_check.png")
    plt.close(fig)
    q_saved = float(np.asarray(b.QowTS)[T])
    record_caption(out[0],
        f"The storm replay is the model: GIS {GIS}, the {y0 + T - 1} storms of the natural 1996-2010 run. "
        "(a) The interior after replaying the year's real storms through Barrier3d.update() from the grid the "
        f"run saved entering the year, minus the grid the run saved at the end of the year: largest difference "
        f"{np.abs(d).max():.2e} m. (b) Overwash volume per storm in the replay (purple where Rhigh reached a dune "
        f"gap); they sum to {sum(qs):.1f} m3/m against the run's recorded {q_saved:.1f} m3/m for the year. "
        "The storm-routing figures replay synthetic storms through the same code from the same grid.")
    print(f"  replay check: max |replay - saved| = {np.abs(d).max():.3e} m; "
          f"overwash {sum(qs):.2f} vs saved {q_saved:.2f} m3/m")
    return out


def fig_storm_routing_ladder():
    b, y0 = saved_domain()
    runs = []
    for label, rh, rl in ladder():
        got, _ = replay(b, T, [(rh, rl, PERIOD_S, DURATION_H)])
        runs.append((label, got[0]))

    nrow = len(runs)
    rows = max(int(np.nonzero((np.abs(s["elevation"][-1] - s["elevation"][0]) > 0.02).any(axis=1))[0].max())
               if (np.abs(s["elevation"][-1] - s["elevation"][0]) > 0.02).any() else 0 for _, s in runs) + 12
    rows = min(rows, runs[0][1]["elevation"].shape[1])
    ncols = runs[0][1]["elevation"].shape[2]
    water_max = max(s["discharge"].sum(0).max() for _, s in runs) or 1.0
    dz_lim = max(0.1, max(np.nanpercentile(np.abs(s["elevation"][-1] - s["elevation"][0]), 99.7) for _, s in runs))

    fig = plt.figure(figsize=figsize("double", height=8.4), constrained_layout=True)
    gs = fig.add_gridspec(nrow + 1, 3, width_ratios=[1, 1, 0.62], height_ratios=[1] * nrow + [0.07])
    im_w = im_z = None
    for i, (label, s) in enumerate(runs):
        water = s["discharge"].sum(0)[:rows] / s["substep"]         # m3 that crossed each cell
        dz = (s["elevation"][-1] - s["elevation"][0])[:rows]
        ext = extent(rows, ncols)
        ax_w = fig.add_subplot(gs[i, 0])
        im_w = ax_w.imshow(np.ma.masked_less_equal(water, 0).T, cmap="Blues",
                           norm=LogNorm(1.0, water_max), origin="lower",
                           extent=ext, aspect="auto", interpolation="nearest")
        ax_w.set_facecolor("0.96")
        ax_w.invert_xaxis()
        ax_w.set_ylabel("alongshore (m)")
        ax_w.set_yticks([0, 250, 500])
        _title(ax_w, 3 * i, label)

        ax_z = fig.add_subplot(gs[i, 1], sharey=ax_w)
        im_z = ax_z.imshow(dz.T, cmap="BrBG_r", norm=TwoSlopeNorm(0, -dz_lim, dz_lim), origin="lower",
                           extent=ext, aspect="auto", interpolation="nearest")
        ax_z.invert_xaxis()
        ax_z.tick_params(labelleft=False)
        _title(ax_z, 3 * i + 1, f"Rhigh {s['rhigh']:.1f} m, {regime(s)}")

        ax_p = fig.add_subplot(gs[i, 2])
        x = np.arange(rows) * CELL_M
        prof = dz.mean(axis=1)
        ax_p.fill_between(x, 0, prof, where=prof > 0, color="#a6761d", lw=0, step="mid")
        ax_p.fill_between(x, 0, prof, where=prof < 0, color="#1b9e77", lw=0, step="mid")
        ax_p.axhline(0, color=INK_MUTED, lw=0.5)
        ax_p.set_xlim(x[-1], x[0])
        ax_p.set_ylabel("mean change (m)")
        open_frame(ax_p)
        _title(ax_p, 3 * i + 2, f"{s['owloss']:.1f} m$^3$/m")
        if i < nrow - 1:
            for a in (ax_w, ax_z, ax_p):
                a.tick_params(labelbottom=False)
        else:
            ax_w.set_xlabel("landward of the dune (m)")
            ax_z.set_xlabel("ocean at the right →")
            ax_p.set_xlabel("landward of the dune (m)")

    cax_w = fig.add_subplot(gs[nrow, 0])
    cb = fig.colorbar(im_w, cax=cax_w, orientation="horizontal")
    cb.set_label("water through the cell, whole storm (m$^3$)")
    cax_z = fig.add_subplot(gs[nrow, 1])
    cb = fig.colorbar(im_z, cax=cax_z, orientation="horizontal")
    cb.set_label("elevation change (m)")

    out = save(fig, OUT / "storm_routing" / "storm_routing_ladder.png")
    plt.close(fig)
    crest = (b.DuneDomain[T - 1].max(axis=1) + b.BermEl) * DAM
    record_caption(out[0],
        f"The same Barrier3D grid (GIS {GIS}, entering {y0 + T - 1}, natural 1996-2010 run) hit by four "
        f"{DURATION_H} h storms of rising strength, each replayed through the model's own update from that "
        f"grid. The dune's lowest crest is {crest.min():.2f} m and its mean crest {crest.mean():.2f} m MHW; "
        "the four storms put Rhigh below every crest, between the two, and above the mean, and the last "
        "also lifts Rlow above the gaps, which switches Barrier3D to its inundation routing. Left: the "
        "water that crossed each cell over the storm (log scale from 1 m3; grey = dry). Middle: elevation change, "
        "deposition brown, erosion green. Right: the change averaged alongshore against distance landward "
        "of the dune, titled with the storm's overwash volume. The ocean is at the right; the dune row sits at 0 m, the bay strip Barrier3D appends "
        "for routing lies beyond the last interior row. What to check: no water and no change when Rhigh is "
        "below every crest; water entering only through the gaps and fanning landward in run-up; broad "
        "sheet flow and deposition reaching further into the interior under inundation; and in every case "
        "erosion at the dune and throats with deposition landward of them, not the reverse. Two things this "
        "shows that are the model's, not the storm's: a storm between the lowest and mean crest routes little "
        "or no water, because Barrier3D drops single-cell and trailing gap cells (storm_routing_response); and "
        "inundation moves far less sand than run-up at a higher water level, because its transport rule "
        "(Ki = 7.5e-6, and a momentum constant the code resets to 0) is much weaker than run-up's (Kr = 7.5e-5).")
    return out


def fig_storm_routing_hours():
    b, y0 = saved_domain()
    label, rh, rl = ladder()[2]
    got, _ = replay(b, T, [(rh, rl, PERIOD_S, DURATION_H)])
    s = got[0]
    nstep = s["elevation"].shape[0]
    moved = (np.abs(s["elevation"][-1] - s["elevation"][0]) > 0.02).any(axis=1)
    rows = min(s["elevation"].shape[1], (int(np.nonzero(moved)[0].max()) if moved.any() else 20) + 12)
    ncols = s["elevation"].shape[2]
    picks = sorted(set([0, nstep // 8, nstep // 3, nstep - 1]))
    qmax = s["discharge"][:, :rows].max() or 1.0
    dz_all = s["elevation"][:, :rows] - s["elevation"][0, :rows]
    dz_lim = max(0.1, np.nanpercentile(np.abs(dz_all[-1]), 99.7))
    ext = extent(rows, ncols)

    fig = plt.figure(figsize=figsize("double", height=7.2), constrained_layout=True)
    gs = fig.add_gridspec(len(picks) + 1, 2, height_ratios=[1] * len(picks) + [0.08])
    for i, k in enumerate(picks):
        hour = (k + 1) / s["substep"]
        ax = fig.add_subplot(gs[i, 0])
        imq = ax.imshow(np.ma.masked_less_equal(s["discharge"][k, :rows], 0).T, cmap="Blues",
                        norm=LogNorm(max(qmax * 1e-4, 1e-2), qmax), origin="lower", extent=ext,
                        aspect="auto", interpolation="nearest")
        ax.set_facecolor("0.96")
        ax.invert_xaxis()
        ax.set_ylabel("alongshore (m)")
        ax.set_yticks([0, 250, 500])
        _title(ax, 2 * i, f"hour {hour:.0f}: water moving")
        ax2 = fig.add_subplot(gs[i, 1], sharey=ax)
        imz = ax2.imshow(dz_all[k].T, cmap="BrBG_r", norm=TwoSlopeNorm(0, -dz_lim, dz_lim), origin="lower",
                         extent=ext, aspect="auto", interpolation="nearest")
        ax2.invert_xaxis()
        ax2.tick_params(labelleft=False)
        _title(ax2, 2 * i + 1, f"hour {hour:.0f}: bed change so far")
        if i < len(picks) - 1:
            ax.tick_params(labelbottom=False)
            ax2.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("landward of the dune (m)  ·  ocean at the right")
            ax2.set_xlabel("landward of the dune (m)  ·  ocean at the right")
    cb = fig.colorbar(imq, cax=fig.add_subplot(gs[-1, 0]), orientation="horizontal")
    cb.set_label("discharge through the cell (m$^3$/hr)")
    cb = fig.colorbar(imz, cax=fig.add_subplot(gs[-1, 1]), orientation="horizontal")
    cb.set_label("elevation change since the storm began (m)")

    out = save(fig, OUT / "storm_routing" / "storm_routing_hours.png")
    plt.close(fig)
    record_caption(out[0],
        f"One storm hour by hour: GIS {GIS} entering {y0 + T - 1}, a {DURATION_H} h storm with Rhigh "
        f"{s['rhigh']:.1f} m MHW ({label}, {regime(s)} routing), replayed through Barrier3D. Left: the "
        "discharge through each cell in that routing step (log scale): water enters at the dune gaps, where "
        "Rhigh exceeds the crest, and is routed landward cell by cell down the steepest descent, spreading "
        "between neighbours in proportion to slope. Right: the cumulative bed change since the storm began: "
        "sediment is picked up where flow is fast and steep (the dune gaps and throats) and dropped where it "
        "slows and spreads. The dune row is lowered linearly over the storm (Goldstein and Moore 2016), "
        f"so later hours route through lower gaps. Overwash for this storm: {s['owloss']:.1f} m3/m.")
    return out


def fig_storm_routing_response():
    """Overwash against storm strength on one grid, as the model is and with
    each suspected defect fixed in memory (storm_replay.FIXES)."""
    b, y0 = saved_domain()
    berm = b.BermEl * DAM
    got0, _ = replay(b, T, [(berm + 0.1, berm + 0.05, PERIOD_S, 1)])
    grown = got0[0]["crest_pre"]                     # the crest the storms meet
    levels = np.round(np.linspace(grown.min() - 0.2, grown.max() + 0.8, 16), 3)
    variants = [((), "as the model is", INK),
                (("gaps", "slice"), "gap cells fixed", C["BASE"]),
                (("gaps", "momentum", "slice"), "gap cells + inundation momentum fixed", C["ACCENT"])]
    cases = [("run-up", lambda rh: berm + 0.05), ("inundation", lambda rh: rh - 0.3)]
    res = {}
    for key, rl_of in cases:
        for fx, _, _ in variants:
            storms = [(rh, max(berm + 0.05, rl_of(rh)), PERIOD_S, DURATION_H) for rh in levels]
            out = []
            for st in storms:
                g, _ = replay(b, T, [st], fixes=fx)
                s0 = g[0]
                wet = int(((s0["discharge"][:, 0, :] > 0).any(axis=0)).sum())
                out.append((s0["owloss"], wet, regime(s0)))
            res[key, fx] = out
    over = np.array([(grown < rh).sum() for rh in levels])

    fig, axes = plt.subplots(2, 2, figsize=figsize("double", height=5.2), constrained_layout=True,
                             sharex=True)
    for j, (key, _) in enumerate(cases):
        ax = axes[0, j]
        for fx, lab, col in variants:
            ax.plot(levels, [o[0] for o in res[key, fx]], "o-", color=col, ms=3, lw=1.2, label=lab)
        ax.axvline(grown.min(), color=C["ADDED"], lw=0.8, ls=(0, (1, 1.5)))
        ax.axvline(grown.mean(), color=C["ADDED"], lw=1.0)
        ax.set_ylabel("overwash (m$^3$/m)")
        open_frame(ax)
        _title(ax, j, f"{key} routing (Rlow {'at the berm' if key == 'run-up' else '= Rhigh − 0.3 m'})")
        ax2 = axes[1, j]
        ax2.step(levels, over, where="mid", color=C["ADDED"], lw=1.2, label="dune cells Rhigh overtops")
        for fx, lab, col in variants:
            ax2.plot(levels, [o[1] for o in res[key, fx]], "o-", color=col, ms=3, lw=1.0)
        ax2.set_ylabel("dune cells (of 50)")
        ax2.set_xlabel("storm Rhigh (m MHW)")
        open_frame(ax2)
        _title(ax2, 2 + j, "overtopped vs given water")
    axes[0, 0].legend(frameon=False, fontsize=7.5, loc="upper left")
    axes[1, 0].legend(frameon=False, fontsize=7.5, loc="upper left")
    out = save(fig, OUT / "storm_routing" / "storm_routing_response.png")
    plt.close(fig)
    record_caption(out[0],
        f"Does overwash grow with the storm? GIS {GIS} entering {y0 + T - 1} (natural 1996-2010 run), one "
        f"{DURATION_H} h storm at a time replayed through Barrier3D, Rhigh stepped from below the lowest dune "
        f"crest to above the highest (crest after the year's dune growth: lowest {grown.min():.2f} m dotted, "
        f"mean {grown.mean():.2f} m solid; berm {berm:.2f} m). Left: Rlow at the berm, run-up routing. "
        "Right: Rlow 0.3 m below Rhigh, which puts the gaps in inundation routing once Rlow clears them. "
        "Top: overwash volume. Bottom: how many of the 50 dune cells Rhigh overtops (amber step) against how "
        "many actually receive water in the routing (points). Black is Barrier3D as it runs in every hindcast; "
        "grey fixes two gap-handling defects (DuneGaps drops the last overtopped cell of the last gap and "
        "any single-cell gap; the gap discharge slice start:stop drops each gap's last cell); purple also "
        "restores the inundation momentum constant C = Cx * AvgSlope, which a 2024 refactor resets to 0 "
        "before routing. The fixes are applied to an in-memory copy of the model for this figure only; "
        "barrier3d.py is unchanged.")
    for key, _ in cases:
        for fx, lab, _ in variants:
            print(f"  {key:10s} {lab:40s}", " ".join(f"{o[0]:5.1f}" for o in res[key, fx]))
    return out


FIGURES = {
    "storm_routing_check": fig_storm_routing_check,
    "storm_routing_ladder": fig_storm_routing_ladder,
    "storm_routing_hours": fig_storm_routing_hours,
    "storm_routing_response": fig_storm_routing_response,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", choices=sorted(FIGURES))
    args = ap.parse_args()
    apply_style()
    for name in args.only or FIGURES:
        out = FIGURES[name]()
        print(f"{name:24s} -> {out[0].relative_to(REPO)}")


if __name__ == "__main__":
    main()
