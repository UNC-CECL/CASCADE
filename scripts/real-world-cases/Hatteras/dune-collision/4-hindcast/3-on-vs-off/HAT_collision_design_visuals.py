"""
What one collision design and Cs looks like in the model (2026-11-19): shoreline, dune rows, and the model grid.

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_design_visuals.py faceretreat_equal 3

(design key, Cs in 1e-6; default face retreat, equal rows, Cs 3e-6, the design the three tests picked:
HAT_collision_three_tests.py). Each against no collision on the same starting rows. Writes, in the design's folder
(collision_designs.py), <key>_cs<Cs>_visuals/:
- position_1996_2009.png: shoreline position change against CoastSat (HAT_collision_on_off_position.figure)
- dune_rows_1996_2009.png: the dune along the island through time, model space (GIS domain x year): no collision's
  crest (the higher row), the design's front row, its back row, and its crest minus no collision's; m above the berm
- model_grid_<zone>.png: plan view of the model's own grid at the hotspot zones (Buxton GIS 9-15, Rodanthe 80-89),
  1996 / after Isabel / 2009, no collision above the design: the two dune rows ahead of the interior, each year
  placed at its own shoreline, so retreat shows as the island moving landward. m NAVD88

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-19
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
_HAT = next(p for p in _HERE.parents if (p / "hatteras_paths.py").exists())
sys.path.insert(0, str(_HAT))
sys.path.insert(0, str(_HERE.parent))
import collision_designs as CD  # noqa: E402
import HAT_collision_design_cs_sweep as DS  # noqa: E402

H = DS.H

# --- CONFIG ------------------------------------------------------------------
GIS = range(1, 91)
STORMS = {1998.7: "Bonnie", 1999.7: "Dennis/Floyd", 2003.7: "Isabel"}
GRID_ZONES = {"buxton": ("Buxton", range(9, 16)), "rodanthe": ("Rodanthe", range(80, 90))}
GRID_YEARS = (0, 8, 13)              # model states: start, 1 Jan after Isabel's year, end
CROSS_M = 500                        # plan view: this far landward of the 1996 shoreline
DUNE_MAX_M = 6.0                     # colour scale top, m above the berm
DIFF_M = 2.0                         # difference scale, +-
GREENS = "Greens"                    # height: one hue, light -> dark
DIVERGING = "RdBu"                   # lower than no collision red, higher blue, grey-white middle
# -----------------------------------------------------------------------------


def load(v):
    from cascade_pipeline.saved_runs import load_saved_cascade
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    rd = H.run_dir(v)
    meta = json.loads(next(rd.glob("*_run_metadata.json")).read_text(encoding="utf-8"))
    c = load_saved_cascade(rd / f"{rd.name}.npz")
    return c, DOM.for_total(len(c.barrier3d)), int(meta["period"]["first_model_year"])


def rows(c, d):
    """(states, 90*50, 2) dune row heights above the berm, m, GIS 1-90 in order."""
    return np.concatenate([np.asarray(c.barrier3d[d.gis_to_pad(g)].DuneDomain, float) * 10 for g in GIS], axis=1)


def plan(c, d, gis, t, x0):
    """Plan view of GIS domains at state t: (cross-shore cells, alongshore cells) m NAVD88, 10 m cells, row 0 at
    the 1996 shoreline x0 (dam); dunes at the shoreline, interior behind them; NaN where nothing is."""
    n = CROSS_M // 10
    out = []
    for g in gis:
        b = c.barrier3d[d.gis_to_pad(g)]
        z = np.full((n, 50), np.nan)
        dune = (np.asarray(b.DuneDomain, float)[t] + b._BermEl).T          # (2, 50) dam MHW
        inter = np.asarray(b.DomainTS[t], float)                           # (cross, 50) dam MHW
        col = np.vstack([dune, inter])
        k = int(round(np.asarray(b.x_s_TS, float)[t] - x0[g]))              # landward shift since 1996, cells
        lo, hi = max(k, 0), min(k + col.shape[0], n)
        if hi > lo:
            z[lo:hi] = col[lo - k:hi - k]
        out.append(z * 10 + H.MHW_NAVD88_M)
    return np.concatenate(out, axis=1)


def main(key="faceretreat_equal", cs=3):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import HAT_collision_on_off_position as P
    import HAT_be_hotspot_zones as Z
    from site_layer import hat_figure_style as S
    S.apply_style()
    design = next(x for x in CD.DESIGNS if f"{x.erosion}_{x.start}" == key)
    ref = DS.reference(design)
    v = H.DESIGN_SWEEP[(key, cs)]
    lab = f"{design.name.rsplit(', Cs', 1)[0]}, Cs {DS.cs_label(cs)}"
    out = design.folder / f"{key}_cs{DS.cs_label(cs)}_visuals"
    (out / "supporting").mkdir(parents=True, exist_ok=True)

    # 1. Shoreline position change against CoastSat
    o, *_ = P.figure(1996, 2009, "HAT_dune_collision",
                     ((ref.run, ref.name + " (reference)", "#7F7F7F", 1.6, ref.linestyle),
                      (v, lab, design.colour, 2.4, "-")),
                     f"Shoreline position change, {lab}, 1996–2009",
                     str((out / "position_1996_2009.png").relative_to(P.OUT)),
                     f"the design {lab} (run {v}; coloured) against {ref.name} (run {ref.run}; grey), the same start "
                     "without the module: " + P.DESIGN_CAPTION, post=P._design_zone_shade)
    print(f"wrote {o}")

    # 2. The dune rows along the island through time
    cm, dm, first = load(v)
    cr, dr, _ = load(ref.run)
    R, Rr = rows(cm, dm), rows(cr, dr)
    years = first + np.arange(R.shape[0])
    crest, crest_ref = R.max(axis=2), Rr.max(axis=2)
    panels = ((crest_ref, f"{ref.name}: dune crest (higher row)", GREENS, 0, DUNE_MAX_M),
              (R[:, :, 0], f"{lab}: front (seaward) row", GREENS, 0, DUNE_MAX_M),
              (R[:, :, 1], f"{lab}: back row", GREENS, 0, DUNE_MAX_M),
              (crest - crest_ref, f"{lab}: crest minus {ref.name.lower()}", DIVERGING, -DIFF_M, DIFF_M))
    fig, axs = plt.subplots(4, 1, figsize=(S.FIG_W_DOUBLE, 9.0), sharex=True, layout="constrained")
    ext = (0.5, 90.5, years[-1] + 0.5, years[0] - 0.5)
    for ax, (a, title, cmap, lo, hi) in zip(axs, panels):
        im = ax.imshow(a, aspect="auto", cmap=cmap, vmin=lo, vmax=hi, extent=ext, interpolation="nearest")
        cb = fig.colorbar(im, ax=ax, pad=0.01, fraction=0.03)
        cb.set_label("m above berm" if lo == 0 else "m (blue higher)", fontsize=8)
        ax.set_title(title, loc="left", fontsize=9, fontweight="bold")
        ax.set_ylabel("Model state\n(1 Jan)", fontsize=8)
        for y, n in STORMS.items():
            ax.axhline(y + 0.3, color="0.2", lw=0.6, ls=":")
        for name, r in Z.ZONES_AVN2.items():
            for x in (min(r) - 0.5, max(r) + 0.5):
                ax.axvline(x, color="0.25", lw=0.7, ls="--")
    for y, n in STORMS.items():
        axs[0].text(90.8, y + 0.3, n, fontsize=6.5, va="center", ha="left", color="0.3", clip_on=False)
    for name, r in Z.ZONES_AVN2.items():
        axs[-1].text(np.mean(list(r)), years[0] - 0.2, name.split(",")[0], ha="center", va="top", fontsize=7,
                     color="0.15", bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    axs[-1].set_xlabel("Alongshore position (GIS domain, south → north)")
    axs[-1].set_xticks(range(0, 91, 10))
    o = out / "dune_rows_1996_2009.png"
    S.save(fig, o, bbox_inches="tight", close=True)
    S.record_caption(o, (
        f"The model's dune through 1996-2009, every 10 m cell of GIS 1-90 (x) at each model state, 1 January (y, "
        f"down = later): {ref.name}'s crest (the higher of its two rows; run {ref.run}), the front and back row of "
        f"{lab} (run {v}), and its crest minus {ref.name.lower()}'s (red lower, blue higher). Heights above the "
        "berm. Dotted: Bonnie, Dennis/Floyd, Isabel; dashed: the hotspot zones Buxton 9-14, Avon south 23-26, "
        "Rodanthe 82-89."))
    print(f"wrote {o}")
    for name, r in Z.ZONES_AVN2.items():
        i = np.concatenate([np.arange((g - 1) * 50, g * 50) for g in r])
        print(f"  {name.split(',')[0]:12s} 2010 state: front {R[-1, i, 0].mean():.2f}, back {R[-1, i, 1].mean():.2f}, "
              f"crest minus ref {(crest - crest_ref)[-1, i].mean():+.2f} m")
    i = np.arange(5 * 50, 90 * 50)
    print(f"  GIS 6-90     2010 state: front {R[-1, i, 0].mean():.2f}, back {R[-1, i, 1].mean():.2f}, "
          f"crest minus ref {(crest - crest_ref)[-1, i].mean():+.2f} m")

    # 3. Plan views of the model grid at the hotspot zones
    for tag, (zname, gis) in GRID_ZONES.items():
        x0 = {g: float(np.asarray(cr.barrier3d[dr.gis_to_pad(g)].x_s_TS, float)[0]) for g in gis}
        fig, axs = plt.subplots(3, len(GRID_YEARS), figsize=(S.FIG_W_DOUBLE, 8.8), sharex=True, sharey=True,
                                layout="constrained")
        ext = (min(gis) - 0.5, max(gis) + 0.5, CROSS_M, 0)
        for col, t in enumerate(GRID_YEARS):
            dz = plan(cm, dm, gis, t, x0) - plan(cr, dr, gis, t, x0)
            ax = axs[2, col]
            ax.set_facecolor("#DCE9F5")
            imd = ax.imshow(np.ma.masked_invalid(dz), cmap=DIVERGING, vmin=-DIFF_M, vmax=DIFF_M, aspect="auto",
                            interpolation="nearest", extent=ext)
            if col == 0:
                ax.set_ylabel("Design minus no collision\nlandward of 1996 shoreline (m)", fontsize=8)
        cbd = fig.colorbar(imd, ax=axs[2], pad=0.01, fraction=0.025)
        cbd.set_label("m (blue higher)", fontsize=8)
        for row, (c, d, who) in enumerate(((cr, dr, ref.name), (cm, dm, lab))):
            for col, t in enumerate(GRID_YEARS):
                ax = axs[row, col]
                z = plan(c, d, gis, t, x0)
                zm = np.ma.masked_where(~np.isfinite(z) | (z < H.MHW_NAVD88_M), z)
                ax.set_facecolor("#DCE9F5")                                # water / nothing: pale blue
                im = ax.imshow(zm, cmap=GREENS, vmin=0, vmax=8, aspect="auto", interpolation="nearest",
                               extent=(min(gis) - 0.5, max(gis) + 0.5, CROSS_M, 0))
                if row == 0:
                    ax.set_title(f"{first + t} (1 Jan)" + (" — after Isabel" if t == 8 else ""), fontsize=9)
                if col == 0:
                    ax.set_ylabel(f"{who}\nlandward of 1996 shoreline (m)", fontsize=8)
                ax.set_xticks(list(gis))
        for ax in axs[-1]:
            ax.set_xlabel("GIS domain")
        cb = fig.colorbar(im, ax=axs[:2], pad=0.01, fraction=0.025)
        cb.set_label("Elevation (m NAVD88)", fontsize=8)
        fig.suptitle(f"Model grid at {zname}: {ref.name} (top), {lab} (middle), difference (bottom)",
                     fontsize=10, fontweight="bold")
        o = out / f"model_grid_{tag}.png"
        S.save(fig, o, bbox_inches="tight", close=True)
        S.record_caption(o, (
            f"Plan view of the model's own 10 m grid at {zname} (GIS {min(gis)}-{max(gis)}), ocean at the top: the "
            "two dune rows (berm plus dune height) ahead of the interior, each year drawn from its own shoreline, "
            "measured landward of the 1996 shoreline, so shoreline retreat moves the island down the panel. "
            f"Model states 1 January {', '.join(str(first + t) for t in GRID_YEARS)}; top {ref.name} (run "
            f"{ref.run}), bottom {lab} (run {v}). Pale blue: below MHW or outside the domain. m NAVD88."))
        print(f"wrote {o}")


SIDE_PLACES = ((12, "Buxton"), (24, "Avon south"), (50, "Mid-island (quiet)"), (85, "Rodanthe"))
SIDE_STATES = ((0, 1996, "start (1996)"), (13, 2009, "end (1 Jan 2010)"))


def side_by_side(cs=3):
    """0-all-designs/designs_side_by_side_cs<Cs>.png: the four erosion designs on equal rows, all at one Cs, side by
    side at four places: each domain's mean front and back row (m above the berm) at the start, after Isabel and at
    the end, the lidar domain-median crest at the matching surveys as markers."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.lines import Line2D
    from site_layer import hat_figure_style as S
    S.apply_style()
    designs = [d for d in CD.DESIGNS if d.start == "equal"]
    runs = {d.erosion: (d.run if d.erosion == "nocollision" else H.DESIGN_SWEEP[(f"{d.erosion}_equal", cs)])
            for d in designs}
    sv = pd.read_csv(H.SURVEY_CRESTS)
    sv = sv[sv.status == "ok"]
    lid = sv.groupby(["gis", "survey"]).crest_m_navd88.median() - H.TOE_AT_BERM_M
    fig, axs = plt.subplots(len(SIDE_PLACES), len(designs), figsize=(S.FIG_W_DOUBLE, 2.0 * len(SIDE_PLACES) + 0.9),
                            sharex=True, sharey=True, layout="constrained")
    shades = {0: ("0.55", 1.3, (0, (3, 2))), 7: (None, 1.8, "-"), 13: (None, 2.6, "-")}
    for j, d in enumerate(designs):
        c, dm, first = load(runs[d.erosion])
        colour = d.colour
        for i, (g, place) in enumerate(SIDE_PLACES):
            ax = axs[i, j]
            dd = np.asarray(c.barrier3d[dm.gis_to_pad(g)].DuneDomain, float) * 10
            for t, yr, _ in SIDE_STATES:
                f, b = dd[t, :, 0].mean(), dd[t, :, 1].mean()
                col, lw, ls = shades[t]
                ax.plot([-6, 0, 0, 10, 10, 20, 20, 26], [0, 0, f, f, b, b, 0, 0], color=col or colour, lw=lw, ls=ls,
                        alpha=0.5 if t == 7 else 1.0, solid_joinstyle="miter")
                if (g, yr) in lid.index:
                    ax.plot(23, lid[(g, yr)], marker={0: "s", 7: "o", 13: "D"}[t], ms=5.5,
                            mfc="white" if t == 7 else ("0.55" if t == 0 else "0.1"), mec="0.1", mew=1.0, ls="none")
            f, b = dd[13, :, 0].mean(), dd[13, :, 1].mean()
            ax.text(0.03, 0.96, f"end: front {f:.1f}, back {b:.1f} m", transform=ax.transAxes, fontsize=7,
                    va="top", color="0.15")
            ax.axhline(0, color="0.6", lw=0.6)
            ax.set_xticks([5, 15])
            ax.set_xticklabels(["front", "back"], fontsize=7.5)
            ax.spines[["top", "right"]].set_visible(False)
            if i == 0:
                ax.set_title(d.name.split(",")[0], fontsize=9.5, fontweight="bold", color=colour)
            if j == 0:
                ax.set_ylabel(f"{place}, GIS {g}\nm above berm", fontsize=8)
    axs[0, 0].set_ylim(-0.2, 7.5)
    axs[0, 0].set_xlim(-6, 26)
    h = [Line2D([], [], color="0.55", lw=1.3, ls=(0, (3, 2)), label="model start (1996, both rows equal)"),
         Line2D([], [], color="0.15", lw=2.6, label="model end (1 Jan 2010)"),
         Line2D([], [], marker="s", mfc="0.55", mec="0.1", ls="none", label="lidar crest 1996"),
         Line2D([], [], marker="D", mfc="0.1", mec="0.1", ls="none", label="lidar crest 2009")]
    fig.legend(handles=h, loc="outside lower center", ncol=2, fontsize=7.5, frameon=False)
    fig.suptitle(f"How each design erodes the two dune rows, 1996 to 2009, equal rows at the start, "
                 f"Cs {DS.cs_label(cs)}", fontsize=10.5, fontweight="bold")
    o = CD.ALL_DIR / f"designs_side_by_side_cs{DS.cs_label(cs)}.png"
    S.save(fig, o, bbox_inches="tight", close=True)
    S.record_caption(o, (
        f"The two model dune rows (front seaward, back landward; each 10 m wide) as a step profile, the mean over "
        f"the domain's 50 cells, m above the berm, for each erosion design on equal rows at Cs {DS.cs_label(cs)} "
        "(columns; no collision has no Cs) at four domains (rows), at the start (both rows equal) and "
        "1 January 2010 (the end). Markers: the lidar domain-median crest at the 1996 and 2009 surveys, "
        "compared with the higher row. Buxton GIS 12 in 2004 and 2009: beside the dune rebuilt along the relocated "
        "NC-12, the picked lidar crest reads about 0.8 m low. Runs: "
        + "; ".join(f"{d.name.split(',')[0]} {runs[d.erosion]}" for d in designs) + "."))
    print(f"wrote {o}")


if __name__ == "__main__":
    a = sys.argv[1:]
    if a and a[0] == "side":
        for n in (a[1:] or ["3", "10"]):
            side_by_side(int(n))
    else:
        main(a[0] if a else "faceretreat_equal", int(a[1]) if len(a) > 1 else 3)
