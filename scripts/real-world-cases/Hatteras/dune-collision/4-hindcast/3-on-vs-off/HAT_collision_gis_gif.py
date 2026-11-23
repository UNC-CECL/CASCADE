"""
A few GIS domains through a collision-design run, year by year (2026-11-20): the model's elevation and the year's
elevation change, in plan view at true scale.

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_gis_gif.py
    python .../HAT_collision_gis_gif.py grid_lowerfront_equal_cs200_vdefault 8 12
    python .../HAT_collision_gis_gif.py cal_step2_growth_noscrape_vdefault 8 12     # no collision

One frame per model state (1 January; the first is the 1996 start). Per frame:
(a) the model's 10 m cells, m NAVD88, terrain colours, ocean at the bottom: the two dune rows (berm plus dune height)
    then the interior, each state placed where its grid lies relative to the 1996 shoreline (each year matched to
    the year before), so shoreline retreat moves the island up the panel; water below MHW in blue;
(b) the elevation change over the year that led to the state (red lower, blue higher; under 5 cm white): dune
    lowering and shoreline retreat at the front, overwash deposition behind.
Both: the 1996 shoreline dashed, GIS domain boundaries along the top.

Writes, in the design's folder (collision_designs.py), <key>_cs<Cs>_visuals/gis<a>-<b>_<y0>_<y1>.gif and the last
frame as a PNG.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-20
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(next(p for p in _HERE.parents if (p / "hatteras_paths.py").exists())))
sys.path.insert(0, str(_HERE.parent))
import collision_designs as CD  # noqa: E402
import HAT_collision_design_cs_sweep as DS  # noqa: E402
import HAT_collision_design_visuals as V  # noqa: E402

H = V.H

# --- CONFIG ------------------------------------------------------------------
DEFAULT_RUN = "grid_lowerfront_equal_cs200_vdefault"
DEFAULT_GIS = (8, 12)
PLACE = {(8, 12): "Buxton"}          # name in the title for a GIS range
CROSS_M = 500                        # drawn this far landward of the 1996 shoreline
PAD_M = 60                           # and this far seaward (ocean)
OCEAN_M, BAY_M = -2.0, -1.0          # fill seaward of the shoreline / landward of the domain, m below MHW
Z_MIN, Z_MAX = -3.0, 6.0             # elevation scale, m NAVD88
CHANGE_M = 1.5                       # +- m, change scale (extends past)
CHANGE_MIN = 0.05                    # changes smaller than this left white
FPS = 0.8
DPI = 150
STORMS = {1998: "Hurricane Bonnie (Aug 1998)", 1999: "Hurricanes Dennis and Floyd (Aug–Sep 1999)",
          2003: "Hurricane Isabel (Sep 2003)", 2006: "Nor'easter (Nov 2006)"}
NO_DATA = "#e9ecef"                  # water in the change panel
# -----------------------------------------------------------------------------


def offsets(b):
    """Per state, how many 10 m cells the domain's interior has moved landward since the start: each year matched
    to the year before by the shift that best lines the two grids up (shoreline position rounded to cells jumps a
    cell when the grid has not moved, which would read as change everywhere)."""
    D = [np.asarray(x, float) for x in b.DomainTS]
    off = [0]
    for t in range(1, len(D)):
        a, p = D[t], D[t - 1]
        best = min(range(-2, 4), key=lambda s: np.nanmean(np.abs(
            a[max(-s, 0):max(-s, 0) + 60] - p[max(s, 0):max(s, 0) + 60])))
        off.append(off[-1] + best)
    return off


def plan(b, off, t):
    """One domain at state t, (cross-shore cells, 50) m NAVD88, row 0 at PAD_M seaward of the 1996 shoreline: ocean,
    the two dune rows (berm plus height), the interior, then bay where the domain ends."""
    pad, n = PAD_M // 10, (PAD_M + CROSS_M) // 10
    z = np.full((n, 50), BAY_M / 10)
    col = np.vstack([(np.asarray(b.DuneDomain, float)[t] + b._BermEl).T, np.asarray(b.DomainTS[t], float)])
    k = pad + off[t]
    z[:max(k, 0)] = OCEAN_M / 10
    lo, hi = max(k, 0), min(k + col.shape[0], n)
    if hi > lo:
        z[lo:hi] = col[lo - k:hi - k]
    return z * 10 + H.MHW_NAVD88_M


def main(run=DEFAULT_RUN, g0=DEFAULT_GIS[0], g1=DEFAULT_GIS[1]):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import PillowWriter
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import ListedColormap, Normalize, TwoSlopeNorm
    from site_layer import hat_figure_style as S
    S.apply_style()

    gis = list(range(g0, g1 + 1))
    sweep = next((k for k, x in H.DESIGN_SWEEP.items() if x == run), None)
    if sweep:                                    # a run of the Cs sweep
        key, cs = sweep
        design = next(x for x in CD.DESIGNS if f"{x.erosion}_{x.start}" == key)
        lab = f"{design.name.rsplit(', Cs', 1)[0]}, Cs {DS.cs_label(cs)}"
        folder = f"{key}_cs{DS.cs_label(cs)}_visuals"
    else:                                        # a design's own run (no collision has no Cs)
        design = next(x for x in CD.DESIGNS if x.run == run)
        key, lab = f"{design.erosion}_{design.start}", design.name
        folder = f"{key}_visuals"
    with_module = (f"with the dune collision module: {lab}" if design.erosion != "nocollision"
                   else f"without the dune collision module ({lab.lower()})")
    c, d, first = V.load(run)
    B = {g: c.barrier3d[d.gis_to_pad(g)] for g in gis}
    n = np.asarray(B[g0].DuneDomain).shape[0]
    offs = {g: offsets(B[g]) for g in gis}
    plans = [np.concatenate([plan(B[g], offs[g], t) for g in gis], axis=1) for t in range(n)]
    mhw = H.MHW_NAVD88_M

    # matplotlib's terrain: its water blues below MHW, its land colours above, meeting at MHW
    terrain = plt.get_cmap("terrain")
    cmap_z = ListedColormap(np.vstack([terrain(np.linspace(0.0, 0.17, 128)), terrain(np.linspace(0.25, 1.0, 128))]))
    norm_z = TwoSlopeNorm(vmin=Z_MIN, vcenter=mhw, vmax=Z_MAX)
    cmap_d = plt.get_cmap("RdBu").copy()
    cmap_d.set_bad(NO_DATA)
    norm_d = Normalize(-CHANGE_M, CHANGE_M)

    L = len(gis) * 0.5                                       # km alongshore
    ext = (0, L, -PAD_M, CROSS_M)                            # km, m

    fig = plt.figure(figsize=(8.0, 5.8))
    ax_z = fig.add_axes((0.08, 0.475, 0.80, 0.31))
    ax_d = fig.add_axes((0.08, 0.125, 0.80, 0.31), sharex=ax_z)
    cax_z = fig.add_axes((0.895, 0.475, 0.013, 0.31))
    cax_d = fig.add_axes((0.895, 0.125, 0.013, 0.31))
    cz = fig.colorbar(ScalarMappable(norm_z, cmap_z), cax=cax_z, extend="both")
    cz.set_ticks([-3, -2, -1, 0, 1, 2, 4, 6])
    cz.set_label("Elevation (m NAVD88)", fontsize=8)
    cd = fig.colorbar(ScalarMappable(norm_d, cmap_d), cax=cax_d, extend="both")
    cd.set_label("Elevation change (m)", fontsize=8)
    for cb in (cz, cd):
        cb.ax.tick_params(labelsize=7)
        cb.outline.set_linewidth(0.6)
    place = PLACE.get((g0, g1))
    fig.text(0.08, 0.96, f"Modelled barrier evolution{', ' + place if place else ''} (GIS {g0}–{g1}), "
             f"{first - 1}–{first + n - 2}", fontsize=11, fontweight="bold", ha="left", va="center")
    fig.text(0.08, 0.922, f"Barrier3D/CASCADE hindcast {with_module}",
             fontsize=8.5, color="0.3", ha="left", va="center")
    fig.text(0.08, 0.035, f"Barrier3D 10 m model cells at true scale. Dashed: 1996 shoreline; dotted: GIS domain "
             f"boundaries. Blue: below mean high water ({mhw:.2f} m NAVD88).\n"
             f"Run {run}.", fontsize=6.5, color="0.4", ha="left", va="center", linespacing=1.4)
    stamp = fig.text(0.88, 0.96, "", fontsize=11, fontweight="bold", ha="right", va="center")
    event = fig.text(0.08, 0.884, "", fontsize=8.5, fontweight="bold", color="#b2182b", ha="left", va="center")

    def frame_axes(ax, letter, title, top_labels):
        ax.set_aspect(0.001)                                  # 1 m cross-shore = 1 m alongshore
        ax.set_xlim(0, L)
        ax.set_ylim(-PAD_M, CROSS_M)
        ax.axhline(0, color="0.1", lw=0.8, ls=(0, (4, 2)))
        for i in range(1, len(gis)):
            ax.axvline(i * 0.5, color="0.15", lw=0.5, ls=":", alpha=0.7)
        ax.set_ylabel("Cross-shore distance\nfrom 1996 shoreline (m)", fontsize=8)
        ax.set_yticks(range(0, CROSS_M + 1, 100))
        ax.tick_params(labelsize=7.5, direction="out", length=3)
        for s in ax.spines.values():
            s.set_linewidth(0.6)
        ax.text(0.006, 0.97, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top", ha="left",
                bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=0.9))
        ax.set_title(title, loc="left", fontsize=9, pad=3 if not top_labels else 14)
        if top_labels:
            for i, g in enumerate(gis):
                ax.text((i + 0.5) * 0.5, 1.015, f"GIS {g}", transform=ax.get_xaxis_transform(), ha="center",
                        va="bottom", fontsize=7, color="0.3")

    def draw(t):
        for ax in (ax_z, ax_d):
            ax.cla()
        yr = first + t - 1
        stamp.set_text(f"1 January {first + t}" if t else f"Start ({first - 1})")
        event.set_text(f"Storm of the year: {STORMS[yr]}" if t and yr in STORMS else "")
        z = plans[t]
        ax_z.imshow(z, cmap=cmap_z, norm=norm_z, origin="lower", extent=ext, interpolation="nearest",
                    aspect="auto")
        ax_z.text(0.5, 0.035, "Atlantic Ocean", transform=ax_z.transAxes, ha="center", va="bottom", fontsize=7.5,
                  style="italic", color="white")
        ax_z.text(0.995, 0.96, "Pamlico Sound", transform=ax_z.transAxes, ha="right", va="top", fontsize=7.5,
                  style="italic", color="white")
        frame_axes(ax_z, "a", "Elevation", True)
        ax_z.tick_params(labelbottom=False)

        if t:
            dz = plans[t] - plans[t - 1]
            land = (plans[t] >= mhw) | (plans[t - 1] >= mhw)
            dz = np.where(np.abs(dz) < CHANGE_MIN, 0.0, dz)
            ax_d.imshow(np.ma.masked_where(~land, dz), cmap=cmap_d, norm=norm_d, origin="lower", extent=ext,
                        interpolation="nearest", aspect="auto")
            ttl = f"Elevation change over {yr}"
        else:
            ax_d.set_facecolor(NO_DATA)
            ax_d.text(0.5, 0.5, "Initial state", transform=ax_d.transAxes, ha="center", va="center",
                      fontsize=8.5, color="0.4")
            ttl = "Elevation change over the year"
        frame_axes(ax_d, "b", ttl, False)
        ax_d.set_xlabel(f"Alongshore distance from the south edge of GIS {g0} (km, south → north)", fontsize=8)

    tag = f"gis{g0:02d}-{g1:02d}_{first - 1}_{first + n - 2}"
    out = design.folder / folder
    out.mkdir(parents=True, exist_ok=True)
    gif = out / f"{tag}.gif"
    w = PillowWriter(fps=FPS)
    with w.saving(fig, str(gif), dpi=DPI):
        for t in range(n):
            draw(t)
            w.grab_frame()
        for _ in range(2):                       # hold the last frame
            w.grab_frame()
    fig.savefig(out / f"{tag}_last_frame.png", dpi=200)
    plt.close(fig)
    S.record_caption(gif, (
        f"GIS {g0}-{g1} through {lab} (run {run}), one frame per model state (1 January; the first the "
        f"{first - 1} start), plan view at true scale, ocean at the bottom. (a) The model's 10 m grid, m NAVD88, "
        "matplotlib terrain colours, each model cell drawn as a block: the two dune rows (berm plus dune height) then the interior, each state "
        "placed where its grid lies relative to the 1996 shoreline (each year matched to the year before), so "
        f"shoreline retreat moves the island up the panel; water below MHW blue (ocean drawn at {OCEAN_M:g} m, bay "
        f"beyond the domain at {BAY_M:g} m, relative to MHW). (b) The elevation change over the year before the "
        f"state, red lower, blue higher, under {CHANGE_MIN * 100:.0f} cm white; grey water. Dashed: the 1996 "
        "shoreline; dotted: GIS domain boundaries."))
    print(f"wrote {gif}")


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0] if a else DEFAULT_RUN, *(int(x) for x in a[1:3]))
