"""
Did the dune line move landward with the shoreline in 1996-2009, inside the villages as outside? The lidar's answer.

    python scripts/real-world-cases/Hatteras/dune-collision/1-observations/8-dune-change-by-survey/HAT_dune_line_position.py [plot]

The question behind it (2026-11-12): CASCADE holds a village's dune line fixed (BeachDuneManager, migration off until
the beach is gone), Barrier3D moves it with the shoreline everywhere else. Which did the real coast do?

Every survey is cut into 10 m cross-shore profiles in the 1996 start's frame (HAT_survey_dune_crests.py: the start's
shear and water cut-off, so a profile index is the same place on the ground in every survey), and the input step's
picker (HAT_dune_toes.crest_and_toe, as the collision runner's toes) finds, per 10 m alongshore cell, the beach start
(the shoreline: the first ground above the beach threshold), the dune toe and the crest, the crest searched in the
survey's dune window (moved landward with the dune, as the crest script did). The change of each from the domain's
first survey to its last, per cell, positive landward, gives how far the shoreline and the dune line moved.

Writes, in the workspace's dune-collision/1-observations/8-dune-change-by-survey/ (D:):
- supporting/dune_line_positions.csv: per survey, domain and cell, the three positions (m from the frame's seaward edge)
- supporting/dune_line_retreat_by_domain.csv: per domain, the median retreat of each, first survey to last
- dune_line_retreat.png: (a) alongshore, the shoreline, toe and crest retreat, villages shaded; (b) the dune line's
  retreat against the shoreline's, villages against the open coast

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-12
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_HAT_DIR = next(p for p in _HERE.parents if (p / "hatteras_paths.py").exists())
for _p in (_HAT_DIR, _HAT_DIR / "input_prep" / "1-barrier3d-domains" / "2-survey-crests",
           _HAT_DIR / "input_prep" / "1-barrier3d-domains" / "3-dune-toes"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
from hatteras_paths import WORKSPACE  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = WORKSPACE / "dune-collision" / "1-observations" / "8-dune-change-by-survey"
START = "1996-start"
CELL_M = 10.0
VILLAGES = {"Buxton": (7, 8), "Avon": (21, 31), "Tri-Village": (68, 83)}    # the site config's community spans
EXCLUDE = (1, 5, 6, 7, 90)                     # the locked ends and the groin, as the zone rules
# -----------------------------------------------------------------------------


def measure():
    import rasterio
    import HAT_survey_dune_crests as S
    import HAT_dune_toes as T
    X, tv = S.X, S.tv
    prod, surveys = S.STARTS[START]
    windows = json.loads(tv.picks_dir(START).joinpath(
        f"HAT_dune_search_windows_{tv.topo_dirs(START)[2]}.json").read_text(encoding="utf-8"))
    frames = {g: S.frame(START, g, windows) for g in range(1, 91)}
    clips = S.product(prod).gapfill_1m
    grids = {g: S.clip_grid(clips / f"clip_domain_{g}_filled.tif") for g in range(1, 91)}
    offsets = S.SURVEY_VERTICAL_OFFSET_M.get(START, {})
    rows = []
    for year, (kind, path, flown, state, note) in surveys.items():
        files = sorted(path.parent.glob(path.name)) if "*" in path.name else [path]
        srcs = [rasterio.open(f) for f in files]
        try:
            unit = "m" if kind == "mosaic" else S.unit_of(srcs, START, grids)[0]
            scale = S.US_FT if unit == "ft" else 1.0
            for g in range(1, 91):
                shear, c0, i0, i1 = frames[g]
                on_clip = S.mosaic_on_clip(srcs[0], grids[g]) if kind == "mosaic" else S.raster_on_clip(srcs, grids[g])
                z = S.profiles(S.downsample(on_clip * scale + offsets.get(year, 0.0)), shear, c0)
                r = T.crest_and_toe(z, *S.survey_windows(START, year, g, i0, i1), X)
                for cell in range(50):
                    rows.append(dict(survey=int(year), gis=g, cell=cell,
                                     shoreline_m=r["beach_start_idx"][cell] * CELL_M,
                                     toe_m=r["toe_idx"][cell] * CELL_M, crest_m=r["crest_idx"][cell] * CELL_M,
                                     crest_m_mhw=r["crest"][cell]))
        finally:
            for s in srcs:
                s.close()
        print(f"{year} done", flush=True)
    t = pd.DataFrame(rows)
    # the crest script's screen (spikes, unsurveyed cells)
    sc = pd.read_csv(S.tv.product_dir(START) / "survey-crests" / "survey_dune_crests.csv")[["survey", "gis", "cell",
                                                                                              "status"]]
    t = t.merge(sc, on=["survey", "gis", "cell"], how="left")
    (OUT / "supporting").mkdir(parents=True, exist_ok=True)
    t.round(2).to_csv(OUT / "supporting" / "dune_line_positions.csv", index=False)
    return t


def village_of(g):
    return next((n for n, (lo, hi) in VILLAGES.items() if lo <= g <= hi), "")


def by_domain(t):
    ok = t[(t.status == "ok") & t[["shoreline_m", "toe_m", "crest_m"]].notna().all(axis=1)]
    rows = []
    for g, d in ok.groupby("gis"):
        ys = sorted(d.survey.unique())
        if len(ys) < 2 or g in EXCLUDE:
            continue
        a, b = d[d.survey == ys[0]].set_index("cell"), d[d.survey == ys[-1]].set_index("cell")
        cells = a.index.intersection(b.index)
        if len(cells) < 10:
            continue
        diff = (b.loc[cells, ["shoreline_m", "toe_m", "crest_m"]] - a.loc[cells, ["shoreline_m", "toe_m", "crest_m"]])
        rows.append(dict(gis=g, village=village_of(g), first=ys[0], last=ys[-1], cells=len(cells),
                         shoreline_retreat_m=diff.shoreline_m.median(), toe_retreat_m=diff.toe_m.median(),
                         crest_retreat_m=diff.crest_m.median(),
                         dune_line_retreat_m=diff[["toe_m", "crest_m"]].median(axis=1).median()))
    m = pd.DataFrame(rows).set_index("gis")
    m.round(1).to_csv(OUT / "supporting" / "dune_line_retreat_by_domain.csv")
    return m


def figure(m):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as S
    S.apply_style()
    fig = plt.figure(figsize=(S.FIG_W_DOUBLE, 6.2), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=(1.0, 1.15), width_ratios=(1.25, 1.0))
    ax = fig.add_subplot(gs[0, :])
    for name, (lo, hi) in VILLAGES.items():
        ax.axvspan(lo - 0.5, hi + 0.5, color="0.90", lw=0, zorder=0)
        ax.text((lo + hi) / 2, 0.97, name, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7.5,
                color=S.INK_MUTED)
    ax.axhline(0, color="0.6", lw=0.6)
    g = m.index.to_numpy()
    ax.plot(g, m.shoreline_retreat_m, color="#D7301F", lw=1.6, ls="--", marker="o", ms=2.5, label="Shoreline (beach start)")
    ax.plot(g, m.toe_retreat_m, color="#2171B5", lw=1.3, marker="o", ms=2.5, label="Dune toe")
    ax.plot(g, m.crest_retreat_m, color="#1A1A1A", lw=1.3, marker="o", ms=2.5, label="Dune crest")
    ax.set_xlim(0.5, 90.5)
    ax.set_xticks([1] + list(range(10, 91, 10)))
    ax.set_xlabel(S.DOMAIN_AXIS_LABEL)
    ax.set_ylabel("Landward movement,\nfirst survey to 2009 (m)")
    ax.grid(axis="y", color="0.92", lw=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper left", ncol=3, frameon=False, fontsize=7.5, bbox_to_anchor=(0, 0.9))
    ax.set_title("(a) How far each moved landward, per domain (median of its 10 m cells)", loc="left", fontsize=9,
                 fontweight="bold")
    axs = fig.add_subplot(gs[1, 0])
    lim = float(np.nanmax(np.abs(m[["shoreline_retreat_m", "dune_line_retreat_m"]].to_numpy()))) + 5
    axs.plot([-lim, lim], [-lim, lim], color="0.6", lw=0.8, ls=":", label="Dune line moves with the shoreline (1:1)")
    axs.axhline(0, color="0.8", lw=0.6)
    axs.axvline(0, color="0.8", lw=0.6)
    for vill, kw in (("", dict(color="0.55", label="Open coast")),
                     ("Avon", dict(color="#1A9850", label="Avon")), ("Tri-Village", dict(color="#2171B5",
                                                                                         label="Tri-Village")),
                     ("Buxton", dict(color="#E66101", label="Buxton village"))):
        s = m[m.village == vill]
        axs.scatter(s.shoreline_retreat_m, s.dune_line_retreat_m, s=18, edgecolor="white", lw=0.4, zorder=3, **kw)
    axs.set_xlim(-lim, lim)
    axs.set_ylim(-lim, lim)
    axs.set_aspect("equal")
    axs.set_xlabel("Shoreline moved landward (m)")
    axs.set_ylabel("Dune line moved landward (m)\n(median of toe and crest)")
    axs.legend(loc="upper left", fontsize=6.5, frameon=False)
    axs.set_title("(b) Dune line against shoreline", loc="left", fontsize=9, fontweight="bold")
    axs.spines[["top", "right"]].set_visible(False)
    axt = fig.add_subplot(gs[1, 1])
    axt.axis("off")
    rows = []
    for name, sel in (("Open coast", m.village == ""), ("Avon", m.village == "Avon"),
                      ("Tri-Village", m.village == "Tri-Village")):
        s = m[sel]
        ratio = (s.dune_line_retreat_m / s.shoreline_retreat_m).where(s.shoreline_retreat_m > 10)
        rows.append([name, len(s), f"{s.shoreline_retreat_m.median():.0f}", f"{s.dune_line_retreat_m.median():.0f}",
                     f"{ratio.median():.2f}" if ratio.notna().any() else "–"])
    tb = axt.table(cellText=rows, colLabels=["", "Domains", "Shoreline\n(m)", "Dune line\n(m)", "Dune ÷ shore\n(>10 m)"],
                   loc="center", cellLoc="center")
    tb.auto_set_font_size(False)
    tb.set_fontsize(7)
    tb.auto_set_column_width(list(range(5)))
    tb.scale(1, 1.6)
    axt.set_title("(c) Medians, first survey to 2009", loc="left", fontsize=9, fontweight="bold")
    out = OUT / "dune_line_retreat.png"
    S.save(fig, out, dpi=300, facecolor="white", close=True)
    S.record_caption(out, (
        "Did the dune line move landward with the shoreline, 1996-2009? From the lidar: each survey cut into 10 m "
        "cross-shore profiles in the 1996 start's frame (the same place on the ground in every survey), the beach "
        "start (the shoreline), dune toe and crest picked per 10 m alongshore cell by the input step's picker "
        "(HAT_dune_toes.crest_and_toe), the crest searched in the survey's dune window; each domain's median change "
        "from its first survey to its last (2009), positive landward, 10 m resolution. (a) Along the island; grey "
        "bands, the villages (Buxton GIS 7-8, Avon 21-31, Salvo to Rodanthe 68-83). (b) The dune line (the median of "
        "the toe's and the crest's movement) against the shoreline's, per domain; dotted, the dune line moving with "
        "the shoreline. (c) Medians by stretch; the ratio over domains whose shoreline moved more than 10 m. GIS 1, "
        "5-7 and 90 left out. In the villages the crest pick can land on a rooftop (Avon south 1996-2004, GIS "
        "82-83), which moves it."))
    print(pd.DataFrame(rows, columns=["stretch", "n", "shoreline", "dune_line", "ratio"]).to_string(index=False))
    print(f"wrote {out}")


def main():
    if "plot" in sys.argv[1:]:
        t = pd.read_csv(OUT / "supporting" / "dune_line_positions.csv")
    else:
        t = measure()
    m = by_domain(t)
    print(m[m.village != ""].round(0).to_string())
    figure(m)


if __name__ == "__main__":
    main()
