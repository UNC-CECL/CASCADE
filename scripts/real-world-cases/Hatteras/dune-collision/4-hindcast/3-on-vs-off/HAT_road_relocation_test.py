"""
The 1999 NC-12 relocation switched on (2026-11-11, test). After Hurricane Dennis, NC-12 north of Buxton was moved
landward, in place by March 2000 (NCDOT 2015); the site config's RelocationEvent moves the road at GIS 9-14 by the
measured 20-80 m (road_relocation_1978_2008.csv, rounded to whole cells). The runs leave it off. Does moving the road
change the modelled dunes?

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_road_relocation_test.py

Runs (HAT_dune_collision.py, 1996-2009): cal_step2_growth_noscrape (module off), cal_toeberm_cs1 and cal_toeberm_cs2
(module on, toe at the berm), each with and without the relocation (the _reloc variants, HAT_RELOCATIONS=1); nothing
else differs. In the model a relocation only moves the road: its setback grows by the displacement, and the dune
rebuild rule (off where the road is 20 m or more behind the dune) and the overwash bulldozed off the road follow from
there. No dune is removed or built, as one was in 2000-2004.

Writes, in the workspace's dune-collision/4-hindcast/5-road/road-relocation-1999/ (D:):
- relocation_crests_over_time_GIS9-15.png: the domain-median crest through time, with and without the relocation
- relocation_difference_by_domain.png: relocation minus none, GIS 5-20: 2009 crest, overwash, shoreline change
- supporting/relocation_road_by_domain.csv: road setback and dune rebuilds per domain and year
- supporting/relocation_metrics.csv: crest RMSE, overwash and shoreline skill per run

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-11
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_collision_on_off_position as P  # noqa: E402
import HAT_collision_row_mode_sweep as R  # noqa: E402

H = R.H

# --- CONFIG ------------------------------------------------------------------
OUT = P.OUT / "5-road" / "road-relocation-1999"
PAIRS = (("cal_step2_growth_noscrape", "Module off", "#8C8C8C"),
         ("cal_toeberm_cs1", "Module on, Cs 1e-6", "#2171B5"),
         ("cal_toeberm_cs2", "Module on, Cs 2e-6", "#7B3F99"))
RELOCATED_GIS = range(9, 15)                        # the 1999 event (HATTERAS_ROAD_EVENTS)
TS_GIS = range(9, 16)
DIFF_GIS = range(5, 21)
EVENT_YEAR = 1999.67                                # Dennis; the road moved by March 2000
# -----------------------------------------------------------------------------


def _load(v):
    from cascade_pipeline.saved_runs import load_saved_cascade
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    rd = H.run_dir(v)
    c = load_saved_cascade(rd / f"{rd.name}.npz")
    return rd, c, DOM.for_total(len(c.barrier3d))


def _crest(b):
    """Each year's crest per cell, m NAVD88 (rows of DuneDomain: years; the higher of the two dune rows)."""
    return (np.asarray(b.DuneDomain, float).max(axis=2) + b._BermEl) * 10 + H.MHW_NAVD88_M


def road_table(runs):
    rows = []
    for v in runs:
        _, c, d = runs[v]
        for g in DIFF_GIS:
            r = c.roadways[d.gis_to_pad(g)]
            if r is None:
                continue
            for k in range(len(r._road_setback_TS) - 1):
                rows.append(dict(variant=v, gis=g, year=H.RUN_START + 1 + k, setback_m=r._road_setback_TS[k],
                                 dunes_rebuilt=int(r._dunes_rebuilt_TS[k]),
                                 rebuild_off=int(r._road_dune_rebuild_disabled_TS[k]),
                                 relocated=int(r._road_relocated_TS[k])))
    return pd.DataFrame(rows)


def crests_over_time(S, plt, runs):
    sv = pd.read_csv(H.SURVEY_CRESTS)
    sv = sv[(sv.status == "ok") & sv.gis.isin(TS_GIS)].copy()
    t = pd.to_datetime(sv.flown.astype(str).str[:7], errors="coerce")
    sv["yr"] = t.dt.year + (t.dt.month - 0.5) / 12
    lid = sv.groupby(["gis", "survey"]).agg(z=("crest_m_navd88", "median"), yr=("yr", "first")).reset_index()
    fig, axs = plt.subplots(len(TS_GIS), 1, figsize=(S.FIG_W_DOUBLE * 0.75, 10.6), sharex=True, sharey=True,
                            layout="constrained")
    for base, lab, colour in PAIRS:
        for v, ls, suffix in ((base, "-", ""), (base + "_reloc", (0, (4, 2)), ", road relocated 1999")):
            _, c, d = runs[v]
            for ax, g in zip(axs, TS_GIS):
                z = _crest(c.barrier3d[d.gis_to_pad(g)])
                ax.plot(H.RUN_START + 1 + np.arange(z.shape[0]), np.median(z, axis=1), color=colour, lw=1.6,
                        ls=ls, label=lab + suffix)
    for ax, g in zip(axs, TS_GIS):
        L = lid[lid.gis == g].sort_values("yr")
        ax.plot(L.yr, L.z, "o", color="#1A1A1A", ms=4.5, label="Lidar (median)", zorder=5)
        ax.axhline(H.TOE_AT_BERM_M, color="0.55", lw=0.8, label="Berm, 1.7 m (the toe)")
        ax.axvline(EVENT_YEAR, color="#AE017E", lw=0.8, ls=":")
        moved = " (road relocated)" if g in RELOCATED_GIS else ""
        ax.set_title(f"GIS {g}{moved}", loc="left", fontsize=9, fontweight="bold")
        ax.set_ylabel("Crest (m)", fontsize=8)
        ax.grid(axis="y", color="0.92", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    axs[0].text(EVENT_YEAR, 1.03, "Dennis; NC-12 relocated", transform=axs[0].get_xaxis_transform(), ha="center",
                fontsize=7, color="#AE017E")
    axs[-1].set_xlabel("Year")
    axs[-1].set_xlim(1996, 2010)
    axs[0].set_ylim(0, 9)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=2, fontsize=7.5, frameon=False)
    out = OUT / "relocation_crests_over_time_GIS9-15.png"
    S.save(fig, out, bbox_inches="tight", close=True)
    S.record_caption(out, (
        "Domain-median dune crest (m NAVD88) through 1996-2009, GIS 9-15, with NC-12 left on its old line (solid) and "
        "with the 1999 relocation (dashed; GIS 9-14 moved 20-80 m landward), for the collision module off and on at "
        "Cs 1e-6 and 2e-6 (toe at the berm), against the lidar domain median at each survey. Model state at 1 January "
        "each year. Magenta dotted: Hurricane Dennis, after which the road was moved (in place by March 2000)."))
    print(f"wrote {out}")


def differences(S, plt, runs):
    rows = []
    for base, lab, _ in PAIRS:
        per = {}
        for v in (base, base + "_reloc"):
            rd, c, d = runs[v]
            rate = pd.read_csv(rd / "tables" / "shoreline_change_rate.csv").set_index("gis_domain").lrr_m_yr
            ow = np.load(next(rd.glob("*_overwash.npy")))
            per[v] = {g: dict(crest=float(np.median(_crest(c.barrier3d[d.gis_to_pad(g)])[-1])),
                              overwash=float(ow[1:, d.gis_to_pad(g)].sum()), lrr=float(rate[g])) for g in DIFF_GIS}
        for g in DIFF_GIS:
            a, b = per[base][g], per[base + "_reloc"][g]
            rows.append(dict(run=lab, gis=g, **{f"{k}_none": a[k] for k in a}, **{f"{k}_reloc": b[k] for k in b},
                             **{f"d_{k}": b[k] - a[k] for k in a}))
    t = pd.DataFrame(rows)
    t.round(4).to_csv(OUT / "supporting" / "relocation_difference_by_domain.csv", index=False)
    panels = (("d_crest", "2009 crest (m)"), ("d_overwash", "Overwash, 1997–2009 (model units)"),
              ("d_lrr", "Shoreline change (m/yr)"))
    fig, axs = plt.subplots(len(panels), 1, figsize=(S.FIG_W_DOUBLE * 0.75, 6.2), sharex=True, layout="constrained")
    for ax, (col, lab) in zip(axs, panels):
        ax.axvspan(RELOCATED_GIS[0] - 0.5, RELOCATED_GIS[-1] + 0.5, color="#AE017E", alpha=0.08, lw=0)
        ax.axhline(0, color="0.3", lw=0.6)
        for _, name, colour in PAIRS:
            s = t[t.run == name]
            ax.plot(s.gis, s[col], "o-", color=colour, lw=1.4, ms=4, label=name)
        ax.set_ylabel(f"Relocated − not\n{lab}", fontsize=8)
        ax.grid(axis="y", color="0.92", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    axs[0].text(np.mean(RELOCATED_GIS), 1.02, "road relocated, GIS 9–14", transform=axs[0].get_xaxis_transform(),
                ha="center", fontsize=7, color="#AE017E")
    axs[-1].set_xlabel("GIS domain")
    axs[-1].set_xticks(list(DIFF_GIS))
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=3, fontsize=7.5, frameon=False)
    out = OUT / "relocation_difference_by_domain.png"
    S.save(fig, out, bbox_inches="tight", close=True)
    S.record_caption(out, (
        "The 1999 NC-12 relocation's effect on the model, GIS 5-20: the run with the relocation minus the run without, "
        "for the collision module off and on (Cs 1e-6, 2e-6; toe at the berm): the domain-median crest in the last "
        "model state (2009), the overwash summed over 1997-2009, and the shoreline change rate. Shaded: the domains "
        "where the road was moved (GIS 9-14)."))
    print(f"wrote {out}")
    return t


NOW_OUT = OUT / "20m-trigger"
NOW_RUNS = (("cal_step2_growth_noscrape_vdefault", "No collision, road not relocated", "#D7301F", 1.8, "-"),
            ("cal_step2_growth_noscrape_vdefault_reloc", "No collision, road relocated 1999", "#D7301F", 1.8,
             (0, (4, 2))))
NOW_GIS = range(8, 16)
# The Buxton dunes flattened to the 1999 lidar at GIS 10-12, without and with the ridge in front of the relocated road
# (HAT_dune_collision.py BUXTON_FLATTEN, 2026-11-20); relocation only for reference. ROYGBIV, CoastSat keeps blue
FLAT_RUNS = (("cal_step2_growth_noscrape_vdefault_reloc", "Road relocated (reference)", "#D7301F", 1.8, (0, (4, 2))),
             ("cal_step2_growth_noscrape_vdefault_reloc_flat", "+ GIS 10-12 flattened to the 1999 lidar", "#F28E2B",
              1.8, "-"),
             ("cal_step2_growth_noscrape_vdefault_reloc_flat_ridge", "+ flattened, ridge before the road 2001",
              "#1A9850", 1.8, "-"))
INDEPENDENT_CRESTS = (P.OUT.parent / "1-observations" / "8-dune-change-by-survey" / "supporting"
                      / "independent_crest_by_domain.csv")


def _zoom(ax):
    """The island-wide position figure zoomed to GIS 8-15: the island's labels outside it clipped away."""
    ax.set_xlim(7.5, 15.5)
    ax.set_ylim(-60, 40)
    ax.set_xticks(list(NOW_GIS))
    ax.axvspan(8.5, 14.5, color="#AE017E", alpha=0.06, lw=0, zorder=0)
    for a in (*ax.texts, *ax.lines, *ax.collections, *ax.patches):
        a.set_clip_on(True)
        a.set_clip_box(ax.bbox)
    for t in ax.texts:
        x, _ = t.get_position()
        if t.get_transform() is not ax.transData or not 7.5 <= x <= 15.5:
            t.set_visible(False)
    ax.text(11.5, 0.97, "road relocated, GIS 9–14", transform=ax.get_xaxis_transform(), ha="center", va="top",
            fontsize=8, color="#AE017E")
    fig = ax.figure
    fig.set_size_inches(7.0, fig.get_size_inches()[1])


NOW_TEXT = dict(
    tag="relocation_20m", title="Shoreline position change north of Buxton, NC-12 relocated or not, 1996–2009",
    runs_text="NC-12 on its 1996 line (solid) or moved 20-80 m landward at GIS 9-14 in 1999 (dashed)",
    rebuild_run=1)
FLAT_TEXT = dict(
    tag="flatten_ridge", title="Shoreline position change north of Buxton, dunes flattened and rebuilt, 1996–2009",
    runs_text="NC-12 moved 20-80 m landward at GIS 9-14 in 1999 in all three: relocation only (red dashed); with the "
              "GIS 10-12 dunes lowered to the 1999 post-Dennis lidar crest after the 1999 storms (orange); and with a "
              "20 m ridge also raised just seaward of the road at the start of 2001 to the 2004 lidar's independent "
              "crest, 5.9-6.5 m (green)",
    rebuild_run=2, indep=True)
# The 1999 relocation as a dune-line shift (HAT_dune_collision.py DUNE_SHIFT, 2026-11-23, advisor's idea): front
# interior rows removed, the dune rebuilt at the 2004 crest, the removed width banked as beach; relocation only for
# reference
SHIFT_RUNS = (("cal_step2_growth_noscrape_vdefault_reloc", "Road relocated (reference)", "#D7301F", 1.8, (0, (4, 2))),
              ("cal_step2_growth_noscrape_vdefault_reloc_shift", "+ dune line shifted landward, GIS 9-12", "#F28E2B",
               1.8, "-"),
              ("cal_step2_growth_noscrape_vdefault_reloc_shift14", "+ dune line shifted landward, GIS 9-14", "#6A3D9A",
               1.8, (0, (1, 1))))
SHIFT_TEXT = dict(
    tag="dune_shift", title="Shoreline position change north of Buxton, dune line moved with the road, 1996–2009",
    runs_text="NC-12 moved 20-80 m landward at GIS 9-14 in 1999 in all three: relocation only (red dashed); with the "
              "dune line moved landward by the road's displacement after the 1999 storms at GIS 9-12 (orange) or "
              "9-14 (purple dotted): the front interior rows removed, both dune rows set to the 2004 lidar's "
              "independent crest, and the removed width added to the beach, which pins the dune until the shoreline "
              "has eroded it",
    rebuild_run=1, indep=True)


# The flattening and the shift together (HAT_dune_collision.py FLAT_SHIFT, 2026-11-23): GIS 10-12 flattened to the 1999
# lidar, the GIS 9-12 dune line moved with the road at the start of 2002; against the shift only and the flatten only
FLATSHIFT_RUNS = (
    ("cal_step2_growth_noscrape_vdefault_reloc", "Road relocated (reference)", "#D7301F", 1.8, (0, (4, 2))),
    ("cal_step2_growth_noscrape_vdefault_reloc_shift", "+ dune line shifted 2000, GIS 9-12", "#F28E2B", 1.8, "-"),
    ("cal_step2_growth_noscrape_vdefault_reloc_flat", "+ GIS 10-12 flattened to the 1999 lidar, never rebuilt",
     "#1A9850", 1.8, (0, (1, 1))),
    ("cal_step2_growth_noscrape_vdefault_reloc_flat_shift", "+ flattened 2000, dune line shifted 2002", "#6A3D9A",
     2.2, "-"))
FLATSHIFT_TEXT = dict(
    tag="flatten_shift", title="Shoreline position change north of Buxton, dunes flattened then moved, 1996–2009",
    runs_text="NC-12 moved 20-80 m landward at GIS 9-14 in 1999 in all four: relocation only (red dashed); the GIS "
              "9-12 dune line moved with the road after the 1999 storms (orange: front interior rows removed, both "
              "dune rows set to the 2004 lidar's independent crest, the removed width added to the beach); the GIS "
              "10-12 dunes lowered to the 1999 post-Dennis lidar crest after the 1999 storms and left (green dotted); "
              "and lowered the same way, then after the 2001 storms the GIS 9-12 dune line moved as in orange, by as "
              "many rows as put the road back at its pre-relocation setback (purple)",
    rebuild_run=3, indep=True)


def shift_designs():
    """The GIS 9-12 dune-line shift on every collision design (HAT_dune_collision.py SHIFT_DESIGNS, 2026-11-23): which
    design matches the CoastSat change best? Writes NOW_OUT/dune_shift_designs_position_GIS8-15.png and _island.png,
    and supporting/dune_shift_designs_score.csv (model minus CoastSat domain median, 1996-2009)."""
    import collision_designs as CD
    from site_layer.hat_observed_rates import net_change_dir
    runs = []
    for d in CD.DESIGNS:
        v = ("cal_step2_growth_noscrape_vdefault_reloc_shift" if d.run == "cal_step2_growth_noscrape_vdefault"
             else d.run + "_reloc_shift")
        runs.append((v, d.name, d.colour, 1.6, d.linestyle))
    caption = ("no other change; NC-12 moved 20-80 m landward at GIS 9-14 in 1999 and, at GIS 9-12, the dune line "
               "moved with it after the 1999 storms (front interior rows removed by the road's displacement, both "
               "dune rows set to the 2004 lidar's independent crest, the removed width added to the beach) in every "
               "run; one line per collision design, each with its own Cs (colour: erosion design; dashed: low back "
               "row start); the road rebuilding dunes only where it is within 20 m behind them.")
    title = "Shoreline position change, dune line moved with the road, every collision design, 1996–2009"
    out, gis, ch = P.figure(1996, 2009, "HAT_dune_collision", runs, title,
                            str((NOW_OUT / "dune_shift_designs_position_GIS8-15.png").relative_to(P.OUT)),
                            caption + " Zoomed to GIS 8-15; shaded: the relocated domains.", post=_zoom)
    print(f"wrote {out}")
    out, *_ = P.figure(1996, 2009, "HAT_dune_collision", runs, title,
                       str((NOW_OUT / "dune_shift_designs_position_island.png").relative_to(P.OUT)), caption)
    print(f"wrote {out}")
    tn = pd.read_csv(net_change_dir(1996, 2009) / "transect_net_change_1996_2009.csv")
    obs = tn.groupby("domain_number").net_change_m.median()
    obs = obs.reindex(gis).to_numpy()
    rows = []
    for v, name, *_ in runs:
        err = np.asarray(ch[v]) - obs
        b = np.isin(gis, list(range(9, 15)))
        ok = np.isfinite(err)
        rows.append(dict(design=name, run=v,
                         **{f"GIS{g}_m": round(float(np.asarray(ch[v])[list(gis).index(g)]), 1) for g in range(9, 15)},
                         buxton_mean_err_m=round(float(np.nanmean(err[b])), 1),
                         buxton_mae_m=round(float(np.nanmean(np.abs(err[b]))), 1),
                         island_rmse_m=round(float(np.sqrt(np.mean(err[ok] ** 2))), 1)))
    rows.append(dict(design="CoastSat (domain median)", run="",
                     **{f"GIS{g}_m": round(float(obs[list(gis).index(g)]), 1) for g in range(9, 15)}))
    sc = pd.DataFrame(rows)
    (NOW_OUT / "supporting").mkdir(parents=True, exist_ok=True)
    sc.to_csv(NOW_OUT / "supporting" / "dune_shift_designs_score.csv", index=False)
    print(sc.drop(columns="run").to_string(index=False))


def now(runs=NOW_RUNS, text=NOW_TEXT):
    """No collision on today's default (village migration, CASCADE's 20 m road-to-dune rebuild trigger), with and
    without the 1999 relocation (2026-11-20): does the road's dune management hold back the Buxton erosion? Writes
    NOW_OUT/<tag>_position_GIS8-15.png, <tag>_crests_GIS8-15.png, supporting/<tag>_road.csv. With FLAT_RUNS/FLAT_TEXT
    (action flat): the dunes flattened, with and without the rebuilt ridge"""
    NOW_RUNS = runs
    tag = text["tag"]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as S
    S.apply_style()
    (NOW_OUT / "supporting").mkdir(parents=True, exist_ok=True)
    # Position change against CoastSat, zoomed to GIS 8-15
    out, gis, ch = P.figure(
        1996, 2009, "HAT_dune_collision", [r[:4] + (r[4],) for r in NOW_RUNS],
        text["title"],
        str((NOW_OUT / f"{tag}_position_GIS8-15.png").relative_to(P.OUT)),
        "no collision module, growth 0.02-0.04, groin 0.16, village dune lines migrating until first nourished, "
        "Rodanthe 40%, the road rebuilding dunes only where it is within 20 m behind them; " + text["runs_text"]
        + ". Zoomed to GIS 8-15; shaded: the relocated domains.",
        post=_zoom)
    print(f"wrote {out}")
    for g in NOW_GIS:
        i = list(gis).index(g)
        print(f"  GIS {g}: " + ", ".join(f"{lab} {np.asarray(ch[v])[i]:+6.1f} m" for v, lab, *_ in NOW_RUNS))
    # Road management and crests through time
    runs = {r[0]: _load(r[0]) for r in NOW_RUNS}
    global DIFF_GIS
    DIFF_GIS = NOW_GIS
    road = road_table(runs)
    road.to_csv(NOW_OUT / "supporting" / f"{tag}_road.csv", index=False)
    print(road.groupby(["variant", "gis"]).agg(setback_1997=("setback_m", "first"), setback_2009=("setback_m", "last"),
                                               rebuilds=("dunes_rebuilt", "sum")).reset_index().to_string(index=False))
    sv = pd.read_csv(H.SURVEY_CRESTS)
    sv = sv[(sv.status == "ok") & sv.gis.isin(NOW_GIS)].copy()
    t = pd.to_datetime(sv.flown.astype(str).str[:7], errors="coerce")
    sv["yr"] = t.dt.year + (t.dt.month - 0.5) / 12
    lid = sv.groupby(["gis", "survey"]).agg(z=("crest_m_navd88", "median"), yr=("yr", "first")).reset_index()
    indep = pd.read_csv(INDEPENDENT_CRESTS)
    fig, axs = plt.subplots(4, 2, figsize=(S.FIG_W_DOUBLE, 7.6), sharex=True, sharey=True, layout="constrained")
    for ax, g in zip(axs.T.ravel(), NOW_GIS):
        for v, lab, colour, lw, ls in NOW_RUNS:
            _, c, d = runs[v]
            z = _crest(c.barrier3d[d.gis_to_pad(g)])
            ax.plot(H.RUN_START + 1 + np.arange(z.shape[0]), np.median(z, axis=1), color=colour, lw=lw, ls=ls,
                    label=lab)
        rb = road[(road.variant == NOW_RUNS[text["rebuild_run"]][0]) & (road.gis == g) & (road.dunes_rebuilt > 0)].year
        for y in rb:
            ax.axvline(y, color="#D7301F", lw=0.6, alpha=0.5)
        L = lid[lid.gis == g].sort_values("yr")
        ax.plot(L.yr, L.z, "o", color="#1A1A1A", ms=4.5, label="Lidar (median)", zorder=5)
        if text.get("indep"):
            ind = indep[(indep.gis == g) & (indep.survey.astype(str) >= "2004")]
            ax.plot(ind.survey.astype(float) + 0.5, ind.independent_crest_m, "o", mfc="none", mec="#1A1A1A", ms=5,
                    label="Lidar, highest point seaward of the road (median)", zorder=5)
        ax.axvline(EVENT_YEAR, color="#AE017E", lw=0.8, ls=":")
        ax.set_title(f"GIS {g}" + (" (road relocated)" if g in RELOCATED_GIS else ""), loc="left", fontsize=9,
                     fontweight="bold")
        ax.grid(axis="y", color="0.92", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axs[:, 0]:
        ax.set_ylabel("Crest (m NAVD88)", fontsize=8)
    for ax in axs[-1]:
        ax.set_xlabel("Year")
    axs[0, 0].set_xlim(1996, 2010)
    axs[0, 0].set_ylim(0, 9)
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=3, fontsize=8, frameon=False)
    out = NOW_OUT / f"{tag}_crests_GIS8-15.png"
    S.save(fig, out, bbox_inches="tight", close=True)
    S.record_caption(out, (
        "Domain-median dune crest (m NAVD88), 1996-2009, GIS 8-15, no collision module, the road rebuilding dunes only "
        "where it is within 20 m behind them: " + text["runs_text"] + ", against the lidar domain median at each "
        "survey" + (" (filled: the picked dune crest; open, from 2004: the highest point seaward of the road, which "
                    "catches the rebuilt dune)" if text.get("indep") else "") + ". The model's dune line only: a "
        "ridge in the interior is not drawn. Model state at 1 January each year. Magenta dotted: Hurricane Dennis, "
        "after which the road was moved; thin red: years the road rebuilt the dunes in the "
        + NOW_RUNS[text["rebuild_run"]][1].lower() + " run."))
    print(f"wrote {out}")


def main():
    if sys.argv[1:] == ["now"]:
        return now()
    if sys.argv[1:] == ["flat"]:
        return now(FLAT_RUNS, FLAT_TEXT)
    if sys.argv[1:] == ["flatshift"]:
        return now(FLATSHIFT_RUNS, FLATSHIFT_TEXT)
    if sys.argv[1:] == ["shiftdesigns"]:
        return shift_designs()
    if sys.argv[1:] == ["shift"]:
        return now(SHIFT_RUNS, SHIFT_TEXT)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as S
    S.apply_style()
    (OUT / "supporting").mkdir(parents=True, exist_ok=True)
    names = [v for base, _, _ in PAIRS for v in (base, base + "_reloc")]
    runs = {v: _load(v) for v in names}
    road = road_table(runs)
    road.to_csv(OUT / "supporting" / "relocation_road_by_domain.csv", index=False)
    summary = road[road.gis.isin(TS_GIS)].groupby(["variant", "gis"]).agg(
        setback_1997=("setback_m", "first"), setback_2009=("setback_m", "last"),
        rebuilds=("dunes_rebuilt", "sum"), years_rebuild_off=("rebuild_off", "sum")).reset_index()
    pd.set_option("display.width", 200)
    print(summary.to_string(index=False))
    m = pd.DataFrame([R.metrics(v, 0.0) for v in names]).drop(columns="cs")
    m.round(4).to_csv(OUT / "supporting" / "relocation_metrics.csv", index=False)
    print(m.round(3).to_string(index=False))
    crests_over_time(S, plt, runs)
    t = differences(S, plt, runs)
    print(t[t.gis.isin(TS_GIS)][["run", "gis", "d_crest", "d_overwash", "d_lrr"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
