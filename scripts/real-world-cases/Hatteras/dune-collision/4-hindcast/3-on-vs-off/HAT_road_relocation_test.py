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

Writes, in the workspace's dune-collision/4-hindcast/3-on-vs-off/6-other-tests/road-relocation-1999/ (D:):
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
OUT = P.OUT / "6-other-tests" / "road-relocation-1999"
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


def main():
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
