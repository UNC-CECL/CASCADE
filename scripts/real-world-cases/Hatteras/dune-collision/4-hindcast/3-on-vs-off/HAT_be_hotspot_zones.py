"""
Background erosion in the missed crest-drop hotspots, step 1 (2026-11-12, Laura's test): where the base run misses the
lidar crest drops, and how those stretches line up with the source/sink residual (smoothed CoastSat LRR minus the
model's) and its zone rules (|residual| > 0.5 m/yr over at least 3 adjacent domains of one sign; GIS 1, 90 and the
groin's 5-7 excluded; input_prep/7-source-sink/README.md, the zone rework step 2).

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_be_hotspot_zones.py [results]

"zones" draws step 1's figure, hotspot_erosion_zones.png, beside the tables. "results" (step 2) reads the 25-100% runs
(HAT_dune_collision.py BE_HOTSPOT) and writes each design's _alongshore.png and _crests_over_time.png in 2-by-design/,
and in 3-compared/ hotspot_erosion_designs_compared.png and supporting/hotspot_erosion_results.csv (per run and zone:
the imposed rate, the extra retreat it produced and the share retained, the gap left to CoastSat, the extra overwash,
the lidar crest drops caught). "crests" and "positions" build the 2 x 2 composites in 3-compared/. The position change
is HAT_collision_on_off_position.py behot.

The 100% rate of each hotspot zone is that residual, fixed from the base run (measured before any background erosion,
so the overwash it later sets off stays a separate term), in Barrier3D's sign ((-) erosion), as BUXTON_SINK was set.

Writes, in the workspace's dune-collision/4-hindcast/6-missed-drops/1-zones/ (D:):
- supporting/be_hotspot_by_gis.csv: per GIS, missed and caught drops, the residual, the zone verdicts
- supporting/be_hotspot_zones.csv: each hotspot zone, its overlap with a source/sink zone, its 100% rates

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-12
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_collision_on_off_position as P  # noqa: E402
import HAT_collision_missed_drops as M  # noqa: E402

H = M.H

# --- CONFIG ------------------------------------------------------------------
OUT = P.OUT / "6-missed-drops"
ZONES_DIR = OUT / "1-zones"                     # step 1: the zones, their tables and figure
DESIGN_DIR = OUT / "2-by-design"                # step 2: each design's own figures
COMPARED_DIR = OUT / "3-compared"               # the designs side by side, the results table
BASE = "cal_toeberm_lat_cs3_maxmed"            # face retreat, low back row, Cs 3e-6: where the hotspots were found
SIG_M_YR, MIN_WIDTH = 0.5, 3                   # the source/sink zone rules
EXCLUDE = {1, 90, 5, 6, 7}                     # locked ends, groin footprint
MAX_GAP = 1                                    # a hotspot may skip one domain without a missed drop
# -----------------------------------------------------------------------------


def runs_of(mask, gis, max_gap=0):
    """Contiguous runs of True over gis (allowing max_gap False domains inside), as (start, end) pairs."""
    on = [g for g, m in zip(gis, mask) if m]
    out = []
    for g in on:
        if out and g - out[-1][1] <= max_gap + 1:
            out[-1][1] = g
        else:
            out.append([g, g])
    return [tuple(r) for r in out]


def target_lrr():
    from cascade_pipeline.coastsat_lowess import build_coastsat_series
    from cascade_pipeline.hindcast import build_target_table
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    sys.path.insert(0, str(H._HAT_DIR / "figure_making" / "model_output"))
    from rerender_run_figures import COASTSAT_DATASETS, LOWESS_CONFIG
    from run_figures import TARGET_WINDOW
    cs = build_coastsat_series(COASTSAT_DATASETS, active_period_start=H.RUN_START, lowess_config=LOWESS_CONFIG,
                               domains=DOM)
    t = build_target_table(next(c for c in cs if c["active"]), LOWESS_CONFIG, DOM, H.TARGET_WINDOW
                           if hasattr(H, "TARGET_WINDOW") else TARGET_WINDOW)
    return t.set_index("gis_domain")["target_lrr_m_yr"]


def main():
    (ZONES_DIR / "supporting").mkdir(parents=True, exist_ok=True)
    d = M.pairs(BASE)
    d["cls"] = d.apply(M.classify, axis=1)
    per = d.groupby("gis").agg(missed=("cls", lambda c: int((c == 0).sum())),
                               caught=("cls", lambda c: int((c == 1).sum())),
                               missed_intervals=("cls", lambda c: ""))
    per["missed_intervals"] = d[d.cls == 0].groupby("gis").apply(
        lambda x: " ".join(f"{a}-{b % 100:02d}" for a, b in zip(x.a, x.b)), include_groups=False)
    per = per.reindex(range(1, 91))
    rate = pd.read_csv(H.run_dir(BASE) / "tables" / "shoreline_change_rate.csv").set_index("gis_domain").lrr_m_yr
    tgt = target_lrr()
    per["model_lrr"], per["coastsat_lrr"] = rate, tgt
    per["residual"] = tgt - rate                              # (-): the model erodes too little; Barrier3D's sign
    eligible = ~per.index.isin(EXCLUDE)
    per["sig_sink"] = eligible & (per.residual <= -SIG_M_YR)
    per["sig_source"] = eligible & (per.residual >= SIG_M_YR)
    gis = list(per.index)
    ss = [(a, b, "sink") for a, b in runs_of(per.sig_sink.to_numpy(), gis) if b - a + 1 >= MIN_WIDTH]
    ss += [(a, b, "source") for a, b in runs_of(per.sig_source.to_numpy(), gis) if b - a + 1 >= MIN_WIDTH]
    per["ss_zone"] = ""
    for a, b, k in ss:
        per.loc[a:b, "ss_zone"] = f"{k} {a}-{b}"
    hot = [(a, b) for a, b in runs_of((per.missed.fillna(0) > 0).to_numpy(), gis, MAX_GAP) if b - a + 1 >= MIN_WIDTH]
    rows = []
    for a, b in hot:
        z = per.loc[a:b]
        rows.append(dict(zone=f"GIS {a}-{b}", start=a, end=b, missed=int(z.missed.sum()), caught=int(z.caught.sum()),
                         mean_residual=round(float(z.residual.mean()), 2),
                         sink_domains=int(z.sig_sink.sum()), source_domains=int(z.sig_source.sum()),
                         ss_zones=", ".join(sorted({s for s in z.ss_zone if s})) or "none",
                         rates_100=" ".join(f"{g}={r:+.2f}" for g, r in z.residual.items() if r < 0)))
    zones = pd.DataFrame(rows)
    per.round(3).to_csv(ZONES_DIR / "supporting" / "be_hotspot_by_gis.csv")
    zones.to_csv(ZONES_DIR / "supporting" / "be_hotspot_zones.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 120)
    print("source/sink zones:", ss)
    print(zones.to_string(index=False))
    print(per[per.missed.fillna(0) > 0][["missed", "caught", "missed_intervals", "model_lrr", "coastsat_lrr",
                                          "residual", "ss_zone"]].round(2).to_string())


# --- step 2: the 25-100% runs on each design (HAT_dune_collision.py BE_HOTSPOT) ---------------------------------------
FRACS = (25, 50, 75, 100)
# Each design: run name -> (plain name, file tag, colour in the comparison figure)
DESIGNS = {
    "cal_step2_growth_noscrape": ("No collision, equal rows", "nocollision_equal", "#7F7F7F"),
    "cal_toeberm_cs2": ("Lower both rows, equal rows, Cs 2e-6", "lowerboth_equal_cs2e-6", "#7B3294"),
    "cal_toeberm_lat_cs10": ("Face retreat, equal rows, Cs 1e-5", "faceretreat_equal_cs1e-5", "#2171B5"),
    "cal_toeberm_lat_cs3_maxmed": ("Face retreat, low back row, Cs 3e-6", "faceretreat_lowback_cs3e-6", "#E66101"),
}
COLOURS = {0: "#1A1A1A", 25: "#CBC9E2", 50: "#9E9AC8", 75: "#756BB1", 100: "#54278F"}
TS_GIS = (11, 12, 13, 23, 84, 86)
EVENTS = ((1998.65, "Bonnie"), (1999.67, "Dennis"), (2003.71, "Isabel"))
ZONES = {"Buxton, GIS 10-14": range(10, 15), "Avon south, GIS 21-24": range(21, 25),
         "Rodanthe, GIS 82-89": range(82, 90)}
ZONE_COLOURS = {"Buxton, GIS 10-14": "#E66101", "Avon south, GIS 21-24": "#1A9850", "Rodanthe, GIS 82-89": "#2171B5"}
SOURCE_HOTSPOTS = ((35, 38), (51, 55))          # missed-drop hotspots on a source/sink source: no erosion added
MARKERS = ("o", "s", "^", "D")
SHARED = "toe at the berm, growth 0.02-0.04, groin 0.16, the road fix, no road rule, Rodanthe 40%, v3's interior"


def runs_for(base):
    return [(base, 0)] + [(f"{base}_be{p}", p) for p in FRACS]


def per_gis(v):
    from cascade_pipeline.saved_runs import load_saved_cascade
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    rd = H.run_dir(v)
    rate = pd.read_csv(rd / "tables" / "shoreline_change_rate.csv").set_index("gis_domain").lrr_m_yr
    ow = np.load(next(rd.glob("*_overwash.npy")))
    c = load_saved_cascade(rd / f"{rd.name}.npz")
    dom = DOM.for_total(len(c.barrier3d))
    crest = {g: (np.asarray(c.barrier3d[dom.gis_to_pad(g)].DuneDomain, float).max(axis=2)
                 + c.barrier3d[dom.gis_to_pad(g)]._BermEl) * 10 + H.MHW_NAVD88_M for g in range(1, 91)}
    t = pd.DataFrame({"lrr": rate, "overwash": [ow[1:, dom.gis_to_pad(g)].sum() for g in range(1, 91)]},
                     index=range(1, 91))
    return t, crest


def drops_of(v):
    d = M.pairs(v)
    d["cls"] = d.apply(M.classify, axis=1)
    return d


def design_results(base, per, R):
    """Per run and zone: imposed rate, extra retreat, share retained, gap left, extra overwash, drops."""
    runs = runs_for(base)
    data = {v: per_gis(v) for v, _ in runs}
    drops = {v: drops_of(v) for v, _ in runs}
    b, _ = data[base]
    rate100 = H.BE_HOTSPOT_100[base]
    all_z = [g for r in ZONES.values() for g in r]
    outside = [g for g in range(2, 90) if g not in all_z]
    rows = []
    for v, p in runs:
        t, _ = data[v]
        m = R.metrics(v, 0.0)
        for zname, zr in [*ZONES.items(), ("All three", all_z), ("Outside them, GIS 2-89", outside)]:
            zr = list(zr)
            imposed = -p / 100 * float(np.mean([rate100.get(g, 0.0) for g in zr])) if zr is not outside else 0.0
            dl = float((b.lrr[zr] - t.lrr[zr]).mean())                           # extra retreat, m/yr (+ landward)
            dd = drops[v][drops[v].gis.isin(zr)]
            rows.append(dict(design=DESIGNS[base][0], run=v, pct=p, zone=zname, imposed_m_yr=imposed,
                             extra_retreat_m_yr=dl, retained=dl / imposed if imposed else np.nan,
                             residual_left_m_yr=float((per.loc[zr, "coastsat_lrr"] - t.lrr[zr]).mean()),
                             overwash_m3_m=float(t.overwash[zr].mean()),
                             extra_overwash_m3_m=float(t.overwash[zr].mean() - b.overwash[zr].mean()),
                             drops_caught=int((dd.cls == 1).sum()), drops_missed=int((dd.cls == 0).sum()),
                             model_only_drops=int((dd.cls == 2).sum()),
                             island_crest_rmse_m=m["crest_rmse_m"], island_false_alarm_pct=m["overwash_false_alarm_pct"],
                             island_shoreline_rmse_m_yr=m["shoreline_rmse_m_yr"]))
    return pd.DataFrame(rows), data, drops


def alongshore_figure(plt, S, base, per, data, res):
    """(a) shoreline rate against CoastSat, (b) overwash added, alongshore; (c-e) per zone against the fraction."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    name, tag, _ = DESIGNS[base]
    g = np.arange(1, 91)
    fig = plt.figure(figsize=(S.FIG_W_DOUBLE, 7.6), layout="constrained")
    gs = fig.add_gridspec(3, 3, height_ratios=(1.25, 0.8, 1.0))
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, :], sharex=ax_a)
    small = [fig.add_subplot(gs[2, k]) for k in range(3)]
    b0, _ = data[base]
    for ax in (ax_a, ax_b):
        for zn, zr in ZONES.items():
            ax.axvspan(zr[0] - 0.5, zr[-1] + 0.5, color=ZONE_COLOURS[zn], alpha=0.10, lw=0, zorder=0)
        for a, b in SOURCE_HOTSPOTS:
            ax.axvspan(a - 0.5, b + 0.5, facecolor="none", edgecolor="0.75", hatch="////", lw=0, zorder=0)
        ax.axhline(0, color="0.55", lw=0.6, zorder=1)
        ax.grid(axis="y", color="0.92", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    for v, p in runs_for(base):
        t, _ = data[v]
        ax_a.plot(g, t.lrr.reindex(g), color=COLOURS[p], lw=1.9 if p == 0 else 1.3, zorder=3 if p == 0 else 2)
        if p:
            ax_b.plot(g, (t.overwash - b0.overwash).reindex(g), color=COLOURS[p], lw=1.3, zorder=2)
    ax_a.plot(g, per.coastsat_lrr.reindex(g), color="#E31A1C", lw=1.6, ls=(0, (5, 2)), zorder=4)
    for zn, zr in ZONES.items():
        ax_a.text((zr[0] + zr[-1]) / 2, 1.01, zn.split(",")[0], transform=ax_a.get_xaxis_transform(),
                  ha="center", va="bottom", fontsize=7.5, color=ZONE_COLOURS[zn], fontweight="bold")
    ax_a.set_ylabel("Shoreline change rate\n1996–2009 (m/yr)")
    ax_b.set_ylabel("Overwash added\n(m³/m)")
    ax_b.set_xlabel(S.DOMAIN_AXIS_LABEL)
    ax_a.set_xlim(0.5, 90.5)
    ax_a.tick_params(labelbottom=False)
    run_h = ([Line2D([], [], color="#E31A1C", lw=1.6, ls=(0, (5, 2)), label="CoastSat (LOWESS)"),
              Line2D([], [], color=COLOURS[0], lw=1.9, label=name)]
             + [Line2D([], [], color=COLOURS[p], lw=1.3, label=f"+{p}% hotspot erosion" if p == 25 else f"+{p}%")
                for p in FRACS]
             + [Patch(facecolor="none", edgecolor="0.6", hatch="////", label="Hotspot on a source (none added)")])
    fig.legend(handles=run_h, loc="outside upper center", ncol=4, fontsize=7.5, frameon=False, handlelength=2.2,
               columnspacing=1.4)
    x = np.array((0, *FRACS))
    for zn in ZONES:
        z = res[res.zone == zn].set_index("pct").reindex(x)
        c = ZONE_COLOURS[zn]
        small[0].plot(x, z.residual_left_m_yr, "-o", color=c, ms=4, lw=1.4)
        small[1].plot(x[1:], z.retained.iloc[1:], "-o", color=c, ms=4, lw=1.4)
        small[2].plot(x, z.drops_caught, "-o", color=c, ms=4, lw=1.4)
        small[2].plot(x, z.model_only_drops, "--", marker="o", mfc="white", color=c, ms=4, lw=1.0)
    small[0].axhline(0, color="0.55", lw=0.6)
    small[1].axhline(1.0, color="0.55", lw=0.6)
    small[0].set_ylabel("Gap to CoastSat (m/yr)\n(−) model erodes too little")
    small[1].set_ylabel("Extra retreat ÷ imposed rate")
    small[2].set_ylabel("Crest drops (count)")
    small[2].yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    for ax in small:
        ax.set_xticks(x)
        ax.set_xlim(-8, 108)
        ax.set_xlabel("Hotspot erosion added (%)")
        ax.grid(axis="y", color="0.92", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    small[1].set_xticks(x[1:])
    small[1].set_xlim(17, 108)
    zone_h = [Line2D([], [], color=c, marker="o", ms=4, lw=1.4, label=n) for n, c in ZONE_COLOURS.items()]
    small[0].legend(handles=zone_h, loc="lower right", fontsize=6.5, frameon=False)
    small[2].legend(handles=[Line2D([], [], color="0.3", marker="o", ms=4, lw=1.4, label="Lidar drop caught"),
                             Line2D([], [], color="0.3", marker="o", mfc="white", ms=4, lw=1.0, ls="--",
                                    label="Model drop, none in lidar")], loc="upper left", fontsize=6.5, frameon=False)
    for k, ax in enumerate((ax_a, ax_b, *small)):
        ax.text(-0.01 if k < 2 else -0.04, 1.02, f"({'abcde'[k]})", transform=ax.transAxes, fontsize=10,
                fontweight="bold", ha="right", va="bottom")
    out = DESIGN_DIR / f"hotspot_erosion_{tag}_alongshore.png"
    S.save(fig, out, bbox_inches="tight", close=True)
    S.record_caption(out, (
        "Background erosion added in the hotspots of missed lidar crest drops that lie on a source/sink sink "
        "(shaded: Buxton GIS 10-14, Avon south GIS 21-24, Rodanthe GIS 82-89), 1996-2009, on the design "
        f"{name} (black; {SHARED}). Purples, light to dark: the same with 25, 50, 75 and 100% of the rate that "
        "closes this design's gap to CoastSat in those domains (smoothed CoastSat LRR minus its own LRR, fixed before "
        "any erosion is added). Hatched: hotspots on a source/sink source (GIS 35-38, 51-55), where no erosion was "
        "added. (a) Shoreline change rate (linear regression) against CoastSat (red dashed, 7-domain LOWESS); (b) "
        "overwash each domain received over 1996-2009 beyond the base design's; (c) zone-mean gap left to CoastSat "
        "(CoastSat minus model; negative, the model erodes too little); (d) the zone's extra retreat over the base "
        "divided by the rate imposed (1, the shoreline moves exactly as imposed); (e) lidar crest drops of at least "
        "0.5 m the run catches (the model drops at least half as much; solid) and model drops with none in the "
        "lidar (dashed, open)."))
    print(f"wrote {out}")


def crests_figure(plt, S, base, data, lid):
    name, tag, _ = DESIGNS[base]
    fig, axs = plt.subplots(3, 2, figsize=(S.FIG_W_DOUBLE, 7.2), sharex=True, sharey=True, layout="constrained")
    for ax, gg in zip(axs.T.ravel(), TS_GIS):
        for v, p in runs_for(base):
            _, crest = data[v]
            yrs = H.RUN_START + 1 + np.arange(crest[gg].shape[0])
            ax.plot(yrs, np.median(crest[gg], axis=1), color=COLOURS[p], lw=1.6,
                    label=name if p == 0 else f"+{p}% hotspot erosion")
        L = lid[lid.gis == gg].sort_values("yr")
        ax.plot(L.yr, L.z, "o", color="#E31A1C", ms=4.5, label="Lidar (domain median)", zorder=5)
        ax.axhline(H.TOE_AT_BERM_M, color="0.55", lw=0.8, label="Berm, 1.7 m (the toe)")
        for x, _ in EVENTS:
            ax.axvline(x, color="0.7", lw=0.6, ls=":")
        ax.set_title(f"GIS {gg}", loc="left", fontsize=9, fontweight="bold")
        ax.grid(axis="y", color="0.92", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axs[:, 0]:
        ax.set_ylabel("Crest (m NAVD88)")
    for ax in axs[-1]:
        ax.set_xlabel("Year")
    for x, n in EVENTS:
        for ax in axs[0]:
            ax.text(x, 0.97, f" {n} " if n != "Bonnie" else f"{n} ", transform=ax.get_xaxis_transform(), fontsize=7,
                    color="0.4", va="top", ha={"Bonnie": "right", "Dennis": "left"}.get(n, "left"))
    axs[0, 0].set_xlim(1996, 2010)
    axs[0, 0].set_ylim(0, 9)
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=4, fontsize=7.5, frameon=False)
    out = DESIGN_DIR / f"hotspot_erosion_{tag}_crests_over_time.png"
    S.save(fig, out, bbox_inches="tight", close=True)
    S.record_caption(out, (
        "Domain-median dune crest (m NAVD88; the higher of Barrier3D's two dune rows) through 1996-2009 in six "
        "hotspot domains (Buxton GIS 11-13, Avon south GIS 23, Rodanthe GIS 84 and 86), against the lidar domain "
        f"median at each survey (red). Black: the design {name} ({SHARED}); purples, light to dark: the same with "
        "25, 50, 75 and 100% hotspot erosion. Model state at 1 January each year; dotted verticals: Bonnie, Dennis, "
        "Isabel."))
    print(f"wrote {out}")


def compare_figure(plt, S, res):
    """Rows: gap left, extra retreat / imposed, overwash added, drops caught, model-only drops; columns: zones."""
    from matplotlib.lines import Line2D
    x = np.array((0, *FRACS))
    rows = (("residual_left_m_yr", "Gap to CoastSat (m/yr)\n(−) erodes too little", 0.0),
            ("retained", "Extra retreat ÷\nimposed rate", 1.0),
            ("extra_overwash_m3_m", "Overwash added\n(m³/m)", 0.0),
            ("drops_caught", "Lidar drops\ncaught", None),
            ("model_only_drops", "Model drops,\nnone in lidar", None))
    fig, axs = plt.subplots(len(rows), 3, figsize=(S.FIG_W_DOUBLE, 8.6), sharex=True, sharey="row",
                            layout="constrained")
    for j, zn in enumerate(ZONES):
        for k, (name, _, c) in enumerate(DESIGNS.values()):
            z = res[(res.design == name) & (res.zone == zn)].set_index("pct").reindex(x)
            if z.run.isna().all():
                continue
            dx = (k - 1.5) * 2.2                         # a small sideways offset, so equal values stay visible
            for i, (col, _, _) in enumerate(rows):
                xx, yy = (x[1:], z[col].iloc[1:]) if col == "retained" else (x, z[col])
                axs[i, j].plot(xx + dx, yy, "-", marker=MARKERS[k], color=c, ms=4, lw=1.3, mfc="white", mew=1.2)
        axs[0, j].set_title(zn, fontsize=9, fontweight="bold", color=ZONE_COLOURS[zn])
        for i, (col, lab, ref) in enumerate(rows):
            ax = axs[i, j]
            if ref is not None:
                ax.axhline(ref, color="0.55", lw=0.6)
            ax.grid(axis="y", color="0.92", lw=0.5)
            ax.spines[["top", "right"]].set_visible(False)
            if j == 0:
                ax.set_ylabel(lab, fontsize=8)
            if col in ("drops_caught", "model_only_drops"):
                ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
        axs[-1, j].set_xticks(x)
        axs[-1, j].set_xlabel("Hotspot erosion added (%)")
    for k, ax in enumerate(axs.ravel()):
        ax.text(-0.04, 1.02, f"({chr(97 + k)})", transform=ax.transAxes, fontsize=8.5, fontweight="bold",
                ha="right", va="bottom")
    h = [Line2D([], [], color=c, marker=MARKERS[k], ms=4, lw=1.3, mfc="white", mew=1.2, label=n)
         for k, (n, _, c) in enumerate(DESIGNS.values())]
    fig.legend(handles=h, loc="outside upper center", ncol=2, fontsize=7.5, frameon=False)
    out = COMPARED_DIR / "hotspot_erosion_designs_compared.png"
    S.save(fig, out, bbox_inches="tight", close=True)
    S.record_caption(out, (
        "Background erosion added in the three hotspot zones of missed lidar crest drops on a source/sink sink "
        "(columns: Buxton GIS 10-14, Avon south GIS 21-24, Rodanthe GIS 82-89), 1996-2009, on four dune collision "
        f"designs (lines; {SHARED}): 0 to 100% of the rate that closes each design's own gap to CoastSat in those "
        "domains. Rows: the zone-mean gap left to CoastSat (CoastSat minus model; negative, the model erodes too "
        "little); the zone's extra retreat over its base divided by the rate imposed; the overwash added over "
        "1996-2009 beyond the base; lidar crest drops of at least 0.5 m caught (the model drops at least half as "
        "much); model drops of at least 0.5 m with none in the lidar. Designs: no collision (module off); lower both "
        "rows (the eroded volume taken from both dune rows in proportion to their height); face retreat (the "
        "seaward row cut back until its sand is gone, then the back row); equal rows (both rows start at the "
        "highest cell of each profile's dune window) or low back row (the back row at the window's median). "
        "The designs are offset slightly along the x axis so that equal values stay visible."))
    print(f"wrote {out}")


def results():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import HAT_collision_row_mode_sweep as R
    from site_layer import hat_figure_style as S
    S.apply_style()
    per = pd.read_csv(ZONES_DIR / "supporting" / "be_hotspot_by_gis.csv", index_col=0)
    sv = pd.read_csv(H.SURVEY_CRESTS)
    sv = sv[(sv.status == "ok") & sv.gis.isin(TS_GIS)].copy()
    tt = pd.to_datetime(sv.flown.astype(str).str[:7], errors="coerce")
    sv["yr"] = tt.dt.year + (tt.dt.month - 0.5) / 12
    lid = sv.groupby(["gis", "survey"]).agg(z=("crest_m_navd88", "median"), yr=("yr", "first")).reset_index()
    every = []
    for base in DESIGNS:
        if any(H.run_dir(v) is None for v, _ in runs_for(base)):
            print(f"  skipping {DESIGNS[base][0]}: runs missing")
            continue
        res, data, _ = design_results(base, per, R)
        every.append(res)
        alongshore_figure(plt, S, base, per, data, res)
        crests_figure(plt, S, base, data, lid)
    res = pd.concat(every)
    res.round(3).to_csv(COMPARED_DIR / "supporting" / "hotspot_erosion_results.csv", index=False)
    pd.set_option("display.width", 250)
    print(res[res.zone.isin(["All three"])].drop(columns=["run", "zone"]).round(2).to_string(index=False))
    compare_figure(plt, S, res)


def crests():
    """The four designs' crests-over-time figures in a 2 x 2 grid."""
    composite("crests_over_time", "hotspot_erosion_crests_designs_compared.png", (
        "Domain-median dune crest through 1996-2009 in six hotspot domains (Buxton GIS 11-13, Avon south GIS 23, "
        "Rodanthe GIS 84 and 86) against the lidar, with 0-100% hotspot erosion, each panel one design's own figure: "
        "no collision, equal rows; lower both rows, equal rows, Cs 2e-6; face retreat, equal rows, Cs 1e-5; face "
        f"retreat, low back row, Cs 3e-6 ({SHARED})."))


def positions():
    """The four designs' position-change figures (HAT_collision_on_off_position.py behot) in a 2 x 2 grid."""
    composite("position_1996_2009", "hotspot_erosion_position_designs_compared.png", (
        "Shoreline position change, 1996-2009, against CoastSat, with 0-100% hotspot erosion (GIS 10-14, 21-24, "
        "82-89), each panel one design's own figure: no collision, equal rows; lower both rows, equal rows, Cs 2e-6; "
        f"face retreat, equal rows, Cs 1e-5; face retreat, low back row, Cs 3e-6 ({SHARED})."))


def composite(kind, name_out, caption):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import record_caption
    imgs = [(mpimg.imread(DESIGN_DIR / f"hotspot_erosion_{tag}_{kind}.png"), name)
            for name, tag, _ in DESIGNS.values()]
    h, w = imgs[0][0].shape[:2]
    scale = 9.0 / w
    fig, axs = plt.subplots(2, 2, figsize=(2 * w * scale, 2 * (h * scale + 0.45)))
    for ax, (img, name) in zip(axs.ravel(), imgs):
        ax.imshow(img)
        ax.set_title(name, loc="left", fontsize=13, fontweight="bold")
        ax.axis("off")
    fig.tight_layout()
    out = COMPARED_DIR / name_out
    fig.savefig(out, dpi=200, facecolor="white")
    plt.close(fig)
    record_caption(out, caption)
    print(f"wrote {out}")


def zones():
    """Step 1 on its own: how the zones were found. (a) every lidar survey interval and domain classed by the base
    run (missed, caught, model only, neither); (b) missed drops per domain, the missed-drop stretches bracketed;
    (c) the source/sink residual with its +-0.5 m/yr rule and the source and sink zones it makes. The three zones
    given background erosion are banded through all three; the two on a source are hatched."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    from site_layer import hat_figure_style as S
    S.apply_style()
    per = pd.read_csv(ZONES_DIR / "supporting" / "be_hotspot_by_gis.csv", index_col=0)
    zt = pd.read_csv(ZONES_DIR / "supporting" / "be_hotspot_zones.csv")
    d = M.pairs(BASE)
    d["cls"] = d.apply(M.classify, axis=1)
    ivs = list(zip(M.SURVEYS[:-1], M.SURVEYS[1:]))
    grid = np.full((len(ivs), 90), np.nan)
    for r in d.itertuples():
        if not np.isnan(r.cls):
            grid[ivs.index((r.a, r.b)), r.gis - 1] = r.cls
    cols = [c for _, c in M.CLASSES] + ["0.80"]                      # 0 missed, 1 caught, 2 model only, 3 none, 4 n/a
    fig, (ax_a, ax_b, ax_c) = plt.subplots(3, 1, figsize=(S.FIG_W_DOUBLE, 7.2), sharex=True, layout="constrained",
                                           gridspec_kw=dict(height_ratios=(1.25, 0.65, 1.0)))
    ax_a.imshow(grid, aspect="auto", interpolation="nearest", cmap=ListedColormap(cols), vmin=-0.5, vmax=4.5,
                extent=(0.5, 90.5, len(ivs) - 0.5, -0.5))
    for k in range(1, len(ivs)):
        ax_a.axhline(k - 0.5, color="white", lw=0.8)
    ax_a.set_yticks(range(len(ivs)), [f"{a}–{b}" + (" (no step)" if (a, b) in M.NO_STEP else "") for a, b in ivs],
                    fontsize=7)
    ax_a.set_ylabel("Lidar survey interval")
    ax_a.legend(handles=[Patch(facecolor=c, edgecolor="0.6", lw=0.3, label=n) for n, c in M.CLASSES]
                + [Patch(facecolor="0.80", label="Same model year, not scored")],
                loc="lower left", bbox_to_anchor=(0, 1.02), ncol=5, frameon=False, fontsize=7, handlelength=1.2,
                columnspacing=1.0)
    g = per.index.to_numpy()
    ax_b.bar(g, per.missed.fillna(0), width=0.8, color=M.CLASSES[0][1], lw=0)
    ymax = per.missed.max() + 1.5
    for r in zt.itertuples():
        ax_b.plot([r.start - 0.4, r.start - 0.4, r.end + 0.4, r.end + 0.4], [ymax - 0.9, ymax - 0.4, ymax - 0.4,
                  ymax - 0.9], color=S.INK, lw=0.8)
        ax_b.text((r.start + r.end) / 2, ymax - 0.25, f"{r.start}–{r.end}", ha="center", va="bottom", fontsize=7,
                  color=S.INK)
    ax_b.set_ylim(0, ymax + 1.9)
    ax_b.set_ylabel("Missed drops\n(count)")
    res = per.residual
    sig = per.sig_sink | per.sig_source
    ax_c.bar(g, res, width=0.8, lw=0, color=np.where(per.index.isin(EXCLUDE), "0.88", np.where(sig, "0.30", "0.70")))
    for y in (-SIG_M_YR, SIG_M_YR):
        ax_c.axhline(y, color=S.INK_MUTED, lw=0.7, ls="--")
    ax_c.axhline(0, color=S.INK, lw=0.6)
    lim = np.nanmax(np.abs(res)) + 0.4
    for zname in sorted({z for z in per.ss_zone.dropna() if z}):
        kind, span = zname.split()
        a, b = map(int, span.split("-"))
        y = -lim + 0.15 if kind == "sink" else lim - 0.15
        ax_c.plot([a - 0.4, b + 0.4], [y, y], color=S.INK, lw=1.6, solid_capstyle="butt")
        ax_c.text((a + b) / 2, y + (0.12 if kind == "sink" else -0.12), f"{kind} {a}–{b}", ha="center",
                  va="bottom" if kind == "sink" else "top", fontsize=6.5, color=S.INK, zorder=6,
                  path_effects=S._halo(2.0))
    ax_c.set_ylim(-lim, lim)
    ax_c.set_ylabel("CoastSat − model shoreline\nrate (m/yr; − sink, + source)")
    ax_c.set_xlim(0.5, 90.5)
    ax_c.set_xticks([1] + list(range(10, 91, 10)))
    ax_c.set_xlabel(S.DOMAIN_AXIS_LABEL)
    for ax in (ax_a, ax_b, ax_c):
        for zn, zr in ZONES.items():
            if ax is ax_a:
                ax.axvspan(zr[0] - 0.5, zr[-1] + 0.5, facecolor="none", edgecolor=ZONE_COLOURS[zn], lw=1.6, zorder=5)
            else:
                ax.axvspan(zr[0] - 0.5, zr[-1] + 0.5, color=ZONE_COLOURS[zn], alpha=0.15, lw=0, zorder=0)
        for a, b in SOURCE_HOTSPOTS:
            if ax is not ax_a:
                ax.axvspan(a - 0.5, b + 0.5, facecolor="none", edgecolor="0.75", hatch="////", lw=0, zorder=0)
        if ax is not ax_a:
            ax.spines[["top", "right"]].set_visible(False)
    for zn, zr in ZONES.items():
        ax_b.text((zr[0] + zr[-1]) / 2, ymax + 1.85, zn.split(",")[0], ha="center", va="top", fontsize=7.5,
                  fontweight="bold", color=ZONE_COLOURS[zn])
    for ax, L in zip((ax_a, ax_b, ax_c), "abc"):
        ax.set_title(f"({L})", loc="left", fontweight="bold", x=-0.085 if ax is ax_a else -0.085)
    out = ZONES_DIR / "hotspot_erosion_zones.png"
    S.save(fig, out, dpi=300, facecolor="white", close=True)
    S.record_caption(out, (
        "How the hotspot erosion zones were found (step 1), from the base run (face retreat, low back row, Cs 3e-6; "
        f"{SHARED}). (a) Each lidar survey interval and domain: red, the lidar crest fell at least "
        f"{M.DROP_M:g} m and the model fell less than half as much (missed); green, caught; blue, the model fell at "
        "least that much with no drop in the lidar; light, neither; white, too few surveyed cells. (b) Missed drops "
        "per domain; bracketed, the missed-drop stretches (three or more domains, a one-domain gap allowed). (c) The "
        "source/sink residual, smoothed CoastSat 1996-2009 shoreline rate minus the model's (negative: the model "
        f"erodes too little); dashed, the ±{SIG_M_YR:g} m/yr rule; dark bars pass it; pale, GIS 1, 5-7 and 90, "
        f"excluded. A sink or source zone is {MIN_WIDTH} or more neighbouring domains past the rule with one sign. "
        "Coloured, the three missed-drop stretches on a sink, which got the background erosion (Buxton GIS 10-14, "
        "Avon south GIS 21-24, Rodanthe GIS 82-89; the sink's extent); hatched, the two on a source (GIS 35-38, "
        "51-55), given none."))
    print(f"wrote {out}")


if __name__ == "__main__":
    {"results": results, "positions": positions, "crests": crests, "zones": zones}.get(sys.argv[1] if len(sys.argv) > 1 else "", main)()
