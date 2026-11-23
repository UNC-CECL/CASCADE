"""
The Cs sensitivity of every collision design (2026-11-13): each module design run at Cs 1, 2, 3, 5, 10, 20, 30 x 1e-6
on the village-migration default (HAT_dune_collision.py DESIGN_SWEEP), scored against the lidar crests, the
overwash imagery and CoastSat, so each design's Cs can be picked by hand from its own curve.

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_design_cs_sweep.py

Writes, in the workspace's dune-collision/4-hindcast/3-collision-design/ (D:; collision_designs.py names the folders):
- <erosion>/<start>/<tag without Cs>_cs_sensitivity_1996_2009.png (+ supporting/ CSV): that design's scores against
  Cs, no collision on the same starting rows dashed grey, the Cs its position plot now uses ringed ("current")
- 0-all-designs/all_designs_cs_sensitivity_1996_2009.png (+ supporting/ CSV): the six module designs together
- <erosion>/<start>/cs-sweep/<tag without Cs>_cs<Cs>_crest_offset_1996_2009.png: each sweep run's own crest figure
  (figures/dune_crests/dune_crest_offset_overview_v2.png: model minus lidar crest by domain and survey, its shoreline
  change against CoastSat on top), copied and renamed; for no collision, <tag>_crest_offset_1996_2009.png in its folder
- <erosion>/<start>/<erosion>_<start>_cs_sweep_position_1996_2009.png: shoreline position change against CoastSat,
  one line per Cs (light to dark in the design's colour, the current Cs thicker), no collision on the same starting
  rows grey underneath ("positions" draws only these)
Runs not made yet are skipped.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-13
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import collision_designs as CD  # noqa: E402
import HAT_collision_row_mode_sweep as R  # noqa: E402  (metrics(), H)
import HAT_be_hotspot_zones as Z  # noqa: E402

H = R.H

# --- CONFIG ------------------------------------------------------------------
ZONES = Z.ZONES_AVN2                       # Buxton 9-14, Avon south 23-26, Rodanthe 82-89
PANELS = (("crest_rmse_m", "Crest RMSE, GIS 6–90 (m)"),
          ("crest_rmse_gis9_15_m", "Crest RMSE, GIS 9–15 (m)"),
          ("crest_rmse_zones_m", "Crest RMSE, hotspot zones (m)"),
          ("crest_offset_zones_m", "Mean crest offset, hotspot zones (m)"),
          ("bias_2009_m", "Mean crest offset, 2009 (m)"),
          ("within_0p5_pct", "Domain-surveys within ±0.5 m (%)"),
          ("overwash_false_alarm_pct", "Overwash false alarms (%)"),
          ("shoreline_rmse_m_yr", "Shoreline RMSE (m/yr)"))
ZERO_LINE = ("bias_2009_m", "crest_offset_zones_m")
# -----------------------------------------------------------------------------


def key(d):
    return f"{d.erosion}_{d.start}"


def score(v, cs):
    """R.metrics' scores, with every crest score recomputed without Buxton's rebuilt-dune cells (GIS 9-14 in 2004,
    2005, 2009: beside the dune rebuilt along the relocated NC-12 the picked lidar crest reads ~0.8 m low;
    HAT_collision_on_off_position.REBUILT_HATCH), plus the hotspot-zone scores."""
    from HAT_collision_on_off_position import REBUILT_HATCH as RB
    m = R.metrics(v, cs)
    from cascade_pipeline.plotting.dune_crest_comparison import TABLE, TABLES
    t = pd.read_csv(H.run_dir(v) / TABLES / TABLE)
    rebuilt = t.gis.between(*RB["gis"]) & t.survey.isin(RB["surveys"])
    t = t[t.gis.between(R.GIS_LO, R.GIS_HI) & (t.survey > 1996) & ~rebuilt]
    e = t.model_minus_survey_m
    f = lambda x: float(np.sqrt((x ** 2).mean()))
    z = t[t.gis.isin([g for r in ZONES.values() for g in r])].model_minus_survey_m
    m.update(crest_rmse_m=f(e), crest_rmse_gis9_15_m=f(t[t.gis.between(9, 15)].model_minus_survey_m),
             bias_2009_m=float(t[t.survey == 2009].model_minus_survey_m.mean()),
             within_0p5_pct=float(100 * (e.abs() <= 0.5).mean()),
             crest_rmse_zones_m=f(z), crest_offset_zones_m=float(z.mean()))
    return m


def table(d):
    """The design's scores at every Cs made so far."""
    rows = []
    for c in H.sweep_cs(key(d)):
        v = H.DESIGN_SWEEP[(key(d), c)]
        if R._done(v):
            rows.append(dict(design=d.name, tag=key(d), current=(v == d.run), **score(v, c * 1e-6)))
    return pd.DataFrame(rows)


def reference(d):
    return next(r for r in CD.DESIGNS if r.erosion == "nocollision" and r.start == d.start)


def draw(S, plt, designs, tables, out, caption, current_apart=False):
    fig, axs = plt.subplots(4, 2, figsize=(S.FIG_W_DOUBLE, 9.6), sharex=True, layout="constrained")
    refs = {}
    for d in designs:
        r = reference(d)
        if r.run not in refs and R._done(r.run):
            refs[r.run] = (r, score(r.run, 0.0))
    for k, (ax, (col, lab)) in enumerate(zip(axs.ravel(), PANELS)):
        for r, m in refs.values():
            ax.axhline(m[col], color="0.5", lw=1.1, ls=(0, (4, 2)) if r.start == "equal" else (0, (1, 1.5)),
                       label=f"{r.name} (no Cs)" if k == 0 else None)
        if col in ZERO_LINE:
            ax.axhline(0, color="0.2", lw=0.7)
        for d in designs:
            t = tables[d.tag].sort_values("cs")
            if t.empty:
                continue
            line = t[~t.current] if current_apart else t      # current_apart: the ringed run is not on the curve
            ax.plot(line.cs, line[col], marker="o", ms=3.5, lw=1.8, color=d.colour, ls=d.linestyle,
                    label=d.name.rsplit(", Cs", 1)[0] if k == 0 else None)
            cur = t[t.current]
            ax.plot(cur.cs, cur[col], "o", ms=9, mfc="none", mec=d.colour, mew=1.4, zorder=6,
                    label=("current Cs (ringed)" if k == 0 and d is designs[0] else None))
        ax.set_xscale("log")
        ticks = sorted({round(c * 1e6) for t in tables.values() for c in t.cs})
        ax.set_xticks([c * 1e-6 for c in ticks])
        ax.set_xticklabels([cs_label(c) for c in ticks], fontsize=7)
        ax.minorticks_off()
        ax.set_ylabel(lab, fontsize=8.5)
        ax.grid(color="0.92", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(-0.15, 1.02, f"({'abcdefgh'[k]})", transform=ax.transAxes, fontsize=10, fontweight="bold")
    for ax in axs[-1]:
        ax.set_xlabel("Cs (log scale)")
    fig.legend(loc="outside lower center", ncol=2, fontsize=8, frameon=False)
    S.save(fig, out, bbox_inches="tight", close=True)
    tt = pd.concat(tables.values())
    tt.round(4).assign(cs=tt.cs).to_csv(out.parent / "supporting" / out.name.replace(".png", ".csv"),
                                               index=False)
    S.record_caption(out, caption + (
        " Crest scores leave out Buxton GIS 9-14 in 2004, 2005 and 2009, where the picked lidar crest reads about "
        "0.8 m low beside the dune rebuilt along the relocated NC-12 (since 2026-11-13)."
        " Scored over 1996-2009: (a-c) dune crest RMSE against the lidar domain medians, surveys after 1996, over "
        "GIS 6-90, GIS 9-15 and the hotspot zones (Buxton 9-14, Avon south 23-26, Rodanthe 82-89); (d) mean model "
        "minus lidar crest in the zones; (e) the same island-wide in 2009; (f) share of domain-surveys within 0.5 m; "
        "(g) overwash false-alarm rate against the imagery; (h) shoreline RMSE against CoastSat. Grey: no collision on "
        "the same starting rows (dashed equal rows, dotted low back row), which has no Cs. Ringed: the Cs the "
        "design's position plot uses now. All runs: toe at the berm (1.7 m NAVD88), growth 0.02-0.04, groin 0.16, "
        "the road fix, no road rule, Rodanthe 40%, village dune lines migrating until the first nourishment; runs "
        "grid_<design>_cs<N>_vdefault (N in 1e-6) and the design runs in 3-collision-design/README.md."))
    print(f"wrote {out}")


# Light to dark per erosion type, one step per Cs in DESIGN_SWEEP_CS; lower both rows yellow-orange-red (an
# orange-only ramp was too alike and turns brown when dark)
CS_RAMPS = {"lowerboth": ("#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026"),
            "lowerfront": ("#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"),
            "faceretreat": ("#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D")}


def ramp(erosion, n):
    """n colours light to dark along the erosion type's ramp (interpolated when a design has more Cs values)."""
    from matplotlib.colors import LinearSegmentedColormap, to_hex
    cmap = LinearSegmentedColormap.from_list(erosion, CS_RAMPS[erosion])
    return [to_hex(cmap(i / max(n - 1, 1))) for i in range(n)]


def positions():
    """Each module design's shoreline position change at every Cs, in one full-style figure."""
    import HAT_collision_on_off_position as P
    for d in CD.DESIGNS:
        if d.erosion == "nocollision":
            continue
        ref = reference(d)
        runs = [(ref.run, ref.name + " (reference)", "#7F7F7F", 1.4, ref.linestyle)]
        name = d.name.rsplit(", Cs", 1)[0]
        cs_all = H.sweep_cs(key(d))
        for c, colour in zip(cs_all, ramp(d.erosion, len(cs_all))):
            v = H.DESIGN_SWEEP[(key(d), c)]
            if R._done(v):
                cur = v == d.run
                runs.append((v, f"Cs {cs_label(c)}" + (" (current)" if cur else ""), colour, 2.6 if cur else 1.5))
        out, gis, ch = P.figure(
            1996, 2009, "HAT_dune_collision", tuple(runs),
            f"Shoreline position change, {name}, Cs {cs_label(cs_all[0])} to {cs_label(cs_all[-1])}, 1996–2009",
            str((d.folder / f"{key(d)}_cs_sweep_position_1996_2009.png").relative_to(P.OUT)),
            f"the design {name} at each Cs of the sweep, light to dark from Cs {cs_label(cs_all[0])} to "
            f"{cs_label(cs_all[-1])} (all solid; thicker: the Cs its own position plot uses now), against "
            f"{ref.name} (grey, underneath; runs grid_{key(d)}_cs<N>_vdefault and {d.run}). " + P.DESIGN_CAPTION,
            post=P._design_zone_shade, legend_rows=5)
        print(f"wrote {out}")
        base = ch[ref.run]
        for v, lab, *_ in runs[1:]:
            dd = np.asarray(ch[v] - base)
            print(f"  {lab:18s} minus reference: mean {dd.mean():+.2f} m, range {dd.min():+.1f} to {dd.max():+.1f} m")


TIMING_REGIONS = {"buxton": ("GIS9-15", range(9, 16)), "avon_south": ("GIS21-27", range(21, 28)),
                  "rodanthe": ("GIS80-89", range(80, 90))}


def timing(region="buxton"):
    """0-all-designs/crest_timing_GIS9-15_faceretreat_vs_lowerfront_1996_2009.png: the domain-median dune crest
    through time at GIS 9-15 against the lidar, no collision against face retreat (Cs 1e-5) and lower front row
    (Cs 2e-5 to 2e-4), all on equal rows: do the drops come when the lidar's do (1998-99, Bonnie and Dennis)?"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from cascade_pipeline.saved_runs import load_saved_cascade
    from site_layer import hat_figure_style as S
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    S.apply_style()
    tag, gis_list = TIMING_REGIONS[region]
    runs = (("cal_step2_growth_noscrape_vdefault", "No collision, equal rows", "#7F7F7F", "-"),
            (H.DESIGN_SWEEP[("faceretreat_equal", 10)], "Face retreat, equal rows, Cs 1e-5", "#6A3D9A", "-"),
            *((H.DESIGN_SWEEP[("lowerfront_equal", c)], f"Lower front row, equal rows, Cs {cs_label(c)}", col, "-")
              for c, col in zip((20, 30, 50, 100, 200), ramp("lowerfront", 7)[2:])))
    sv = pd.read_csv(H.SURVEY_CRESTS)
    sv = sv[(sv.status == "ok") & sv.gis.isin(gis_list)].copy()
    t = pd.to_datetime(sv.flown.astype(str).str[:7], errors="coerce")
    sv["yr"] = t.dt.year + (t.dt.month - 0.5) / 12
    lid = sv.groupby(["gis", "survey"]).agg(z=("crest_m_navd88", "median"), yr=("yr", "first")).reset_index()
    fig, axs = plt.subplots(len(gis_list), 1, figsize=(S.FIG_W_DOUBLE * 0.8, 1.55 * len(gis_list) + 0.2),
                            sharex=True, sharey=True,
                            layout="constrained")
    rows = []
    for v, lab, colour, ls in runs:
        rd = H.run_dir(v)
        c = load_saved_cascade(rd / f"{rd.name}.npz")
        d = DOM.for_total(len(c.barrier3d))
        for ax, g in zip(axs, gis_list):
            b = c.barrier3d[d.gis_to_pad(g)]
            z = np.median((np.asarray(b.DuneDomain, float).max(axis=2) + b._BermEl) * 10 + H.MHW_NAVD88_M, axis=1)
            yrs = H.RUN_START + 1 + np.arange(z.size)
            ax.plot(yrs, z, color=colour, ls=ls, lw=2.2 if "No collision" in lab else 1.9, label=lab,
                    zorder=3 if "No collision" in lab else 4)
            rows += [dict(gis=g, run=v, design=lab, year=int(y), crest_m_navd88=float(zz)) for y, zz in zip(yrs, z)]
    for ax, g in zip(axs, gis_list):
        L = lid[lid.gis == g].sort_values("yr")
        ax.plot(L.yr, L.z, "o", color="#1A1A1A", ms=5, label="Lidar (domain median)", zorder=6)
        for x in (1998.65, 1999.67, 2003.71):
            ax.axvline(x, color="0.6", lw=0.6, ls=":")
        zn = next((n.split(",")[0] for n, r in Z.ZONES_AVN2.items() if g in r), None)
        if zn:
            ax.text(0.995, 0.06, f"{zn} zone", transform=ax.transAxes, ha="right", fontsize=7, color="0.35")
        ax.set_title(f"GIS {g}", loc="left", fontsize=9, fontweight="bold")
        ax.set_ylabel("Crest (m NAVD88)", fontsize=8)
        ax.grid(axis="y", color="0.92", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    for x, n in ((1998.65, "Bonnie"), (1999.67, "Dennis"), (2003.71, "Isabel")):
        axs[0].text(x, 1.04, n, transform=axs[0].get_xaxis_transform(), ha="center", fontsize=7, color="0.4")
    axs[-1].set_xlabel("Year")
    axs[-1].set_xlim(1996, 2010)
    axs[0].set_ylim(0, 9)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=2, fontsize=7.5, frameon=False)
    out = CD.ALL_DIR / f"crest_timing_{tag}_faceretreat_vs_lowerfront_1996_2009.png"
    S.save(fig, out, bbox_inches="tight", close=True)
    pd.DataFrame(rows).to_csv(out.parent / "supporting" / out.name.replace(".png", ".csv"), index=False)
    S.record_caption(out, (
        "Domain-median dune crest (the higher dune row) through 1996-2009 at GIS " + f"{gis_list[0]}-{gis_list[-1]}" + ", model state at 1 January each "
        "year, against the lidar domain median at each survey (black dots): no collision (grey), face retreat (Cs "
        "1e-5, violet) and lower front row (Cs 2e-5 to 2e-4, light to dark green), all on equal rows, village dune "
        "lines migrating. Dotted: Bonnie, Dennis, Isabel. Runs: " + "; ".join(f"{lab} ({v})" for v, lab, *_ in runs)
        + "."))
    print(f"wrote {out}")
    tt = pd.DataFrame(rows)
    for g in gis_list:
        L = lid[lid.gis == g].sort_values("yr")
        drop = L.z.diff().idxmin() if len(L) > 1 else None
        msg = f"  GIS {g}: lidar biggest drop {L.z.diff().min():+.1f} m by {L.loc[drop, 'survey']}" if drop is not None else f"  GIS {g}"
        for v, lab, *_ in runs[1:]:
            z = tt[(tt.gis == g) & (tt.run == v)].set_index("year").crest_m_navd88
            dz = z.diff()
            msg += f" | {lab.split(',')[0][:12]} {lab.split('Cs ')[1]}: {dz.min():+.1f} m in {int(dz.idxmin())}"
        print(msg)


def road():
    """0-all-designs/road_rule_check_1996_2009.png (+ CSV): the front-runners (face retreat, lower front row, lower
    both rows, equal rows) at their best Cs and either side, without (dashed) and with (solid) the 198 m road rule,
    no collision both ways as grey lines; Rodanthe at 40% in all."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as S
    S.apply_style()
    rows = []
    for (k, c), v in H.ROAD_TEST.items():
        base = v[: -len("_road")]
        for rule, run in (("off", base), ("on", v)):
            if R._done(run):
                rows.append(dict(key=k, road_rule=rule, **score(run, (c or 0) * 1e-6)))
    t = pd.DataFrame(rows)
    designs = {key(d): d for d in CD.DESIGNS}
    panels = [p for p in PANELS if p[0] in ("crest_rmse_m", "crest_rmse_gis9_15_m", "crest_rmse_zones_m",
                                            "bias_2009_m", "overwash_false_alarm_pct", "shoreline_rmse_m_yr")]
    fig, axs = plt.subplots(3, 2, figsize=(S.FIG_W_DOUBLE, 7.8), sharex=True, layout="constrained")
    for i, (ax, (col, lab)) in enumerate(zip(axs.ravel(), panels)):
        for rule, ls in (("off", (0, (4, 2))), ("on", "-")):
            nc = t[(t.key == "nocollision_equal") & (t.road_rule == rule)]
            if len(nc):
                ax.axhline(nc[col].iloc[0], color="0.5", lw=1.1, ls=ls,
                           label=f"No collision, road rule {rule}" if i == 0 else None)
            for k in ("lowerboth_equal", "lowerfront_equal", "faceretreat_equal"):
                m = t[(t.key == k) & (t.road_rule == rule)].sort_values("cs")
                if len(m):
                    d = designs[k]
                    ax.plot(m.cs, m[col], marker="o", ms=3.5, lw=1.8, ls=ls, color=d.colour,
                            label=f"{d.name.rsplit(', Cs', 1)[0]}, road rule {rule}" if i == 0 else None)
        if col in ZERO_LINE:
            ax.axhline(0, color="0.2", lw=0.7)
        ax.set_xscale("log")
        ticks = sorted({round(c * 1e6) for c in t.cs if c > 0})
        ax.set_xticks([c * 1e-6 for c in ticks])
        ax.set_xticklabels([cs_label(c) for c in ticks], fontsize=7)
        ax.minorticks_off()
        ax.set_ylabel(lab, fontsize=8.5)
        ax.grid(color="0.92", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(-0.15, 1.02, f"({'abcdef'[i]})", transform=ax.transAxes, fontsize=10, fontweight="bold")
    for ax in axs[-1]:
        ax.set_xlabel("Cs (log scale)")
    fig.legend(loc="outside lower center", ncol=2, fontsize=7.5, frameon=False)
    out = CD.ALL_DIR / "road_rule_check_1996_2009.png"
    S.save(fig, out, bbox_inches="tight", close=True)
    t.round(4).assign(cs=t.cs).to_csv(out.parent / "supporting" / out.name.replace(".png", ".csv"), index=False)
    S.record_caption(out, (
        "The three front-running collision designs on equal rows (lower both rows orange, lower front row green, face "
        "retreat violet) at their best Cs and the values either side, without (dashed) and with (solid) the NCDOT "
        "198 m road-dune rule, against no collision both ways (grey); everything else as the Cs sweep (toe at the "
        "berm, growth 0.02-0.04, groin 0.16, the road fix, Rodanthe 40%, village dune lines migrating until the first "
        "nourishment). Crest scores leave out Buxton GIS 9-14 in 2004, 2005 and 2009 (rebuilt dune). (a-c) Crest RMSE "
        "against the lidar, GIS 6-90, GIS 9-15 and the hotspot zones; (d) mean 2009 crest offset; (e) overwash "
        "false alarms; (f) shoreline RMSE against CoastSat. Runs: <run>_road (HAT_dune_collision.py ROAD_TEST)."))
    print(f"wrote {out}")
    cols = ["crest_rmse_m", "crest_rmse_gis9_15_m", "crest_rmse_zones_m", "bias_2009_m", "overwash_false_alarm_pct",
            "shoreline_rmse_m_yr"]
    print(t.assign(N=(t.cs * 1e6).round().astype(int)).pivot_table(index=["key", "N"], columns="road_rule",
                                                                     values=cols).round(3).to_string())


CREST_FIG = Path("figures") / "dune_crests" / "dune_crest_offset_overview_v2.png"


def cs_label(c):
    """1 -> '1e-6', 10 -> '1e-5', 20 -> '2e-5'."""
    m, e = f"{c * 1e-6:.0e}".split("e")
    return f"{m}e{int(e)}"


def stamp(src, out, label, run):
    """Copy a run's figure with a band added on top naming its design and Cs (the run's own figure does not)."""
    from PIL import Image, ImageDraw, ImageFont
    im = Image.open(src).convert("RGB")
    h = max(60, im.height // 18)
    canvas = Image.new("RGB", (im.width, im.height + h), "white")
    canvas.paste(im, (0, h))
    draw = ImageDraw.Draw(canvas)
    fonts = Path(r"C:\Windows\Fonts")
    try:
        bold = ImageFont.truetype(str(fonts / "arialbd.ttf"), int(h * 0.55))
        small = ImageFont.truetype(str(fonts / "arial.ttf"), int(h * 0.36))
    except OSError:
        bold = small = ImageFont.load_default()
    draw.text((int(h * 0.4), h // 2), label, fill="#1A1A1A", font=bold, anchor="lm")
    draw.text((im.width - int(h * 0.4), h // 2), f"run {run}", fill="#666666", font=small, anchor="rm")
    canvas.save(out, dpi=im.info.get("dpi", (300, 300)))


def crests(S):
    """Copy each run's crest-offset overview into its design folder, renamed by design and Cs and labelled with
    them on a band across the top, with a caption."""
    note = (" Each run's own figure (CASCADE run figures/dune_crests/dune_crest_offset_overview_v2.png), copied: "
            "model minus lidar median dune crest by GIS domain and lidar survey, the run's shoreline position change "
            "1996-2009 against CoastSat on top. Toe at the berm, growth 0.02-0.04, groin 0.16, the road fix, no road "
            "rule, Rodanthe 40%, village dune lines migrating until the first nourishment.")
    jobs = [(d.run, d.folder / f"{d.tag}_crest_offset_1996_2009.png", f"{d.name} (run {d.run}), no Cs.",
             f"{d.name} (no collision module)") for d in CD.DESIGNS if d.erosion == "nocollision"]
    for d in CD.DESIGNS:
        if d.erosion == "nocollision":
            continue
        for c in H.sweep_cs(key(d)):
            v = H.DESIGN_SWEEP[(key(d), c)]
            name = d.name.rsplit(", Cs", 1)[0]
            jobs.append((v, d.folder / "cs-sweep" / f"{key(d)}_cs{cs_label(c)}_crest_offset_1996_2009.png",
                         f"{name}, Cs {cs_label(c)} (run {v})" + (", the Cs its position plot uses now." if v == d.run
                                                                     else "."),
                         f"{name}, Cs {cs_label(c)}" + (" (current)" if v == d.run else "")))
    n = 0
    for v, out, what, label in jobs:
        rd = H.run_dir(v)
        if rd is None or not (rd / CREST_FIG).is_file():
            print(f"  missing crest figure: {v}")
            continue
        out.parent.mkdir(parents=True, exist_ok=True)
        try:
            stamp(rd / CREST_FIG, out, label, v)
        except OSError as e:          # open in a viewer on Windows
            print(f"  could not write {out.name} (open elsewhere?): {e}")
            continue
        S.record_caption(out, "Dune crest offset and shoreline change, " + what + note)
        n += 1
    print(f"copied {n} of {len(jobs)} crest figures")


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as S
    S.apply_style()
    designs = [d for d in CD.DESIGNS if d.erosion != "nocollision"]
    tables = {d.tag: table(d) for d in designs}
    for d in designs:
        if tables[d.tag].empty:
            print(f"  skipping {d.name}: no runs yet")
            continue
        out = d.folder / f"{key(d)}_cs_sensitivity_1996_2009.png"
        (out.parent / "supporting").mkdir(parents=True, exist_ok=True)
        draw(S, plt, [d], {d.tag: tables[d.tag]}, out,
             f"Sensitivity of the design {d.name.rsplit(', Cs', 1)[0]} to the collision module's Cs (Larson et al. "
             f"2004), {len(tables[d.tag])} of {len(H.sweep_cs(key(d)))} Cs values run.")
    out = CD.ALL_DIR / "all_designs_cs_sensitivity_1996_2009.png"
    (out.parent / "supporting").mkdir(parents=True, exist_ok=True)
    draw(S, plt, designs, tables, out,
         "Sensitivity of the six collision designs with the module on to Cs (Larson et al. 2004), each in its "
         "colour: orange lower both rows, green lower front row, violet face retreat; solid equal rows, dashed low "
         "back row.")
    crests(S)
    t = pd.concat(tables.values())
    print(t.pivot_table(index="cs", columns="tag", values="crest_rmse_m").round(3).to_string())


if __name__ == "__main__":
    if sys.argv[1:] == ["positions"]:
        positions()
    elif sys.argv[1:] == ["road"]:
        road()
    elif sys.argv[1:2] == ["timing"]:
        for r in (sys.argv[2:] or ["buxton", "avon_south", "rodanthe"]):
            timing(r)
    else:
        main()
        positions()
