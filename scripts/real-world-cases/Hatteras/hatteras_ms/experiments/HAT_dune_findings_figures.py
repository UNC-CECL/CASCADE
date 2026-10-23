"""
Three summary figures of the dune ceiling (Dmaxel) and growth-rate studies.

    python scripts/hatteras_ms/experiments/HAT_dune_findings_figures.py

1. dune_ceiling_reached_by_year_3: under the adopted settings (5.5 m, r 0.55-0.95)
   the dunes reach the ceiling within about 3 years of each start, the low dunes
   with them, so far fewer cells sit below each year's highest run-up than on the
   starting DEM.
2. dune_ceiling_at_measured_maximum: ceilings set from the measured dunes (DEM P95
   8.4 m, DEM maximum 12.1 m, 10-13; starting-DEM P99 9.03 m with slow growth,
   10-19) against the adopted 5.5 m: overwash skill, crest trajectory, 2009 crest.
3. dune_growth_rate_vs_2009_lidar: every growth rate run from the 1996 start
   (no growth to 0.55-0.95, ceilings 5.0 and 5.5 m), the modelled 2009 crest
   against the 2009 lidar on the same 10 m cells, and the overwash skill.

Read-only on the runs of the 2026-10-13 and 2026-10-19 dunes
experiments (managed edgeBE, groin off). The 2009 lidar is the survey on the
1996 dune windows (HAT_survey_dune_crests.py), set against the model state of
1 January 2010, as HAT_dune_crest_vs_surveys.py does; so the crest offsets here
are cell-matched and differ slightly from the grid study's domain-median offsets.
Crests are cached per run in the experiment folder (cache/), since each run file
is ~220 MB.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-20
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_crest_dem_vs_model as D  # noqa: E402  (also puts the overwash scripts on the path)
import HAT_dune_crest_vs_surveys as V  # noqa: E402
import HAT_dune_ceiling_uniform as U  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-20-dune-findings-figures"
EXPS = U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / "dunes"
EXP_DIR = EXPS / TAG.split("/")[1]
RANGE = "2026-10-13-dune-ceiling-range"
DEMG = "2026-10-19-dune-ceiling-dem-growth"
GRID = "2026-10-19-dune-ceiling-growth-grid"
LOW = "2026-10-19-dune-growth-low-rates"
END = {1996: 2009, 2009: 2025}
FIRST = {1996: 1997, 2009: 2010}        # 1 January of the first saved state (= the DEM)
BERM_NAVD88 = 1.7
MHW_NAVD88 = 0.36
DESIGN_NAVD88 = 4.3
AT_CEILING_TOL = 0.05                   # m below the ceiling still counted as at it
CONTROL = (GRID, "control_s0")
# Growth-rate ladder from the 1996 start: label -> run, per ceiling (the same runs the study tables score)
RATES = ["0", "0.01-0.03", "0.03-0.07", "0.05-0.15", "0.10-0.20", "0.15-0.30", "0.25-0.45", "0.35-0.55",
         "0.55-0.95"]
LADDER = {
    5.5: {"0": (LOW, "c550_r000_000_s0"), "0.01-0.03": (LOW, "c550_r001_003_s0"),
          "0.03-0.07": (LOW, "c550_r003_007_s0"), "0.05-0.15": (LOW, "c550_r005_015_s0"),
          "0.10-0.20": (LOW, "c550_r010_020_s0"), "0.15-0.30": (GRID, "c550_r015_030_s0"),
          "0.25-0.45": (GRID, "c550_r025_045_s0"), "0.35-0.55": (GRID, "c550_r035_055_s0"),
          "0.55-0.95": CONTROL},
    5.0: {"0": (LOW, "c500_r000_000_s0"), "0.01-0.03": (LOW, "c500_r001_003_s0"),
          "0.03-0.07": (LOW, "c500_r003_007_s0"), "0.05-0.15": (LOW, "c500_r005_015_s0"),
          "0.10-0.20": (LOW, "c500_r010_020_s0"), "0.15-0.30": (GRID, "c500_r015_030_s0"),
          "0.25-0.45": (GRID, "c500_r025_045_s0"), "0.35-0.55": (GRID, "c500_r035_055_s0")},
}
# Measured-maximum ceilings, 1996 start
MEASURED = {"5.5 m, r 0.55-0.95 (adopted)": CONTROL,
            "8.4 m (DEM P95), r 0.55-0.95": (RANGE, "c840"),
            "12.1 m (DEM max), r 0.55-0.95": (RANGE, "c1210"),
            "9.03 m (DEM P99), r 0.05-0.15": (DEMG, "p99_r005_015"),
            "9.03 m (DEM P99), r 0.15-0.30": (DEMG, "p99_r015_030")}
MEASURED_COLOURS = dict(zip(MEASURED, [fs.INK, "#f4a582", fs.C_1984, "#c2a5cf", fs.C["ACCENT"]]))
START_COLOURS = {1996: fs.C_1984, 2009: fs.C_1997}   # calibration red, test blue
# -----------------------------------------------------------------------------


def run_dir(study, variant, start):
    base = EXPS / study / "runs" / variant / f"{start}_{END[start]}" / "edgeBE"
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json")))
    if not hits:
        raise FileNotFoundError(base)
    return hits[-1]


# Crest per saved state and 10 m cell (m NAVD88) and the run's ceiling, cached
def crests(study, variant, start):
    cache = EXP_DIR / "cache" / f"{study}__{variant}__{start}.npz"
    if cache.exists():
        z = np.load(cache)
        return z["crest"], float(z["ceiling"])
    c = U.load(run_dir(study, variant, start))
    cc = U.cell_crests(c, DOM)
    crest = np.stack(cc, axis=1).reshape(cc[0].shape[0], -1) + BERM_NAVD88
    ceiling = c.barrier3d[DOM.gis_to_pad(1)]._Dmaxel * 10 + MHW_NAVD88
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, crest=crest.astype(np.float32), ceiling=ceiling)
    print(f"cached {cache.name}", flush=True)
    return crest, ceiling


def survey_points(start):
    s = V.surveys(start)
    rows = []
    for year, g in s.groupby("survey"):
        g = g.sort_values(["gis", "cell"])
        rows.append(dict(survey=year, year=int(g.model_state.iloc[0][:4]), kind=g.kind.iloc[0],
                         median=g.crest.median(), p10=g.crest.quantile(0.1), obs=g.crest.to_numpy()))
    return pd.DataFrame(rows)


def lidar_2009():
    p = survey_points(1996)
    return p[p.survey == 2009].iloc[0]


def yearly_runup(start):
    forcing = D.SV.sf.load_forcing()
    storms = D.SV.sf.classify(D.SV.sf.load_record(D.SV.env.DEFAULT_STORM_VARIANT, forcing), forcing)
    s = storms[storms.pi == list(D.SV.ovm.WINDOWS).index((start, END[start]))]
    return s.groupby("calendar_year").rhigh_m.max() + MHW_NAVD88


def domain_median(h):
    return np.median(h.reshape(90, 50), axis=1)


def survey_markers(ax, pts, key, label=True):
    for kind, marker, lab in (("mosaic", "o", "Surveyed, 1 m DEM"), ("raster", "s", "Surveyed, coarse 3-9 m grid")):
        p = pts[pts.kind == kind]
        if len(p):
            ax.plot(p.year, p[key], ls="none", marker=marker, ms=4.5, mfc="white", mec=fs.INK, mew=1.1,
                    label=lab if label else None, zorder=6)


def finish(ax):
    ax.grid(axis="y", color=fs.GRID_C, lw=0.5)
    fs.open_frame(ax)


# --- Figure 1: the adopted settings reach the ceiling within ~3 years ------------
def figure_reached(out):
    fs.apply_style()
    runs = {s: crests(GRID, "control_s0", s) for s in (1996, 2009)}
    fig, axes = plt.subplots(2, 2, figsize=fs.figsize("double", height=5.2), layout="constrained")
    (a, b), (c, d) = axes
    rows = []
    for start, (cr, ceil) in runs.items():
        col, yrs = START_COLOURS[start], FIRST[start] + np.arange(cr.shape[0])
        at = (cr >= ceil - AT_CEILING_TOL).mean(axis=1) * 100
        a.plot(yrs - FIRST[start], at, color=col, lw=1.4, marker="o", ms=3, label=f"{start} start")
        b.plot(yrs, np.percentile(cr, 10, axis=1), color=col, lw=1.4, label=f"Model, {start} start")
        survey_markers(b, survey_points(start), "p10", label=start == 1996)
        ru = yearly_runup(start)
        k = np.array([y - FIRST[start] for y in ru.index if 0 <= y - FIRST[start] < cr.shape[0]])
        ys = FIRST[start] + k
        mod = [(cr[i] < ru[y]).mean() * 100 for i, y in zip(k, ys)]
        dem = [(cr[0] < ru[y]).mean() * 100 for y in ys]
        c.plot(ys, dem, color=col, lw=1.0, ls="--", marker="o", ms=3, mfc="white",
               label=f"{start} DEM held fixed")
        c.plot(ys, mod, color=col, lw=1.4, marker="o", ms=3, label=f"Model, {start} start")
        for t in (0, 3):
            d.hist(cr[t], bins=np.arange(1.5, 10.01, 0.1), histtype="step", density=True, color=col,
                   lw=1.0 if t == 0 else 1.6, ls=":" if t == 0 else "-",
                   label=f"{start} start, {'DEM (year 0)' if t == 0 else 'model year 3'}")
        for t in range(cr.shape[0]):
            rows.append(dict(start=start, year=FIRST[start] + t, model_year=t, ceiling_m_navd88=round(ceil, 2),
                             share_at_ceiling_pct=round(at[t], 1), p10_m_navd88=round(np.percentile(cr[t], 10), 2),
                             median_m_navd88=round(np.median(cr[t]), 2)))
    ceil = runs[1996][1]
    a.set(xlabel="Years since the start DEM", ylabel=f"Cells at the ceiling (%, within {AT_CEILING_TOL} m)",
          ylim=(0, 100))
    a.axvline(3, color=fs.INK_MUTED, lw=0.7, ls=":")
    a.legend(frameon=False, fontsize=7, loc="lower right")
    b.axhline(ceil, color=fs.INK, lw=0.8, ls="--", label=f"Ceiling, Dmaxel {ceil:.1f} m")
    b.axhline(DESIGN_NAVD88, color=fs.C["REF"], lw=0.8, ls=":", label="NC-12 rebuild height")
    b.set(xlabel="1 January of", ylabel="Low dunes, 10th percentile crest (m NAVD88)")
    b.legend(frameon=False, fontsize=6.5, loc="upper left", bbox_to_anchor=(0.0, -0.2), ncol=3)
    c.set(xlabel="Calendar year (that year's highest storm run-up)", ylabel="Cells below the year's highest run-up (%)")
    c.set_ylim(top=62)
    c.legend(frameon=False, fontsize=6.5, loc="upper right", ncol=2)
    d.axvline(ceil, color=fs.INK, lw=0.8, ls="--")
    d.set(xlabel="Dune crest (m NAVD88), 10 m cells", ylabel="Density", xlim=(1.5, 10))
    d.legend(frameon=False, fontsize=6.5, loc="upper right")
    for i, ax in enumerate(axes.flat):
        finish(ax)
        fs._letter_inside(ax, i)
    fs.save(fig, out, close=True)
    fs.record_caption(out, (
        "Under the adopted dune settings (Dmaxel 5.5 m NAVD88, growth 0.55-0.95; managed edgeBE runs, groin off; "
        "1996 start red, 2009 start blue) the dunes reach the ceiling within about three years. (a) share of the "
        "4,500 10 m cells within 0.05 m of the ceiling, by years since the start DEM. (b) 10th percentile crest "
        "(the low dunes) each 1 January, model lines and surveys (circles 1 m DEMs, squares coarse 1997-2005 lidar "
        "grids), on each start's dune windows. (c) share of cells whose crest on 1 January sits below that year's "
        "highest storm run-up: model (solid) and the start DEM held fixed (dashed). (d) crest distribution at the "
        "start (dotted) and at model year 3 (solid); dashed line the ceiling. Producer: HAT_dune_findings_figures.py."))
    return pd.DataFrame(rows)


# --- Figure 2: ceilings from the measured dunes -----------------------------------
def figure_measured(out):
    fs.apply_style()
    rng = pd.read_csv(EXPS / RANGE / "tables" / "scores.csv")
    rng = rng[rng.period == "1996_2009"].sort_values("ceiling_m_navd88")
    p99 = pd.read_csv(EXPS / DEMG / "tables" / "scores.csv")
    p99 = p99[(p99.period == "1996_2009") & p99.variant.str.startswith("p99")]
    runup = yearly_runup(1996)
    isabel = float(runup.loc[2003])
    lid, pts = lidar_2009(), survey_points(1996)

    fig = plt.figure(figsize=fs.figsize("double", height=5.6), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1])
    a, b, c = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])
    a.plot(rng.ceiling_m_navd88, rng.skill, color=fs.INK, lw=1.2, marker="o", ms=3.5,
           label="Growth 0.55-0.95 (10-13 ceiling range)")
    for _, r in p99.iterrows():
        lab = r.variant.replace("p99_r", "r ").replace("_", "-").replace("r 0", "r 0.").replace("-0", "-0.")
        a.plot(9.03, r.skill, ls="none", marker="D", ms=4.5, color=fs.C["ACCENT"], mfc="white", mew=1.1)
        a.annotate(lab.replace("..", "."), (9.03, r.skill), xytext=(5, 0), textcoords="offset points",
                   fontsize=6, va="center", color=fs.C["ACCENT"])
    a.plot([], [], ls="none", marker="D", ms=4.5, color=fs.C["ACCENT"], mfc="white", mew=1.1,
           label="9.03 m (DEM P99), slow growth (10-19)")
    a.axvline(isabel, color=fs.C_1997, lw=0.8, ls=":")
    a.text(isabel, 0.72, f" Isabel run-up {isabel:.2f} m", fontsize=6.5, color=fs.C_1997, va="top")
    for x, t in ((8.4, "DEM P95"), (12.1, "DEM max")):
        a.annotate(t, (x, rng.set_index("ceiling_m_navd88").skill.get(x, 0.02)), xytext=(0, 7),
                   textcoords="offset points", ha="center", fontsize=6.5)
    a.set(xlabel="Dune ceiling, Dmaxel (m NAVD88)", ylabel="1996-2009 overwash skill (hit rate − false-alarm rate)",
          ylim=(0, 0.75))
    a.legend(frameon=False, fontsize=6.5, loc="upper right", bbox_to_anchor=(1.0, 0.9))

    final = {}
    for label, (study, variant) in MEASURED.items():
        cr, _ = crests(study, variant, 1996)
        yrs = FIRST[1996] + np.arange(cr.shape[0])
        b.plot(yrs, np.median(cr, axis=1), color=MEASURED_COLOURS[label], lw=1.4, label=label)
        k = lid.year - FIRST[1996]
        final[label] = cr[k]
        c.step(np.arange(1, 91), domain_median(cr[k]), where="mid", color=MEASURED_COLOURS[label], lw=1.2)
    survey_markers(b, pts, "median")
    b.set(xlabel="1 January of", ylabel="Median dune crest (m NAVD88)")

    obs = lid.obs.reshape(90, 50)
    full = np.isfinite(obs).sum(1) >= 25
    od = np.where(full, np.nanmedian(np.where(np.isfinite(obs), obs, np.nan), axis=1), np.nan)
    c.step(np.arange(1, 91), od, where="mid", color=fs.INK, lw=2.2, alpha=0.35, label="2009 lidar", zorder=1)
    c.axhline(DESIGN_NAVD88, color=fs.C["REF"], lw=0.8, ls=":")
    c.set(xlim=(0.5, 90.5), xlabel=fs.DOMAIN_AXIS_LABEL + "; 500 m domain medians",
          ylabel="Crest, 1 Jan 2010 (m NAVD88)")
    fs.town_bands(c)
    h, l = b.get_legend_handles_labels()
    h2, l2 = c.get_legend_handles_labels()
    fig.legend(h + h2, l + l2, loc="outside lower center", ncol=4, frameon=False, fontsize=6.5)
    rows = []
    for label, m in final.items():
        cmp = V.compare(lid.obs, m)
        rows.append(dict(run=label, median_offset_vs_2009_lidar_m=round(cmp["median_offset_m"], 2),
                         domain_median_r=round(cmp["domain_median_r"], 2)))
    for i, ax in enumerate((a, b, c)):
        finish(ax)
        fs._letter_inside(ax, i)
    fs.save(fig, out, close=True)
    fs.record_caption(out, (
        "Dune ceilings set from the measured dunes, 1996-2009, managed edgeBE runs, groin off. (a) overwash skill "
        "against Dmaxel with the adopted growth 0.55-0.95 (2026-10-13 range: 4.5 m to the 1996 DEM maximum 12.1 m; "
        "8.4 m is the DEM 95th percentile) and with the ceiling at the starting-DEM 99th percentile (9.03 m) and "
        "slower growth (2026-10-19, diamonds). Skill collapses once the ceiling passes Isabel's run-up (dotted): "
        "the dunes grow above the only storm able to overtop them. (b) median crest each 1 January, model lines "
        "and surveys (circles 1 m DEMs, squares coarse 1997-2005 lidar grids). (c) 500 m domain-median crest on "
        "1 January 2010 against the 2009 lidar on the same 10 m cells (grey; domains surveyed in at least half "
        "their cells); dotted line the NC-12 rebuild height. Producer: HAT_dune_findings_figures.py."))
    return pd.DataFrame(rows)


# --- Figure 3: growth rate against the 2009 lidar ---------------------------------
def figure_growth(out):
    fs.apply_style()
    lid = lidar_2009()
    k = lid.year - FIRST[1996]
    scores = pd.concat([pd.read_csv(EXPS / LOW / "tables" / "stage1_1996.csv"),
                        pd.read_csv(EXPS / GRID / "tables" / "stage1_grid_1996.csv")]).drop_duplicates("variant")
    scores = scores.set_index("variant")
    o = lid.obs
    ok = np.isfinite(o)
    rows = []
    for ceil, ladder in LADDER.items():
        for rate, (study, variant) in ladder.items():
            cr, _ = crests(study, variant, 1996)
            m = cr[k]
            cmp = V.compare(o, m)
            rows.append(dict(ceiling_m_navd88=ceil, growth=rate, run=f"{study}/{variant}",
                             median_offset_m=cmp["median_offset_m"], domain_median_r=cmp["domain_median_r"],
                             p10_error_m=np.percentile(m[ok], 10) - np.percentile(o[ok], 10),
                             mean_abs_error_m=cmp["mean_abs_error_m"],
                             skill=scores.skill.get(variant, np.nan), hits=scores.hits.get(variant, np.nan),
                             false_alarms=scores.false_alarms.get(variant, np.nan)))
    t = pd.DataFrame(rows)

    fig = plt.figure(figsize=fs.figsize("double", height=6.0), layout="constrained")
    gs = fig.add_gridspec(2, 4, height_ratios=[1.15, 1])
    top = fig.add_subplot(gs[0, :])
    axs = [fig.add_subplot(gs[1, i]) for i in range(4)]
    obs = o.reshape(90, 50)
    full = np.isfinite(obs).sum(1) >= 25
    od = np.where(full, np.nanmedian(np.where(np.isfinite(obs), obs, np.nan), axis=1), np.nan)
    top.step(np.arange(1, 91), od, where="mid", color=fs.INK, lw=2.4, alpha=0.35, label="2009 lidar", zorder=1)
    show = ["0", "0.10-0.20", "0.25-0.45", "0.55-0.95"]
    ramp = plt.cm.Purples(np.linspace(0.45, 1.0, len(show)))
    for col, rate in zip(ramp, show):
        cr, _ = crests(*LADDER[5.5][rate], 1996)
        lab = "no growth" if rate == "0" else f"r {rate}" + (" (adopted)" if rate == "0.55-0.95" else "")
        top.step(np.arange(1, 91), domain_median(cr[k]), where="mid", color=col, lw=1.1, label=f"5.5 m, {lab}")
    top.axhline(5.5, color=fs.INK, lw=0.7, ls="--")
    top.set(xlim=(0.5, 90.5), xlabel=fs.DOMAIN_AXIS_LABEL + "; 500 m domain medians",
            ylabel="Crest, 1 Jan 2010 (m NAVD88)")
    top.legend(frameon=False, fontsize=6.5, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=5)
    fs.town_bands(top)

    x = np.arange(len(RATES))
    for ax, key, ylab, ref in ((axs[0], "median_offset_m", "Median crest offset,\nmodel − lidar (m)", 0),
                               (axs[1], "p10_error_m", "Low dunes (P10),\nmodel − lidar (m)", 0),
                               (axs[2], "domain_median_r", "Alongshore pattern r\n(domain medians)", None),
                               (axs[3], "skill", "1996-2009 overwash skill", None)):
        for ceil, ls, mk in ((5.5, "-", "o"), (5.0, "--", "s")):
            s = t[t.ceiling_m_navd88 == ceil].set_index("growth").reindex(RATES)
            ax.plot(x, s[key], ls=ls, marker=mk, ms=3.5, color=fs.C["ACCENT"] if ceil == 5.0 else fs.INK,
                    mfc="white", lw=1.1, label=f"Dmaxel {ceil} m")
        if ref is not None:
            ax.axhline(ref, color=fs.INK_MUTED, lw=0.7)
        ax.set_xticks(x, ["none" if r == "0" else r for r in RATES], rotation=60, ha="right", fontsize=6)
        ax.set_ylabel(ylab, fontsize=7)
        ax.set_xlabel("Growth rate, rmin-rmax", fontsize=7)
    axs[3].plot(x[-1], t[(t.ceiling_m_navd88 == 5.5) & (t.growth == "0.55-0.95")].skill, marker="*", ms=9,
                color=fs.C_1984, ls="none", label="Adopted")
    axs[0].legend(frameon=False, fontsize=6.5, loc="lower right")
    axs[3].legend(frameon=False, fontsize=6.5, loc="lower right")
    for i, ax in enumerate([top] + axs):
        finish(ax)
        fs._letter_inside(ax, i)
    fs.save(fig, out, close=True)
    fs.record_caption(out, (
        "Which growth rate from the 1996 start reproduces the 2009 dunes? Managed edgeBE runs, groin off, "
        "Barrier3D's default growth-rate draw; the 2009 lidar is taken on the 1996 dune windows and set against "
        "the model state of 1 January 2010, cell by cell. (a) 500 m domain-median crest: lidar (grey) and the "
        "5.5 m ceiling with four growth rates (light to dark: slow to the adopted 0.55-0.95). (b) median of "
        "model minus lidar over the cells; (c) the same for the 10th percentile (the low dunes); (d) alongshore "
        "correlation of the domain medians; (e) 1996-2009 overwash skill from the study tables. Solid: Dmaxel "
        "5.5 m; dashed: 5.0 m. Slower growth follows the 2009 dunes better, but below about 0.25 the extra "
        "overwash is mostly false alarms. Producer: HAT_dune_findings_figures.py."))
    return t


def main():
    fig_dir, tab_dir = EXP_DIR / "figures" / "5-detailed", EXP_DIR / "tables"
    tab_dir.mkdir(parents=True, exist_ok=True)
    jobs = [("dune_ceiling_reached_by_year_3_1996_2009_starts", figure_reached),
            ("dune_ceiling_at_measured_maximum_1996_2009", figure_measured),
            ("dune_growth_rate_vs_2009_lidar_1996_start", figure_growth)]
    for name, f in jobs:
        t = f(fig_dir / f"{name}.png")
        t.round(3).to_csv(tab_dir / f"{name}.csv", index=False)
        with pd.option_context("display.width", 220, "display.max_columns", 20):
            print(t.round(3).to_string(index=False))
        print(f"figure: {fig_dir / name}.png", flush=True)


if __name__ == "__main__":
    main()
