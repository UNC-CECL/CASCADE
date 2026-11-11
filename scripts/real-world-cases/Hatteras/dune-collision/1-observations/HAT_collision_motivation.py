"""
Why the dune collision module is needed: one figure from the 1996-2009 lidar and the model without it.

    python scripts/real-world-cases/Hatteras/dune-collision/1-observations/HAT_collision_motivation.py profiles
    python scripts/real-world-cases/Hatteras/dune-collision/1-observations/HAT_collision_motivation.py plot

profiles   the 1 m lidar profiles of one Buxton cell before and after Dennis (Sep 1998, Sep 1999), from the
           survey rasters as HAT_survey_dune_crests reads them; written once to profiles_dennis.csv (slow)
plot       the four panels:
           (a) the mechanism: Dennis's peak water level reaches the dune face but not the crest, and the
               lidar crest drops while the model's, with no overtopping, does not move
           (b) when, at Buxton: GIS 10-12's lidar crests against the model's, 1996-2009
           (c) when, island-wide: the share of domains whose modelled crest stands more than 1 m above
               (or below) the lidar's, survey by survey
           (d) what sets it: crest loss in collision against beach width, Dennis and 2004-05

The model is the calibration's control (cal_control_g16: the 2026-10-29 defaults, Hs 3.0, StormStart 1, groin
0.16, the module off). Lidar crests: survey_dune_crests.csv (status ok). Toes, beach widths and crest loss in
collision: the 2026-10-26-dune-collision tables (dune_toes.csv, collision_calibration_cells.csv). The beach
width in (d) is each cell's median over the surveys, 1998 left out, as the module's: the 1998 grid's own
widths put most Buxton cells behind beaches wider than 60 m.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-03
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in sys.path:
    sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
# One folder per topic in 1-observations/, each with its PNGs and a supporting/ (2026-11-03); this script draws 3-mechanism
ROOT_DIR = PROJECT_ROOT / "dune-collision" / "1-observations"
OUT_DIR = ROOT_DIR / "3-mechanism"
RUNS = PROJECT_ROOT / "output" / "raw_runs" / "dune-storm-response"
COLLISION_DIR = RUNS / "2026-10-26-dune-collision"
CONTROL_NPZ = (COLLISION_DIR / "runs" / "cal_control_g16" / "1996_2009" / "edgeBE" / "HAT_1996-2009_managed_groin_edgeBE"
               / "HAT_1996-2009_managed_groin_edgeBE.npz")
START = "1996-start"
SURVEY_CRESTS = (PROJECT_ROOT / "data" / "hatteras_init" / "1-barrier3d-domains" / START
                 / "survey-crests" / "survey_dune_crests.csv")
MHW_NAVD88_M = 0.36
FIRST_MODEL_YEAR = 1997          # DuneDomain row 0 is 1 Jan 1997
# (a) the example cell: narrow beach, never overtopped by Dennis, about 2.5 m of crest lost; the 1999 survey
# was flown after Dennis and before Floyd
PROFILE_GIS, PROFILE_CELL = 12, 22
PROFILE_SURVEYS = (1998, 1999)
SPIKE_M = 20.0                   # 1 m returns above this are lidar spikes (the 1998 grid reaches 50 m)
PROFILE_HALF_M = 60              # cross-shore half-window around the 1998 crest
# (b) the Buxton hotspot; storm landfalls as decimal years
BUXTON_GIS = (10, 12)
STORMS = {"Bonnie": 1998.65, "Dennis": 1999.67, "Isabel": 2003.71}
# (c) a domain counts if this many of its cells have a lidar crest
MIN_CELLS = 10
MISMATCH_M = 1.0
# (d) the storm intervals, their peak water level from collision_calibration.csv
D_INTERVALS = ("1998-1999", "2004-2005")
D_BINS_M = np.arange(0, 101, 10)
D_MIN_CELLS = 15
# Colours: lidar and model as the crest figures; the storm water and the 2004-05 interval on the orange ramp;
# model-too-high from the crest-offset map's red pole (checked for colour-blind separation, 2026-11-03)
LIDAR_C, MODEL_C, WATER_C, LATER_C, HIGH_C = "#2166AC", "#8C8C8C", "#D94801", "#FDAE6B", "#B2182B"
# The two-dune figure (two_dunes): one storm, two dunes alike in crest and toe, behind a narrow and a wide beach.
# Chosen 2026-11-03 from collision_calibration_cells.csv, Dennis (1998-99): crest 5.87 / 5.86 m, toe 2.32 / 2.28 m
# NAVD88 (10 m cell values), dry beach 10 / 40 m (to the crest 25 / 70 m), crest loss 1.80 / 0.01 m; profiles free of the
# 1998 grid's spikes. GIS 12 lies where NC-12 was relocated in 1999, but its change is on the seaward face only
TWO_DUNES = {"narrow": (12, 22), "wide": (21, 33)}       # GIS, cell
# Each dune's toe and dry beach width (beach start to toe) as the model reads them, and the module's limit on it
# Each cell's own median over the surveys (1998 left out): the per-cell table the module read until 2026-11-04; it
# now reads one toe and one beach factor per domain, so the live dune_toe_cells.csv no longer holds these values
TOE_CELLS = (PROJECT_ROOT / "data" / "hatteras_init" / "9-dune-toes" / "0-previous-methods" / START
             / "per-cell-median-until-2026-11-04" / "dune_toe_cells.csv")
W0_M = 25.0                      # hat_run.yaml dune_collision beach_width_limit_m, on the dry beach since 2026-11-03
TWO_DUNES_X = (-110, 40)         # cross-shore window, m from each dune's 1998 crest
# -----------------------------------------------------------------------------

SUPPORT_DIR = OUT_DIR / "supporting"     # the house rule: PNGs at the top, PDFs and tables under supporting/
PROFILES_CSV = SUPPORT_DIR / "profiles_dennis.csv"
TWO_DUNES_CSV = SUPPORT_DIR / "profiles_two_dunes.csv"


def cell_profiles(cells):
    """1 m cross-shore profiles of each (GIS, cell) in the 1998 and 1999 surveys, ocean first, m NAVD88: the
    median of the cell's ten 1 m rows, as HAT_survey_dune_crests reads the surveys. Columns z_<gis>_<cell>_<year>."""
    import warnings
    import rasterio
    warnings.filterwarnings("ignore")
    sys.path.insert(0, str(_HAT_DIR / "input_prep" / "1-barrier3d-domains" / "2-survey-crests"))
    import HAT_survey_dune_crests as S
    X = S.X
    prod, surveys = S.STARTS[START]
    clips = S.product(prod).gapfill_1m
    # unit_of reads GIS 10, 30, 50 and 70 to tell feet from metres
    grids = {g: S.clip_grid(clips / f"clip_domain_{g}_filled.tif") for g in {g for g, _ in cells} | {10, 30, 50, 70}}
    cols = {}
    for year in PROFILE_SURVEYS:
        kind, path, flown, state, note = surveys[year]
        files = sorted(path.parent.glob(path.name)) if "*" in path.name else [path]
        srcs = [rasterio.open(f) for f in files]
        try:
            unit = "m" if kind == "mosaic" else S.unit_of(srcs, START, grids)[0]
            scale = S.US_FT if unit == "ft" else 1.0
            for g, cell in cells:
                a = (S.raster_on_clip(srcs, grids[g]) * scale + S.SURVEY_VERTICAL_OFFSET_M[START].get(year, 0.0))
                a = X.orient_ocean_right(a, X.OCEAN_LOC)[:, ::-1]       # rows alongshore, ocean first, 1 m
                rows = a[10 * cell:10 * cell + 10]
                cols[f"z_{g}_{cell}_{year}"] = np.nanmedian(np.where(rows < SPIKE_M, rows, np.nan), axis=0)
        finally:
            for s in srcs:
                s.close()
        print(f"{year}: {unit}", flush=True)
    t = pd.DataFrame(cols)
    t.insert(0, "x_m", np.arange(len(t)))
    return t


def profiles():
    """The example cell's profiles for panel (a), written once to PROFILES_CSV (slow: reads the survey rasters)."""
    t = cell_profiles([(PROFILE_GIS, PROFILE_CELL)])
    t = t.rename(columns={f"z_{PROFILE_GIS}_{PROFILE_CELL}_{y}": f"z_{y}_m_navd88" for y in PROFILE_SURVEYS})
    SUPPORT_DIR.mkdir(parents=True, exist_ok=True)
    t.round(3).to_csv(PROFILES_CSV, index=False)
    print(f"wrote {PROFILES_CSV}")


def model_crests():
    """The control's crest per year, domain and cell (m NAVD88), keyed by GIS."""
    from cascade_pipeline.saved_runs import load_saved_cascade
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS
    c = load_saved_cascade(CONTROL_NPZ)
    dom = HATTERAS_DOMAINS.for_total(len(c.barrier3d))
    out = {}
    for g in range(1, 91):
        b = c.barrier3d[dom.gis_to_pad(g)]
        out[g] = (np.asarray(b.DuneDomain, float).max(axis=2) + b._BermEl) * 10 + MHW_NAVD88_M
    return out


def lidar_crests():
    s = pd.read_csv(SURVEY_CRESTS)
    s = s[s.status == "ok"].copy()
    s["row"] = pd.to_datetime(s.model_state).dt.year - FIRST_MODEL_YEAR
    span = s.flown.str.split(" to ")
    first, last = pd.to_datetime(span.str[0]), pd.to_datetime(span.str[-1])
    s["when"] = (first.dt.year + (first.dt.month + last.dt.month) / 2 / 12
                 + (last.dt.year - first.dt.year) * last.dt.month / 24)
    return s


def module_beach_width():
    """Each cell's median beach width over the surveys, 1998 left out, as the module's (m)."""
    t = pd.read_csv(COLLISION_DIR / "dune_toes.csv")
    t = t[(t.status == "ok") & ~t.survey.astype(str).isin(["start", "1998"])]
    return t.groupby(["gis", "cell"]).beach_width_m.median().rename("w_module_m")


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FormatStrFormatter, MultipleLocator, PercentFormatter
    from site_layer.hat_figure_style import INK, INK_MUTED, apply_style, record_caption
    apply_style()

    if not PROFILES_CSV.exists():
        raise SystemExit(f"run the profiles action first ({PROFILES_CSV} missing)")
    prof = pd.read_csv(PROFILES_CSV)
    s = lidar_crests()
    model = model_crests()
    storm_peak = pd.read_csv(COLLISION_DIR / "collision_calibration.csv").set_index("interval").twl_max_m_navd88
    dennis = float(storm_peak["1998-1999"])
    cells = pd.read_csv(COLLISION_DIR / "collision_calibration_cells.csv")

    fig, axs = plt.subplots(2, 2, figsize=(7.4, 6.4))

    def title(ax, letter, text):
        ax.set_title(f"({letter})", loc="left", fontweight="bold")
        ax.set_title(text, loc="center")

    # (a) the mechanism, one cell
    ax = axs[0, 0]
    z98, z99 = prof[f"z_{PROFILE_SURVEYS[0]}_m_navd88"].to_numpy(), prof[f"z_{PROFILE_SURVEYS[1]}_m_navd88"].to_numpy()
    c98 = int(np.nanargmax(z98[: prof.x_m.size // 2]))        # the foredune crest, in the seaward half
    x = prof.x_m.to_numpy() - c98                             # m from the 1998 crest, ocean negative
    keep = np.abs(x) <= PROFILE_HALF_M
    toe = float(cells[(cells.gis == PROFILE_GIS) & (cells.cell == PROFILE_CELL)].toe_cell_m_navd88.iloc[0])
    # the water reaches the 1998 dune face where the face first climbs above it
    face = np.flatnonzero((x < 0) & (z98 >= dennis))
    x_reach = x[face[0]] if face.size else 0
    wet = keep & (x <= x_reach)
    ax.fill_between(x[wet], z98[wet], dennis, where=z98[wet] < dennis, color=WATER_C, alpha=0.12, lw=0)
    ax.plot([-PROFILE_HALF_M, x_reach], [dennis, dennis], color=WATER_C, lw=1.6)
    ax.text(-PROFILE_HALF_M + 2, dennis + 0.15, f"Dennis peak water level {dennis:.2f} m", color=INK,
            fontsize=7.5, va="bottom")
    ax.plot(x[keep], z98[keep], color=LIDAR_C, lw=1.4, ls=(0, (4, 2)), label="Lidar, Sep 1998 (before)")
    ax.plot(x[keep], z99[keep], color=LIDAR_C, lw=2.0, label="Lidar, Sep 1999 (after)")
    ax.axhline(toe, color=INK_MUTED, lw=0.7, ls=":")
    ax.text(-PROFILE_HALF_M + 2, toe + 0.08, f"dune toe {toe:.1f} m", color=INK_MUTED, fontsize=7, ha="left",
            va="bottom")
    m98, m99 = (float(model[PROFILE_GIS][1999 - FIRST_MODEL_YEAR, PROFILE_CELL]),
                float(model[PROFILE_GIS][2000 - FIRST_MODEL_YEAR, PROFILE_CELL]))
    ax.plot([8, 26], [m99, m99], color=MODEL_C, lw=2.2, solid_capstyle="round")
    ax.text(28, m99, f"model crest\n{m98:.1f} → {m99:.1f} m", color=INK, fontsize=7, va="center")
    ax.annotate("", xy=(3, np.nanmax(z99[keep]) + 0.1), xytext=(3, np.nanmax(z98[keep]) - 0.1),
                arrowprops=dict(arrowstyle="->", color=INK, lw=0.9))
    ax.text(8, (np.nanmax(z98[keep]) + np.nanmax(z99[keep])) / 2,
            f"−{np.nanmax(z98[keep]) - np.nanmax(z99[keep]):.1f} m", fontsize=7.5, color=INK, va="center")
    ax.set_xlim(-PROFILE_HALF_M, PROFILE_HALF_M)
    ax.set_ylim(-1.0, 8.6)
    ax.set_xlabel("Cross-shore distance from the 1998 crest (m, ocean left)")
    ax.set_ylabel("Elevation (m NAVD88)")
    ax.legend(loc="upper left", fontsize=7)
    title(ax, "a", f"Dennis at GIS {PROFILE_GIS}")

    # (b) when, at Buxton
    ax = axs[0, 1]
    lo_g, hi_g = BUXTON_GIS
    b = s[s.gis.between(lo_g, hi_g)]
    mod_rows = np.concatenate([model[g] for g in range(lo_g, hi_g + 1)], axis=1)   # (years, cells)
    years = FIRST_MODEL_YEAR + np.arange(mod_rows.shape[0])
    q = np.percentile(mod_rows, [25, 50, 75], axis=1)
    ax.fill_between(years, q[0], q[2], color=MODEL_C, alpha=0.22, lw=0)
    ax.plot(years, q[1], color=MODEL_C, lw=2.0, label="Model without the module")
    for sv, grp in b.groupby("survey"):
        qq = grp.crest_m_navd88.quantile([0.25, 0.5, 0.75]).to_numpy()
        ax.errorbar(grp.when.iloc[0], qq[1], yerr=[[qq[1] - qq[0]], [qq[2] - qq[1]]], fmt="o", ms=4.5,
                    color=LIDAR_C, mfc=LIDAR_C if sv != 2001 else "white", elinewidth=1, capsize=0, zorder=5,
                    label="Lidar" if sv == b.survey.min() else None)
    ax.axhline(dennis, color=WATER_C, lw=1.0, ls=(0, (4, 2)))
    ax.text(2009.4, dennis - 0.12, "Dennis peak", color=INK, fontsize=7, ha="right", va="top")
    for name, t in STORMS.items():
        ax.axvline(t, color=INK_MUTED, lw=0.6, ls=":")
        ax.text(t - 0.1, 8.45, name, rotation=90, ha="right", va="top", fontsize=6.5, color=INK_MUTED)
    ax.set_xlim(1996.3, 2010.2)
    ax.xaxis.set_major_locator(FixedLocator(range(1996, 2011, 2)))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
    ax.set_ylim(1.0, 8.6)
    ax.set_ylabel("Dune crest (m NAVD88)")
    ax.legend(loc="lower left", fontsize=7)
    title(ax, "b", f"Dune crests, GIS {lo_g}–{hi_g}")

    # (c) when, island-wide
    ax = axs[1, 0]
    rows = []
    for (sv, g), grp in s.groupby(["survey", "gis"]):
        if len(grp) < MIN_CELLS:
            continue
        r = int(grp.row.iloc[0])
        rows.append(dict(survey=sv, when=grp.when.iloc[0], gis=g,
                         offset_m=float(np.median(model[g][r][grp.cell.to_numpy()]) - grp.crest_m_navd88.median())))
    d = pd.DataFrame(rows)
    share = d.groupby(["survey", "when"]).offset_m.agg(
        high=lambda v: (v > MISMATCH_M).mean(), low=lambda v: (v < -MISMATCH_M).mean(), n="size").reset_index()
    ax.plot(share.when, share.high, "o-", color=HIGH_C, lw=1.8, ms=4.5, label=f"Model > {MISMATCH_M:g} m above lidar")
    ax.plot(share.when, share.low, "s--", color=INK, lw=1.2, ms=3.5, mfc="white",
            label=f"Model > {MISMATCH_M:g} m below lidar")
    for name, t in STORMS.items():
        ax.axvline(t, color=INK_MUTED, lw=0.6, ls=":")
    ax.set_xlim(1996.3, 2010.2)
    ax.xaxis.set_major_locator(FixedLocator(range(1996, 2011, 2)))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
    ax.set_ylim(0, max(0.3, share[["high", "low"]].max().max() * 1.15))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax.set_ylabel("Share of domains")
    ax.set_xlabel("Lidar survey")
    ax.legend(loc="upper left", fontsize=7)
    title(ax, "c", "Model against lidar, island-wide")

    # (d) what sets the loss
    ax = axs[1, 1]
    c = cells.join(module_beach_width(), on=["gis", "cell"])
    c["r_high"] = c.interval.map(storm_peak)
    c = c[c.interval.isin(D_INTERVALS) & (c.toe_cell_m_navd88 < c.r_high) & (c.crest_m_navd88 > c.r_high)]
    mids = (D_BINS_M[:-1] + D_BINS_M[1:]) / 2
    d_rows = []
    for iv, col, lw, band in ((D_INTERVALS[0], WATER_C, 2.0, True), (D_INTERVALS[1], LATER_C, 1.6, False)):
        x_iv = c[c.interval == iv]
        grp = x_iv.groupby(pd.cut(x_iv.w_module_m, D_BINS_M, right=False), observed=False).crest_loss_m
        n = grp.size().to_numpy()
        qs = [grp.quantile(p).to_numpy() for p in (0.25, 0.5, 0.75)]
        ok = n >= D_MIN_CELLS
        if band:
            ax.fill_between(mids[ok], qs[0][ok], qs[2][ok], color=col, alpha=0.15, lw=0)
        name = "Dennis, 1998–99" if iv == "1998-1999" else "2004–05"
        ax.plot(mids[ok], qs[1][ok], "o-", color=col, lw=lw, ms=4.5,
                label=f"{name} (peak {storm_peak[iv]:.2f} m)")
        d_rows += [dict(interval=iv, beach_width_bin_m=f"{lo}-{lo + 10}", cells=k, loss_p25_m=a, loss_median_m=m,
                        loss_p75_m=hh) for lo, k, a, m, hh in zip(D_BINS_M[:-1], n, *qs)]
    ax.axhline(0, color=INK_MUTED, lw=0.7)
    ax.axvline(45, color=INK_MUTED, lw=0.7, ls=":")
    ax.text(46, ax.get_ylim()[0] + 0.05, "45 m", color=INK_MUTED, fontsize=7, va="bottom")
    ax.set_xlim(0, 100)
    ax.set_xlabel("Beach width in front of the dune (m)")
    ax.set_ylabel("Crest lost between surveys (m)")
    ax.legend(loc="upper right", fontsize=7)
    title(ax, "d", "Crest loss and the beach")

    fig.tight_layout(h_pad=1.6, w_pad=1.4)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "collision_motivation.png"
    fig.savefig(out, dpi=300, facecolor="white")
    SUPPORT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(SUPPORT_DIR / "collision_motivation.pdf", facecolor="white")
    plt.close(fig)

    share.round(3).to_csv(SUPPORT_DIR / "collision_motivation_share_by_survey.csv", index=False)
    d.round(3).to_csv(SUPPORT_DIR / "collision_motivation_model_minus_lidar_by_domain.csv", index=False)
    pd.DataFrame(d_rows).round(3).to_csv(SUPPORT_DIR / "collision_motivation_loss_by_beach_width.csv", index=False)
    print(share.round(2).to_string(index=False))
    print(f"(a) GIS {PROFILE_GIS} cell {PROFILE_CELL}: lidar crest {np.nanmax(z98[keep]):.2f} -> "
          f"{np.nanmax(z99[keep]):.2f} m, model {m98:.2f} -> {m99:.2f} m, toe {toe:.2f}, Dennis {dennis:.2f}")

    record_caption(out, (
        "Why Barrier3D needs dune erosion in the collision regime (Sallenger 2000), 1996–2009. "
        f"(a) Lidar profile of GIS {PROFILE_GIS} (cell {PROFILE_CELL}, median of its ten 1 m transects) before "
        f"(Sep 1998) and after (Sep 1999, before Floyd) Hurricane Dennis, whose peak total water level in the "
        f"model's storm series ({dennis:.2f} m NAVD88) reached above the dune toe but stayed below the crest. "
        "The lidar crest dropped while the model's, with no overtopping, did not (grey bar: the control's crest "
        "on 1 Jan 1999 and 1 Jan 2000). "
        f"(b) Dune crests at GIS {lo_g}–{hi_g}: the model without the collision module (line, median over the "
        "cells; band, 25th–75th percentile) and each lidar survey (dot, median; bar, 25th–75th percentile; "
        "hollow, 2001 unclassified lidar). Dotted: Bonnie, Dennis, Isabel. "
        f"(c) Share of the domains with at least {MIN_CELLS} surveyed cells whose median modelled crest stands "
        f"more than {MISMATCH_M:g} m above (red) or below (black) the lidar's, per survey. "
        "(d) Crest lost between the surveys around Dennis (1998–99) and 2004–05 by the cells the storm water "
        "reached but did not overtop (toe < peak water level < crest), against each cell's median beach "
        "width over the surveys (1998 left out, as the module's), in 10 m bins with at least "
        f"{D_MIN_CELLS} cells (line, median; band, 25th–75th percentile for Dennis). Dotted, the 45 m beach "
        "width beyond which the module's beach factor is zero. "
        "Model: the calibration control (Hs 3.0, StormStart 1, groin 0.16, collision module off)."))
    print(f"wrote {out}")


def two_dunes():
    """One storm, two dunes alike in crest and toe: Dennis cuts the one behind a narrow beach and leaves the one
    behind a wide beach, which is what the beach factor encodes. Profiles read once into TWO_DUNES_CSV."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from site_layer.hat_figure_style import INK, INK_MUTED, apply_style, record_caption
    apply_style()

    if not TWO_DUNES_CSV.exists():
        t = cell_profiles(list(TWO_DUNES.values()))
        SUPPORT_DIR.mkdir(parents=True, exist_ok=True)
        t.round(3).to_csv(TWO_DUNES_CSV, index=False)
        print(f"wrote {TWO_DUNES_CSV}")
    prof = pd.read_csv(TWO_DUNES_CSV)
    dennis = float(pd.read_csv(COLLISION_DIR / "collision_calibration.csv").set_index("interval")
                   .twl_max_m_navd88["1998-1999"])
    cells = pd.read_csv(COLLISION_DIR / "collision_calibration_cells.csv")
    cells = cells[cells.interval == "1998-1999"].set_index(["gis", "cell"])

    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.3), sharey=True)
    rows = []
    inputs = pd.read_csv(TOE_CELLS).set_index(["gis", "cell"])
    for ax, (letter, (key, (g, c))) in zip(axs, enumerate(TWO_DUNES.items())):
        w = float(inputs.loc[(g, c), "beach_width_m"])
        z98, z99 = prof[f"z_{g}_{c}_1998"].to_numpy(), prof[f"z_{g}_{c}_1999"].to_numpy()
        # the foredune: the highest point within 40 m landward of where the profile first rises above the storm
        # water (the profile's own maximum can be higher ground behind it)
        first = int(np.flatnonzero(np.nan_to_num(z98, nan=-9) >= dennis)[0])
        c98 = first + int(np.nanargmax(z98[first:first + 40]))
        x = prof.x_m.to_numpy() - c98
        keep = (x >= TWO_DUNES_X[0]) & (x <= TWO_DUNES_X[1])
        near = keep & (np.abs(x) <= 15)
        crest98, crest99 = np.nanmax(z98[near]), np.nanmax(z99[near])
        toe = float(inputs.loc[(g, c), "toe_m_navd88"])
        factor = max(0.0, 1 - w / W0_M)
        # the storm water over the beach, up to where it meets the 1998 dune face
        face = np.flatnonzero(keep & (x < 0) & (z98 >= dennis))
        x_reach = x[face[0]] if face.size else 0
        wet = keep & (x <= x_reach)
        ax.fill_between(x[wet], z98[wet], dennis, where=z98[wet] < dennis, color=WATER_C, alpha=0.12, lw=0)
        ax.plot([TWO_DUNES_X[0], x_reach], [dennis, dennis], color=WATER_C, lw=1.6)
        ax.plot(x[keep], z98[keep], color=LIDAR_C, lw=1.4, ls=(0, (4, 2)))
        ax.plot(x[keep], z99[keep], color=LIDAR_C, lw=2.0)
        ax.axhline(toe, color=INK_MUTED, lw=0.7, ls=":")
        ax.text(TWO_DUNES_X[0] + 3, toe + 0.08, f"dune toe {toe:.1f} m", color=INK_MUTED, fontsize=7, va="bottom")
        # the beach in front of the dune, as the module measures it: the cell's median over the surveys
        # the dry beach ends at the toe: where the 1998 face, walked seaward from the crest, drops to the toe
        below = np.flatnonzero(keep & (x < 0) & (z98 <= toe))
        x_toe = x[below[-1]] if below.size else 0
        yb = -0.55
        if w > 0:
            ax.annotate("", xy=(x_toe - w, yb), xytext=(x_toe, yb), arrowprops=dict(arrowstyle="<->", color=INK, lw=0.9))
        ax.text(x_toe - w / 2, yb - 0.12, f"dry beach {w:.0f} m", ha="center", va="top", fontsize=8, color=INK)
        # the crest change
        lost = crest98 - crest99
        if lost > 0.3:
            ax.annotate("", xy=(4, crest99 + 0.08), xytext=(4, crest98 - 0.08),
                        arrowprops=dict(arrowstyle="->", color=INK, lw=0.9))
            ax.text(7, (crest98 + crest99) / 2, f"crest −{lost:.1f} m", fontsize=8, color=INK, va="center")
        else:
            ax.text(7, crest98 + 0.15, "crest unchanged", fontsize=8, color=INK, va="bottom")
        ax.text(0.03, 0.97, f"beach factor = max(0, 1 − {w:.0f}/{W0_M:.0f}) = {factor:.2f}", transform=ax.transAxes,
                ha="left", va="top", fontsize=8.5, color=INK,
                bbox=dict(facecolor="white", edgecolor="0.8", boxstyle="round,pad=0.3"))
        ax.set_title(f"({'ab'[letter]})", loc="left", fontweight="bold")
        ax.set_title(f"{key.capitalize()} beach: GIS {g}", loc="center")
        ax.set_xlim(*TWO_DUNES_X)
        ax.set_xlabel("Distance from the 1998 crest (m, ocean left)")
        rows.append(dict(dune=key, gis=g, cell=c, beach_width_m=w, beach_factor=round(factor, 3),
                         toe_m_navd88=toe, crest_1998_m=round(crest98, 2), crest_1999_m=round(crest99, 2),
                         crest_change_m=round(crest99 - crest98, 2),
                         cell_crest_loss_m=float(cells.loc[(g, c), "crest_loss_m"])))
    axs[0].set_ylabel("Elevation (m NAVD88)")
    axs[0].set_ylim(-1.3, 8.0)
    fig.legend(handles=[Line2D([], [], color=LIDAR_C, lw=1.4, ls=(0, (4, 2)), label="Lidar, Sep 1998 (before Dennis)"),
                        Line2D([], [], color=LIDAR_C, lw=2.0, label="Lidar, Sep 1999 (after Dennis)"),
                        Line2D([], [], color=WATER_C, lw=1.6, label=f"Dennis peak water level, {dennis:.2f} m")],
               loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    SUPPORT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "two_dunes.png"
    fig.savefig(out, dpi=300, facecolor="white")
    fig.savefig(SUPPORT_DIR / "two_dunes.pdf", facecolor="white")
    plt.close(fig)
    pd.DataFrame(rows).to_csv(SUPPORT_DIR / "two_dunes.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    n, wd = TWO_DUNES["narrow"], TWO_DUNES["wide"]
    record_caption(out, (
        "One storm, two dunes. Lidar profiles (median of each 10 m cell's ten 1 m transects) before (Sep 1998) and "
        f"after (Sep 1999, before Floyd) Hurricane Dennis at (a) GIS {n[0]} cell {n[1]} and (b) GIS {wd[0]} cell {wd[1]}. "
        "The two dunes stood about as high and had about the same toe, and Dennis's peak water level "
        f"({dennis:.2f} m NAVD88, the model's storm series) reached above both toes and below both crests: both were "
        "in the collision regime (Sallenger 2000). The dune behind the narrow beach was cut back and lowered; the "
        "one behind the wide beach did not change. Beach width as the collision module uses it, each cell's median "
        "over the surveys (1998 left out), the dry beach from the beach start (0.5 m above MHW) to the toe; the beach factor, "
        f"max(0, 1 − W/{W0_M:.0f} m), is the share of Larson et al.'s (2004) erosion each dune receives in the model."))
    print(f"wrote {out}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("action", choices=("profiles", "plot", "two_dunes"))
    a = p.parse_args()
    {"profiles": profiles, "plot": plot, "two_dunes": two_dunes}[a.action]()


if __name__ == "__main__":
    main()
