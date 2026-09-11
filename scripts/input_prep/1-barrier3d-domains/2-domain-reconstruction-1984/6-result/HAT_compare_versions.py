#!/usr/bin/env python3
r"""
HAT_compare_versions.py
==============================================================================
v2 against v3 under the same hindcast, side by side: what the 1984
reconstruction changes in what the model DOES. One figure and one table per
run pair, read from the two runs' saved state (the .npz) and their metadata;
nothing is re-run and nothing is re-scored.

THE PAIRS (Hannah's advisor, 2026-09-09: "a comparison between v2 and v3
under the modules' automatic behaviour, full management and calibrated")
    emergent    HAT_1984_2004_calibBE_road_bdm_groin
                full management, calibBE, groin on, the roadway and
                beach-dune modules acting on their own (no prescribed
                relocations). Both versions in output/raw_runs/version-pair/<v>/,
                run by HAT_run_version_pair.py on the same code the same day.
                (The earlier pair - v2 in the calibration tree of 2026-09-07,
                v3 in behindroad-copy of 09-08 - sat on different commits, with
                the pipeline and the live setback CSV changed between them, so
                it was not a clean pair and was re-run.)
    prescribed  HAT_1984_2004_calibBE_road_reloc_bdm_groin
                the same with the recorded 1989 (GIS 84-87) and 1999
                (GIS 9-14) relocations imposed: the control, showing what
                the v3 setbacks change when the module is not deciding.
                Both versions in output/raw_runs/version-pair/<v>/, run by
                HAT_run_version_pair.py --relocations 1.

THE FIGURES, per pair (two double-column figures since 2026-09-10, in the
house style of hat_figure_style: v2 grey C["BASE"], v3 purple C["ACCENT"],
the recorded events C["REF"]; no title sentence on the canvas, the run name
and the pair go to CAPTIONS.md beside the PNGs)
  HAT_compare_v2_v3_<pair>_relocations.png
    (a) the NC-12 setback the model starts with, per road domain, two thin
        lines with markers on a symlog axis
    (b) every year the model relocated NC-12, per domain, the recorded
        events outlined; counts in the legend. In the prescribed pair the
        1989/1999 rows are inputs, so only the OTHER relocations are the
        module's own.
    (c) the number of relocations per year, island-wide, through time
  HAT_compare_v2_v3_<pair>_geometry.png
    (a) interior width per domain at 1984 (dashed) and 2004 (solid), both
        versions: what the footprint added or removed, and what the run
        then did with it
    (b) island-mean interior width through time, (c) the difference
        v3 - v2 on its own axis, (d) island-total cumulative overwash
    (e) mean interior elevation of the land cells at 2004 per domain, both
        versions, and (f) the difference v3 - v2 on its own axis
    The island-wide shoreline skill of both runs (they are near-identical by
    construction: the shoreline offset does not read the topography) is in
    the table and the report, not on the figure.

THE TABLE  version_compare_<pair>.csv, one row per domain: initial setback,
    relocation years and count, drowned, interior rows and width at 1984
    and 2004, mean land elevation at 1984 and 2004, cumulative overwash -
    for v2, for v3, and the difference.
THE REPORT  HAT_compare_versions.txt: the skill of the four runs, the
    relocation counts and their timing against the recorded events (mean
    error over the event blocks; a domain that never relocated is censored
    and counted, not averaged), and the island-wide geometry medians.

UNITS. Barrier3D stores decametres: widths and elevations x 10 -> m;
QowTS is dam^3 per dam of shoreline per year, x 100 -> m^3/m. The buffer
is 15 domains: GIS g is index g + 14.

USAGE
    python HAT_compare_versions.py                 # both pairs, whatever exists
    python HAT_compare_versions.py --pairs emergent
==============================================================================
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "0-elevation" / "3-figures"))
from cascade_pipeline.run_layout import resolve as resolve_run_file  # noqa: E402
from hat_topo_version import insert_figures_dir, insert_scope_step  # noqa: E402
from hat_figure_style import (  # noqa: E402
    C as STYLE_C, INK, DOMAIN_AXIS_LABEL, apply_style, figsize, open_frame, record_caption, save, town_bands, _title,
)

RAW = REPO / "output" / "raw_runs"
PAIRS = {
    "emergent": {
        "name": "HAT_1984_2004_calibBE_road_bdm_groin",
        "v2": RAW / "version-pair" / "v2" / "1984_2004" / "calibBE" / "HAT_1984_2004_calibBE_road_bdm_groin",
        "v3": RAW / "version-pair" / "v3" / "1984_2004" / "calibBE" / "HAT_1984_2004_calibBE_road_bdm_groin",
        "what": "full management, calibBE, groin; the modules act on their own",
    },
    "prescribed": {
        "name": "HAT_1984_2004_calibBE_road_reloc_bdm_groin",
        "v2": RAW / "version-pair" / "v2" / "1984_2004" / "calibBE" / "HAT_1984_2004_calibBE_road_reloc_bdm_groin",
        "v3": RAW / "version-pair" / "v3" / "1984_2004" / "calibBE" / "HAT_1984_2004_calibBE_road_reloc_bdm_groin",
        "what": "the same with the recorded 1989 and 1999 relocations prescribed (the control)",
    },
}
BUFFER = 15
START, YEARS = 1984, 20
DAM = 10.0
EVENTS = {1989: (84, 87), 1999: (9, 14)}
VS = ("v2", "v3")
C = {"v2": STYLE_C["BASE"], "v3": STYLE_C["ACCENT"]}       # the unmodified input, the modification under test
C_REF = STYLE_C["REF"]                                     # the recorded events
MK = {"v2": "o", "v3": "s"}
LABEL = {"v2": "as extracted (1996 surface)", "v3": "1984 reconstruction"}
SHORT = {"v2": "as extracted", "v3": "reconstruction"}
FIG_DIR = insert_figures_dir("1984-start", "6-result", "island")
TAB_DIR = insert_scope_step("1984-start", "6-result")


# =============================================================================
# READING A RUN
# =============================================================================

def load_run(d: Path, name: str) -> dict:
    # RESOLVED, NOT JOINED: run_layout knows where a run folder keeps each of
    # its files, in the new layout and the old flat one alike.
    meta = json.load(open(resolve_run_file(d, "metadata_json", name)))
    c = np.load(resolve_run_file(d, "archive", name, must_exist=True),
                allow_pickle=True)["cascade"][0]
    per = {}
    for g in range(1, 91):
        i = g + BUFFER - 1
        b = c.barrier3d[i]
        rec = {"domain": g}
        w = np.asarray(b.InteriorWidth_AvgTS, float) * DAM
        rec["width_1984_m"], rec["width_2004_m"] = float(w[0]), float(w[-1])
        rec["width_TS"] = w
        q = np.asarray(b.QowTS, float) * 100.0            # dam^3/dam -> m^3/m
        rec["qow_cum_m3_m"] = float(np.nansum(q))
        rec["qow_TS"] = q
        for lab, k in (("1984", 0), ("2004", -1)):
            z = np.asarray(b.DomainTS[k], float) * DAM
            land = z[z > -2.99]
            rec[f"rows_{lab}"] = int(z.shape[0])
            rec[f"zmean_{lab}_m"] = float(land.mean()) if land.size else np.nan
        rec["dune_1984_m"] = float(np.mean(np.asarray(b.DuneDomain[0], float)[:, 0]) * DAM)
        rec["dune_2004_m"] = float(np.mean(np.asarray(b.DuneDomain[-1], float)[:, 0]) * DAM)
        mgr = c.roadways[i]
        if mgr is not None:
            sb = np.asarray(getattr(mgr, "_road_setback_TS", []), float).ravel()
            rel = np.asarray(getattr(mgr, "_road_relocated_TS", []), float).ravel()
            rec["setback_1984_m"] = float(sb[0]) if sb.size else np.nan
            rec["setback_2004_m"] = float(sb[-1]) if sb.size else np.nan
            rec["reloc_years"] = [START + int(k) for k in np.flatnonzero(np.nan_to_num(rel) > 0)]
            # a PRESCRIBED event is applied as a displacement of the setback and
            # does not raise the relocated flag: read it off the setback series
            # as the jump at the event year in the event's block
            # (displacements run 17-108 m; setbacks quantise to 10 m; an emergent
            # relocation the same year would carry the flag and is excluded)
            flagged = set(np.flatnonzero(np.nan_to_num(rel) > 0).tolist())
            rec["prescribed_years"] = [yr for yr, (lo, hi) in EVENTS.items()
                                       if lo <= g <= hi and sb.size > yr - START
                                       and sb[yr - START] - sb[yr - START - 1] >= 10.0
                                       and (yr - START) not in flagged]
            rec["setback_TS"] = sb
            rec["drowned"] = bool(getattr(mgr, "drown_break", False))
        else:
            rec["setback_1984_m"] = rec["setback_2004_m"] = np.nan
            rec["reloc_years"] = []
            rec["prescribed_years"] = []
            rec["setback_TS"] = np.array([])
            rec["drowned"] = False
        rec["n_reloc"] = len(rec["reloc_years"])
        per[g] = rec
    return {"meta": meta, "per": per, "dir": d}


def timing_score(per: dict) -> dict:
    """Mean error of the first relocation against the recorded event year, per
    block; domains that never relocate are censored, counted, not averaged.
    Meaningful only when the events are NOT prescribed (scoring a run against
    its own input is circular); the prescribed pair reports counts instead."""
    out = {}
    for yr, (lo, hi) in EVENTS.items():
        errs, cens = [], 0
        for g in range(lo, hi + 1):
            ys = per[g]["reloc_years"]
            if ys:
                errs.append(min(ys) - yr)
            else:
                cens += 1
        out[yr] = {"mean_err_yr": float(np.mean(errs)) if errs else np.nan, "n_scored": len(errs), "n_censored": cens}
    return out


def own_relocations(per: dict) -> int:
    """Relocations that are the module's own in the prescribed pair: everything
    outside the (event year, event block) pairs."""
    n = 0
    for g, r in per.items():
        for y in r["reloc_years"]:
            if not any(y == yr and lo <= g <= hi for yr, (lo, hi) in EVENTS.items()):
                n += 1
    return n


# =============================================================================
# THE FIGURES
# =============================================================================
# Two double-column figures per pair (2026-09-10, the house style): the road
# (setback, relocations, relocations per year) and the geometry (widths,
# overwash, elevation, each difference on its own axis). The run name and the
# pair description go to CAPTIONS.md beside the PNGs, not onto the canvas.

def _events(ax, label=True):
    """The recorded relocation blocks, outlined, with the event year as a bar."""
    for yr, (lo, hi) in EVENTS.items():
        ax.axvspan(lo - .5, hi + .5, facecolor="none", edgecolor=C_REF, lw=0.8, ls=(0, (3, 2)), zorder=1)
        ax.hlines(yr, lo - .5, hi + .5, color=C_REF, lw=1.6, zorder=3)
        if label:
            right = hi < 60
            ax.text(hi + 1.6 if right else lo - 1.6, yr, f"recorded {yr}", va="center",
                    ha="left" if right else "right", fontsize=7, color=C_REF, zorder=6)


def _domain_axis(ax, label=True):
    ax.set_xlim(0.2, 90.8)
    ax.set_xticks([1] + list(range(10, 91, 10)))
    if label:
        ax.set_xlabel(DOMAIN_AXIS_LABEL)
    else:
        ax.tick_params(labelbottom=False)


def fig_relocations(key: str, runs: dict, out: Path) -> Path:
    apply_style()
    per = {v: runs[v]["per"] for v in VS}
    doms = np.arange(1, 91)
    years = np.arange(START, START + YEARS + 1)

    fig = plt.figure(figsize=figsize("double", height=6.6), constrained_layout=True)
    gs = fig.add_gridspec(3, 1, height_ratios=[0.85, 1.35, 0.8])
    a = fig.add_subplot(gs[0])
    b = fig.add_subplot(gs[1], sharex=a)
    c = fig.add_subplot(gs[2])

    # (a) the setback the model starts with: two thin lines with markers
    town_bands(a)                      # names here, at the top, where nothing is drawn
    for v in VS:
        y = np.array([per[v][g]["setback_1984_m"] for g in doms], float)
        a.plot(doms, y, ls="-", lw=0.9, marker=MK[v], ms=2.6, mew=0, color=C[v], zorder=3 + (v == "v3"),
               label=LABEL[v])
    a.axhline(0, color=INK, lw=0.5, zorder=2)
    a.set_yscale("symlog", linthresh=20, linscale=0.8)
    a.set_yticks([0, 10, 20, 50, 100, 200, 500])
    a.set_yticklabels(["0", "10", "20", "50", "100", "200", "500"])
    a.set_ylim(-2, 1500)
    a.set_ylabel("NC-12 setback at 1984\n(m landward of row 0)")
    a.grid(axis="y")
    a.set_axisbelow(True)
    open_frame(a)
    a.legend(loc="lower center", ncol=2, frameon=False, handlelength=2.2)
    _domain_axis(a, label=False)
    _title(a, 0, "the setback the model starts with")

    # (b) every relocation, by year
    town_bands(b, label=False)
    _events(b)
    handles = []
    for v, dx in (("v2", -0.18), ("v3", 0.18)):
        xs, ys = [], []
        for g in doms:
            for y in per[v][g]["reloc_years"]:
                xs.append(g + dx); ys.append(y)
        b.plot(xs, ys, MK[v], ms=4.0, color=C[v], mec="white", mew=0.5, ls="none", zorder=4)
        handles.append(Line2D([0], [0], marker=MK[v], ms=4.0, color=C[v], mec="white", mew=0.5, ls="none",
                              label=f"{LABEL[v]}: {len(xs)} relocations"))
        # a prescribed event is an input: a tick in the run's colour across the
        # recorded bar where it was applied, one per version, side by side
        px = [(g + 1.5 * dx, y) for g in doms for y in per[v][g]["prescribed_years"]]
        if px:
            b.plot([x for x, _ in px], [y for _, y in px], "|", ms=7, color=C[v], mew=1.2, ls="none", zorder=5)
        dr = [g for g in doms if per[v][g]["drowned"]]
        if dr:
            b.plot([g + 1.8 * dx for g in dr], [START + YEARS + 0.8] * len(dr), "x", color=C[v], ms=4.5, mew=1.3,
                   ls="none", zorder=5)
    npre = {v: sum(len(per[v][g]["prescribed_years"]) for g in doms) for v in VS}
    if any(npre.values()):
        handles.append(Line2D([0], [0], marker="|", ms=7, color=INK, mew=1.2, ls="none",
                              label=f"recorded event applied as prescribed ({npre['v2']}, {npre['v3']})"))
    ndr = {v: sum(1 for g in doms if per[v][g]["drowned"]) for v in VS}
    if any(ndr.values()):
        handles.append(Line2D([0], [0], marker="x", ms=4.5, mew=1.3, color=INK, ls="none",
                              label=f"road drowned ({ndr['v2']}, {ndr['v3']})"))
    handles.append(Line2D([0], [0], color=C_REF, lw=1.6, label="recorded relocation, its block outlined"))
    b.set_ylim(START - 0.6, START + YEARS + 1.6)
    b.set_yticks(range(START, START + YEARS + 1, 4))
    b.set_ylabel("year of a relocation")
    b.grid(axis="y")
    b.set_axisbelow(True)
    open_frame(b)
    b.legend(handles=handles, loc="upper center", ncol=2, frameon=False, handletextpad=0.5)
    _domain_axis(b)
    _title(b, 1, "when the model relocated NC-12")

    # (c) relocations per year, island-wide
    for v, dx in (("v2", -0.2), ("v3", 0.2)):
        cnt = np.zeros(len(years))
        for g in doms:
            for y in per[v][g]["reloc_years"]:
                cnt[y - START] += 1
        c.step(years, np.cumsum(cnt), where="post", color=C[v], lw=1.4, label=f"{LABEL[v]}: {int(cnt.sum())} in all")
        c.bar(years + dx, cnt, 0.4, color=C[v], alpha=0.45, lw=0)
    c.set_xlim(START - 0.6, START + YEARS + 0.6)
    c.set_xticks(range(START, START + YEARS + 1, 2))
    c.set_ylabel("relocations")
    c.set_xlabel("year")
    c.grid(axis="y")
    c.set_axisbelow(True)
    open_frame(c)
    c.legend(loc="upper left", frameon=False)
    _title(c, 2, "relocations per year (bars) and cumulative (steps)")

    save(fig, out)
    plt.close(fig)
    return out.with_suffix(".png")


def fig_geometry(key: str, runs: dict, out: Path) -> Path:
    apply_style()
    per = {v: runs[v]["per"] for v in VS}
    doms = np.arange(1, 91)
    years = np.arange(START, START + YEARS + 1)
    yticks = range(START, START + YEARS + 1, 4)

    fig = plt.figure(figsize=figsize("double", height=8.6), constrained_layout=True)
    gs = fig.add_gridspec(4, 3, height_ratios=[1.15, 1.0, 1.0, 0.7])
    a = fig.add_subplot(gs[0, :])
    b = fig.add_subplot(gs[1, 0])
    c = fig.add_subplot(gs[1, 1], sharex=b)
    d = fig.add_subplot(gs[1, 2], sharex=b)
    e = fig.add_subplot(gs[2, :], sharex=a)
    f = fig.add_subplot(gs[3, :], sharex=a)

    # (a) interior width at 1984 and 2004
    town_bands(a, where="bottom")
    for v in VS:
        a.plot(doms, [per[v][g]["width_1984_m"] for g in doms], ls=(0, (3, 2)), lw=1.0, color=C[v])
        a.plot(doms, [per[v][g]["width_2004_m"] for g in doms], "-", lw=1.3, color=C[v])
    handles = [Line2D([0], [0], color=C[v], lw=1.3, label=LABEL[v]) for v in VS]
    handles += [Line2D([0], [0], color=INK, lw=1.0, ls=(0, (3, 2)), label="at 1984 (start)"),
                Line2D([0], [0], color=INK, lw=1.3, label="at 2004 (end)")]
    a.set_ylim(0, None)
    a.set_ylabel("interior width (m)")
    a.grid(axis="y")
    a.set_axisbelow(True)
    open_frame(a)
    a.legend(handles=handles, loc="upper right", ncol=2, frameon=False)
    _domain_axis(a, label=False)
    _title(a, 0, "interior width per domain at the start and the end of the run")

    # (b, c, d) through time, island-wide
    W = {v: np.array([per[v][g]["width_TS"] for g in doms]) for v in VS}
    Q = {v: np.array([per[v][g]["qow_TS"] for g in doms]) for v in VS}
    wmean = {v: np.nanmean(W[v], axis=0) for v in VS}
    qcum = {v: np.cumsum(np.nansum(Q[v], axis=0)) / 1e3 for v in VS}
    for v in VS:
        b.plot(years, wmean[v], color=C[v], lw=1.4)
        d.plot(years, qcum[v], color=C[v], lw=1.4)
    c.plot(years, wmean["v3"] - wmean["v2"], color=C["v3"], lw=1.4)
    c.axhline(0, color=INK, lw=0.5)
    b.set_ylabel("island-mean interior width (m)")
    c.set_ylabel(f"{SHORT['v3']} − {SHORT['v2']} (m)")
    d.set_ylabel("island total (10³ m³/m)")
    for ax in (b, c, d):
        ax.set_xticks(yticks)
        ax.set_xlabel("year")
        ax.grid(axis="y")
        ax.set_axisbelow(True)
        open_frame(ax)
    _title(b, 1, "island-mean width")
    _title(c, 2, "difference in (b)")
    _title(d, 3, "cumulative overwash")

    # (e) mean land elevation at 2004, (f) the difference on its own axis
    town_bands(e, label=False)
    for v in VS:
        e.plot(doms, [per[v][g]["zmean_2004_m"] for g in doms], "-", lw=1.2, color=C[v], label=LABEL[v])
    e.set_ylabel("mean elevation of\nland cells (m MHW)")
    e.grid(axis="y")
    e.set_axisbelow(True)
    open_frame(e)
    e.legend(loc="upper right", ncol=2, frameon=False)
    _domain_axis(e, label=False)
    _title(e, 4, "mean interior elevation at 2004")

    town_bands(f, label=False)
    dz = np.array([per["v3"][g]["zmean_2004_m"] - per["v2"][g]["zmean_2004_m"] for g in doms])
    f.bar(doms, dz, 0.7, color=C["v3"], lw=0, zorder=3)
    f.axhline(0, color=INK, lw=0.5)
    lim = 1.15 * np.nanmax(np.abs(dz)) if np.isfinite(dz).any() and np.nanmax(np.abs(dz)) > 0 else 0.1
    f.set_ylim(-lim, lim)
    f.set_ylabel("difference (m)")
    f.grid(axis="y")
    f.set_axisbelow(True)
    open_frame(f)
    _domain_axis(f)
    _title(f, 5, f"elevation at 2004, {SHORT['v3']} − {SHORT['v2']}")

    save(fig, out)
    plt.close(fig)
    return out.with_suffix(".png")


def write_captions(key: str, runs: dict, rel: Path, geo: Path) -> None:
    info = PAIRS[key]
    per = {v: runs[v]["per"] for v in VS}
    n = {v: sum(per[v][g]["n_reloc"] for g in per[v]) for v in VS}
    nd = {v: sum(1 for g in per[v] if per[v][g]["n_reloc"]) for v in VS}
    dr = {v: sum(1 for g in per[v] if per[v][g]["drowned"]) for v in VS}
    npre = {v: sum(len(per[v][g]["prescribed_years"]) for g in per[v]) for v in VS}
    sk = {v: runs[v]["meta"]["skill"] for v in VS}
    head = (f"Run `{info['name']}` ({key}: {info['what']}), the 1984–2004 hindcast on the domains as extracted "
            f"from the 1996 surface (grey; dune-topo v2) and on the 1984 reconstruction (purple; v3, rows behind "
            f"NC-12 filled by copy, 1984 setbacks), both run by `HAT_run_version_pair.py` on the same code. ")
    if key == "prescribed":
        rel_text = (f"(b) Every year the model relocated NC-12, per domain; the recorded 1989 (GIS 84–87) and 1999 "
                    f"(GIS 9–14) relocations are green bars across their outlined blocks and, prescribed here, are "
                    f"applied as displacements of the setback (a tick across the bar in the run's colour where applied: "
                    f"{npre['v2']} and {npre['v3']}); "
                    f"the filled markers are the module's own, {n['v2']} in {nd['v2']} domains as extracted and "
                    f"{n['v3']} in {nd['v3']} on the reconstruction; a cross at the top marks a drowned road "
                    f"({dr['v2']} and {dr['v3']}). ")
    else:
        rel_text = (f"(b) Every year the model relocated NC-12, per domain, every one the module's own: {n['v2']} "
                    f"in {nd['v2']} domains as extracted, {n['v3']} in {nd['v3']} on the reconstruction; the recorded "
                    f"1989 (GIS 84–87) and 1999 (GIS 9–14) relocations are green bars across their outlined blocks; "
                    f"a cross at the top would mark a drowned road ({dr['v2']} and {dr['v3']}). ")
    record_caption(rel, head +
                   "(a) The NC-12 setback the model starts with, per domain, on a symmetric-log axis (linear below "
                   "20 m): measured on the 1996 surface and floored at 0, or the 1984 measurement, unfloored. " +
                   rel_text +
                   "(c) Relocations per year (bars) and cumulative (steps), island-wide. Domain 1 is Cape Point, 90 "
                   "Pea Island; villages banded. Interior shoreline RMSE "
                   f"{sk['v2']['rmse_interior_m_yr']} (as extracted) and {sk['v3']['rmse_interior_m_yr']} "
                   "(reconstruction) m/yr, near-identical by construction: the shoreline offset does not read the "
                   "topography. Skill, counts and timing in `6-result/HAT_compare_versions.txt`; the per-domain "
                   f"table in `version_compare_{key}.csv`.")
    doms = range(1, 91)
    w84 = np.array([per["v3"][g]["width_1984_m"] - per["v2"][g]["width_1984_m"] for g in doms])
    dz = np.array([per["v3"][g]["zmean_2004_m"] - per["v2"][g]["zmean_2004_m"] for g in doms])
    record_caption(geo, head +
                   "(a) Interior width per domain at 1984 (dashed: what the reconstruction added or removed, "
                   f"{int((w84 > 0).sum())} domains wider and {int((w84 < 0).sum())} narrower) and at 2004 (solid: "
                   "what the run left). (b) Island-mean interior width through time and (c) its difference, "
                   "reconstruction minus as extracted, on its own axis. (d) Cumulative overwash flux summed over the "
                   "90 domains. (e) Mean elevation of the land cells (above −3 m) at 2004 per domain and (f) the "
                   f"difference, reconstruction minus as extracted (largest |difference| {np.nanmax(np.abs(dz)):.2f} m "
                   f"at GIS {int(np.nanargmax(np.abs(dz))) + 1}). Domain 1 is Cape Point, 90 Pea Island; villages "
                   "banded. Island medians in `6-result/HAT_compare_versions.txt`; the per-domain table in "
                   f"`version_compare_{key}.csv`.")


# =============================================================================
# TABLE AND REPORT
# =============================================================================

def table(key: str, runs: dict) -> Path:
    rows = []
    for g in range(1, 91):
        r = {"domain": g}
        for v in ("v2", "v3"):
            p = runs[v]["per"][g]
            for k in ("setback_1984_m", "setback_2004_m", "n_reloc", "drowned", "rows_1984", "rows_2004",
                      "width_1984_m", "width_2004_m", "zmean_1984_m", "zmean_2004_m", "dune_1984_m", "dune_2004_m",
                      "qow_cum_m3_m"):
                r[f"{k}_{v}"] = p[k]
            r[f"reloc_years_{v}"] = " ".join(str(y) for y in p["reloc_years"])
        for k in ("setback_1984_m", "n_reloc", "rows_1984", "width_1984_m", "width_2004_m", "zmean_2004_m", "qow_cum_m3_m"):
            r[f"d_{k}"] = (r[f"{k}_v3"] - r[f"{k}_v2"]) if np.isfinite(_f(r[f"{k}_v3"])) and np.isfinite(_f(r[f"{k}_v2"])) else np.nan
        rows.append(r)
    p = TAB_DIR / f"version_compare_{key}.csv"
    pd.DataFrame(rows).round(2).to_csv(p, index=False)
    return p


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def report(results: dict, figs: list[Path], tabs: list[Path]) -> Path:
    L = [f"HAT_compare_versions.txt - v2 against v3 under the same hindcast ({datetime.now():%Y-%m-%d %H:%M})", ""]
    for key, runs in results.items():
        info = PAIRS[key]
        L.append(f"PAIR {key}: {info['name']}")
        L.append(f"  {info['what']}")
        for v in ("v2", "v3"):
            m = runs[v]["meta"]
            L.append(f"  {v}: {runs[v]['dir'].relative_to(REPO)}")
            L.append(f"      run {m['identity']['timestamp']}  topo {m['identity']['topo_dune_version']}  "
                     f"commit {m['identity']['git_commit'][:8]}  relocations_enabled {m['dunes and roadway']['relocations_enabled']}")
            L.append(f"      skill: interior RMSE {m['skill']['rmse_interior_m_yr']} m/yr, bias {m['skill']['mean_bias_interior_m_yr']}; "
                     f"roads drowned {m['verification']['roads_drowned']}")
        per = {v: runs[v]["per"] for v in ("v2", "v3")}
        for v in ("v2", "v3"):
            n = sum(per[v][g]["n_reloc"] for g in per[v])
            nd = sum(1 for g in per[v] if per[v][g]["n_reloc"])
            if key == "prescribed":
                npre = sum(len(per[v][g]["prescribed_years"]) for g in per[v])
                L.append(f"  {v}: {npre} prescribed relocations applied (the 1989 and 1999 events, read off the setback "
                         f"jumps; no timing score) and {n} of the module's own in {nd} domains")
                continue
            L.append(f"  {v}: {n} relocations in {nd} domains")
            for yr, s_ in timing_score(per[v]).items():
                L.append(f"      first relocation vs the recorded {yr} event (GIS {EVENTS[yr][0]}-{EVENTS[yr][1]}): "
                         f"mean error {s_['mean_err_yr']:+.1f} yr over {s_['n_scored']} domains, {s_['n_censored']} never relocated")
        doms = range(1, 91)
        for lab, k in (("interior width 1984 (m)", "width_1984_m"), ("interior width 2004 (m)", "width_2004_m"),
                       ("mean land elevation 2004 (m MHW)", "zmean_2004_m"), ("cumulative overwash (m3/m)", "qow_cum_m3_m")):
            a2 = np.array([per["v2"][g][k] for g in doms]); a3 = np.array([per["v3"][g][k] for g in doms])
            L.append(f"  {lab}: island median v2 {np.nanmedian(a2):.1f}, v3 {np.nanmedian(a3):.1f}; "
                     f"median difference {np.nanmedian(a3 - a2):+.1f}, max |diff| {np.nanmax(np.abs(a3 - a2)):.1f} at GIS "
                     f"{int(np.nanargmax(np.abs(a3 - a2))) + 1}")
        L.append("")
    L.append("figures: " + ", ".join(str(p.relative_to(REPO)) for p in figs))
    L.append("tables:  " + ", ".join(str(p.relative_to(REPO)) for p in tabs))
    p = TAB_DIR / "HAT_compare_versions.txt"
    p.write_text("\n".join(L) + "\n", encoding="utf-8")
    return p


def main() -> None:
    ap = argparse.ArgumentParser(description="v2 against v3 under the same hindcast")
    ap.add_argument("--pairs", default=",".join(PAIRS))
    a = ap.parse_args()
    results, figs, tabs = {}, [], []
    for key in [k.strip() for k in a.pairs.split(",") if k.strip()]:
        info = PAIRS[key]
        missing = [v for v in ("v2", "v3")
                   if not resolve_run_file(info[v], "metadata_json",
                                           info["name"]).is_file()]
        if missing:
            print(f"{key}: no run for {missing} at {[str(info[v].relative_to(REPO)) for v in missing]} - skipped")
            continue
        print(f"{key}: loading both runs ...")
        runs = {v: load_run(info[v], info["name"]) for v in ("v2", "v3")}
        results[key] = runs
        rel = fig_relocations(key, runs, FIG_DIR / f"HAT_compare_v2_v3_{key}_relocations.png")
        geo = fig_geometry(key, runs, FIG_DIR / f"HAT_compare_v2_v3_{key}_geometry.png")
        write_captions(key, runs, rel, geo)
        figs += [rel, geo]
        tabs.append(table(key, runs))
        print(f"  wrote {rel.name}, {geo.name}, {tabs[-1].name}")
    if results:
        rep = report(results, figs, tabs)
        print(f"wrote {rep}")


if __name__ == "__main__":
    main()
