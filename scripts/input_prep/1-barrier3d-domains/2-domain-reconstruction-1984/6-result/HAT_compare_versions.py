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

THE FIGURE, per pair (rows; v2 grey, v3 red on every panel)
    (a) the NC-12 setback the model starts with, per road domain
    (b) every year the model relocated NC-12, per domain, the recorded
        events outlined; counts in the legend. In the prescribed pair the
        1989/1999 rows are inputs, so only the OTHER relocations are the
        module's own.
    (c) the number of relocations per year, island-wide, through time
    (d) interior width per domain at 1984 (dashed) and 2004 (solid), both
        versions: what the footprint added or removed, and what the run
        then did with it
    (e) island-mean interior width and island-total cumulative overwash
        through time, both versions
    (f) mean interior elevation of the land cells at 2004 per domain, both
        versions, and the difference v3 - v2
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
from matplotlib.patches import Patch


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "0-elevation" / "3-figures"))
from hat_topo_version import insert_figures_dir, insert_scope_step  # noqa: E402
import HAT_plot_duneline_offset as off  # noqa: E402

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
C = {"v2": "0.45", "v3": off.C_1984}
MK = {"v2": "o", "v3": "s"}
INK = off.INK
FIG_DIR = insert_figures_dir("1984-start", "6-result", "island")
TAB_DIR = insert_scope_step("1984-start", "6-result")


# =============================================================================
# READING A RUN
# =============================================================================

def load_run(d: Path, name: str) -> dict:
    meta = json.load(open(d / f"{name}_run_metadata.json"))
    c = np.load(next(d.glob("*.npz")), allow_pickle=True)["cascade"][0]
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
# THE FIGURE
# =============================================================================

def _towns(ax, label=False):
    ann = off.HATTERAS_ANNOTATIONS
    for name, (lo, hi) in ann.town_spans.items():
        ax.axvspan(lo - .5, hi + .5, color="0.93", zorder=0)
        if label:
            ax.text((lo + hi) / 2, 0.985, name, transform=ax.get_xaxis_transform(), ha="center", va="top",
                    fontsize=7, color=off.INK_MUTED)


def fig_pair(key: str, runs: dict, out: Path) -> Path:
    off.apply_style()
    info = PAIRS[key]
    per = {v: runs[v]["per"] for v in ("v2", "v3")}
    doms = np.arange(1, 91)
    road = [g for g in doms if np.isfinite(per["v2"][g]["setback_1984_m"]) or np.isfinite(per["v3"][g]["setback_1984_m"])]
    years = np.arange(START, START + YEARS + 1)

    fig = plt.figure(figsize=(14.5, 17.0), constrained_layout=True)
    gs = fig.add_gridspec(6, 2, height_ratios=[0.8, 1.2, 0.7, 1.0, 0.9, 1.0])
    a = fig.add_subplot(gs[0, :])
    b = fig.add_subplot(gs[1, :], sharex=a)
    c = fig.add_subplot(gs[2, :])
    d = fig.add_subplot(gs[3, :], sharex=a)
    e1 = fig.add_subplot(gs[4, 0])
    e2 = fig.add_subplot(gs[4, 1])
    f = fig.add_subplot(gs[5, :], sharex=a)

    # (a) the setback the model starts with
    _towns(a, label=True)
    for v, dx in (("v2", -0.2), ("v3", 0.2)):
        a.bar([g + dx for g in road], [per[v][g]["setback_1984_m"] for g in road], 0.4, color=C[v], zorder=3,
              label=f"{v}: setback at 1984")
    a.axhline(0, color=INK, lw=0.6)
    a.set_yscale("symlog", linthresh=50, linscale=1.0)
    a.set_yticks([0, 10, 20, 50, 100, 200, 500])
    a.set_yticklabels(["0", "10", "20", "50", "100", "200", "500"])
    a.set_ylabel("NC-12 setback at 1984\n(m landward of row 0)")
    a.legend(loc="upper center", ncol=2, fontsize=8)
    off._title(a, 0, "the setback the model starts with: as measured on the 1996 surface (v2) and the 1984 one (v3)")

    # (b) every relocation, by year
    _towns(b)
    for yr, (lo, hi) in EVENTS.items():
        b.axvspan(lo - .5, hi + .5, facecolor="none", edgecolor="#2c6e49", lw=0.9, ls=(0, (3, 2)), zorder=1)
        b.hlines(yr, lo - .5, hi + .5, color="#2c6e49", lw=1.6, zorder=3)
        b.text(hi + 0.7, yr, f"recorded {yr}", va="center", fontsize=7, color="#2c6e49")
    for v, dx in (("v2", -0.18), ("v3", 0.18)):
        xs, ys = [], []
        for g in doms:
            for y in per[v][g]["reloc_years"]:
                xs.append(g + dx); ys.append(y)
        b.plot(xs, ys, MK[v], ms=4.2, color=C[v], mec="white", mew=0.5, ls="none", zorder=4,
               label=f"{v}: {len(xs)} relocations")
        px = [(g + dx, y) for g in doms for y in per[v][g]["prescribed_years"]]
        if px:
            b.plot([x for x, _ in px], [y for _, y in px], MK[v], ms=6.5, mfc="none", mec="#2c6e49", mew=1.4,
                   ls="none", zorder=5, label=f"{v}: prescribed event applied ({len(px)})")
        dr = [g for g in doms if per[v][g]["drowned"]]
        if dr:
            b.plot([g + dx for g in dr], [START + YEARS + 0.8] * len(dr), "x", color=C[v], ms=6, mew=1.5, ls="none",
                   zorder=5, label=f"{v}: road drowned ({len(dr)})")
    b.set_ylim(START - 0.5, START + YEARS + 1.5)
    b.set_yticks(range(START, START + YEARS + 1, 4))
    b.set_ylabel("model year of a relocation")
    b.legend(loc="upper left", ncol=4, fontsize=7.5)
    off._title(b, 1, "when the model relocated NC-12" + (" (1989 and 1999 are prescribed here; the rest is the module's own)"
                                                          if key == "prescribed" else " (every one the module's own)"))

    # (c) relocations per year, island-wide
    for v in ("v2", "v3"):
        cnt = np.zeros(len(years))
        for g in doms:
            for y in per[v][g]["reloc_years"]:
                cnt[y - START] += 1
        c.step(years, np.cumsum(cnt), where="post", color=C[v], lw=1.6, label=f"{v}: cumulative ({int(cnt.sum())})")
        c.bar(years + (-0.2 if v == "v2" else 0.2), cnt, 0.4, color=C[v], alpha=0.5)
    c.set_xlim(START - 0.5, START + YEARS + 0.5)
    c.set_xticks(range(START, START + YEARS + 1, 2))
    c.set_ylabel("relocations")
    c.set_xlabel("year")
    c.legend(loc="upper left", fontsize=8)
    off._title(c, 2, "relocations per year (bars) and cumulative (steps), island-wide")

    # (d) interior width at 1984 and 2004
    _towns(d)
    for v in ("v2", "v3"):
        d.plot(doms, [per[v][g]["width_1984_m"] for g in doms], ls=(0, (3, 2)), lw=1.1, color=C[v], label=f"{v} at 1984")
        d.plot(doms, [per[v][g]["width_2004_m"] for g in doms], "-", lw=1.5, color=C[v], label=f"{v} at 2004")
    d.set_ylabel("interior width (m)")
    d.legend(loc="upper right", ncol=4, fontsize=8)
    off._title(d, 3, "interior width per domain: what the footprint changed at 1984 (dashed) and what the run left at 2004 (solid)")

    # (e) through time, island-wide
    for v in ("v2", "v3"):
        W = np.array([per[v][g]["width_TS"] for g in doms])
        e1.plot(years, np.nanmean(W, axis=0), color=C[v], lw=1.6, label=f"{v}")
        Q = np.array([per[v][g]["qow_TS"] for g in doms])
        e2.plot(years, np.cumsum(np.nansum(Q, axis=0)), color=C[v], lw=1.6, label=f"{v}")
    e1.set_ylabel("island-mean interior width (m)")
    e1.set_xlabel("year")
    e1.set_xticks(range(START, START + YEARS + 1, 4))
    e2.set_xticks(range(START, START + YEARS + 1, 4))
    e1.legend(loc="upper right", fontsize=8)
    off._title(e1, 4, "island-mean interior width through time")
    e2.set_ylabel("cumulative overwash, island total (m³/m)")
    e2.set_xlabel("year")
    e2.legend(loc="upper left", fontsize=8)
    off._title(e2, 5, "overwash through time, summed over the 90 domains")

    # (f) mean land elevation at 2004 and the difference
    _towns(f)
    for v in ("v2", "v3"):
        f.plot(doms, [per[v][g]["zmean_2004_m"] for g in doms], "-", lw=1.4, color=C[v], label=f"{v} at 2004")
    dz = np.array([per["v3"][g]["zmean_2004_m"] - per["v2"][g]["zmean_2004_m"] for g in doms])
    f2 = f.twinx()
    f2.bar(doms, dz, 0.6, color=off.C_1984, alpha=0.25, zorder=1, label="v3 − v2")
    f2.axhline(0, color=INK, lw=0.5)
    f2.set_ylabel("v3 − v2 (m)", color=off.C_1984)
    f.set_ylabel("mean elevation of land cells (m MHW)")
    f.set_xlabel("domain (1 = south, Cape Hatteras)")
    f.set_xlim(0.2, 90.8)
    f.set_xticks([1] + list(range(10, 91, 10)))
    h1, l1 = f.get_legend_handles_labels()
    h2, l2 = f2.get_legend_handles_labels()
    f.legend(h1 + h2, l1 + l2, loc="upper right", ncol=3, fontsize=8)
    off._title(f, 6, "mean interior elevation at 2004, and the difference")

    fig.suptitle(f"{info['name']}: v2 against v3 — {info['what']}", fontsize=11, x=0.01, ha="left")
    fig.savefig(out, dpi=180, facecolor="white")
    plt.close(fig)
    return out


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
        missing = [v for v in ("v2", "v3") if not (info[v] / f"{info['name']}_run_metadata.json").is_file()]
        if missing:
            print(f"{key}: no run for {missing} at {[str(info[v].relative_to(REPO)) for v in missing]} - skipped")
            continue
        print(f"{key}: loading both runs ...")
        runs = {v: load_run(info[v], info["name"]) for v in ("v2", "v3")}
        results[key] = runs
        figs.append(fig_pair(key, runs, FIG_DIR / f"HAT_compare_v2_v3_{key}.png"))
        tabs.append(table(key, runs))
        print(f"  wrote {figs[-1].name}, {tabs[-1].name}")
    if results:
        rep = report(results, figs, tabs)
        print(f"wrote {rep}")


if __name__ == "__main__":
    main()
