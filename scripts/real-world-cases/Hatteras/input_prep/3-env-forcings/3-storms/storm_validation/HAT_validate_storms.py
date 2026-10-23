"""
Validate a CASCADE storm file against the documented named-storm record, per hindcast period.

    python scripts/input_prep/3-env-forcings/3-storms/storm_validation/HAT_validate_storms.py

Matches model storms to named storms by date window, reports capture and
false positives, and writes a figure, table and report per period. Details: scripts/input_prep/3-env-forcings/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-20
"""
import contextlib
import io
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Anchored 2026-09-14: these literals named trees renamed twice over.
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(_p for _p in Path(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as _REPO_LAYOUT  # noqa: E402

try:
    from HAT_storm_catalog import HISTORICAL_STORMS, storms_in_period
except ImportError as _e:
    raise SystemExit(
        f"\nCould not import HAT_storm_catalog.py  ({_e})\n"
        f"  It must sit in the SAME folder as this script:\n"
        f"    {Path(__file__).parent}\n"
        f"  It is the shared named-storm catalog. Both this script and\n"
        f"  HAT_storm_record_1984_2021.py should import from it, so the list\n"
        f"  cannot drift between them again.\n")

# User configuration

# Periods

# Each period writes its own folder; the 2004 boundary fix is decided, not applied (README)

# Paths are DERIVED from each period's name to match the folder layout

import sys as _envsys
from pathlib import Path as _EnvP
_envsys.path.insert(0, str(next(_q for _q in _EnvP(__file__).resolve().parents
                                if (_q / "hatteras_paths.py").exists())))
from site_layer import hat_env_forcings as _env  # noqa: E402
from site_layer.hat_figure_style import (C_1984, C_1997, INK, INK_MUTED, apply_style, caption,  # noqa: E402
                                         figsize, open_frame, save)
PANEL_COLOURS = (C_1984, C_1997)   # calibration red, test blue
# --- CONFIG ------------------------------------------------------------------
STORM_ROOT = _env.HINDCAST_STORMS

SUMMARY_TEMPLATE     = "{name}_storms_" + _env.DEFAULT_STORM_VARIANT + "_summary.csv"   # the adopted series
VALIDATION_SUBFOLDER = "validation"

PERIODS = [
    {"name": "1997_2009", "label": "Calibration", "begin": 1997, "end": 2009},
    {"name": "2010_2025", "label": "Test",        "begin": 2010, "end": 2025},
]

# Which to run: "all", or a list of names e.g. ["1984_2004"]
RUN_PERIODS = "all"

# Datum

# MHW conversion for the Duck gauge [m]
MHW = 0.36

# Input format: auto, v3 (dam MHW) or readable (m NAVD88)
INPUT_FORMAT = "auto"

# Matching

# Match mode: overlap of the storm window, or the peak time inside it
MATCH_MODE = "overlap"

# The series' storm rule, for the figure (historical_storm_creation_v3_HAT.py)
BERM_NAVD88 = 1.7    # m NAVD88
GROUPING_H = 24      # above-berm hours less than this apart are one storm
MIN_STORM_H = 8      # shorter storms are not kept
# The figure reads each named storm's water within this many days of its closest approach
APPROACH_WINDOW_DAYS = 2
CAPE_HATTERAS = (35.25, -75.53)

# Days of slack either side of the named storm's HURDAT2 active period
MATCH_WINDOW_DAYS = 3

# Sweep these windows to show how sensitive the capture rate is to the choice.
SENSITIVITY_WINDOWS = [0, 1, 2, 3, 5, 7]

# Output

# Each period writes into <STORM_ROOT>/<name>/validation/ 
SAVE_FIGURE = True
SAVE_TABLE  = True
SAVE_REPORT = True
# -----------------------------------------------------------------------------


# Step 1 — load, detect format, normalise

# Which storm-summary schema a table uses
def detect_format(df):
    cols = set(df.columns)
    if {"Storm_Start", "Storm_End", "Rhigh_m"} <= cols:
        return "readable"
    if {"StartTime", "EndTime", "Rhigh"} <= cols:
        return "v3"
    if {"time", "Rhigh", "Rlow", "period", "duration"} == cols:
        raise SystemExit(
            "This is the bare CASCADE storm file (time/Rhigh/Rlow/period/duration).\n"
            "It has no dates, only a model-year index, so it cannot be matched\n"
            "against the historical record. Use the *_summary.csv from the same\n"
            "run instead (v3 writes one alongside)."
        )
    raise SystemExit(f"Unrecognised storm file schema: {sorted(cols)}")


# Normalise any supported schema to
def load_storms(path, fmt, mhw, begin_year, end_year):
    df = pd.read_csv(path)
    fmt = detect_format(df) if fmt == "auto" else fmt

    out = pd.DataFrame()
    if fmt == "v3":
        out["start_ts"] = pd.to_datetime(df["StartTime"])
        out["end_ts"]   = pd.to_datetime(df["EndTime"])
        out["peak_ts"]  = pd.NaT
        out["year"]     = df["calendar_year"].astype(int)
        # v3 stores (TWL - MHW)/10 dam -> back to m NAVD88
        out["Rhigh_m"]  = df["Rhigh"] * 10.0 + mhw
        out["Rlow_m"]   = df["Rlow"] * 10.0 + mhw
        out["duration_h"] = df["duration"]
        out["period_s"]   = df["period"]
    elif fmt == "readable":
        out["start_ts"] = pd.to_datetime(df["Storm_Start"])
        out["end_ts"]   = pd.to_datetime(df["Storm_End"])
        out["peak_ts"]  = (pd.to_datetime(df["Peak_TWL_Time"])
                           if "Peak_TWL_Time" in df.columns else pd.NaT)
        out["year"]     = df["Calendar_Year"].astype(int)
        # already m NAVD88
        out["Rhigh_m"]  = df["Rhigh_m"]
        out["Rlow_m"]   = df["Rlow_m"]
        out["duration_h"] = df.get("Duration_h", np.nan)
        out["period_s"]   = df.get("Wave_Period_s", np.nan)
    else:
        raise SystemExit(f"Unknown INPUT_FORMAT: {fmt}")

    out = out[(out["year"] >= begin_year) & (out["year"] <= end_year)]
    out = out.sort_values("start_ts").reset_index(drop=True)

    print(f"Storm file      : {path.name}")
    print(f"Format detected : {fmt}")
    print(f"Events loaded   : {len(out)}  ({begin_year}-{end_year})")
    print(f"Rhigh range     : {out.Rhigh_m.min():.2f} - {out.Rhigh_m.max():.2f} m NAVD88")
    if fmt == "v3":
        print(f"                  (converted from dam rel. MHW using MHW={mhw} m)")
    return out, fmt


# Step 2 — match

# For each named storm, find the best-matching model event
def match_storms(model, historical, window_days, mode):
    win = pd.Timedelta(days=window_days)
    model = model.copy()
    model["matched_to"] = None
    model["match_type"] = "unmatched"
    model["offset_h"]   = np.nan

    if mode == "peak" and model["peak_ts"].isna().all():
        print("\n  NOTE: MATCH_MODE='peak' but this file has no Peak_TWL_Time."
              "\n        Falling back to 'overlap'.")
        mode = "overlap"

    results = []
    for s in historical:
        lo, hi = s["start_ts"] - win, s["end_ts"] + win
        if mode == "peak":
            hit = (model["peak_ts"] >= lo) & (model["peak_ts"] <= hi)
        else:  # overlap: model window intersects the padded named window
            hit = (model["start_ts"] <= hi) & (model["end_ts"] >= lo)
        cand = model[hit]

        label = f"{s['name']} {s['start_ts'].year}"
        rec = {"label": label, "name": s["name"], "year": s["start_ts"].year,
               "cat": s["cat"], "start": s["start_ts"], "end": s["end_ts"],
               "landfall": s.get("landfall"), "source": s.get("source", "noaa"),
               "note": s.get("note", ""), "matched": False, "matched_idx": None,
               "matched_rhigh": np.nan, "offset_h": np.nan, "shared": False,
               "n_candidates": len(cand)}

        if len(cand):
            idx = cand["Rhigh_m"].idxmax()
            row = model.loc[idx]
            ref = row["peak_ts"] if pd.notna(row["peak_ts"]) else \
                  row["start_ts"] + (row["end_ts"] - row["start_ts"]) / 2
            # hours outside the named storm's active window; 0 = inside
            if ref < s["start_ts"]:
                off = (ref - s["start_ts"]).total_seconds() / 3600
            elif ref > s["end_ts"]:
                off = (ref - s["end_ts"]).total_seconds() / 3600
            else:
                off = 0.0

            already = model.loc[idx, "match_type"] == "matched"
            rec.update(matched=True, matched_idx=idx,
                       matched_rhigh=row["Rhigh_m"], offset_h=off,
                       shared=already)
            if not already:
                model.loc[idx, ["matched_to", "match_type", "offset_h"]] = \
                    [label, "matched", off]
            else:
                model.loc[idx, "matched_to"] = \
                    f"{model.loc[idx, 'matched_to']} + {label}"
        results.append(rec)
    return model, results


# (matched, total) named storms at one window
def capture_rate(model, historical, window_days, mode):
    _, res = match_storms(model, historical, window_days, mode)
    return sum(r["matched"] for r in res), len(res)


# Step 3 — report

# The period's matching report to the console
def print_report(model, results, window_days, mode):
    matched = [r for r in results if r["matched"]]
    missed  = [r for r in results if not r["matched"]]
    shared  = [r for r in results if r["shared"]]
    n_un    = int((model["match_type"] == "unmatched").sum())

    print("\n" + "=" * 74)
    print(f"STORM VALIDATION SUMMARY   (mode={mode}, window=+/-{window_days} d)")
    print("=" * 74)
    print(f"Named storms in period :  {len(results)}")
    print(f"  captured             :  {len(matched)}  "
          f"({100*len(matched)/len(results):.0f}%)")
    print(f"  missed               :  {len(missed)}")
    print(f"Model events           :  {len(model)}")
    print(f"  matched to a name    :  {int((model['match_type']=='matched').sum())}")
    print(f"  unmatched            :  {n_un}  <- expected: nor'easters/unnamed")
    if shared:
        print(f"\n  {len(shared)} named storm(s) share a model event with another storm.")
        print("  A merged mega-event can span several names — the model cannot")
        print("  distinguish them, so 'captured' overstates resolution here:")
        for r in shared:
            print(f"    {r['label']:<18} shares event #{r['matched_idx']}")

    print("\n" + "=" * 74)
    print("CAPTURED")
    print("=" * 74)
    print(f"  {'storm':<18}{'cat':<5}{'Rhigh_m':>9}{'offset_h':>10}   flags")
    for r in sorted(matched, key=lambda x: x["start"]):
        flags = []
        if r["landfall"]:            flags.append("landfall")
        if r["source"] == "local":   flags.append("source=local")
        if r["shared"]:              flags.append("SHARED EVENT")
        if r["n_candidates"] > 1:    flags.append(f"{r['n_candidates']} candidates")
        print(f"  {r['label']:<18}{r['cat']:<5}{r['matched_rhigh']:>9.3f}"
              f"{r['offset_h']:>+10.0f}   {', '.join(flags)}")

    print("\n" + "=" * 74)
    print("MISSED")
    print("=" * 74)
    if not missed:
        print("  none")
    else:
        for r in sorted(missed, key=lambda x: x["start"]):
            flags = []
            if r["landfall"]:          flags.append("LANDFALL — investigate")
            if r["landfall"] is False: flags.append("passed by, minor impact")
            if r["landfall"] is None:  flags.append("no local-impact narrative")
            if r["source"] == "local": flags.append("source=local")
            print(f"  {r['label']:<18}{r['cat']:<5} {r['start']:%Y-%m-%d} -> "
                  f"{r['end']:%Y-%m-%d}   {', '.join(flags)}")
        hard = [r for r in missed if r["landfall"]]
        if hard:
            print(f"\n  {len(hard)} missed storm(s) have DOCUMENTED Hatteras impact.")
            print("  These are the ones that matter — a documented landfall with no")
            print("  model event is either a water-level gap or a threshold problem:")
            for r in hard:
                print(f"    {r['label']}: {r['note']}")


# Capture rate against the match window
def print_sensitivity(model, historical, mode, windows):
    print("\n" + "=" * 74)
    print("MATCH-WINDOW SENSITIVITY")
    print("=" * 74)
    print("  Capture rate vs window. A rate that climbs steeply with the window")
    print("  means matches are being bought with slack, not skill.")
    print(f"\n  {'window (d)':<12}{'captured':<12}{'rate'}")
    for w in windows:
        n, tot = capture_rate(model, historical, w, mode)
        bar = "#" * int(round(20 * n / tot))
        print(f"  {w:<12}{f'{n}/{tot}':<12}{100*n/tot:>3.0f}%  {bar}")


# Step 4 — figure

# When each named storm was closest to Cape Hatteras (HURDAT2, any status): {(name, year): (time, km)}
def closest_approach(hurdat2_file):
    rows, name = [], None
    for line in open(hurdat2_file):
        p = [x.strip() for x in line.split(",")]
        if p[0][:2].isalpha():
            name = p[1].title()
            continue
        if int(p[0][:4]) < 1984:                       # before the record; some old lines are malformed
            continue
        lat = float(p[4][:-1]) * (1 if p[4][-1] == "N" else -1)
        lon = float(p[5][:-1]) * (-1 if p[5][-1] == "W" else 1)
        rows.append((name, pd.Timestamp(p[0] + p[1]), lat, lon))
    h = pd.DataFrame(rows, columns=["name", "t", "lat", "lon"])
    la, lo = np.radians(h.lat), np.radians(h.lon)
    la0, lo0 = np.radians(CAPE_HATTERAS[0]), np.radians(CAPE_HATTERAS[1])
    h["km"] = 6371.0 * 2 * np.arcsin(np.sqrt(np.sin((la - la0) / 2) ** 2
                                             + np.cos(la) * np.cos(la0) * np.sin((lo - lo0) / 2) ** 2))
    h["year"] = h.t.dt.year
    near = h.loc[h.groupby(["name", "year"]).km.idxmin()]
    return {(r.name, r.year): (r.t, r.km) for r in near.itertuples()}


# Each named storm's water: its highest above-berm spell (grouped as the series is) around its closest approach
def storm_water(results):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    with contextlib.redirect_stdout(io.StringIO()):
        import storm_figures as sf
    twl = sf.load_forcing().TWL
    near = closest_approach(sf.HURDAT2_FILE)
    up = twl[twl > BERM_NAVD88]
    gap = np.diff(up.index.values).astype("timedelta64[h]").astype(int)
    gid = np.concatenate([[0], np.cumsum(gap >= GROUPING_H)])
    spells = (pd.DataFrame({"t": up.index, "twl": up.values, "g": gid}).groupby("g")
              .agg(start=("t", "first"), end=("t", "last"), hours=("t", "size"), peak=("twl", "max")))
    win = pd.Timedelta(days=APPROACH_WINDOW_DAYS)
    for r in results:
        t, km = near.get((r["name"], r["year"]), (None, np.nan))
        # not in HURDAT2 (the 2017 unnamed system): its catalog dates
        lo, hi = (t - win, t + win) if t is not None else (r["start"], r["end"] + pd.Timedelta(hours=23))
        r["closest_approach"], r["closest_km"] = t, round(float(km), 0) if np.isfinite(km) else np.nan
        c = spells[(spells.start <= hi) & (spells.end >= lo)]
        if len(c):
            s = c.loc[c.peak.idxmax()]
            r["hours_above_berm"], r["peak_twl_m_navd88"] = int(s.hours), round(float(s.peak), 3)
        else:
            r["hours_above_berm"] = 0
            r["peak_twl_m_navd88"] = round(float(twl.loc[lo:hi].max()), 3)
        # The series keeps a spell of MIN_STORM_H hours or more, so this is exact; "matched" is the looser date match
        r["in_series"] = r["hours_above_berm"] >= MIN_STORM_H
        r["reason"] = ("in the model" if r["in_series"] else
                       "never above the berm" if r["hours_above_berm"] == 0 else
                       f"above the berm {r['hours_above_berm']} h, under the {MIN_STORM_H} h minimum")
    return results


# Place labels next to their dots, clear of each other and of the axes edge
def _labels(ax, xs, ys, texts):
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    box, placed = ax.get_window_extent(rend), []
    offsets = [(6, 0), (-6, 0), (0, 7), (0, -7), (8, 8), (-8, 8), (8, -8), (-8, -8), (16, 0), (-16, 0),
               (0, 15), (0, -15), (18, 12), (-18, 12), (18, -12), (-18, -12)]
    for x, y, t in sorted(zip(xs, ys, texts), key=lambda z: -z[1]):
        best = None
        for dx, dy in offsets:
            ha = "left" if dx > 0 else "right" if dx < 0 else "center"
            va = "bottom" if dy > 0 else "top" if dy < 0 else "center"
            a = ax.annotate(t, (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, va=va, fontsize=6.8)
            a.update_positions(rend)
            bb = a.get_window_extent(rend).expanded(1.04, 1.1)
            a.remove()
            inside = box.x0 <= bb.x0 and bb.x1 <= box.x1 and box.y0 <= bb.y0 and bb.y1 <= box.y1
            cost = (0 if inside else 1e6) + 1e3 * sum(bb.overlaps(o) for o in placed) + abs(dx) + abs(dy)
            if best is None or cost < best[0]:
                best = (cost, dx, dy, ha, va, bb)
        _, dx, dy, ha, va, bb = best
        placed.append(bb)
        ax.annotate(t, (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, va=va, fontsize=6.8,
                    color=INK, arrowprops=(dict(arrowstyle="-", lw=0.4, color=INK_MUTED, shrinkA=0, shrinkB=2)
                                           if abs(dx) + abs(dy) > 12 else None))


# Peak water level against hours above the berm, one dot per named storm, calibration and test side by side
def make_figure(runs, path):
    from matplotlib.lines import Line2D
    apply_style()
    fig, axes = plt.subplots(1, len(runs), figsize=figsize("double", height=3.9), sharey=True,
                             constrained_layout=True)
    axes = np.atleast_1d(axes)
    zero_x = 0.75                                    # where a storm with no hour above the berm sits
    allres = pd.concat([pd.DataFrame(r["results"]) for r in runs])
    ymin, ymax = min(1.0, allres.peak_twl_m_navd88.min() - 0.15), allres.peak_twl_m_navd88.max() + 0.35
    for ax, run, colour in zip(axes, runs, PANEL_COLOURS):
        res = pd.DataFrame(run["results"])
        x = res.hours_above_berm.where(res.hours_above_berm > 0, zero_x)
        ax.axvspan(MIN_STORM_H, 400, color=colour, alpha=0.06, lw=0, zorder=0)
        ax.axhline(BERM_NAVD88, color=INK_MUTED, lw=0.8, ls=(0, (4, 2)), zorder=1)
        ax.axvline(MIN_STORM_H, color=INK_MUTED, lw=0.8, ls=(0, (4, 2)), zorder=1)
        hit = res.in_series.astype(bool)
        ax.scatter(x[hit], res.peak_twl_m_navd88[hit], s=34, color=colour, edgecolor="white", lw=0.5, zorder=3)
        ax.scatter(x[~hit], res.peak_twl_m_navd88[~hit], s=34, facecolor="white", edgecolor=colour, lw=1.3,
                   zorder=3)
        ax.set_xscale("log")
        ax.set_xlim(0.55, 330)
        ax.set_ylim(ymin, ymax)
        ax.set_yticks([v for v in range(1, int(ymax) + 1) if v >= ymin])
        ax.set_xticks([zero_x, 2, 4, MIN_STORM_H, 24, 48, 96, 192])
        ax.set_xticklabels(["0", "2", "4", str(MIN_STORM_H), "24", "48", "96", "192"])
        ax.minorticks_off()
        ax.set_xlabel("Hours above the berm")
        c = run["cfg"]
        ax.text(0.98, 0.02, f"{c['label'].lower()} {c['begin']}–{c['end']}\n{int(hit.sum())} of {len(res)} in the model",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color=INK)
        ax.text(MIN_STORM_H * 1.06, ymax - 0.05, "kept as a storm", fontsize=7, color=INK_MUTED, va="top")
        ax.text(320, BERM_NAVD88 + 0.03, f"berm {BERM_NAVD88} m", fontsize=7, color=INK_MUTED, va="bottom", ha="right")
        open_frame(ax)
        up = res.hours_above_berm > 0
        _labels(ax, list(x[up]), list(res.peak_twl_m_navd88[up]), list(res.label[up]))
        if (~up).any():
            ax.text(0.03, 0.82, "never above the berm:\n" + "\n".join(res.label[~up]), transform=ax.transAxes,
                    ha="left", va="top", fontsize=6.8, color=INK_MUTED)
    axes[0].set_ylabel("Peak total water level (m NAVD88)")
    axes[0].legend(handles=[Line2D([], [], marker="o", ls="", color=INK, ms=6, label="in the model"),
                            Line2D([], [], marker="o", ls="", mfc="white", mec=INK, mew=1.2, ms=6, label="missed")],
                   loc="upper left", frameon=False, fontsize=7.5)
    tot = [pd.DataFrame(r["results"]) for r in runs]
    caption(fig, (
        f"Named storms whose track passed within 65 nmi of Cape Hatteras (HAT_storm_catalog.py) against the "
        f"model's storm series ({_env.DEFAULT_STORM_VARIANT}), "
        + "; ".join(f"{r['cfg']['label'].lower()} {r['cfg']['begin']}–{r['cfg']['end']}" for r in runs)
        + f". Each dot is one named storm: its highest spell of total water level (Duck gauge + Stockdon R2% from "
        f"WIS waves) above the {BERM_NAVD88} m NAVD88 berm within ±{APPROACH_WINDOW_DAYS} days of its closest approach "
        f"to Cape Hatteras (HURDAT2), "
        f"with hours less than {GROUPING_H} h apart grouped as in the series. x: the spell's hours above the berm "
        f"(log; storms never above it sit at 0, at their highest water, and are listed); y: its peak. Filled: "
        f"the spell is in the model's series (it lasted at least {MIN_STORM_H} h); open: missed. "
        f"A storm is kept only if it is above the berm for at least {MIN_STORM_H} h (shaded). Missed storms left "
        f"of the line topped the berm too briefly; those below the berm line never reached it. "
        + " ".join(f"{r['cfg']['label']}: {int(t.in_series.sum())} of {len(t)} in the model." for r, t in zip(runs, tot))
        + " The report's date match (any model event within ±" + str(MATCH_WINDOW_DAYS) + " days of a storm's "
        "catalog dates) is looser and can credit a storm with a neighbouring event; the figure does not use it."))
    path.parent.mkdir(parents=True, exist_ok=True)
    out = save(fig, path, bbox_inches="tight")
    print(f"\nFigure: {out[0]}")
    plt.close(fig)


# Path resolution

# Derive (summary_file, output_dir) from the period name + folder layout
def resolve_paths(cfg):
    folder = Path(cfg["folder"]) if "folder" in cfg else STORM_ROOT / cfg["name"]
    file   = (Path(cfg["file"]) if "file" in cfg
              else folder / SUMMARY_TEMPLATE.format(name=cfg["name"]))
    # Beside the storm_check output for the same window
    outdir = (Path(cfg["outdir"]) if "outdir" in cfg
              else _env.STORM_VALIDATION / cfg["name"])
    return folder, file, outdir


# Print every resolved path and flag missing ones BEFORE any work starts
def preflight(periods):
    print("\n" + "=" * 74)
    print("PATHS")
    print("=" * 74)
    print(f"  STORM_ROOT : {STORM_ROOT}")
    print(f"               {'exists' if STORM_ROOT.exists() else 'DOES NOT EXIST <-- fix this first'}")
    ok = STORM_ROOT.exists()
    for cfg in periods:
        folder, file, outdir = resolve_paths(cfg)
        print(f"\n  [{cfg['name']}]")
        print(f"    input  : {file}")
        print(f"             {'found' if file.exists() else 'NOT FOUND'}")
        print(f"    output : {outdir}")
        if not file.exists():
            ok = False
            if folder.exists():
                found = sorted(p.name for p in folder.glob("*.csv"))
                print(f"    -> folder exists. CSVs in it:")
                for f in found or ["(none)"]:
                    print(f"         {f}")
                print(f"    -> expected '{SUMMARY_TEMPLATE.format(name=cfg['name'])}'."
                      f" Adjust SUMMARY_TEMPLATE, or add an explicit"
                      f' "file" key for this period.')
            else:
                print(f"    -> folder does not exist: {folder}")
    print("\n" + "=" * 74)
    return ok


# Run one period

# Validate one period and write its outputs
def run_period(cfg):
    name, begin, end = cfg["name"], cfg["begin"], cfg["end"]
    folder, file, outdir = resolve_paths(cfg)

    print("\n" + "#" * 74)
    print(f"# {cfg['label']}   ({begin}-{end})")
    print(f"#   in  : {file}")
    print(f"#   out : {outdir}")
    print("#" * 74)

    if not file.exists():
        print(f"\n  SKIPPED — summary file not found.")
        print("  v3 writes <save_name>_summary.csv alongside <save_name>.csv;")
        print("  the bare .csv/.npy have only a model-year index and cannot be")
        print("  date-matched. Re-run storm creation with save_dfs = True if the")
        print("  summary is missing.")
        return None

    model, fmt = load_storms(file, INPUT_FORMAT, MHW, begin, end)
    historical = storms_in_period(begin, end)
    n_local = sum(1 for s in historical if s.get("source") == "local")
    print(f"Named storms    : {len(historical)}  ({n_local} source='local')")

    model, results = match_storms(model, historical, MATCH_WINDOW_DAYS, MATCH_MODE)
    print_report(model, results, MATCH_WINDOW_DAYS, MATCH_MODE)
    print_sensitivity(model, historical, MATCH_MODE, SENSITIVITY_WINDOWS)

    results = storm_water(results)
    outdir.mkdir(parents=True, exist_ok=True)
    if SAVE_TABLE:
        tp = outdir / f"match_table_{name}.csv"
        pd.DataFrame(results).to_csv(tp, index=False)
        print(f"\nMatch table : {tp}")

    return {"cfg": cfg, "model": model, "results": results, "outdir": outdir}


# Cross-period checks

# Compare the periods, flagging any overlapping year
def cross_period_report(runs):
    if len(runs) < 2:
        return
    print("\n" + "=" * 74)
    print("CROSS-PERIOD SUMMARY")
    print("=" * 74)
    print(f"  {'period':<14}{'years':<13}{'events':>8}{'captured':>11}{'max Rhigh':>12}")
    for r in runs:
        c = r["cfg"]
        n_cap = sum(x["matched"] for x in r["results"])
        n_tot = len(r["results"])
        yrs = f"{c['begin']}-{c['end']}"
        cap = f"{n_cap}/{n_tot}"
        print(f"  {c['name']:<14}{yrs:<13}{len(r['model']):>8}{cap:>11}"
              f"{r['model'].Rhigh_m.max():>11.2f} m")

    # Boundary overlap
    print("\n  Boundary check:")
    for a, b in zip(runs, runs[1:]):
        ca, cb = a["cfg"], b["cfg"]
        lo, hi = max(ca["begin"], cb["begin"]), min(ca["end"], cb["end"])
        if lo > hi:
            print(f"    {ca['name']} / {cb['name']}: no overlap — clean.")
            continue
        ov_a = a["model"][(a["model"].year >= lo) & (a["model"].year <= hi)]
        ov_b = b["model"][(b["model"].year >= lo) & (b["model"].year <= hi)]
        print(f"    {ca['name']} / {cb['name']}: years {lo}-{hi} appear in BOTH")
        print(f"      {len(ov_a)} events in {ca['name']}, {len(ov_b)} in {cb['name']}")
        key = ["Rhigh_m", "Rlow_m", "duration_h", "period_s"]
        same = (len(ov_a) == len(ov_b) and np.allclose(
            ov_a[key].values, ov_b[key].values, atol=1e-6, equal_nan=True))
        print(f"      identical: {same}"
              + ("  (deterministic — but they will be forced twice)" if same
                 else "  <- DIFFER: the boundary is cutting events differently"))
        if len(ov_a):
            print(f"      -> if Period 2 initialises from Period 1's final state,")
            print(f"         these {len(ov_a)} events are applied twice. v1 ended")
            print(f"         Period 1 at 2003-12-31 23:00 to avoid this.")


# Run: every chosen period, then the cross-period summary
def main():
    print("\nHAT_validate_storms.py")
    print(f"Catalog: {len(HISTORICAL_STORMS)} named storms (HAT_storm_catalog.py)")

    wanted = ([p for p in PERIODS] if RUN_PERIODS == "all"
              else [p for p in PERIODS if p["name"] in RUN_PERIODS])
    if not wanted:
        raise SystemExit(f"RUN_PERIODS={RUN_PERIODS!r} matched no entry in PERIODS.")

    if not preflight(wanted):
        print("\n  One or more paths above are missing. Periods with a missing")
        print("  input are skipped; the rest still run.\n")

    runs = []
    for cfg in wanted:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            r = run_period(cfg)
        text = buf.getvalue()
        print(text, end="")
        if r is not None:
            if SAVE_REPORT:
                rp = r["outdir"] / f"report_{cfg['name']}.txt"
                rp.parent.mkdir(parents=True, exist_ok=True)
                rp.write_text(text, encoding="utf-8")
                print(f"Report      : {rp}")
            runs.append(r)

    cross_period_report(runs)
    if SAVE_FIGURE and runs:
        make_figure(runs, _env.STORM_VALIDATION / f"storm_validation_{'_'.join(r['cfg']['name'] for r in runs)}.png")
    print("\n" + "=" * 74)
    for r in runs:
        print(f"  {r['cfg']['name']}: {r['outdir']}")
    print("=" * 74)
    return runs


if __name__ == "__main__":
    all_runs = main()
