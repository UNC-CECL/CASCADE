# =============================================================================
# HAT_validate_storms.py
# Unified storm-record validation — Hatteras Island
# -----------------------------------------------------------------------------
# Compares a CASCADE storm file against the documented record of named tropical
# and extratropical storms near Hatteras Island, NC.
#
# REPLACES: HAT_validate_storms.py + HAT_validate_storms_windowchange.py
#           (which differed only in period config and raw_offset convention)
#
# -----------------------------------------------------------------------------
# WHAT CHANGED, AND WHY:
#
#   Catalog moved out          The list lived in three files and had drifted:
#                              BERTHA/FRAN 1996 and ISAIAS 2020 were on the
#                              record chart but tested by nothing. Now imported
#                              from HAT_storm_catalog.py. Do not paste it back.
#
#   Format auto-detection      Reads BOTH schemas:
#                                - v3 summary (Lexi's): StartTime, EndTime,
#                                  calendar_year, Rhigh [dam rel. MHW]
#                                - HAT_create_storms readable: Storm_Start,
#                                  Storm_End, Peak_TWL_Time, Rhigh_m [m NAVD88]
#                              The original validator required Peak_TWL_Time, so it
#                              could not read v3 comparison at all.
#
#   Datum normalisation        v3 stores (TWL - MHW)/10 dam; HAT_create_storms
#                              stores TWL/10 dam. Everything is converted to
#                              m NAVD88 internally so numbers are comparable
#                              across formats. Get MHW right or the axis lies.
#
#   Overlap matching           Default MATCH_MODE="overlap": the model storm
#                              WINDOW must intersect the named storm window
#                              (+/- MATCH_WINDOW_DAYS). Works for every format,
#                              needs no Peak_TWL_Time, and is the better test:
#                              a storm is captured if the model has an event
#                              running at the same time, not if a single
#                              instant lands inside a fuzzy box.
#                              MATCH_MODE="peak" reproduces the original behaviour
#                              where Peak_TWL_Time exists.
#
#   Shared-match reporting     A merged mega-event can span several named
#                              storms (the 2003 event runs 09-09 -> 09-19 and
#                              could claim more than Isabel). The original code
#                              silently let the first claimant win while still
#                              reporting every storm as matched. Shared matches
#                              are now counted and printed.
#
# Usage:
#   Set STORM_FILE / BEGIN_YEAR / END_YEAR below and run.
#
# Author: Hannah Henry
# =============================================================================

import contextlib
import io
from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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

# =============================================================================
# USER CONFIGURATION
# =============================================================================

# --- Periods ------------------------------------------------------------------
# Each period validates independently and writes to its own subfolder under
# validation/ subfolder, so the two never overwrite each other and each can
# be checked on its own.
#
# BOUNDARY — DECIDED BUT NOT YET APPLIED.
# 2004 currently appears in BOTH storm files, with the same 10 events. If
# Period 2 initialises from Period 1's final morphological state, those 10 are
# forced twice. Agreed fix: Period 1 ends 2003-12-31 23:00, Period 2 begins
# 2004-01-01 — matching historical_storm_creation_HAT.py (v1), which ended
# Period 1 at 2003-12-31 23:00 for exactly this reason.
#
# Config below matches the files AS THEY EXIST TODAY (Period 1 -> 2004), so the
# script runs against them unchanged. The CROSS-PERIOD SUMMARY will keep
# reporting the 2004 overlap every run — that is deliberate, and it stops when
# the input is actually regenerated.
#
# WHEN YOU REGENERATE, change in TWO places:
#   1. historical_storm_creation_v3_HAT.py -> end_time = '2003-12-31 23:00:00'
#                                             save_name = "1984_2003_storms_v3"
#   2. here -> "file": "1984_2003_storms_v3_summary.csv",  "end": 2003
# Note this also drops Alex/Bonnie 2004 from Period 1's catalog slice, so its
# denominator goes 22 -> 20 and the capture rate shifts for that reason alone.
#
# Point these at the *_summary.csv files, NOT the bare CASCADE files. The bare
# files have only a model-year index and cannot be date-matched.

# Paths are DERIVED from each period's name to match the folder layout:
#
#   <STORM_ROOT>/
#     1984_2004/
#       1984_2004_storms_v3_summary.csv      <- input  (SUMMARY_TEMPLATE)
#       validation/                          <- comparison (VALIDATION_SUBFOLDER)
#     2004_2024/
#       2004_2024_storms_v3_summary.csv
#       validation/
#
# So renaming a period (e.g. 1984_2004 -> 1984_2003) updates the input path AND
# the comparison folder together. Override per period with explicit "file" /
# "outdir" keys if a run ever sits outside this convention.

STORM_ROOT = Path(r"/data/hatteras_init/storms/hindcast_storms")

SUMMARY_TEMPLATE     = "{name}_storms_v3_summary.csv"
VALIDATION_SUBFOLDER = "validation"

PERIODS = [
    {"name": "1984_2004", "label": "Period 1 - calibration", "begin": 1984, "end": 2004},
    {"name": "2004_2024", "label": "Period 2 - test",  "begin": 2004, "end": 2024},
]

# Which to run: "all", or a list of names e.g. ["1984_2004"]
RUN_PERIODS = "all"

# --- Datum --------------------------------------------------------------------
# MHW conversion for the Duck gauge [m]: 0 m NAVD88 = MHW m MHW.
# Used to convert v3's (TWL - MHW)/10 dam back to m NAVD88. Must match the MHW
# in the storm-creation config or every Rhigh here is shifted.
MHW = 0.36

# "auto"  -> infer from columns (recommended)
# "v3"    -> StartTime/EndTime/calendar_year, Rhigh in dam relative to MHW
# "readable" -> HAT_create_storms.py readable CSV, Rhigh_m in m NAVD88
INPUT_FORMAT = "auto"

# --- Matching -----------------------------------------------------------------
# "overlap" -> model storm window intersects [named_start - W, named_end + W]
# "peak"    -> Peak_TWL_Time falls in that range (needs the readable format)
MATCH_MODE = "overlap"

# Days of slack either side of the named storm's HURDAT2 active period.
#   (1) Duck gauge (8651370) is ~40 km north of Hatteras centre — surge timing
#       differs by hours depending on approach direction.
#   (2) HURDAT2 reports full lifetime at sea, not the Hatteras influence window.
# Recommended 3-5. Larger windows inflate the capture rate; see the sensitivity
# sweep printed at the end.
MATCH_WINDOW_DAYS = 3

# Sweep these windows to show how sensitive the capture rate is to the choice.
SENSITIVITY_WINDOWS = [0, 1, 2, 3, 5, 7]

# --- Output -------------------------------------------------------------------
# Each period writes into <STORM_ROOT>/<name>/validation/ :
#     storm_check_<name>.png     figure
#     match_table_<name>.csv     per-named-storm match detail
#     report_<name>.txt          the full console report
SAVE_FIGURE = True
SAVE_TABLE  = True
SAVE_REPORT = True
SHOW_FIGURE = True   # False when running both periods, or two windows open

# =============================================================================
# COLOUR SCHEME
# =============================================================================

CAT_COLORS = {"H5": "#7b0000", "H4": "#b22222", "H3": "#e05c1a",
              "H2": "#e8a020", "H1": "#f0d040", "TS": "#5590d0", "ET": "#8888aa"}
CAT_ORDER  = ["H5", "H4", "H3", "H2", "H1", "TS", "ET"]

MATCHED_COLOR   = "#2ca02c"   # model event matched a named storm
UNMATCHED_COLOR = "#b0b0b0"   # model event is a nor'easter / unnamed
MISSED_COLOR    = "#d62728"   # named storm not captured


# =============================================================================
# STEP 1 — LOAD, DETECT FORMAT, NORMALISE
# =============================================================================

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


def load_storms(path, fmt, mhw, begin_year, end_year):
    """
    Normalise any supported schema to:
        start_ts, end_ts, peak_ts (NaT if absent), year,
        Rhigh_m, Rlow_m [m NAVD88], duration_h, period_s
    """
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


# =============================================================================
# STEP 2 — MATCH
# =============================================================================

def match_storms(model, historical, window_days, mode):
    """
    For each named storm, find the best-matching model event.
    Best = highest Rhigh among candidates. Returns (model, results).
    """
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


def capture_rate(model, historical, window_days, mode):
    _, res = match_storms(model, historical, window_days, mode)
    return sum(r["matched"] for r in res), len(res)


# =============================================================================
# STEP 3 — REPORT
# =============================================================================

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


# =============================================================================
# STEP 4 — FIGURE
# =============================================================================

def make_figure(model, results, window_days, mode, begin_year, end_year, path):
    fig = plt.figure(figsize=(16, 11))
    gs  = fig.add_gridspec(3, 1, height_ratios=[1, 1.5, 1.1], hspace=0.32)

    # --- Panel 1: annual counts, stacked matched/unmatched --------------------
    ax1 = fig.add_subplot(gs[0])
    years = np.arange(begin_year, end_year + 1)
    m_cnt = [int(((model.year == y) & (model.match_type == "matched")).sum()) for y in years]
    u_cnt = [int(((model.year == y) & (model.match_type == "unmatched")).sum()) for y in years]
    ax1.bar(years, u_cnt, color=UNMATCHED_COLOR, label="Unnamed / nor'easter")
    ax1.bar(years, m_cnt, bottom=u_cnt, color=MATCHED_COLOR, label="Matched to named storm")
    ax1.axhline(5, color="#B71C1C", ls="--", lw=1.2, zorder=5)
    ax1.text(0.995, 5, " Barrier3D ~5/yr ", transform=ax1.get_yaxis_transform(),
             ha="right", va="bottom", fontsize=9, color="#B71C1C",
             bbox=dict(fc="white", ec="none", pad=1))
    ax1.set_ylabel("Storms per year")
    ax1.set_title(f"Model storm events vs named-storm record, {begin_year}-{end_year}"
                  f"   (mode={mode}, window=+/-{window_days} d)",
                  fontsize=13, fontweight="bold", loc="left")
    ax1.legend(loc="upper left", fontsize=9, framealpha=0.9)
    ax1.set_xlim(begin_year - 0.6, end_year + 0.6)

    # --- Panel 2: timeline, Rhigh vs date ------------------------------------
    ax2 = fig.add_subplot(gs[1])
    for r in results:
        ax2.axvspan(r["start"], r["end"], color=CAT_COLORS.get(r["cat"], "#999"),
                    alpha=0.16, lw=0)
    mm = model[model.match_type == "matched"]
    uu = model[model.match_type == "unmatched"]
    ax2.scatter(uu.start_ts, uu.Rhigh_m, s=26, c=UNMATCHED_COLOR,
                edgecolor="none", label="Unnamed / nor'easter", zorder=3)
    ax2.scatter(mm.start_ts, mm.Rhigh_m, s=52, c=MATCHED_COLOR,
                edgecolor="#1a1a1a", linewidth=0.4, label="Matched", zorder=4)
    for r in results:
        if not r["matched"]:
            ax2.axvline(r["start"], color=MISSED_COLOR, lw=1.6, alpha=0.85, zorder=2)
            ax2.annotate(r["label"], (r["start"], ax2.get_ylim()[1]),
                         rotation=90, fontsize=7.5, color=MISSED_COLOR,
                         ha="right", va="top")
    top = model.nlargest(6, "Rhigh_m")
    for _, row in top.iterrows():
        if row.matched_to:
            ax2.annotate(str(row.matched_to), (row.start_ts, row.Rhigh_m),
                         xytext=(0, 9), textcoords="raw_offset points",
                         ha="center", fontsize=8.5, fontweight="bold")
    ax2.set_ylabel("Rhigh [m NAVD88]")
    ax2.legend(loc="upper left", fontsize=9, framealpha=0.9)
    ax2.set_title("Shaded bands = named-storm active periods; red lines = missed",
                  fontsize=10, loc="left", color="#555")

    # --- Panel 3: capture table ----------------------------------------------
    ax3 = fig.add_subplot(gs[2]); ax3.axis("off")
    rows, colors = [], []
    for r in sorted(results, key=lambda x: x["start"]):
        rows.append([r["label"], r["cat"],
                     "yes" if r["matched"] else "NO",
                     f"{r['matched_rhigh']:.2f}" if r["matched"] else "-",
                     f"{r['offset_h']:+.0f}" if r["matched"] else "-",
                     ("landfall" if r["landfall"] else
                      "minor" if r["landfall"] is False else "-")
                     + (" | local" if r["source"] == "local" else "")
                     + (" | shared" if r["shared"] else "")])
        colors.append(["#eaf7ea" if r["matched"] else "#fdeaea"] * 6)
    tbl = ax3.table(cellText=rows,
                    colLabels=["Named storm", "Cat", "Captured", "Rhigh_m",
                               "Offset (h)", "Notes"],
                    cellColours=colors, loc="upper center",
                    colWidths=[0.18, 0.06, 0.09, 0.09, 0.09, 0.30])
    tbl.auto_set_font_size(False); tbl.set_fontsize(7.6); tbl.scale(1, 1.16)
    for j in range(6):
        tbl[0, j].set_facecolor("#1c2b39")
        tbl[0, j].set_text_props(color="white", fontweight="bold")

    handles = [mpatches.Patch(color=CAT_COLORS[c], label=c) for c in CAT_ORDER]
    handles += [mlines.Line2D([], [], color=MISSED_COLOR, lw=2, label="Missed")]
    ax1.legend(handles=handles + ax1.get_legend_handles_labels()[0],
               loc="upper left", fontsize=8, ncol=5, framealpha=0.9)

    if SAVE_FIGURE:
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=180, bbox_inches="tight")
        print(f"\nFigure: {path}")
    return fig


# =============================================================================
# PATH RESOLUTION
# =============================================================================

def resolve_paths(cfg):
    """Derive (summary_file, output_dir) from the period name + folder layout."""
    folder = Path(cfg["folder"]) if "folder" in cfg else STORM_ROOT / cfg["name"]
    file   = (Path(cfg["file"]) if "file" in cfg
              else folder / SUMMARY_TEMPLATE.format(name=cfg["name"]))
    outdir = (Path(cfg["outdir"]) if "outdir" in cfg
              else folder / VALIDATION_SUBFOLDER)
    return folder, file, outdir


def preflight(periods):
    """Print every resolved path and flag missing ones BEFORE any work starts."""
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


# =============================================================================
# RUN ONE PERIOD
# =============================================================================

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

    outdir.mkdir(parents=True, exist_ok=True)
    if SAVE_TABLE:
        tp = outdir / f"match_table_{name}.csv"
        pd.DataFrame(results).to_csv(tp, index=False)
        print(f"\nMatch table : {tp}")

    fig = make_figure(model, results, MATCH_WINDOW_DAYS, MATCH_MODE, begin, end,
                      outdir / f"storm_check_{name}.png")
    if SHOW_FIGURE:
        plt.show()
    else:
        plt.close(fig)

    return {"cfg": cfg, "model": model, "results": results, "outdir": outdir}


# =============================================================================
# CROSS-PERIOD CHECKS
# =============================================================================

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

    # --- boundary overlap -----------------------------------------------------
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


# =============================================================================
# MAIN
# =============================================================================

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
    print("\n" + "=" * 74)
    for r in runs:
        print(f"  {r['cfg']['name']}: {r['outdir']}")
    print("=" * 74)
    return runs


if __name__ == "__main__":
    all_runs = main()
