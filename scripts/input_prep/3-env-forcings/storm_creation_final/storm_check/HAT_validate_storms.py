# =============================================================================
# HAT_validate_storms.py
# Historical Storm Record Validation — Hatteras Island
# -----------------------------------------------------------------------------
# Description:
#   Compares a CASCADE-format storm readable CSV (comparison of HAT_create_storms.py)
#   against the documented historical record of named tropical and extratropical
#   storms that passed within 60 nautical miles of Hatteras Island, NC. The
#   historical catalog is sourced from NOAA's Hurricane Research Division
#   Atlantic hurricane database (HURDAT2), filtered to storms reaching at least
#   tropical storm strength.
#
#   The comparison uses a configurable time window to account for two sources
#   of lag between the historical record and the Duck, NC gauge-derived TWL:
#     (1) The Duck gauge (Station 8651370) is ~40 km north of Hatteras Island's
#         center. For storms approaching from the south, peak surge at Hatteras
#         can precede peak surge at Duck by several hours. For nor'easters
#         approaching from the northeast, the Duck gauge leads Hatteras.
#     (2) The NOAA historical catalog reports the storm's full lifetime, not
#         only the period of direct Hatteras influence. A storm active for two
#         weeks may only affect Hatteras for 24–48 hours within that window.
#   A ±MATCH_WINDOW_DAYS window centered on Hannah's Peak_TWL_Time is applied
#   to the historical storm's active date range to determine a match.
#
#   Importantly, the historical catalog contains ONLY named storms (tropical
#   and extratropical). Hannah's TWL-exceedance method captures ALL storms,
#   including winter nor'easters, which are the dominant morphological forcing
#   at Hatteras. Unmatched events in Hannah's record are expected and not
#   indicative of error — they represent nor'easters and unnamed coastal storms.
#
#   Figure comparison:
#     Panel 1 — Annual storm count bar chart (Hannah's comparison), bars colored
#                by fraction of matched/unmatched events.
#     Panel 2 — Timeline scatter: Hannah's storms plotted by year vs Rhigh,
#                colored by match status. Historical named storms overlaid as
#                labeled markers.
#     Panel 3 — Named storm capture table: for each historical named storm,
#                shows whether it was captured, which event matched it, the
#                timing raw_offset, and the matched Rhigh.
#
# Usage:
#   Set READABLE_CSV, BEGIN_YEAR, END_YEAR at the top of USER CONFIGURATION,
#   then run. Works for both 1984–2004 and 2004–2024 periods.
#
# Author: Hannah Henry
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from pathlib import Path
from datetime import timedelta

# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Path to the readable CSV from HAT_create_storms.py
READABLE_CSV = Path(r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\storms\hindcast_storms\base_storms\storms_2004_2024_readable.csv")

# Period to validate — controls which slice of HISTORICAL_STORMS is used
BEGIN_YEAR = 2004
END_YEAR   = 2024

# Matching window (days). A named storm is considered matched if Hannah's
# Peak_TWL_Time falls within [named_start - window, named_end + window].
# The Duck gauge position, storm tracking uncertainty, and the fact that HURDAT2
# reports full storm lifetime (not Hatteras influence period) all justify a
# generous window. Recommended: 3–5 days.
MATCH_WINDOW_DAYS = 5

# Output
OUTPUT_DIR  = Path(r"/scripts/input_prep/storm_creation_final/storm_check/comparison/storm_check_2004_2024.png")
SAVE_FIGURE = True
SHOW_FIGURE = True

# =============================================================================
# HISTORICAL STORM CATALOG
# All named tropical and extratropical storms passing within 60 nautical miles
# of Hatteras Island, NC. Source: NOAA HURDAT2 / Hurricane Research Division.
# Date ranges are the storm's full HURDAT2 active period — not Hatteras-only.
#
# Category codes: H5/H4/H3/H2/H1 = hurricane, TS = tropical storm, ET = extratropical
# =============================================================================

HISTORICAL_STORMS = [
    # -------------------------------------------------------------------------
    # 1984–2004
    # -------------------------------------------------------------------------
    {"name": "DIANA",    "start": "1984-09-08", "end": "1984-09-16", "cat": "H4"},
    {"name": "GLORIA",   "start": "1985-09-16", "end": "1985-10-02", "cat": "H4"},
    {"name": "KATE",     "start": "1985-11-15", "end": "1985-11-23", "cat": "H3"},
    {"name": "CHARLEY",  "start": "1986-08-13", "end": "1986-08-30", "cat": "H1"},
    {"name": "ALBERTO",  "start": "1988-08-05", "end": "1988-08-08", "cat": "TS"},
    {"name": "BOB",      "start": "1991-08-16", "end": "1991-08-29", "cat": "H3"},
    {"name": "DANIELLE", "start": "1992-09-22", "end": "1992-09-26", "cat": "TS"},
    {"name": "EMILY",    "start": "1993-08-22", "end": "1993-09-06", "cat": "H3"},
    {"name": "ALLISON",  "start": "1995-06-03", "end": "1995-06-11", "cat": "H1"},
    {"name": "ARTHUR",   "start": "1996-06-17", "end": "1996-06-23", "cat": "ET"},
    {"name": "JOSEPHINE","start": "1996-10-04", "end": "1996-10-16", "cat": "TS"},
    {"name": "BONNIE",   "start": "1998-08-19", "end": "1998-08-31", "cat": "H3"},
    {"name": "DENNIS",   "start": "1999-08-24", "end": "1999-09-08", "cat": "H2"},
    {"name": "IRENE",    "start": "1999-10-12", "end": "1999-10-19", "cat": "H2"},
    {"name": "ARTHUR",   "start": "2002-07-14", "end": "2002-07-19", "cat": "TS"},
    {"name": "GUSTAV",   "start": "2002-09-08", "end": "2002-09-15", "cat": "H2"},
    {"name": "KYLE",     "start": "2002-09-20", "end": "2002-10-12", "cat": "H1"},
    {"name": "ISABEL",   "start": "2003-09-06", "end": "2003-09-20", "cat": "H5"},
    # -------------------------------------------------------------------------
    # 2004–2024
    # -------------------------------------------------------------------------
    {"name": "ALEX",      "start": "2004-07-31", "end": "2004-08-06", "cat": "H3"},
    {"name": "BONNIE",    "start": "2004-08-03", "end": "2004-08-14", "cat": "TS"},
    {"name": "OPHELIA",   "start": "2005-09-06", "end": "2005-09-23", "cat": "H1"},
    {"name": "ALBERTO",   "start": "2006-06-10", "end": "2006-06-19", "cat": "TS"},
    {"name": "BARRY",     "start": "2007-05-31", "end": "2007-06-05", "cat": "TS"},
    {"name": "GABRIELLE", "start": "2007-09-08", "end": "2007-09-11", "cat": "TS"},
    {"name": "CRISTOBAL", "start": "2008-07-19", "end": "2008-07-23", "cat": "TS"},
    {"name": "IRENE",     "start": "2011-08-21", "end": "2011-08-30", "cat": "H3"},
    {"name": "BERYL",     "start": "2012-05-25", "end": "2012-06-02", "cat": "TS"},
    {"name": "ARTHUR",    "start": "2014-06-28", "end": "2014-07-09", "cat": "H2"},
    {"name": "CLAUDETTE", "start": "2015-07-12", "end": "2015-07-15", "cat": "TS"},
    {"name": "BONNIE",    "start": "2016-05-27", "end": "2016-06-09", "cat": "TS"},
    {"name": "COLIN",     "start": "2016-06-05", "end": "2016-06-08", "cat": "ET"},
    {"name": "HERMINE",   "start": "2016-08-28", "end": "2016-09-08", "cat": "H1"},
    {"name": "JULIA",     "start": "2016-09-13", "end": "2016-09-21", "cat": "TS"},
    {"name": "MATTHEW",   "start": "2016-09-28", "end": "2016-10-10", "cat": "H5"},
    {"name": "UNNAMED",   "start": "2017-08-27", "end": "2017-08-29", "cat": "ET"},
    {"name": "DORIAN",    "start": "2019-08-24", "end": "2019-09-09", "cat": "H5"},
    {"name": "ARTHUR",    "start": "2020-05-16", "end": "2020-05-21", "cat": "ET"},
    {"name": "FAY",       "start": "2020-07-05", "end": "2020-07-11", "cat": "TS"},
    {"name": "CLAUDETTE", "start": "2021-06-17", "end": "2021-06-23", "cat": "TS"},
]

# =============================================================================
# COLOUR SCHEME
# =============================================================================

CAT_COLORS = {
    "H5": "#7b0000",
    "H4": "#b22222",
    "H3": "#e05c1a",
    "H2": "#e8a020",
    "H1": "#f0d040",
    "TS": "#5590d0",
    "ET": "#8888aa",
}

CAT_ORDER = ["H5", "H4", "H3", "H2", "H1", "TS", "ET"]

MATCHED_COLOR   = "#2ca02c"   # green — Hannah event matched a named storm
UNMATCHED_COLOR = "#b0b0b0"   # grey  — Hannah event is a nor'easter or unnamed event
MISSED_COLOR    = "#d62728"   # red   — named storm not captured by Hannah's method


# =============================================================================
# STEP 1 — LOAD AND PREPARE DATA
# =============================================================================

def load_and_filter(csv_path: Path, begin_year: int, end_year: int) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    for col in ['Storm_Start', 'Storm_End', 'Peak_TWL_Time']:
        df[col] = pd.to_datetime(df[col])
    df = df[(df['Calendar_Year'] >= begin_year) & (df['Calendar_Year'] <= end_year)]
    df = df.reset_index(drop=True)
    print(f"Loaded {len(df)} storm events for {begin_year}–{end_year}")
    return df


def filter_historical(begin_year: int, end_year: int) -> list:
    storms = []
    for s in HISTORICAL_STORMS:
        start = pd.Timestamp(s['start'])
        end   = pd.Timestamp(s['end'])
        if start.year >= begin_year and start.year <= end_year:
            storms.append({**s, 'start_ts': start, 'end_ts': end})
    print(f"Historical named storms in period: {len(storms)}")
    return storms


# =============================================================================
# STEP 2 — MATCH EVENTS
# =============================================================================

def match_storms(hannah_df: pd.DataFrame, historical: list,
                 window_days: int) -> tuple:
    """
    For each historical named storm, find the best-matching Hannah event
    within the time window. Best match = highest Rhigh within the window.

    Returns:
        hannah_df    : original DataFrame with added 'matched_to' and 'match_type' columns
        match_results: list of dicts with match details per historical storm
    """
    window = timedelta(days=window_days)
    hannah_df = hannah_df.copy()
    hannah_df['matched_to']    = None    # name of historical storm matched
    hannah_df['match_type']    = 'unmatched'
    hannah_df['timing_offset_h'] = np.nan

    match_results = []

    for storm in historical:
        search_start = storm['start_ts'] - window
        search_end   = storm['end_ts']   + window

        candidates = hannah_df[
            (hannah_df['Peak_TWL_Time'] >= search_start) &
            (hannah_df['Peak_TWL_Time'] <= search_end)
        ]

        label = f"{storm['name']} {storm['start_ts'].year}"

        if len(candidates) == 0:
            match_results.append({
                'label':         label,
                'name':          storm['name'],
                'year':          storm['start_ts'].year,
                'cat':           storm['cat'],
                'start':         storm['start_ts'],
                'end':           storm['end_ts'],
                'matched':       False,
                'matched_idx':   None,
                'matched_rhigh': np.nan,
                'timing_offset_h': np.nan,
            })
        else:
            # Best match = highest Rhigh within window
            best_idx  = candidates['Rhigh_m'].idxmax()
            best_row  = hannah_df.loc[best_idx]
            # Timing raw_offset: hours between peak TWL and nearest edge of storm window
            storm_mid = storm['start_ts'] + (storm['end_ts'] - storm['start_ts']) / 2
            offset_h  = (best_row['Peak_TWL_Time'] - storm_mid).total_seconds() / 3600

            # Only mark matched if not already claimed by a closer storm
            if hannah_df.loc[best_idx, 'match_type'] == 'unmatched':
                hannah_df.loc[best_idx, 'matched_to']      = label
                hannah_df.loc[best_idx, 'match_type']      = 'matched'
                hannah_df.loc[best_idx, 'timing_offset_h'] = offset_h

            match_results.append({
                'label':           label,
                'name':            storm['name'],
                'year':            storm['start_ts'].year,
                'cat':             storm['cat'],
                'start':           storm['start_ts'],
                'end':             storm['end_ts'],
                'matched':         True,
                'matched_idx':     best_idx,
                'matched_rhigh':   best_row['Rhigh_m'],
                'timing_offset_h': offset_h,
            })

    return hannah_df, match_results


# =============================================================================
# STEP 3 — PRINT SUMMARY
# =============================================================================

def print_summary(hannah_df: pd.DataFrame, match_results: list,
                  window_days: int):
    matched   = [r for r in match_results if r['matched']]
    missed    = [r for r in match_results if not r['matched']]
    n_total   = len(match_results)
    n_matched = len(matched)
    n_missed  = len(missed)
    n_noreasters = (hannah_df['match_type'] == 'unmatched').sum()

    print(f"\n{'='*70}")
    print(f"STORM VALIDATION SUMMARY  (window = ±{window_days} days)")
    print(f"{'='*70}")
    print(f"Named storms in period:       {n_total}")
    print(f"  Captured (matched):         {n_matched}  ({100*n_matched/n_total:.0f}%)")
    print(f"  Missed (not captured):      {n_missed}   ({100*n_missed/n_total:.0f}%)")
    print(f"Hannah events total:          {len(hannah_df)}")
    print(f"  Matched to named storm:     {(hannah_df['match_type']=='matched').sum()}")
    print(f"  Unmatched (nor'easters):    {n_noreasters}")

    print(f"\n{'='*70}")
    print("MATCHED NAMED STORMS")
    print(f"{'='*70}")
    for r in sorted(matched, key=lambda x: x['start']):
        offset_str = f"{r['timing_offset_h']:+.0f} h from storm midpoint"
        print(f"  ✓ {r['label']:<18} ({r['cat']})  →  "
              f"Rhigh={r['matched_rhigh']:.3f} m  |  {offset_str}")

    print(f"\n{'='*70}")
    print(f"MISSED NAMED STORMS  (no Hannah event within ±{window_days}-day window)")
    print(f"{'='*70}")
    if missed:
        for r in sorted(missed, key=lambda x: x['start']):
            print(f"  ✗ {r['label']:<18} ({r['cat']})  "
                  f"{r['start'].strftime('%b %d')}–{r['end'].strftime('%b %d, %Y')}")
        print()
        print("  Likely reasons for missing:")
        print("  • Storm's wind field was too far offshore → Duck TWL stayed below threshold")
        print("  • Tropical storm made landfall south of Hatteras → surge attenuated")
        print("    before reaching Duck (consider increasing SURGE_MULTIPLIER)")
        print("  • Short-duration event that didn't exceed threshold in hourly data")
    else:
        print("  All named storms captured!")

    print(f"\n{'='*70}")
    print(f"UNMATCHED HANNAH EVENTS (likely nor'easters / unnamed coastal storms)")
    print(f"{'='*70}")
    noreasters = hannah_df[hannah_df['match_type'] == 'unmatched'].sort_values('Peak_TWL_Time')
    for _, row in noreasters.iterrows():
        season = 'Winter/Spring' if row['Peak_TWL_Time'].month in [11,12,1,2,3,4] else 'Summer/Fall'
        print(f"  {row['Peak_TWL_Time'].strftime('%Y-%m-%d')}  "
              f"Rhigh={row['Rhigh_m']:.3f} m  Hs={row['Peak_Hs_m']:.2f} m  [{season}]")


# =============================================================================
# STEP 4 — FIGURE
# =============================================================================

def make_validation_figure(hannah_df: pd.DataFrame, match_results: list,
                            begin_year: int, end_year: int,
                            window_days: int, save_dir: Path = None):
    """
    Three-panel validation figure.

    Panel 1 — Annual storm counts, stacked bars (matched / unmatched)
    Panel 2 — Rhigh scatter: Hannah events + named storm markers on timeline
    Panel 3 — Named storm capture table
    """
    all_years = list(range(begin_year, end_year + 1))
    matched_df   = hannah_df[hannah_df['match_type'] == 'matched']
    unmatched_df = hannah_df[hannah_df['match_type'] == 'unmatched']

    missed_storms   = [r for r in match_results if not r['matched']]
    captured_storms = [r for r in match_results if r['matched']]

    fig = plt.figure(figsize=(20, 15))
    gs  = fig.add_gridspec(3, 1, height_ratios=[1.2, 2.5, 1.8], hspace=0.38)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2])

    fig.suptitle(
        f"Storm Validation — Hatteras Island ({begin_year}–{end_year})\n"
        f"Duck, NC gauge TWL vs NOAA Named Storm Catalog (±{window_days}-day match window)",
        fontsize=14, fontweight='bold', y=0.98
    )

    # ---- PANEL 1: Annual storm counts ----
    matched_counts   = matched_df.groupby('Calendar_Year').size().reindex(all_years, fill_value=0)
    unmatched_counts = unmatched_df.groupby('Calendar_Year').size().reindex(all_years, fill_value=0)

    ax1.bar(all_years, matched_counts.values,
            color=MATCHED_COLOR, alpha=0.85, label='Matched to named storm',
            edgecolor='white', lw=0.5)
    ax1.bar(all_years, unmatched_counts.values,
            bottom=matched_counts.values,
            color=UNMATCHED_COLOR, alpha=0.75, label='Unmatched (nor\'easter / unnamed)',
            edgecolor='white', lw=0.5)

    # Mark years with missed named storms
    missed_years = set(r['year'] for r in missed_storms)
    for yr in missed_years:
        total = matched_counts[yr] + unmatched_counts[yr]
        ax1.annotate('✗', xy=(yr, total + 0.05),
                     ha='center', va='bottom', fontsize=11,
                     color=MISSED_COLOR, fontweight='bold')

    ax1.set_ylabel('Storms / year', fontsize=10)
    ax1.set_ylim(0, max((matched_counts + unmatched_counts).max() + 2, 6))
    ax1.set_xticks(all_years)
    ax1.set_xticklabels(all_years, rotation=45, ha='right', fontsize=8)
    ax1.grid(axis='y', alpha=0.3)
    ax1.legend(fontsize=9, loc='upper right')
    ax1.set_title('Annual storm counts  (red ✗ = named storm not captured by TWL method)',
                  fontsize=10, pad=4)

    # ---- PANEL 2: Rhigh timeline ----
    # Hannah unmatched events (nor'easters) — small grey circles
    ax2.scatter(unmatched_df['Calendar_Year'] +
                (unmatched_df['Peak_TWL_Time'].dt.dayofyear / 365),
                unmatched_df['Rhigh_m'],
                s=unmatched_df['Duration_h'] * 2,
                color=UNMATCHED_COLOR, alpha=0.55, edgecolors='none', zorder=2,
                label="Nor'easter / unmatched")

    # Hannah matched events — green circles, sized by duration
    ax2.scatter(matched_df['Calendar_Year'] +
                (matched_df['Peak_TWL_Time'].dt.dayofyear / 365),
                matched_df['Rhigh_m'],
                s=matched_df['Duration_h'] * 2,
                color=MATCHED_COLOR, alpha=0.85, edgecolors='k', linewidths=0.5,
                zorder=4, label='Matched (Hannah event)')

    # Named storms — coloured diamonds at estimated mid-storm timing
    for r in match_results:
        storm_mid_frac = r['start'].dayofyear / 365
        x_pos = r['year'] + storm_mid_frac
        cat   = r['cat']
        color = CAT_COLORS.get(cat, '#888888')

        if r['matched']:
            ax2.scatter(x_pos, r['matched_rhigh'],
                        marker='D', s=90, color=color,
                        edgecolors='k', linewidths=0.8, zorder=5)
            # Short name label above diamond
            short = r['name'][:4] + f"\n{r['year']}"
            ax2.annotate(short, xy=(x_pos, r['matched_rhigh']),
                         xytext=(0, 8), textcoords='raw_offset points',
                         ha='center', fontsize=6.5, color='#222222',
                         fontweight='bold')
        else:
            # Missed storms shown as red X on the berm line
            ax2.scatter(x_pos, 1.7,
                        marker='x', s=100, color=MISSED_COLOR,
                        linewidths=2, zorder=6)
            short = r['name'][:4] + f"\n{r['year']}"
            ax2.annotate(short, xy=(x_pos, 1.7),
                         xytext=(0, -16), textcoords='raw_offset points',
                         ha='center', fontsize=6.5, color=MISSED_COLOR,
                         fontweight='bold')

    ax2.axhline(1.7, color='steelblue', ls='--', lw=1.2,
                label='Berm crest (1.7 m)', zorder=1)
    ax2.axhline(1.5, color='orange', ls=':', lw=1.0,
                label='Storm threshold (1.5 m)', zorder=1)

    ax2.set_ylabel('Rhigh (m NAVD88)', fontsize=10)
    ax2.set_xlabel('Year', fontsize=10)
    ax2.set_xlim(begin_year - 0.3, end_year + 1.3)
    ax2.set_xticks(all_years)
    ax2.set_xticklabels(all_years, rotation=45, ha='right', fontsize=8)
    ax2.grid(alpha=0.25)
    ax2.set_title('Storm timeline: Hannah events (circles) + named storm catalog (diamonds = matched, ✗ = missed)',
                  fontsize=10, pad=4)

    # Category legend for diamonds
    cat_handles = [
        mlines.Line2D([], [], marker='D', color='w',
                      markerfacecolor=CAT_COLORS[c], markeredgecolor='k',
                      markersize=8, label=c)
        for c in CAT_ORDER if any(r['cat'] == c for r in match_results)
    ]
    other_handles = [
        mlines.Line2D([], [], marker='o', color='w',
                      markerfacecolor=MATCHED_COLOR, markeredgecolor='k',
                      markersize=8, label='Hannah matched'),
        mlines.Line2D([], [], marker='o', color='w',
                      markerfacecolor=UNMATCHED_COLOR,
                      markersize=8, label="Hannah nor'easter"),
        mlines.Line2D([], [], marker='x', color=MISSED_COLOR,
                      markersize=9, linewidth=2, label='Named storm missed'),
        mlines.Line2D([], [], linestyle='--', color='steelblue',
                      label='Berm crest'),
    ]
    ax2.legend(handles=cat_handles + other_handles, fontsize=8,
               loc='upper left', ncol=3, framealpha=0.85)

    # ---- PANEL 3: Named storm capture table ----
    ax3.axis('off')

    n_cols  = 5
    headers = ['Storm', 'Category', 'Date Range', 'Status', 'Matched Rhigh (m)']
    col_x   = [0.01, 0.18, 0.32, 0.62, 0.82]

    # Header row
    for hdr, x in zip(headers, col_x):
        ax3.text(x, 0.97, hdr, transform=ax3.transAxes,
                 fontsize=9, fontweight='bold', va='top', color='#222222')
    ax3.plot([0, 1], [0.93, 0.93], color='#cccccc', lw=1, transform=ax3.transAxes)

    all_results = sorted(match_results, key=lambda r: r['start'])
    row_h = 0.88 / max(len(all_results), 1)

    for i, r in enumerate(all_results):
        y     = 0.92 - i * row_h
        cat   = r['cat']
        color = CAT_COLORS.get(cat, '#888888')
        bg    = MATCHED_COLOR if r['matched'] else MISSED_COLOR
        status_text = (f"✓ Captured  |  raw_offset {r['timing_offset_h']:+.0f} h"
                       if r['matched'] else "✗ Not captured")

        label_yr   = f"{r['name']} {r['year']}"
        date_range = (f"{r['start'].strftime('%b %d')}–"
                      f"{r['end'].strftime('%b %d')}")
        rhigh_str  = (f"{r['matched_rhigh']:.3f}"
                      if r['matched'] else "—")

        ax3.text(col_x[0], y, label_yr,   transform=ax3.transAxes, fontsize=8, va='top')
        ax3.text(col_x[1], y, cat,        transform=ax3.transAxes, fontsize=8, va='top',
                 color=color, fontweight='bold')
        ax3.text(col_x[2], y, date_range, transform=ax3.transAxes, fontsize=8, va='top',
                 color='#444444')
        ax3.text(col_x[3], y, status_text, transform=ax3.transAxes, fontsize=8, va='top',
                 color=bg, fontweight='bold')
        ax3.text(col_x[4], y, rhigh_str,  transform=ax3.transAxes, fontsize=8, va='top',
                 color='#222222')

        if i % 2 == 0:
            ax3.axhspan(y - row_h * 0.15, y + row_h * 0.85,
                        xmin=0, xmax=1,
                        color='#f5f5f5', transform=ax3.transAxes, zorder=0)

    n_cap = len(captured_storms)
    n_tot = len(all_results)
    ax3.set_title(
        f"Named storm capture rate: {n_cap}/{n_tot} ({100*n_cap/n_tot:.0f}%)  |  "
        f"Window = ±{window_days} days  |  "
        f"Unmatched Hannah events (nor'easters): {(hannah_df['match_type']=='unmatched').sum()}",
        fontsize=10, pad=8
    )

    plt.tight_layout(rect=[0, 0, 1, 0.97])

    if save_dir and SAVE_FIGURE:
        save_dir.mkdir(parents=True, exist_ok=True)
        out_path = save_dir / f"storm_validation_{begin_year}_{end_year}.png"
        fig.savefig(out_path, dpi=250, bbox_inches='tight')
        print(f"\n✓ Saved figure: {out_path}")

    if SHOW_FIGURE:
        plt.show()
    else:
        plt.close()


# =============================================================================
# STEP 5 — SAVE MATCH RESULTS CSV
# =============================================================================

def save_match_table(match_results: list, hannah_df: pd.DataFrame,
                     begin_year: int, end_year: int, save_dir: Path):
    """Save a CSV table of match results for use in the dissertation."""
    rows = []
    for r in sorted(match_results, key=lambda x: x['start']):
        matched_event = None
        if r['matched'] and r['matched_idx'] is not None:
            matched_event = hannah_df.loc[r['matched_idx']]

        rows.append({
            'Storm_Name'       : r['name'],
            'Year'             : r['year'],
            'Category'         : r['cat'],
            'HURDAT_Start'     : r['start'].strftime('%Y-%m-%d'),
            'HURDAT_End'       : r['end'].strftime('%Y-%m-%d'),
            'Captured'         : r['matched'],
            'Matched_Event_Peak': matched_event['Peak_TWL_Time'].strftime('%Y-%m-%d %H:%M') if matched_event is not None else '',
            'Timing_Offset_h'  : round(r['timing_offset_h'], 1) if r['matched'] else '',
            'Matched_Rhigh_m'  : round(r['matched_rhigh'], 3) if r['matched'] else '',
            'Matched_Duration_h': matched_event['Duration_h'] if matched_event is not None else '',
        })

    out_df   = pd.DataFrame(rows)
    csv_path = save_dir / f"storm_validation_{begin_year}_{end_year}.csv"
    save_dir.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(csv_path, index=False)
    print(f"✓ Saved match table: {csv_path}")

    # Also save an annotated Hannah CSV with match status
    annotated_path = save_dir / f"storms_{begin_year}_{end_year}_annotated.csv"
    hannah_df.to_csv(annotated_path, index=False)
    print(f"✓ Saved annotated Hannah CSV: {annotated_path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("\n" + "="*70)
    print("HAT_validate_storms.py — Historical Storm Record Validation")
    print(f"Period:  {BEGIN_YEAR}–{END_YEAR}  |  Match window: ±{MATCH_WINDOW_DAYS} days")
    print(f"CSV:     {READABLE_CSV.name}")
    print("="*70)

    # 1. Load
    hannah_df  = load_and_filter(READABLE_CSV, BEGIN_YEAR, END_YEAR)
    historical = filter_historical(BEGIN_YEAR, END_YEAR)

    # 2. Match
    hannah_df, match_results = match_storms(hannah_df, historical, MATCH_WINDOW_DAYS)

    # 3. Summary
    print_summary(hannah_df, match_results, MATCH_WINDOW_DAYS)

    # 4. Figure
    make_validation_figure(
        hannah_df, match_results,
        BEGIN_YEAR, END_YEAR, MATCH_WINDOW_DAYS,
        save_dir=OUTPUT_DIR
    )

    # 5. Save tables
    if SAVE_FIGURE:
        save_match_table(match_results, hannah_df, BEGIN_YEAR, END_YEAR, OUTPUT_DIR)


if __name__ == "__main__":
    main()
