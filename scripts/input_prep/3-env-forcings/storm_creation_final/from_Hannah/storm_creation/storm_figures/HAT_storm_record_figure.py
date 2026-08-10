# =============================================================================
# HAT_storm_record_figure.py
# Presentation-quality Storm Record Figure — Hatteras Island
# -----------------------------------------------------------------------------
# Produces a clean, two-panel figure suitable for PowerPoint / dissertation
# figures. Emphasises the matched named storm record without calling out misses.
#
# Panel 1 (top) — Compact annual bar chart: storms/year, bars split into
#   summer/fall (hurricane season, Jun–Nov) and winter/spring (nor'easter
#   season, Dec–May), with matched named storms annotated above bars.
#
# Panel 2 (bottom) — Main timeline: all storm events as small grey circles
#   (the full nor'easter + tropical record), named storm matches as large
#   coloured diamonds labelled by name, Isabel called out with a special
#   annotation. Berm crest shown as a clean reference line.
#
# Usage:
#   1. Point READABLE_CSV to your _readable.csv comparison from HAT_create_storms.py
#   2. Set BEGIN_YEAR / END_YEAR to match your hindcast period
#   3. Run — figure is displayed and optionally saved to OUTPUT_PATH
#
# Author: Hannah Henry
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from pathlib import Path
from datetime import timedelta

# =============================================================================
# USER CONFIGURATION
# =============================================================================

READABLE_CSV = Path(r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\storms\hindcast_storms\fixed_storms\storms_2004_2024_readable.csv")

BEGIN_YEAR = 2004
END_YEAR   = 2024

MATCH_WINDOW_DAYS = 5   # must match HAT_validate_storms.py setting
BERM_CREST_M      = 1.7 # m NAVD88

OUTPUT_PATH = Path(
    r"/scripts/input_prep/storm_creation_final/from_Hannah/storm_creation/storm_figures/comparison/storm_record_2004_2024_fixed.png")
SAVE_FIGURE = True
SHOW_FIGURE = True

# =============================================================================
# HISTORICAL STORM CATALOG  (1984–2024 full record)
# All named tropical and extratropical storms passing within 60 nautical miles
# of Hatteras Island, NC. Source: NOAA HURDAT2 / Hurricane Research Division.
# Date ranges are the storm's full HURDAT2 active period — not Hatteras-only.
# match_storms() filters to BEGIN_YEAR–END_YEAR automatically.
#
# Category codes: H5/H4/H3/H2/H1 = hurricane, TS = tropical storm, ET = extratropical
# =============================================================================

HISTORICAL_STORMS = [
    # -------------------------------------------------------------------------
    # 1984–2004
    # -------------------------------------------------------------------------
    {"name": "Diana",    "start": "1984-09-08", "end": "1984-09-16", "cat": "H4"},
    {"name": "Gloria",   "start": "1985-09-16", "end": "1985-10-02", "cat": "H4"},
    {"name": "Kate",     "start": "1985-11-15", "end": "1985-11-23", "cat": "H3"},
    {"name": "Charley",  "start": "1986-08-13", "end": "1986-08-30", "cat": "H1"},
    {"name": "Alberto",  "start": "1988-08-05", "end": "1988-08-08", "cat": "TS"},
    {"name": "Bob",      "start": "1991-08-16", "end": "1991-08-29", "cat": "H3"},
    {"name": "Danielle", "start": "1992-09-22", "end": "1992-09-26", "cat": "TS"},
    {"name": "Emily",    "start": "1993-08-22", "end": "1993-09-06", "cat": "H3"},
    {"name": "Allison",  "start": "1995-06-03", "end": "1995-06-11", "cat": "H1"},
    {"name": "Arthur",   "start": "1996-06-17", "end": "1996-06-23", "cat": "ET"},
    {"name": "Josephine","start": "1996-10-04", "end": "1996-10-16", "cat": "TS"},
    {"name": "Bonnie",   "start": "1998-08-19", "end": "1998-08-31", "cat": "H3"},
    {"name": "Dennis",   "start": "1999-08-24", "end": "1999-09-08", "cat": "H2"},
    {"name": "Irene",    "start": "1999-10-12", "end": "1999-10-19", "cat": "H2"},
    {"name": "Arthur",   "start": "2002-07-14", "end": "2002-07-19", "cat": "TS"},
    {"name": "Gustav",   "start": "2002-09-08", "end": "2002-09-15", "cat": "H2"},
    {"name": "Kyle",     "start": "2002-09-20", "end": "2002-10-12", "cat": "H1"},
    {"name": "Isabel",   "start": "2003-09-06", "end": "2003-09-20", "cat": "H5"},
    # -------------------------------------------------------------------------
    # 2004–2024
    # -------------------------------------------------------------------------
    {"name": "Alex",      "start": "2004-07-31", "end": "2004-08-06", "cat": "H3"},
    {"name": "Bonnie",    "start": "2004-08-03", "end": "2004-08-14", "cat": "TS"},
    {"name": "Ophelia",   "start": "2005-09-06", "end": "2005-09-23", "cat": "H1"},
    {"name": "Alberto",   "start": "2006-06-10", "end": "2006-06-19", "cat": "TS"},
    {"name": "Barry",     "start": "2007-05-31", "end": "2007-06-05", "cat": "TS"},
    {"name": "Gabrielle", "start": "2007-09-08", "end": "2007-09-11", "cat": "TS"},
    {"name": "Cristobal", "start": "2008-07-19", "end": "2008-07-23", "cat": "TS"},
    {"name": "Irene",     "start": "2011-08-21", "end": "2011-08-30", "cat": "H3"},
    {"name": "Beryl",     "start": "2012-05-25", "end": "2012-06-02", "cat": "TS"},
    {"name": "Arthur",    "start": "2014-06-28", "end": "2014-07-09", "cat": "H2"},
    {"name": "Claudette", "start": "2015-07-12", "end": "2015-07-15", "cat": "TS"},
    {"name": "Bonnie",    "start": "2016-05-27", "end": "2016-06-09", "cat": "TS"},
    {"name": "Colin",     "start": "2016-06-05", "end": "2016-06-08", "cat": "ET"},
    {"name": "Hermine",   "start": "2016-08-28", "end": "2016-09-08", "cat": "H1"},
    {"name": "Julia",     "start": "2016-09-13", "end": "2016-09-21", "cat": "TS"},
    {"name": "Matthew",   "start": "2016-09-28", "end": "2016-10-10", "cat": "H5"},
    {"name": "Unnamed",   "start": "2017-08-27", "end": "2017-08-29", "cat": "ET"},
    {"name": "Dorian",    "start": "2019-08-24", "end": "2019-09-09", "cat": "H5"},
    {"name": "Arthur",    "start": "2020-05-16", "end": "2020-05-21", "cat": "ET"},
    {"name": "Fay",       "start": "2020-07-05", "end": "2020-07-11", "cat": "TS"},
    {"name": "Claudette", "start": "2021-06-17", "end": "2021-06-23", "cat": "TS"},
]

# =============================================================================
# SAFFIR-SIMPSON COLOUR PALETTE  — muted, print-safe, colourblind-considerate
# =============================================================================

CAT_COLORS = {
    "H5": "#6b0000",   # deep burgundy
    "H4": "#c0392b",   # strong red
    "H3": "#e67e22",   # burnt orange
    "H2": "#f39c12",   # warm amber
    "H1": "#f1c40f",   # golden yellow
    "TS": "#2980b9",   # steel blue
    "ET": "#7f8c8d",   # slate grey
}

CAT_LABELS = {
    "H5": "Category 5",
    "H4": "Category 4",
    "H3": "Category 3",
    "H2": "Category 2",
    "H1": "Category 1",
    "TS": "Tropical Storm",
    "ET": "Extratropical",
}

# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def load_data(csv_path, begin_year, end_year):
    df = pd.read_csv(csv_path)
    for col in ['Storm_Start', 'Storm_End', 'Peak_TWL_Time']:
        df[col] = pd.to_datetime(df[col])
    df = df[(df['Calendar_Year'] >= begin_year) &
            (df['Calendar_Year'] <= end_year)].reset_index(drop=True)
    df['month']  = df['Peak_TWL_Time'].dt.month
    df['season'] = df['month'].apply(
        lambda m: 'hurricane' if m in [6, 7, 8, 9, 10, 11] else 'winter'
    )
    return df


def match_storms(df, window_days, begin_year=None, end_year=None):
    window = timedelta(days=window_days)
    df = df.copy()
    df['matched_to']  = None
    df['match_type']  = 'unmatched'
    df['match_cat']   = None
    match_results     = []

    for s in HISTORICAL_STORMS:
        start_ts = pd.Timestamp(s['start'])
        end_ts   = pd.Timestamp(s['end'])
        # Skip storms outside the active hindcast period
        if start_ts.year < begin_year or start_ts.year > end_year:
            continue
        cands    = df[(df['Peak_TWL_Time'] >= start_ts - window) &
                      (df['Peak_TWL_Time'] <= end_ts   + window)]
        label    = f"{s['name']} {start_ts.year}"

        if len(cands) == 0:
            match_results.append({'label': label, 'matched': False, **s,
                                  'start_ts': start_ts, 'end_ts': end_ts,
                                  'matched_rhigh': np.nan, 'matched_idx': None})
        else:
            best_idx = cands['Rhigh_m'].idxmax()
            if df.loc[best_idx, 'match_type'] == 'unmatched':
                df.loc[best_idx, 'matched_to'] = label
                df.loc[best_idx, 'match_type'] = 'matched'
                df.loc[best_idx, 'match_cat']  = s['cat']
            match_results.append({'label': label, 'matched': True, **s,
                                  'start_ts': start_ts, 'end_ts': end_ts,
                                  'matched_rhigh': df.loc[best_idx, 'Rhigh_m'],
                                  'matched_idx': best_idx})
    return df, match_results


# Bar chart season colours — defined at module level so legend patches can reference them
WINTER_C = '#a8c4e0'   # soft steel blue  (Dec–May nor'easter season)
HURR_C   = '#f4a261'   # warm sand orange (Jun–Nov hurricane season)

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE
# ─────────────────────────────────────────────────────────────────────────────

def make_figure(df, match_results, begin_year, end_year):

    matched_df   = df[df['match_type'] == 'matched']
    unmatched_df = df[df['match_type'] == 'unmatched']
    all_years    = list(range(begin_year, end_year + 1))

    # ── per-year season split ─────────────────────────────────────────────────
    hurr_counts  = (df[df['season'] == 'hurricane']
                    .groupby('Calendar_Year').size()
                    .reindex(all_years, fill_value=0))
    wint_counts  = (df[df['season'] == 'winter']
                    .groupby('Calendar_Year').size()
                    .reindex(all_years, fill_value=0))

    # ── layout ────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(16, 8.5), facecolor='white')
    gs  = fig.add_gridspec(2, 1, height_ratios=[1, 2.8],
                           hspace=0.08, top=0.91, bottom=0.09,
                           left=0.065, right=0.97)
    ax1 = fig.add_subplot(gs[0])   # bar chart
    ax2 = fig.add_subplot(gs[1])   # main timeline

    # ── figure title ──────────────────────────────────────────────────────────
    fig.text(0.5, 0.965,
             f'Hatteras Island Storm Record  ·  {begin_year}–{end_year}',
             ha='center', va='top', fontsize=15, fontweight='bold',
             color='#1a1a2e', fontfamily='DejaVu Sans')
    fig.text(0.5, 0.937,
             'Duck, NC tide gauge  ·  WIS hindcast waves  ·  Stockdon (2006) runup  ·  '
             f'TWL threshold {BERM_CREST_M} m (berm crest)',
             ha='center', va='top', fontsize=8.5, color='#555577',
             fontfamily='DejaVu Sans')

    # ══════════════════════════════════════════════════════════════════════════
    # PANEL 1 — Annual bar chart
    # ══════════════════════════════════════════════════════════════════════════

    ax1.bar(all_years, wint_counts.values,
            color=WINTER_C, width=0.72, label='Winter/Spring (nor\'easters)',
            zorder=2)
    ax1.bar(all_years, hurr_counts.values,
            bottom=wint_counts.values,
            color=HURR_C, width=0.72, label='Hurricane Season (Jun–Nov)',
            zorder=2)

    # Annotate matched named storms above bars with storm name
    match_yr_labels = {}
    for r in match_results:
        if r['matched']:
            yr = r['start_ts'].year
            name = r['name']
            if yr not in match_yr_labels:
                match_yr_labels[yr] = []
            match_yr_labels[yr].append(name)

    for yr, names in match_yr_labels.items():
        total = hurr_counts[yr] + wint_counts[yr]
        label_str = ', '.join(names)
        ax1.text(yr, total + 0.15, label_str,
                 ha='center', va='bottom', fontsize=6.8,
                 color='#2c3e50', style='italic',
                 fontfamily='DejaVu Sans')

    ax1.set_xlim(begin_year - 0.6, end_year + 0.6)
    ax1.set_ylim(0, max((hurr_counts + wint_counts).max() + 1.8, 6))
    ax1.set_xticks([])
    ax1.set_ylabel('Storms / yr', fontsize=9, color='#444455')
    ax1.tick_params(axis='y', labelsize=8, colors='#444455')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_color('#ccccdd')
    ax1.yaxis.grid(True, color='#eeeeee', lw=0.8, zorder=0)
    ax1.set_axisbelow(True)

    # Bar legend removed — season colours explained in main panel legend below

    # ══════════════════════════════════════════════════════════════════════════
    # PANEL 2 — Rhigh timeline
    # ══════════════════════════════════════════════════════════════════════════

    # Subtle hurricane season shading every year
    for yr in all_years:
        ax2.axvspan(yr + 6/12, yr + 11/12,
                    color='#fff3e0', alpha=0.55, zorder=0, lw=0)

    # Berm crest reference
    ax2.axhline(BERM_CREST_M, color='#8dafc8', lw=1.2, ls='--',
                zorder=1, alpha=0.85)
    ax2.text(end_year + 0.45, BERM_CREST_M,
             f'Berm crest\n{BERM_CREST_M} m', va='center',
             fontsize=7.5, color='#6899b8', fontfamily='DejaVu Sans')

    # All storms — grey circles, sized lightly by Hs
    x_all = unmatched_df['Calendar_Year'] + \
            unmatched_df['Peak_TWL_Time'].dt.dayofyear / 365
    ax2.scatter(x_all, unmatched_df['Rhigh_m'],
                s=unmatched_df['Peak_Hs_m'] * 7,
                color='#c8ccd8', alpha=0.55, edgecolors='none',
                zorder=2, label='Nor\'easter / unnamed event')

    # ── Named storm matches — coloured diamonds ───────────────────────────────
    # Collect categories actually present for legend
    cats_present = set()

    for r in match_results:
        if not r['matched']:
            continue

        idx      = r['matched_idx']
        row      = df.loc[idx]
        cat      = r['cat']
        color    = CAT_COLORS[cat]
        cats_present.add(cat)

        x_pos = row['Calendar_Year'] + row['Peak_TWL_Time'].dayofyear / 365
        rhigh = row['Rhigh_m']

        # Diamond marker
        ax2.scatter(x_pos, rhigh,
                    marker='D', s=220, color=color,
                    edgecolors='white', linewidths=1.2,
                    zorder=5, alpha=0.95)

        # Label — name only (year is implicit from x position)
        name = r['name']
        yr   = r['start_ts'].year

        # Count how many matched storms share the same calendar year
        # so we can raw_offset labels horizontally to avoid overlap
        same_yr_matched = [rr for rr in match_results
                           if rr['matched'] and rr['start_ts'].year == yr]
        n_same = len(same_yr_matched)
        yr_rank = [i for i, rr in enumerate(same_yr_matched)
                   if rr['name'] == r['name'] and
                   rr['start_ts'].strftime('%Y-%m-%d') == r['start_ts'].strftime('%Y-%m-%d')]
        yr_rank = yr_rank[0] if yr_rank else 0

        # Horizontal nudge: spread labels when multiple storms share a year
        if n_same > 1:
            x_nudge = (yr_rank - (n_same - 1) / 2) * 22  # points
        else:
            x_nudge = 0

        # Offset labels to avoid overlap; Isabel gets special treatment
        if name == 'Isabel':
            ax2.annotate(
                name,
                xy=(x_pos, rhigh),
                xytext=(x_nudge, 11), textcoords='raw_offset points',
                ha='center', fontsize=8.5, fontweight='bold',
                color=color, fontfamily='DejaVu Sans',
                arrowprops=dict(arrowstyle='-', color=color,
                                lw=0.8, shrinkA=0, shrinkB=3)
            )
        else:
            y_offset = 10 if rhigh < 2.8 else -14
            ax2.annotate(
                name,
                xy=(x_pos, rhigh),
                xytext=(x_nudge, y_offset), textcoords='raw_offset points',
                ha='center', fontsize=7.8, color='#2c3e50',
                fontfamily='DejaVu Sans',
                arrowprops=dict(arrowstyle='-', color='#aaaaaa',
                                lw=0.5, shrinkA=0, shrinkB=3)
            )

    # ── axes cosmetics ────────────────────────────────────────────────────────
    ax2.set_xlim(begin_year - 0.6, end_year + 0.85)
    ax2.set_ylim(1.2, df['Rhigh_m'].max() * 1.08)
    ax2.set_xticks(all_years)
    ax2.set_xticklabels(all_years, rotation=45, ha='right',
                        fontsize=8.5, color='#333344')
    ax2.set_ylabel('Total Water Level — Rhigh  (m NAVD88)',
                   fontsize=9.5, color='#444455', labelpad=8)
    ax2.set_xlabel('Year', fontsize=9.5, color='#444455', labelpad=8)
    ax2.tick_params(axis='y', labelsize=8.5, colors='#444455')

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_color('#ccccdd')
    ax2.spines['bottom'].set_color('#ccccdd')
    ax2.yaxis.grid(True, color='#eeeeee', lw=0.8, zorder=0)
    ax2.set_axisbelow(True)

    # ── legend — all items consolidated here ─────────────────────────────────
    cat_order = ["H5", "H4", "H3", "H2", "H1", "TS", "ET"]
    cat_handles = [
        mlines.Line2D([], [], marker='D', color='w',
                      markerfacecolor=CAT_COLORS[c],
                      markeredgecolor='white',
                      markersize=9,
                      label=CAT_LABELS[c])
        for c in cat_order if c in cats_present
    ]
    norea_handle = mlines.Line2D([], [], marker='o', color='w',
                                 markerfacecolor='#c8ccd8',
                                 markersize=7,
                                 label="Nor'easter / unnamed event")
    hurr_patch  = mpatches.Patch(facecolor=HURR_C,  alpha=0.85,
                                 label='Hurricane season bars (Jun–Nov)')
    wint_patch  = mpatches.Patch(facecolor=WINTER_C, alpha=0.85,
                                 label='Winter / Spring bars (Dec–May)')

    ax2.legend(
        handles=cat_handles + [norea_handle, hurr_patch, wint_patch],
        fontsize=8, loc='upper left',
        framealpha=0.92, edgecolor='#ddddee',
        ncol=1, handlelength=1.5,
        title='Storm record legend',
        title_fontsize=8.5,
    )

    # ── hurricane season shading label (one annotation, not per year) ─────────
    ax2.text(begin_year + 0.5, 1.24,
             'Hurricane season →', fontsize=7, color='#c8965a',
             fontstyle='italic', va='bottom')

    # ─────────────────────────────────────────────────────────────────────────
    # Save / show
    # ─────────────────────────────────────────────────────────────────────────
    if SAVE_FIGURE:
        OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(OUTPUT_PATH, dpi=300, bbox_inches='tight',
                    facecolor='white')
        print(f'✓ Saved: {OUTPUT_PATH}')

    if SHOW_FIGURE:
        plt.show()
    else:
        plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    df             = load_data(READABLE_CSV, BEGIN_YEAR, END_YEAR)
    df, match_res  = match_storms(df, MATCH_WINDOW_DAYS, BEGIN_YEAR, END_YEAR)

    n_matched = sum(1 for r in match_res if r['matched'])
    n_total   = len(match_res)
    print(f'Loaded {len(df)} events  |  Named storm matches: {n_matched}/{n_total}')

    make_figure(df, match_res, BEGIN_YEAR, END_YEAR)


if __name__ == '__main__':
    main()
