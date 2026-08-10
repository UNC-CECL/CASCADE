# =============================================================================
# HAT_storm_panel.py
# Simple square storm record panel
# -----------------------------------------------------------------------------
# Single-panel square figure: Rhigh by year, nor'easters as grey background,
# matched named storms as large coloured diamonds with name labels.
# Designed to sit alongside another panel on a 16:9 slide.
#
# Author: Hannah Henry
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from pathlib import Path
from datetime import timedelta

# =============================================================================
# USER CONFIGURATION
# =============================================================================

READABLE_CSV = Path(r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\storms\hindcast_storms\base_storms\storms_1984_2004_base_readable.csv")
BEGIN_YEAR        = 1984
END_YEAR          = 2004
MATCH_WINDOW_DAYS = 5
BERM_CREST_M      = 1.7

OUTPUT_PATH = Path(
    r"/scripts/input_prep/storm_creation_final/from_Hannah/storm_creation/storm_figures/comparison/storm_panel_1984_2004.png")
SAVE_FIGURE = True
SHOW_FIGURE = True

# =============================================================================
# HISTORICAL CATALOG  (1984–2004)
# =============================================================================

HISTORICAL_STORMS = [
    {"name": "Diana",     "start": "1984-09-08", "end": "1984-09-16", "cat": "H4"},
    {"name": "Gloria",    "start": "1985-09-16", "end": "1985-10-02", "cat": "H4"},
    {"name": "Kate",      "start": "1985-11-15", "end": "1985-11-23", "cat": "H3"},
    {"name": "Charley",   "start": "1986-08-13", "end": "1986-08-30", "cat": "H1"},
    {"name": "Alberto",   "start": "1988-08-05", "end": "1988-08-08", "cat": "TS"},
    {"name": "Bob",       "start": "1991-08-16", "end": "1991-08-29", "cat": "H3"},
    {"name": "Danielle",  "start": "1992-09-22", "end": "1992-09-26", "cat": "TS"},
    {"name": "Emily",     "start": "1993-08-22", "end": "1993-09-06", "cat": "H3"},
    {"name": "Allison",   "start": "1995-06-03", "end": "1995-06-11", "cat": "H1"},
    {"name": "Arthur",    "start": "1996-06-17", "end": "1996-06-23", "cat": "ET"},
    {"name": "Josephine", "start": "1996-10-04", "end": "1996-10-16", "cat": "TS"},
    {"name": "Bonnie",    "start": "1998-08-19", "end": "1998-08-31", "cat": "H3"},
    {"name": "Dennis",    "start": "1999-08-24", "end": "1999-09-08", "cat": "H2"},
    {"name": "Irene",     "start": "1999-10-12", "end": "1999-10-19", "cat": "H2"},
    {"name": "Arthur",    "start": "2002-07-14", "end": "2002-07-19", "cat": "TS"},
    {"name": "Gustav",    "start": "2002-09-08", "end": "2002-09-15", "cat": "H2"},
    {"name": "Kyle",      "start": "2002-09-20", "end": "2002-10-12", "cat": "H1"},
    {"name": "Isabel",    "start": "2003-09-06", "end": "2003-09-20", "cat": "H5"},
    {"name": "Alex",      "start": "2004-07-31", "end": "2004-08-06", "cat": "H3"},
    {"name": "Bonnie",    "start": "2004-08-03", "end": "2004-08-14", "cat": "TS"},
]

CAT_COLORS = {
    "H5": "#6b0000",
    "H4": "#c0392b",
    "H3": "#e67e22",
    "H2": "#f39c12",
    "H1": "#f1c40f",
    "TS": "#2980b9",
    "ET": "#7f8c8d",
}

CAT_LABELS = {
    "H5": "Cat. 5", "H4": "Cat. 4", "H3": "Cat. 3",
    "H2": "Cat. 2", "H1": "Cat. 1", "TS": "Trop. Storm", "ET": "Extratropical",
}

# =============================================================================
# HELPERS
# =============================================================================

def load_and_match(csv_path, begin_year, end_year, window_days):
    df = pd.read_csv(csv_path)
    for col in ['Storm_Start', 'Storm_End', 'Peak_TWL_Time']:
        df[col] = pd.to_datetime(df[col])
    df = df[(df['Calendar_Year'] >= begin_year) &
            (df['Calendar_Year'] <= end_year)].reset_index(drop=True)

    window = timedelta(days=window_days)
    df['match_type'] = 'unmatched'
    df['match_cat']  = None
    df['match_name'] = None
    match_results = []

    for s in HISTORICAL_STORMS:
        start_ts = pd.Timestamp(s['start'])
        end_ts   = pd.Timestamp(s['end'])
        if not (begin_year <= start_ts.year <= end_year):
            continue
        cands = df[(df['Peak_TWL_Time'] >= start_ts - window) &
                   (df['Peak_TWL_Time'] <= end_ts   + window)]
        if len(cands) > 0:
            best = cands['Rhigh_m'].idxmax()
            if df.loc[best, 'match_type'] == 'unmatched':
                df.loc[best, 'match_type'] = 'matched'
                df.loc[best, 'match_cat']  = s['cat']
                df.loc[best, 'match_name'] = s['name']
            match_results.append({**s, 'start_ts': start_ts, 'end_ts': end_ts,
                                   'matched_idx': best,
                                   'matched_rhigh': df.loc[best, 'Rhigh_m']})
    return df, match_results


# =============================================================================
# FIGURE
# =============================================================================

def make_panel(df, match_results, begin_year, end_year):

    unmatched = df[df['match_type'] == 'unmatched']
    all_years = list(range(begin_year, end_year + 1))
    cats_present = set(df[df['match_type'] == 'matched']['match_cat'].dropna())

    # ── square canvas ─────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 7), facecolor='white')
    ax.set_facecolor('white')

    # ── hurricane season shading ───────────────────────────────────────────────
    for yr in all_years:
        ax.axvspan(yr + 6/12, yr + 11/12,
                   color='#fff4e6', alpha=0.7, zorder=0, lw=0)

    # ── berm crest ────────────────────────────────────────────────────────────
    ax.axhline(BERM_CREST_M, color='#9ab8cc', lw=1.1, ls='--',
               zorder=1, alpha=0.9)
    ax.text(end_year + 0.75, BERM_CREST_M, 'Berm\ncrest',
            va='center', ha='left', fontsize=7.5, color='#7a9db8',
            linespacing=1.3)

    # ── nor'easters — grey circles ────────────────────────────────────────────
    x_u = (unmatched['Calendar_Year'] +
            unmatched['Peak_TWL_Time'].dt.dayofyear / 365)
    ax.scatter(x_u, unmatched['Rhigh_m'],
               s=unmatched['Peak_Hs_m'] * 5,
               color='#ced4da', alpha=0.5, edgecolors='none', zorder=2)

    # ── matched named storms — diamonds ───────────────────────────────────────
    # Pre-compute per-year counts for horizontal label spreading
    yr_groups = {}
    for r in match_results:
        yr = r['start_ts'].year
        yr_groups.setdefault(yr, []).append(r)

    for r in match_results:
        idx   = r['matched_idx']
        row   = df.loc[idx]
        cat   = r['cat']
        color = CAT_COLORS[cat]
        name  = r['name']
        yr    = r['start_ts'].year

        x_pos = row['Calendar_Year'] + row['Peak_TWL_Time'].dayofyear / 365
        rhigh = row['Rhigh_m']

        # Diamond
        ax.scatter(x_pos, rhigh,
                   marker='D', s=160, color=color,
                   edgecolors='white', linewidths=0.9,
                   zorder=5, alpha=0.96)

        # Label horizontal spread for same-year storms
        siblings = yr_groups[yr]
        n        = len(siblings)
        rank     = next(i for i, rr in enumerate(siblings)
                        if rr['name'] == name and
                        rr['start'] == r['start'])
        x_nudge  = (rank - (n - 1) / 2) * 20 if n > 1 else 0

        # Isabel: bold coloured label above
        if name == 'Isabel':
            ax.annotate(name, xy=(x_pos, rhigh),
                        xytext=(x_nudge, 10), textcoords='raw_offset points',
                        ha='center', fontsize=9, fontweight='bold',
                        color=color,
                        arrowprops=dict(arrowstyle='-', color=color,
                                        lw=0.8, shrinkA=0, shrinkB=3))
        else:
            y_nudge = 9 if rhigh < 2.8 else -13
            ax.annotate(name, xy=(x_pos, rhigh),
                        xytext=(x_nudge, y_nudge), textcoords='raw_offset points',
                        ha='center', fontsize=7.5, color='#2c3e50',
                        arrowprops=dict(arrowstyle='-', color='#bbbbbb',
                                        lw=0.5, shrinkA=0, shrinkB=3))

    # ── axes ──────────────────────────────────────────────────────────────────
    ax.set_xlim(begin_year - 0.4, end_year + 1.4)
    ax.set_ylim(1.3, df['Rhigh_m'].max() * 1.07)

    # Show only every other year to keep x-axis uncluttered
    tick_years = [y for y in all_years if y % 2 == 0]
    ax.set_xticks(tick_years)
    ax.set_xticklabels([str(y) for y in tick_years],
                       rotation=45, ha='right', fontsize=9, color='#333344')
    ax.set_yticks([2.0, 2.5, 3.0, 3.5, 4.0, 4.5])
    ax.set_yticklabels(['2.0', '2.5', '3.0', '3.5', '4.0', '4.5'],
                       fontsize=9, color='#333344')

    ax.set_ylabel('Total Water Level — R$_{high}$  (m NAVD88)',
                  fontsize=10, color='#333344', labelpad=6)
    ax.set_xlabel('Year', fontsize=10, color='#333344', labelpad=6)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#ddddee')
    ax.spines['bottom'].set_color('#ddddee')
    ax.yaxis.grid(True, color='#eeeeee', lw=0.8, zorder=0)
    ax.set_axisbelow(True)

    # ── title ─────────────────────────────────────────────────────────────────
    ax.set_title(f'Storm Record  ·  Hatteras Island  ·  {begin_year}–{end_year}',
                 fontsize=11, fontweight='bold', color='#1a1a2e', pad=10)

    # ── compact legend ────────────────────────────────────────────────────────
    cat_order = ["H5", "H4", "H3", "H2", "H1", "TS", "ET"]
    cat_handles = [
        mlines.Line2D([], [], marker='D', color='w',
                      markerfacecolor=CAT_COLORS[c],
                      markeredgecolor='white', markersize=8,
                      label=CAT_LABELS[c])
        for c in cat_order if c in cats_present
    ]
    norea = mlines.Line2D([], [], marker='o', color='w',
                          markerfacecolor='#ced4da', markersize=7,
                          label="Nor'easter / unnamed")

    ax.legend(handles=cat_handles + [norea],
              fontsize=7.5, loc='upper left',
              framealpha=0.92, edgecolor='#ddddee',
              handlelength=1.2, labelspacing=0.4,
              title='Named storms\n(matched to record)',
              title_fontsize=7.5)

    plt.tight_layout()

    if SAVE_FIGURE:
        OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(OUTPUT_PATH, dpi=300, bbox_inches='tight', facecolor='white')
        print(f'✓ Saved: {OUTPUT_PATH}')

    if SHOW_FIGURE:
        plt.show()
    else:
        plt.close()


# =============================================================================
# MAIN
# =============================================================================

def main():
    df, match_results = load_and_match(
        READABLE_CSV, BEGIN_YEAR, END_YEAR, MATCH_WINDOW_DAYS
    )
    n = sum(1 for r in match_results)
    print(f'{len(df)} events loaded  |  {n} named storms matched')
    make_panel(df, match_results, BEGIN_YEAR, END_YEAR)


if __name__ == '__main__':
    main()
