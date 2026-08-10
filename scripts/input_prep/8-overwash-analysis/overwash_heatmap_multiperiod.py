"""
Hatteras Island Overwash Heatmap — Multi-Period
================================================

PERIOD SELECTION  (change PERIOD at the top of the config block):
    PERIOD = (1984, 2004)  →  Period 1 only
    PERIOD = (2004, 2024)  →  Period 2 only
    PERIOD = (1984, 2024)  →  Combined (both periods)

DISPLAY MODE:
    MODE = 'observed'     →  Observed imagery heatmap only (default)
    MODE = 'comparison'   →  Observed (top panel) + Modeled CASCADE (bottom panel)

MODELED DATA:
    Set MODELED_FILE to a .npy path when CASCADE results are available.
    If MODELED_FILE = None, the script falls back to 'observed' mode
    automatically (no crash, just a warning).

    Expected format for MODELED_FILE:
        np.save('cascade_overwash.npy', array)
        where array.shape == (n_display_rows, 90), float dtype.
        Values: 1.0 = overwash, 0.0 = no overwash, np.nan = gap year / no data.
        Rows must align 1-to-1 with the all_rows list printed at runtime.
        Run this script once in 'observed' mode to see the row order printed
        in the console; use that to build your modeled matrix from CASCADE comparison.

MULTI-OBS YEARS:
    Years with two imagery dates (2009, 2014, 2023) appear as two rows
    labeled with season abbreviations (e.g., '2009 Spr' / '2009 Fall').
    The storm timeline correctly routes each storm to the matching season row.

STORM TIMELINE TOGGLE:
    SHOW_STORM_TIMELINE = True / False
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import warnings
warnings.filterwarnings('ignore')

# ══════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit here
# ══════════════════════════════════════════════════════════════════════
FILE_PATH = r'C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\overwash_analysis\Hatteras_Overwash_Data.xlsx'
SHEET     = 'Overwash_Matrix'
DPI       = 200

# --- Period ---
PERIOD = (1984, 2004)          # (1984,2004) | (2004,2024) | (1984,2024)

# --- Display mode ---
MODE = 'observed'              # 'observed' | 'comparison'

# --- Modeled data (comparison mode only) ---
# Path to .npy binary matrix (n_display_rows × 90).  Set None = not yet available.
MODELED_FILE = None            # e.g. r'...\cascade_overwash_period1.npy'

# --- Storm timeline ---
SHOW_STORM_TIMELINE = False    # True = show timeline panel to the right of heatmap

# ══════════════════════════════════════════════════════════════════════
# OUTPUT FILE — auto-selected from PERIOD and MODE
# ══════════════════════════════════════════════════════════════════════
_BASE  = r'C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\overwash_analysis'
_PTAG  = {(1984, 2004): 'period1', (2004, 2024): 'period2', (1984, 2024): 'combined'}
_MTAG  = {'observed': 'obs', 'comparison': 'vs_model'}
OUTPUT_FILE = rf'{_BASE}\overwash_heatmap_{_PTAG[PERIOD]}_{_MTAG[MODE]}.png'

# ══════════════════════════════════════════════════════════════════════
# ISLAND SECTIONS  (matches ANN_TOWN_SPANS in CASCADE code)
# ══════════════════════════════════════════════════════════════════════
SECTIONS = [
    ("Cape Point",                              1,   6,  "inter"),
    ("Buxton",                                  7,   8,  "village"),
    ("Inter-village:\nBuxton–Avon",             9,  20,  "inter"),
    ("Avon",                                   21,  31,  "village"),
    ("Inter-village:\nAvon–Tri-Village\n(Wimble Shoals)", 32, 67, "inter"),
    ("Tri-Village",                            68,  83,  "village"),
    ("Pea Island /\nN. Rodanthe",              84,  90,  "inter"),
]

# ══════════════════════════════════════════════════════════════════════
# COLORS
# ══════════════════════════════════════════════════════════════════════
CLR_ZERO        = '#FFFFFF'   # assessed — no overwash
CLR_OBS         = '#C0392B'   # overwash observed  (red)
CLR_MOD         = '#2471A3'   # overwash modeled   (blue — visually distinct from observed)
CLR_HATCH_FACE  = '#F2F2F2'   # no imagery / gap background
CLR_HATCH_EDGE  = '#BBBBBB'   # hatch line color
CLR_BAR_VILLAGE = '#C4A882'   # section bar — village
CLR_BAR_INTER   = '#7FA8C4'   # section bar — inter-village
CLR_PERIOD_DIV  = '#34495E'   # period boundary dashed line (combined plot only)

SEASON_ABBR = {'Winter': 'Win', 'Spring': 'Spr', 'Summer': 'Sum', 'Fall': 'Fall'}

# ══════════════════════════════════════════════════════════════════════
# POOR-QUALITY / PARTIAL-COVERAGE OBSERVATIONS
# Set of (year, season) tuples → y-axis label gets a '*' marker
# ══════════════════════════════════════════════════════════════════════
POOR_QUALITY_OBS = {
    (1996, 'Fall'),     # poor image quality
    (2006, 'Summer'),   # medium-to-poor image quality
    (2023, 'Fall'),     # coverage stops at ~domain 17
    (2024, 'Spring'),   # coverage starts at ~domain 17
}

# ══════════════════════════════════════════════════════════════════════
# STORM DATA
# ──────────────────────────────────────────────────────────────────────
# Format per entry:  (display_name, category, imagery_matched, season_hint)
#
#   imagery_matched: True  → storm predates imagery date (effects are visible)
#                   False → storm postdates imagery, or no imagery that year
#
#   season_hint:    None  → single-obs year (no ambiguity needed)
#                   str   → for multi-obs years, which seasonal observation
#                           is associated with this storm entry
#                           (e.g., 'Fall' routes Ida to the 2009 Fall row)
# ══════════════════════════════════════════════════════════════════════
STORMS_P1 = {
    1984: [('DIANA',                'H4',   True,  None)],
    1985: [('GLORIA',               'H4',   True,  None),
           ('KATE',                 'H3',   True,  None)],
    1986: [('CHARLEY',              'H1',   False, None)],
    1988: [('ALBERTO',              'TS',   False, None)],
    1991: [('Perfect Storm',        'NE 5', True,  None),
           ('BOB',                  'H3',   True,  None)],
    1992: [("Dec '92 Nor'easter",   'NE 4', True,  None)],
    1993: [('Storm of the Century', 'NE 5', True,  None),
           ('EMILY',                'H3',   True,  None)],
    1996: [("Jan '96 Nor'easter",   'NE 4', True,  None)],
    1998: [('BONNIE',               'H3',   True,  None)],
    1999: [('DENNIS',               'H2',   False, None),
           ('IRENE',                'H2',   False, None)],
    2002: [('KYLE',                 'H1',   False, None)],
    2003: [('ISABEL',               'H5',   False, None)],
    # Spring 2004 imagery (May 25): captures Sep 2003 Isabel damage
    # Alex (Jul 2004) post-dates the imagery
    2004: [('ISABEL (effects)',      'H5',   True,  None),
           ('ALEX',                  'H3',   False, None)],
}

STORMS_P2 = {
    # Spring 2004 imagery (May 25): captures Sep 2003 Isabel damage
    2004: [('ISABEL (effects)',       'H5',   True,  None),
           ('ALEX',                   'H3',   False, None)],
    # Fall 2005 imagery (Sep 1): Ophelia hit Sep 14 — imagery predates storm
    2005: [('OPHELIA',               'H1',   False, None)],
    # Summer 2008 imagery (Jun 25): Hanna & Ike both Sep 2008 — imagery predates both
    2008: [('HANNA',                 'H1',   False, None),
           ('IKE',                   'H4',   False, None)],
    # 2009 has Spring (May 30) and Fall (Nov 15) imagery
    # Ida made extratropical landfall Nov 12 — captured in Fall obs only
    2009: [('Ida (extratrop.)',       'TS',   True,  'Fall')],
    # Summer 2011 imagery (Aug 27): Irene landfall was Aug 27 — same day, effects captured
    2011: [('IRENE',                 'H1',   True,  None)],
    # Sandy Oct 2012: no imagery that year
    2012: [('SANDY',                 'H3',   False, None)],
    # Spring 2013 imagery (Apr 6): captures Oct 2012 Sandy damage
    2013: [('SANDY (effects)',        'H3',   True,  None)],
    # Matthew Oct 2016: no imagery that year
    2016: [('MATTHEW',               'H1',   False, None)],
    # Summer 2018 imagery (Aug 22): Feb 2018 nor'easter preceded imagery;
    # Florence (Sep 2018) postdates
    2018: [("Feb '18 Nor'easter",    'NE 3', True,  None),
           ('FLORENCE',              'H4',   False, None)],
    # Fall 2019 imagery (Sep 10): Dorian landfall Sep 6 — just before imagery
    2019: [('DORIAN',                'H2',   True,  None)],
    # Winter 2022 imagery (Feb 1): May 2022 nor'easter postdates imagery
    2022: [("May '22 Nor'easter",    'NE 3', False, None)],
    # 2023 has Spring (Apr 24) and Fall (Oct 8) imagery
    # Fall obs captures May 2022 nor'easter effects (first good imagery since)
    2023: [("May '22 Nor'easter",    'NE 3', True,  'Fall')],
    # Spring 2024 imagery (Mar 14): Debby (Aug 2024) postdates imagery
    2024: [('DEBBY',                 'H1',   False, None)],
}

# Observed years with imagery but no catalog storm → gray dot in timeline
# Format: set of (year, season) tuples
NO_EVENT_OBS_P1 = {(1987,'Summer'), (1989,'Summer'), (1995,'Fall'), (1997,'Fall')}
NO_EVENT_OBS_P2 = {
    (2006, 'Summer'), (2009, 'Spring'),
    (2014, 'Winter'), (2014, 'Summer'),
    (2017, 'Winter'), (2022, 'Winter'),
    (2023, 'Spring'), (2024, 'Spring'),
}

# Merge STORMS_P1 + STORMS_P2 for the combined period, avoiding duplicate names
def _merge_storms(p1, p2):
    merged = {}
    all_years = sorted(set(list(p1.keys()) + list(p2.keys())))
    for yr in all_years:
        entries = list(p1.get(yr, []))
        existing_names = {e[0] for e in entries}
        entries += [e for e in p2.get(yr, []) if e[0] not in existing_names]
        merged[yr] = entries
    return merged

STORMS_ALL = _merge_storms(STORMS_P1, STORMS_P2)

# ── Storm category → (marker_size, matched_color, unmatched_color) ──
CAT_STYLE = {
    'H5':   (9,  '#7B0000', '#D4AAAA'),
    'H4':   (8,  '#A00000', '#DDBBBB'),
    'H3':   (7,  '#C0392B', '#E5CCCC'),
    'H2':   (6,  '#C0570B', '#E2CCBB'),
    'H1':   (5,  '#CC7030', '#E8D8CC'),
    'TS':   (4,  '#777777', '#BBBBBB'),
    'NE 5': (9,  '#1A237E', '#AABBD8'),
    'NE 4': (8,  '#2B5BAA', '#BBCCE5'),
    'NE 3': (7,  '#3D7DBF', '#BBCCDD'),
}

def cat_style(cat, matched):
    sz, cm, cu = CAT_STYLE.get(cat, (5, '#555555', '#AAAAAA'))
    return sz, (cm if matched else cu)

# ══════════════════════════════════════════════════════════════════════
# PERIOD-DEPENDENT CONFIG  (auto-selected)
# ══════════════════════════════════════════════════════════════════════
_TITLES = {
    (1984, 2004): 'Hatteras Island — Overwash Observations 1984–2004',
    (2004, 2024): 'Hatteras Island — Overwash Observations 2004–2024',
    (1984, 2024): 'Hatteras Island — Overwash Observations 1984–2024',
}
PLOT_TITLE = _TITLES[PERIOD]

if PERIOD == (1984, 2004):
    STORMS_BY_YEAR = STORMS_P1
    NO_EVENT_OBS   = NO_EVENT_OBS_P1
elif PERIOD == (2004, 2024):
    STORMS_BY_YEAR = STORMS_P2
    NO_EVENT_OBS   = NO_EVENT_OBS_P2
else:
    STORMS_BY_YEAR = STORMS_ALL
    NO_EVENT_OBS   = NO_EVENT_OBS_P1 | NO_EVENT_OBS_P2

SHOW_PERIOD_DIVIDER = (PERIOD == (1984, 2024))

# ══════════════════════════════════════════════════════════════════════
# LOAD OBSERVED DATA
# ══════════════════════════════════════════════════════════════════════
df_raw = pd.read_excel(FILE_PATH, sheet_name=SHEET, skiprows=3, header=0)
df_raw['Year'] = pd.to_numeric(df_raw['Year'], errors='coerce')
df_obs = df_raw[df_raw['Year'].between(*PERIOD)].copy()
df_obs['Year'] = df_obs['Year'].astype(int)

domain_cols = [c for c in df_raw.columns if str(c).replace('.0', '').strip().isdigit()]
domain_ints = [int(float(c)) for c in domain_cols]
n_d = len(domain_ints)

# ══════════════════════════════════════════════════════════════════════
# BUILD DISPLAY ROW LIST
# ──────────────────────────────────────────────────────────────────────
# Each entry: (row_label, year, season_str_or_None, obs_datarow_or_None)
#   Gap years       → obs_datarow = None
#   Multi-obs years → one entry per observation, labeled with season abbrev.
#                     e.g., '2009 Spr' and '2009 Fall' as separate rows
# ══════════════════════════════════════════════════════════════════════
all_rows = []
year_obs_count = df_obs.groupby('Year').size().to_dict()

for yr in range(PERIOD[0], PERIOD[1] + 1):
    if yr not in year_obs_count:
        # Gap year — no imagery
        all_rows.append((str(yr), yr, None, None))
    else:
        yr_data = df_obs[df_obs['Year'] == yr].sort_values('Imagery_Date')
        multi   = (year_obs_count[yr] > 1)
        for _, row in yr_data.iterrows():
            season = str(row.get('Season', ''))
            label  = (f"{yr} {SEASON_ABBR.get(season, season)}" if multi
                      else str(yr))
            all_rows.append((label, yr, season, row))

n_rows = len(all_rows)

# ── Build observed matrix (n_rows × n_d) ──
obs_matrix  = np.full((n_rows, n_d), np.nan)
obs_row_set = set()

for ri, (lbl, yr, season, row) in enumerate(all_rows):
    if row is not None:
        obs_row_set.add(ri)
        obs_matrix[ri, :] = pd.to_numeric(row[domain_cols].values, errors='coerce')

# ── Load modeled matrix if available (comparison mode) ──
mod_matrix = None
if MODE == 'comparison':
    if MODELED_FILE is None:
        print(
            "\nWARNING: MODE='comparison' requested but MODELED_FILE=None.\n"
            "Falling back to MODE='observed'. Set MODELED_FILE to a .npy\n"
            "binary matrix (n_display_rows × 90) to enable comparison.\n"
        )
        MODE = 'observed'
    else:
        mod_matrix = np.load(MODELED_FILE).astype(float)
        if mod_matrix.shape != (n_rows, n_d):
            raise ValueError(
                f"Modeled matrix shape {mod_matrix.shape} does not match "
                f"expected ({n_rows}, {n_d}).\n"
                "Rows must correspond 1-to-1 with all_rows (see console comparison).\n"
                "Use NaN rows for gap years."
            )

# ══════════════════════════════════════════════════════════════════════
# STORM → DISPLAY ROW MAPPING
# ──────────────────────────────────────────────────────────────────────
# row_storms[ri] = [(name, cat, matched), ...]
# Accounts for season_hint in multi-obs years: Ida → 2009 Fall row, not Spring.
# Also handles storm-only gap years (e.g. Sandy 2012 — no imagery that year).
# ══════════════════════════════════════════════════════════════════════
row_storms = {}

for yr, storm_list in STORMS_BY_YEAR.items():
    if not (PERIOD[0] <= yr <= PERIOD[1]):
        continue
    year_row_info = [
        (ri, season)
        for ri, (lbl, yyr, season, row) in enumerate(all_rows)
        if yyr == yr
    ]
    if not year_row_info:
        continue

    for (name, cat, matched, season_hint) in storm_list:
        if season_hint is None:
            # Single-obs year or no season ambiguity → first row for this year
            target_ri = year_row_info[0][0]
        else:
            # Multi-obs year: match by season string
            target_ri = next(
                (ri for ri, s in year_row_info if s == season_hint),
                year_row_info[0][0]  # fallback to first if season not found
            )
        row_storms.setdefault(target_ri, []).append((name, cat, matched))

# ══════════════════════════════════════════════════════════════════════
# FIGURE LAYOUT
# ──────────────────────────────────────────────────────────────────────
# Figure height scales with number of display rows.
# Comparison mode stacks two panels so figure is ~1.7× taller.
# ══════════════════════════════════════════════════════════════════════
fig_h_single = max(9.0, n_rows * 0.40 + 2.0)
fig_h        = fig_h_single * 1.7 if MODE == 'comparison' else fig_h_single

# Fixed-inch margins → fractions of figure height
TOP_IN  = 0.5   # top white space
BOT_IN  = 1.6   # space for legend + section bar + footnote
SECB_IN = 0.45  # section bar height
GAP_IN  = 0.40  # gap between panels (comparison mode)

top_f  = TOP_IN  / fig_h
bot_f  = BOT_IN  / fig_h
secb_f = SECB_IN / fig_h
gap_f  = GAP_IN  / fig_h

if SHOW_STORM_TIMELINE:
    fig_w      = 17.0
    hm_x       = 0.055
    hm_w       = 0.62
    tl_x       = 0.705
    tl_w       = 0.27
    legend_anc = 0.39
else:
    fig_w      = 14.0
    hm_x       = 0.07
    hm_w       = 0.88
    tl_x       = None
    tl_w       = None
    legend_anc = 0.50

fig = plt.figure(figsize=(fig_w, fig_h))

if MODE == 'comparison':
    panel_h = (1.0 - top_f - bot_f - secb_f - gap_f) / 2
    obs_bot = bot_f + secb_f + panel_h + gap_f
    mod_bot = bot_f + secb_f
    ax_obs  = fig.add_axes([hm_x, obs_bot, hm_w, panel_h])
    ax_mod  = fig.add_axes([hm_x, mod_bot, hm_w, panel_h])
    axb     = fig.add_axes([hm_x, bot_f,   hm_w, secb_f])
    axt     = (fig.add_axes([tl_x, obs_bot, tl_w, panel_h])
               if SHOW_STORM_TIMELINE else None)
else:
    hm_h   = 1.0 - top_f - bot_f - secb_f
    hm_bot = bot_f + secb_f
    ax_obs = fig.add_axes([hm_x, hm_bot, hm_w, hm_h])
    ax_mod = None
    axb    = fig.add_axes([hm_x, bot_f,  hm_w, secb_f])
    axt    = (fig.add_axes([tl_x, hm_bot, tl_w, hm_h])
              if SHOW_STORM_TIMELINE else None)

# ══════════════════════════════════════════════════════════════════════
# HEATMAP DRAWING FUNCTION
# ──────────────────────────────────────────────────────────────────────
# Draws one binary heatmap panel.  Called once for 'observed' mode and
# twice (observed + modeled) for 'comparison' mode.
# ══════════════════════════════════════════════════════════════════════
def draw_heatmap(ax, matrix, fill_color, panel_label, show_xlabel):
    """
    ax          : matplotlib Axes to draw on
    matrix      : (n_rows × n_d) float array; NaN = no data / gap year
    fill_color  : color for overwash cells (CLR_OBS or CLR_MOD)
    panel_label : string appended to title, e.g. 'Observed' or 'Modeled (CASCADE)'
                  pass '' for single-panel mode
    show_xlabel : bool — draw x-axis labels and label (False for top panel in comparison)
    """
    ax.set_facecolor(CLR_HATCH_FACE)

    # Determine which rows have any data in this matrix
    data_row_set = {ri for ri in range(n_rows)
                    if not np.all(np.isnan(matrix[ri, :]))}

    # Hatch entire rows with no data (gap years or un-modeled years)
    for ri in range(n_rows):
        if ri not in data_row_set:
            ax.add_patch(mpatches.Rectangle(
                (-0.5, ri - 0.5), n_d, 1,
                facecolor=CLR_HATCH_FACE, edgecolor=CLR_HATCH_EDGE,
                hatch='////', linewidth=0, zorder=1))

    # Hatch individual NaN cells within data rows
    for ri in data_row_set:
        for xi, v in enumerate(matrix[ri, :]):
            if np.isnan(v):
                ax.add_patch(mpatches.Rectangle(
                    (xi - 0.5, ri - 0.5), 1, 1,
                    facecolor=CLR_HATCH_FACE, edgecolor=CLR_HATCH_EDGE,
                    hatch='////', linewidth=0, zorder=4))

    # Data cells
    masked = np.ma.masked_where(np.isnan(matrix), matrix)
    cmap   = mcolors.ListedColormap([CLR_ZERO, fill_color])
    norm   = mcolors.BoundaryNorm([0, 0.5, 1.0], cmap.N)
    cmap.set_bad(color='none')
    ax.imshow(masked, aspect='auto', cmap=cmap, norm=norm,
              interpolation='none',
              extent=[-0.5, n_d - 0.5, n_rows - 0.5, -0.5], zorder=2)

    # Row separators
    for i in range(n_rows + 1):
        ax.axhline(y=i - 0.5, color='white', linewidth=0.35, zorder=3)

    # Section vertical dividers
    for _, lo, _, _ in SECTIONS:
        if lo in domain_ints:
            ax.axvline(x=domain_ints.index(lo) - 0.5,
                       color='white', linewidth=1.2, zorder=3)

    # Period 1 / Period 2 boundary (combined plot only)
    if SHOW_PERIOD_DIVIDER:
        for ri, (lbl, yr, season, row) in enumerate(all_rows):
            if yr == 2005:
                ax.axhline(y=ri - 0.5, color=CLR_PERIOD_DIV,
                           linewidth=2.0, linestyle='--', zorder=5)
                ax.text(n_d - 1, ri - 0.65,
                        '  Period 1  ↑       ↓  Period 2  ',
                        ha='right', va='bottom', fontsize=7.5,
                        color=CLR_PERIOD_DIV, style='italic',
                        fontweight='bold', zorder=6)
                break

    # ── Y-axis ──
    ax.set_yticks(range(n_rows))
    ylabels = []
    for (lbl, yr, season, row) in all_rows:
        pq = (yr, season) in POOR_QUALITY_OBS
        ylabels.append(lbl + ('*' if pq else ''))
    ax.set_yticklabels(ylabels, fontsize=8.5)

    # Bold labels for rows with imagery
    for ri, (lbl, yr, season, row) in enumerate(all_rows):
        if row is not None:
            ax.get_yticklabels()[ri].set_fontweight('bold')

    ylabel_str = 'Year' + (f'  [{panel_label}]' if panel_label else '')
    ax.set_ylabel(ylabel_str, fontsize=11, labelpad=6)

    # ── X-axis ──
    if show_xlabel:
        tick_pos = [i for i, d in enumerate(domain_ints)
                    if d == 1 or d % 10 == 0]
        ax.set_xticks(tick_pos)
        ax.set_xticklabels([str(domain_ints[i]) for i in tick_pos],
                           fontsize=8)
    else:
        ax.set_xticks([])

    # ── Title ──
    title = PLOT_TITLE + (f'  —  {panel_label}' if panel_label else '')
    ax.set_title(title, fontsize=13, fontweight='bold', pad=10)


# ══════════════════════════════════════════════════════════════════════
# DRAW HEATMAP(S)
# ══════════════════════════════════════════════════════════════════════
if MODE == 'comparison':
    draw_heatmap(ax_obs, obs_matrix, CLR_OBS, 'Observed',          show_xlabel=False)
    draw_heatmap(ax_mod, mod_matrix, CLR_MOD, 'Modeled (CASCADE)', show_xlabel=True)
else:
    draw_heatmap(ax_obs, obs_matrix, CLR_OBS, '', show_xlabel=True)


# ══════════════════════════════════════════════════════════════════════
# STORM TIMELINE PANEL  (only when SHOW_STORM_TIMELINE = True)
# ══════════════════════════════════════════════════════════════════════
if SHOW_STORM_TIMELINE and axt is not None:
    axt.set_xlim(0, 1)
    axt.set_ylim(n_rows - 0.5, -0.5)
    axt.set_facecolor('white')
    axt.axis('off')

    axt.text(0.5, -0.95, 'Storm Events', va='center', ha='center',
             fontsize=9, fontweight='bold', color='#333333',
             transform=axt.transData)

    # Subtle year guide lines
    for i, (lbl, yr, season, row) in enumerate(all_rows):
        lc = '#EEEEEE' if row is None else '#E0E0E0'
        axt.axhline(y=i, color=lc, linewidth=0.5, zorder=0)

    AXIS_X   = 0.10
    CONN_END = AXIS_X + 0.10
    LABEL_X  = AXIS_X + 0.13
    axt.axvline(x=AXIS_X, color='#CCCCCC', linewidth=1.2, zorder=1)

    # Build label positions — spread multi-storm rows vertically
    storm_positions = []
    for ri, storms in row_storms.items():
        n = len(storms)
        offsets = ([0.0]         if n == 1 else
                   [-0.28, 0.28] if n == 2 else
                   [-(n - 1) * 0.25 + i * 0.25 for i in range(n)])
        for (name, cat, matched), offset in zip(storms, offsets):
            storm_positions.append((ri, ri + offset, name, cat, matched))
    storm_positions.sort(key=lambda x: x[1])

    # Draw storm dots, connectors, and labels
    for (ri_exact, ri_label, name, cat, matched) in storm_positions:
        sz, color = cat_style(cat, matched)
        alpha = 1.0 if matched else 0.65
        # Dot at exact year row
        axt.plot(AXIS_X, ri_exact, 'o', color=color, markersize=sz, zorder=4,
                 markeredgecolor='white', markeredgewidth=0.5, alpha=alpha)
        # Vertical jog connector if label is raw_offset
        if abs(ri_label - ri_exact) > 0.05:
            axt.plot([AXIS_X, AXIS_X], [ri_exact, ri_label],
                     '-', color=color, linewidth=0.6, alpha=0.35, zorder=2)
        axt.plot([AXIS_X, CONN_END], [ri_label, ri_label],
                 '-', color=color, linewidth=0.7, alpha=0.45, zorder=2)
        axt.text(LABEL_X, ri_label, f'{name}  ({cat})',
                 va='center', ha='left', fontsize=7,
                 fontweight='bold' if matched else 'normal',
                 fontstyle='normal' if matched else 'italic',
                 color=color, alpha=alpha)

    # No-storm obs rows and quiet gap-year rows
    for ri, (lbl, yr, season, row) in enumerate(all_rows):
        if ri in row_storms:
            continue
        if row is not None:
            # Observed row with no catalog storm
            if (yr, season) in NO_EVENT_OBS:
                axt.plot(AXIS_X, ri, 'o', color='#BBBBBB', markersize=4,
                         zorder=3, markeredgecolor='white', markeredgewidth=0.5)
                axt.text(LABEL_X, ri, '(no named storm)', va='center',
                         ha='left', fontsize=6.5, color='#BBBBBB',
                         fontstyle='italic')
        else:
            # Gap year — open circle
            axt.plot(AXIS_X, ri, 'o', color='#DDDDDD', markersize=3.5,
                     zorder=3, markeredgecolor='#CCCCCC', markeredgewidth=0.5,
                     fillstyle='none')

    # Intensity mini-legend
    leg_y = n_rows + 0.2
    axt.text(0.0, leg_y, 'Intensity scale:', va='top', ha='left',
             fontsize=6, color='#555555', fontweight='bold')
    leg_y += 0.7
    for cat_lbl in [('H5', 'H5'), ('H3', 'H3'), ('H1', 'H1'),
                    ('NE 5', 'NE 5'), ('NE 4', 'NE 4'), ('NE 3', 'NE 3'), ('TS', 'TS')]:
        sz, color = cat_style(cat_lbl[0], True)
        axt.plot(0.05, leg_y, 'o', color=color, markersize=sz,
                 markeredgecolor='white', markeredgewidth=0.4)
        axt.text(0.16, leg_y, cat_lbl[1], va='center', ha='left',
                 fontsize=6, color=color)
        leg_y += 0.55
    axt.plot(0.05, leg_y, 'o', color='#DDDDDD', markersize=3.5,
             markeredgecolor='#CCCCCC', markeredgewidth=0.5, fillstyle='none')
    axt.text(0.16, leg_y, 'Quiet / no imagery', va='center', ha='left',
             fontsize=6, color='#AAAAAA')


# ══════════════════════════════════════════════════════════════════════
# SECTION ANNOTATION BAR
# ══════════════════════════════════════════════════════════════════════
axb.set_xlim(-0.5, n_d - 0.5)
axb.set_ylim(0, 1)
axb.axis('off')

for sec_name, lo, hi, stype in SECTIONS:
    lo_idx = domain_ints.index(lo) if lo in domain_ints else 0
    hi_idx = domain_ints.index(hi) if hi in domain_ints else n_d - 1
    mid    = (lo_idx + hi_idx) / 2
    bc = CLR_BAR_VILLAGE if stype == 'village' else CLR_BAR_INTER
    lc = '#3D2B1F'       if stype == 'village' else '#1A3545'
    axb.add_patch(mpatches.FancyBboxPatch(
        (lo_idx - 0.5, 0.52), hi_idx - lo_idx + 1, 0.48,
        boxstyle='square,pad=0', lw=0.8, edgecolor='white', facecolor=bc))
    kw = dict(ha='center', va='top', color=lc, multialignment='center')
    axb.text(mid, 0.38, sec_name,
             fontsize=6 if sec_name == 'Buxton' else (6.2 if '\n' in sec_name else 7.2),
             rotation=90 if sec_name == 'Buxton' else 0, **kw)

# ── X-axis description ──
# Placed as fig.text() below the section bar to avoid overlapping section labels
fig.text(
    hm_x + hm_w / 2,
    bot_f - (0.3 / fig_h),
    'CASCADE Domain   (1 = Cape Point  →  90 = Pea Island / N. Rodanthe)',
    ha='center', va='top', fontsize=10, color='black')
# ══════════════════════════════════════════════════════════════════════
# LEGEND
# ══════════════════════════════════════════════════════════════════════
hatch_patch = mpatches.Patch(
    facecolor=CLR_HATCH_FACE, edgecolor=CLR_HATCH_EDGE,
    hatch='////', label='No imagery / not assessed')

legend_handles = [
    mpatches.Patch(fc=CLR_OBS,  ec='#888', lw=0.5, label='Overwash observed'),
    mpatches.Patch(fc=CLR_ZERO, ec='#888', lw=0.5, label='No overwash (assessed)'),
    hatch_patch,
]
if MODE == 'comparison':
    legend_handles.insert(1, mpatches.Patch(
        fc=CLR_MOD, ec='#888', lw=0.5, label='Overwash modeled (CASCADE)'))

fig.legend(
    handles=legend_handles,
    loc='lower center',
    bbox_to_anchor=(legend_anc, 0.02),
    ncol=len(legend_handles), fontsize=8.5, framealpha=0.95,
    edgecolor='#BBBBBB', columnspacing=1.4,
    handlelength=1.6, handleheight=1.2)


# ══════════════════════════════════════════════════════════════════════
# FOOTNOTE
# ══════════════════════════════════════════════════════════════════════
footnote_lines = [
    'Bold years = imagery available  |  * poor image quality or partial domain coverage  '
    '|  Sources: Hapke & Henderson (2007); Google Earth'
]
if PERIOD in [(1984, 2004), (1984, 2024)]:
    footnote_lines.append(
        'P1 — ALEX (H3, Jul 2004) postdates May 2004 imagery; '
        'ISABEL (H5, Sep 2003) effects visible in Spring 2004 imagery'
    )
if PERIOD in [(2004, 2024), (1984, 2024)]:
    footnote_lines.append(
        'P2 — OPHELIA (2005): imagery Sep 1 predates storm (Sep 14);  '
        "IRENE (2011): imagery Aug 27 = landfall day;  "
        "2023 Fall imagery captures May 2022 nor'easter effects"
    )
if MODE == 'comparison':
    footnote_lines.append(
        'Modeled panel: CASCADE binary overwash '
        '(1 = overwash flux > threshold at domain in model year)'
    )
if SHOW_STORM_TIMELINE:
    footnote_lines.append(
        'Timeline: bold label = storm visible in imagery; '
        'italic = no imagery that year or imagery predates storm  |  '
        'Dot size scaled by peak intensity'
    )

footnote = '\n'.join(footnote_lines)
fig.text(hm_x, 0.003, footnote,
         fontsize=6.5, color='#555555', va='bottom', linespacing=1.5)


# ══════════════════════════════════════════════════════════════════════
# SAVE
# ══════════════════════════════════════════════════════════════════════
plt.savefig(OUTPUT_FILE, dpi=DPI, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {OUTPUT_FILE}")
print(f"  Period: {PERIOD[0]}–{PERIOD[1]}  |  Mode: {MODE}  |  "
      f"Rows: {n_rows} total ({len(obs_row_set)} with imagery)")

# Print row order for reference when building MODELED_FILE
print("\nDisplay row order (use this to align your MODELED_FILE):")
for ri, (lbl, yr, season, row) in enumerate(all_rows):
    flag = 'OBS' if row is not None else 'gap'
    print(f"  Row {ri:2d}: {lbl:<12}  ({flag})")

plt.show()
