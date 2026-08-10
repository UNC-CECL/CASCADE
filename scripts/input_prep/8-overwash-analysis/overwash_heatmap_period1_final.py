"""
Hatteras Island Overwash Heatmap — Period 1 (1984–2004)

TOGGLE:
    SHOW_STORM_TIMELINE = True   → heatmap + storm timeline panel
    SHOW_STORM_TIMELINE = False  → heatmap only (wider)
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
FILE_PATH   = r'C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\overwash_analysis\Hatteras_Overwash_Data.xlsx'
SHEET       = 'Overwash_Matrix'
PERIOD      = (1984, 2004)
OUTPUT_FILE = r'C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\overwash_analysis\overwash_heatmap_period1.png'
DPI         = 200

SHOW_STORM_TIMELINE = False   # ← toggle: True = show timeline, False = heatmap only

# ══════════════════════════════════════════════════════════════════════
# ISLAND SECTIONS (matches ANN_TOWN_SPANS in CASCADE code)
# ══════════════════════════════════════════════════════════════════════
SECTIONS = [
    ("Cape Point",                              1,  6,  "inter"),
    ("Buxton",                                  7,  8,  "village"),
    ("Inter-village:\nBuxton–Avon",             9,  20, "inter"),
    ("Avon",                                   21,  31, "village"),
    ("Inter-village:\nAvon–Tri-Village\n(Wimble Shoals)", 32, 67, "inter"),
    ("Tri-Village",                            68,  83, "village"),
    ("Pea Island /\nN. Rodanthe",              84,  90, "inter"),
]

# ══════════════════════════════════════════════════════════════════════
# COLORS
# ══════════════════════════════════════════════════════════════════════
CLR_ZERO        = '#FFFFFF'   # assessed — no overwash
CLR_ONE         = '#C0392B'   # overwash observed
CLR_HATCH_FACE  = '#F2F2F2'   # no imagery background
CLR_HATCH_EDGE  = '#BBBBBB'   # no imagery hatch lines
CLR_BAR_VILLAGE = '#C4A882'   # section bar — village
CLR_BAR_INTER   = '#7FA8C4'   # section bar — inter-village

POOR_QUALITY_YEARS = {1996}

# ══════════════════════════════════════════════════════════════════════
# STORM DATA
# ══════════════════════════════════════════════════════════════════════
# Format per entry: (display_name, category, imagery_year_match)
#   imagery_year_match = True  → bold, vivid color (storm in imagery year)
#   imagery_year_match = False → italic, muted color (no imagery that year,
#                                or storm occurred after imagery date)
STORMS_BY_YEAR = {
    1984: [('DIANA',                'H4',   True)],
    1985: [('GLORIA',               'H4',   True),
           ('KATE',                 'H3',   True)],
    1986: [('CHARLEY',              'H1',   False)],
    1988: [('ALBERTO',              'TS',   False)],
    1991: [('Perfect Storm',        'NE 5', True),
           ('BOB',                  'H3',   True)],
    1992: [("Dec '92 Nor'easter",   'NE 4', True)],
    1993: [('Storm of the Century', 'NE 5', True),
           ('EMILY',                'H3',   True)],
    1996: [("Jan '96 Nor'easter",   'NE 4', True)],
    1998: [('BONNIE',               'H3',   True)],
    1999: [('DENNIS',               'H2',   False),
           ('IRENE',                'H2',   False)],
    2002: [('KYLE',                 'H1',   False)],
    2003: [('ISABEL',               'H5',   False)],
    2004: [('ALEX',                 'H3',   False)],   # post-dates May 2004 imagery
}
# Imagery years with no named storm in catalog
NO_EVENT_OBS_YEARS = {1987, 1989, 1995, 1997}

# Storm category → (marker_size, color_matched, color_unmatched)
CAT_STYLE = {
    'H5':   (9,  '#7B0000', '#D4AAAA'),
    'H4':   (8,  '#A00000', '#DDBBBB'),
    'H3':   (7,  '#C0392B', '#E5CCCC'),
    'H2':   (6,  '#C0570B', '#E2CCBB'),
    'H1':   (5,  '#CC7030', '#E8D8CC'),
    'TS':   (4,  '#777777', '#BBBBBB'),
    'NE 5': (9,  '#1A237E', '#AABBD8'),
    'NE 4': (8,  '#2B5BAA', '#BBCCE5'),
}
def cat_style(cat, matched):
    sz, cm, cu = CAT_STYLE.get(cat, (5, '#555555', '#AAAAAA'))
    return sz, (cm if matched else cu)

# ══════════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════════
df_raw = pd.read_excel(FILE_PATH, sheet_name=SHEET, skiprows=3, header=0)
df_raw['Year'] = pd.to_numeric(df_raw['Year'], errors='coerce')
df_obs = df_raw[df_raw['Year'].between(*PERIOD)].copy()
df_obs['Year'] = df_obs['Year'].astype(int)

domain_cols = [c for c in df_raw.columns if str(c).replace('.0','').strip().isdigit()]
domain_ints = [int(float(c)) for c in domain_cols]
n_d = len(domain_ints)

all_years = list(range(PERIOD[0], PERIOD[1]+1))
n_y       = len(all_years)
obs_years = set(df_obs['Year'].tolist())

orig_matrix = np.full((n_y, n_d), np.nan)
for yi, yr in enumerate(all_years):
    row = df_obs[df_obs['Year'] == yr]
    if not row.empty:
        orig_matrix[yi, :] = pd.to_numeric(row.iloc[0][domain_cols].values, errors='coerce')

# ══════════════════════════════════════════════════════════════════════
# FIGURE LAYOUT — adapts to toggle
# ══════════════════════════════════════════════════════════════════════
if SHOW_STORM_TIMELINE:
    fig_w        = 17.0
    heatmap_x    = 0.055
    heatmap_w    = 0.62
    timeline_x   = 0.705
    timeline_w   = 0.27
    legend_anchor = 0.39
else:
    fig_w        = 14.0
    heatmap_x    = 0.07
    heatmap_w    = 0.88
    timeline_x   = None
    timeline_w   = None
    legend_anchor = 0.50

fig = plt.figure(figsize=(fig_w, 10))
ax  = fig.add_axes([heatmap_x, 0.26, heatmap_w, 0.64])
axb = fig.add_axes([heatmap_x, 0.18, heatmap_w, 0.06])
axt = fig.add_axes([timeline_x, 0.26, timeline_w, 0.64]) if SHOW_STORM_TIMELINE else None

# ══════════════════════════════════════════════════════════════════════
# DRAW HEATMAP
# ══════════════════════════════════════════════════════════════════════
ax.set_facecolor(CLR_HATCH_FACE)

# Hatched rows for unobserved years
for yi, yr in enumerate(all_years):
    if yr not in obs_years:
        ax.add_patch(mpatches.Rectangle(
            (-0.5, yi-0.5), n_d, 1,
            facecolor=CLR_HATCH_FACE, edgecolor=CLR_HATCH_EDGE,
            hatch='////', linewidth=0, zorder=1))

# Hatched cells for individual NaN within observed rows
for yi, yr in enumerate(all_years):
    row = df_obs[df_obs['Year'] == yr]
    if not row.empty:
        vals = pd.to_numeric(row.iloc[0][domain_cols].values, errors='coerce')
        for xi, v in enumerate(vals):
            if np.isnan(v):
                ax.add_patch(mpatches.Rectangle(
                    (xi-0.5, yi-0.5), 1, 1,
                    facecolor=CLR_HATCH_FACE, edgecolor=CLR_HATCH_EDGE,
                    hatch='////', linewidth=0, zorder=4))

# Observed data (NaN transparent so hatching shows through)
masked = np.ma.masked_where(np.isnan(orig_matrix), orig_matrix)
cmap2  = mcolors.ListedColormap([CLR_ZERO, CLR_ONE])
norm2  = mcolors.BoundaryNorm([0, 0.5, 1.0], cmap2.N)
cmap2.set_bad(color='none')
ax.imshow(masked, aspect='auto', cmap=cmap2, norm=norm2, interpolation='none',
          extent=[-0.5, n_d-0.5, n_y-0.5, -0.5], zorder=2)

# Row separators and section dividers
for i in range(n_y+1):
    ax.axhline(y=i-0.5, color='white', linewidth=0.35, zorder=3)
for _, lo, _, _ in SECTIONS:
    if lo in domain_ints:
        ax.axvline(x=domain_ints.index(lo)-0.5, color='white', linewidth=1.2, zorder=3)

# Y-axis
ax.set_yticks(range(n_y))
y_labels = [str(yr) + ('*' if yr in POOR_QUALITY_YEARS else '') for yr in all_years]
ax.set_yticklabels(y_labels, fontsize=8.5)
for i, yr in enumerate(all_years):
    if yr in obs_years:
        ax.get_yticklabels()[i].set_fontweight('bold')
ax.set_ylabel('Year', fontsize=11, labelpad=6)

# X-axis
tick_pos = [i for i, d in enumerate(domain_ints) if d == 1 or d % 10 == 0]
ax.set_xticks(tick_pos)
ax.set_xticklabels([str(domain_ints[i]) for i in tick_pos], fontsize=8)
ax.set_xlabel('CASCADE Domain   (1 = Cape Point  →  90 = Pea Island / N. Rodanthe)',
              fontsize=10, labelpad=4)
ax.set_title('Hatteras Island — Overwash Observations 1984–2004',
             fontsize=13, fontweight='bold', pad=10)

# ══════════════════════════════════════════════════════════════════════
# STORM TIMELINE PANEL (only when SHOW_STORM_TIMELINE = True)
# ══════════════════════════════════════════════════════════════════════
if SHOW_STORM_TIMELINE:
    axt.set_xlim(0, 1)
    axt.set_ylim(n_y-0.5, -0.5)
    axt.set_facecolor('white')
    axt.axis('off')

    # Panel title
    axt.text(0.5, -0.95, 'Storm Events', va='center', ha='center',
             fontsize=9, fontweight='bold', color='#333333',
             transform=axt.transData)

    # Subtle year guides
    for i, yr in enumerate(all_years):
        lc = '#EEEEEE' if yr not in obs_years else '#E0E0E0'
        axt.axhline(y=i, color=lc, linewidth=0.5, zorder=0)

    AXIS_X   = 0.10
    CONN_END = AXIS_X + 0.10
    LABEL_X  = AXIS_X + 0.13
    axt.axvline(x=AXIS_X, color='#CCCCCC', linewidth=1.2, zorder=1)

    # Build label positions — spread multi-storm years vertically
    storm_positions = []
    for yr, storm_list in STORMS_BY_YEAR.items():
        if yr < PERIOD[0] or yr > PERIOD[1]: continue
        yi = all_years.index(yr)
        n  = len(storm_list)
        offsets = ([0.0]             if n == 1 else
                   [-0.28, 0.28]     if n == 2 else
                   [-(n-1)*0.25 + i*0.25 for i in range(n)])
        for (name, cat, matched), offset in zip(storm_list, offsets):
            storm_positions.append((yi, yi+offset, name, cat, matched))
    storm_positions.sort(key=lambda x: x[1])

    # Draw storm dots, connectors, and labels
    for (yi_exact, yi_label, name, cat, matched) in storm_positions:
        sz, color = cat_style(cat, matched)
        alpha = 1.0 if matched else 0.65
        # Dot at exact year
        axt.plot(AXIS_X, yi_exact, 'o', color=color, markersize=sz, zorder=4,
                 markeredgecolor='white', markeredgewidth=0.5, alpha=alpha)
        # Connector: vertical jog if label raw_offset, then horizontal to label
        if abs(yi_label - yi_exact) > 0.05:
            axt.plot([AXIS_X, AXIS_X], [yi_exact, yi_label],
                     '-', color=color, linewidth=0.6, alpha=0.35, zorder=2)
        axt.plot([AXIS_X, CONN_END], [yi_label, yi_label],
                 '-', color=color, linewidth=0.7, alpha=0.45, zorder=2)
        # Label
        axt.text(LABEL_X, yi_label, f"{name}  ({cat})",
                 va='center', ha='left', fontsize=7,
                 fontweight='bold' if matched else 'normal',
                 fontstyle='normal' if matched else 'italic',
                 color=color, alpha=alpha)

    # Imagery years with no named storm
    for yr in all_years:
        yi = all_years.index(yr)
        if yr in STORMS_BY_YEAR:
            continue
        if yr in NO_EVENT_OBS_YEARS:
            axt.plot(AXIS_X, yi, 'o', color='#BBBBBB', markersize=4, zorder=3,
                     markeredgecolor='white', markeredgewidth=0.5)
            axt.text(LABEL_X, yi, '(no named storm)', va='center', ha='left',
                     fontsize=6.5, color='#BBBBBB', fontstyle='italic')
        else:
            # Quiet year — open circle
            axt.plot(AXIS_X, yi, 'o', color='#DDDDDD', markersize=3.5, zorder=3,
                     markeredgecolor='#CCCCCC', markeredgewidth=0.5, fillstyle='none')

    # 2004 note
    axt.text(LABEL_X, all_years.index(2004)+0.55,
             '† imagery predates Alex', va='top', ha='left',
             fontsize=5.5, color='#BBBBBB', fontstyle='italic')

    # Intensity mini-legend
    leg_y = n_y + 0.2
    axt.text(0.0, leg_y, 'Intensity scale:', va='top', ha='left',
             fontsize=6, color='#555555', fontweight='bold')
    leg_y += 0.7
    for cat, label in [('H5','H5'),('H3','H3'),('H1','H1'),('NE 5','NE 5'),('NE 4','NE 4'),('TS','TS')]:
        sz, color = cat_style(cat, True)
        axt.plot(0.05, leg_y, 'o', color=color, markersize=sz,
                 markeredgecolor='white', markeredgewidth=0.4)
        axt.text(0.16, leg_y, label, va='center', ha='left', fontsize=6, color=color)
        leg_y += 0.55
    axt.plot(0.05, leg_y, 'o', color='#DDDDDD', markersize=3.5,
             markeredgecolor='#CCCCCC', markeredgewidth=0.5, fillstyle='none')
    axt.text(0.16, leg_y, 'No storm / no imagery', va='center', ha='left',
             fontsize=6, color='#AAAAAA')

# ══════════════════════════════════════════════════════════════════════
# SECTION ANNOTATION BAR
# ══════════════════════════════════════════════════════════════════════
axb.set_xlim(-0.5, n_d-0.5); axb.set_ylim(0, 1); axb.axis('off')
for sec_name, lo, hi, stype in SECTIONS:
    lo_idx = domain_ints.index(lo) if lo in domain_ints else 0
    hi_idx = domain_ints.index(hi) if hi in domain_ints else n_d-1
    mid    = (lo_idx + hi_idx) / 2
    bc = CLR_BAR_VILLAGE if stype == "village" else CLR_BAR_INTER
    lc = '#3D2B1F'       if stype == "village" else '#1A3545'
    axb.add_patch(mpatches.FancyBboxPatch(
        (lo_idx-0.5, 0.52), hi_idx-lo_idx+1, 0.48,
        boxstyle="square,pad=0", lw=0.8, edgecolor='white', facecolor=bc))
    kw = dict(ha='center', va='top', color=lc, multialignment='center')
    axb.text(mid, 0.38, sec_name,
             fontsize=6 if sec_name == "Buxton" else (6.2 if '\n' in sec_name else 7.2),
             rotation=90 if sec_name == "Buxton" else 0, **kw)

# ══════════════════════════════════════════════════════════════════════
# LEGEND
# ══════════════════════════════════════════════════════════════════════
hatch_patch = mpatches.Patch(
    facecolor=CLR_HATCH_FACE, edgecolor=CLR_HATCH_EDGE,
    hatch='////', label='No imagery / not assessed')
fig.legend(
    handles=[
        mpatches.Patch(fc=CLR_ONE,  ec='#888', lw=0.5, label='Overwash observed'),
        mpatches.Patch(fc=CLR_ZERO, ec='#888', lw=0.5, label='No overwash (assessed)'),
        hatch_patch,
    ],
    loc='lower center', bbox_to_anchor=(legend_anchor, 0.05),
    ncol=3, fontsize=8.5, framealpha=0.95, edgecolor='#BBBBBB',
    columnspacing=1.4, handlelength=1.6, handleheight=1.2)

# ══════════════════════════════════════════════════════════════════════
# FOOTNOTE
# ══════════════════════════════════════════════════════════════════════
footnote = (
    'Bold years = imagery available  |  * poor image quality  |  '
    '† Alex (H3) Jul–Aug 2004 post-dates May 2004 imagery; Isabel (H5, Sep 2003) '
    'most likely accounts for overwash visible in 2004 imagery\n'
    'Sources: Hapke & Henderson (2007); Google Earth'
)
if SHOW_STORM_TIMELINE:
    footnote += (
        '\nTimeline: bold label = storm in imagery year; italic = no imagery that year  |  '
        'Dot size scaled by peak intensity'
    )
fig.text(heatmap_x, 0.02, footnote,
         fontsize=6.5, color='#555555', va='bottom', linespacing=1.5)

# ══════════════════════════════════════════════════════════════════════
# SAVE
# ══════════════════════════════════════════════════════════════════════
plt.savefig(OUTPUT_FILE, dpi=DPI, bbox_inches='tight', facecolor='white')
print(f"Saved: {OUTPUT_FILE}  (timeline={'ON' if SHOW_STORM_TIMELINE else 'OFF'})")
plt.show()
