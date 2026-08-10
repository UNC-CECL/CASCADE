#!/usr/bin/env python3
"""
HAT_management_rules_table.py
Single combined management rules table as a publication/presentation figure.
"""
import textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib import rcParams

PNG = r'C:\Users\hanna\PycharmProjects\CASCADE\scripts\management_plot\tables\HAT_management_rules_table.png'
PDF = r'C:\Users\hanna\PycharmProjects\CASCADE\scripts\management_plot\tables\HAT_management_rules_table.pdf'

rcParams.update({'font.family': 'DejaVu Serif', 'font.size': 9})

# ── COLOURS ───────────────────────────────────────────────────────────────────
C_TITLE_BG  = '#1c3a5e'
C_COL_BG    = '#2c5f9e'
C_NOUR_SEC  = '#9e6a1a'    # amber — nourishment section header
C_ROAD_SEC  = '#404040'    # charcoal — road section header
C_NOUR_A    = '#fef8ee'    # light amber rows
C_NOUR_B    = '#fef2d8'
C_ROAD_A    = '#f6f6f6'    # light grey rows
C_ROAD_B    = '#ececec'
WHITE       = 'white'

# ── FIGURE SETUP ──────────────────────────────────────────────────────────────
FW, FH = 15, 7.6
fig = plt.figure(figsize=(FW, FH), facecolor='white')
ax  = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, FW); ax.set_ylim(0, FH); ax.axis('off')

# ── TABLE GEOMETRY ────────────────────────────────────────────────────────────
ML, MR = 0.28, 0.20          # left / right margin (inches)
TL = ML
TR = FW - MR
TW = TR - TL                  # total table width

# Column widths (inches); must sum to TW
CW = [1.90, 1.10, 1.20, TW - 1.90 - 1.10 - 1.20]   # Location | Domains | Period | Description

CX = [TL + sum(CW[:i]) for i in range(len(CW))]   # left x of each column

# Row heights
H_TITLE = 0.52
H_COL   = 0.46
H_SEC   = 0.36
H_ROW   = 0.50

# Build y positions top-to-bottom
rows_meta = (
    ['title']
    + ['col_hdr']
    + ['sec_nour']
    + ['nour'] * 3
    + ['sec_road']
    + ['road'] * 5
)
heights = {
    'title':    H_TITLE,
    'col_hdr':  H_COL,
    'sec_nour': H_SEC,
    'nour':     H_ROW,
    'sec_road': H_SEC,
    'road':     H_ROW,
}
total_h = sum(heights[r] for r in rows_meta)
y_top   = FH - (FH - total_h) / 2   # vertically centred

row_ys = []
y = y_top
for r in rows_meta:
    row_ys.append(y - heights[r])
    y -= heights[r]

# ── DRAW HELPERS ──────────────────────────────────────────────────────────────
PAD_L, PAD_R = 0.10, 0.08

def fill(x, y, w, h, color):
    ax.add_patch(Rectangle((x, y), w, h, fc=color, ec='none', zorder=2))

def divider(y, color='#cccccc', lw=0.5):
    ax.plot([TL, TR], [y, y], color=color, lw=lw, zorder=5)

def txt(x, y, w, h, text, fs=9, bold=False, italic=False,
        color='#111111', ha='left', va='center', wrap_w=None):
    if wrap_w:
        text = textwrap.fill(text, wrap_w)
    tx = x + PAD_L if ha == 'left' else (x + w - PAD_R if ha == 'right' else x + w/2)
    ty = y + h/2 if va == 'center' else (y + h - 0.08 if va == 'top' else y + 0.08)
    ax.text(tx, ty, text, fontsize=fs, color=color,
            ha=ha, va='center' if va=='center' else va,
            fontweight='bold' if bold else 'normal',
            fontstyle='italic' if italic else 'normal',
            family='DejaVu Serif', zorder=6, clip_on=False)

def section_hdr(y, h, label, bg, accent):
    fill(TL, y, TW, h, bg + '33')          # very subtle bg tint
    fill(TL, y, 0.18, h, accent)            # left accent bar
    ax.text(TL + 0.25, y + h/2, label,
            fontsize=10, color=WHITE, ha='left', va='center',
            fontweight='bold', family='DejaVu Serif', zorder=6)

# ── ROWS: TITLE ───────────────────────────────────────────────────────────────
ri = 0
y_r = row_ys[ri]; h_r = H_TITLE
fill(TL, y_r, TW, h_r, C_TITLE_BG)
ax.text(TL + TW/2, y_r + h_r/2,
        'Hatteras Island CASCADE Hindcast: Management Rules Summary (1984–2024)',
        fontsize=11.5, color=WHITE, ha='center', va='center',
        fontweight='bold', family='DejaVu Serif', zorder=6)

# ── ROWS: COLUMN HEADERS ──────────────────────────────────────────────────────
ri = 1
y_r = row_ys[ri]; h_r = H_COL
fill(TL, y_r, TW, h_r, C_COL_BG)
for i, lbl in enumerate(['Location', 'CASCADE\nDomains', 'Period\nActive', 'Description']):
    ha_i = 'center' if i in (1, 2) else 'left'
    ax.text(
        CX[i] + (CW[i]/2 if ha_i == 'center' else PAD_L),
        y_r + h_r/2, lbl,
        fontsize=9, color=WHITE, ha=ha_i, va='center',
        fontweight='bold', family='DejaVu Serif', zorder=6)
divider(y_r, color=C_COL_BG, lw=1.2)

# ── ROWS: NOURISHMENT SECTION HEADER ─────────────────────────────────────────
ri = 2
y_r = row_ys[ri]; h_r = H_SEC
fill(TL, y_r, TW, h_r, '#f5e8d0')
fill(TL, y_r, 0.18, h_r, C_NOUR_SEC)
ax.text(TL + 0.25, y_r + h_r/2,
        'Beach Nourishment Events  —  sediment added to specified domains',
        fontsize=9.5, color='#5a3a00', ha='left', va='center',
        fontweight='bold', family='DejaVu Serif', zorder=6)

# ── NOURISHMENT DATA ROWS ─────────────────────────────────────────────────────
NOUR_DATA = [
    ('Rodanthe', 'D85–88', '2014',
     'Emergency fill following storm erosion; material placed at southern Pea Island NWR boundary '
     'adjacent to Rodanthe  (1,620,000 cy)'),
    ('Buxton',   'D6–15',  '2022',
     'USACE emergency shore-protection nourishment project'),
    ('Avon',     'D23–26', '2022',
     'USACE emergency shore-protection nourishment project'),
]

for k, (loc, dom, per, desc) in enumerate(NOUR_DATA):
    ri = 3 + k
    y_r = row_ys[ri]; h_r = H_ROW
    bg  = C_NOUR_A if k % 2 == 0 else C_NOUR_B
    fill(TL, y_r, TW, h_r, bg)
    txt(CX[0], y_r, CW[0], h_r, loc,  fs=9)
    txt(CX[1], y_r, CW[1], h_r, dom,  fs=9, ha='center')
    txt(CX[2], y_r, CW[2], h_r, per,  fs=9, ha='center')
    txt(CX[3], y_r, CW[3], h_r, desc, fs=8.8, wrap_w=80)
    divider(y_r, color='#e0cba8', lw=0.4)

# ── ROAD SECTION HEADER ───────────────────────────────────────────────────────
ri = 6
y_r = row_ys[ri]; h_r = H_SEC
fill(TL, y_r, TW, h_r, '#e8e8e8')
fill(TL, y_r, 0.18, h_r, C_ROAD_SEC)
ax.text(TL + 0.25, y_r + h_r/2,
        'Road Relocation — Disabled Zones  —  NC-12 position held fixed',
        fontsize=9.5, color='#1a1a1a', ha='left', va='center',
        fontweight='bold', family='DejaVu Serif', zorder=6)

# ── ROAD DATA ROWS ────────────────────────────────────────────────────────────
ROAD_DATA = [
    ('Cape Hatteras',                    'D1–6',    '1984–2024',
     'No NC-12 road present; cape terminus with no through-road infrastructure'),
    ('Buxton',                           'D7–8',    '1984–2024',
     'Developed community — dense built infrastructure prevents landward road relocation'),
    ('Avon',                             'D21–31',  '1984–2024',
     'Developed community — dense built infrastructure prevents landward road relocation'),
    ('Tri-Village\n(Salvo / Waves / Rodanthe)', 'D68–83', '1984–2024',
     'Developed community — dense built infrastructure prevents landward road relocation'),
    ('N. Rodanthe /\nPea Island NWR',    'D82–88',  '2022–2024',
     'Jug Handle Bridge (July 2022): NC-12 permanently elevated off surface — '
     'road relocation physically impossible while bridge is in place'),
]

for k, (loc, dom, per, desc) in enumerate(ROAD_DATA):
    ri = 7 + k
    y_r = row_ys[ri]; h_r = H_ROW
    bg  = C_ROAD_A if k % 2 == 0 else C_ROAD_B
    fill(TL, y_r, TW, h_r, bg)
    txt(CX[0], y_r, CW[0], h_r, loc,  fs=9)
    txt(CX[1], y_r, CW[1], h_r, dom,  fs=9, ha='center')
    txt(CX[2], y_r, CW[2], h_r, per,  fs=9, ha='center')
    txt(CX[3], y_r, CW[3], h_r, desc, fs=8.8, wrap_w=80)
    divider(y_r, color='#cccccc', lw=0.4)

# ── OUTER BORDER ─────────────────────────────────────────────────────────────
# Top and bottom thick rules; left/right thin
y_bot = row_ys[-1]
y_top_tbl = y_top
for yy, lw in [(y_top_tbl, 1.4), (y_bot, 1.4)]:
    ax.plot([TL, TR], [yy, yy], color='#333333', lw=lw, zorder=7)
ax.plot([TL, TL], [y_bot, y_top_tbl], color='#aaaaaa', lw=0.6, zorder=7)
ax.plot([TR, TR], [y_bot, y_top_tbl], color='#aaaaaa', lw=0.6, zorder=7)

# Column separator lines (very subtle)
for cx in CX[1:]:
    ax.plot([cx, cx], [y_bot, row_ys[1]], color='#dddddd', lw=0.4, zorder=5)

# ── FOOTER ────────────────────────────────────────────────────────────────────
ax.text(TL, y_bot - 0.12,
        'Note: Inter-village zones (D9–20, D32–67) retain road relocation capability throughout. '
        'Domain 1 = Cape Hatteras (south); Domain 90 = Oregon Inlet (north).',
        fontsize=7.8, color='#555555', ha='left', va='top',
        family='DejaVu Serif', fontstyle='italic', zorder=6)

# ── SAVE ─────────────────────────────────────────────────────────────────────
fig.savefig(PNG, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(PDF,           bbox_inches='tight', facecolor='white')
print(f'Saved: {PNG}')
plt.close(fig)
