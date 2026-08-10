#!/usr/bin/env python3
"""
HAT_management_timeline_figure.py  (v8)
Changes:
  - Cape Point → Cape Hatteras, colored as inter-village
  - N. Rodanthe label → Pea Island NWR
  - Period 1/2 shown as geological-epoch bar below x-axis (much clearer)
  - Background period tints removed (cleaner data area)
  - Point-event width reduced to 0.5 yr (true year slice)
  - Buxton corrected to D7–8; Buxton–Avon inter-village D9–20
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import matplotlib.ticker as ticker
from matplotlib import rcParams

PNG = r'C:\Users\hanna\PycharmProjects\CASCADE\scripts\figure_making\management_plot\figures\HAT_management_1984_2024_timeline.png'
PDF = r'C:\Users\hanna\PycharmProjects\CASCADE\scripts\figure_making\management_plot\figures\HAT_management_1984_2024_timeline.pdf'

Y0, Y1  = 1984, 2024
PBREAK  = 2004
PW      = 0.5    # point-event width — true year slice

# ── COLOURS ───────────────────────────────────────────────────────────────────
C_N  = '#c8912a'   # amber/sand    → nourishment
C_R  = '#5a5a5a'   # charcoal      → road
C_B  = '#2c5f9e'   # steel blue    → bridge

C_VILLAGE = '#e4eef7'   # soft blue  — community zones
C_INTERV  = '#efefeb'   # warm grey  — inter-village / undeveloped

# Period epoch-bar colours (more saturated than the original background tints)
C_P1_BAR = '#aecde8'   # clear blue
C_P2_BAR = '#f0dc90'   # clear amber-yellow

# ── COMMUNITY ZONES ───────────────────────────────────────────────────────────
COMMUNITIES = [
    # label                                  d0  d1   facecolor
    ('Cape Hatteras',                         1,  6,  C_INTERV),   # not a village → inter-village colour
    ('Buxton',                                7,  8,  C_VILLAGE),  # D7–8 only
    ('Buxton–Avon\n(inter-village)',          9, 20,  C_INTERV),
    ('Avon',                                 21, 31,  C_VILLAGE),
    ('Avon–Tri-Village\n(inter-village)',    32, 67,  C_INTERV),
    ('Tri-Village: Salvo /\nWaves / Rodanthe', 68, 83, C_VILLAGE),
    ('Pea Island NWR',                       84, 90,  C_INTERV),   # not a village community
]

# ── EVENT CATALOGUE ───────────────────────────────────────────────────────────
EVENTS = [
    # Nourishment (point events)
    (2014, 2014, 85, 88, 'nourish', ''),
    (2022, 2022,  6, 15, 'nourish', ''),
    (2022, 2022, 23, 26, 'nourish', ''),
    # Road relocations (point events)
    (1989, 1989, 84, 87, 'road',    ''),
    (1999, 1999,  9, 15, 'road',    ''),
    # Bridge (persistent)
    (2022, 2024, 82, 88, 'bridge',  ''),
]

STYLE = {
    'nourish': (C_N, '#8a620e', 0.90, 5, 0.7),
    'road':    (C_R, '#333333', 0.88, 4, 0.7),
    'bridge':  (C_B, '#1a3f70', 0.88, 4, 0.7),
}

# ── GLOBAL RCPARAMS ───────────────────────────────────────────────────────────
rcParams.update({
    'font.family':      'DejaVu Sans',
    'font.size':        9,
    'axes.linewidth':   0.8,
    'figure.facecolor': 'white',
    'axes.facecolor':   'white',
})

fig, ax = plt.subplots(figsize=(14, 9.0))

# ── COMMUNITY BANDS ───────────────────────────────────────────────────────────
for _, d0, d1, fc in COMMUNITIES:
    ax.axhspan(d0-0.5, d1+0.5, color=fc, alpha=1.0, zorder=1, lw=0)

for _, d0, _, _ in COMMUNITIES[1:]:
    ax.axhline(d0-0.5, color='#c0c0c0', lw=0.6, zorder=2)

# ── PERIOD BOUNDARY LINE ─────────────────────────────────────────────────────
ax.axvline(PBREAK, color='#444444', lw=1.5, ls='--', zorder=9, dashes=(5,3))
ax.text(PBREAK+0.3, 91.8, '2004',
        fontsize=7.5, color='#555555', ha='left', va='center',
        style='italic', zorder=10)

# ── EPOCH BAR — period indicator below x-axis ─────────────────────────────────
# Sits in the extended y range below domain 1
BAR_Y0, BAR_Y1 = -5.5, -2.5   # domain-unit coordinates for the bar

ax.add_patch(Rectangle((Y0,     BAR_Y0), PBREAK-Y0, BAR_Y1-BAR_Y0,
                        fc=C_P1_BAR, ec='#6a9fc0', lw=0.9, zorder=5,
                        clip_on=False))
ax.add_patch(Rectangle((PBREAK, BAR_Y0), Y1-PBREAK, BAR_Y1-BAR_Y0,
                        fc=C_P2_BAR, ec='#b09020', lw=0.9, zorder=5,
                        clip_on=False))

ax.text((Y0+PBREAK)/2, (BAR_Y0+BAR_Y1)/2,
        'Period 1  (1984–2004)',
        fontsize=8.5, fontweight='semibold', color='#1a3f6e',
        ha='center', va='center', zorder=6, clip_on=False)
ax.text((PBREAK+Y1)/2, (BAR_Y0+BAR_Y1)/2,
        'Period 2  (2004–2024)',
        fontsize=8.5, fontweight='semibold', color='#5a4000',
        ha='center', va='center', zorder=6, clip_on=False)

# Divider tick at 2004
ax.plot([PBREAK, PBREAK], [BAR_Y0, BAR_Y1],
        color='#888888', lw=1.0, zorder=7, clip_on=False)

# ── DRAW EVENT BARS ───────────────────────────────────────────────────────────
for (yr0, yr1, d0, d1, etype, _) in EVENTS:
    fc, ec, alpha, zo, lw = STYLE[etype]
    is_point = (yr0 == yr1)
    x0 = yr0 - PW/2 if is_point else yr0
    x1 = yr0 + PW/2 if is_point else yr1
    y0_, y1_ = d0-0.5, d1+0.5

    ax.add_patch(Rectangle((x0, y0_), x1-x0, y1_-y0_,
                            fc=fc, ec=ec, alpha=alpha,
                            linewidth=lw, zorder=zo))
    if not is_point:
        ax.plot([yr0, yr0], [y0_, y1_], color=ec, lw=2.0,
                solid_capstyle='butt', zorder=zo+1)

# ── ARROW CALLOUT ANNOTATIONS ─────────────────────────────────────────────────
def callout(text, xy, xytext, color, rad=0.0):
    ax.annotate(
        text, xy=xy, xytext=xytext,
        fontsize=8.0, color='#111111', ha='center', va='center',
        fontweight='normal', zorder=11,
        arrowprops=dict(arrowstyle='->', color=color, lw=1.0,
                        connectionstyle=f'arc3,rad={rad}'),
        bbox=dict(boxstyle='round,pad=0.35', lw=0.9,
                  fc='white', ec=color, alpha=0.97))

callout('Rodanthe emergency\nnourishment  (2014)',
        xy=(2014, 86.5), xytext=(2009.5, 74.5), color=C_N)
callout('Buxton beach\nnourishment  (2022)',
        xy=(2022, 10.5), xytext=(2015.5, 20), color=C_N)
callout('Avon beach\nnourishment  (2022)',
        xy=(2022, 24.5), xytext=(2015.5, 36), color=C_N)
callout('S-curve relocation  (1989)',
        xy=(1989, 85.5), xytext=(1992, 79), color=C_R)
callout('Post-Dennis relocation  (1999)',
        xy=(1999, 12), xytext=(1993, 22.5), color=C_R)
callout('Jug Handle Bridge  (2022–);\nNC-12 no longer managed',
        xy=(2023, 85), xytext=(2016.5, 62), color=C_B)

# ── RIGHT-AXIS COMMUNITY LABELS ───────────────────────────────────────────────
for lbl, d0, d1, fc in COMMUNITIES:
    mid   = (d0 + d1) / 2
    is_iv = (fc == C_INTERV)
    ax.text(Y1 + 1.2, mid, lbl,
            fontsize=7.8 if not is_iv else 7.2,
            color='#333333' if not is_iv else '#666666',
            va='center', ha='left',
            style='normal' if not is_iv else 'italic',
            clip_on=False, zorder=12)

# ── AXES ─────────────────────────────────────────────────────────────────────
ax.set_xlim(Y0-0.5, Y1+0.8)
ax.set_ylim(-6.5, 95.5)   # extended bottom for epoch bar

ax.set_xlabel('Year', fontsize=10, labelpad=36, color='#222222')
ax.set_ylabel('CASCADE Domain\n(1 = S / Cape Hatteras  →  90 = N / Oregon Inlet)',
              fontsize=9.5, labelpad=8, color='#222222')

ax.xaxis.set_major_locator(ticker.MultipleLocator(4))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.yaxis.set_major_locator(ticker.FixedLocator([0,10,20,30,40,50,60,70,80,90]))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(5))

ax.grid(axis='x', which='major', color='#d8d8d8', lw=0.6, zorder=1.5)
ax.grid(axis='x', which='minor', color='#efefef', lw=0.35, zorder=1.4)
ax.grid(axis='y', which='major', color='#e5e5e5', lw=0.35, zorder=1.4)

for spine in ('top', 'right'):
    ax.spines[spine].set_visible(False)
ax.spines['left'].set_color('#aaaaaa')
ax.spines['bottom'].set_color('#aaaaaa')
ax.tick_params(axis='both', labelsize=8.5, color='#888888',
               labelcolor='#333333')

# Hide tick labels in the negative y (epoch bar) region
ax.set_yticks([t for t in ax.get_yticks() if t >= 0])

# ── TITLE ─────────────────────────────────────────────────────────────────────
ax.set_title('Hatteras Island: Management Action History  (1984–2024)',
             fontsize=12, fontweight='bold', pad=10,
             color='#111111', loc='left', x=0.01)

# ── LEGEND ────────────────────────────────────────────────────────────────────
_blank = mpatches.Patch(fc='none', ec='none', label=' ')
legend_handles = [
    # Row 1 — management action types
    mpatches.Patch(fc=C_N,       ec='#8a620e', alpha=0.90, label='Beach nourishment'),
    mpatches.Patch(fc=C_R,       ec='#333333', alpha=0.88, label='NC-12 road relocation'),
    mpatches.Patch(fc=C_B,       ec='#1a3f70', alpha=0.88, label='Bridge / infrastructure (persistent)'),
    _blank,
    # Row 2 — background zones
    mpatches.Patch(fc=C_VILLAGE, ec='#aaaaaa', alpha=1.0, lw=0.6, label='Village / community zone'),
    mpatches.Patch(fc=C_INTERV,  ec='#aaaaaa', alpha=1.0, lw=0.6, label='Inter-village / undeveloped zone'),
    mpatches.Patch(fc=C_P1_BAR,  ec='#6a9fc0', alpha=1.0, lw=0.8, label='Period 1  (1984–2004)'),
    mpatches.Patch(fc=C_P2_BAR,  ec='#b09020', alpha=1.0, lw=0.8, label='Period 2  (2004–2024)'),
]
leg = ax.legend(
    handles=legend_handles,
    loc='upper center', bbox_to_anchor=(0.42, -0.175),
    ncol=4, fontsize=8.2, frameon=True, framealpha=0.97,
    edgecolor='#cccccc', title='Management action',
    title_fontsize=8.8,
    handlelength=1.6, handleheight=0.9,
    borderpad=0.7, columnspacing=1.2,
)
leg.get_frame().set_linewidth(0.8)

plt.tight_layout(rect=[0, 0.10, 0.87, 0.97])
fig.savefig(PNG, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(PDF,           bbox_inches='tight', facecolor='white')
print(f'Saved: {PNG}')
plt.close(fig)
