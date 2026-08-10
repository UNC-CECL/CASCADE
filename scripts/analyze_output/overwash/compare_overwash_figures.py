"""
Hatteras Island — CASCADE vs Observed Overwash Comparison
==========================================================

Figure 1  (PLOT_STACKED = True)
    Top:    CASCADE Qow heatmap, imagery years highlighted
    Bottom: Observed overwash (binary)
    Shared domain x-axis + section annotations

Figure 2  (PLOT_CONTINGENCY = True)
    Imagery years only, equal-height rows
    Hit / Miss / False Alarm / Correct Rejection per domain
    Summary stats printed to console

USAGE
-----
Set NPZ_PATH, OBS_XLSX_PATH, OUT_DIR, and QOW_THRESHOLD, then run.
"""

import os, io, pickle, zipfile, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
warnings.filterwarnings('ignore')

# ══════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════
NPZ_PATH      = r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_1984_2004_basestorms_Hs2p0\HAT_1984_2004_basestorms_Hs2p0.npz"
OBS_XLSX_PATH = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\overwash_analysis\Hatteras_Overwash_Data.xlsx"
OUT_DIR       = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\overwash_analysis"

START_YEAR          = 1984
END_YEAR            = 2004    # model has END_YEAR - START_YEAR timesteps
NUM_REAL_DOMAINS    = 90
NUM_BUFFER_DOMAINS  = 15
FIRST_GIS_DOMAIN_ID = 1

PLOT_STACKED     = True
PLOT_CONTINGENCY = True
PLOT_SPATIAL     = True    # spatial profile comparison (two versions)

# Model threshold for "predicted overwash"
# Any cell where Qow > QOW_THRESHOLD is counted as a model overwash prediction.
# Start at 0 (any non-zero flux), raise if too many false alarms.
QOW_THRESHOLD = 0.0    # dam³/yr

# Imagery year remapping
# Some imagery dates post-date the storm they capture. Map them to the model
# year that best represents the event visible in the image.
# e.g. 2004 imagery (May 2004) shows Isabel 2003 deposits → remap to 2003.
YEAR_REMAP = {
    2004: 2003,    # May 2004 imagery captures Hurricane Isabel (Sep 2003) deposits
}

DPI = 200

# ══════════════════════════════════════════════════════════════════════
# ISLAND SECTIONS
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
POOR_QUALITY_YEARS = {1996}

# ══════════════════════════════════════════════════════════════════════
# UNIFIED COLOUR PALETTE
# ──────────────────────────────────────────────────────────────────────
# Concept: muted coastal earth & ocean tones, consistent warmth across
# all elements.  Section bar uses desaturated versions of the same
# warm/cool families used in the contingency legend.
# Inspired by AGU / Nature Geoscience figure aesthetics.
# ══════════════════════════════════════════════════════════════════════

# ── Geography (section bar) — linen/maritime, very desaturated ────────
CLR_BAR_VILLAGE = '#CABB9E'   # warm linen / sand
CLR_BAR_INTER   = '#8DAFC2'   # dusty maritime blue

# ── Contingency categories ─────────────────────────────────────────────
CLR_HIT  = '#2E7D5A'   # deep maritime green   — model & obs agree on overwash
CLR_MISS = '#2A5F8F'   # deep maritime blue    — obs yes, model missed
CLR_FA   = '#A85840'   # warm brick            — model yes, obs no (over-prediction)
CLR_CR   = '#F2EFE9'   # warm cream/parchment  — both show no overwash
CLR_NA   = '#BEB9B4'   # warm greige           — no imagery this year

# ── Observation binary (stacked figure) ────────────────────────────────
CLR_OBS_ON   = '#A85840'   # brick — overwash observed (same family as FA)
CLR_OBS_OFF  = '#FFFFFF'   # white — no overwash (assessed)
CLR_HATCH_F  = '#F2F2EE'   # warm off-white — no imagery background
CLR_HATCH_E  = '#C0BBB5'   # warm mid-grey  — hatch lines
CLR_IMG_BAND = '#EDF3F8'   # very light maritime blue — imagery row highlight

# ══════════════════════════════════════════════════════════════════════
# CASCADE NPZ LOADER
# ══════════════════════════════════════════════════════════════════════
_cls_reg = {}

def _dummy(module, name):
    key = (module, name)
    if key not in _cls_reg:
        def _ss(self, s):
            if isinstance(s, dict): self.__dict__.update(s)
            else: self._s = s
        _cls_reg[key] = type(f"{module}.{name}", (),
            {"__init__": lambda self, *a, **k: None, "__setstate__": _ss})
    return _cls_reg[key]

class _FU(pickle.Unpickler):
    def find_class(self, m, n):
        try: return super().find_class(m, n)
        except: return _dummy(m, n)

def load_qow_matrix(npz_path):
    n_yrs  = END_YEAR - START_YEAR
    s_idx  = NUM_BUFFER_DOMAINS
    e_idx  = s_idx + NUM_REAL_DOMAINS
    gis_ids = np.arange(FIRST_GIS_DOMAIN_ID,
                         FIRST_GIS_DOMAIN_ID + NUM_REAL_DOMAINS, dtype=int)
    with zipfile.ZipFile(npz_path) as zf:
        npy_m = next(m for m in zf.namelist() if m.endswith(".npy"))
        raw   = zf.read(npy_m)
    casc = _FU(io.BytesIO(raw[128:])).load()
    casc = casc.flat[0] if isinstance(casc, np.ndarray) else casc
    Q_rows = []
    for b in casc._barrier3d[s_idx:e_idx]:
        for attr in ["_QowTS","QowTS","_Qow_TS","Qow_TS","_Qow","Qow"]:
            if hasattr(b, attr):
                ts = np.asarray(getattr(b, attr), float).squeeze()
                if ts.ndim == 1 and ts.size > 1:
                    Q_rows.append(ts[:n_yrs]); break
    Q     = np.column_stack(Q_rows)
    years = np.arange(START_YEAR, START_YEAR + Q.shape[0], dtype=int)
    return years, Q, gis_ids

# ══════════════════════════════════════════════════════════════════════
# OBSERVATION LOADER
# ══════════════════════════════════════════════════════════════════════
def load_obs_matrix(xlsx_path, years_model):
    df = pd.read_excel(xlsx_path, sheet_name='Overwash_Matrix', skiprows=3, header=0)
    df['Year'] = pd.to_numeric(df['Year'], errors='coerce')
    df = df[df['Year'].notna()].copy()
    df['Year'] = df['Year'].astype(int)
    domain_cols = [c for c in df.columns
                   if str(c).replace('.0','').strip().isdigit()]
    domain_ints = [int(float(c)) for c in domain_cols]
    n_y = len(years_model)
    n_d = len(domain_ints)
    mat = np.full((n_y, n_d), np.nan)
    obs_years = []
    for yi, yr in enumerate(years_model):
        # Check if any imagery year remaps to this model year
        imagery_yr = yr
        for img_yr, model_yr in YEAR_REMAP.items():
            if model_yr == yr:
                imagery_yr = img_yr
                break

        row = df[df['Year'] == imagery_yr]
        if not row.empty:
            mat[yi, :] = pd.to_numeric(row.iloc[0][domain_cols].values,
                                        errors='coerce')
            obs_years.append(yr)   # store as MODEL year for alignment
    return mat, np.array(obs_years, dtype=int), domain_ints

# ══════════════════════════════════════════════════════════════════════
# SHARED DRAWING HELPERS
# ══════════════════════════════════════════════════════════════════════
def section_dividers(ax, n_domains, zorder=3):
    """White vertical lines at section boundaries. x = 0-indexed domain."""
    section_starts = [s[1] - 1 for s in SECTIONS]   # 0-indexed
    for x in section_starts:
        ax.axvline(x=x - 0.5, color='white', linewidth=1.1, zorder=zorder)

def section_bar(axb, n_domains):
    axb.set_xlim(-0.5, n_domains - 0.5)
    axb.set_ylim(0, 1); axb.axis('off')
    for sec_name, lo, hi, stype in SECTIONS:
        lo_i = lo - 1; hi_i = hi - 1   # 0-indexed
        mid  = (lo_i + hi_i) / 2
        bc   = CLR_BAR_VILLAGE if stype == "village" else CLR_BAR_INTER
        lc   = '#4A3728'       if stype == "village" else '#1E3D52'
        axb.add_patch(mpatches.FancyBboxPatch(
            (lo_i - 0.5, 0.52), hi_i - lo_i + 1, 0.48,
            boxstyle="square,pad=0", lw=0.8, edgecolor='white', facecolor=bc))
        kw = dict(ha='center', va='top', color=lc, multialignment='center')
        if sec_name == "Buxton":
            axb.text(mid, 0.38, 'Buxton', ha='center', va='top',
                     fontsize=6, color=lc, rotation=90)
        else:
            axb.text(mid, 0.38, sec_name,
                     fontsize=6.2 if '\n' in sec_name else 7, **kw)

# ══════════════════════════════════════════════════════════════════════
# FIGURE 1 — STACKED DUAL-PANEL
# ══════════════════════════════════════════════════════════════════════
def plot_stacked(years_m, Q, gis_ids, obs_mat, obs_years, domain_ints):
    n_y = len(years_m)
    n_d = len(domain_ints)

    # --- Layout ---
    fig  = plt.figure(figsize=(15, 12))
    ax_m = fig.add_axes([0.07, 0.42, 0.82, 0.50])   # model  (top 50%)
    ax_o = fig.add_axes([0.07, 0.19, 0.82, 0.21])   # obs    (middle 21%)
    ax_b = fig.add_axes([0.07, 0.12, 0.82, 0.05])   # section bar

    # ── MODEL HEATMAP ─────────────────────────────────────────────────
    # Highlight imagery year rows on model panel
    for yi, yr in enumerate(years_m):
        if yr in obs_years:
            ax_m.axhspan(yi - 0.5, yi + 0.5,
                         facecolor=CLR_IMG_BAND, alpha=0.55, zorder=0)

    vmax = max(np.nanpercentile(Q[Q > 0], 99), 1.0) if np.any(Q > 0) else 1.0
    im = ax_m.imshow(
        Q, aspect='auto', cmap='YlOrRd',
        norm=mcolors.Normalize(vmin=0, vmax=vmax),
        extent=[-0.5, n_d-0.5, n_y-0.5, -0.5],
        origin='upper', interpolation='nearest', zorder=1)

    section_dividers(ax_m, n_d)

    cb = fig.colorbar(im, ax=ax_m, fraction=0.022, pad=0.01, aspect=35)
    cb.set_label('Overwash Flux  Qow  (dam³/yr)', fontsize=9)
    cb.ax.tick_params(labelsize=8)

    ax_m.set_yticks(range(n_y))
    y_labels = [str(yr) + ('*' if yr in POOR_QUALITY_YEARS else '') for yr in years_m]
    ax_m.set_yticklabels(y_labels, fontsize=7.5)
    for i, yr in enumerate(years_m):
        if yr in obs_years:
            ax_m.get_yticklabels()[i].set_fontweight('bold')
            ax_m.get_yticklabels()[i].set_color('#1A4080')
    ax_m.set_ylabel('Year', fontsize=10, labelpad=6)
    ax_m.set_xticks([])
    ax_m.set_title(
        '(a)  CASCADE Model — Annual Overwash Flux Qow\n'
        'Blue bold years = aerial imagery available for that year',
        fontsize=11, fontweight='bold', loc='left', pad=6)

    # ── OBSERVATION HEATMAP ───────────────────────────────────────────
    ax_o.set_facecolor(CLR_HATCH_F)

    # Hatch unobserved years
    for yi, yr in enumerate(years_m):
        if yr not in obs_years:
            ax_o.add_patch(mpatches.Rectangle(
                (-0.5, yi-0.5), n_d, 1,
                facecolor=CLR_HATCH_F, edgecolor=CLR_HATCH_E,
                hatch='////', linewidth=0, zorder=1))

    # Individual NaN cells
    for yi, yr in enumerate(years_m):
        if yr in obs_years:
            for xi, v in enumerate(obs_mat[yi, :]):
                if np.isnan(v):
                    ax_o.add_patch(mpatches.Rectangle(
                        (xi-0.5, yi-0.5), 1, 1,
                        facecolor=CLR_HATCH_F, edgecolor=CLR_HATCH_E,
                        hatch='////', linewidth=0, zorder=4))

    # Observed data
    masked = np.ma.masked_where(np.isnan(obs_mat), obs_mat)
    cmap_o = mcolors.ListedColormap([CLR_OBS_OFF, CLR_OBS_ON])
    cmap_o.set_bad(color='none')
    ax_o.imshow(masked, aspect='auto',
                cmap=cmap_o, norm=mcolors.BoundaryNorm([0, 0.5, 1.0], 2),
                extent=[-0.5, n_d-0.5, n_y-0.5, -0.5],
                origin='upper', interpolation='none', zorder=2)

    section_dividers(ax_o, n_d)
    for i in range(n_y+1):
        ax_o.axhline(y=i-0.5, color='white', linewidth=0.3, zorder=3)

    ax_o.set_yticks(range(n_y))
    ax_o.set_yticklabels(y_labels, fontsize=7.5)
    for i, yr in enumerate(years_m):
        if yr in obs_years:
            ax_o.get_yticklabels()[i].set_fontweight('bold')
    ax_o.set_ylabel('Year', fontsize=10, labelpad=6)
    ax_o.set_xticks([])

    # Obs legend
    hh = mpatches.Patch(fc=CLR_HATCH_F, ec=CLR_HATCH_E, hatch='////', label='No imagery')
    ax_o.legend(handles=[
        mpatches.Patch(fc=CLR_OBS_ON,  ec='#888', lw=0.5, label='Overwash observed'),
        mpatches.Patch(fc=CLR_OBS_OFF, ec='#888', lw=0.5, label='No overwash (assessed)'),
        hh], loc='upper right', fontsize=7.5, framealpha=0.92,
        edgecolor='#BBBBBB', ncol=3, bbox_to_anchor=(1.0, 1.04))
    ax_o.set_title('(b)  Observed Overwash — Aerial Imagery',
                   fontsize=11, fontweight='bold', loc='left', pad=6)

    section_bar(ax_b, n_d)
    fig.text(0.48, 0.075,
             'CASCADE Domain   (1 = Cape Point  →  90 = Pea Island / N. Rodanthe)',
             ha='center', va='top', fontsize=10)
    fig.text(0.48, 0.975,
             'CASCADE Model vs Observed Overwash — Hatteras Island 1984–2003',
             ha='center', va='top', fontsize=12, fontweight='bold')
    fig.text(0.07, 0.03,
             '* poor image quality  |  '
             + ('  '.join([f'† {img}→{mdl}: imagery from {img} assigned to model year {mdl} (deposits from that storm season)'
                           for img, mdl in YEAR_REMAP.items()])
                + '  |  ' if YEAR_REMAP else '')
             + 'Sources: Hapke & Henderson (2007); Google Earth  |  '
             f'Model: {os.path.basename(os.path.dirname(NPZ_PATH))}',
             fontsize=6.5, color='#555555', va='top')

    out = os.path.join(OUT_DIR, 'compare_stacked.png')
    plt.savefig(out, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f"Saved: {out}")
    plt.show()

# ══════════════════════════════════════════════════════════════════════
# FIGURE 2 — CONTINGENCY (imagery years only, equal-height rows)
# ══════════════════════════════════════════════════════════════════════
def plot_contingency(years_m, Q, gis_ids, obs_mat, obs_years, domain_ints):
    n_y  = len(years_m)    # ALL model years (e.g. 20 for 1984-2003)
    n_d  = len(domain_ints)
    obs_set = set(obs_years.tolist())

    # Build full n_y × n_d contingency matrix — all model years as rows
    # 0 = Correct Rejection, 1 = Hit, 2 = False Alarm, 3 = Miss, 4 = Not Assessed
    # Non-imagery years → all 4 (grey "not assessed")
    cont = np.full((n_y, n_d), 4, dtype=float)

    for yi, yr in enumerate(years_m):
        if yr not in obs_set:
            continue    # leave as 4 (grey) — no imagery this year
        q_row = Q[yi, :]
        o_row = obs_mat[yi, :]
        for xi in range(n_d):
            if np.isnan(o_row[xi]):
                cont[yi, xi] = 4; continue
            pred = q_row[xi] > QOW_THRESHOLD
            obs  = o_row[xi] == 1
            if   pred and     obs: cont[yi, xi] = 1   # Hit
            elif pred and not obs: cont[yi, xi] = 2   # False Alarm
            elif not pred and obs: cont[yi, xi] = 3   # Miss
            else:                  cont[yi, xi] = 0   # Correct Rejection

    # Summary stats (imagery years only)
    h  = np.sum(cont == 1); m  = np.sum(cont == 3)
    fa = np.sum(cont == 2); cr = np.sum(cont == 0)
    tot = h + m + fa + cr
    print(f"\nContingency (threshold = {QOW_THRESHOLD} dam³/yr):")
    print(f"  Hits             {h:5d}  ({100*h/tot:.1f}%)")
    print(f"  Misses           {m:5d}  ({100*m/tot:.1f}%)")
    print(f"  False Alarms    {fa:5d}  ({100*fa/tot:.1f}%)")
    print(f"  Correct Rej.    {cr:5d}  ({100*cr/tot:.1f}%)")
    if h+m   > 0: print(f"  POD = {h/(h+m):.3f}")
    if h+m+fa> 0: print(f"  CSI = {h/(h+m+fa):.3f}")

    # --- Figure ---
    row_h = 0.40
    fig_h = max(6.0, n_y * row_h + 4.5)   # extra space at bottom for legend
    fig   = plt.figure(figsize=(15, fig_h))
    # More room at bottom for legend
    ax_c  = fig.add_axes([0.07, 0.26, 0.82, 0.65])
    ax_b  = fig.add_axes([0.07, 0.18, 0.82, 0.05])

    # ── Colormap using unified palette ────────────────────────────────
    # 0=CR, 1=Hit, 2=FA, 3=Miss, 4=NA
    cmap_c = mcolors.ListedColormap([CLR_CR, CLR_HIT, CLR_FA, CLR_MISS, CLR_NA])
    norm_c = mcolors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], cmap_c.N)

    # Extent uses ROW INDICES → every year gets identical row height
    ax_c.imshow(cont, aspect='auto', cmap=cmap_c, norm=norm_c,
                extent=[-0.5, n_d-0.5, n_y-0.5, -0.5],
                origin='upper', interpolation='nearest', zorder=2)

    section_dividers(ax_c, n_d)
    for i in range(n_y + 1):
        ax_c.axhline(y=i - 0.5, color='white', linewidth=0.35, zorder=3)

    # Y-axis: one tick per model year, bold blue for imagery years
    ax_c.set_yticks(range(n_y))
    ylabels = []
    for yr in years_m:
        lbl = str(yr)
        if yr in POOR_QUALITY_YEARS: lbl += '*'
        ylabels.append(lbl)
    ax_c.set_yticklabels(ylabels, fontsize=8)
    for i, yr in enumerate(years_m):
        if yr in obs_set:
            ax_c.get_yticklabels()[i].set_fontweight('bold')
            ax_c.get_yticklabels()[i].set_color('#1A4080')

    ax_c.set_ylabel('Year', fontsize=10, labelpad=6)
    ax_c.set_xticks([])
    pod_str = f"{h/(h+m):.2f}" if h+m > 0 else "—"
    csi_str = f"{h/(h+m+fa):.2f}" if h+m+fa > 0 else "—"
    far_str = f"{fa/(fa+cr):.2f}" if fa+cr > 0 else "—"
    ax_c.set_title(
        f'CASCADE vs Observed — Contingency  |  Qow threshold = {QOW_THRESHOLD} dam³/yr\n'
        f'Blue bold years = imagery available  |  '
        f'POD = {pod_str}   CSI = {csi_str}   False Alarm Rate = {far_str}',
        fontsize=11, fontweight='bold', loc='left', pad=6)

    # Legend: placed OUTSIDE axes in figure space below section bar
    # so it never covers any data
    leg_handles = [
        mpatches.Patch(fc=CLR_HIT,  ec='#1A5C3A', lw=0.5,
                       label='Hit — model & obs both show overwash'),
        mpatches.Patch(fc=CLR_MISS, ec='#1A3F62', lw=0.5,
                       label='Miss — observed overwash, model did not predict'),
        mpatches.Patch(fc=CLR_FA,   ec='#7A3A28', lw=0.5,
                       label='False Alarm — model predicted, obs shows none'),
        mpatches.Patch(fc=CLR_CR,   ec='#AAA49E', lw=0.5,
                       label='Correct Rejection — neither shows overwash'),
        mpatches.Patch(fc=CLR_NA,   ec='#8A847F', lw=0.5,
                       label='Not assessed — no imagery this year'),
    ]
    fig.legend(handles=leg_handles, fontsize=8, framealpha=0.95,
               edgecolor='#BBBBBB', ncol=3,
               loc='lower center',
               bbox_to_anchor=(0.48, 0.09))

    section_bar(ax_b, n_d)
    fig.text(0.48, 0.145,
             'CASCADE Domain   (1 = Cape Point  →  90 = Pea Island / N. Rodanthe)',
             ha='center', va='top', fontsize=10)
    fig.text(0.48, 0.975,
             'CASCADE vs Observed Overwash — Contingency Analysis  |  Hatteras Island',
             ha='center', va='top', fontsize=12, fontweight='bold')
    fig.text(0.07, 0.02,
             f'* poor image quality  |  Threshold = {QOW_THRESHOLD} dam³/yr  |  '
             'POD = Hits/(Hits+Misses)   CSI = Hits/(Hits+Misses+False Alarms)  |  '
             + ('  '.join([f'† {img}→{mdl}: imagery from {img} assigned to model year {mdl}'
                           for img, mdl in YEAR_REMAP.items()])
                + '  |  ' if YEAR_REMAP else '')
             + f'Model: {os.path.basename(os.path.dirname(NPZ_PATH))}',
             fontsize=6.5, color='#555555', va='bottom')

    out = os.path.join(OUT_DIR, f'compare_contingency_T{QOW_THRESHOLD}.png')
    plt.savefig(out, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f"Saved: {out}")
    plt.show()

# ══════════════════════════════════════════════════════════════════════
# SPATIAL PROFILE COMPARISON — shared data prep helper
# ══════════════════════════════════════════════════════════════════════
def compute_spatial_data(Q, obs_mat, domain_ints):
    """
    Returns per-domain summary arrays used by both spatial profile figures.

    qow_mean  : mean annual Qow per domain across all model years (dam³/yr)
    qow_p25   : 25th percentile of annual Qow per domain
    qow_p75   : 75th percentile of annual Qow per domain
    obs_freq  : fraction of assessed imagery years with overwash per domain
                (NaN where no assessed imagery exists for that domain)
    bar_colors: colour per domain based on section type (village vs inter-village)
    """
    n_d = len(domain_ints)

    qow_mean = Q.mean(axis=0)
    qow_p25  = np.percentile(Q, 25, axis=0)
    qow_p75  = np.percentile(Q, 75, axis=0)

    obs_freq = np.full(n_d, np.nan)
    for xi in range(n_d):
        col      = obs_mat[:, xi]
        assessed = col[~np.isnan(col)]
        if len(assessed) > 0:
            obs_freq[xi] = np.sum(assessed == 1) / len(assessed)

    # Colour each domain bar by section type for geographic context
    bar_colors = []
    village_doms = set()
    for _, lo, hi, stype in SECTIONS:
        if stype == "village":
            village_doms.update(range(lo, hi+1))
    for d in domain_ints:
        bar_colors.append(CLR_BAR_VILLAGE if d in village_doms else CLR_BAR_INTER)

    return qow_mean, qow_p25, qow_p75, obs_freq, bar_colors


def _draw_section_background(ax, domain_ints, alpha=0.08):
    """Very subtle alternating section shading on an axes."""
    village_doms = set()
    for _, lo, hi, stype in SECTIONS:
        if stype == "village": village_doms.update(range(lo, hi+1))
    for i, d in enumerate(domain_ints):
        if d in village_doms:
            ax.axvspan(i-0.5, i+0.5, facecolor=CLR_BAR_VILLAGE,
                       alpha=alpha, zorder=0, linewidth=0)


# ══════════════════════════════════════════════════════════════════════
# OPTION 1 — NORMALISED SINGLE-PANEL
# Both model and observations scaled to 0–1 on the same y-axis.
# ══════════════════════════════════════════════════════════════════════
def plot_spatial_normalised(Q, obs_mat, domain_ints, obs_years, years_m):
    qow_mean, qow_p25, qow_p75, obs_freq, bar_colors = \
        compute_spatial_data(Q, obs_mat, domain_ints)

    n_d    = len(domain_ints)
    x      = np.arange(n_d)
    n_img  = len(obs_years)

    # Normalise both to 0–1
    qow_norm   = qow_mean / np.nanmax(qow_mean)
    qow_p25_n  = qow_p25  / np.nanmax(qow_mean)
    qow_p75_n  = qow_p75  / np.nanmax(qow_mean)
    freq_norm  = obs_freq   # already 0–1

    fig = plt.figure(figsize=(16, 7))
    ax  = fig.add_axes([0.07, 0.32, 0.88, 0.57])   # more bottom space for 2-line footnote
    axb = fig.add_axes([0.07, 0.18, 0.88, 0.06])

    _draw_section_background(ax, domain_ints)

    # ── Observation frequency bars (behind curve) ─────────────────────
    bars = ax.bar(x, freq_norm, width=0.85, color=bar_colors,
                  alpha=0.75, zorder=2, linewidth=0,
                  label=f'Observed overwash frequency\n(fraction of {n_img} imagery years)')

    # ── Model IQR shading ─────────────────────────────────────────────
    ax.fill_between(x, qow_p25_n, qow_p75_n,
                    color=CLR_MISS, alpha=0.18, zorder=3,
                    label='Modelled Qow IQR (25th–75th pct.)')

    # ── Model mean curve ──────────────────────────────────────────────
    ax.plot(x, qow_norm, color=CLR_MISS, linewidth=1.8,
            zorder=4, label='Modelled mean annual Qow (normalised)')
    ax.fill_between(x, 0, qow_norm, color=CLR_MISS, alpha=0.10, zorder=3)

    # ── Section dividers ──────────────────────────────────────────────
    for _, lo, _, _ in SECTIONS:
        ax.axvline(x=lo-1-0.5, color='#AAAAAA', linewidth=0.7,
                   linestyle='--', alpha=0.6, zorder=5)

    ax.set_xlim(-0.5, n_d-0.5)
    ax.set_ylim(0, 1.08)
    ax.set_xticks([i for i,d in enumerate(domain_ints) if d==1 or d%10==0])
    ax.set_xticklabels([str(domain_ints[i])
                        for i,d in enumerate(domain_ints) if d==1 or d%10==0],
                       fontsize=8)
    ax.set_ylabel('Normalised value  (0 – 1)', fontsize=10)
    ax.set_xlabel('CASCADE Domain   (1 = Cape Point  →  90 = Pea Island / N. Rodanthe)',
                  fontsize=10, labelpad=8)
    ax.yaxis.grid(True, color='#CCCCCC', linewidth=0.5, alpha=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title('Modelled Overwash Intensity vs Observed Frequency — Hatteras Island '
                 f'{years_m[0]}–{years_m[-1]}',
                 fontsize=11, fontweight='bold', loc='left', pad=6)

    # Single consolidated legend — upper centre
    all_handles = [
        mpatches.Patch(fc=CLR_MISS, alpha=0.18, ec=CLR_MISS, lw=0,
                       label='Modelled Qow IQR (25th–75th percentile)'),
        plt.Line2D([0], [0], color=CLR_MISS, linewidth=2.0,
                   label='Modelled mean annual Qow (normalised)'),
        mpatches.Patch(fc=CLR_BAR_INTER, ec='none', alpha=0.75,
                       label='Observed overwash frequency — inter-village domain'),
        mpatches.Patch(fc=CLR_BAR_VILLAGE, ec='none', alpha=0.75,
                       label='Observed overwash frequency — village domain'),
    ]
    ax.legend(handles=all_handles, fontsize=8.2, framealpha=0.93,
              edgecolor='#BBBBBB', loc='upper center',
              bbox_to_anchor=(0.5, 0.99), ncol=2,
              handlelength=1.8, columnspacing=1.2)

    section_bar(axb, n_d)

    # Professional footnote — two lines to prevent horizontal overflow
    norm_line = (
        f'Both series normalised by division by their respective maxima '
        f'(mean Qow$_{{max}}$ = {qow_mean.max():.1f} dam\u00b3/yr; '
        f'obs. freq.$_{{max}}$ = {np.nanmax(obs_freq):.2f}, '
        f'fraction of {n_img} imagery years), '
        f'enabling direct shape comparison independent of physical dimensions.')
    source_line = (
        f'Model: {os.path.basename(os.path.dirname(NPZ_PATH))}  |  '
        f'Obs sources: Hapke \u0026 Henderson (2007); Google Earth'
        + (f'  |  ' + '  '.join([f'{img}\u2192{mdl} remap applied'
           for img, mdl in YEAR_REMAP.items()]) if YEAR_REMAP else ''))
    fig.text(0.07, 0.075, norm_line + '\n' + source_line,
             fontsize=6.5, color='#555555', va='top', linespacing=1.4)

    out = os.path.join(OUT_DIR, 'spatial_profile_normalised.png')
    plt.savefig(out, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f"Saved: {out}")
    plt.show()


# ══════════════════════════════════════════════════════════════════════
# OPTION 2 — DUAL PANEL, SHARED X-AXIS
# Model and observations on separate panels, each with its own y-axis.
# ══════════════════════════════════════════════════════════════════════
def plot_spatial_dual_panel(Q, obs_mat, domain_ints, obs_years, years_m):
    qow_mean, qow_p25, qow_p75, obs_freq, bar_colors = \
        compute_spatial_data(Q, obs_mat, domain_ints)

    n_d   = len(domain_ints)
    x     = np.arange(n_d)
    n_img = len(obs_years)

    fig  = plt.figure(figsize=(16, 9))
    ax_m = fig.add_axes([0.07, 0.50, 0.88, 0.41])   # model (top)  — shorter, higher up
    ax_o = fig.add_axes([0.07, 0.25, 0.88, 0.19])   # obs   (bottom) — gap of 0.06 above
    axb  = fig.add_axes([0.07, 0.13, 0.88, 0.06])   # section bar

    # ── TOP: modelled Qow curve ───────────────────────────────────────
    _draw_section_background(ax_m, domain_ints)

    ax_m.fill_between(x, qow_p25, qow_p75,
                      color=CLR_MISS, alpha=0.20, zorder=2,
                      label='IQR across years (25th–75th percentile)')
    ax_m.fill_between(x, 0, qow_mean,
                      color=CLR_MISS, alpha=0.12, zorder=3)
    ax_m.plot(x, qow_mean, color=CLR_MISS, linewidth=2.0,
              zorder=4, label='Mean annual Qow')

    for _, lo, _, _ in SECTIONS:
        ax_m.axvline(x=lo-1-0.5, color='#AAAAAA', linewidth=0.7,
                     linestyle='--', alpha=0.6, zorder=5)

    ax_m.set_xlim(-0.5, n_d-0.5)
    ax_m.set_ylim(bottom=0)
    ax_m.set_xticks([])
    ax_m.yaxis.grid(True, color='#CCCCCC', linewidth=0.5, alpha=0.5, zorder=0)
    ax_m.set_axisbelow(True)
    ax_m.set_ylabel('Mean annual Qow  (dam³/yr)', fontsize=10)
    ax_m.legend(fontsize=8.5, framealpha=0.93, edgecolor='#BBBBBB',
                loc='upper center')
    ax_m.set_title('(a)  CASCADE Model — Mean Annual Overwash Flux per Domain  '
                   f'({years_m[0]}–{years_m[-1]})',
                   fontsize=10, fontweight='bold', loc='left', pad=6)

    # ── BOTTOM: observation frequency bars ───────────────────────────
    _draw_section_background(ax_o, domain_ints)

    ax_o.bar(x, obs_freq, width=0.85, color=bar_colors,
             alpha=0.85, zorder=2, linewidth=0)

    # Count labels only for bars where count >= 2 to reduce clutter
    for xi, (f, d) in enumerate(zip(obs_freq, domain_ints)):
        if not np.isnan(f) and f > 0:
            count = round(f * n_img)
            if count >= 2:
                ax_o.text(xi, f + 0.012, str(count),
                          ha='center', va='bottom', fontsize=6,
                          color='#333333', zorder=6)

    for _, lo, _, _ in SECTIONS:
        ax_o.axvline(x=lo-1-0.5, color='#AAAAAA', linewidth=0.7,
                     linestyle='--', alpha=0.6, zorder=5)

    ax_o.set_xlim(-0.5, n_d-0.5)
    ax_o.set_ylim(0, min(1.15, np.nanmax(obs_freq)*1.3 + 0.05)
                  if np.nanmax(obs_freq) > 0 else 1.0)
    ax_o.yaxis.grid(True, color='#CCCCCC', linewidth=0.5, alpha=0.5, zorder=0)
    ax_o.set_axisbelow(True)
    ax_o.set_xticks([i for i,d in enumerate(domain_ints) if d==1 or d%10==0])
    ax_o.set_xticklabels([str(domain_ints[i])
                          for i,d in enumerate(domain_ints) if d==1 or d%10==0],
                         fontsize=8)
    ax_o.set_ylabel('Fraction of imagery\nyears with overwash', fontsize=9)
    ax_o.set_xlabel('CASCADE Domain   (1 = Cape Point  →  90 = Pea Island / N. Rodanthe)',
                    fontsize=10, labelpad=4)

    # Secondary right-axis tick labels showing raw counts
    ax_o_r = ax_o.twinx()
    ax_o_r.set_ylim(ax_o.get_ylim())
    ax_o_r.set_yticks([i/n_img for i in range(n_img+1)])
    ax_o_r.set_yticklabels([str(i) for i in range(n_img+1)], fontsize=7)
    ax_o_r.set_ylabel(f'Count of imagery\nyears  (max {n_img})', fontsize=8)
    ax_o_r.spines['right'].set_color('#AAAAAA')
    ax_o_r.tick_params(colors='#AAAAAA')
    ax_o_r.yaxis.label.set_color('#AAAAAA')

    # Section colour bar legend
    village_patch = mpatches.Patch(fc=CLR_BAR_VILLAGE, ec='none',
                                   label='Village domain')
    inter_patch   = mpatches.Patch(fc=CLR_BAR_INTER,   ec='none',
                                   label='Inter-village domain')
    ax_o.legend(handles=[village_patch, inter_patch], fontsize=8,
                framealpha=0.9, edgecolor='#BBBBBB',
                loc='upper right', bbox_to_anchor=(0.99, 0.99))
    ax_o.set_title(f'(b)  Observed Overwash Frequency per Domain\n'
                   f'Numbers \u2265\u20092 shown above bars (count of imagery years '
                   f'with overwash, out of {n_img} assessed)',
                   fontsize=10, fontweight='bold', loc='left', pad=6)

    section_bar(axb, n_d)

    fig.text(0.48, 0.975,
             'Modelled Overwash Intensity vs Observed Frequency — Hatteras Island',
             ha='center', va='top', fontsize=12, fontweight='bold')
    fig.text(0.07, 0.04,
             f'Model: {os.path.basename(os.path.dirname(NPZ_PATH))}  |  '
             f'Observation sources: Hapke & Henderson (2007); Google Earth  |  '
             + ('  '.join([f'{img}→{mdl} remap applied'
                           for img, mdl in YEAR_REMAP.items()])
                if YEAR_REMAP else ''),
             fontsize=6.5, color='#555555', va='top')

    out = os.path.join(OUT_DIR, 'spatial_profile_dual_panel.png')
    plt.savefig(out, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f"Saved: {out}")
    plt.show()


# ══════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════
def main():
    print("Loading CASCADE model ...")
    years_m, Q, gis_ids = load_qow_matrix(NPZ_PATH)
    print(f"  {years_m[0]}–{years_m[-1]}, {Q.shape[1]} domains, max Qow = {Q.max():.2f}")

    print("Loading observations ...")
    obs_mat, obs_years, domain_ints = load_obs_matrix(OBS_XLSX_PATH, years_m)
    print(f"  {len(obs_years)} imagery years in model range: {list(obs_years)}")

    if PLOT_STACKED:
        print("\nGenerating stacked comparison ...")
        plot_stacked(years_m, Q, gis_ids, obs_mat, obs_years, domain_ints)

    if PLOT_CONTINGENCY:
        print("\nGenerating contingency figure ...")
        plot_contingency(years_m, Q, gis_ids, obs_mat, obs_years, domain_ints)

    if PLOT_SPATIAL:
        print("\nGenerating spatial profile — normalised ...")
        plot_spatial_normalised(Q, obs_mat, domain_ints, obs_years, years_m)
        print("Generating spatial profile — dual panel ...")
        plot_spatial_dual_panel(Q, obs_mat, domain_ints, obs_years, years_m)

    print("\nDone.")

if __name__ == "__main__":
    main()
