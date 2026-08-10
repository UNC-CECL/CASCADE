"""
Hatteras Island — CASCADE vs Observed Overwash Comparison
==========================================================
Produces two figures:

  Figure 1 — Stacked dual-panel (PLOT_STACKED = True)
      Top panel  : CASCADE Qow heatmap (continuous, all model years)
      Bottom panel: Observed overwash  (binary, imagery years + hatching)
      Both share domain x-axis and island section annotations.
      Imagery years are flagged on the model panel.

  Figure 2 — Contingency heatmap (PLOT_CONTINGENCY = True)
      For each imagery year × domain where observations exist, classifies
      each cell as Hit / Miss / False Alarm / Correct Rejection.
      Requires a Qow threshold to binarise model comparison.

USAGE
-----
1. Set NPZ_PATH and OBS_XLSX_PATH below.
2. Adjust PERIOD, domain constants, and SECTIONS to match your run.
3. Set QOW_THRESHOLD for the contingency analysis.
4. python compare_overwash.py
"""

import os, io, pickle, zipfile, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
warnings.filterwarnings('ignore')

# ══════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit here
# ══════════════════════════════════════════════════════════════════════

# --- Paths ---
NPZ_PATH      = r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_1984_2004_basestorms_Hs2p0\HAT_1984_2004_basestorms_Hs2p0.npz"
OBS_XLSX_PATH = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\overwash_analysis\Hatteras_Overwash_Data.xlsx"
OUT_DIR       = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\overwash_analysis"

# --- Model run ---
START_YEAR          = 1984
END_YEAR            = 2004    # exclusive upper bound (model has END_YEAR - START_YEAR rows)
NUM_REAL_DOMAINS    = 90
NUM_BUFFER_DOMAINS  = 15
FIRST_GIS_DOMAIN_ID = 1

# --- Which figures to produce ---
PLOT_STACKED     = True    # stacked dual-panel comparison
PLOT_CONTINGENCY = True    # contingency heatmap (Hit / Miss / FA / CR)

# --- Contingency threshold ---
# Model cell is "predicted overwash" if Qow > QOW_THRESHOLD
# Start with 0 (any non-zero Qow counts); raise if too many false alarms.
QOW_THRESHOLD = 0.0        # dam³/yr

# --- Figure comparison ---
DPI = 200

# ══════════════════════════════════════════════════════════════════════
# ISLAND SECTIONS — matches CASCADE ANN_TOWN_SPANS
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
CLR_BAR_VILLAGE = '#C4A882'
CLR_BAR_INTER   = '#7FA8C4'

# ══════════════════════════════════════════════════════════════════════
# COLOURS
# ══════════════════════════════════════════════════════════════════════
CLR_OBS_ONE        = '#C0392B'   # observed overwash
CLR_OBS_ZERO       = '#FFFFFF'   # assessed, no overwash
CLR_HATCH_FACE     = '#F2F2F2'   # no imagery background
CLR_HATCH_EDGE     = '#BBBBBB'   # hatch lines
CLR_IMAGERY_BAND   = '#E8F0FA'   # highlight band on model panel for imagery years

# Contingency colours
CLR_HIT  = '#2E7D32'   # green  — model & obs both show overwash
CLR_MISS = '#1565C0'   # blue   — obs yes, model no
CLR_FA   = '#E65100'   # orange — model yes, obs no (false alarm)
CLR_CR   = '#F5F5F5'   # light grey — both show no overwash
CLR_NA   = '#DDDDDD'   # medium grey — not assessed (NaN in obs)

POOR_QUALITY_YEARS = {1996}

# ══════════════════════════════════════════════════════════════════════
# CASCADE NPZ LOADER (from plot_overwash.py — duplicated for portability)
# ══════════════════════════════════════════════════════════════════════
_class_registry = {}

def _make_dummy_class(module, name):
    key = (module, name)
    if key not in _class_registry:
        def _ss(self, state):
            if isinstance(state, dict): self.__dict__.update(state)
            else: self._unpickled_state = state
        _class_registry[key] = type(f"{module}.{name}", (),
            {"__init__": lambda self, *a, **k: None, "__setstate__": _ss})
    return _class_registry[key]

class _FlexUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        try: return super().find_class(module, name)
        except (ModuleNotFoundError, AttributeError):
            return _make_dummy_class(module, name)

def load_qow_matrix(npz_path):
    """Load CASCADE NPZ and return (years, Q, gis_domain_ids)."""
    total  = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS
    s_idx  = NUM_BUFFER_DOMAINS
    e_idx  = s_idx + NUM_REAL_DOMAINS
    n_yrs  = END_YEAR - START_YEAR
    gis_ids = np.arange(FIRST_GIS_DOMAIN_ID,
                        FIRST_GIS_DOMAIN_ID + NUM_REAL_DOMAINS, dtype=int)

    with zipfile.ZipFile(npz_path) as zf:
        npy_m = next(m for m in zf.namelist() if m.endswith(".npy"))
        raw   = zf.read(npy_m)
    buf      = io.BytesIO(raw[128:])
    wrapper  = _FlexUnpickler(buf).load()
    casc     = wrapper.flat[0] if isinstance(wrapper, np.ndarray) else wrapper

    b3d_list = casc._barrier3d
    Q_rows   = []
    for b in b3d_list[s_idx:e_idx]:
        for attr in ["_QowTS","QowTS","_Qow_TS","Qow_TS","_Qow","Qow"]:
            if hasattr(b, attr):
                ts = np.asarray(getattr(b, attr), float).squeeze()
                if ts.ndim == 1 and ts.size > 1:
                    Q_rows.append(ts[:n_yrs]); break

    Q     = np.column_stack(Q_rows)          # shape (n_yrs, n_real_domains)
    years = np.arange(START_YEAR, START_YEAR + Q.shape[0], dtype=int)
    return years, Q, gis_ids

# ══════════════════════════════════════════════════════════════════════
# OBSERVATION LOADER
# ══════════════════════════════════════════════════════════════════════
def load_obs_matrix(xlsx_path, years_model):
    """Load observation Excel and return (obs_matrix, obs_years_in_model)."""
    df = pd.read_excel(xlsx_path, sheet_name='Overwash_Matrix', skiprows=3, header=0)
    df['Year'] = pd.to_numeric(df['Year'], errors='coerce')
    df = df[df['Year'].notna()].copy()
    df['Year'] = df['Year'].astype(int)

    domain_cols = [c for c in df.columns
                   if str(c).replace('.0','').strip().isdigit()]
    domain_ints = [int(float(c)) for c in domain_cols]

    # Build matrix aligned to model years
    n_y   = len(years_model)
    n_d   = len(domain_ints)
    mat   = np.full((n_y, n_d), np.nan)
    obs_years = []
    for yi, yr in enumerate(years_model):
        row = df[df['Year'] == yr]
        if not row.empty:
            mat[yi, :] = pd.to_numeric(row.iloc[0][domain_cols].values,
                                        errors='coerce')
            obs_years.append(yr)
    return mat, np.array(obs_years), domain_ints

# ══════════════════════════════════════════════════════════════════════
# SHARED HELPERS
# ══════════════════════════════════════════════════════════════════════
def draw_section_dividers(ax, domain_ids, color='white', lw=1.2, zorder=3):
    """Draw vertical section boundary lines on a heatmap axes."""
    for _, lo, _, _ in SECTIONS:
        if lo in domain_ids:
            ax.axvline(x=lo - 0.5, color=color, linewidth=lw, zorder=zorder)

def draw_section_bar(axb, domain_ids):
    """Draw the island section colour bar on axb."""
    n_d = len(domain_ids)
    axb.set_xlim(domain_ids[0] - 0.5, domain_ids[-1] + 0.5)
    axb.set_ylim(0, 1); axb.axis('off')
    for sec_name, lo, hi, stype in SECTIONS:
        lo_px = lo - 0.5; hi_px = hi + 0.5
        mid   = (lo + hi) / 2
        bc    = CLR_BAR_VILLAGE if stype == "village" else CLR_BAR_INTER
        lc    = '#3D2B1F'       if stype == "village" else '#1A3545'
        axb.add_patch(mpatches.FancyBboxPatch(
            (lo_px, 0.52), hi_px - lo_px, 0.48,
            boxstyle="square,pad=0", lw=0.8, edgecolor='white', facecolor=bc))
        kw = dict(ha='center', va='top', color=lc, multialignment='center')
        axb.text(mid, 0.38, sec_name,
                 fontsize=6 if sec_name=="Buxton" else (6.2 if '\n' in sec_name else 7),
                 rotation=90 if sec_name=="Buxton" else 0, **kw)

# ══════════════════════════════════════════════════════════════════════
# FIGURE 1 — STACKED DUAL-PANEL
# ══════════════════════════════════════════════════════════════════════
def plot_stacked(years_m, Q, gis_ids, obs_mat, obs_years, domain_ints):
    n_y_m = len(years_m)
    n_y_o = n_y_m           # both panels show the same year range
    n_d   = len(gis_ids)

    fig = plt.figure(figsize=(15, 12))

    # Axes: model top (55%), obs bottom (30%), section bar (8%)
    ax_m  = fig.add_axes([0.07, 0.44, 0.82, 0.50])   # model heatmap
    ax_o  = fig.add_axes([0.07, 0.19, 0.82, 0.23])   # observation heatmap
    ax_b  = fig.add_axes([0.07, 0.12, 0.82, 0.05])   # section bar

    # ── MODEL PANEL ───────────────────────────────────────────────────
    # Highlight imagery year rows (subtle background band)
    for yi, yr in enumerate(years_m):
        if yr in obs_years:
            ax_m.axhspan(yr - 0.5, yr + 0.5,
                         facecolor=CLR_IMAGERY_BAND, alpha=0.6,
                         zorder=0, linewidth=0)

    vmax = np.nanpercentile(Q[Q > 0], 99) if np.any(Q > 0) else 1.0
    cmap_m = plt.get_cmap('YlOrRd')
    norm_m = mcolors.Normalize(vmin=0, vmax=vmax)

    extent_m = [gis_ids[0]-0.5, gis_ids[-1]+0.5, years_m[-1]+0.5, years_m[0]-0.5]
    im = ax_m.imshow(Q, aspect='auto', cmap=cmap_m, norm=norm_m,
                     extent=extent_m, origin='upper', interpolation='nearest',
                     zorder=1)

    draw_section_dividers(ax_m, list(gis_ids), color='white', lw=1.0)

    # Left-edge triangles marking imagery year rows
    for yr in obs_years:
        ax_m.annotate('▶', xy=(gis_ids[0]-0.5, yr), xycoords='data',
                      fontsize=7, color='#1F4E79', va='center', ha='right',
                      fontweight='bold', zorder=4)

    # Colorbar
    cb = fig.colorbar(im, ax=ax_m, fraction=0.025, pad=0.01, aspect=30)
    cb.set_label('Overwash Flux  Qow  (dam³/yr)', fontsize=9)
    cb.ax.tick_params(labelsize=8)

    ax_m.set_yticks(years_m)
    ax_m.set_yticklabels(
        [str(yr) + ('▶' if yr in obs_years else '') for yr in years_m],
        fontsize=7.5)
    for i, yr in enumerate(years_m):
        if yr in obs_years:
            ax_m.get_yticklabels()[i].set_fontweight('bold')
            ax_m.get_yticklabels()[i].set_color('#1F4E79')
    ax_m.set_ylabel('Year', fontsize=10, labelpad=6)
    ax_m.set_xticks([]); ax_m.set_xlabel('')
    ax_m.set_title(
        f'(a)  CASCADE Model — Overwash Flux Qow\n'
        f'Blue ▶ years = imagery available for comparison',
        fontsize=11, fontweight='bold', loc='left', pad=6)

    # ── OBSERVATION PANEL ─────────────────────────────────────────────
    ax_o.set_facecolor(CLR_HATCH_FACE)

    # Hatched rows — unobserved years
    for yi, yr in enumerate(years_m):
        if yr not in obs_years:
            ax_o.add_patch(mpatches.Rectangle(
                (gis_ids[0]-0.5, yr-0.5),
                n_d, 1,
                facecolor=CLR_HATCH_FACE, edgecolor=CLR_HATCH_EDGE,
                hatch='////', linewidth=0, zorder=1))

    # Individual NaN cells in observed rows
    for yi, yr in enumerate(years_m):
        if yr in obs_years:
            oi = np.where(years_m == yr)[0][0]
            for xi, v in enumerate(obs_mat[oi, :]):
                if np.isnan(v):
                    d = domain_ints[xi]
                    ax_o.add_patch(mpatches.Rectangle(
                        (d-0.5, yr-0.5), 1, 1,
                        facecolor=CLR_HATCH_FACE, edgecolor=CLR_HATCH_EDGE,
                        hatch='////', linewidth=0, zorder=4))

    # Binary data
    obs_plot = np.full((n_y_m, n_d), np.nan)
    for yi, yr in enumerate(years_m):
        if yr in obs_years:
            oi = np.where(years_m == yr)[0][0]
            obs_plot[yi, :] = obs_mat[oi, :]

    masked_obs = np.ma.masked_where(np.isnan(obs_plot), obs_plot)
    cmap_o = mcolors.ListedColormap([CLR_OBS_ZERO, CLR_OBS_ONE])
    norm_o = mcolors.BoundaryNorm([0, 0.5, 1.0], cmap_o.N)
    cmap_o.set_bad(color='none')
    ax_o.imshow(masked_obs, aspect='auto', cmap=cmap_o, norm=norm_o,
                interpolation='none',
                extent=extent_m, origin='upper', zorder=2)

    draw_section_dividers(ax_o, list(gis_ids), color='white', lw=1.0)
    for i in range(n_y_m):
        yr = years_m[i]
        ax_o.axhline(y=yr-0.5, color='white', linewidth=0.3, zorder=3)

    ax_o.set_yticks(years_m)
    ax_o.set_yticklabels([str(yr) for yr in years_m], fontsize=7.5)
    for i, yr in enumerate(years_m):
        if yr in obs_years:
            ax_o.get_yticklabels()[i].set_fontweight('bold')
    ax_o.set_ylabel('Year', fontsize=10, labelpad=6)
    ax_o.set_xticks([])
    ax_o.set_xlim(gis_ids[0]-0.5, gis_ids[-1]+0.5)

    # Obs legend (inline, top right of obs panel)
    hatch_h = mpatches.Patch(fc=CLR_HATCH_FACE, ec=CLR_HATCH_EDGE,
                              hatch='////', label='No imagery')
    obs_handles = [
        mpatches.Patch(fc=CLR_OBS_ONE,  ec='#888', lw=0.5, label='Overwash observed'),
        mpatches.Patch(fc=CLR_OBS_ZERO, ec='#888', lw=0.5, label='No overwash (assessed)'),
        hatch_h,
    ]
    ax_o.legend(handles=obs_handles, loc='upper right', fontsize=7.5,
                framealpha=0.92, edgecolor='#BBBBBB', ncol=3,
                bbox_to_anchor=(1.0, 1.02))
    ax_o.set_title('(b)  Observed Overwash — Aerial Imagery\n'
                   'Bold years = imagery available',
                   fontsize=11, fontweight='bold', loc='left', pad=6)

    # ── SECTION BAR ───────────────────────────────────────────────────
    draw_section_bar(ax_b, list(gis_ids))

    # Shared x-label below section bar
    fig.text(0.48, 0.07, 'CASCADE Domain   (1 = Cape Point  →  90 = Pea Island / N. Rodanthe)',
             ha='center', va='top', fontsize=10)

    # Overall title
    fig.text(0.48, 0.97,
             'CASCADE Model Output vs Aerial Imagery Observations — Hatteras Island 1984–2003',
             ha='center', va='top', fontsize=12, fontweight='bold')

    # Footnote
    fig.text(0.07, 0.04,
             f'▶ on y-axis and blue year labels = imagery year  |  '
             f'Model: {os.path.basename(os.path.dirname(NPZ_PATH))}  |  '
             f'Obs sources: Hapke & Henderson (2007); Google Earth',
             fontsize=6.5, color='#555555', va='top')

    out = os.path.join(OUT_DIR, 'compare_stacked.png')
    plt.savefig(out, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f"Saved: {out}")
    plt.show()

# ══════════════════════════════════════════════════════════════════════
# FIGURE 2 — CONTINGENCY HEATMAP
# ══════════════════════════════════════════════════════════════════════
def plot_contingency(years_m, Q, gis_ids, obs_mat, obs_years, domain_ints):
    """
    For each imagery year × domain where obs data exists, classify as:
      Hit (H)             : Qow > T  AND  obs = 1   (green)
      Miss (M)            : Qow <= T AND  obs = 1   (blue)
      False Alarm (FA)    : Qow > T  AND  obs = 0   (orange)
      Correct Rejection(C): Qow <= T AND  obs = 0   (light grey)
      Not Assessed (NA)   : obs = NaN               (medium grey)
    Rows: imagery years only. Columns: all 90 domains.
    """
    n_obs = len(obs_years)
    n_d   = len(domain_ints)

    # Build contingency matrix (rows = imagery years, cols = domains)
    # Encoding: 0=CR, 1=Hit, 2=FA, 3=Miss, 4=NA
    cont = np.full((n_obs, n_d), 4, dtype=float)   # default NA

    for ri, yr in enumerate(obs_years):
        yi = np.where(years_m == yr)[0][0]
        q_row  = Q[yi, :]                  # shape (n_d,)
        o_row  = obs_mat[yi, :]            # shape (n_d,) — NaN where unassessed

        for xi in range(n_d):
            obs_val = o_row[xi]
            if np.isnan(obs_val):
                cont[ri, xi] = 4           # not assessed
                continue
            predicted = q_row[xi] > QOW_THRESHOLD
            observed  = obs_val == 1
            if   predicted and     observed: cont[ri, xi] = 1   # Hit
            elif predicted and not observed: cont[ri, xi] = 2   # False Alarm
            elif not predicted and observed: cont[ri, xi] = 3   # Miss
            else:                            cont[ri, xi] = 0   # Correct Rejection

    # Domain-level summary stats
    hits = np.sum(cont == 1, axis=0)
    miss = np.sum(cont == 3, axis=0)
    fa   = np.sum(cont == 2, axis=0)
    cr   = np.sum(cont == 0, axis=0)
    pod  = np.where(hits+miss > 0, hits/(hits+miss), np.nan)   # Prob of Detection
    csi  = np.where(hits+miss+fa > 0, hits/(hits+miss+fa), np.nan) # Critical Success Index

    # Overall stats
    h_tot = np.nansum(hits); m_tot = np.nansum(miss)
    f_tot = np.nansum(fa);   c_tot = np.nansum(cr)
    tot   = h_tot + m_tot + f_tot + c_tot
    print(f"\nContingency summary (threshold = {QOW_THRESHOLD} dam³/yr):")
    print(f"  Hits:              {h_tot:4d}  ({100*h_tot/tot:.1f}%)")
    print(f"  Misses:            {m_tot:4d}  ({100*m_tot/tot:.1f}%)")
    print(f"  False Alarms:      {f_tot:4d}  ({100*f_tot/tot:.1f}%)")
    print(f"  Correct Rej.:      {c_tot:4d}  ({100*c_tot/tot:.1f}%)")
    print(f"  Overall POD:       {h_tot/(h_tot+m_tot):.3f}" if h_tot+m_tot>0 else "  POD: N/A")
    print(f"  Overall CSI:       {h_tot/(h_tot+m_tot+f_tot):.3f}" if h_tot+m_tot+f_tot>0 else "  CSI: N/A")

    # ── FIGURE ────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(15, max(5, n_obs * 0.55 + 3.5)))
    ax_c = fig.add_axes([0.07, 0.30, 0.82, 0.58])
    ax_b = fig.add_axes([0.07, 0.20, 0.82, 0.05])
    ax_s = fig.add_axes([0.07, 0.05, 0.82, 0.10])  # summary stats strip

    # Colormap: 0=CR, 1=Hit, 2=FA, 3=Miss, 4=NA
    cmap_c = mcolors.ListedColormap([CLR_CR, CLR_HIT, CLR_FA, CLR_MISS, CLR_NA])
    norm_c = mcolors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], cmap_c.N)

    extent_c = [gis_ids[0]-0.5, gis_ids[-1]+0.5,
                obs_years[-1]+0.5, obs_years[0]-0.5]
    ax_c.imshow(cont, aspect='auto', cmap=cmap_c, norm=norm_c,
                extent=extent_c, origin='upper', interpolation='nearest', zorder=2)

    draw_section_dividers(ax_c, list(gis_ids), color='white', lw=1.0)
    for yr in obs_years:
        ax_c.axhline(y=yr-0.5, color='white', linewidth=0.3, zorder=3)

    ax_c.set_yticks(obs_years)
    ax_c.set_yticklabels([str(yr) + ('*' if yr in POOR_QUALITY_YEARS else '')
                          for yr in obs_years], fontsize=8.5, fontweight='bold')
    ax_c.set_ylabel('Imagery Year', fontsize=10, labelpad=6)
    ax_c.set_xticks([])
    ax_c.set_xlim(gis_ids[0]-0.5, gis_ids[-1]+0.5)
    ax_c.set_title(
        f'(c)  Model vs Observation Contingency  |  Qow threshold = {QOW_THRESHOLD} dam³/yr\n'
        f'Imagery years only ({n_obs} years with data)',
        fontsize=11, fontweight='bold', loc='left', pad=6)

    # Legend
    leg_handles = [
        mpatches.Patch(fc=CLR_HIT,  ec='#555', lw=0.5,
                       label='Hit — model & obs both show overwash'),
        mpatches.Patch(fc=CLR_MISS, ec='#555', lw=0.5,
                       label='Miss — obs overwash, model did not predict'),
        mpatches.Patch(fc=CLR_FA,   ec='#555', lw=0.5,
                       label='False Alarm — model predicted, obs shows none'),
        mpatches.Patch(fc=CLR_CR,   ec='#999', lw=0.5,
                       label='Correct Rejection — both show no overwash'),
        mpatches.Patch(fc=CLR_NA,   ec='#999', lw=0.5,
                       label='Not assessed (no imagery for that domain)'),
    ]
    ax_c.legend(handles=leg_handles, loc='upper right', fontsize=7.5,
                framealpha=0.93, edgecolor='#BBBBBB', ncol=2,
                bbox_to_anchor=(1.0, 1.02))

    # Section bar
    draw_section_bar(ax_b, list(gis_ids))
    fig.text(0.48, 0.175,
             'CASCADE Domain   (1 = Cape Point  →  90 = Pea Island / N. Rodanthe)',
             ha='center', va='top', fontsize=10)

    # Summary stats strip — POD and CSI by domain
    ax_s.set_facecolor('#FAFAFA')
    ax_s.set_xlim(gis_ids[0]-0.5, gis_ids[-1]+0.5)
    ax_s.set_ylim(0, 1); ax_s.set_yticks([0, 0.5, 1])
    ax_s.set_yticklabels(['0', '0.5', '1'], fontsize=7)
    ax_s.set_xticks([]); ax_s.spines['top'].set_visible(False)
    ax_s.set_ylabel('Score\nby domain', fontsize=7, labelpad=4)

    ax_s.bar(gis_ids, np.nan_to_num(pod, nan=0),
             width=0.8, color='#2E7D32', alpha=0.7, label='POD', zorder=2)
    ax_s.plot(gis_ids, np.nan_to_num(csi, nan=0),
              'o-', color='#C0392B', markersize=2.5, linewidth=0.8,
              label='CSI', zorder=3, alpha=0.85)
    ax_s.axhline(0.5, color='#AAAAAA', lw=0.6, ls='--', zorder=1)
    ax_s.legend(fontsize=6.5, loc='upper right', framealpha=0.9, ncol=2)

    # Overall title
    fig.text(0.48, 0.97,
             'CASCADE vs Observed — Contingency Analysis  |  Hatteras Island 1984–2003',
             ha='center', va='top', fontsize=12, fontweight='bold')
    fig.text(0.07, 0.015,
             f'* poor image quality  |  POD = Probability of Detection = Hits/(Hits+Misses)  |  '
             f'CSI = Critical Success Index = Hits/(Hits+Misses+False Alarms)  |  '
             f'Model: {os.path.basename(os.path.dirname(NPZ_PATH))}',
             fontsize=6.5, color='#555555', va='bottom')

    out = os.path.join(OUT_DIR, f'compare_contingency_T{QOW_THRESHOLD}.png')
    plt.savefig(out, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f"Saved: {out}")
    plt.show()

# ══════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════
def main():
    print("Loading CASCADE model comparison ...")
    years_m, Q, gis_ids = load_qow_matrix(NPZ_PATH)
    print(f"  Model: {years_m[0]}–{years_m[-1]}, {Q.shape[1]} domains, "
          f"Qow max = {Q.max():.1f} dam³/yr")

    print("Loading observation data ...")
    obs_mat, obs_years, domain_ints = load_obs_matrix(OBS_XLSX_PATH, years_m)
    print(f"  Observations: {len(obs_years)} imagery years within model range: "
          f"{list(obs_years)}")

    if PLOT_STACKED:
        print("\nGenerating stacked comparison figure ...")
        plot_stacked(years_m, Q, gis_ids, obs_mat, obs_years, domain_ints)

    if PLOT_CONTINGENCY:
        print("\nGenerating contingency figure ...")
        plot_contingency(years_m, Q, gis_ids, obs_mat, obs_years, domain_ints)

    print("\nDone.")

if __name__ == "__main__":
    main()
