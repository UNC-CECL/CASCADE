"""
Automorph (Itzkin, Moore, Ruggiero & Hacker 2020, "The effect of sand fencing on the morphology of natural dune
systems", Geomorphology 352, 106995), ported from Michael Itzkin's MATLAB (automorph_matlab/, a CECL tool, added with
his permission) and run on our 1996 transect profiles and shoreline (HAT_toe_method_line.run). Until 2026-11-08 this
script was a rebuild from the paper alone; its tables are kept in 2-automorph/rebuilt-from-paper-until-2026-11-08/.

    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_method_automorph.py

The versions:
    as sent (no prefix)  the MATLAB's automated pass, automorph.m with its first-year crest finder
    paper_               the same, plus what the paper describes and the code does not apply: the 35 deg building
                         check (Stockdon et al. 2009; the code computes the slope and never uses it) and a 3 m crest
                         floor (2.5 m in the code)
    adapted_             adapted for Hatteras (2026-11-08), as sent with three changes; step1_ and step12_ hold the
                         first and the first two, to see what each does:
                         1. the crest is your reviewed crest (crest_picks_1996.json), moved to the detected peak
                            within AM_SNAP_M (crest_finder_previous.m), no toe where you rejected it, in place of
                            his automated crest; his workflow reviews crests by hand too
                         2. only a peak with the backshore drop behind it, a real dune, moves the toe (step 4 below)
                         6. no crest of 12 m or more (his cap), and the toe screen of the other methods: at least
                            TOE_MIN_M_NAVD88, and TOE_BELOW_CREST_M below the crest

Steps, each transect (his file in brackets):
1. Smooth with a moving mean over about 3 m (automorph.m: 11 points at his 3.36 points per metre; 3 points here).
2. Mean high water: our shoreline (0.36 m NAVD88), not his closest point to the contour, since our profiles run across
   the whole island. The search stops CREST_SEARCH_M landward of it, near his 85-135 m profiles.
3. The crest (crest_finder_mull_ruggiero_2014.m): peaks above 2.5 m and below 12 m, seaward first; the first with a
   0.6 m drop before the next taller peak (Mull & Ruggiero 2014); else the first with any drop; else the highest
   point. The most landward peak is never tested (his loop stops one short). Then a seaward peak within 1% of the
   crest's height takes over, repeatedly. The code intends a 1.0 m drop from 2010 on but its test
   (num2str(year) < 2010) is always true, so 0.6 m is what it runs.
4. The toe (toe_finder.m): the point furthest below the straight line from the crest to MHW; then downhill landward
   to the local low; then, if a peak of step 3 stands between crest and toe, the lowest point from the crest to the
   most seaward such peak.
5. The berm (berm_finder.m): the tallest peak of the profile above the straight line from MHW to the toe (or to the
   2 m contour, where the toe is landward of it).
6. The toe again, from the crest-to-berm line (toe_finder.m), where a berm was found.
Paper version, between steps 4 and 5: where the crest-to-toe slope is steeper than 35 deg the crest moves to the next
peak seaward and the toe is found again; no peak seaward, no toe.
Left out: his year-to-year crest tracking (crest_finder_previous.m, for the years after the first), the fence steps,
and the manual review (profileReview.m, automorphCorrections*.m).
A toe needs the profile complete from MHW to the crest and a point below the line. Per domain, the median of the kept
transects when at least 3.

Writes, in data/hatteras_init/9-dune-toes/1-transects-1996/3-methods/2-automorph/ (D:):
automorph_transects_1996.csv, automorph_domains_1996.csv and automorph_method_1996.png (the method figure).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-08
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

sys.path.insert(0, str(Path(__file__).resolve().parent))
from HAT_toe_common import (CREST_SEARCH_M, EXAMPLES, FILL, MHW, MIN_KEPT, M_AUTOMORPH, PROFILE,  # noqa: E402
                            STEP_M, TOE, TOE_BELOW_CREST_M, TOE_MIN_M_NAVD88, WATER, YEAR, crest_picks)
from HAT_toe_method_line import run  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
AM_SMOOTH_POINTS = 3       # moving mean: his 11 points at 3.36 points per metre, about 3 m, at our 1 m spacing
AM_DROP_M = 0.6            # backshore drop behind the crest (Mull & Ruggiero 2014)
AM_CREST_MIN_M = 2.5       # peaks above this
AM_CREST_MAX_M = 12.0      # ... and below this (not a building)
AM_BUFFER = 0.01           # a seaward peak within this share of the crest's height takes over
AM_SNAP_M = 5.0            # adapted: your crest moves to the detected peak within this distance
AM_TWO_M = 2.0             # the contour that limits the berm search
PAPER_CREST_MIN_M = 3.0    # paper version: the crest floor
PAPER_MAX_SLOPE_DEG = 35.0  # paper version: a steeper crest-to-toe slope is a building (Stockdon et al. 2009)
AM = "#7B3F99"             # Automorph's colour in every figure
AM_PAPER = "#B48CC8"       # the paper version
AM_ADAPTED = "#3B1566"     # the version adapted for Hatteras
VERSIONS = {               # column prefix: the options of automorph() and its label
    "": (dict(), "as sent"),
    "paper_": (dict(paper=True), "with the paper's checks"),
    "step1_": (dict(picks=True), "+ 1 your crests"),
    "step12_": (dict(picks=True, real_peaks=True), "+ 2 only real dunes move the toe"),
    "adapted_": (dict(picks=True, real_peaks=True, screen=True), "+ 6 cap and toe screen (adapted)"),
}
# -----------------------------------------------------------------------------


def smooth(z):
    """The moving mean over surveyed points; gaps stay gaps."""
    m = pd.Series(z).rolling(AM_SMOOTH_POINTS, center=True, min_periods=1).mean().to_numpy()
    return np.where(np.isfinite(z), m, np.nan)


def local_peaks(a):
    return find_peaks(np.where(np.isfinite(a), a, -9999.0))[0]


def _peak_loop(z, p, end, threshold):
    """The first peak, seaward first, with threshold of drop before the next taller peak (or the profile's end)."""
    for k, i in enumerate(p[:-1]):
        taller = next((j for j in p[k + 1:] if z[j] > z[i]), end - 1)
        if z[i] - np.nanmin(z[i:taller + 1]) >= threshold:
            return i
    return None


def find_crest(z, shore, end, floor):
    """The crest (index) and the peaks that passed the height limits, seaward first."""
    p = shore + local_peaks(z[shore:end])
    p = np.array([i for i in p if floor < z[i] < AM_CREST_MAX_M], dtype=int)
    high = shore + int(np.nanargmax(z[shore:end]))
    if not len(p):
        return high, p
    c = _peak_loop(z, p, end, AM_DROP_M)
    if c is None:
        c = _peak_loop(z, p, end, 0.0)
    if c is None:
        return high, p
    while True:
        sea = [i for i in p if i < c and (1 - AM_BUFFER) * z[c] <= z[i] < z[c]]
        if not sea:
            return int(c), p
        c = sea[-1]


def find_toe(z, crest, anchor, p):
    """The toe on the crest-to-anchor line (anchor seaward of the crest); None when nothing lies below the line."""
    sag = np.linspace(z[anchor], z[crest], crest - anchor + 1) - z[anchor:crest + 1]
    if not np.nanmax(sag) > 0:
        return None
    t = anchor + int(np.nanargmax(sag))
    while t + 1 < len(z) and z[t + 1] < z[t]:                    # downhill landward: not on a high point
        t += 1
    between = [i for i in p if t < i < crest]                    # a peak between: the low behind the most seaward
    if between:
        t = between[0] + int(np.nanargmin(z[between[0]:crest + 1]))
    return t


def find_berm(z, shore, end, toe):
    """The tallest peak above the MHW-to-toe line (to the 2 m contour where the toe is landward of it), or None."""
    two = shore
    while two + 1 < end and not (z[two + 1] > AM_TWO_M and z[two] < AM_TWO_M):
        two += 1
    to = toe if toe < two else two
    if to <= shore:
        return None
    above = z[shore:to + 1] - np.linspace(z[shore], z[to], to - shore + 1)
    pk = local_peaks(above)
    return int(shore + pk[np.argmax(above[pk])]) if len(pk) else None


def auto_crest(z_raw, shore):
    """His automated crest, with the 3 m floor, on one profile (index): what HAT_toe_automorph_review.py suggests."""
    z = smooth(z_raw)
    end = min(len(z), shore + int(CREST_SEARCH_M / STEP_M) + 1)
    return None if not np.isfinite(z[shore:end]).any() else find_crest(z, shore, end, PAPER_CREST_MIN_M)[0]


def real_dunes(z, p, end):
    """The peaks with the backshore drop before the next taller peak (or the profile's end)."""
    keep = []
    for k, i in enumerate(p):
        taller = next((j for j in p[k + 1:] if z[j] > z[i]), end - 1)
        if z[i] - np.nanmin(z[i:taller + 1]) >= AM_DROP_M:
            keep.append(i)
    return np.array(keep, dtype=int)


def automorph(s, z_raw, shore, paper=False, pick=None, picks=False, real_peaks=False, screen=False):
    """Crest, first toe, berm and final toe (indices) on one profile, the peaks, the crest-to-toe slope (deg), whether
    the building check moved the crest, and the reason when there is no toe. pick: the transect's entry in
    crest_picks_1996.json, used with picks=True."""
    out = dict(crest=None, first_toe=None, berm=None, toe=None, peaks=[], moved=False, slope=np.nan, reason="")
    if shore is None:
        return out | dict(reason="no shoreline")
    if picks and (pick or {}).get("crest_x") is None:
        return out | dict(reason="crest rejected in the picker" if (pick or {}).get("status") == "rejected"
                          else "no crest pick")
    z = smooth(z_raw)
    end = min(len(z), shore + int(CREST_SEARCH_M / STEP_M) + 1)
    if picks:
        hand = min(len(z) - 1, shore + int(round(pick["crest_x"] / STEP_M)))
        end = max(end, min(len(z), hand + 1))
    if not np.isfinite(z[shore:end]).any():
        return out | dict(reason="no data behind the shoreline")
    crest, p = find_crest(z, shore, end, PAPER_CREST_MIN_M if paper else AM_CREST_MIN_M)
    if picks:                                                    # the detected peak nearest your crest, as his
        near = shore + local_peaks(z[shore:end])                 # crest_finder_previous does with a known crest
        near = near[np.abs(near - hand) <= AM_SNAP_M / STEP_M]
        crest = int(near[np.argmin(np.abs(near - hand))]) if len(near) else hand
    if screen and z[crest] >= AM_CREST_MAX_M:
        return out | dict(crest=crest, reason=f"crest at or above {AM_CREST_MAX_M:g} m")
    if real_peaks:
        p = real_dunes(z, p, end)
    out |= dict(crest=crest, peaks=[int(i) for i in p])
    while True:
        if crest - shore < 2:
            return out | dict(crest=crest, reason="crest at the shoreline")
        if not np.isfinite(z[shore:crest + 1]).all():
            return out | dict(crest=crest, reason="gap between shoreline and crest")
        toe = find_toe(z, crest, shore, p)
        if toe is None:
            return out | dict(crest=crest, reason="nothing below the crest-to-MHW line")
        slope = float(np.degrees(np.arctan2(z[crest] - z[toe], abs(s[crest] - s[toe]) or STEP_M)))
        if not paper or slope <= PAPER_MAX_SLOPE_DEG:
            break
        sea = [i for i in p if i < crest]
        if not sea:
            return out | dict(crest=crest, slope=slope, reason="steeper than 35 deg, no peak seaward")
        crest, out["moved"] = int(sea[-1]), True
    out |= dict(crest=crest, first_toe=toe, slope=slope)
    berm = find_berm(z, shore, end, toe)
    final = toe if berm is None or berm >= crest else find_toe(z, crest, berm, p)
    if final is None:
        return out | dict(berm=berm, reason="nothing below the crest-to-berm line")
    slope = float(np.degrees(np.arctan2(z[crest] - z[final], abs(s[crest] - s[final]) or STEP_M)))
    if screen and z[final] < TOE_MIN_M_NAVD88:
        return out | dict(berm=berm, slope=slope, reason=f"toe below {TOE_MIN_M_NAVD88:g} m")
    if screen and z[crest] - z[final] < TOE_BELOW_CREST_M:
        return out | dict(berm=berm, slope=slope, reason=f"crest under {TOE_BELOW_CREST_M:g} m above toe")
    return out | dict(berm=berm, toe=final, slope=slope)


def toes(profiles=None):
    """Automorph on every transect, in every version; elevations from the smoothed profile."""
    if profiles is None:
        profiles = run(write=False)[2]
    cp = crest_picks()
    rows = []
    for tid, (s, z, p, *_) in profiles.items():
        sh = p["shore"]
        zs = smooth(z)
        row = dict(transect=tid)
        for pre, (opts, _) in VERSIONS.items():
            paper = opts.get("paper", False)
            a = automorph(s, z, sh, pick=cp.get(str(tid)), **opts)
            ok, hc = a["toe"] is not None, a["crest"] is not None
            row |= {f"{pre}toe_m_navd88": float(zs[a["toe"]]) if ok else np.nan,
                    f"{pre}toe_from_shoreline_m": float(s[a["toe"]] - s[sh]) if ok else np.nan,
                    f"{pre}crest_m_navd88": float(zs[a["crest"]]) if hc else np.nan,
                    f"{pre}crest_from_shoreline_m": float(s[a["crest"]] - s[sh]) if hc else np.nan,
                    f"{pre}berm_m_navd88": float(zs[a["berm"]]) if a["berm"] is not None else np.nan,
                    f"{pre}slope_deg": a["slope"], f"{pre}reason": a["reason"] or "kept"}
            if paper:
                row["paper_building_check_moved_crest"] = a["moved"]
        rows.append(row)
    return pd.DataFrame(rows)


def domains(t):
    med = lambda col: t.groupby("gis")[col].median().where(t.groupby("gis")[col].count() >= MIN_KEPT)
    d = pd.DataFrame({f"{pre}toe_m_navd88": med(f"{pre}toe_m_navd88") for pre in VERSIONS})
    for pre in VERSIONS:
        d[f"{pre}n_kept"] = t.groupby("gis")[f"{pre}toe_m_navd88"].count()
    d.index.name = "gis"
    return d.reindex(range(1, 91))


def method_figure(profiles, gis_of):
    """The method in three steps on the two example transects (columns): the crest, the first toe and the berm, and
    the final toes of the two versions."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    fig, axes = plt.subplots(3, 2, figsize=(fs.figsize("double")[0], 7.4), layout="constrained", sharey="col")
    for col, tid in enumerate(EXAMPLES):
        s, z_raw, p = profiles[tid][:3]
        sh = p["shore"]
        z = smooth(z_raw)
        a, b = automorph(s, z_raw, sh), automorph(s, z_raw, sh, paper=True)
        c = automorph(s, z_raw, sh, pick=crest_picks().get(str(tid)), **VERSIONS["adapted_"][0])
        x = s - s[sh]
        xmax = max(140.0, (x[a["crest"]] if a["crest"] is not None else 0) + 40)
        v = (x >= -20) & (x <= xmax)
        for row, ax in enumerate(axes[:, col]):
            ax.fill_between(x[v], z[v], -1, color=FILL, lw=0)
            ax.plot(x[v], z[v], color=PROFILE, lw=1.0)
            ax.axhline(MHW, color=WATER, lw=0.6)
            ax.plot(0, z[sh], "o", ms=5, color=WATER, zorder=6, label="MHW" if row == 0 else None)
            ax.set_xlim(-20, xmax)
            ax.set_ylim(-1, np.nanmax(z[v]) + 1.2)
            fs.open_frame(ax)
            ax.text(-0.02, 1.02, "abcdef"[2 * row + col], transform=ax.transAxes, fontsize=10, fontweight="bold",
                    ha="right", va="bottom")
            if col == 0:
                ax.set_ylabel("Elevation (m NAVD88)")
            if row == 2:
                ax.set_xlabel("Distance landward of the shoreline (m)")
        ci, t1, bi, ti = a["crest"], a["first_toe"], a["berm"], a["toe"]
        ax = axes[0, col]                                        # 1: the crest
        for i in a["peaks"]:
            ax.plot(x[i], z[i], "^", ms=7 if i == ci else 5.5, mfc=AM if i == ci else "white", mec=AM, mew=0.9, zorder=6)
        ax.text(0.02, 0.97, f"Transect {tid} (GIS {gis_of[tid]})", transform=ax.transAxes, fontsize=8, va="top",
                color=PROFILE)
        ax = axes[1, col]                                        # 2: the first toe and the berm
        if ci is not None and t1 is not None:
            ax.plot([0, x[ci]], [z[sh], z[ci]], color=AM, lw=0.9, ls=(0, (4, 2)))
            ax.plot(x[ci], z[ci], "^", ms=7, color=AM, zorder=6)
            ax.plot(x[t1], z[t1], "s", ms=6, mfc="white", mec=AM, mew=1.0, zorder=7)
            if bi is not None:
                ax.plot(x[bi], z[bi], "D", ms=5, color=WATER, zorder=7)
        ax = axes[2, col]                                        # 3: the final toes
        if ci is not None and ti is not None:
            anchor = sh if bi is None else bi
            ax.plot([x[anchor], x[ci]], [z[anchor], z[ci]], color=AM, lw=0.9, ls=(0, (4, 2)))
            ax.plot(x[ci], z[ci], "^", ms=7, color=AM, zorder=6)
            ax.plot(x[ti], z[ti], "s", ms=7, color=AM, mec="white", mew=0.6, zorder=7)
            ax.text(x[ti] + 3, z[ti] - 0.45, f"toe {z[ti]:.2f} m", fontsize=7, color=AM)
        if b["toe"] is not None:
            ax.plot(x[b["toe"]], z[b["toe"]], "o", ms=9, mfc="none", mec=AM_PAPER, mew=1.3, zorder=8)
        if c["toe"] is not None:
            ax.plot(x[c["toe"]], z[c["toe"]], "P", ms=8, color=AM_ADAPTED, mec="white", mew=0.6, zorder=9)
        if p["toe"] is not None:                              # ours, for reference
            ax.plot(x[p["toe"]], z_raw[p["toe"]], "o", ms=6, mfc="none", mec=TOE, mew=1.0, zorder=8)
    for row, txt in enumerate(("1. Crest: most seaward peak, 2.5-12 m, with a 0.6 m drop behind",
                               "2. First toe under the crest-to-MHW line, then the berm",
                               "3. Final toe under the crest-to-berm line")):
        axes[row, 0].text(0.0, 1.12, txt, transform=axes[row, 0].transAxes, fontsize=8, fontweight="bold")
    h = [plt.Line2D([], [], lw=0, marker="^", ms=6, mfc="white", mec=AM, label="Other peak within the height limits"),
         plt.Line2D([], [], lw=0, marker="^", ms=7, color=AM, label="Automorph crest"),
         plt.Line2D([], [], lw=0, marker="s", ms=6, mfc="white", mec=AM, label="First toe"),
         plt.Line2D([], [], lw=0, marker="D", ms=5, color=WATER, label="Berm"),
         plt.Line2D([], [], lw=0, marker="s", ms=6, color=AM, label="Final toe, code as sent"),
         plt.Line2D([], [], lw=0, marker="o", ms=8, mfc="none", mec=AM_PAPER, label="Final toe, with the paper's checks"),
         plt.Line2D([], [], lw=0, marker="P", ms=8, color=AM_ADAPTED, label="Final toe, adapted for Hatteras"),
         plt.Line2D([], [], lw=0, marker="o", ms=6, mfc="none", mec=TOE, label="Our toe")]
    fig.legend(handles=h, loc="outside lower center", ncol=4, fontsize=7, frameon=False)
    png = M_AUTOMORPH / f"automorph_method_{YEAR}.png"
    fs.save(fig, png, close=True)
    fs.record_caption(png, (
        f"Automorph (Itzkin et al. 2020), ported from the authors' MATLAB code, on two {YEAR} lidar transects: "
        f"{EXAMPLES[0]} (left, a single foredune) and {EXAMPLES[1]} (right, a dune field), profiles smoothed over "
        f"about 3 m. (a, b) The crest: the most seaward peak between {AM_CREST_MIN_M:g} and {AM_CREST_MAX_M:g} m with "
        f"a {AM_DROP_M:g} m drop behind it (Mull and Ruggiero 2014). (c, d) The first toe, the point furthest below "
        f"the straight line from the crest to mean high water (after Mitasova et al. 2011), and the berm. (e, f) The "
        f"final toe, furthest below the line from the crest to the berm; open circles, the same with the paper's "
        f"{PAPER_MAX_SLOPE_DEG:.0f}° building check and {PAPER_CREST_MIN_M:g} m crest floor, and our toe; cross, the "
        f"version adapted for Hatteras (hand-reviewed crests, only peaks with the backshore drop move the toe, the "
        f"{AM_CREST_MAX_M:g} m cap and the toe screen)."))
    print(f"wrote {png}")


def main():
    t0, dom0, profiles = run(write=False)[:3]
    t = toes(profiles)
    t["gis"] = t.transect.map(t0.set_index("transect").gis)
    t = t[["transect", "gis"] + [c for c in t.columns if c not in ("transect", "gis")]]
    d = domains(t)
    M_AUTOMORPH.mkdir(parents=True, exist_ok=True)
    t.round(3).to_csv(M_AUTOMORPH / f"automorph_transects_{YEAR}.csv", index=False)
    d.round(3).to_csv(M_AUTOMORPH / f"automorph_domains_{YEAR}.csv")
    ours = dom0.set_index("gis").toe_m_navd88
    current = dom0.set_index("gis").current_toe_m_navd88
    rmse = lambda a, b: float(np.sqrt(((a - b).dropna() ** 2).mean()))
    for pre, (_, name) in VERSIONS.items():
        tz, dz = t[f"{pre}toe_m_navd88"], d[f"{pre}toe_m_navd88"]
        print(f"Automorph {name}: {int(tz.notna().sum())}/{len(t)} transects with a toe, {int(dz.notna().sum())} "
              f"domain medians, median toe {tz.median():.2f} m, lowest {tz.min():.2f} m, below 1.7 m {int((tz < 1.7).sum())}, "
              f"highest {tz.max():.2f} m; domains vs ours RMSE {rmse(dz, ours):.2f} m, > 0.5 m "
              f"{int(((dz - ours).abs() > 0.5).sum())}; vs current RMSE {rmse(dz, current):.2f} m")
        print("   " + t[f"{pre}reason"].value_counts().to_string().replace("\n", "\n   "))
    print(f"building check moved {int(t.paper_building_check_moved_crest.sum())} crests")
    print(f"wrote {M_AUTOMORPH}")
    method_figure(profiles, t.set_index("transect").gis.to_dict())


if __name__ == "__main__":
    main()
