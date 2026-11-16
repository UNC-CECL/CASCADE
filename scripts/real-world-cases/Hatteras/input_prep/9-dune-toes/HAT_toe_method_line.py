"""
Our dune toe, the line method (2026-11-05): one toe per GIS domain from the 100 m shore-normal transects and the
1996 lidar alone, the survey the model starts from. A simpler alternative to HAT_dune_toes.py (10 m model cells, 1996
and 1997 averaged). The toe rule (the point furthest below the shoreline-to-crest line) is the perpendicular-distance
method of pybeach (Beuzen 2019) and Automorph's (Itzkin et al. 2020, after Mitasova et al. 2011); ours adds the crest
inside the model's dune search window and the berm floor.

    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_method_line.py [--gis 22]
    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_method_line.py --overview

Steps, each transect:
1. The transects: the 100 m transects of the BRIE offset step (2-brie-offset/transects/transects_100m.geojson):
   east-west lines 10 km long from the ocean across the island, on the model domains' frame, 5 per 500 m domain.
2. The profile: the 1996 lidar (NOAA, 9 ft grid, ft NAVD88) read every 1 m along it (bilinear), in m NAVD88, with
   the survey's vertical correction to the start DEM (-0.24 m, HAT_survey_dune_crests.SURVEY_VERTICAL_OFFSET_M).
3. The shoreline: the most seaward point at or above mean high water (0.36 m NAVD88) with the next
   SHORE_RUN_M all above it.
4. The front dune's crest, within CREST_SEARCH_M landward of the shoreline: the most seaward peak with a dip of
   FOREDUNE_DROP_M (1 m) behind it, as HAT_dune_toes.py; the highest point where no peak has one. Without it, on a
   field of two or more ridges the back ridge is picked and the toe falls in the swale between them.
   Second version, compared (window_ columns, 2026-11-05): the crest is the highest point of the model's dune search
   window on the transect (Windows: the window the model's crest is taken from, mapped onto the transect through
   the domain's 1 m clip grid; checked against survey_dune_crests.csv). Where a window spans two ridges it takes
   the back one (GIS 22, 35, 65-67, 70: the two toes > 0.5 m apart).
5. The toe: the point furthest below the straight line from the shoreline to the crest.
6. The screen: a toe is kept when the profile is complete from the shoreline to the crest, the toe is at least
   TOE_MIN_M_NAVD88 and at least TOE_BELOW_CREST_M below the crest (as HAT_dune_toes.py).
7. The domain toe: the median of its kept transects, when at least MIN_KEPT; otherwise the mean of the nearest
   domains each side that have one (flagged "filled").

Writes, in data/hatteras_init/9-dune-toes/1-transects-1996/ (D:): 1-toes/ transect_toes_1996.csv (one row per
transect) and domain_toes_1996.csv (one row per domain); 3-methods/1-line/ one figure per step for an example domain
(the method figure), all steps as one lettered figure, the alongshore result, and with --overview the island map.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-05
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parent))
from HAT_toe_common import (BEND_M, BEND_SHARE, BERM_M, CHORD, CREST, CREST_SEARCH_M, CURRENT, DROPPED,
    EXAMPLE_GIS, FILL, FOREDUNE_DROP_M, MHW, MIN_KEPT, M_LINE, PEAK_HALF_M, PROFILE, S, SHORE_RUN_M, START, STEP_M,
    T, TOE, TOES, TOE_BELOW_CREST_M, TOE_MIN_M_NAVD88, TRANSECTS, VIEW_M, WATER, WINDOW, Windows, YEAR, _legend,
    apply_qc, crest_picks, first_dune, line_points, sample, survey, transects)  # noqa: E402



def pick(s, z, window=None, drop=None, mode="window_front", crest_x=None):
    """Shoreline, crest and toe on one profile (indices), and why a toe is not kept.
    mode "window_front" (used): the crest is the front dune (first_dune) searched inside the model's dune search
    window (window: its s_seaward, s_landward), landward of the shoreline; with no window, the front dune within
    CREST_SEARCH_M ("crest_fallback"). mode "front": the front dune within CREST_SEARCH_M. mode "window_max": the
    window's highest point, the model's own crest. crest_x: a crest you picked on this transect (m landward of the
    shoreline, crest_picker()), used in place of the search.
    The toe is the point furthest below the shoreline-to-crest line, no lower than the berm: where that point is
    below BERM_M it moves landward to where the profile first reaches the berm ("floored"). "alt": the other bend,
    where a second point almost as far below the line (BEND_SHARE of the deepest) gives a toe BEND_M or more apart."""
    run = int(SHORE_RUN_M / STEP_M)
    above = np.isfinite(z) & (z >= MHW)
    shore = next((i for i in range(len(z) - run) if above[i:i + run + 1].all()), None)
    out = dict(shore=shore, crest=None, toe=None, floored=False, alt=None, crest_fallback=False)
    if shore is None:
        return out | dict(reason="no shoreline")
    if crest_x is not None:
        crest = min(len(z) - 1, shore + int(round(crest_x / STEP_M)))
    elif mode != "front" and window is not None:
        lo = max(shore, int(np.searchsorted(s, window[0])))
        hi = min(len(z), int(np.searchsorted(s, window[1], side="right")))
        if hi <= lo:
            return out | dict(reason="window seaward of the shoreline")
        seg = z[lo:hi]
        if not np.isfinite(seg).all():
            return out | dict(reason="window not surveyed")
        crest = lo + (int(np.argmax(seg)) if mode == "window_max" else first_dune(seg, drop))
    elif mode == "window_max":
        return out | dict(reason="no window")
    else:
        out["crest_fallback"] = mode == "window_front"
        hi = min(len(z), shore + int(CREST_SEARCH_M / STEP_M) + 1)
        seg = z[shore:hi]
        if not np.isfinite(seg).any():
            return out | dict(reason="no data behind the shoreline")
        crest = shore + first_dune(seg, drop)
    out["crest"] = crest
    if crest - shore < 2:
        return out | dict(reason="crest at the shoreline")
    prof = z[shore:crest + 1]
    sag = np.linspace(z[shore], z[crest], len(prof)) - prof
    toe = shore + int(np.nanargmax(sag))
    floor = lambda i: i if z[i] >= BERM_M else next((j for j in range(i, crest + 1) if z[j] >= BERM_M), crest)
    w = int(PEAK_HALF_M / STEP_M)                                # the two-bends check, on the floored toes
    sg = np.where(np.isfinite(sag), sag, -np.inf)
    for i in range(len(sg)):
        if sg[i] < BEND_SHARE * sg.max() or sg[i] < sg[max(0, i - w):i + w + 1].max():
            continue
        if abs(z[floor(shore + i)] - z[floor(toe)]) >= BEND_M:
            out["alt"] = floor(shore + i)
            break
    out["floored"] = bool(z[toe] < BERM_M)
    toe = floor(toe)
    reason = ""
    if not np.isfinite(prof).all():
        reason = "gap between shoreline and crest"
    elif z[toe] < TOE_MIN_M_NAVD88:
        reason = f"toe below {TOE_MIN_M_NAVD88:g} m"
    elif z[crest] - z[toe] < TOE_BELOW_CREST_M:
        reason = f"crest under {TOE_BELOW_CREST_M:g} m above toe"
    return out | dict(toe=toe, reason=reason)


def _row(s, z, p, prefix=""):
    g = lambda k: (float(s[p[k]]), float(z[p[k]])) if p[k] is not None else (np.nan, np.nan)
    (ss, _), (sc, zc), (st, zt) = g("shore"), g("crest"), g("toe")
    return {f"{prefix}crest_s_m": sc, f"{prefix}crest_m_navd88": zc, f"{prefix}toe_s_m": st,
            f"{prefix}toe_m_navd88": zt, f"{prefix}toe_from_shoreline_m": st - ss,
            f"{prefix}kept": p["reason"] == "", f"{prefix}reason": p["reason"] or "kept",
            f"{prefix}floored_to_berm": bool(p.get("floored")),
            f"{prefix}other_bend_toe_m_navd88": float(z[p["alt"]]) if p.get("alt") is not None else np.nan,
            f"{prefix}crest_fallback": bool(p.get("crest_fallback"))}


def _medians(t, prefix=""):
    """One toe per domain: the median of the kept transects (MIN_KEPT or more), else the nearest measured domains'
    mean each side."""
    k = t[t[f"{prefix}kept"]]
    agg = k.groupby("gis")[f"{prefix}toe_m_navd88"].agg(["median", "count"])
    d = pd.DataFrame({"gis": range(1, 91)}).join(agg, on="gis")
    d["count"] = d["count"].fillna(0).astype(int)
    d.loc[d["count"] < MIN_KEPT, "median"] = np.nan
    src = np.where(d["median"].notna(), "median", "filled")
    have = d.dropna(subset=["median"]).set_index("gis")["median"]
    for i in d.index[d["median"].isna()]:
        g = d.at[i, "gis"]
        lo, hi = have[have.index < g], have[have.index > g]
        d.at[i, "median"] = float(np.mean(([lo.iloc[-1]] if len(lo) else []) + ([hi.iloc[0]] if len(hi) else [])))
    return d["median"].to_numpy(), d["count"].to_numpy(), src


def run(write=True):
    src, scale, offset = survey()
    tr, domains = transects()
    win = Windows()
    cp = crest_picks()                                           # your crests, where you picked them
    rows, profiles = [], {}
    for r in tr.itertuples():
        s, xy, p0, u = line_points(r.geometry)
        z = sample(src, scale, offset, xy, tr.crs)
        wa = win.along(r.gis, p0, u, float(np.mean(np.asarray(r.geometry.coords)[:, 1])))
        e = cp.get(str(r.ObjectID), {})
        p = pick(s, z, wa[1:] if wa else None, crest_x=e.get("crest_x"))   # used: your crest, else the front dune
        if e.get("status") == "rejected":                                  # inside the model window
            p = p | dict(reason="rejected in the crest picker")
        pw = pick(s, z, wa[1:] if wa else None, mode="window_max")
        pf = pick(s, z, mode="front")
        profiles[r.ObjectID] = (s, z, p, xy, pw, wa, pf)
        ss = float(s[p["shore"]]) if p["shore"] is not None else np.nan
        rows.append(dict(transect=r.ObjectID, gis=r.gis, cell=wa[0] if wa else -1, shoreline_s_m=ss,
                         window_sea_s_m=wa[1] if wa else np.nan, window_land_s_m=wa[2] if wa else np.nan,
                         crest_source=("picked" if e.get("status") == "picked" else
                                       "rejected" if e.get("status") == "rejected" else
                                       "accepted" if e.get("status") == "accepted" else "automatic"),
                         **_row(s, z, p), **_row(s, z, pw, "window_"), **_row(s, z, pf, "free_")))
    t = apply_qc(pd.DataFrame(rows))                         # your picker checks, where there are any

    # Check the windows land on the model's dune: the transect's highest point in the window against the 1996
    # crest of the same 10 m cell in survey_dune_crests.csv (10 m cells average the centre 2 x 2 m, so a little lower)
    sc = pd.read_csv(S.tv.product_dir(START) / "survey-crests" / "survey_dune_crests.csv")
    sc = sc[(sc.survey == YEAR) & (sc.status == "ok")].set_index(["gis", "cell"]).crest_m_navd88
    t["model_cell_crest_m_navd88"] = [sc.get((g, c), np.nan) for g, c in zip(t.gis, t.cell)]
    dc = (t.window_crest_m_navd88 - t.model_cell_crest_m_navd88).dropna()
    print(f"window check, transect crest minus the cell's crest in survey_dune_crests.csv ({len(dc)} transects): "
          f"median {dc.median():+.2f} m, within 0.5 m {100 * (dc.abs() <= 0.5).mean():.0f}%")

    dom = pd.DataFrame({"gis": range(1, 91)})
    dom["toe_m_navd88"], dom["n_kept"], dom["source"] = _medians(t)
    dom["toe_window_m_navd88"], dom["n_kept_window"], dom["source_window"] = _medians(t, "window_")
    dom["toe_free_m_navd88"], dom["n_kept_free"], dom["source_free"] = _medians(t, "free_")
    dom["toe_qc_m_navd88"], dom["n_kept_qc"], dom["source_qc"] = _medians(t, "qc_")
    cur = T._tables()[START][1].groupby("gis").toe_m_navd88.first()
    dom["current_toe_m_navd88"] = dom.gis.map(cur)
    dom["toe_m_mhw"] = dom.toe_m_navd88 - MHW
    dom["toe_window_m_mhw"] = dom.toe_window_m_navd88 - MHW
    # where the window's crest is a different ridge from the front dune (the window picks a crest > 10 m behind it)
    t["window_behind_front_dune_m"] = t.window_crest_s_m - t.crest_s_m
    dom["n_window_behind_front_dune"] = dom.gis.map(
        t[t.window_behind_front_dune_m > 10].groupby("gis").size()).fillna(0).astype(int)
    if not write:                                           # the other methods read the profiles only
        return t, dom, profiles, tr, domains, src, scale, offset
    TOES.mkdir(parents=True, exist_ok=True)
    t.round(3).to_csv(TOES / f"transect_toes_{YEAR}.csv", index=False)
    dom.round(3).to_csv(TOES / f"domain_toes_{YEAR}.csv", index=False)
    print(f"wrote {TOES}")
    print(f"used method: {int(t.floored_to_berm.sum())} toes floored to the berm, {int(t.other_bend_toe_m_navd88.notna().sum())} "
          f"with a second bend >= {BEND_M:g} m apart, {int(t.crest_fallback.sum())} without a window (front dune "
          f"within {CREST_SEARCH_M:.0f} m)")
    for name, pre, col in (("front dune in the window (used)", "", "toe_m_navd88"),
                           ("window highest point", "window_", "toe_window_m_navd88"),
                           ("front dune within 120 m", "free_", "toe_free_m_navd88")):
        d = dom[col] - dom.current_toe_m_navd88
        print(f"{name}: kept {int(t[f'{pre}kept'].sum())}/{len(t)}; minus current: mean {d.mean():+.2f} m, "
              f"RMSE {np.sqrt((d ** 2).mean()):.2f} m")
        print("   " + t[f"{pre}reason"].value_counts().to_string().replace("\n", "\n   "))
    d = dom.toe_window_m_navd88 - dom.toe_m_navd88
    print(f"window minus front dune: mean {d.mean():+.2f} m, |diff| > 0.5 m in {int((d.abs() > 0.5).sum())} domains; "
          f"transects whose window crest is a ridge behind the front dune: "
          f"{int((t.window_behind_front_dune_m > 10).sum())}")
    return t, dom, profiles, tr, domains, src, scale, offset


class _Ctx:
    """Everything the panels share for the example domain."""

    def __init__(self, t, dom, profiles, tr, domains, src, scale, offset, gis):
        self.t, self.dom, self.profiles, self.tr, self.domains = t, dom, profiles, tr, domains
        self.src, self.scale, self.offset, self.gis = src, scale, offset, gis
        self.ids = list(tr[tr.gis == gis].ObjectID)
        self.ex = self.ids[len(self.ids) // 2]                   # the middle transect
        s, z, p, _ = profiles[self.ex][:4]
        self.s, self.z, self.p = s, z, p
        self.x = s - s[p["shore"]]                               # m landward of the shoreline
        v = (self.x >= VIEW_M[0]) & (self.x <= VIEW_M[1])
        self.v = v
        self.top = float(np.ceil(np.nanmax(z[v]) + 0.8))
        self.rows = t[t.gis == gis].set_index("transect")


def _profile(ax, c, fs):
    ax.fill_between(c.x[c.v], c.z[c.v], -1, color=FILL, lw=0, zorder=1)
    ax.plot(c.x[c.v], c.z[c.v], color=PROFILE, lw=1.1, zorder=4, label=f"Profile, transect {c.ex}")
    ax.axhline(MHW, color=WATER, lw=0.7, zorder=2)
    ax.text(-3, MHW + 0.12, "MHW", fontsize=7, color=WATER, ha="right")
    ax.axhline(BERM_M, color=CREST, lw=0.6, ls=(0, (2, 2)), zorder=2)
    ax.text(VIEW_M[0] + 2, BERM_M + 0.12, f"Berm ({BERM_M:g} m)", fontsize=7, color=CREST, ha="left")
    ax.set_xlim(*VIEW_M)
    ax.set_ylim(-0.5, c.top)
    ax.set_xlabel("Distance landward of the shoreline (m)")
    ax.set_ylabel("Elevation (m NAVD88)")
    fs.open_frame(ax)


def panel_map(ax, c, fs, fig):
    from matplotlib.colors import LightSource
    poly = c.domains[c.domains.domain_id == c.gis].geometry.iloc[0]
    x0, y0, x1, y1 = poly.bounds
    pts = np.vstack([c.profiles[i][3] for i in c.ids])
    ok = np.isfinite(np.concatenate([c.profiles[i][1] for i in c.ids]))
    e0, e1 = pts[ok, 0].min() - 60, pts[ok, 0].max() + 60
    gx, gy = np.meshgrid(np.arange(e0, e1, 3.0), np.arange(y1 + 40, y0 - 40, -3.0))
    gz = sample(c.src, c.scale, c.offset, np.c_[gx.ravel(), gy.ravel()], c.tr.crs).reshape(gx.shape)
    ext = (0, e1 - e0, y0 - 40 - y0, y1 + 40 - y0)
    hs = LightSource(azdeg=315, altdeg=45).hillshade(np.nan_to_num(gz, nan=0.0), vert_exag=3, dx=3, dy=3)
    ax.imshow(np.where(np.isfinite(gz), hs, np.nan), cmap="Greys_r", extent=ext, origin="upper", zorder=1,
              interpolation="bilinear")
    im = ax.imshow(gz, cmap="terrain", vmin=-1, vmax=8, alpha=0.8, extent=ext, origin="upper", zorder=2,
                   interpolation="bilinear")
    ax.plot(*(np.array(poly.exterior.coords) - [e0, y0]).T, color=fs.INK, lw=0.8, ls=(0, (4, 2)), zorder=3,
            label=f"GIS {c.gis} domain")
    for i in c.ids:
        _, _, pi, p, pw, wa = c.profiles[i][:6]
        ax.plot(p[:, 0] - e0, p[:, 1] - y0, color=TOE, lw=0.9, zorder=4)
        ax.text(6, p[0, 1] - y0 + 7, str(i), fontsize=6.5, color=TOE, zorder=5)
        if wa:                                                   # the model's dune search window on it
            k = (c.profiles[i][0] >= wa[1]) & (c.profiles[i][0] <= wa[2])
            ax.plot(p[k, 0] - e0, p[k, 1] - y0, color=WINDOW, lw=3.2, alpha=0.8, solid_capstyle="butt", zorder=4)
        if pi["toe"] is not None:
            ax.plot(p[pi["toe"], 0] - e0, p[pi["toe"], 1] - y0, "o", ms=4.5, mfc="white", mec=TOE, mew=1.0,
                    zorder=5)
    ax.set_xlim(0, e1 - e0)
    ax.set_ylim(ext[2], ext[3])
    ax.set_aspect("equal")
    ax.set_xlabel("Easting (m, local)")
    ax.set_ylabel("Northing (m, local)")
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02, aspect=25)
    cb.set_label("Elevation (m NAVD88)", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    cb.outline.set_linewidth(0.5)
    ax.text(0.985, 0.5, "Ocean", transform=ax.transAxes, ha="right", va="center", rotation=90, fontsize=7,
            style="italic", color="white")
    ax.text(6, y1 - y0 - 8, f"GIS {c.gis}", fontsize=7, ha="left", va="top", color=fs.INK)


def panel_profiles(ax, c, fs, fig=None):
    shades = __import__("matplotlib").pyplot.get_cmap("Blues")(np.linspace(0.45, 0.95, len(c.ids)))
    o = c.s[c.p["shore"]]
    for i, col in zip(c.ids, shades):
        si, zi, _, _ = c.profiles[i][:4]
        ax.plot(si - o, zi, color=col, lw=1.0, label=f"Transect {i}")
    ax.axhline(MHW, color=WATER, lw=0.7)
    ax.set_xlim(VIEW_M[0] - 20, VIEW_M[1] + 120)
    ax.set_ylim(-0.8, None)
    ax.set_xlabel(f"Distance landward of the shoreline of transect {c.ex} (m)")
    ax.set_ylabel("Elevation (m NAVD88)")
    fs.open_frame(ax)
    _legend(ax)


def panel_shoreline(ax, c, fs, fig=None):
    _profile(ax, c, fs)
    ax.axvline(0, color=WATER, lw=0.6, ls=(0, (2, 2)))
    ax.plot(0, c.z[c.p["shore"]], "o", ms=6, color=WATER, zorder=6, label="Shoreline (MHW)")
    _legend(ax)


def panel_crest(ax, c, fs, fig=None):
    _profile(ax, c, fs)
    s, wa = c.profiles[c.ex][0], c.profiles[c.ex][5]
    if wa:
        o = s[c.p["shore"]]
        ax.axvspan(wa[1] - o, wa[2] - o, color=WINDOW, alpha=0.12, lw=0, zorder=0, label="Model dune window")
    ax.plot(0, c.z[c.p["shore"]], "o", ms=5, color=WATER, zorder=6, label="Shoreline")
    cx = c.x[c.p["crest"]]
    ax.plot(cx, c.z[c.p["crest"]], "^", ms=7, color=CREST, zorder=6,
            label=f"Foredune crest ({c.z[c.p['crest']]:.2f} m)")
    _legend(ax)


def panel_toe(ax, c, fs, fig=None):
    _profile(ax, c, fs)
    zs, zc, zt = c.z[c.p["shore"]], c.z[c.p["crest"]], c.z[c.p["toe"]]
    cx, tx = c.x[c.p["crest"]], c.x[c.p["toe"]]
    ax.plot([0, cx], [zs, zc], color=CHORD, lw=0.8, ls=(0, (4, 2)), zorder=3, label="Shoreline–crest chord")
    at = zs + (zc - zs) * tx / cx
    ax.annotate("", (tx, zt), (tx, at), arrowprops=dict(arrowstyle="<->", color=TOE, lw=0.8,
                                                        shrinkA=0, shrinkB=0))
    ax.plot(0, zs, "o", ms=5, color=WATER, zorder=6, label="Shoreline")
    ax.plot(cx, zc, "^", ms=7, color=CREST, zorder=6, label="Foredune crest")
    i0, i1 = c.p["shore"], c.p["crest"]
    raw = i0 + int(np.nanargmax(np.linspace(zs, zc, i1 - i0 + 1) - c.z[i0:i1 + 1]))
    if c.p.get("floored"):
        rx, rz = c.x[raw], c.z[raw]
        ax.plot(rx, rz, "o", ms=6, mfc="white", mec=TOE, mew=1.0, zorder=7,
                label=f"Deepest below line ({rz:.2f} m)")
        ax.annotate("", (tx, zt), (rx, rz), arrowprops=dict(arrowstyle="->", color=TOE, lw=0.8))
    ax.plot(tx, zt, "o", ms=7, color=TOE, zorder=7,
            label=(f"Toe at berm ({zt:.2f} m)" if c.p.get("floored") else f"Dune toe ({zt:.2f} m)"))
    _legend(ax)


def panel_screen(ax, c, fs, fig=None):
    shades = __import__("matplotlib").pyplot.get_cmap("Blues")(np.linspace(0.45, 0.95, len(c.ids)))
    for i, col in zip(c.ids, shades):
        si, zi, pi, _ = c.profiles[i][:4]
        if pi["shore"] is None:
            continue
        r = c.rows.loc[i]
        x = si - si[pi["shore"]]
        ax.plot(x, zi, color=col if r.kept else DROPPED, lw=0.9, zorder=3)
        if pi["toe"] is not None:
            ax.plot(x[pi["toe"]], zi[pi["toe"]], "o" if r.kept else "x", ms=6 if r.kept else 7,
                    color=col if r.kept else "#B03A2E", mec="white" if r.kept else None, mew=0.6 if r.kept else 1.5,
                    zorder=6, label=f"{i}: {r.toe_m_navd88:.2f} m" + ("" if r.kept else f" ({r.reason})"))
    ax.axhline(BERM_M, color=CREST, lw=0.6, ls=(0, (2, 2)))
    ax.text(VIEW_M[0] + 2, BERM_M + 0.12, f"Berm, lowest toe ({BERM_M:g} m)", fontsize=7,
            color=CREST, ha="left")
    ax.axhline(MHW, color=WATER, lw=0.7)
    ax.set_xlim(*VIEW_M)
    ax.set_ylim(-0.5, c.top + 3)
    ax.set_xlabel("Distance landward of the shoreline (m)")
    ax.set_ylabel("Elevation (m NAVD88)")
    fs.open_frame(ax)
    _legend(ax, title="Transect: toe", title_fontsize=7)


def panel_domain(ax, c, fs, fig=None):
    r, d = c.rows, c.dom.set_index("gis").loc[c.gis]
    xs = np.arange(len(c.ids))
    k = r[r.kept]
    ax.scatter([xs[c.ids.index(i)] for i in k.index], k.toe_m_navd88, s=28, color=TOE, zorder=5,
               label="Transect toe (kept)")
    dr = r[~r.kept & r.toe_m_navd88.notna()]
    if len(dr):
        ax.scatter([xs[c.ids.index(i)] for i in dr.index], dr.toe_m_navd88, s=30, marker="x", color=DROPPED,
                   zorder=4, label="Transect toe (rejected)")
    ax.axhline(d.toe_m_navd88, color=TOE, lw=1.4, label=f"Median ({d.toe_m_navd88:.2f} m)")
    mean = float(k.toe_m_navd88.mean())
    ax.axhline(mean, color="#6BAED6", lw=1.2, ls=(0, (4, 2)), label=f"Mean ({mean:.2f} m)")
    vals = r.toe_m_navd88.dropna().tolist()
    ax.set_ylim(min(vals) - 0.4, max(vals) + 0.9)
    ax.set_xlim(-0.5, len(c.ids) - 0.5)
    ax.set_xticks(xs, [str(i) for i in c.ids])
    ax.set_xlabel(f"Transect (GIS {c.gis}, south → north)")
    ax.set_ylabel("Dune toe (m NAVD88)")
    fs.open_frame(ax)
    _legend(ax, ncol=2)


def panel_alongshore(ax, c, fs, fig=None):
    t, dom = c.t, c.dom
    kk = t[t.kept]
    jit = (kk.groupby("gis").cumcount() - 2) * 0.12
    ax.scatter(kk.gis + jit, kk.toe_m_navd88, s=4, color=TOE, alpha=0.22, lw=0, zorder=2,
               label="Transect toes (kept)")
    ax.plot(dom.gis, dom.current_toe_m_navd88, color=CURRENT, lw=1.0, ls=(0, (4, 2)), zorder=3,
            label="Current (10 m cells, 1996–97)")
    ax.plot(dom.gis, dom.toe_m_navd88, color=TOE, lw=1.4, zorder=4, label="Transects, 1996 (domain median)")
    fil = dom[dom.source == "filled"]
    if len(fil):
        ax.plot(fil.gis, fil.toe_m_navd88, "o", ms=4, mfc="white", mec=TOE, zorder=5,
                label="Filled from neighbours")
    ax.set_xlim(0.5, 90.5)
    q = kk.toe_m_navd88.quantile([0.005, 0.995])
    ax.set_ylim(np.floor(q.iloc[0] - 0.3), np.ceil(q.iloc[1] + 0.8))
    ax.set_xlabel(getattr(fs, "DOMAIN_AXIS_LABEL", "GIS domain (south → north)"))
    ax.set_ylabel("Dune toe (m NAVD88)")
    fs.open_frame(ax)
    _legend(ax, ncol=2)


PANELS = (  # (file, draw, caption)
    ("step1_transects", panel_map,
     "Shore-normal transects (100 m spacing; the BRIE offset transects, on the model domains' frame) over the 1996 "
     "NOAA lidar (9 ft grid, hillshade with elevation) at GIS {gis}; five transects cross each 500 m domain (dashed "
     "outline), labelled by transect number; open circles mark the dune toes picked on each transect, thick teal segments the model's dune search windows. The ocean is to the right; the 1996 survey "
     "covers only the beach and foredune strip."),
    ("step2_profiles", panel_profiles,
     "Elevation profiles of the five transects of GIS {gis}, sampled every 1 m (bilinear) from the 1996 lidar, in m "
     "NAVD88 with the survey's {offset:+.2f} m correction to the start DEM; distance measured from the shoreline of "
     "transect {ex}. Gaps are areas the survey did not reach."),
    ("step3_shoreline", panel_shoreline,
     "Shoreline on transect {ex} (GIS {gis}): the most seaward point at or above mean high water (MHW, {mhw:g} m "
     "NAVD88) with the following {run:.0f} m also above it. Distances in later panels are measured landward of it."),
    ("step4_crest", panel_crest,
     "Foredune crest on transect {ex} (GIS {gis}): inside the model's dune search window (shaded; the window the "
     "model's dune crest is taken from) and landward of the shoreline, the most seaward peak with a dip of at least "
     "{drop:g} m behind it, or the window's highest point where no peak has one. Where the window spans two ridges "
     "this selects the seaward one."),
    ("step5_toe", panel_toe,
     "Dune toe on transect {ex} (GIS {gis}): the profile point with the greatest vertical distance below the chord "
     "from the shoreline to the foredune crest, i.e. the break in slope between the beach and the dune face, no "
     "lower than the model's berm ({berm:g} m NAVD88; dotted): the dune stands on the berm in the model, so a lower "
     "toe is moved landward to where the profile first reaches the berm (open circle to filled)."),
    ("step6_screen", panel_screen,
     "The five transects of GIS {gis}, aligned at their shorelines, with their dune toes. A toe is retained where "
     "the profile is continuous from the shoreline to the crest and the crest is at least {below:g} m above the "
     "toe; no toe is below the berm ({berm:g} m NAVD88, dotted). Rejected transects in grey."),
    ("step7_domain_toe", panel_domain,
     "Dune toe of GIS {gis}: the median of its retained transect toes (at least {nmin} required; otherwise the "
     "mean of the nearest measured domains on each side; solid), with the mean of the same toes (dashed) for "
     "comparison. The median is used: one high or low transect moves it less. The model takes the toe in m above "
     "MHW (minus {mhw:g} m)."),
    ("result_alongshore", panel_alongshore,
     "Domain dune toes alongshore from the 1996 transects (median of each domain's retained transects; {nmed} "
     "domains measured, {nfill} filled from neighbours) compared with the toes currently used by the model (10 m "
     "cells, mean of the 1996 and 1997 surveys); points are individual retained transects. Transect minus current: "
     "mean {dmean:+.2f} m, RMSE {drmse:.2f} m."),
)


def figures(t, dom, profiles, tr, domains, src, scale, offset, gis):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    out = M_LINE
    c = _Ctx(t, dom, profiles, tr, domains, src, scale, offset, gis)
    dd = dom.toe_m_navd88 - dom.current_toe_m_navd88
    fmt = dict(gis=gis, ex=c.ex, offset=offset, mhw=MHW, run=SHORE_RUN_M, win=CREST_SEARCH_M, berm=BERM_M,
               drop=FOREDUNE_DROP_M, tmin=TOE_MIN_M_NAVD88, below=TOE_BELOW_CREST_M, nmin=MIN_KEPT,
               nmed=int((dom.source == "median").sum()), nfill=int((dom.source == "filled").sum()),
               dmean=dd.mean(), drmse=float(np.sqrt((dd ** 2).mean())))

    # One file per panel, all the same size
    for name, draw, cap in PANELS:
        size = (4.6, 3.6) if draw is panel_map else fs.figsize("double", height=2.9)   # the map is near square
        fig, ax = plt.subplots(figsize=size, layout="constrained")
        draw(ax, c, fs, fig)
        png = out / f"{name}.png"
        fs.save(fig, png, close=True)
        fs.record_caption(png, cap.format(**fmt))
        print(f"wrote {png}")

    # The same panels as one lettered figure
    fig = plt.figure(figsize=(fs.figsize("double")[0], 10.2), layout="constrained")
    gs = fig.add_gridspec(4, 2)
    letters = "abcdefgh"
    caps = []
    for k, (name, draw, cap) in enumerate(PANELS):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        draw(ax, c, fs, fig)
        ax.text(-0.02, 1.02, letters[k], transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right",
                va="bottom")
        caps.append(f"({letters[k]}) " + cap.format(**fmt))
    png = out / f"dune_toe_method_{YEAR}.png"
    fs.save(fig, png, close=True)
    fs.record_caption(png, f"Dune toe from the {YEAR} lidar along the 100 m transects. " + " ".join(caps))
    print(f"wrote {png}")


def overview():
    """Every transect along the island over the 1996 lidar (resampled to 20 m), plan view, north up."""
    import geopandas as gpd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from rasterio.warp import Resampling, reproject
    from rasterio.transform import from_origin
    from matplotlib.colors import LightSource
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    src, scale, offset = survey()
    tr, domains = transects()
    allt = gpd.read_file(TRANSECTS)
    x0, y0, x1, y1 = domains.total_bounds
    res = 20.0
    w, h = int((x1 - x0) / res) + 1, int((y1 - y0) / res) + 1
    z = np.full((h, w), np.nan)
    reproject(S.clean(src.read(1), src.nodata), z, src_transform=src.transform, src_crs=src.crs,
              dst_transform=from_origin(x0, y1, res, res), dst_crs=tr.crs, src_nodata=np.nan, dst_nodata=np.nan,
              resampling=Resampling.average)
    z = z * scale + offset
    km = lambda v, o: (v - o) / 1000.0
    ext = (0, km(x0 + w * res, x0), 0, km(y1, y0))
    fig, ax = plt.subplots(figsize=(4.2, 9.0), layout="constrained")
    hs = LightSource(azdeg=315, altdeg=45).hillshade(np.nan_to_num(z, nan=0.0), vert_exag=4, dx=res, dy=res)
    ax.imshow(np.where(np.isfinite(z), hs, np.nan), cmap="Greys_r", extent=(0, km(x0 + w * res, x0),
              km(y1 - h * res, y0), km(y1, y0)), zorder=1)
    im = ax.imshow(z, cmap="terrain", vmin=-1, vmax=8, alpha=0.8, zorder=2,
                   extent=(0, km(x0 + w * res, x0), km(y1 - h * res, y0), km(y1, y0)))
    used = set(tr.ObjectID)
    poly = domains.set_index("domain_id").geometry
    for r in tr.itertuples():                                    # each transect across its own domain's 2 km box
        seg = r.geometry.intersection(poly[r.gis])
        for part in getattr(seg, "geoms", [seg]):
            c = np.asarray(part.coords)
            if len(c):
                ax.plot(km(c[:, 0], x0), km(c[:, 1], y0), color=TOE, lw=0.4, zorder=3)
    for g in range(10, 91, 10):                                  # every tenth domain named beside it
        b = poly[g].bounds
        ax.text(km(b[0], x0) - 0.1, km(b[1], y0) + 0.25, f"GIS {g}", fontsize=6.5, color=fs.INK, ha="right",
                va="center", zorder=4)
    ax.set_xlim(0, km(x1, x0))
    ax.set_ylim(0, km(y1, y0))
    ax.set_aspect("equal")
    ax.set_xlabel("Easting (km)")
    ax.set_ylabel("Northing (km)")
    cb = fig.colorbar(im, ax=ax, shrink=0.4, pad=0.02, aspect=20)
    cb.set_label("1996 elevation (m NAVD88)", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    fs.open_frame(ax)
    png = M_LINE / "transects_overview.png"
    fs.save(fig, png, close=True)
    fs.record_caption(png, (
        f"The {len(used)} 100 m transects (the BRIE offset step's, 5 per 500 m GIS domain) along Hatteras Island, "
        "each drawn across its domain's 2 km box, over the 1996 NOAA lidar (resampled to 20 m; the survey covers "
        "the beach and foredune strip). Every tenth domain is named at its south edge. Axes in km from the "
        "domains' south-west corner (UTM 18N)."))
    print(f"wrote {png}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gis", type=int, default=EXAMPLE_GIS, help="the example domain for the step figures")
    ap.add_argument("--overview", action="store_true", help="only the island-wide transect map")
    a = ap.parse_args()
    if a.overview:
        return overview()
    figures(*run(), a.gis)


if __name__ == "__main__":
    main()
