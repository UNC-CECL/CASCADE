"""
The pickers for checking the toe by hand (2026-11-05), in a window per domain like the dune search window picker.

    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_pickers.py --crests [--gis 22] [--flagged]
    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_pickers.py --qc [--gis 22] [--flagged]
    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_pickers.py --viewer-data

Crest picker (--crests): each domain's transect profiles; a left click sets the crest used to find the toe (snaps to
the highest point within 3 m) and both toes, the line method's and maximum curvature's, are recalculated on the
spot. Right click rejects a transect, middle click resets it; enter accepts the automatic crest on the rest and goes
north. Saved at once to 2-picks/crest_picks_1996.json, and the window for each domain to
2-picks/crest_picker_views/. Toe picker (--qc): the same for the toe itself, saved to 2-picks/toe_qc_1996.json.
--flagged: only the transects that can change their domain's toe (flags()). Every method's run reads the picks.
--viewer-data rewrites the profiles the pickers draw (1-toes/viewer_data_1996.json).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-05
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parent))
from HAT_toe_common import (BERM_M, CHORD, CREST, CREST_FILE, FILL, MHW, MIN_KEPT, PICKS, PROFILE, QC_FILE, TOE,
    TOES, TOE_BELOW_CREST_M, WATER, WINDOW, YEAR, crest_picks, qc_edits)  # noqa: E402
from HAT_toe_method_line import (pick, run)  # noqa: E402

VIEW_DATA_M = (-40, 400)          # the pickers' profiles: m from the shoreline, every 1 m (dune fields reach
                                  # 250 m and more: GIS 32-35, 64-67)
VIEW_PAD_M = 30                   # each domain's view ends this far past its windows and crests
MAP_STEP_M = 2.0                  # the plan-view map's lidar grid
MAP_BEYOND_M = 150                # the map shows this much of each neighbouring domain
DUNELINE_YEAR = {1996: 1997, 2009: 2009}[YEAR]   # the traced dune line nearest the survey (1996: 12 Oct 1997)
DUNELINE_COLOR = "#E7298A"
PHOTO_STEP_M = 1.0                # the aerial photo on the map frame (1996: USGS scan of 14 Oct 1996, ~0.5 m)
PHOTO_LABEL = {1996: "14 Oct 1996 aerial photo", 2009: "2008 aerial mosaic (NOAA IOCM)"}[YEAR]
# The quick review (--qc --flagged) shows only the transects that can change their DOMAIN's toe (the median is what
# the model uses, so a transect that cannot move it needs no look):
FLAG_MEDIAN_M = 0.25              # "swing": leaving this transect out moves the domain median this much
FLAG_FEW = 4                      # "no_toe": a transect without a toe, where the domain has fewer kept than this
FLAG_WINDOW_M = 0.5               # "window": the domain's window and front-dune medians differ this much, and
                                  # this transect's two toes do
FLAG_DIP_M = 0.25                 # "dip": the domain median moves this much with a 0.5 or 1.5 m foredune dip, and
                                  # this transect's toe does
FLAG_NAMES = {"swing": "moves the domain median", "no_toe": "no toe, domain short of transects",
              "window": "window and front-dune toes differ", "dip": "changes with the dip rule",
              "bend": "two bends, toe could be the other one"}


def flags(t, dom, profiles):
    """Why a transect is worth a look: a list of FLAG_NAMES keys, empty where it cannot change its domain's toe."""
    d = dom.set_index("gis")
    toe = {r.transect: r.toe_m_navd88 for r in t.itertuples() if r.kept}
    dip_toe = {}                                                 # transect -> {drop: toe or None}
    for r in t.itertuples():
        s, z = profiles[r.transect][:2]
        wa = profiles[r.transect][5]
        dip_toe[r.transect] = {}
        for drop in (0.5, 1.5):
            q = pick(s, z, wa[1:] if wa else None, drop=drop)
            dip_toe[r.transect][drop] = None if q["toe"] is None or q["reason"] else float(z[q["toe"]])
    out = {}
    for g, grp in t.groupby("gis"):
        ids = list(grp.transect)
        v = {i: toe[i] for i in ids if i in toe}
        med = float(np.median(list(v.values()))) if len(v) >= MIN_KEPT else np.nan
        dip_med = {}
        for drop in (0.5, 1.5):
            w = [dip_toe[i][drop] for i in ids if dip_toe[i][drop] is not None]
            dip_med[drop] = float(np.median(w)) if len(w) >= MIN_KEPT else np.nan
        dip_domain = any(not np.isfinite(m) or abs(m - med) > FLAG_DIP_M for m in dip_med.values())
        win_domain = abs(d.toe_window_m_navd88[g] - d.toe_m_navd88[g]) > FLAG_WINDOW_M
        for r in grp.itertuples():
            f = []
            if r.transect not in v:
                if len(v) < FLAG_FEW:
                    f.append("no_toe")
            else:
                rest = [x for i, x in v.items() if i != r.transect]
                if len(rest) < MIN_KEPT or abs(float(np.median(rest)) - med) > FLAG_MEDIAN_M:
                    f.append("swing")
                if win_domain and r.window_kept and abs(r.window_toe_m_navd88 - r.toe_m_navd88) > FLAG_WINDOW_M:
                    f.append("window")
                if dip_domain and any(x is None or abs(x - r.toe_m_navd88) > FLAG_DIP_M
                                      for x in dip_toe[r.transect].values()):
                    f.append("dip")
            if np.isfinite(r.other_bend_toe_m_navd88):
                f.append("bend")
            out[r.transect] = f
    return out


def view_end(rs):
    """Where a domain's profile view ends (m from the shoreline): past its windows and crests, at least 150 m."""
    far = [v for r in rs for v in ((r.get("w") or [None, None])[1], (r.get("fc") or [None])[0]) if v is not None]
    return int(min(VIEW_DATA_M[1], max([150.0] + [v + VIEW_PAD_M for v in far])))


def viewer_data(t, dom, profiles, *_, path=None):
    """Every transect's profile and picks as JSON for the toe QC viewer (an artifact page): distances in m from
    the transect's shoreline, elevations m NAVD88, 2 decimals."""
    import json
    lo, hi = VIEW_DATA_M
    r2 = lambda v: None if v is None or not np.isfinite(v) else round(float(v), 2)
    fl = flags(t, dom, profiles)
    recs = []
    for r in t.itertuples():
        s, z, p, _, pw, wa = profiles[r.transect][:6]
        rec = dict(id=int(r.transect), gis=int(r.gis), kept=bool(r.kept), reason=r.reason,
                   wkept=bool(r.window_kept), wreason=r.window_reason, flags=fl[r.transect])
        if p["shore"] is not None:
            o = s[p["shore"]]
            xs = np.arange(lo, hi + 1)
            zz = np.interp(xs, s - o, z, left=np.nan, right=np.nan)
            rec["z"] = [r2(v) for v in zz]
            pt = lambda q, k: [r2(s[q[k]] - o), r2(z[q[k]])] if q.get(k) is not None else None
            xy = profiles[r.transect][3]
            u = (xy[-1] - xy[0]) / np.hypot(*(xy[-1] - xy[0]))         # ocean -> sound
            rec |= dict(sz=r2(z[p["shore"]]), fc=pt(p, "crest"), ft=pt(p, "toe"), wc=pt(pw, "crest"),
                        wt=pt(pw, "toe"), w=[r2(wa[1] - o), r2(wa[2] - o)] if wa else None,
                        o=[round(float(xy[p["shore"]][0]), 1), round(float(xy[p["shore"]][1]), 1)],
                        u=[round(float(u[0]), 5), round(float(u[1]), 5)])
        recs.append(rec)
    doms = [dict(gis=int(d.gis), current=r2(d.current_toe_m_navd88), auto=r2(d.toe_m_navd88),
                 window=r2(d.toe_window_m_navd88)) for d in dom.itertuples()]
    out = dict(year=YEAR, mhw=MHW, x0=lo, toe_min=BERM_M, below=TOE_BELOW_CREST_M, min_kept=MIN_KEPT,
               transects=recs, domains=doms)
    path = Path(path or TOES / f"viewer_data_{YEAR}.json")
    path.write_text(json.dumps(out, separators=(",", ":")), encoding="utf-8")
    print(f"wrote {path} ({path.stat().st_size / 1e3:.0f} kB)")
    n = pd.Series([k for v in fl.values() for k in v]).value_counts()
    nf = sum(bool(v) for v in fl.values())
    nd = len({r.gis for r in t.itertuples() if fl[r.transect]})
    print(f"flagged {nf}/{len(fl)} transects in {nd} domains: " + ", ".join(f"{FLAG_NAMES[k]} {v}" for k, v in n.items()))
# crest_x is m landward of the shoreline. run() uses a picked or accepted crest in place of the automatic one, and
# leaves a rejected transect out. Status: picked (clicked), accepted (the automatic crest, confirmed), rejected
CREST_HELP = ("left click = the crest here (snaps to the highest point within 3 m)   right click = reject / un-reject"
              "   middle click = reset this transect\n"
              "enter = accept the rest and go north   n / → = north   p / ← = south   r = reset domain   q = save and quit")
CREST_SNAP_M = 3
PICKER_VIEWS = PICKS / "crest_picker_views"   # the picker window per domain as last seen, saved as you go


def toes_from_crest(z, x0, crest_x, smooth_m=5):
    """Both toes on a viewer profile (1 m, from x0 m) for a crest at crest_x m from the shoreline: the line method
    (furthest below the shoreline-to-crest line) and maximum curvature (Stockdon), each no lower than the berm.
    Returns ((x, z) line, (x, z) curvature), None where there is no answer."""
    z = np.asarray([np.nan if v is None else v for v in z], float)
    i0, i1 = -x0, int(round(crest_x)) - x0
    if i1 - i0 < 3 or not np.isfinite(z[i0:i1 + 1]).all():
        return None, None
    seg = z[i0:i1 + 1]
    floor = lambda k: next((j for j in range(k, len(seg)) if seg[j] >= BERM_M), len(seg) - 1)
    kl = int(np.argmax(np.linspace(seg[0], seg[-1], len(seg)) - seg))
    sm = pd.Series(seg).rolling(smooth_m, center=True, min_periods=1).mean().to_numpy()
    d2 = np.gradient(np.gradient(sm))
    d2[:2], d2[-2:] = -np.inf, -np.inf
    kc = int(np.argmax(d2))
    out = []
    for k in (kl, kc):
        k = k if seg[k] >= BERM_M else floor(k)
        out.append((float(k), float(seg[k])))
    return out[0], out[1]


def domain_map(rs, x0, xend):
    """The domain's 1996 lidar in plan view on the profiles' cross-shore axis: m landward of the domain's most
    seaward shoreline across (ocean left, as on the profiles), m north along the domain up. Sampled every
    MAP_STEP_M along each transect's direction from the shoreline line."""
    from HAT_toe_common import S, sample, survey, transects
    rs = [r for r in rs if r.get("o")]
    if not rs:
        return None
    src, scale, offset = survey()
    crs = transects()[0].crs
    u = np.mean([r["u"] for r in rs], axis=0)
    u /= np.hypot(*u)
    nrm = np.array([-u[1], u[0]])                                 # alongshore
    if nrm[1] < 0:
        nrm = -nrm                                               # pointing north
    O = np.array(min((r["o"] for r in rs), key=lambda o: np.dot(o, u)))   # the most seaward shoreline
    al = [float(np.dot(np.array(r["o"]) - O, nrm)) for r in rs]
    far = max(float(np.dot(np.array(r["o"]) - O, u)) for r in rs)   # the most landward shoreline
    a = np.arange(min(al) - 50 - MAP_BEYOND_M, max(al) + 50 + MAP_BEYOND_M, MAP_STEP_M)
    c = np.arange(x0, far + xend + MAP_STEP_M, MAP_STEP_M)
    A, C = np.meshgrid(a, c, indexing="ij")
    xy = O + C.ravel()[:, None] * u + A.ravel()[:, None] * nrm
    z = sample(src, scale, offset, xy, crs).reshape(A.shape)
    # where each transect's shoreline sits on the common cross-shore axis, and how far north it is
    off = {r["id"]: float(np.dot(np.array(r["o"]) - O, u)) for r in rs}
    north = {r["id"]: float(np.dot(np.array(r["o"]) - O, nrm)) for r in rs}
    # the traced dune line: the map's frame, and where it crosses each transect (the transect's own row)
    dl, dlx = [], {}
    try:
        import geopandas as gpd
        from shapely.geometry import LineString, box
        f = S.REPO / "data" / "hatteras_init" / "2-brie-offset" / "dunelines" / f"duneline_{DUNELINE_YEAR}.geojson"
        line = gpd.read_file(f).to_crs(crs).geometry.union_all()
        corners = [O + cc * u + aa * nrm for cc in (c[0], c[-1]) for aa in (a[0], a[-1])]
        xs_, ys_ = zip(*corners)
        part = line.intersection(box(min(xs_), min(ys_), max(xs_), max(ys_)))
        for g in getattr(part, "geoms", [part]):
            if g.is_empty or not hasattr(g, "coords"):
                continue
            q = np.asarray(g.coords)[:, :2] - O
            dl.append((q @ u, q @ nrm))
        for r in rs:
            o_ = np.array(r["o"])
            row = LineString([o_ + c[0] * u - (o_ - O) @ u * u, o_ + (c[-1] - (o_ - O) @ u) * u])
            hit = row.intersection(line)
            pts = [np.asarray(h.coords)[0] for h in getattr(hit, "geoms", [hit]) if not h.is_empty]
            dlx[r["id"]] = [float((pt - O) @ u) for pt in pts]
    except Exception as e:                                       # the picker works without it
        print(f"[warn] no dune line on the map: {e}")
    rgb = None
    try:                                                         # the 1996 photo on the same frame
        from scipy.ndimage import map_coordinates
        import HAT_domain_imagery_duneline as IM
        if YEAR == 2009:
            IM.PHOTO_GLOB = S.GIS_ROOT / "aerial_imagery" / "2008" / "2008_full_aerial.tif"
        pa = np.arange(a[0], a[-1] + PHOTO_STEP_M, PHOTO_STEP_M)
        pc = np.arange(c[0], c[-1] + PHOTO_STEP_M, PHOTO_STEP_M)
        PA, PC = np.meshgrid(pa, pc, indexing="ij")
        pxy = O + PC.ravel()[:, None] * u + PA.ravel()[:, None] * nrm
        bx0, by0 = pxy.min(axis=0) - 5
        bx1, by1 = pxy.max(axis=0) + 5
        img = IM.photo(crs, (bx0, by0, bx1, by1))                # (3, h, w), north up, IM.RES_M pixels
        rows, cols = (by1 - pxy[:, 1]) / IM.RES_M, (pxy[:, 0] - bx0) / IM.RES_M
        rgb = np.stack([map_coordinates(img[k].astype(float), [rows, cols], order=1).reshape(PA.shape)
                        for k in range(3)], axis=-1)
        rgb = None if not rgb.any() else np.clip(rgb / 255.0, 0, 1)
    except Exception as e:                                       # the picker works without it
        print(f"[warn] no photo on the map: {e}")
    return dict(z=z, a=a, c=c, off=off, north=north, edges=(min(al) - 50, max(al) + 50), dl=dl, dlx=dlx, rgb=rgb)


def draw_map(ax, mp, rs, crest_of, toes, picks, xend, photo=False):
    """A map panel at true scale: the lidar (or, photo=True, the 14 Oct 1996 aerial photo), each transect with its
    dune window, crest (automatic brown, yours amber) and line toe (blue), and the traced dune line."""
    from matplotlib.colors import LightSource
    ax.clear()
    if mp is None:
        ax.text(0.5, 0.5, "no map (rerun --viewer-data)", ha="center", transform=ax.transAxes)
        return
    z, a, c = mp["z"], mp["a"], mp["c"]
    ext = (c[0], c[-1], a[0], a[-1])
    hs = LightSource(azdeg=315, altdeg=45).hillshade(np.nan_to_num(z, nan=0.0), vert_exag=3, dx=MAP_STEP_M,
                                                      dy=MAP_STEP_M)
    if photo and mp.get("rgb") is not None:
        ax.imshow(mp["rgb"], extent=ext, origin="lower")
    elif photo:
        ax.text(0.5, 0.5, "no photo here", ha="center", transform=ax.transAxes)
    else:
        ax.imshow(np.where(np.isfinite(z), hs, np.nan), cmap="Greys_r", extent=ext, origin="lower")
        ax.imshow(z, cmap="terrain", vmin=-1, vmax=8, alpha=0.75, extent=ext, origin="lower")
    for r in rs:
        if r["id"] not in mp["off"]:
            continue
        o, y = mp["off"][r["id"]], mp["north"][r["id"]]
        ax.plot([c[0], c[-1]], [y, y], color="white", lw=0.6, alpha=0.8)
        if r.get("w"):
            ax.plot([o + r["w"][0], o + r["w"][1]], [y, y], color=WINDOW, lw=5, alpha=0.6, solid_capstyle="butt")
        ax.plot(o, y, "o", ms=4, color=WATER)
        st = picks.get(str(r["id"]), {}).get("status")
        if r.get("fc"):
            ax.plot(o + r["fc"][0], y, "^", ms=7, color=CREST, mec="white", mew=0.5,
                    alpha=0.4 if st == "picked" else 1)
        cr = crest_of(r)
        if st == "picked" and cr:
            ax.plot(o + cr[0], y, "^", ms=9, color="#D9822B", mec="white", mew=0.8)
        tl = toes(r)[0]
        if tl:
            ax.plot(o + tl[0], y, "o", ms=6, color=TOE, mec="white", mew=0.5)
        ax.text(c[-1] - 3, y + 6, str(r["id"]), fontsize=8, color="0.1", fontweight="bold")
    for cc, aa in mp.get("dl", []):                              # the traced dune line
        ax.plot(cc, aa, color=DUNELINE_COLOR, lw=1.4, ls=(0, (5, 3)), zorder=4)
    for e in mp["edges"]:                                        # the domain's south and north edges
        ax.axhline(e, color="white", lw=1.0, ls=(0, (5, 3)))
    ax.set_xlim(c[-1], c[0])                                     # as a map, north up: the ocean on the right
    ax.set_aspect("equal")                                       # true scale: angles as on the ground
    ax.set_ylim(a[0], a[-1])
    ax.set_xlabel("m landward of the most seaward shoreline", fontsize=8)
    ax.tick_params(labelsize=8)
    if photo:
        ax.set_title(PHOTO_LABEL, loc="left", fontsize=9)
        ax.tick_params(labelleft=False)
    else:
        ax.set_ylabel("m north along the domain", fontsize=8)
        ax.set_title(f"{YEAR} lidar", loc="left", fontsize=9)


def crest_picker(gis=None, flagged=False):
    """Pick the crest Barrier3D's dune should have, transect by transect; both toes (line method and maximum
    curvature) are recomputed from it as you click. One window per domain with its five profiles (flagged: only
    the flagged transects). Every change is saved at once to CREST_FILE; the next run of the script uses it."""
    import json
    from datetime import datetime
    import matplotlib
    try:
        import tkinter  # noqa: F401
        matplotlib.use("TkAgg", force=True)
    except Exception:
        print("[warn] TkAgg unavailable; the picker needs an interactive backend.")
    import matplotlib.pyplot as plt
    plt.switch_backend(matplotlib.get_backend())
    vd = TOES / f"viewer_data_{YEAR}.json"
    if not vd.exists() or '"flags"' not in vd.read_text(encoding="utf-8")[:4000]:
        viewer_data(*run())
    D = json.loads(vd.read_text(encoding="utf-8"))
    by = {}
    for r in D["transects"]:
        by.setdefault(r["gis"], []).append(r)
    picks = crest_picks()
    x0 = D["x0"]
    ORANGE, AMBER = "#E07B39", "#D9822B"

    def save():
        CREST_FILE.write_text(json.dumps(picks, indent=1), encoding="utf-8")

    def crest_of(r):
        e = picks.get(str(r["id"]), {})
        if e.get("status") == "rejected":
            return None
        if e.get("crest_x") is not None:
            return e["crest_x"], e["crest_z"]
        return tuple(r["fc"]) if r.get("fc") else None

    def toes(r):
        c = crest_of(r)
        return toes_from_crest(r["z"], x0, c[0]) if c and r.get("z") else (None, None)

    def medians(rs):
        out = []
        for j in (0, 1):
            v = [toes(r)[j][1] for r in rs if toes(r)[j] is not None]
            out.append(float(np.median(v)) if len(v) >= D["min_kept"] else None)
        return out

    shown = lambda k: [r for r in by.get(k, []) if r.get("flags")] if flagged else by.get(k, [])
    order = [k for k in range(1, 91) if shown(k)]
    start = gis if gis in order else next((k for k in order if any(str(r["id"]) not in picks for r in shown(k))),
                                          order[0] if order else None)
    i = order.index(start) if start in order else len(order)
    while 0 <= i < len(order):
        g = order[i]
        rs, alls = shown(g), by.get(g, [])
        rs = sorted(rs, key=lambda r: -(r.get("o") or [0, 0])[1])     # north at the top, as on the map
        xend = view_end(rs)
        n = len(next((r["z"] for r in rs if r.get("z")), [None] * (VIEW_DATA_M[1] - x0 + 1)))
        xs = np.arange(x0, x0 + n)
        vis = xs <= xend
        zs = [z for r in rs for z, v in zip(r.get("z") or [], vis) if z is not None and v]
        ylo, yhi = np.floor(min([-0.5] + zs)), np.ceil(max([4.0] + zs) + 0.5)
        fig = plt.figure(figsize=(18, 2.1 * len(rs) + 1.6))
        gsp = fig.add_gridspec(len(rs), 3, width_ratios=(1, 1, 2.3))
        ax_map = fig.add_subplot(gsp[:, 0])
        ax_photo = fig.add_subplot(gsp[:, 1], sharex=ax_map, sharey=ax_map)   # zoom one, both follow
        axes = [fig.add_subplot(gsp[k, 2]) for k in range(len(rs))]
        state = {"go": None}
        mp = domain_map(rs, x0, xend)
        sh = lambda r: mp["off"].get(r["id"], 0.0) if mp else 0.0     # its shoreline on the map's axis
        xlim = (mp["c"][0], mp["c"][-1]) if mp else (xs[0], xend)

        def draw():
            draw_map(ax_map, mp, rs, crest_of, toes, picks, xend)
            draw_map(ax_photo, mp, rs, crest_of, toes, picks, xend, photo=True)
            for ax, r in zip(axes, rs):
                ax.clear()
                e = picks.get(str(r["id"]), {})
                s = e.get("status", "unchecked")
                o = sh(r)
                at = lambda q: (q[0] + o, q[1])
                if r.get("z"):
                    z = np.array([np.nan if v is None else v for v in r["z"]])
                    ax.fill_between(xs + o, z, ylo, color=FILL, lw=0)
                    ax.plot(xs + o, z, color=PROFILE, lw=1.2)
                    if r.get("w"):
                        ax.axvspan(r["w"][0] + o, r["w"][1] + o, color=WINDOW, alpha=0.12, lw=0)
                    ax.plot(o, z[-x0], "o", ms=5, color=WATER, zorder=5)            # its shoreline
                    for xd in (mp["dlx"].get(r["id"], []) if mp else []):         # the traced dune line
                        ax.axvline(xd, color=DUNELINE_COLOR, lw=1.2, ls=(0, (5, 3)), zorder=4)
                    if r.get("fc"):
                        ax.plot(*at(r["fc"]), "^", ms=8, color=CREST, alpha=0.35 if s == "picked" else 1)
                    c = crest_of(r)
                    if s == "picked" and c:
                        ax.plot(*at(c), "^", ms=10, color=AMBER, mec="white", mew=1)
                    tl, tc = toes(r)
                    if c and tl:
                        ax.plot([o, c[0] + o], [z[-x0], c[1]], color=CHORD, lw=0.7, ls=(0, (4, 2)))
                        ax.plot(*at(tl), "o", ms=8, color=TOE)
                    if tc:
                        ax.plot(*at(tc), "D", ms=6, color=ORANGE, mec="white", mew=0.6)
                ax.axhline(D["mhw"], color=WATER, lw=0.7)
                ax.axhline(BERM_M, color=CREST, lw=0.6, ls=(0, (2, 2)))
                ax.set_ylim(ylo, yhi)
                ax.set_xlim(xlim[1], xlim[0])                    # the ocean on the right, as on the map
                c, (tl, tc) = crest_of(r), toes(r)
                f = lambda p: "—" if p is None else f"{p[1]:.2f} m"
                ax.set_title(f"Transect {r['id']}   crest {f(c)}   line toe {f(tl)}   Stockdon toe {f(tc)}",
                             loc="left", fontsize=9.5)
                if r.get("flags"):
                    ax.text(0.99, 0.94, "; ".join(FLAG_NAMES.get(k, k) for k in r["flags"]), transform=ax.transAxes,
                            ha="right", va="top", fontsize=8, color="#B03A2E")
                col = {"picked": AMBER, "rejected": "#B03A2E", "accepted": TOE}.get(s, "0.4")
                ax.set_title(s.upper(), loc="right", fontsize=9, color=col, fontweight="bold")
                if s == "rejected":
                    ax.set_facecolor("#FBEBEB")
                ax.set_ylabel("m NAVD88", fontsize=8)
                ax.tick_params(labelsize=8)
            for ax in axes[:-1]:
                ax.tick_params(labelbottom=False)
            axes[-1].set_xlabel("m landward of the most seaward shoreline (as the maps; dot = this transect's shoreline)",
                                fontsize=9)
            ml, mc = medians(alls)
            f = lambda v: "—" if v is None else f"{v:.2f}"
            n_done = sum(str(r["id"]) in picks for r in rs)
            fig.suptitle(f"GIS {g}   domain toe with these crests: line method {f(ml)}  ·  Stockdon {f(mc)} m   "
                         f"({n_done}/{len(rs)} checked; domain {i + 1} of {len(order)})\n{CREST_HELP}\n"
                         "brown triangle = automatic crest, amber = yours    circle = line toe    diamond = Stockdon "
                         "toe    teal band = the model's dune window    dotted = berm    "
                         f"pink dashed = {DUNELINE_YEAR} dune line traced on the aerial photo    maps: white dashed = "
                         "domain edges, true scale",
                         fontsize=9.5, x=0.01, ha="left")
            fig.canvas.draw_idle()
            if state.get("laid_out"):                            # the window as you see it, kept per domain
                PICKER_VIEWS.mkdir(parents=True, exist_ok=True)
                fig.savefig(PICKER_VIEWS / f"crest_picker_gis{g:02d}.png", dpi=110)

        def stamp(r, **kw):
            picks[str(r["id"])] = dict(gis=g, updated=datetime.now().isoformat(timespec="seconds"), **kw)
            save()
            draw()

        def on_click(ev):
            if ev.inaxes not in list(axes) or ev.xdata is None:
                return
            r = rs[list(axes).index(ev.inaxes)]
            k = str(r["id"])
            if ev.button == 1 and r.get("z"):
                z = np.array([np.nan if v is None else v for v in r["z"]])
                x = int(round(ev.xdata - sh(r)))                # back to m from this transect's shoreline
                lo, hi = max(1, x - CREST_SNAP_M) - x0, min(xs[-1], x + CREST_SNAP_M) - x0 + 1
                if hi > lo and np.isfinite(z[lo:hi]).any():
                    j = lo + int(np.nanargmax(z[lo:hi]))
                    stamp(r, status="picked", crest_x=float(xs[j]), crest_z=float(z[j]))
            elif ev.button == 3:
                if picks.get(k, {}).get("status") == "rejected":
                    picks.pop(k)
                    save()
                    draw()
                else:
                    stamp(r, status="rejected", crest_x=None, crest_z=None)
            elif ev.button == 2 and k in picks:
                picks.pop(k)
                save()
                draw()

        def on_key(ev):
            if ev.key in ("enter", "return"):
                for r in rs:
                    if str(r["id"]) not in picks and r.get("fc"):
                        picks[str(r["id"])] = dict(gis=g, status="accepted", crest_x=r["fc"][0], crest_z=r["fc"][1],
                                                   updated=datetime.now().isoformat(timespec="seconds"))
                save()
                state["go"] = +1
                plt.close(fig)
            elif ev.key in ("n", "right"):
                state["go"] = +1
                plt.close(fig)
            elif ev.key in ("p", "left"):
                state["go"] = -1
                plt.close(fig)
            elif ev.key == "r":
                for r in rs:
                    picks.pop(str(r["id"]), None)
                save()
                draw()
            elif ev.key in ("q", "escape"):
                state["go"] = 0
                plt.close(fig)

        draw()
        fig.tight_layout(rect=(0, 0, 1, 1 - 0.95 / fig.get_size_inches()[1]))
        state["laid_out"] = True
        draw()                                                   # and saves the view as it opens
        fig.canvas.mpl_connect("button_press_event", on_click)
        fig.canvas.mpl_connect("key_press_event", on_key)
        plt.show()
        if not state["go"]:
            break
        i += state["go"]
    save()
    n = {}
    for e in picks.values():
        n[e["status"]] = n.get(e["status"], 0) + 1
    print(f"saved {CREST_FILE}: {len(picks)} transects {n}")
QC_HELP = ("left click = put the toe here   right click = reject / un-reject   middle click = reset this transect\n"
           "enter = accept the rest and go north   n / → = north   p / ← = south   r = reset domain   q = save and quit")


def qc(gis=None, flagged=False):
    """The toe picker: one window per domain with its five transect profiles (as the figures), like the dune
    search window picker. Click to move a toe, right click to reject a transect; every change is saved to QC_FILE
    at once, and run() takes it into the qc_ columns and the domain's toe_qc_m_navd88. flagged: the quick review,
    only the domains with flagged transects (flags()) and in each only those transects."""
    import json
    from datetime import datetime
    import matplotlib
    try:
        import tkinter  # noqa: F401
        matplotlib.use("TkAgg", force=True)
    except Exception:
        print("[warn] TkAgg unavailable; the picker needs an interactive backend.")
    import matplotlib.pyplot as plt
    plt.switch_backend(matplotlib.get_backend())
    vd = TOES / f"viewer_data_{YEAR}.json"
    if not vd.exists() or '"flags"' not in vd.read_text(encoding="utf-8")[:4000]:
        viewer_data(*run())
    D = json.loads(vd.read_text(encoding="utf-8"))
    by = {}
    for r in D["transects"]:
        by.setdefault(r["gis"], []).append(r)
    dom = {d["gis"]: d for d in D["domains"]}
    edits = qc_edits()
    x0 = D["x0"]

    def save():
        QC_FILE.write_text(json.dumps(edits, indent=1), encoding="utf-8")

    def z_at(r, x):
        i = int(round(x)) - x0
        return r["z"][i] if r.get("z") and 0 <= i < len(r["z"]) else None

    def toe_of(r):
        e = edits.get(str(r["id"]), {})
        if e.get("status") == "rejected":
            return None
        if e.get("status") == "moved":
            return e["toe_z"]
        return r["ft"][1] if r["kept"] and r.get("ft") else None

    def median(rs, use_edits):
        v = [toe_of(r) if use_edits else (r["ft"][1] if r["kept"] and r.get("ft") else None) for r in rs]
        v = [x for x in v if x is not None]
        return float(np.median(v)) if len(v) >= D["min_kept"] else None

    shown = lambda k: [r for r in by.get(k, []) if r.get("flags")] if flagged else by.get(k, [])
    order = [k for k in range(1, 91) if shown(k)]                # the domains this review visits
    if flagged:
        print(f"quick review: {sum(len(shown(k)) for k in order)} flagged transects in {len(order)} domains")
    start = gis if gis in order else next((k for k in order if any(str(r["id"]) not in edits for r in shown(k))),
                                          order[0] if order else None)
    i = order.index(start) if start in order else len(order)
    while 0 <= i < len(order):
        g = order[i]
        rs, alls = shown(g), by.get(g, [])
        zs = [z for r in rs for z in (r.get("z") or []) if z is not None]
        ylo, yhi = np.floor(min([-0.5] + zs)), np.ceil(max([4.0] + zs) + 0.5)
        fig, axes = plt.subplots(len(rs), 1, figsize=(11, 2.1 * len(rs) + 1.3), sharex=True, squeeze=False)
        axes = axes[:, 0]
        state = {"go": None}
        xs = np.arange(x0, x0 + len(rs[0]["z"]) if rs and rs[0].get("z") else VIEW_DATA_M[1] + 1)
        xend = view_end(rs)

        def draw():
            for ax, r in zip(axes, rs):
                ax.clear()
                e = edits.get(str(r["id"]), {})
                s = e.get("status", "unchecked")
                if r.get("z"):
                    z = np.array([np.nan if v is None else v for v in r["z"]])
                    ax.fill_between(xs, z, ylo, color=FILL, lw=0)
                    ax.plot(xs, z, color=PROFILE, lw=1.2)
                    if r.get("w"):
                        ax.axvspan(*r["w"], color=WINDOW, alpha=0.12, lw=0)
                    if r.get("fc"):
                        ax.plot(*r["fc"], "^", ms=8, color=CREST)
                    if r.get("wt"):
                        ax.plot(*r["wt"], "o", ms=7, mfc="none", mec=WINDOW, mew=1.4)
                    if r.get("ft"):
                        ax.plot(*r["ft"], "o", ms=8, color=TOE, alpha=0.4 if s == "moved" else 1)
                    if s == "moved":
                        ax.plot(e["toe_x"], e["toe_z"], "o", ms=9, color="#D9822B", mec="white", mew=1)
                ax.axhline(D["mhw"], color=WATER, lw=0.7)
                ax.axhline(D["toe_min"], color="#B03A2E", lw=0.6, ls=(0, (2, 2)))
                ax.set_ylim(ylo, yhi)
                ax.set_xlim(xs[0], xend)
                auto = f"{r['ft'][1]:.2f} m" if r.get("ft") else "none"
                extra = "" if r["kept"] else f" ({r['reason']})"
                mine = f" · yours {e['toe_z']:.2f} m" if s == "moved" else ""
                col = {"moved": "#D9822B", "rejected": "#B03A2E", "accepted": TOE}.get(s, "0.4")
                ax.set_title(f"Transect {r['id']}   auto toe {auto}{extra}{mine}", loc="left", fontsize=9.5)
                if r.get("flags"):
                    ax.set_title("; ".join(FLAG_NAMES.get(k, k) for k in r["flags"]), loc="center", fontsize=9,
                                 color="#B03A2E")
                ax.set_title(s.upper(), loc="right", fontsize=9, color=col, fontweight="bold")
                if s == "rejected":
                    ax.set_facecolor("#FBEBEB")
                ax.set_ylabel("m NAVD88", fontsize=8)
                ax.tick_params(labelsize=8)
            axes[-1].set_xlabel("Distance landward of the shoreline (m)")
            a, b = median(alls, False), median(alls, True)
            d = dom[g]
            f = lambda v: "—" if v is None else f"{v:.2f}"
            n_done = sum(str(r["id"]) in edits for r in rs)
            fig.suptitle(f"GIS {g}   domain toe: auto {f(a)}  ·  with your checks {f(b)}  ·  window {f(d['window'])}"
                         f"  ·  current {f(d['current'])} m   ({n_done}/{len(rs)} checked)"
                         + (f"   quick review: {len(rs)} of {len(alls)} transects shown, domain {i + 1} of "
                            f"{len(order)}" if flagged else "") + f"\n{QC_HELP}",
                         fontsize=9.5, x=0.01, ha="left")
            fig.canvas.draw_idle()

        def stamp(r, **kw):
            edits[str(r["id"])] = dict(gis=g, updated=datetime.now().isoformat(timespec="seconds"), **kw)
            save()
            draw()

        def on_click(ev):
            if ev.inaxes not in list(axes) or ev.xdata is None:
                return
            r = rs[list(axes).index(ev.inaxes)]
            k = str(r["id"])
            if ev.button == 1:
                x = int(round(ev.xdata))
                z = z_at(r, x)
                if z is not None:
                    stamp(r, status="moved", toe_x=x, toe_z=z)
            elif ev.button == 3:
                if edits.get(k, {}).get("status") == "rejected":
                    edits.pop(k)
                    save()
                    draw()
                else:
                    stamp(r, status="rejected", toe_x=None, toe_z=None)
            elif ev.button == 2 and k in edits:
                edits.pop(k)
                save()
                draw()

        def on_key(ev):
            if ev.key in ("enter", "return"):
                for r in rs:
                    if str(r["id"]) not in edits and r.get("ft"):
                        edits[str(r["id"])] = dict(gis=g, status="accepted", toe_x=r["ft"][0], toe_z=r["ft"][1],
                                                   updated=datetime.now().isoformat(timespec="seconds"))
                save()
                state["go"] = +1
                plt.close(fig)
            elif ev.key in ("n", "right"):
                state["go"] = +1
                plt.close(fig)
            elif ev.key in ("p", "left"):
                state["go"] = -1
                plt.close(fig)
            elif ev.key == "r":
                for r in rs:
                    edits.pop(str(r["id"]), None)
                save()
                draw()
            elif ev.key in ("q", "escape"):
                state["go"] = 0
                plt.close(fig)

        draw()
        fig.tight_layout(rect=(0, 0, 1, 1 - 0.75 / fig.get_size_inches()[1]))
        fig.canvas.mpl_connect("button_press_event", on_click)
        fig.canvas.mpl_connect("key_press_event", on_key)
        plt.show()                                               # blocks until the window closes
        if not state["go"]:                                      # q, escape or the window's X
            break
        i += state["go"]
    save()
    n = {}
    for e in edits.values():
        n[e["status"]] = n.get(e["status"], 0) + 1
    print(f"saved {QC_FILE}: {len(edits)} transects checked {n}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gis", type=int, default=None, help="the domain to start at")
    ap.add_argument("--crests", action="store_true", help="the crest picker")
    ap.add_argument("--qc", action="store_true", help="the toe picker")
    ap.add_argument("--flagged", action="store_true", help="only the flagged transects (quick review)")
    ap.add_argument("--viewer-data", action="store_true", help="only rewrite the pickers' JSON")
    a = ap.parse_args()
    if a.viewer_data:
        return viewer_data(*run())
    if a.qc:
        return qc(a.gis, flagged=a.flagged)
    return crest_picker(a.gis, flagged=a.flagged)


if __name__ == "__main__":
    main()
