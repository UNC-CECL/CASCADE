"""
Domain-averaged toes on every collision design (2026-11-23, Laura's suggestion): each module design with one toe per
domain (the mean of its screened 1996 and 1997 toe cells, island median 2.89 m NAVD88) instead of the berm (1.7 m),
at Cs 5, 10, 20, 50, 100 x 1e-6 (HAT_dune_collision.py TOEDOM), against the same design on the berm toe.

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_domain_toes.py

Writes, in the workspace's dune-collision/4-hindcast/4-cs-and-toes/domain-toes/ (D:):
- domain_toes_cs_sensitivity_1996_2009.png (+ supporting/ CSV): the six designs' scores against Cs on the domain
  toes (as HAT_collision_design_cs_sweep), each design's berm-toe run at its current Cs as a ringed point
- <erosion>_<start>_domain_toes_position_1996_2009.png: shoreline position change against CoastSat, the design on
  the berm toe (its current Cs, dashed, the design's colour) and on the domain toes at every Cs (light to dark),
  no collision on the same starting rows grey underneath
- supporting/domain_toes_scores.csv: every run's scores plus the 1998-99 crest drop at GIS 10-12 (lidar: about -2 m)
Runs not made yet are skipped.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-23
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import collision_designs as CD  # noqa: E402
import HAT_collision_design_cs_sweep as W  # noqa: E402  (score, draw, ramp, cs_label, reference)
import HAT_collision_row_mode_sweep as R  # noqa: E402

H = R.H
OUT = CD.OUT / "4-cs-and-toes" / "domain-toes"
DROP_GIS = (10, 11, 12)                    # the lidar's 1998-99 drop: about -1.9, -1.9, -2.0 m


def drop_1999(v):
    """Domain-median crest change over the 1999 storm year (after 1998 -> after 1999) at DROP_GIS, m."""
    import HAT_road_relocation_test as T
    _, c, d = T._load(v)
    out = {}
    for g in DROP_GIS:
        z = T._crest(c.barrier3d[d.gis_to_pad(g)])
        out[f"drop_1999_gis{g}_m"] = float(np.median(z[3]) - np.median(z[2]))
    return out


def tables():
    rows = []
    for d in CD.DESIGNS:
        if d.erosion == "nocollision":
            continue
        k = W.key(d)
        berm_cs = H.VARIANTS[d.run][1]
        if R._done(d.run):
            rows.append(dict(design=d.name.rsplit(", Cs", 1)[0], tag=d.tag, toe="berm", current=True, run=d.run,
                             **W.score(d.run, berm_cs), **drop_1999(d.run)))
        for c in H.TOEDOM_CS:
            v = H.TOEDOM[(k, c)]
            if R._done(v):
                rows.append(dict(design=d.name.rsplit(", Cs", 1)[0], tag=d.tag, toe="domain", current=False, run=v,
                                 **W.score(v, c * 1e-6), **drop_1999(v)))
    return pd.DataFrame(rows)


def positions():
    import HAT_collision_on_off_position as P
    for d in CD.DESIGNS:
        if d.erosion == "nocollision":
            continue
        k = W.key(d)
        ref = W.reference(d)
        name = d.name.rsplit(", Cs", 1)[0]
        runs = [(ref.run, ref.name + " (reference)", "#7F7F7F", 1.4, ref.linestyle),
                (d.run, f"Berm toe, Cs {W.cs_label(round(H.VARIANTS[d.run][1] * 1e6))} (current)", d.colour, 2.6,
                 (0, (4, 2)))]
        for c, colour in zip(H.TOEDOM_CS, W.ramp(d.erosion, len(H.TOEDOM_CS))):
            v = H.TOEDOM[(k, c)]
            if R._done(v):
                runs.append((v, f"Domain toes, Cs {W.cs_label(c)}", colour, 1.5))
        if len(runs) < 3:
            print(f"  skipping {name}: no domain-toe runs yet")
            continue
        out, gis, ch = P.figure(
            1996, 2009, "HAT_dune_collision", tuple(runs),
            f"Shoreline position change, {name}, berm toe vs domain toes, 1996–2009",
            str((OUT / f"{k}_domain_toes_position_1996_2009.png").relative_to(P.OUT)),
            f"the design {name} on the berm toe (1.7 m NAVD88 everywhere; dashed, its current Cs) and on one toe per "
            f"domain, the mean of its screened 1996 and 1997 toe cells (island median 2.89 m), at Cs "
            f"{W.cs_label(H.TOEDOM_CS[0])} to {W.cs_label(H.TOEDOM_CS[-1])} (light to dark), against {ref.name} "
            f"(grey, underneath; runs toedom_{k}_cs<N>_vdefault and {d.run}). " + P.DESIGN_CAPTION,
            post=P._design_zone_shade, legend_rows=5)
        print(f"wrote {out}")


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as S
    S.apply_style()
    (OUT / "supporting").mkdir(parents=True, exist_ok=True)
    t = tables()
    t.round(4).assign(cs=t.cs).to_csv(OUT / "supporting" / "domain_toes_scores.csv", index=False)
    designs = [d for d in CD.DESIGNS if d.erosion != "nocollision"]
    # draw(): one curve per design over Cs (the domain-toe runs), the ringed "current" point its berm-toe run
    tabs = {d.tag: t[t.tag == d.tag] for d in designs}
    W.draw(S, plt, designs, tabs, OUT / "domain_toes_cs_sensitivity_1996_2009.png",
           "Domain-averaged toes on the six collision designs: each design's scores against Cs with one toe per "
           "domain (the mean of its screened 1996 and 1997 toe cells, island median 2.89 m NAVD88; runs "
           "toedom_<design>_cs<N>_vdefault), and ringed, the same design on the berm toe (1.7 m NAVD88) at its "
           "current Cs.", current_apart=True)
    cols = ["design", "toe", "cs", "crest_rmse_m", "crest_rmse_gis9_15_m", "bias_2009_m", "within_0p5_pct",
            "overwash_false_alarm_pct", "shoreline_rmse_m_yr", *[f"drop_1999_gis{g}_m" for g in DROP_GIS]]
    pd.set_option("display.width", 250)
    print(t[[c for c in cols if c in t.columns]].round(3).to_string(index=False))


if __name__ == "__main__":
    if sys.argv[1:] == ["positions"]:
        positions()
    else:
        main()
        positions()
