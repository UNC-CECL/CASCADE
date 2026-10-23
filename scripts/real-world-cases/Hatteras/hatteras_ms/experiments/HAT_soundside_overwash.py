"""
How much of the observed overwash could not have come from the ocean side, and so is likely sound-side?

    python scripts/hatteras_ms/experiments/HAT_soundside_overwash.py

The overwash record has no direction. For every observed new overwash
(image, domain) this compares the highest ocean storm run-up in the image's
window with the domain's lowest surveyed dune crest: run-up below the crest
makes ocean overwash implausible. NOAA sound-side gauge peaks for the named
storms show which windows had high sound water. Scores the adopted run on the
two groups. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-13
"""
from __future__ import annotations

import io
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_ceiling_range as R  # noqa: E402
import HAT_dune_ceiling_uniform as U  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / "storms-overwash" / \
    "2026-10-13-soundside-overwash"
RUN = U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / "dunes" / \
    "2026-10-13-dune-ceiling-range" / "runs" / "c550"          # the adopted single 5.5 m ceiling
MHW_NAVD88_M = 0.36
SOUND = {"8652587": "Oregon Inlet Marina", "8654467": "USCG Hatteras", "8653215": "Rodanthe (sound)"}
EVENTS = {"Bonnie 1998": "1998-08-27", "Dennis 1999": "1999-09-04", "Floyd 1999": "1999-09-16",
          "Isabel 2003": "2003-09-18", "Ophelia 2005": "2005-09-14", "Nov 2006 nor'easter": "2006-11-22",
          "Nor'Ida 2009": "2009-11-13", "Earl 2010": "2010-09-03", "Irene 2011": "2011-08-27",
          "Sandy 2012": "2012-10-29", "Arthur 2014": "2014-07-04", "Joaquin 2015": "2015-10-03",
          "Matthew 2016": "2016-10-09", "Mar 2018 nor'easter": "2018-03-04", "Florence 2018": "2018-09-14",
          "Dorian 2019": "2019-09-06", "Isaias 2020": "2020-08-04", "Ian 2022": "2022-09-30",
          "Lee 2023": "2023-09-16", "Erin 2025": "2025-08-20"}
SOUND_HIGH_M = 1.0           # sound gauge peak, m above MHW, counted as high sound water
# -----------------------------------------------------------------------------


# Peak water level (m above the gauge's MHW) within 3 days of a date, from NOAA; cached
def gauge_peak(station, day):
    cache = OUT / "gauges" / f"{station}_{day}.csv"
    if not cache.exists():
        cache.parent.mkdir(parents=True, exist_ok=True)
        d = pd.Timestamp(day)
        b, e = (d - pd.Timedelta(days=3)).strftime("%Y%m%d"), (d + pd.Timedelta(days=3)).strftime("%Y%m%d")
        raw = ""
        for product in ("hourly_height", "water_level"):
            url = ("https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?"
                   f"product={product}&application=CASCADE_Hatteras&begin_date={b}&end_date={e}"
                   f"&datum=MHW&station={station}&time_zone=GMT&units=metric&format=csv")
            try:
                raw = urllib.request.urlopen(url, timeout=60).read().decode()
            except Exception:
                raw = ""
            if raw.lower().startswith("date time"):
                break
        cache.write_text(raw if raw.lower().startswith("date time") else "date time,v\n", encoding="utf-8")
    df = pd.read_csv(cache)
    v = pd.to_numeric(df.iloc[:, 1], errors="coerce").dropna() if len(df.columns) > 1 else pd.Series(dtype=float)
    return float(v.max()) if len(v) else np.nan


def main():
    import overwash_vs_model as ovm
    from site_layer import hat_env_forcings as env
    from site_layer import hat_overwash as ow
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    ev = pd.DataFrame([dict(event=k, date=pd.Timestamp(v), **{s: gauge_peak(s, v) for s in SOUND})
                       for k, v in EVENTS.items()])
    ev["sound_peak_m_mhw"] = ev[list(SOUND)].max(axis=1)
    rows = []
    for start, (end, storm_years) in U.PERIODS.items():
        rd = sorted((RUN / f"{start}_{end}" / U.PRESET).glob("HAT_*"))[-1]
        c = U.load(rd)
        lowest = {}
        for gis in range(1, 91):                                   # lowest surveyed crest, m NAVD88
            b = c.barrier3d[DOM.gis_to_pad(gis)]
            crest = np.asarray(b.DuneDomain[0], float).max(axis=1)
            lowest[gis] = (crest.min() + b._BermEl) * 10 + MHW_NAVD88_M
        s = pd.read_csv(env.storm_summary_file(*storm_years), parse_dates=["StartTime", "EndTime"])
        s["dated"] = s.EndTime - ovm.GRACE
        s["rhigh_m"] = s.Rhigh * 10 + MHW_NAVD88_M
        d = R.overwash_cells(c, obs, DOM, storm_years)
        imgs = ovm.images_in((storm_years[0], storm_years[1] + 1), obs).set_index("Obs_ID")
        for _, r in d.iterrows():
            im = imgs.loc[r.obs_id]
            win = s[(s.dated > im["from"]) & (s.dated <= im.date)]
            evw = ev[(ev.date > im["from"]) & (ev.date <= im.date + ovm.GRACE)]
            rows.append(dict(period=f"{start}_{end}", image=im.date.date(), gis=r.gis, zone=r.zone,
                             observed=int(r.observed), model=int(r.model),
                             ocean_runup_max_m=float(win.rhigh_m.max()) if len(win) else np.nan,
                             lowest_crest_m=lowest[r.gis],
                             sound_peak_m_mhw=float(evw.sound_peak_m_mhw.max()) if len(evw) else np.nan,
                             sound_events=", ".join(evw.event)))
    t = pd.DataFrame(rows)
    t["ocean_reaches_dune"] = t.ocean_runup_max_m >= t.lowest_crest_m
    t["sound_high"] = t.sound_peak_m_mhw >= SOUND_HIGH_M
    OUT.mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT / "observed_overwash_by_source.csv", index=False)
    ev.to_csv(OUT / "sound_gauge_peaks.csv", index=False)

    print("Sound-side gauge peaks, m above MHW (blank = no data):")
    print(ev[["event", *SOUND, "sound_peak_m_mhw"]].rename(columns=SOUND).round(2).to_string(index=False))
    o = t[t.observed == 1]
    for per, g in o.groupby("period"):
        n = len(g)
        imp = g[~g.ocean_reaches_dune]
        print(f"\n{per}: {n} observed new-overwash domain-images")
        print(f"  ocean run-up below the lowest surveyed crest: {len(imp)} ({len(imp) / n:.0%}); "
              f"of those, in a window with high sound water: {int(imp.sound_high.sum())}")
        for lab, gg in (("ocean could reach the dune", g[g.ocean_reaches_dune]), ("ocean could not", imp)):
            if len(gg):
                print(f"  model hit rate where {lab}: {gg.model.mean():.2f} ({int(gg.model.sum())}/{len(gg)})")
        by = g.assign(implausible=~g.ocean_reaches_dune).groupby("image").agg(
            observed=("gis", "size"), ocean_could_not=("implausible", "sum"), model_hits=("model", "sum"),
            ocean_runup_max_m=("ocean_runup_max_m", "first"), sound_peak_m_mhw=("sound_peak_m_mhw", "first"),
            sound_events=("sound_events", "first"))
        with pd.option_context("display.width", 250, "display.max_colwidth", 60):
            print(by.round(2).to_string())


if __name__ == "__main__":
    main()
