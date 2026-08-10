# =============================================================================
# HAT_storm_catalog.py
# Named-storm record for Hatteras Island, NC — SINGLE SOURCE OF TRUTH
# -----------------------------------------------------------------------------
# WHY THIS FILE EXISTS:
#   This catalog previously lived in three places: HAT_storm_record_1984_2021.py
#   (42 storms, with wind/pressure/landfall/note) and both validator scripts
#   (39 storms, name/start/end/cat only). They had drifted — BERTHA 1996,
#   FRAN 1996 and ISAIAS 2020 appeared on the record chart but were never
#   tested against any storm file. Those three are exactly the source="local"
#   entries added by hand from the hurricane-history doc, i.e. the ones most
#   in need of testing.
#
#   Import from here everywhere. Do not copy the list.
#
#       from HAT_storm_catalog import HISTORICAL_STORMS, storms_in_period
#
#   HAT_storm_record_1984_2021.py should also be edited to import from here
#   and delete its own copy.
#
# -----------------------------------------------------------------------------
# FIELDS:
#   name      storm name ("Unnamed" for the 2017 ET event)
#   start/end HURDAT2 full active period — NOT the Hatteras influence window.
#             A storm active two weeks may affect Hatteras for 24–48 h.
#   cat       max Saffir-Simpson category over full lifetime (not at landfall)
#   wind      max sustained wind (mph) corresponding to `cat`
#   pressure  min central pressure (mb) corresponding to `cat`
#   landfall  True  -> documented direct Hatteras/Dare Co. impact or evacuation
#             False -> passed by, minor/no local impact
#             None  -> no local-impact narrative available (NOAA-only)
#   note      local-impact detail, where documented
#   source    "noaa"  -> NOAA Historical Hurricane Tracks, 60 nm buffer of
#                        Hatteras, Dare County NC
#             "local" -> hurricane-history doc only; wind/pressure cross-checked
#                        against NHC/NWS tropical cyclone reports. Track fell
#                        outside the 60 nm buffer.
#
# KNOWN GAPS (deliberate, pending a decision):
#   Earl 2010, Sandy 2012, Maria 2017, Florence 2018, Michael 2018 appear in the
#   hurricane-history doc but are absent from the 60 nm buffer search, and their
#   NC landfalls/tracks were far enough from Hatteras that local impact was
#   minor/indirect per the doc's own text.
#
#   Nor'easters and non-tropical winter storms are NOT included — no catalog
#   available. At Hatteras these are the dominant morphological forcing, so
#   unmatched events in a TWL-derived storm file are EXPECTED, not errors.
# =============================================================================

import pandas as pd

PERIOD_BOUNDARY = 2004  # Period 1 (calibration) / Period 2 (validation)


HISTORICAL_STORMS = [
    # -------------------------------------------------------------------------
    # 1984-2004  (Period 1: calibration)
    # -------------------------------------------------------------------------
    {"name": "Diana",     "start": "1984-09-08", "end": "1984-09-16", "cat": "H4", "wind": 115, "pressure": 949,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Gloria",    "start": "1985-09-16", "end": "1985-10-02", "cat": "H4", "wind": 125, "pressure": 920,  "landfall": True,  "note": "Direct hit \u2014 Cat 2 at landfall, 6\u20138 ft surge", "source": "noaa"},
    {"name": "Kate",      "start": "1985-11-15", "end": "1985-11-23", "cat": "H3", "wind": 105, "pressure": 954,  "landfall": False, "note": "",                                                  "source": "noaa"},
    {"name": "Charley",   "start": "1986-08-13", "end": "1986-08-30", "cat": "H1", "wind": 70,  "pressure": 980,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Alberto",   "start": "1988-08-05", "end": "1988-08-08", "cat": "TS", "wind": 35,  "pressure": 1002, "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Bob",       "start": "1991-08-16", "end": "1991-08-29", "cat": "H3", "wind": 100, "pressure": 950,  "landfall": False, "note": "",                                                  "source": "noaa"},
    {"name": "Danielle",  "start": "1992-09-22", "end": "1992-09-26", "cat": "TS", "wind": 55,  "pressure": 1001, "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Emily",     "start": "1993-08-22", "end": "1993-09-06", "cat": "H3", "wind": 100, "pressure": 960,  "landfall": False, "note": "Dare Co. evacuated, $12M damage",                    "source": "noaa"},
    {"name": "Allison",   "start": "1995-06-03", "end": "1995-06-11", "cat": "H1", "wind": 65,  "pressure": 982,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Arthur",    "start": "1996-06-17", "end": "1996-06-23", "cat": "ET", "wind": 45,  "pressure": 992,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Bertha",    "start": "1996-07-05", "end": "1996-07-17", "cat": "H3", "wind": 115, "pressure": 960,  "landfall": True,  "note": "Landfall nr. Topsail; Cat 2, 104 mph, 5 ft surge",   "source": "local"},
    {"name": "Fran",      "start": "1996-08-23", "end": "1996-09-10", "cat": "H3", "wind": 120, "pressure": 946,  "landfall": False, "note": "Landfall Cape Fear; no Dare Co. evacuation",         "source": "local"},
    {"name": "Josephine", "start": "1996-10-04", "end": "1996-10-16", "cat": "TS", "wind": 60,  "pressure": 970,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Bonnie",    "start": "1998-08-19", "end": "1998-08-31", "cat": "H3", "wind": 100, "pressure": 954,  "landfall": True,  "note": "Dare Co. evacuated, 6\u20138 ft surge",              "source": "noaa"},
    {"name": "Dennis",    "start": "1999-08-24", "end": "1999-09-08", "cat": "H2", "wind": 90,  "pressure": 962,  "landfall": True,  "note": "Dare Co. evacuated, $10M damage",                    "source": "noaa"},
    {"name": "Irene",     "start": "1999-10-12", "end": "1999-10-19", "cat": "H2", "wind": 95,  "pressure": 958,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Arthur",    "start": "2002-07-14", "end": "2002-07-19", "cat": "TS", "wind": 50,  "pressure": 992,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Gustav",    "start": "2002-09-08", "end": "2002-09-15", "cat": "H2", "wind": 85,  "pressure": 960,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Kyle",      "start": "2002-09-20", "end": "2002-10-12", "cat": "H1", "wind": 75,  "pressure": 980,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Isabel",    "start": "2003-09-06", "end": "2003-09-20", "cat": "H5", "wind": 145, "pressure": 915,  "landfall": True,  "note": "Breached island; $167M damage",                      "source": "noaa"},
    # -------------------------------------------------------------------------
    # 2004-2024  (Period 2: validation; data currently through 2021)
    # -------------------------------------------------------------------------
    {"name": "Alex",      "start": "2004-07-31", "end": "2004-08-06", "cat": "H3", "wind": 105, "pressure": 957,  "landfall": True,  "note": "Sound-side flooding, $2.4M damage",                  "source": "noaa"},
    {"name": "Bonnie",    "start": "2004-08-03", "end": "2004-08-14", "cat": "TS", "wind": 55,  "pressure": 1001, "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Ophelia",   "start": "2005-09-06", "end": "2005-09-23", "cat": "H1", "wind": 75,  "pressure": 976,  "landfall": True,  "note": "Hatteras Island evacuated",                          "source": "noaa"},
    {"name": "Alberto",   "start": "2006-06-10", "end": "2006-06-19", "cat": "TS", "wind": 60,  "pressure": 969,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Barry",     "start": "2007-05-31", "end": "2007-06-05", "cat": "TS", "wind": 50,  "pressure": 990,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Gabrielle", "start": "2007-09-08", "end": "2007-09-11", "cat": "TS", "wind": 50,  "pressure": 1004, "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Cristobal", "start": "2008-07-19", "end": "2008-07-23", "cat": "TS", "wind": 55,  "pressure": 998,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Irene",     "start": "2011-08-21", "end": "2011-08-30", "cat": "H3", "wind": 105, "pressure": 942,  "landfall": True,  "note": "Dare Co. evacuated; significant flooding, $54M damage", "source": "noaa"},
    {"name": "Beryl",     "start": "2012-05-25", "end": "2012-06-02", "cat": "TS", "wind": 60,  "pressure": 992,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Arthur",    "start": "2014-06-28", "end": "2014-07-09", "cat": "H2", "wind": 85,  "pressure": 973,  "landfall": True,  "note": "Earliest hurricane on record; Hatteras evacuated",   "source": "noaa"},
    {"name": "Claudette", "start": "2015-07-12", "end": "2015-07-15", "cat": "TS", "wind": 45,  "pressure": 1003, "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Bonnie",    "start": "2016-05-27", "end": "2016-06-09", "cat": "TS", "wind": 40,  "pressure": 1006, "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Colin",     "start": "2016-06-05", "end": "2016-06-08", "cat": "ET", "wind": 50,  "pressure": 987,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Hermine",   "start": "2016-08-28", "end": "2016-09-08", "cat": "H1", "wind": 70,  "pressure": 981,  "landfall": True,  "note": "Villages flooded; 4 ft surge, $5.4M damage",         "source": "noaa"},
    {"name": "Julia",     "start": "2016-09-13", "end": "2016-09-21", "cat": "TS", "wind": 45,  "pressure": 1007, "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Matthew",   "start": "2016-09-28", "end": "2016-10-10", "cat": "H5", "wind": 145, "pressure": 934,  "landfall": False, "note": "Landfall SC; Dare Co. gusts to 94 mph",              "source": "noaa"},
    {"name": "Unnamed",   "start": "2017-08-27", "end": "2017-08-29", "cat": "ET", "wind": 40,  "pressure": 1004, "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Dorian",    "start": "2019-08-24", "end": "2019-09-09", "cat": "H5", "wind": 160, "pressure": 910,  "landfall": True,  "note": "Cat 1 landfall on Hatteras; 4\u20137 ft surge, $14.8M damage", "source": "noaa"},
    {"name": "Arthur",    "start": "2020-05-16", "end": "2020-05-21", "cat": "ET", "wind": 55,  "pressure": 989,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Fay",       "start": "2020-07-05", "end": "2020-07-11", "cat": "TS", "wind": 50,  "pressure": 998,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Isaias",    "start": "2020-07-30", "end": "2020-08-05", "cat": "H1", "wind": 90,  "pressure": 986,  "landfall": True,  "note": "72 mph in Avon; 1\u20133 ft surge, Zone A evacuated", "source": "local"},
    {"name": "Claudette", "start": "2021-06-17", "end": "2021-06-23", "cat": "TS", "wind": 40,  "pressure": 1003, "landfall": None,  "note": "",                                                  "source": "noaa"},
]

# --- convenience -------------------------------------------------------------

def storms_in_period(begin_year, end_year, catalog=None):
    """
    Return catalog entries whose START year falls in [begin_year, end_year],
    with 'start_ts'/'end_ts' Timestamps added.

    NOTE: filtering on start year (not overlap) matches the original validator
    behaviour. A storm starting 2003-12-28 and ending 2004-01-03 belongs to the
    year it started in, so it cannot be double-counted across the two periods.
    """
    cat = catalog if catalog is not None else HISTORICAL_STORMS
    out = []
    for s in cat:
        start = pd.Timestamp(s["start"])
        end   = pd.Timestamp(s["end"])
        if begin_year <= start.year <= end_year:
            out.append({**s, "start_ts": start, "end_ts": end})
    return sorted(out, key=lambda s: s["start_ts"])


if __name__ == "__main__":
    print(f"{len(HISTORICAL_STORMS)} storms in catalog")
    for tag, a, b in [("Period 1", 1984, 2004), ("Period 2", 2004, 2024)]:
        sub = storms_in_period(a, b)
        loc = [s["name"] for s in sub if s["source"] == "local"]
        print(f"  {tag} ({a}-{b}): {len(sub)} storms"
              + (f"  | source='local': {', '.join(loc)}" if loc else ""))
