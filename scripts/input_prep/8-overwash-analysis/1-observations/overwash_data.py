"""
The overwash observation record, read off Hatteras_Overwash_Data.xlsx, in the shape the figures need.

    from overwash_data import load_observations, load_storms, assign_capture, SECTIONS

One row per image and one column per domain; a storm is matched to the first image
taken after it (7-day grace). Nothing is plotted here. Details: scripts/input_prep/8-overwash-analysis/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-22
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = next(
    _p for _p in HERE.parents
    if (_p / "pyproject.toml").exists())
import sys  # noqa: E402
sys.path.insert(0, str(REPO / "scripts"))
from site_layer import hat_overwash as ow  # noqa: E402

# Folders resolved by site_layer/hat_overwash.py; OUT_DIR is the root the siblings import
OUT_DIR = ow.OVERWASH_ROOT
# The hand-digitised workbook lives with the data, not the code
XLSX = ow.WORKBOOK

# --- CONFIG ------------------------------------------------------------------
PERIODS = {
    "period1": (1984, 2004),
    "period2": (2004, 2024),
    "combined": (1984, 2024),
}

# Matches ANN_TOWN_SPANS in the site config; Tri-Village is Rodanthe, Waves and Salvo
SECTIONS = [
    ("Cape Point",                 1,  6, "inter"),
    ("Buxton",                     7,  8, "village"),
    ("Buxton–Avon",                9, 20, "inter"),
    ("Avon",                      21, 31, "village"),
    ("Avon–Rodanthe\n(Wimble Shoals)", 32, 67, "inter"),
    ("Rodanthe–Waves–Salvo",      68, 83, "village"),
    ("Pea Island",                84, 90, "inter"),
]

# Images flagged poor or partial, by (year, season); Image_Quality is free text, so not parsed
POOR_QUALITY = {
    (1996, "Fall"),     # poor image quality
    (2006, "Summer"),   # medium-to-poor image quality
    (2023, "Fall"),     # coverage stops at ~domain 17
    (2024, "Spring"),   # coverage starts at ~domain 17
}

# Storm id -> (year, season) of the image that shows it, where the date rule is overruled (README)
CAPTURE_OVERRIDE = {
    "NOR-2022-052": (2023, "Fall"),
}

# Storms missing from the reference sheet but kept from the old figure script
EXTRA_STORMS = [
    dict(id="STM-2009-IDA", name="Ida (extratropical)", cat="TS",
         end=pd.Timestamp("2009-11-13"), approx=False),
    dict(id="STM-2024-DEB", name="DEBBY", cat="H1",
         end=pd.Timestamp("2024-08-09"), approx=False),
]

# Days of grace between a storm's listed last day and an image (README)
CAPTURE_GRACE_DAYS = 7
# -----------------------------------------------------------------------------


# The record: (image table, domain ids 1..90, image x domain matrix of 1 / 0 / NaN)
def load_observations():
    df = pd.read_excel(XLSX, sheet_name="Overwash_Matrix", skiprows=3, header=0)
    df["Year"] = pd.to_numeric(df["Year"], errors="coerce")
    df = df.dropna(subset=["Year"]).copy()
    df["Year"] = df["Year"].astype(int)
    df["Imagery_Date"] = pd.to_datetime(df["Imagery_Date"])
    df = df.sort_values("Imagery_Date").reset_index(drop=True)

    dcols = [c for c in df.columns
             if str(c).replace(".0", "").strip().isdigit()]
    dnum = np.array([int(float(c)) for c in dcols])
    order = np.argsort(dnum)
    domains = dnum[order]
    matrix = (df[dcols].apply(pd.to_numeric, errors="coerce")
              .to_numpy(dtype=float)[:, order])

    obs = df[["Obs_ID", "Imagery_Date", "Year", "Season", "Image_Quality",
              "Source", "Linked_Event_ID"]].copy()
    obs["Season"] = obs["Season"].astype(str)
    obs["poor"] = [(y, s) in POOR_QUALITY
                   for y, s in zip(obs["Year"], obs["Season"])]
    return obs, domains, matrix


# Storm category from the sheet: 'H4' stays, 'Class 5' -> 'NE 5', '~Class 3' -> 'NE 3'
def _category(raw: str) -> str:
    s = str(raw).strip()
    if s.startswith("H") and s[1:2].isdigit():
        return s[:2]
    if "Class" in s:
        digits = [ch for ch in s if ch.isdigit()]
        return f"NE {digits[0]}" if digits else "NE"
    if s.upper().startswith("TS"):
        return "TS"
    if s.upper().startswith("ET"):
        return "ET"
    return s


# Every storm in the reference sheet plus EXTRA_STORMS, sorted by end date
def load_storms():
    df = pd.read_excel(XLSX, sheet_name="Storm_Reference", skiprows=2, header=0)
    df = df.dropna(subset=["Event_ID"])
    out = []
    for _, r in df.iterrows():
        raw = str(r["Search_GE_After"]).strip()
        end = pd.to_datetime(raw, errors="coerce")
        if pd.isna(end):
            continue
        # "Nov 2021" parses to the 1st; the day is unknown, say so.
        approx = not any(ch.isdigit() for ch in raw.split(",")[0].split(" ")[-1:]) \
            or len(raw.split()) < 3
        out.append(dict(id=str(r["Event_ID"]), name=str(r["Storm_Name"]),
                        cat=_category(r["Max_Category"]), end=end, approx=approx))
    out += [dict(s) for s in EXTRA_STORMS]
    for s in out:
        s["year"] = int(s["end"].year)
        s["month"] = s["end"].strftime("%b")
    return sorted(out, key=lambda s: s["end"])


# For each storm, the image that shows it: the first on or after its end less the grace
def assign_capture(storms, obs):
    dates = obs["Imagery_Date"].to_numpy()
    grace = np.timedelta64(CAPTURE_GRACE_DAYS, "D")
    for s in storms:
        if s["id"] in CAPTURE_OVERRIDE:
            y, sea = CAPTURE_OVERRIDE[s["id"]]
            hit = obs.index[(obs["Year"] == y) & (obs["Season"] == sea)]
            s["capture"] = int(hit[0]) if len(hit) else None
            continue
        later = np.nonzero(dates >= np.datetime64(s["end"]) - grace)[0]
        s["capture"] = int(later[0]) if later.size else None

    by_id = {s["id"]: s for s in storms}
    for i, r in obs.iterrows():
        link = r["Linked_Event_ID"]
        if isinstance(link, str) and link in by_id:
            got = by_id[link]["capture"]
            if got != i:
                shown = ("no image" if got is None else
                         obs.loc[got, "Imagery_Date"].strftime("%Y-%m-%d"))
                print(f"  note: {r['Obs_ID']} ({r['Imagery_Date'].date()}) is "
                      f"linked to {link} in the sheet; by date that storm is "
                      f"first shown by {shown}")
    return storms


# Long table: one row per image and domain
def observations_long(obs, domains, matrix) -> pd.DataFrame:
    recs = []
    for i, r in obs.iterrows():
        for j, d in enumerate(domains):
            v = matrix[i, j]
            recs.append((r["Obs_ID"], r["Imagery_Date"].date(), r["Year"],
                         r["Season"], int(d),
                         "" if np.isnan(v) else int(v)))
    return pd.DataFrame(recs, columns=["Obs_ID", "Imagery_Date", "Year",
                                       "Season", "domain", "overwash"])


# One row per storm: which image the figures decided shows it
def storms_table(storms, obs) -> pd.DataFrame:
    recs = []
    for s in storms:
        cap = s.get("capture")
        recs.append(dict(
            storm_id=s["id"], storm=s["name"], peak_category=s["cat"],
            last_day=s["end"].date(), last_day_approx=bool(s["approx"]),
            first_image_after=(obs.loc[cap, "Imagery_Date"].date()
                               if cap is not None else ""),
            first_image_obs_id=(obs.loc[cap, "Obs_ID"] if cap is not None else ""),
            days_to_image=((obs.loc[cap, "Imagery_Date"] - s["end"]).days
                           if cap is not None else ""),
            override=s["id"] in CAPTURE_OVERRIDE))
    return pd.DataFrame(recs)


# Captions
CAPTIONS_HEADER = (
    "# Figure captions\n\n"
    "Written by the scripts in `scripts/input_prep/8-overwash-analysis/`. "
    "Each heading names the figure and the folder it is in. The "
    "figures carry no in-image titles or footnotes on purpose; use these "
    "under them.\n")
# Folder order in CAPTIONS.md, whatever order the scripts ran in
_FOLDER_RANK = {label: i for i, label in enumerate(ow.CAPTION_FOLDERS.values())}


# Replace or add the entry for <name> in CAPTIONS.md
def upsert_caption(name: str, folder: str, text: str) -> Path:
    p = ow.CAPTIONS
    old = p.read_text(encoding="utf-8") if p.exists() else CAPTIONS_HEADER
    parts = re.split(r"(?m)^(?=## )", old)
    header, sections = parts[0], parts[1:]
    key = f"## `{name}`"
    sections = [x for x in sections if not x.startswith(key)]
    sections.append(f"{key} ({ow.CAPTION_FOLDERS[folder]})\n\n{text}\n")

    def rank(sec):
        m = re.search(r"\(([^)]+)\)", sec.splitlines()[0])
        return _FOLDER_RANK.get(m.group(1) if m else "", 9)

    sections.sort(key=rank)
    p.write_text(header.rstrip() + "\n\n"
                 + "\n".join(x.rstrip() + "\n" for x in sections),
                 encoding="utf-8")
    return p


# Drop the entry for a figure that no longer exists
def remove_caption(name: str) -> None:
    p = ow.CAPTIONS
    if not p.exists():
        return
    old = p.read_text(encoding="utf-8")
    parts = re.split(r"(?m)^(?=## )", old)
    key = f"## `{name}`"
    keep = [x for x in parts[1:] if not x.startswith(key)]
    p.write_text(parts[0].rstrip() + "\n\n" + "\n".join(x.rstrip() + "\n" for x in keep),
                 encoding="utf-8")
