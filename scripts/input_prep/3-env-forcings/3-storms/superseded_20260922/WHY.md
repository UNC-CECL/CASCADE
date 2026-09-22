# storm_check/ — retired 2026-09-22

The two storm validators that `storm_validation/HAT_validate_storms.py`
replaced. Its own header says so:

    REPLACES: HAT_validate_storms.py + HAT_validate_storms_windowchange.py
              (which differed only in period config and raw_offset convention)

**Why it was moved rather than left beside the replacement.** Both folders held
a file called `HAT_validate_storms.py`, and they were the only duplicate
basename in `scripts/`. Two files with one name and different behaviour, in
sibling folders, one superseded by the other — opening the wrong one is a
mistake nothing would catch, and `input_prep/README.md` described the pair as
"both validators", which made it sound like a deliberate division of labour.

What the replacement does that these did not, from its own header: the storm
catalog moved out into `HAT_storm_catalog.py`, because the list lived in three
files and had drifted (BERTHA/FRAN 1996 and ISAIAS 2020 were on the record
chart but tested by nothing); and it auto-detects both the v3 summary schema
and the `HAT_create_storms` readable schema, which is why two scripts were
needed before.

Nothing here is maintained or expected to run. Rule 4 of `ORGANIZATION.md`:
retirement is a dated folder with a note, and deletes nothing.

Note `storm_check_<window>.png` is still a LIVE output filename, written by
`storm_validation/HAT_validate_storms.py`. The name outlived the folder; do
not grep for "storm_check" and assume every hit is retired.
