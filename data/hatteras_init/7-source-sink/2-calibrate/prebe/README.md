# prebe/ — the BE field before each calibration pass

`be_apply_fit_to_config.py` copies `scripts/site_layer/hatteras_site_config.py` here
**before** it overwrites the `HATTERAS_BE_RATES_CALIBRATED` block. Each file is
therefore the field as it stood going *into* that pass, not coming out of it —
the output of pass *n* is the input of the file stamped *n+1*.

## The 2026-09-14 lineage

| File | Holds |
|---|---|
| `..._175853.py` | the 2026-08-24 field, retired (46 nonzero P1, 65 P2) |
| `..._180700.py` | **pass 0, the one-shot solve** (43 nonzero P1, 63 P2) |
| `..._181336.py` | pass 1 |
| `..._182007.py` | pass 2 |

`..._180700.py` is not an archive. `plot_be_zones.py` reads it as the pass-0
field, and the `be_pass0_*` / `iteration_added_*` columns of the exported CSV
are the difference between it and the live config. **Do not delete these four.**

The reason that warning is here: the equivalent pass-0 backup of the retired
lineage, `..._20260824_223143.py`, was never committed and is gone. That is why
the BE zones figure could not be drawn and why those CSV columns came out empty
until the iteration was re-run on 09-14.

## Why they are here and not in `scripts/`

They are data — a snapshot of a solved field — that happens to carry a `.py`
extension, which is also why `hat_layout_check.py`'s Rule 1 never flagged them
(it matches on data *suffixes*). Sitting loose at the root of `scripts/` they
read as stray copies of the site config, which is how the 08-24 one came to be
discarded. Moved here 2026-09-18, beside the rest of the calibrate step's
output; the apply script now writes new ones here too.

Backups are always written to this production folder, even when a what-if pass
sets `HAT_BE_OUTPUT_DIR`. A what-if pass still overwrites the real config, so
its backup is a real backup of the real lineage.
