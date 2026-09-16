# 2026-09-14-recode

**Question.** The whole 1984-2004 relocation arm (12 runs: three presets x
bdm/nobdm x groin/nogroin, relocations on) re-run under the 2026-09-14 code
on the same v1 topography, so the difference against the stored runs is CODE
only -- chiefly the whole-cell rounding of prescribed displacements, which
moved GIS 11 over the drowning line by one cell.

**Made by.** `tools/HAT_rerun_arm.py`; compared with `tools/HAT_compare_rerun.py`.

**Answer lives in.** `hatteras_site_config.py` (the note above `CELL_M`) and
the compare tool's output.

**Runs deletable?** Yes once the stored v1 matrix is itself archived: both
sides of the comparison go together.
