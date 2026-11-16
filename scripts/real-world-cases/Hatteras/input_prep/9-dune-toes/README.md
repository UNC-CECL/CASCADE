# 9-dune-toes — one dune toe per domain from the 100 m transects

Prototype (2026-11-05), for review with Laura: a simpler toe than `1-barrier3d-domains/3-dune-toes/HAT_dune_toes.py`
(10 m model cells, 1996 and 1997 surveys averaged), from the survey the model starts from alone, set against the
published toe methods. Not used by the model yet; the calibration variants `cal_transect_*`
(`dune-collision/4-hindcast/1-calibration-1996-2009/`) read `1-toes/domain_toes_1996.csv`.

Outputs go to the workspace's `data/hatteras_init/9-dune-toes/` (D:), whose README maps every folder.

## One script per method

Every method runs on the same profiles (the 1996 lidar every 1 m along the 100 m transects) and the same shoreline
(mean high water, 0.36 m NAVD88), so they differ only in their own rules. Each script's docstring gives the
reference, each step and any choice the paper leaves open, and each draws a figure explaining its method on the same
two transects: 119 (GIS 22, a single foredune) and 326 (GIS 63, a dune field).

| Script | Method (reference) | Writes (D:, `1-transects-1996/`) |
|---|---|---|
| `HAT_toe_method_line.py [--gis 22] [--overview]` | Ours: the front dune inside the model's dune search window (or your crest), the point furthest below the shoreline-to-crest line, no lower than the berm (1.7 m), the median of 5 transects | `1-toes/` (the toes), `3-methods/1-line/` (`step1`-`step7`, the method figure) |
| `HAT_toe_method_automorph.py` | Automorph (Itzkin, Moore, Ruggiero & Hacker 2020), ported from his MATLAB (`automorph_matlab/`, a CECL tool, added with his permission): the code as sent, with what the paper adds (35° building check, 3 m crest floor), and adapted for Hatteras (your reviewed crests, only real dunes move the toe, 12 m cap and toe screen) | `3-methods/2-automorph/` |
| `HAT_toe_automorph_review.py [--gis 22] [--sample 50] [--all] [--report]` | Check Automorph's automated crest (3 m floor) and the adapted toe by hand, one window per domain | `2-picks/automorph_review_1996.json` |
| `HAT_toe_automorph_final.py` | The hand-checked Automorph toes, one per transect and one per domain (read by the `cal_automorph_*` hindcast variants) | `1-toes/automorph_transect_toes_1996.csv`, `1-toes/automorph_domain_toes_1996.csv` |
| `HAT_toe_automorph_cells.py [--all]` | Automorph's automatic toe on every 10 m model cell, unreviewed (GIS 9-15, or `--all`), a test of toes per cell against per domain (2026-11-09) | `1-toes/automorph_cell_toes[_all]_1996.csv` |

**2026-11-09 to 11.** What these toes did in the model, and why the hindcast moved to the toe at the berm instead
(2026-11-10): `dune-collision/4-hindcast/TESTING_2026-11-09_to_11.md`. The 2009 lidar's datum was checked the same day
(section 9 there): consistent with the model's DEM and the earlier surveys, no correction.
| `HAT_toe_method_curvature.py [--check] [--why]` | Maximum curvature (Stockdon et al. 2007) | `3-methods/3-curvature/` |
| `HAT_toe_method_pybeach.py [--figure]` | pybeach (Beuzen 2019): perpendicular distance, maximum curvature, relative relief (Wernette et al. 2016) and machine learning (Beuzen et al. 2019; three pre-trained classifiers). Its own environment: see the docstring | `3-methods/4-pybeach/` |
| `HAT_toe_methods_compare.py` | All of them side by side: per transect, per domain, against ours and the current toes | `3-methods/5-comparison/` |

Order: `line` (also writes `1-toes/`), `HAT_toe_pickers.py --viewer-data` (pybeach reads it), then the others in any
order, `compare` last.

## Shared code, checks and pickers

- `HAT_toe_common.py`: the settings, the lidar and the transects, the model's dune windows on a transect, the figure
  palette, and the files of your picks. Not run on its own.
- `HAT_toe_checks.py --crest-choice | --dip`: which crest (front dune in the window, the window's highest point, or the
  front dune within 120 m) and the foredune dip test (0.5 / 1 / 1.5 m). Writes `4-checks/`.
- `HAT_toe_pickers.py`: the pickers, in a window per domain like the dune search window picker.
  - Crest picker, `--crests [--gis 22] [--flagged]`: a left click sets the crest used to find the toe (snaps to the
    highest point within 3 m), and both toes, the line method's and maximum curvature's, are recalculated on the
    spot, with the domain medians in the header. Right click rejects a transect, middle click resets it; enter
    accepts the automatic crest on the rest and goes north. Saved at once to `2-picks/crest_picks_1996.json`; the
    window for each domain is saved as you go to `2-picks/crest_picker_views/`. These crests are for the toe only:
    Barrier3D's dune crests and search windows are not changed.
  - Toe picker, `--qc [--gis 22] [--flagged]`: the same for the toe itself. Left click puts the toe there, right
    click rejects or restores a transect, middle click resets one; enter accepts the rest and goes north, n/p or the
    arrows move without marking, r resets the domain, q saves and quits. Saved to `2-picks/toe_qc_1996.json`.
  - Quick review, `--flagged`: only the transects that can change their domain's toe: leaving it out moves the median
    > 0.25 m; it has no toe where the domain has fewer than 4; in a domain whose median moves > 0.25 m with a 0.5/1.5
    m dip or differs > 0.5 m between the window and front-dune crests, its own toe does; or it has a second bend
    (FLAG_* in the script; 116 transects in 47 domains for 1996).
  - `--viewer-data`: rewrites `1-toes/viewer_data_1996.json`, the profiles the pickers draw.
  Every method's run reads the picks.

## Other scripts

| Script | Writes (in `1-transects-1996/`) |
|---|---|
| `HAT_crest_picks_paper.py [--gis N] [--gis 66 --wide 250]` | `5-crest-pick-figures/`: the final crest per transect, one figure per domain |
| `HAT_domain_imagery_duneline.py [--gis 85]` | `6-imagery/`: a domain on the 1996 aerial photo with the traced dune line |
| `HAT_cs_refit_transect_toes.py` | `7-cs-refit/`: offline Larson Cs refit with the transect toes |
