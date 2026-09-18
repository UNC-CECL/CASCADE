# 7-source-sink — scripts

The background-erosion (BE) calibration: the per-domain source/sink field that
CASCADE carries as `DOMAIN_BE_RATES`. It is derived from what the modules could
NOT explain — the residual between the LOESS-smoothed CoastSat rate and the
model's own LRR — so every value has a named physical zone behind it.

`scripts/site_layer/hatteras_site_config.py` is the source of truth for the field; the
copies under `data/hatteras_init/7-source-sink/` are exported FROM it by stage
4, never maintained alongside it.

## The stages

| stage | file | what it does |
|---|---|---|
| `1-prepare/` | `HAT_backfill_run_lrr.py` | Adds `lrr_m_yr` / `lrr_r2` to run rate CSVs written before the LRR estimator existed, recovering them exactly from each run's own `*_shoreline_matrix.npy`. A precondition, not a fit: stage 2 reads the LRR column and stops if a run lacks it. All 191 current runs already carry it — this is for a restored archive. |
| `2-calibrate/` | `HAT_be_zone_residual_fit.py` | The calibration. LOESS-smooths the observed rate, differences it against the base-run LRR, identifies spatially coherent zones with a named mechanism, and writes `DOMAIN_BE_RATES*.txt` plus the metrics tables. |
| | `HAT_be_apply_fit_to_config.py` | Writes those rates into `hatteras_site_config.py`, preserving the two locked ends and the zone labels that a bare paste would destroy. `--add` accumulates instead of replacing. |
| | `HAT_be_edge_domain_solve.py` | The two locked end domains (GIS 1, 90), which are solved separately by Newton steps on a secant rather than fitted from a residual. Reads finished runs and prints the next probe; it does not run the model. |
| `3-figures/` | `HAT_plot_be_convergence.py` | Did the fixed-point solve converge, and was the zone set fixed before it ran. |
| | `HAT_plot_be_zones.py` | Which domains were eligible at all, and how much of each rate came from the one-shot solve versus the iteration. |
| | `HAT_plot_groin_reserved_residual.py` | Why the largest residual in the hindcast (D6) is deliberately left uncorrected. |
| `4-export/` | `HAT_export_be_calibration.py` | Publishes the converged field to `data/hatteras_init/7-source-sink/` — dicts, per-domain CSV, figures, provenance README. Refuses to run on figures older than the newest calibBE run. |

## Stage 2 is a loop, not three steps in a row

The numbering says stage 2 comes after stage 1 and before stage 3. It does not
say the three files inside it run once each in listed order. The documented
method (`hatteras_site_config.py`, METHOD) interleaves the interior fit with the
edge solve:

```
pass 0     interior from the edgeBE base runs          replace
pass 1-3   interior from the calibBE base runs         --add
GIS 90     re-solved after the interior settled        Newton, +3.0 probe
pass 4     final interior pass at the final edge values --add
```

Each pass is fit → apply → re-run the model → fit again. The additive form is
the point: imposing X m/yr of background erosion at a domain moves that domain's
rate by well under X once BRIE has diffused it alongshore, so each pass closes a
fraction of whatever misfit remains and no estimate of that fraction is ever
needed.

Zone MEMBERSHIP is identified once and frozen (`FROZEN_ZONE_DOMAINS`). Only
magnitude iterates. Re-deriving zones each pass would let the arithmetic rewrite
the science — later passes would start correcting the spillover of earlier ones,
which never terminates.

## Run

```
cd scripts/input_prep/7-source-sink

python 1-prepare/HAT_backfill_run_lrr.py --check       # only after restoring an archive

python 2-calibrate/HAT_be_zone_residual_fit.py                      # pass 0
python 2-calibrate/HAT_be_apply_fit_to_config.py --check
python 2-calibrate/HAT_be_apply_fit_to_config.py
#   re-run the model at calibBE, then for each further pass:
HAT_BE_BASE_PRESET=calibBE python 2-calibrate/HAT_be_zone_residual_fit.py
python 2-calibrate/HAT_be_apply_fit_to_config.py --add
#   and for the ends:
python 2-calibrate/HAT_be_edge_domain_solve.py --period 1996 --run <run_name>

python 3-figures/HAT_plot_be_convergence.py
python 3-figures/HAT_plot_be_zones.py
python 3-figures/HAT_plot_groin_reserved_residual.py

python 4-export/HAT_export_be_calibration.py --check
```

`--check` writes nothing anywhere it is offered. Use it first: stage 2 edits the
site config in place and stage 4 moves the previous export to `superseded_<date>/`.

## HAT_BE_OUTPUT_DIR

An exploratory pass — a different Hs, a trial base run — MUST set
`HAT_BE_OUTPUT_DIR`. The production directory holds a converged calibration
whose stopping point is a recorded scientific claim, and a what-if pass
overwriting it would destroy the provenance silently.

Both `HAT_be_zone_residual_fit.py` and `HAT_be_apply_fit_to_config.py` honour it,
so a redirected pass is applied from the directory it actually wrote. The apply
step prints the path it read before it reads it, so which calibration is being
applied is never left to be inferred from the environment.

## Where things land

Tables go to `data/hatteras_init/7-source-sink/2-calibrate/<pair>/` (mirroring the
code folder that writes them), figures to `.../7-source-sink/3-figures/<pair>/`,
and the stage 4 export to `.../7-source-sink/4-export/`, with its README at the
top of `7-source-sink/`. A pair is `<p1start>_<p1end>__<p2start>_<p2end>`, and
every pair has one, the default `1984_2004__2004_2024` included (2026-09-18;
before that the default wrote to the unlabelled root of both folders). Config
backups go to `2-calibrate/prebe/`, shared by every pair. Every script resolves
these through `scripts/site_layer/hat_source_sink.py`; do not type them.

Style: `scripts/site_layer/hat_figure_style.py`. No in-image titles or footnotes; the words
are in `CAPTIONS.md`.
