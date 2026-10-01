# 6-scr-smooth - LOWESS smoothing of the CoastSat shoreline change rates

Takes the per-transect and per-domain LRR tables that `5-scr` produces and
smooths them along the coast, so the noisy transect signal becomes a curve the
model can be scored against.

    5-scr/3-rates/coastsat/lrr/<window>/transect_lrr_full.csv     906 rows, one per transect (~50 m spacing)
    5-scr/3-rates/coastsat/lrr/<window>/domain_lrr_summary.csv    90 rows, one per GIS domain

(Both sat under `5-scr/CoastSat/<period>/` until 2026-09-18. A rate fit spans
an interval, so the folder is a WINDOW, `<start>_<end>` - rule 2 of
`ORGANIZATION.md`.)

**Nothing in this stage writes a model input.** It produced a *decision* -
smooth transects first, then average to domains, at a 7-10 domain window - and
the hindcast implements that decision itself in
`scripts/cascade_pipeline/coastsat_lowess.py`. No CSV crosses from here into a
run. See "What the hindcast actually runs" below.

## The stage

```
lowess_method_comparison.py   CURRENT. Transect-first vs domain-first smoothing.
lowess_dsas_vs_coastsat.py    Side question: DSAS vs CoastSat as sources.
```

Both scripts run from anywhere - every path is anchored on the `pyproject.toml`
at the repo root, not typed as an absolute literal.

Every script here is `lowess_<what it compares>`. The `HAT_` prefix they carried
until 2026-09-22 was dropped to match `5-scr`, the stage this one reads from:
the two shoreline stages are read together, and a reader should not have to
remember that one spells its scripts differently. That does leave these two as
the only `input_prep` stages without the prefix - `0-elevation`,
`3-env-forcings`, `4-mgmt-forcings` and `7-source-sink` still carry it on 46 of
their 50 active scripts. Each script writes its products to a folder under
`data/hatteras_init/6-scr-smooth/` named for what it compares
(`method_comparison/`, `dsas_vs_coastsat/`; they were `<script name>_output/`
until 2026-09-18), resolved through `site_layer/hat_observed_rates.py` - scripts live under
`scripts/`, products under `data/`, the same split every other stage uses, and
those folders are gitignored because one run regenerates them. The one exception
is the snapshot, which keeps the live module's filename `coastsat_lowess.py` on
purpose - it is a mirror, and the name is what makes it one.

Both reconfigure stdout to UTF-8 at import. Their progress lines use arrows and
en-dashes, and a cp1252 Windows console cannot encode those: the DSAS script
used to die in its closing summary *after* writing every figure, which looks
like a failed run that had in fact finished.

### lowess_method_comparison.py

The one that settled the method. Runs **both** smoothers on every call:

- **transect-first** - LOWESS over 906 transects using along-coast metres as
  x, then average the smoothed signal down to 90 domains;
- **domain-first** - LOWESS over the 90 pre-averaged domain means directly.

`LOWESS_WINDOW_KM = 3.5` (7 domains) is the primary window;
`COMPARE_WINDOWS_KM = [2.5, 3.5, 5.0]` (5, 7, 10 domains) is the sensitivity set.
Writes to `data/hatteras_init/6-scr-smooth/method_comparison/`
in four numbered subfolders:
`01_transect_based`, `02_domain_averaged`, `03_cascade_inputs`,
`04_method_comparison`.

The two methods agree across most of the island and separate only where the
gradient is steepest - Rodanthe around domains 76-80 in period 1, Buxton around
7-8 in period 2, and the first few domains. That is the argument for
transect-first: domain-first pre-averages away the structure exactly where it
matters.

Verified against the module the runs use: at a matched 7-domain window this
script and `cascade_pipeline/coastsat_lowess.py` agree to 8.9e-16 m/yr. The
decision recorded here and the code implementing it are the same method, not
two lookalikes.

Its transect-space figures carry the same geographic marks as the domain-space
ones. The `T_*` positions are in metres, derived from the domain-space values
using this script's own conventions - domain *d* occupies `[(d-1)*500, d*500)`
and plots at its band centre `(d-0.5)*500` - so an inclusive span `lo..hi`
becomes `((lo-1)*500, hi*500)` and a point feature becomes `(d-0.5)*500`.
Re-derive them if `DOMAIN_SPACING_M` changes. One mismatch is deliberately left
and commented: `GROINS` puts the Buxton groin at domain 6, band centre 2750 m,
while the hindcast places it at GIS 5.5 - the domain 5/6 boundary, 2500 m.

`SKIP_SOUTHERN_DOMAINS = 10` withholds the smoothed series across domains
1-10 here too - the same guard as the DSAS script and the hindcast, so all
three now cut at the same place. It is applied at the two smoothing functions'
returns, after the fit, which is the whole implementation: every figure and the
export inherit it from there, and `aggregate_to_domains` picks it up for free
because a per-domain mean of an all-NaN group is NaN. The band is drawn in both
x-spaces - domains 1-10 in domain space, 0-5000 m in transect space, which is
the same cut since domain *d* occupies `[(d-1)*500, d*500)`.

Figures that draw raw data anyway simply show it across the zone. The
smoothed-only ones would otherwise be blank there, so they draw the raw domain
means as a dotted line - which is what the hindcast shows in that zone too.

It bites less here than in the domain-space scripts: smoothing runs at transect
resolution, ~10 points per domain, so the edge fit has far more local support.
This script never showed the -6.21 m/yr excursion the DSAS one did. The guard
is applied for consistency with what the model is actually scored against, not
because the artifact appeared.

Note the folder name `03_cascade_inputs` overstates it: nothing reads that CSV,
and it is smoothed at this script's primary window (7 domains), whereas the
hindcast scores against the 10-domain curve. Its `cs_lrr_smooth_*` columns are
blank across domains 1-10 like everything else; the raw columns are complete.
Treat it as a table to read, not as an input to a run.

### lowess_dsas_vs_coastsat.py

A different question - not method, but **source**: does DSAS agree with CoastSat?
It runs on the older **1978-1997 / 1997-2019** period pair (from
`5-scr/archive/coastsat_lrr_superseded_20260810/`), domain-level only, at a fixed
`LOWESS_FRAC = 0.15`. Not part of the 1984-2024 lineage; keep the period
difference in mind before reading its figures next to the others. The legends
name the periods and nothing else - they used to say "Calibration Period" and
"Validation Period", which name the hindcast's split, not this pair - and both
two-panel figures carry a footnote saying which pair they are.

**What it shows.** DSAS reads systematically more accretional than CoastSat in
1978-1997: +1.74 m/yr averaged over all 90 domains, r = 0.64. The offset
largely closes in 1997-2019 (-0.69 m/yr, r = 0.90). The disagreement is
structured, not noise.

Three ways these figures used to mislead, each now marked:

- **DSAS has no 1997-2019 data for domains 1-7** (83 of 90). Its line simply
  started late beside a full-span CoastSat line, which reads as agreement
  rather than absence. `shade_missing_dsas()` hatches and labels the gap on
  every panel that draws DSAS. It finds the runs from the data, so it stays
  correct if coverage changes, and draws nothing when a period is complete. It
  matters most on the difference panels, where an absent bar reads as *zero
  difference* - a meaningful value on that axis.

- **The smoothed scatter panel is not an agreement statistic.** Both series are
  LOWESS-smoothed along the same domain order, so neighbouring smoothed values
  come from overlapping windows and are not independent samples: r goes
  0.64 -> 0.87 in 1978-1997 without either source having moved. The bias does
  *not* shrink under smoothing (+1.74 -> +1.95 m/yr), which is the tell - the
  offset is real and only the apparent scatter is manufactured. Both panels now
  carry their r; the smoothed one also carries the raw r and a count of
  effectively independent spans (a window holds `frac*n` domains, so there are
  about `1/frac` of them, ~7 here, not the 90 points plotted).

- **LOWESS extrapolates at the southern edge.** It is a local *linear* fit, so
  at the very edge it runs past the data: CoastSat 1978-1997 smoothed to
  -6.21 m/yr at domain 1, where the raw value is -0.59 and the whole local
  spread over domains 1-6 is -7.01..-0.59. `SKIP_SOUTHERN_DOMAINS = 10`
  withholds the smoothed curves there - the same guard, at the same width, as
  the hindcast's `LowessConfig(skip_southern_domains=10)`, and for the same
  reason: Oregon Inlet dominates that zone. Display only, as it is there. The
  fit still runs over all 90 domains, so the southern data still pulls the
  values just north of the cut; only the result is withheld, and the raw series
  still plots across the whole island. The smoothed-only figures, which would
  otherwise be blank there, draw the raw domain means as a dotted line. Set it
  to `0` to smooth everywhere.

`comparison_table_smoothed.csv` follows the figures rather than the fit: its
three `*_smooth` columns are blank across domains 1-10, so the table cannot be
read as evidence the plots withheld. Raw columns are complete at 90 of 90 -
except DSAS 1997-2019 at 83, which is that source's real coverage gap and not
the guard.

## What the hindcast actually runs

`scripts/cascade_pipeline/coastsat_lowess.py`, imported by
`scripts/hatteras_ms/HAT_hindcast_1984_2024.py` (section 8) and by the groin
sweep, the sensitivity plotter and `7-source-sink/2-calibrate/`. It reads the
same raw `transect_lrr_full.csv` and applies the same transect-first method
this stage chose, configured as:

    LowessConfig(window_domains=(7,), skip_southern_domains=10)
    TARGET_WINDOW = 7           # rate_comparison uses max(window_domains); 10 until 2026-09-28

`skip_southern_domains=10` is **display-only**: LOWESS still fits over all
transects and only the result is truncated across GIS 1-10, so the southern
transects still pull the values just north of the cut.

**Read `scripts/cascade_pipeline/coastsat_lowess.py` for what the runs do.**
There is no copy of it in this folder any more.

`hindcast_lowess_snapshot/coastsat_lowess.py` held one until 2026-09-22. It
called itself a read-only mirror, and keeping it accurate depended on someone
refreshing it by hand after every edit to the live module. Nobody did: by the
time it was removed it was 70 lines and twenty days behind, and nothing in the
tree could tell. A copy that drifts silently is worse than no copy, because it
answers the question wrongly instead of sending you to the source.

The `method_comparison/` figures were produced against the 2026-09-02 version,
if you need to know exactly what they ran on:

    git show f0b64cf1:scripts/input_prep/6-scr-smooth/hindcast_lowess_snapshot/coastsat_lowess.py

## How the method got here (v1-v4, deleted 2026-10-01)

The lineage was kept in `superseded_20260902/` until 2026-10-01, when it was
deleted (Hannah: delete, git keeps them); none of it ran against the current
tree (its paths still named `input_preperation`). To read one:

```
git log --diff-filter=D --oneline -- scripts/input_prep/6-scr-smooth/superseded_20260902/<file>
git show <commit>^:scripts/input_prep/6-scr-smooth/superseded_20260902/<file>
```

Oldest first:

The `v1`-`v4` order is reconstructed from file dates (2026-04-07, 05-14, 07-20,
08-06), not from anything the files themselves declare - so read the ordering as
likely, and the description as the reliable part.

| file | was | what it added |
|---|---|---|
| `loess_v1_domain_only.py` | `coastsat_smoothed_single LOWESS.py` | The start. Domain-averaged only, one hardcoded `frac = 0.111`, no sensitivity. |
| `loess_v2_transect_mode_switch.py` | `coastsat_smoothed_transect_domain_dottedline.py` | First transect-first implementation, behind a hand-set `SMOOTHING_MODE` switch. |
| `loess_v3_window_sensitivity.py` | `coastsat_smoothed_final.py` | Window sizes expressed in *domains* rather than a raw frac, plus `window_comparison.png`. Still domain-only. Named "final" - it was not. |
| `loess_v4_both_modes.py` | `coastsat_smoothed_transect_domain.py` | Drops the switch - both modes every run, windows in km. Tested 5/6/7/8 domains. No method-comparison figure yet. |

The current script is that last one plus the method-comparison plots and the
5/7/10 domain window set that the hindcast's `(7, 10)` came from.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### lowess_dsas_vs_coastsat.py

DSAS against CoastSat LRR per domain, raw and LOWESS-smoothed, for the two DSAS windows.

From the script's original header:

```text
DSAS vs CoastSat LRR Comparison — Hatteras Island  (SMOOTHED VERSION)
Extends the original comparison script with LOWESS smoothing and
multiple visualization options for collaborator review.

New outputs (saved to OUTPUT_DIR)
  overview_smoothed.png          – 2-panel both periods, raw + LOWESS overlay
  smoothed_only_comparison.png   – 2-panel both periods, smoothed lines only
  smoothing_sensitivity.png      – 3-panel showing frac=0.10, 0.15, 0.20 side by side
  combined_sources.png           – single panel combining both periods + both sources
  scatter_smoothed_1978_1997.png – scatter of smoothed values, period 1
  scatter_smoothed_1997_2019.png – scatter of smoothed values, period 2

All original outputs are also regenerated.

comparison_table_smoothed.csv carries the raw values for every domain, but its
LOWESS columns are blank across the southern boundary zone (domains
1..SKIP_SOUTHERN_DOMAINS), matching what the figures draw and why.

Smoothing method: LOWESS (locally weighted scatterplot smoothing)
  - Applied independently to each series (DSAS and CoastSat)
  - Preserves large-scale spatial patterns while removing per-domain noise
```

Notes that were in the code:

```text
pathlib must be imported before the CONFIG block because every path below is
built from PROJECT_BASE_DIR at module level.
```

```text
Anchored 2026-09-14: this named a home directory, or a tree renamed since.
Rule 5 of ORGANIZATION.md.
```

```text
ANCHORED, NOT TYPED. Every path below used to be an absolute literal: the
output one had lost its drive (str(_PATH_REPO / "scripts" / "input_prep" / "...")) and so wrote its
figures to C:\scripts\ instead of into the repository, and the input ones
still spelled the folder "input_preperation" and pointed at a CoastSat tree
that has since moved under 5-scr. Anchoring on the pyproject.toml at the repo
root makes all of them follow the checkout and survive this file changing
depth.
```

```text
--- CoastSat CSVs ---
The retired windows. old_time_periods/ became superseded_20260810/ and
these two had not resolved since.
```

```text
--- LOWESS bandwidth (fraction of data used per local fit) ---
0.10 = ~9 domains  → more local, preserves more variation
0.167 = ~15 domains → recommended default (1.5km smoothing window)
0.20 = ~18 domains → smoother, loses finer spatial patterns
```

```text
--- Southern boundary guard ---
Domains 1..N are dropped from the SMOOTHED curves. Oregon Inlet dominates
that zone, and LOWESS is a local linear fit, so at the very edge it
extrapolates: CoastSat 1978-1997 smooths to -6.21 m/yr at domain 1 where
the raw value is -0.59 and the local raw spread is -7.01..-0.59. The
smoothed value lands outside the data it claims to summarise.

Same guard, same width as the hindcast's
cascade_pipeline/coastsat_lowess.py: LowessConfig(skip_southern_domains=10).
DISPLAY ONLY, as it is there -- the LOWESS still fits over all 90 domains,
so the southern data still pulls the values just north of the cut; only the
result is withheld. Set to 0 to show LOWESS everywhere.
```

```text
--- Output ---
Products live under data/hatteras_init/<stage>/, beside every other
input_prep stage's output; only the scripts live under scripts/. Resolved
through hat_observed_rates.py since 2026-09-18, when the folder was renamed
from lowess_dsas_vs_coastsat_output/.
```

```text
Windows consoles default to cp1252, which cannot encode the arrows and
en-dashes in the closing summary -- the run died there after writing every
figure. UTF-8 here, matching lowess_method_comparison.py.
```

```text
HOUSE STYLE: one typeface and one palette across every figure in this
project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
found by searching upward (ORGANIZATION.md rule 5). This file drew in
matplotlib's defaults until 2026-09-17 -- it never called apply_style().
```

```text
Town labels sit in this band, measured down from the top of the axes; the
stats boxes start below it. Both are in axes fractions so they cannot drift
into each other when ylim changes.
```

```text
Backed in white like the town labels: on the difference panels the
mean-difference line runs straight through this text otherwise.
```

```text
txt = (f"Raw:     R²={s_raw['r2']:.2f}  RMSE={s_raw['rmse']:.2f} m/yr  "
f"Bias={s_raw['bias']:+.2f} m/yr\n"
f"Smoothed: R²={s_smooth['r2']:.2f}  RMSE={s_smooth['rmse']:.2f} m/yr  "
f"Bias={s_smooth['bias']:+.2f} m/yr")
ax.text(0.01, 0.97, txt, transform=ax.transAxes, fontsize=8.5,
va="top", family="monospace",
bbox=dict(boxstyle="round", fc="white", alpha=0.88, ec="0.7"))
```

```text
DSAS only: this is the presentation figure, and two dotted series
across the withheld zone read as clutter rather than as context.
```

```text
DSAS only, one line per period. Four dotted series across the
withheld zone was clutter on a panel already carrying four curves.
```

```text
Both panels carry their r. The smoothed one also carries the raw r
and a count of independent spans, because its own r is not a
like-for-like improvement: LOWESS runs along domain order, so
neighbouring smoothed values are built from overlapping windows and
are not independent samples. A window holds frac*n domains, so the
number of effectively independent spans is about 1/frac -- roughly
6 here, not the 90 points plotted.
```

```text
Same guard as the figures, so the table cannot be read as evidence
the plots withheld. The three *_smooth columns are therefore EMPTY
across domains 1..SKIP_SOUTHERN_DOMAINS; the raw columns are not,
and still cover the whole island.
```

<details><summary>Function notes (the original docstrings)</summary>

**`mask_southern_smoothed()`**

```text
Blank the smoothed columns across the southern boundary zone.

Applied after the fit, never before: the LOWESS still sees every domain,
matching how the hindcast splices this zone out. Raw columns are
untouched, so the raw series still plots across the whole island.
```

**`draw_raw_in_guard_zone()`**

```text
Raw domain values across the withheld zone, so it is not simply blank.

Matches the hindcast, which omits the LOWESS line across the southern
domains and shows the raw values there rather than nothing. Only the
smoothed-only figures need this; the others already draw raw everywhere.
```

**`shade_missing_dsas()`**

```text
Hatch the domains where DSAS has no data, so absence is not read as agreement.

The 1997-2019 DSAS export covers domains 8-90; without this the DSAS
series simply starts late while CoastSat runs the full span, which reads
as the two sources agreeing rather than as one of them being absent.

Runs are found from the data, so this stays correct if the DSAS coverage
changes. Only the first run is labelled, to keep one legend entry.
```

**`plot_overview_smoothed()`**

```text
2-panel figure: raw lines (faded) + LOWESS overlay (bold).
This is the recommended figure for collaborator review.
```

**`plot_smoothed_only()`**

```text
2-panel: LOWESS smoothed lines only, no raw data.
Cleanest version for presentations or dissertation figures.
```

**`plot_smoothing_sensitivity()`**

```text
3-panel showing effect of different LOWESS bandwidths.
Helps collaborators understand smoothing choice.
```

**`plot_combined_sources()`**

```text
Single panel: all 4 smoothed series together.
Good for seeing overall pattern and period differences simultaneously.
```

</details>

### lowess_method_comparison.py

CoastSat LRR smoothing along the island, transect-based and domain-averaged LOWESS side by side.

From the script's original header:

```text
CoastSat LRR Smoothing — Hatteras Island
Always runs both transect-based and domain-averaged smoothing.

  "domain"    Smooth pre-averaged domain LRR values (original approach).
              Input : domain_lrr_summary.csv — one row per CASCADE domain.

  "transect"  Smooth individual transect LRR values first, then aggregate
              the smoothed signal back to domain resolution for CASCADE.
              Input : transect_lrr_full.csv  — one row per CoastSat transect.

Within "transect" mode, TRANSECT_X_AXIS controls what the smoother uses as x:
  "transect_id"   sequential integer derived from sort order
  "along_coast_m" cumulative along-coast distance derived from domain position

along_coast_m is derived automatically from domain number if not present in
the CSV: each domain's transects are spread evenly across its 500 m band.
Physical spacing for the LOWESS frac is always estimated from along_coast_m.

Outputs — domain-space figures (both modes)
  overview_smoothed.png              raw + LOWESS overlay, both periods
  smoothed_only_comparison.png       clean version for presentations
  combined_periods.png               both periods on one panel
  smoothing_sensitivity_*.png        3-panel bandwidth sensitivity
  window_comparison.png              all window sizes overlaid
  coastsat_<mode>_smoothed_table.csv domain-level raw + smoothed values

Additional outputs — transect mode only
  transect_smoothed_overview.png     raw transect scatter + LOWESS in transect space
  transect_window_comparison.png     window sensitivity in transect space
  coastsat_transect_lrr_*.csv        full transect-level table with lrr_smooth
```

Notes that were in the code:

```text
pathlib must be imported before the CONFIG block because every path below is
built from PROJECT_BASE_DIR at module level.
```

```text
Anchored 2026-09-14: this named a home directory, or a tree renamed since.
Rule 5 of ORGANIZATION.md.
```

```text
ANCHORED, NOT TYPED. Every path below used to be an absolute literal: the
output one had lost its drive (str(_PATH_REPO / "scripts" / "input_prep" / "...")) and so wrote its
figures to C:\scripts\ instead of into the repository, and the input ones
still spelled the folder "input_preperation" and pointed at a CoastSat tree
that has since moved under 5-scr. Anchoring on the pyproject.toml at the repo
root makes all of them follow the checkout and survive this file changing
depth.
```

```text
── Smoothing x-axis ─────────────────────────────────────────
"along_coast_m" is recommended — keeps the physical window consistent.
"transect_id" is available but causes non-uniform spacing artefacts.
```

```text
── Domain-mode inputs ───────────────────────────────────────
Resolved through hat_observed_rates.py (2026-09-18), not typed.
```

```text
── Transect-mode inputs ─────────────────────────────────────
Point to your transect_lrr_full.csv files for each period
```

```text
── LOWESS window ─────────────────────────────────────────────
Physical window width in km — applies to both modes.
Converted to a frac automatically based on data resolution.
2.5 km = 5 domains | 3.5 km = 7 domains | 4.0 km = 8 domains
```

```text
── Southern boundary guard ──────────────────────────────────
Domains 1..N are dropped from the SMOOTHED series. LOWESS is a local linear
fit, so at the edge of the reach it extrapolates rather than smooths, and
Oregon Inlet dominates that zone anyway.

Same guard, same width as the hindcast's cascade_pipeline/coastsat_lowess.py:
LowessConfig(skip_southern_domains=10). Applied AFTER the fit, never before,
so the southern data still pulls the values just north of the cut - only the
result is withheld. Raw series are untouched and still cover the whole
island. Set to 0 to smooth everywhere.

This bites less here than in the domain-space scripts: smoothing runs at
transect resolution, ~10 points per domain, so the edge fit has far more
local support. It is applied for consistency with what the model is scored
against, not because this script showed the same excursion.
```

```text
── Geographic annotations ──────────────────────────────
This block used to hold a copy of the town spans, the village centres, the
piers, the groins and the Wimble Shoals zone, in domain units and again in
metres, with its own colours -- a second description of the island to keep
in step with scripts/site_layer/hatteras_site_config.py. It is gone: the village
shading now comes from the house helper town_bands() and everything else is
read from HATTERAS_ANNOTATIONS. See ANNOTATION HELPERS below.

The colours those figures use are set after the imports, with the house
style, since they are taken from it.
```

```text
── Output ───────────────────────────────────────────────────
Products live under data/hatteras_init/<stage>/, beside every other
input_prep stage's output; only the scripts live under scripts/. Resolved
through hat_observed_rates.py since 2026-09-18, when the folder was renamed
from lowess_method_comparison_output/.
```

```text
Windows consoles default to cp1252, which cannot encode the arrows and
en-dashes in the progress output -- the script died on its first status
line. UTF-8 here so it runs the same from PyCharm, a terminal or a
scheduled call.
```

```text
The house figure style and the site's annotation config are siblings in
scripts/, which is not on sys.path when this file is run from its own folder.
```

```text
One typographic and colour standard for every Hatteras figure, in
scripts/site_layer/hat_figure_style.py. This script used to set its own rcParams and
name its own hex colours; both are gone.
```

```text
Two periods drawn together are the house vintage pair: the earlier one red,
the later one blue.
```

```text
The band over the domains whose LOWESS is withheld, and the shade
town_bands() uses for a village span (repeated here only so the legend can
show a patch that matches it).
```

```text
Three LOWESS windows are compared on the sweep figures: grey, purple and
green from the house palette, three hues that also separate on luminance.
```

```text
Strip any filename prefix accidentally prepended to column names
e.g. "domain_lrr_summary.csvdomain_number" -> "domain_number"
```

```text
Derive along_coast_m from domain position if not present in CSV.
Each domain's transects are spread evenly across its 500 m band so that
physical spacing can be estimated for the LOWESS frac calculation.
```

```text
Masked on domain, not on along_coast_m: identical cut, and it carries
through aggregate_to_domains, whose per-domain mean of an all-NaN group
is NaN. So every domain-space figure and the CSV export inherit the
guard without a second mask.
```

```text
Where the villages, piers, groins and shoal zones are is settled in
scripts/site_layer/hatteras_site_config.py (HATTERAS_ANNOTATIONS). This script used to
carry its own copy of the spans, the village centres and their colours, so
the island had two descriptions of itself that had to be kept in step. The
village shading now comes from the house helper town_bands(); only the marks
town_bands does not draw -- shoal zones, piers, groins -- are added here, at
the positions and in the colours the site config gives them.
```

```text
DOMAIN-SPACE FIGURES
Works identically for both modes — receives a domain-level DataFrame
regardless of whether it came from load_domain_csv or aggregate_to_domains.
```

```text
TRANSECT-SMOOTHED WINDOWS IN DOMAIN SPACE
Smooths at transect level for each window, aggregates to domain
resolution, then plots with domain number and geographic annotations.
```

```text
METHOD COMPARISON — transect-based vs domain-averaged LOWESS
Both smoothing approaches overlaid for each window size.
```

```text
The two methods are the point of these four figures, so they carry the
colour: grey C["BASE"] is smoothing the domain means, the original
approach, and purple C["ACCENT"] is smoothing the individual transects,
the one under test. The window, where more than one is shown, is the line
style.
```

```text
MAIN
Always runs both transect-based and domain-averaged smoothing.
Outputs are organised into clearly labelled subfolders.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_domain_csv()`**

```text
Load domain-averaged LRR summary CSV. Returns standardised DataFrame or None.
Tolerates common column name variations and strips filename-prefix corruption
(e.g. column named "domain_lrr_summary.csvdomain_number" instead of "domain_number").
```

**`load_transect_csv()`**

```text
Load transect-level LRR CSV. Returns standardised DataFrame or None.

Handles:
  - String transect IDs (e.g. 'usa_NC_0032_0021') — sorted by domain then
    ID string, then replaced with a sequential integer (1, 2, 3 …)
  - Missing along_coast_m — derived by spreading each domain's transects
    evenly across its 500 m band (domain 1 → 0–500 m, domain 2 → 500–1000 m …)
  - Physical spacing for the LOWESS frac always estimated from along_coast_m
```

**`smooth_transect_df()`**

```text
Apply LOWESS in transect space. Returns copy of df with lrr_smooth column.

Physical spacing for the frac is always estimated from along_coast_m
(derived or real), so the window is correct in physical kilometres
regardless of whether transect_id or along_coast_m is the plot x-axis.
```

**`aggregate_to_domains()`**

```text
Average smoothed (and raw) transect values within each CASCADE domain.

Returns a domain-level DataFrame matching the domain CSV schema so all
domain-space plot functions work unchanged:
  domain        — CASCADE domain number
  cs_lrr        — mean of raw transect LRRs within the domain
  cs_std        — std  of raw transect LRRs within the domain
  cs_lrr_smooth — mean of smoothed transect LRRs within the domain
```

**`_domain_to_m()`**

```text
A GIS domain number as along-coast metres. Domain d occupies
[(d-1)*500, d*500) m -- the convention load_transect_csv uses -- so its
centre is (d-0.5)*500 and the edges of a span lo..hi fall out of the same
call. Nothing about the transect-space figures is measured separately.
```

**`_reference_marks()`**

```text
Shoal zones, piers and groins from the site config. `to_x` maps a GIS
domain number onto this panel's x units, so the same positions serve the
domain-space and the along-coast figures.
```

**`draw_raw_in_guard_zone()`**

```text
Raw domain means across the withheld zone, so it is not simply blank.

Matches what the hindcast does there: splice_lowess_with_raw_south omits
the LOWESS line across the southern domains and the raw values are shown
instead. On figures that already draw raw everywhere this adds nothing, so
it is called only from the smoothed-only ones.
```

**`plot_transect_overview()`**

```text
Raw transect scatter + LOWESS smoothed curve in along-coast space,
with domain-averaged LRR overlaid as a dashed line with open markers.
Shows how much variability domain averaging collapses.
```

**`plot_transect_windows_domain_space()`**

```text
For each window in COMPARE_WINDOWS_KM:
  1. Apply LOWESS at transect resolution
  2. Aggregate smoothed values to domain means
  3. Plot against CASCADE domain number with geographic annotations

Replaces plot_domain_window_comparison in transect mode so all curves
shown are transect-based — no domain-averaged smoothing is mixed in.
```

**`plot_transect_sensitivity()`**

```text
3-panel bandwidth sensitivity — transect mode.
Smooths at transect level for each window then aggregates to domains,
matching plot_transect_windows_domain_space exactly.
```

**`plot_method_comparison()`**

```text
For each window in COMPARE_WINDOWS_KM, plots both:
  — purple : transect-based LOWESS (smooth transects → aggregate to domains)
  — grey   : domain-averaged LOWESS (smooth domain means directly)
the window carried by the line style. Both in domain space.
```

**`plot_method_comparison_single()`**

```text
Single-window method comparison: transect-based vs domain-averaged LOWESS.
Shows one window size only so the two curves can be read clearly.
```

</details>
