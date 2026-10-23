# NC-12 relocation comparison — results, 2026-08-31

Six runs of `HAT_relocation_comparison.py` on the rebuilt matrix: three
source/sink presets × groin off and on.

This file exists because the comparison writes to
`output/comparisons/relocation_1984_2004/`, which `.gitignore` excludes — the
same reason `CALIBRATION_FIGURES.md` exists. The artefacts do not survive a
clone; the numbers and the reasoning should.

**Layout since 2026-09-17:** one tree, `output/comparisons/relocation/`:
`<start>_<end>/<version>/<preset>[_groin]/` for a set, `events/1999/<version>/`
for the cross-window report, `versions/v2_vs_v3/` for the cross-version one,
`standard_setback/` for the GIS 11 drowning figure. Paths quoted below in
dated sections are the paths of their day. The `v1/` sets (the six results
below, 2026-09-01, superseded topography) were deleted the same day; this
file is what survives of them, and five of the six regenerate from the
calibration tree. Every set's `report.txt` and `tables/*.csv` are tracked
since then, so the numbers no longer depend on this file alone.

**Regenerate with:**

```
python scripts/hatteras_ms/experiments/HAT_relocation_comparison.py --preset <name>
# groin arms need explicit arms, since ARM_SCENARIO_TOKENS pins "nogroin".
# Pass --preset too: it only labels the report header, and omitting it
# writes "preset zeroBE" above a pair of calibBE arms.
python scripts/hatteras_ms/experiments/HAT_relocation_comparison.py \
  --preset <preset> \
  --arm-a output/raw_runs/1984_2004/<preset>/HAT_1984_2004_<preset>_road_bdm_groin \
  --arm-b output/raw_runs/1984_2004/<preset>/HAT_1984_2004_<preset>_road_reloc_bdm_groin \
  --out  output/comparisons/relocation_1984_2004/<preset>_groin
```

---

## The six results

| preset | groin | ±2 yr | ±5 yr | false positives |
|---|---|---|---|---|
| zeroBE | off | 0.10 | **0.30** | 0/45 |
| zeroBE | on | 0.10 | 0.20 | 0/45 |
| edgeBE | off | 0.10 | 0.20 | 0/45 |
| edgeBE | on | 0.10 | 0.20 | 0/45 |
| calibBE | off | 0.30 | **0.60** | 0/45 |
| calibBE | on | **0.40** | 0.50 | 0/45 |

Recall is hits out of the 10 historical relocation domains (1989, GIS 84–87;
1999, GIS 9–14). **False positives stayed at 0/45 control domains in every
one of the six**, with a control margin of median 170 m / min 40 m — so no
improvement here comes from the trigger simply firing more often.

**The number to quote is calibBE, groin off, 0.60 at ±5 yr.** It is the best of
the six and the like-for-like comparison against the 2026-08-27 baseline.

---

## Against the 2026-08-27 archive

`output/superseded_20260828/comparisons/relocation_1984_2004/`, no-groin
arms, so directly comparable:

| preset | ±2 then → now | ±5 then → now |
|---|---|---|
| zeroBE | 0.10 → 0.10 | 0.20 → **0.30** |
| edgeBE | 0.00 → **0.10** | 0.00 → **0.20** |
| calibBE | 0.30 → 0.30 | 0.50 → **0.60** |

Every preset improved or held; none regressed. edgeBE moved off zero, having
previously matched nothing at any tolerance.

**What changed between the two dates is the source/sink field, not the groin** —
the staleness catch-up in `8f8851b` plus three convergence passes. The groin's
own contribution to the field was ≤ 0.5 m/yr and confined to D9–D17
(`daf5372`).

---

## Only calibBE responds, and that is predicted

`HAT_relocation_comparison.py`'s own docstring says so: *"edgeBE carries rates
on GIS 1 and 90 ONLY, so at every domain under test here edgeBE and zeroBE are
the same forcing... Only calibBE puts a background-erosion term on the
relocation domains, which makes a calibBE re-run the natural sensitivity test
once those source/sink terms are updated."*

That is what happened. calibBE is the only preset that moves materially, and it
moves in the direction the mechanism predicts.

There is also a direct link to the groin work: **the 1999 event is at GIS 9–14,
and the groin-aware recalibration moved the calibBE field at D9–D17.** The
groin sits at D5/D6, nowhere near the road, but its correction propagates to
the edge of the excluded groin zone — which is where the 1999 relocation
domains begin.

---

## The groin does not systematically change relocation recall

Across the three presets it moves recall by **at most one domain in either
direction**, with no consistent sign:

| preset | effect of turning the groin on |
|---|---|
| calibBE | **gains** GIS 10 at ±2; **loses** GIS 12 at ±5 |
| zeroBE | **loses** GIS 11 at ±5 |
| edgeBE | nothing — identical hits *and* identical miss lists |

On a 10-domain metric each domain is worth 0.10, so these are single-domain
threshold crossings, not demonstrated effects. Reporting calibBE's ±2 rise to
0.40 without its ±5 fall to 0.50 would be picking the tolerance.

### A prediction that failed, and what it implies

Before running the zeroBE groin arms the expectation recorded was that they
would be **identical** to zeroBE groin-off, because zeroBE puts no
background-erosion term at the relocation domains and the groin therefore has
no path to them.

They are not identical: GIS 11 drops out of ±5. There is a route that does not
go through the BE field — the dipole sits at D5/D6 and BRIE diffuses roughly
six domains over twenty years, putting D11 at the edge of its reach through
alongshore transport.

What makes that reading coherent rather than convenient is that **edgeBE does
not move at all**. edgeBE and zeroBE differ only at GIS 1 and 90, and edgeBE's
GIS 1 carries −42.6 m/yr. Under edgeBE that edge forcing already propagates
north through the same domains and dominates; under zeroBE it is absent, so the
groin is the only local perturbation and D11 becomes sensitive to it.

Hold it loosely: it is one domain crossing a threshold, which is the metric's
resolution rather than a measured effect.

---

## What the comparison can still not say

Unchanged from the script's own docstring, and worth repeating before any of
these numbers are quoted: CASCADE's trigger is purely geometric — the dune line
overruns the road. There is no storm damage, no cost and no maintenance
decision in it. A match means *"the modelled physics would have overrun NC-12
near that year"*, **not** *"NCDOT would have moved the road then."*

Provenance: all six reports carry a header naming both arms with their run
times, topography product and git commit. The no-groin arms ran 2026-08-31
10:57–10:59, the groin arms 11:30–11:32 (edgeBE, calibBE) and 15:12–15:17
(zeroBE), all on `1984-start/v1`.


---

## 2026-09-04 — by interior (the row-insert set), groin on, calibBE

Seven off/on pairs under `output/raw_runs/row-insert/<arm>/`, reports under
`output/comparisons/relocation_1984_2004/row-insert/<arm>/`, digest in
`output/experiments/row_insert_set/relocation/`. Versions per
`1984-start/2-domain-reconstruction-1984/DUNE_TOPO_VERSION_GUIDE.md`.

**All three of those locations were deleted on 2026-09-07**, together with the
layers v4–v8 themselves (Hannah's decision: keep only unmodified topography).
This table is the surviving record of the comparison; it cannot be re-run.

| arm | version | ±2 yr | ±5 yr | hits at ±5 | false positives |
|---|---|---|---|---|---|
| original (v1 + v1-era setbacks) | v1 | 0.30 | 0.40 | 84 85 86 10 | 0/45 |
| none (re-pick base) | v2 | 0.30 | 0.40 | 84 85 86 10 | 0/45 |
| measured-floor | v4 | 0.00 | 0.30 | 86 10 11 | 0/45 |
| median | v5 | 0.00 | 0.30 | 86 10 11 | 0/45 |
| platform | v6 | 0.10 | 0.30 | 86 10 11 | 0/45 |
| matched-crest | v7 | 0.00 | 0.30 | 86 10 11 | 0/45 |
| matched-nocrest | v8 | 0.10 | 0.30 | 86 10 11 | 0/45 |

`original` reproduces the calibBE groin-on row above (0.30 / 0.40) to the
domain, so the 08-31 numbers stand on the renumbered tree. The 1984 row insert
(v4-v8) moves every emergent relocation 3-10 years LATER (GIS 85 1985→1995,
84 1987→1999, 86 1985→1992, 10 1998→2003, 11 1992→2004): the control's early
hits at 84-86 were the setback-0 artefact, and the inserted interiors miss
honestly. The FILL changes timing by at most one year. Mean absolute error
under the inserts 5.6 yr, signed +5 (late). Full reading in the digest README.

---

## Layout change, 2026-09-09

`output/comparisons/relocation_1984_2004/` is now **version-first**:
`<dune-topo version>/<preset>[_groin]/`, the version read from the runs'
metadata by the script. The six sets above sit under `v1/`, marked superseded
(v1 stopped being CURRENT on 2026-09-04). `v2/calibBE_groin/` and
`v3/calibBE_groin/` are the like-for-like pair run on 2026-09-09 from
`output/raw_runs/version-pair/`; their numbers are in each set's `report.txt`
and `tables/confusion.csv` (every set: `report.txt`, `tables/`, one folder per place with `topography.gif` and `dune-and-road.gif`), and the broader v2-against-v3 figure is in the
reconstruction's `6-result/`. The folder README maps it.

---

## 2026-09-15: the 1999 event from a 1996 start

Four new sets, all on dune-topo `v2`, all groin off, and a cross-window
report. The 1984-2004 pair was RE-RUN on v2 for this (the calibration tree
is still v1), filed under `output/raw_runs/arms/version-pair/v2/1984_2004/`;
the 1996-2010 pair is the calibration tree's own, on island offsets v2.

```
python scripts/hatteras_ms/experiments/HAT_relocation_comparison.py --period 1996 --preset <zeroBE|edgeBE>
python scripts/hatteras_ms/experiments/HAT_relocation_comparison.py --preset <p> \
  --arm-a output/raw_runs/arms/version-pair/v2/1984_2004/<p>/HAT_1984_2004_<p>_road_bdm_nogroin \
  --arm-b output/raw_runs/arms/version-pair/v2/1984_2004/<p>/HAT_1984_2004_<p>_road_reloc_bdm_nogroin
python scripts/hatteras_ms/experiments/HAT_relocation_period_compare.py --presets zeroBE edgeBE --version v2
```

Outputs: `output/comparisons/relocation_1996_2010/v2/`, `relocation_1984_2004/v2/`
(zeroBE, edgeBE added beside calibBE_groin), `relocation_periods/1999_event/v2/`.

### The 1996 window never fires

| window | preset | 1999 domains relocated unaided | recall ±2 / ±5 (six 1999 domains) | control false positives |
|---|---|---|---|---|
| 1984-2004 | zeroBE | GIS 11 in 2004 (+5 yr) | 0.00 / 0.17 | 0/45 |
| 1984-2004 | edgeBE | GIS 10 in 1999 (0), GIS 11 in 1994 (−5) | 0.17 / 0.33 | 0/45 |
| 1996-2010 | zeroBE | none | 0.00 / 0.00 | 0/49 |
| 1996-2010 | edgeBE | none | 0.00 / 0.00 | 0/49 |

(The per-window reports score all ten historical domains for 1984-2004: zeroBE
0.10 / 0.30, edgeBE 0.10 / 0.20, unchanged from the v1 numbers above.)

**Why, in one number: retreat accumulated at the road before 1999.** Both
windows start GIS 9-14 at the same setbacks (40, 20, 10, 20, 20, 60 m; the 1996
file is the 1984 file with the 1989 event applied, which does not touch these
domains). From 1984 the free arm has 15 model years and erodes a mean 7-12 m of
it before 1999; from 1996 it has 3 years and erodes 0 m at every domain. The
closest 1996 approach is GIS 11, which loses its one remaining cell in 2001
and then sits at 0 m for nine years -- the trigger is strictly `< 0` and the
setback moves in whole cells, so a road on the dune line has not fired.

**The dune line stops because the modelled shoreline does.** Over 1996-2010
the free run moves the GIS 9-14 shorelines by 0.2-1.2 m landward in total
(rates −0.02 to −0.09 m/yr), where CoastSat has −1.1 to −2.4 m/yr, i.e. 15-33 m
or one to three cells. Under zeroBE/edgeBE nothing puts that erosion on the
interior domains, and the roadway manager removes every overwash and rebuilds
the dunes (12 rebuilds across the six domains), so the barrier does not
migrate either. This is the same mechanism as the 08-31 finding that only
calibBE moves recall -- and no calibBE field exists for 1996-2010 yet.

### The position check at 2004 cuts the other way

Mean |modelled − surveyed 2004 setback| over GIS 9-14, m:

| window | free arm | prescribed arm |
|---|---|---|
| 1984-2004 (2004 is the end) | 21 (zeroBE) / 32 (edgeBE) | 22 / 19 |
| 1996-2010 (2004 is year 8 of 14) | 19 / 16 | 29 / 32 |

The 1996 free arm is CLOSER to the surveyed 2004 road than the 1984 free arm,
because it has not eroded the setback that history did not erode either. But
the 1996 prescribed arm is WORSE than the 1984 one: the prescribed
displacement is added to the current setback, and from 1996 the current
setback is still the full 1984 value, so the road lands 10-60 m behind the
surveyed position (GIS 11: 80 m against 20 m). From 1984 the 15 years of
modelled retreat had consumed part of the setback first, so the same
displacement landed nearer the truth. The prescribed arm's accuracy is
therefore a property of the window as much as of the displacement -- and the
1984-2004 edgeBE prescribed arm on v2 drowns GIS 11 (the 20 m target clears
the drowning threshold by one cell; see hat_run.yaml).

### What this does and does not say

The 1996 window cannot reproduce 1999 emergently because the model needs more
than three years of its own (too-slow) retreat to reach a road 1-6 cells
behind the dune. That is a statement about the start year AND about the
missing interior erosion, and the two are separable only with a 1996 calibBE
field. The 1984 window is the one on which relocation skill should be quoted;
the 1996 window's value for the road is the position check, where its free
arm is the best of the four.
