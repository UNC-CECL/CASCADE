# NC-12 relocation comparison — results, 2026-08-31

Six runs of `HAT_relocation_comparison.py` on the rebuilt matrix: three
source/sink presets × groin off and on.

This file exists because the comparison writes to
`output/comparisons/relocation_1984_2004/`, which `.gitignore` excludes — the
same reason `CALIBRATION_FIGURES.md` exists. The artefacts do not survive a
clone; the numbers and the reasoning should.

**Regenerate with:**

```
python scripts/hatteras_ms/HAT_relocation_comparison.py --preset <name>
# groin arms need explicit arms, since ARM_SCENARIO_TOKENS pins "nogroin":
python scripts/hatteras_ms/HAT_relocation_comparison.py \
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

`output/archived_output_20260828/comparisons/relocation_1984_2004/`, no-groin
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
