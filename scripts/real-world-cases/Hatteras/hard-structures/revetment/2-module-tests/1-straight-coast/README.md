# 1-straight-coast — the simple revetment, wall vs no wall

Run 2026-10-12. Does the simple wall do what it says, on a coast with nothing
else going on? The wall is **repaired, so it never fails** (Hannah, 10-12):
once the shoreline reaches it, the shoreline is put back at the wall each year.
A first version that let the wall fail at gap 0 was replaced before commit; its
results are summarised at the end.

    python straight_coast_revetment.py run       # both cases, ~90 s each
    python straight_coast_revetment.py figures

- `revetment_schedule.py`: the module (`RevetmentSchedule`), applied around
  each `cascade.update()`. It will move to `cascade/` if it's adopted.
- `straight_coast_revetment.py`: the two runs, the checks and the figure.
- `tables/`: every domain every year per case, the wall's own record, `checks.csv`.
- Runs: `output/raw_runs/experiments/revetment/2026-10-12-straight-coast/<case>/`.
- How the module works, as schematics: `../figures/`.

## Setup

Full CASCADE on the test-suite inputs (`tests/test_coupling`): 15 domains of
500 m, 60 years, RSLR 4 mm/yr, background erosion −1 m/yr everywhere so the
shoreline reaches the wall. Wall on domains 6–8 from year 1, crest floor 1.2 m
above the berm (Dmax here is 1.5 m), starting gap = Barrier3D's implicit beach
width, 30 m. Waves are the CASCADE test defaults (Hs 1 m, Tp 7 s, asymmetry 0.8,
high-angle fraction 0.2), not option A: option A's alongshore diffusivity is
negative at 0° and gets clipped to zero, so a straight coast would freeze and
the neighbours couldn't respond. Storms and dune growth are seeded, so the two
runs differ only by the wall.

## Results

| check | result |
|---|---|
| dune line held | yes: 0 dune-cell migrations in the wall domains (no wall: 15 in domain 7) |
| crest floor held | yes: dunes grow to Dmax (1.5 m); two big storms knock them to the 1.2 m floor (years ~28, ~45), and they regrow |
| gap closes, shoreline held | gap closes at ~1.7 m/yr, reaches 0 in years 17–18, then held at 0 |
| overwash behind | 355 m³/m in the wall domains vs 2,244 without (−84%), only from storms that top the crest |
| sand to hold the shoreline | ~2,000 m³/m per wall domain over 60 yr, 3.0 M m³ for the 1.5 km wall, and rising: 7, 40, 41, 56, 58 m³/m/yr per domain by decade |
| neighbours | no end erosion: neighbours end 13–84 m *seaward* of no-wall |

**The held shoreline feeds the coast.** The wall domains end 126 m seaward of
no-wall, sticking out from a coast that keeps retreating. BRIE's alongshore
transport works like diffusion, so a domain that sticks out sheds sand to its
neighbours, and the pin refills it every year. So the hidden nourishment ends
up feeding the neighbours, which is why the demand keeps rising and why the
two edge domains need more than twice the middle one's sand (2,400–2,500 vs
1,100 m³/m): they're the ones shedding sideways. This is the "protection, not
end erosion" outcome predicted before the run, and it's large.

**Not like a real wall.** A real seawall adds no sand: the beach in front lowers
and the neighbours erode (flanking). Here the sand is supplied and the
neighbours gain. Two ways to look at it later, not decided:
- block the alongshore export from the wall domains (the groin's blocking
  idea), so the pin only replaces cross-shore loss;
- represent the lowering in front through storms (open question 4) and report
  the hidden nourishment as the cost of holding the line.

**Shoreline sums aren't a sand budget here.** The pin put back 1,254 m of
shoreline (summed over domains and years); the coast's summed change is
+838 m against no-wall. The difference is not explained yet: pinning `x_s`
steepens the shoreface, and BRIE's row-scaled solve isn't conservative (see the
groin's straight-coast test). Volumes weren't measured.

**Starting gap.** In Barrier3D the dune sits at the shoreline and the beach is
implicit, always in front. Holding the dune while `x_s` retreats uses that
implicit beach up, so the gap starts at its width (30 m).

## First version: the wall failed at gap 0

Same setup, superseded 10-12. The wall failed in years 17–18 and the domain went
natural. While standing it slowed the wall domains' retreat (21 m ahead of
no-wall at year 17; less overwash, so less sand drawn off the shoreface; the
interior ended 18 m narrower). After failure BRIE spread that lead along the
coast: every domain ended 5–6 m seaward. The raised dune outlived the wall until
a storm in year 24.
