# revetment — a shore-parallel wall

Scoped 2026-10-12. This note holds the simple first version and the questions
still open. Keep it simple first: no pinned shoreline, no storm feedback, no
end erosion.

```
2-module-tests/            what the module does on its own; no Hatteras fit
    1-straight-coast/          the simple wall vs no wall, full CASCADE, 60 yr
    figures/                   how the module works, as schematics
```

**Where it stands (2026-10-12):** the simple version works on a straight coast
(`2-module-tests/1-straight-coast/README.md`). The wall is **repaired, never
fails** (Hannah, 10-12): once the gap closes, the shoreline is held at the wall
through CASCADE's nourishment path and the sand is tallied. Dune line, crest
floor and overwash all behave. But the held shoreline sticks out and BRIE
spreads its sand to the neighbours, so the wall needs a growing supply (~58
m³/m/yr per domain by the last decade) and the neighbours gain instead of
eroding. How to treat that is open (question 1 below).

## What a revetment does, and what the model can show

A revetment fixes the line it protects. The beach in front keeps retreating
until the shoreline reaches the wall, the sand behind it stops feeding the
coast, and overwash behind it stops, so the island can't roll over.

Barrier3D has no beach: the grid starts at the dune, and the beach is an
implicit fixed width (`BermEl / beta`) used only in storm runup. So the
"beach" in front of the wall is bookkeeping, the distance from the shoreline
to the wall (the **gap**), the same kind of number `BeachDuneManager` keeps
for nourished beaches.

## The simple first version

A schedule object applied by the runner each year before `cascade.update()`,
like the nourishment schedule. No change to CASCADE, BRIE or Barrier3D.

Inputs: armoured domains, build year, crest height above the berm (m),
starting gap (m).

Each year, in each armoured domain:

1. **Hold the dune line:** `barrier3d.dune_migration_on = False`, the switch
   `BeachDuneManager` already uses (`barrier3d.py:1988`).
2. **Floor the dune height at the crest:** set the domain's `DuneRestart` to the
   crest height (dam). A storm can't lower the dune below it; overwash still
   happens when runup tops it.
3. **Track the gap:** wall position minus `x_s`, recorded per year.
4. **Hold the shoreline at the wall once the gap closes** (the wall is
   repaired, never fails): put `x_s` back at the wall with
   `shoreface_nourishment` and sync BRIE with
   `update_brie_for_human_modifications`, as a nourishment does, and record
   the sand that takes.

Output: one table per run (year, domain, gap, shoreline put back, sand added, crest).

Not in this version, on purpose: the
narrowing beach strengthening storms, end erosion, blocking overwash entirely,
and a wall and nourishment in the same domain (both set `dune_migration_on`,
so setup refuses it).

Things to check while building: the `DuneRestart` floor is capped at `Dmax`
(`barrier3d.py:1480`), and it's one value per domain, so the wall covers the
whole 500 m domain.

## Open questions

1. **What happens when the shoreline reaches the wall?** DECIDED 10-12: the
   wall is repaired, so the shoreline is held at the wall (the first test let it
   fail; results kept in the test README). OPEN: the held shoreline sheds sand
   to its neighbours through BRIE's diffusion, so the sand to hold it grows and
   the neighbours gain instead of eroding. Options, not chosen: block the
   alongshore export from the wall domains (the groin's blocking idea), or keep
   it and report the hidden nourishment as the cost of holding the line.
2. **Straight coast first, or Hatteras?** DONE: straight coast first. The groin's problems only showed on a
   straight coast. The revetment acts on Barrier3D, so the test needs full
   CASCADE on a small straight coast (10–20 domains, 50 yr), not the BRIE
   emulator. Checks: the dune line holds, the gap falls at the background
   rate, overwash behind stops, total shoreline change across all domains
   matches a no-wall run, how the neighbours respond, and what happens at gap 0.
   Hatteras has no revetment to fit (NC hardening ban), and its cape, fills
   and source/sink rates would confound the result. **Lean: straight coast.**
3. **Starting gap.** Sets how many years before the shoreline reaches the wall. Options:
   Barrier3D's implicit width `BermEl / beta` (consistent with the storm
   runup), `BeachDuneManager`'s default 30 m, or dune line minus the CoastSat
   shoreline at the build year (real, but Hatteras only). Used the implicit
   width (30 m) in the test.
4. **Deferred: the narrowing beach strengthening storms.** Wanted eventually.
   Storm runup comes from the fixed storm series (Stockdon, slope 0.06), so the
   gap doesn't reach the storms yet. Hook: the per-domain `StormSeries` setter
   (`barrier3d.py:1939`), as the parked runup-by-slope experiment used. Rules:
   - A. Stockdon with a changing slope: no basis.
   - B. Lower the toe in front of the wall: add up to (berm − MHW) to Rhigh
     and Rlow as the gap goes to 0. No new parameters. **Lean.**
   - C. Runup on the structure once the beach is gone (EurOtop): needs Hs at
     the toe, which the storm series doesn't carry.
   - D. A tuned ramp: nothing at Hatteras to fit it to.

## Why not the existing sandbag code

`check_sandbag_need` (`cascade/roadway_manager.py:439`) is a dune-height floor
near the road, not a fixed line, and it can't run here: it writes a per-cell
`DuneRestart` that only Barrier3D's upstream `Benton-split-dune-update` branch
has, so on `hatteras/adopted` it would fail in the first managed year. Its reset
branch also writes 0.075 dam (0.75 m), 10x the default. Sandbags are off in
`HAT_hindcast_config`. Step 2 above is the working part of its idea.
