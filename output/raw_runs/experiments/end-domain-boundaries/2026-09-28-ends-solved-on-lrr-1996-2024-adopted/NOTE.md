# end-domain-boundaries/2026-09-28-ends-solved-on-lrr-1996-2024-adopted

**Question.** On the adopted model, what end rates at GIS 1 and 90 make the model match the full-period 1996–2024 CoastSat LRR in each window?

- Adopted model (2026-09-28): Barrier3D `hatteras/adopted` (the three overwash fixes and per-cell dune ceilings), storms `v3_trim24`. Otherwise option A: metres offset, Hs 2.0 / Tp 7.5 / asymmetry 0.6 / high-angle 0.5, full management, edgeBE, no groin.
- This is the 09-27 solve (`../2026-09-27-ends-solved-on-lrr-1996-2024-option-a/`), redone because the model changed (Hannah, 2026-09-28: "re-solve the LRR 1996-2024 ends on the adopted setup too").
- It supplies the `ends_solved_on_coastsat` set of `output/comparisons/target_comparison/projected/`. `target_comparison.py` reads it through `FULL_SOLVE_DIR`.

**Method.**

- Driver: `be_dune_edgesolve_loop.py --exp end-domain-boundaries/2026-09-28-ends-solved-on-lrr-1996-2024-adopted --windows 1996 2010 --target coastsat --coastsat-window 1996_2024 --tol 0.02 --max-steps 8`.
- The solve works as the matrix ends do: the model's OLS rate against the target table, GIS 1 against the raw domain mean, and GIS 90 against the LOWESS value. That LOWESS is **7 domains** now (`be_edge_domain_solve.py`, switched 2026-09-28); on 09-27 it was 10.
- The target is the 1996–2024 table, not each window's own.
- Start: the rebuilt matrix's zeroBE and edgeBE runs, on the adopted model.

**Answer** (`loop_log.csv`; `target_comparison.py` reads the last step):

| window | ends, GIS 1 / 90 (m/yr) | residuals (m/yr) | 09-27, pre-adoption |
|---|---|---|---|
| 1996–2010 | +3.9 / +27.7 | −0.047 / +0.001 | +4.5 / +27.5 |
| 2010–2024 | +3.6 / +18.6 | +0.038 / +0.010 | +4.8 / +20.5 |

**Both chains were accepted just outside the 0.02 tolerance.**

- The solver writes probes to 0.1 m/yr.
- From step 3 (1996) and step 4 (2010), it printed no further step: steps 4–8 re-ran the same probe.
- The misfits are the size of those accepted on 09-27 (+0.022 and +0.046).
- On the adopted model, 2010 GIS 1 responded at once: +4.06 → +0.33 in one step, where on 09-27 it needed seven.

**Runs.** `coastsat/step<k>/<window>/edgeBE/<run>/` stay on disk only (`.gitignore`): the paths pass git's 260-character limit. Logs are in `logs/`.
