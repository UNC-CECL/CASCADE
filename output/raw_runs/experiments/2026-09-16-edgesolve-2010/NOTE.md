# 2026-09-16-edgesolve-2010

**Question.** What should the two locked end domains carry in 2010-2024?
The period was wired 2026-09-11 with zeroBE only; this is its Newton solve,
done exactly as 1996 was (HAT_be_edge_domain_solve.py, base geometry, ends
GIS 1 and 90, model lrr_m_yr against target_lrr_m_yr, GIS 1 raw and GIS 90
LOESS-10). full_management, nourishment on, no groin, relocations off,
Hs 2.5. The 2010 matrix zeroBE runs preceded this and are the same run as
`base/` here; base is re-run under the experiment only so the solve script
can read stage 0 and each probe from one kind.

**Layout.**

```
base/2010_2024/zeroBE/<run>/     stage 0: nothing imposed at either end
step<k>/2010_2024/edgeBE/<run>/  the k-th probe, both ends through HAT_BE_OVERRIDE
SOLVED                           the step whose values went into HATTERAS_BE_EDGE_ONLY
logs/                            one log per run
```

**Runner change made for this.** An edgeBE run on a period with no edge
entry used to refuse before the override could act. Since 2026-09-16 it
starts from an empty mapping when HAT_BE_OVERRIDE is set, and the existing
guard still refuses a run with an end left at zero.
