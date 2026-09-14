# The 1984 seaward row-insert set — retired 2026-09-07, filed 2026-09-13

These four scripts drove and read an experiment that **cannot be re-run as
written**, because its subject was deleted.

    HAT_run_row_insert_set.py            the driver, one run per fill
    HAT_plot_row_insert_set.py           the in-model comparison
    HAT_digest_relocation_by_interior.py one table across the interiors
    HAT_gif_domain_by_interior.py        one domain, every interior, per year

## Why they are here rather than in experiments/

The set ran the 1984-2004 hindcast once per FILL of the seaward row insert,
on topography layers v3 to v8. Those layers, and every run built on them, were
deleted on 2026-09-07 when the decision was made to keep unmodified
extractions only. The driver's own arm table records it: seven arms cut to
two, with the note that a run on the others can no longer be built.

The evidence that nothing survives:

* `output/experiments/row_insert_set/` exists and holds zero files, as do the
  `relocation/` and `gifs/` subfolders the digest and the gif write into;
* the driver writes arms under the set name `row-insert`, and no such arm
  exists in the run tree.

## Why they were not deleted

They are the only record of how that experiment was driven, and each docstring
carries the reasoning and the retirement note. The conclusions it reached are
summarised where the decision was taken, in
`data/hatteras_init/1-barrier3d-domains/1984-start/2-domain-reconstruction-1984/`.

Re-running any of this means rebuilding the layers first, which was
deliberately given up.
