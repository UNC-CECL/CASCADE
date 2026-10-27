# Changelog for CASCADE

## Unreleased


- Added the groin as the ``GroinManager`` module (``cascade/groin.py``), configured
  through ``Cascade``'s ``groin_module`` and ``groin_*`` arguments and stepped in
  ``Cascade.update()`` before BRIE's alongshore transport solve. It replaces the
  groin object previously built by the Hatteras runner and attached to a Cascade
  from outside.


## 2.0.0 (2026-06-09)


- Added outwash processes. This is the version of *cascade* used for the simulations in
  "Barrier vulnerability following outwash: A balance of overwash and dune gap recovery"
  by Van Blunk et al., (in review).


## 1.0.0 (2024-02-09)


- Finalized for initial model publication. No changes to the code, but finalized to version 1.0.0.


## 0.1.1 (2023-03-21)

- Added additional tests for the ``rebuild_dunes`` function that uses that new
  *numpy* interpolator, ``RegularGridInterpolator``.
  (`#36 <https://github.com/UNC-CECL/cascade/issues/36>`_)


## 0.1.0 (2023-03-17)

- This is the version of *cascade* used for the simulations in "The Future of
  Developed Barrier Systems: Pathways Toward Uninhabitability, Drowning, and
  Rebound" by Anarde et al., (in review).
