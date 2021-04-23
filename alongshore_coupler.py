"""Alongshore Coupler

This module couples barrier3d with brie via alongshore sediment transport

References
----------

.. [1] Jaap H. Nienhuis, Jorge Lorenzo Trueba; Simulating barrier island response to sea level rise with the barrier
    island and inlet environment (BRIE) model v1.0 ; Geosci. Model Dev., 12, 4013–4030, 2019;
    https://doi.org/10.5194/gmd-12-4013-2019


Notes
---------

"""


class AlongshoreCoupler:
    #
    # """Couple brie and barrier3d for AST
    #
    # Examples
    # --------
    # >>> from cascade.alongshore_coupler import AlongshoreCoupler
    # >>> ast_coupler = AlongshoreCoupler()
    # >>> ast_coupler.update(brie, barrier3d)
    # """

    def __init__(self, ny=3):
        """The AlongshoreCoupler module.

        Parameters
        ----------
        ny: int, optional
            The number of alongshore Barrier3D domains for simulation in BRIE
        """

        self._ny = ny

    def update(self, brie, barrier3d, x_t_dt, x_s_dt, h_b_dt):

        # pass shoreline and shoreface values from B3D subdomains to brie for use in second time step
        brie.x_t_dt = x_t_dt
        brie.x_s_dt = x_s_dt
        brie.x_b_dt = 0  # dummy variable, will set x_b below
        brie.h_b_dt = h_b_dt

        # update brie one time step (this is time_index = 2 at start of loop)
        brie.update()

        for iB3D in range(self._ny):
            # pass shoreline position back to B3D from Brie (convert from m to dam)
            barrier3d[iB3D].x_s = brie.x_s[iB3D] / 10
            barrier3d[iB3D].x_s_TS[-1] = brie.x_s[iB3D] / 10

            # update dune domain in B3D (erode/prograde) based on shoreline change from Brie
            barrier3d[iB3D].update_dune_domain()

            # update back-barrier shoreline location in BRIE based on new shoreline + average interior width in B3D
            brie.x_b[iB3D] = barrier3d[iB3D].x_b_TS[-1] * 10
            brie.x_b_save[iB3D, brie.time_index - 1] = brie.x_b[iB3D]
