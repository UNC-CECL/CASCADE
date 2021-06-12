import numpy as np

from CASCADE import Cascade


def test_initialize():
    cascade = Cascade(
        name="test",
        datadir="/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/",
        roadway_management_module=False,
        alongshore_transport_module=False,
        beach_nourishment_module=False,
        community_dynamics_module=False,
    )

    # check that the shoreface toe and shoreline are correct between the two models
    x_t_TS, x_s_TS, x_b_TS, h_b_TS = [[] for _ in range(4)]
    for iB3D in range(cascade.brie.ny):
        x_t_TS.append(cascade.barrier3d[iB3D].x_t_TS[0] * 10)
        x_s_TS.append(cascade.barrier3d[iB3D].x_s_TS[0] * 10)
        x_b_TS.append(cascade.barrier3d[iB3D].x_b_TS[0] * 10)
        h_b_TS.append(cascade.barrier3d[iB3D].h_b_TS[0] * 10)

    dt = cascade.brie.x_t - x_t_TS
    ds = cascade.brie.x_s - x_s_TS  # this isn't always zero; rounding error
    db = cascade.brie.x_b - x_b_TS
    dh = cascade.brie.h_b - h_b_TS

    return dt, ds, db, dh


def test_shoreline_variable_exchange_AST():

    cascade = Cascade(
        name="test",
        datadir="/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/",
        roadway_management_module=False,
        alongshore_transport_module=True,
        beach_nourishment_module=False,
        community_dynamics_module=False,
        time_step_count=3,
    )

    cascade.update()
    cascade.update()

    x_t_TS, x_s_TS, x_b_TS, h_b_TS = [[] for _ in range(4)]
    for iB3D in range(cascade.brie.ny):
        x_t_TS.append(np.array(cascade.barrier3d[iB3D].x_t_TS) * 10)
        x_s_TS.append(np.array(cascade.barrier3d[iB3D].x_s_TS) * 10)
        x_b_TS.append(np.array(cascade.barrier3d[iB3D].x_b_TS) * 10)
        h_b_TS.append(np.array(cascade.barrier3d[iB3D].h_b_TS) * 10)

    dt = cascade.brie._x_t_save - np.array(x_t_TS).astype(int)
    db = cascade.brie._x_b_save - np.array(x_b_TS).astype(int)
    ds = cascade.brie._x_s_save - np.array(x_s_TS).astype(int)
    dh = cascade.brie._h_b_save - np.array(h_b_TS)  # this isn't always zero; rounding error

    return dt, ds, db, dh


# other checks on the AST model?????