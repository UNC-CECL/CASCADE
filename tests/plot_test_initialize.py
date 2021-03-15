import numpy as np

from cascade import CASCADE as CASCADE


def test_initialize():
    brie, barrier3d = CASCADE.initialize(
        name="test",
        datadir="/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/",
    )

    # debugging: check that the shoreface toe and shoreline are correct between the two models
    x_t_TS, x_s_TS, x_b_TS, h_b_TS = [[] for _ in range(4)]
    for iB3D in range(brie.ny):
        x_t_TS.append(barrier3d[iB3D].x_t_TS[0] * 10)
        x_s_TS.append(barrier3d[iB3D].x_s_TS[0] * 10)
        x_b_TS.append(barrier3d[iB3D].x_b_TS[0] * 10)
        h_b_TS.append(barrier3d[iB3D].h_b_TS[0] * 10)

    dt = brie.x_t - x_t_TS
    ds = brie.x_s - x_s_TS  # this isn't always zero; rounding error
    db = brie.x_b - x_b_TS
    dh = brie.h_b - h_b_TS

    return dt, ds, db, dh


def test_shoreline_variable_exchange():

    brie, barrier3d = CASCADE.initialize(
        name="test",
        datadir="/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/",
        nt=30,
        ny=3,
    )

    # short time loop
    brie, barrier3d = CASCADE.time_loop(brie, barrier3d, num_cores=3)

    x_t_TS, x_s_TS, x_b_TS, h_b_TS = [[] for _ in range(4)]
    for iB3D in range(brie.ny):
        x_t_TS.append(np.array(barrier3d[iB3D].x_t_TS) * 10)
        x_s_TS.append(np.array(barrier3d[iB3D].x_s_TS) * 10)
        x_b_TS.append(np.array(barrier3d[iB3D].x_b_TS) * 10)
        h_b_TS.append(np.array(barrier3d[iB3D].h_b_TS) * 10)

    dt = brie._x_t_save - np.array(x_t_TS).astype(int)
    db = brie._x_b_save - np.array(x_b_TS).astype(int)
    ds = brie._x_s_save - np.array(x_s_TS).astype(int)
    dh = brie._h_b_save - np.array(h_b_TS)

    return dt, ds, db, dh
