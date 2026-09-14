###############################################################################
# run brie with LTA for comparison to CASCADE
###############################################################################

from brie import Brie


def LTA(
    name,
    wave_height,
    wave_period,
    asym_frac,
    high_ang_frac,
    slr,
    ny,
    nt,
    w_b_crit,
    h_b_crit,
    Qow_max,
):
    # update the initial conditions
    ast_model = True  # shoreface formulations on
    barrier_model = True  # LTA14 overwash model on
    inlet_model = False  # inlet model off
    b3d = False  # B3d overwash model on

    # barrier model parameters
    s_background = 0.001  # background slope (for shoreface toe position, back-barrier & inlet calculations)
    z = 10.0  # initial sea level (for tracking SL, Eulerian reference frame)
    bb_depth = 3.0  # back-barrier depth

    # inlet parameters (use default; these are here to remind me later that they are important and I can change)
    Jmin = 10000  # minimum inlet spacing [m]
    a0 = 0.5  # amplitude of tide [m]
    marsh_cover = 0.5  # % of backbarrier covered by marsh and therefore does not contribute to tidal prism

    # model setup
    dy = 100  # m, length of alongshore section (NOT the same as B3D, but overwash model performs better with dy=100 m)
    ny = ny * int(
        500 / dy
    )  # number of alongshore sections (NOTE, currently hard-coded for B3D dy = 500 m)
    dt = 0.05  # yr, timestep (NOT the same as B3D, but again, LTA14 performs better with dt = 0.05 yr)
    nt = int(nt / dt)  # equivalent timesteps to B3D
    dtsave = int(1 / dt)  # save spacing (equivalent of yearly for 0.05 time step)

    brieLTA = Brie(
        name=name,
        ast_model=ast_model,
        barrier_model=barrier_model,
        inlet_model=inlet_model,
        b3d=b3d,
        wave_height=wave_height,
        wave_period=wave_period,
        wave_asymmetry=asym_frac,
        wave_angle_high_fraction=high_ang_frac,
        sea_level_rise_rate=slr,
        sea_level_initial=z,
        barrier_height_critical=h_b_crit,
        barrier_width_critical=w_b_crit,
        max_overwash_flux=Qow_max,
        tide_amplitude=a0,
        back_barrier_marsh_fraction=marsh_cover,
        back_barrier_depth=bb_depth,
        xshore_slope=s_background,
        inlet_min_spacing=Jmin,
        alongshore_section_length=dy,
        alongshore_section_count=ny,
        time_step=dt,
        time_step_count=nt,
        save_spacing=dtsave,
    )  # initialize class

    for time_step in range(int(brieLTA.nt) - 1):
        brieLTA.update()

    return brieLTA
