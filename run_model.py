# Model runfile to couple barrier3d and brie
# Written by K.Anarde

from yaml import full_load, dump
import numpy as np
import matplotlib.pyplot as plt
import time

from barrier3d import Barrier3dBmi
from brie import Brie
import Barrier3D_Functions as B3Dplt
import CASCADE_plotters as CASCADEplt


# temporary placement - I want this added to the BMI but keep getting an object error
def set_yaml(var_name, new_vals, file_name):
    with open(file_name) as f:
        doc = full_load(f)

    doc[var_name] = new_vals

    with open(file_name, 'w') as f:
        dump(doc, f)


###############################################################################
# initial conditions for Brie
###############################################################################

# start by initializing brie because it has physical parameters related to wave climate that we will use to calculate
# params for B3D
brie = Brie()  # initialize class

# update the initial conditions
brie._name = 'bmiB3D'
brie._plot_on = False
brie._make_gif = False
brie._ast_model_on = True  # shoreface formulations on
brie._inlet_model_on = False  # inlet model off
brie._barrier_model_on = False  # overwash model off
brie._b3d_barrier_model_on = True  # B3d overwash model on

# wave climate parameters
brie._wave_height = 1  # m (lowered from 1 m to reduce k_sf)
brie._wave_period = 10  # s (lowered from 10 s to reduce k_sf)
brie._wave_asym = 0.8  # fraction approaching from left
brie._wave_high = 0.2  # fraction of waves approaching from higher than 45 degrees

# barrier model parameters (the following are needed for other calculations even if the barrier model is off)
brie._slr = 0.002  # m/yr
brie._s_background = 0.001  # background slope (for shoreface toe position, back-barrier & inlet calculations)
brie._z = 10  # initial sea level (for tracking SL, Eulerian reference frame)
brie._bb_depth = 3  # back-barrier depth

# model setup
brie._dy = 500  # m, length of alongshore section (same as B3D)
brie._ny = 6  # number of alongshore sections (6=3 km for testing AST, make 30 for inlets=15 km)
brie._dt = 1  # yr, timestep (same as B3D)
brie._nt = 25  # timesteps for 200 morphologic years
brie._dtsave = 1  # save spacing (every year)

# inlet parameters (use default)
brie._Jmin = 10000  # minimum inlet spacing [m]
brie._a0 = 0.5  # amplitude of tide [m]
brie._marsh_cover = 0.5  # % of backbarrier covered by marsh and therefore does not contribute to tidal prism

# get dependent variables
Brie.dependent(brie)

###############################################################################
# initial conditions for Barrier3D
###############################################################################

# when using the BMI, all values currently have to be updated in the yaml
#datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/B3D_Inputs/barrier3d-parameters.yaml"
datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/barrier3d-parameters.yaml"

# for each B3D subgrid, set the initial shoreface geometry equal to what is set in brie (some random perturbations);
# all other B3D variables are set equal
barrier3d = []

for iB3D in range(brie._ny):
    barrier3d.append(Barrier3dBmi())  # initialize each class of the BMI

    # update yaml file (these are the only variables that I'm like to change from default)
    set_yaml('TMAX', brie._nt, datadir)  # [yrs] Duration of simulation (if brie._dt = 1 yr, set to ._nt)
    set_yaml('BarrierLength', brie._dy, datadir)  # [m] Static length of island segment (comprised of 10x10 cells)
    set_yaml('DShoreface', brie._d_sf, datadir)  # [m] Depth of shoreface (set to brie depth, function of wave height)
    set_yaml('LShoreface', float(brie._x_s[iB3D] - brie._x_t[iB3D]),
             datadir)  # [m] Length of shoreface (calculate from brie variables, shoreline - shoreface toe)
    set_yaml('ShorefaceToe', float(brie._x_t[iB3D]), datadir)  # [m] Start location of shoreface toe
    set_yaml('BermEl', 2.0 , datadir) # [m] Static elevation of berm; berm elevation + dune height = dune elevation
    set_yaml('BayDepth', brie._bb_depth, datadir)  # [m] Depth of bay behind island segment (set to brie bay depth)
    set_yaml('MHW', brie._a0, datadir)  # [m] Elevation of Mean High Water (set to brie tidal amplitude?) ????????????
    set_yaml('DuneParamStart', True, datadir)  # Dune height will come from external file
    set_yaml('GrowthParamStart', True, datadir)  # Dune growth parameter will come from external file
    set_yaml('Dmaxel', 3.4, datadir)  # [m] Maximum elevation of dunes
    set_yaml('Rat', 0.0,
             datadir)  # [m / y] corresponds to Qat in LTA formulations (!!! must set to 0 because Qs is calculated in brie !!!)
    set_yaml('RSLR_Constant', True,
             datadir)  # Relative sea-level rise rate will be constant, otherwise logistic growth function used for time series
    set_yaml('RSLR_const', brie._slr, datadir)  # [m / y] Relative sea-level rise rate
    set_yaml('beta', 0.04, datadir)  # Beach slope for runup calculations
    set_yaml('k_sf', float(brie._k_sf),
             datadir)  # [m^3 / m / y] Shoreface flux rate constant (function of wave parameters from brie)
    set_yaml('s_sf_eq', float(brie._s_sf_eq),
             datadir)  # Equilibrium shoreface slope (function of wave and sediment parameters from brie)

    barrier3d[iB3D].initialize(datadir)

    # debugging: check that the shoreface toe, shoreline, back-barrier, and h_b are correct between the two models
    # brie._x_t[iB3D] - (barrier3d[iB3D]._model._x_t_TS[0] * 10)  # this isn't zero; rounding error
    # brie._x_s[iB3D] - (barrier3d[iB3D]._model._x_s_TS[0] * 10)

    # now update brie x_b, x_b_save[:,0], h_b, h_b_save[:,0] from B3D so all the initial conditions are the same
    brie._x_b[iB3D] = (barrier3d[iB3D]._model._x_b_TS[0] * 10)
    brie._h_b[iB3D] = (barrier3d[iB3D]._model._h_b_TS[0] * 10)
    brie._x_b_save[iB3D, 0] = brie._x_b[iB3D]
    brie._h_b_save[iB3D, 0] = brie._h_b[iB3D]

###############################################################################
# time loop
###############################################################################

# preallocate arrays for shoreline and barrier height change used in time loops
x_t_dt, x_s_dt, x_b_dt, h_b_dt = [
    np.zeros(np.size(barrier3d)).astype(float) for _ in range(4)
]

Time = time.time()

for time_step in range(brie._nt-1):

    # Print time step to screen
    print("\r", 'Time Step: ', time_step, end="")

    for iB3D in range(brie._ny):
        # advance B3D by one time step (this is time_index = 2 at start of loop)
        barrier3d[iB3D].update()

        # get values for passing to brie (all in dam)
        x_t_TS = barrier3d[iB3D]._values["shoreface_toe_position"]()
        x_s_TS = barrier3d[iB3D]._values["shoreline_position"]()
        x_b_TS = barrier3d[iB3D]._values["back_barrier_shoreline_position"]()
        h_b_TS = barrier3d[iB3D]._values["height_of_barrier"]()

        # calculate the diff in shoreface toe, shorelines (dam), height of barrier and convert to m
        x_t_dt[iB3D] = (x_t_TS[-1] - x_t_TS[-2]) * 10
        x_s_dt[iB3D] = (x_s_TS[-1] - x_s_TS[-2]) * 10
        x_b_dt[iB3D] = (x_b_TS[-1] - x_b_TS[-2]) * 10
        h_b_dt[iB3D] = (h_b_TS[-1] - h_b_TS[-2]) * 10

    # pass values from B3D subdomains to brie for use in second timestep
    # (there has to be a better way to do this with the BMI)
    brie._x_t_dt = x_t_dt
    brie._x_s_dt = x_s_dt
    brie._x_b_dt = x_b_dt
    brie._h_b_dt = h_b_dt

    # update brie one time step (this is time index = 2 at start of loop)
    brie.update()

    # loop to pass x_s and x_b (maybe this will be important for inlets with the x_b_fldt) back to B3D from Brie
    # (convert back to dam)
    for iB3D in range(brie._ny):
        barrier3d[iB3D]._model._x_s = brie._x_s[iB3D] / 10
        barrier3d[iB3D]._model._x_s_TS[-1] = brie._x_s[iB3D] / 10
        # barrier3d[iB3D]._model._x_b = brie._x_b[iB3D] / 10   # maybe will need this if tidal inlets on?
        # barrier3d[iB3D]._model._x_b_TS[-1] = brie._x_b[iB3D] / 10

SimDuration = time.time() - Time
print()
print('Elapsed Time: ', SimDuration, 'sec')  # Print elapsed time of simulation

###############################################################################
# run brie without B3D
###############################################################################

brieLTA = Brie()  # initialize class

# update the initial conditions
brieLTA._name = 'LTA'
brieLTA._plot_on = False
brieLTA._make_gif = False
brieLTA._ast_model_on = True  # shoreface formulations on
brieLTA._inlet_model_on = False  # inlet model off
brieLTA._barrier_model_on = True  # overwash model on
brieLTA._b3d_barrier_model_on = False  # B3d overwash model off

# wave climate parameters
brieLTA._wave_height = 1  # m (lowered from 1 m to reduce k_sf)
brieLTA._wave_period = 10  # s (lowered from 10 s to reduce k_sf)
brieLTA._wave_asym = 0.8  # fraction approaching from left
brieLTA._wave_high = 0.2  # fraction of waves approaching from higher than 45 degrees

# barrier model parameters (the following are needed for other calculations even if the barrier model is off)
brieLTA._slr = 0.002  # m/yr
brieLTA._s_background = 0.001  # background slope (for shoreface toe position & inlet calculations)
brieLTA._z = 10  # initial sea level (for tracking SL, Eulerian reference frame)
brieLTA._bb_depth = 3  # back-barrier depth

# also need these to be as similar as possible to "storm conditions" for B3D
brieLTA._w_b_crit = 450  # critical barrier width [m]
brieLTA._h_b_crit = 2 #(barrier3d[iB3D]._model._BermEl + barrier3d[iB3D]._model._MHW) * 10  # critical barrier height [m] (should equal B3D original BermEl above)
brieLTA._Qow_max = 20  # max overwash flux [m3/m/yr]

# model setup
brieLTA._dy = 100  # m, length of alongshore section (same as B3D)
brieLTA._ny = 60  # number of alongshore sections (10=6 km for testing AST, make 30 for inlets=15 km)
brieLTA._dt = 0.05  # yr, timestep (same as B3D)
brieLTA._nt = 100000  # timesteps for 200 morphologic years
brieLTA._dtsave = 1000  # save spacing (every year)

# inlet parameters (use default)
brieLTA._Jmin = 10000  # minimum inlet spacing [m]
brieLTA._a0 = 0.5  # amplitude of tide [m]
brieLTA._marsh_cover = 0.5  # % of backbarrier covered by marsh and therefore does not contribute to tidal prism

# get dependent variables
Brie.dependent(brieLTA)

# NOTE: the LTA model does not work well when you set the width and height from B3D
# now update brie x_b, x_b_save[:,0], h_b, h_b_save[:,0] from B3D so all the initial conditions are the same
#for iB3D in range(brie._ny):
#    brieLTA._x_b[iB3D] = (barrier3d[iB3D]._model._x_b_TS[0] * 10)
#    brieLTA._h_b[iB3D] = (barrier3d[iB3D]._model._h_b_TS[0] * 10)
#    brieLTA._x_b_save[iB3D, 0] = brieLTA._x_b[iB3D]
#    brieLTA._h_b_save[iB3D, 0] = brieLTA._h_b[iB3D]

    # debugging
    # brieLTA._h_b_save[:, 0] - brie._h_b_save[:, 0]
    # brieLTA._x_b_save[:, 0] - brie._x_b_save[:, 0]
    # brieLTA._x_t_save[:, 0] - brie._x_t_save[:, 0]
    # brieLTA._x_s_save[:, 0] - brie._x_s_save[:, 0] # need to turn off the random number generator

for time_step in range(int(brieLTA._nt) - 1):
    brieLTA.update()

###############################################################################
# plot
###############################################################################

# 4: Animation Frames of Barrier and Dune Elevation

#def plot_ElevAnimation(InteriorWidth_AvgTS, ShorelineChange, DomainTS, DuneDomain, SL, x_s_TS, PercentCoverTS,
#                       TMAX, DeadPercentCoverTS):

InteriorWidth_AvgTS = []
ShorelineChange = []
DomainTS = np.array()
DuneDomain = []
x_s_TS = []
SL = []
TMAX = []

for iB3D in range(brie._ny):
    InteriorWidth_AvgTS.append(barrier3d[iB3D]._model._InteriorWidth_AvgTS)
    ShorelineChange.append(barrier3d[iB3D]._model._ShorelineChange)
    DomainTS.append(barrier3d[iB3D]._model._DomainTS)
    DuneDomain.append(barrier3d[iB3D]._model._DuneDomain)
    x_s_TS.append(barrier3d[iB3D]._model._x_s_TS)
    SL.append(barrier3d[iB3D]._model._SL)
    TMAX.append(barrier3d[iB3D]._model._TMAX)

InteriorWidth_AvgTS = np.vstack(InteriorWidth_AvgTS).astype(float)
ShorelineChange = np.vstack(ShorelineChange).astype(float)
DomainTS = np.vstack(DomainTS).astype(float)

BeachWidth = 6
#OriginY = 10
OriginY = [10]
AniDomainWidth = int(max(InteriorWidth_AvgTS) + BeachWidth + abs(ShorelineChange) + OriginY + 35)  # was +15  # KA, just set for the first subgrid?

for t in range(TMAX):

    for iB3D in range(brie._ny):
        # Build beach elevation domain
        BeachDomain = np.zeros([barrier3d[iB3D]._model._BeachWidth, barrier3d[iB3D]._model._BarrierLength])
        berm = math.ceil(barrier3d[iB3D]._model._BeachWidth * 0.65)
        BeachDomain[berm:barrier3d[iB3D]._model._BeachWidth + 1, :] = barrier3d[iB3D]._model._BermEl
        add = (barrier3d[iB3D]._model._BermEl - barrier3d[iB3D]._model._SL) / berm
        for i in range(berm):
            BeachDomain[i, :] = barrier3d[iB3D]._model._SL + add * i

        # Make animation frame domain
        Domain = barrier3d[iB3D]._model._DomainTS[t] * 10
        Dunes = (barrier3d[iB3D]._model._DuneDomain[t, :, :] + barrier3d[iB3D]._model._BermEl) * 10
        Dunes = np.rot90(Dunes)
        Dunes = np.flipud(Dunes)
        Beach = BeachDomain * 10
        Domain = np.vstack([Beach, Dunes, Domain])
        Domain[Domain < 0] = -1
        AnimateDomain = np.ones([AniDomainWidth + 1, barrier3d[iB3D]._model._BarrierLength]) * -1
        widthTS = len(Domain)
        scts = [(x - barrier3d[iB3D]._model._x_s_TS[0]) for x in barrier3d[iB3D]._model._x_s_TS]
        if scts[t] >= 0:
            OriginTstart = OriginY + math.floor(scts[t])
        else:
            OriginTstart = OriginY + math.ceil(scts[t])
        OriginTstop = OriginTstart + widthTS
        #AnimateDomain[OriginTstart:OriginTstop, 0: barrier3d[iB3D]._model._BarrierLength] = Domain



    # Plot and save
    elevFig1 = plt.figure(figsize=(7, 12))
    ax = elevFig1.add_subplot(111)
    cax = ax.matshow(AnimateDomain, origin='lower', cmap='terrain', vmin=-1.1,
                     vmax=4.0)  # , interpolation='gaussian') # analysis:ignore
    ax.xaxis.set_ticks_position('bottom')
    # cbar = elevFig1.colorbar(cax)
    plt.xlabel('Alongshore Distance (dam)')
    plt.ylabel('Cross-Shore Diatance (dam)')
    plt.title('Interior Elevation')
    plt.tight_layout()
    timestr = 'Time = ' + str(t) + ' yrs'
    newpath = 'Output/SimFrames/'
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    plt.text(1, 1, timestr)
    name = 'Output/SimFrames/elev_' + str(t)
    elevFig1.savefig(name)  # dpi=200
    plt.close(elevFig1)

frames = []
for filenum in range(TMAX):
    filename = 'Output/SimFrames/elev_' + str(filenum) + '.png'
    frames.append(imageio.imread(filename))
imageio.mimsave('Output/SimFrames/elev.gif', frames, 'GIF-FI')
print()
print('[ * GIF successfully generated * ]')

# plot_ShorelinePositions(x_s_TS, x_b_TS):
B3Dplt.plot_ShorelinePositions(barrier3d[5]._model._x_s_TS, barrier3d[5]._model._x_b_TS)

#===================================================
# 1: CASCADE: shoreline transects

# plot_XShoreTransects(InteriorDomain, DuneDomain, SL, TMAX):
B3Dplt.plot_XShoreTransects(barrier3d[5]._model._InteriorDomain, barrier3d[5]._model._DuneDomain, barrier3d[5]._model._SL, barrier3d[5]._model._TMAX)

# plot_LTATransects(SL, TMAX, x_b_TS, x_t_TS, x_s_TS):
# B3Dplt.plot_LTATransects(barrier3d[5]._model._SL,
#                          barrier3d[5]._model._TMAX,
#                          barrier3d[5]._model._x_b_TS,
#                          barrier3d[5]._model._x_t_TS,
#                          barrier3d[5]._model._x_s_TS,
#                          barrier3d[5]._model._RSLR,
#                          barrier3d[5]._model._DShoreface,
#                          barrier3d[5]._model._BayDepth,
#                          barrier3d[5]._model._BermEl)

xmax = barrier3d[5]._model._x_b_TS[barrier3d[5]._model._TMAX - 1] + 20

SFfig = plt.figure(figsize=(20, 5))
colors = plt.cm.jet(np.linspace(0, 1, barrier3d[5]._model._TMAX))

for t in range(0, barrier3d[5]._model._TMAX, 5):  # Plots one transect every 5 years
    # Create data points
    Tx = barrier3d[5]._model._x_t_TS[t]
    Ty = ((barrier3d[5]._model._SL + (t * barrier3d[5]._model._RSLR[t])) - barrier3d[5]._model._DShoreface) * 10
    Sx = barrier3d[5]._model._x_s_TS[t]
    Sy = (barrier3d[5]._model._SL + (t * barrier3d[5]._model._RSLR[t])) * 10
    Bx = barrier3d[5]._model._x_b_TS[t]
    By = ((barrier3d[5]._model._SL + (t * barrier3d[5]._model._RSLR[t])) - barrier3d[5]._model._BayDepth) * 10
    Hx1 = Sx
    Hy1 = ((t * barrier3d[5]._model._RSLR[t]) + barrier3d[5]._model._BermEl) * 10
    Hx2 = Bx
    Hy2 = Hy1
    Mx = xmax
    My = By

    x = [Tx, Sx, Hx1, Hx2, Bx, Mx]
    y = [Ty, Sy, Hy1, Hy2, By, My]

    # Plot
    plt.plot(x, y, color=colors[t])

plt.xlabel('Alongshore Distance (dam)')
plt.ylabel('Elevation (m)')
plt.title('Shoreface Evolution')
plt.show()


#===================================================
# 1: Dune Height Over Time for CASCADE

# plot dune domain for all sub-grids
DuneCrest = []
Dmax = []

for iB3D in range(brie._ny):
    DuneCrest.append(barrier3d[iB3D]._model._DuneDomain.max(axis=2))
    Dmax.append(barrier3d[iB3D]._model._Dmax)

DuneCrest = np.hstack(DuneCrest).astype(float)
Dmax = np.max(Dmax)
duneFig = plt.figure(figsize=(14, 8))
plt.rcParams.update({'font.size': 13})
ax = duneFig.add_subplot(111)
ax.matshow((DuneCrest) * 10, origin='lower', cmap='bwr', aspect='auto', vmin=0, vmax=Dmax * 10)
cax = ax.xaxis.set_ticks_position('bottom')  # analysis:ignore
# cbar = duneFig.colorbar(cax)
# cbar.set_label('Dune Height Above Berm Elevation (m)', rotation=270)
plt.xlabel('Alongshore Distance (dam)')
plt.ylabel('Year')
plt.title('Dune Height (m)')
name = 'Output/Dunes'
# duneFig.savefig(name)

#===================================================
# CASCADE vs brie (LTA)

# sum Qoverwash for entire B3D grid
Qoverwash = np.zeros(np.size(barrier3d[iB3D]._model._QowTS)) # m^3/m/yr

for iB3D in range(brie._ny):
    Qoverwash = Qoverwash + (np.array(barrier3d[iB3D]._model._QowTS) * (barrier3d[iB3D]._model._BarrierLength * 10)) # m^3/yr

Qoverwash = Qoverwash / (brie._ny * brie._dy)
QoverwashLTA = brieLTA._Qoverwash / (brieLTA._ny * brieLTA._dy) # from brie in m^3/yr

plt.figure()
#plt.plot(Qoverwash)
plt.plot(QoverwashLTA)
fig = plt.gcf()
fig.set_size_inches(14, 5)
plt.xlabel('Time (yrs)')
plt.ylabel('Qow (m^3/m/yr)')
plt.title('Overwash Flux')
plt.legend(['CASCADE (Barrier3D)', 'BRIE (LTA14)'])
plt.show()
name = 'Output/Overwash'

# mean barrier width from brie+LTA (self._x_b - self._x_s)
plt.figure()
plt.plot(brieLTA._x_b_save[5, :]-brieLTA._x_s_save[5, :])

plt.figure()
plt.plot(np.array(barrier3d[5]._model._x_b_TS[:]) - np.array(barrier3d[5]._model._x_s_TS[:]))

# plot_BRIE_LTA


# plot_BRIE_frames
# CASCADEplt.plot_BRIE_frames(brieLTA._y,
#                             brieLTA._x_t_save,
#                             brieLTA._x_s_save,
#                             brieLTA._x_b_save,
#                             brieLTA._nt,
#                             brieLTA._ny,
#                             brieLTA._dy,
#                             brieLTA._dtsave,
#                             brieLTA._inlet_idx,
#                             'True',
#                             'test')

y = brieLTA._y
x_t = brieLTA._x_t_save
x_s = brieLTA._x_s_save
x_b = brieLTA._x_b_save
nt = brieLTA._nt
ny = brieLTA._ny
dy = brieLTA._dy
dtsave = brieLTA._dtsave
inlet_idx = brieLTA._inlet_idx
make_gif = 'True'
file_name = 'test'
directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"

# create time array
t = np.arange(0, nt, dtsave)

if make_gif:
    os.chdir(directory)
    os.mkdir(directory+'/GIF')

fig, axes = plt.subplots()
frames = []
ymax = int(math.ceil(np.max(x_b[:, -1]) / 100.0)) * 100
ymin = int(math.floor(np.min(x_t[:, 0]) / 100.0)) * 100

for i in np.arange(0, np.size(t)):

    # plot initial conditions
    axes.plot(y / 1000, x_t[:, i], color="cornflowerblue")
    axes.plot(y / 1000, x_s[:, i], color="gold")
    axes.plot(y / 1000, x_b[:, i], color="teal")
    axes.fill_between(y / 1000, ymin, x_t[:, i], color="royalblue")
    axes.fill_between(y / 1000, x_b[:, i], ymax, color="teal", alpha=0.6)
    axes.fill_between(y / 1000, x_t[:, i], x_s[:, i], color="cornflowerblue", alpha=0.6)
    axes.fill_between(y / 1000, x_s[:, i], x_b[:, i], color="gold", alpha=0.6)
    axes.legend(['x_t', 'x_s', 'x_b'])
    plt.xlabel('Alongshore (km)')
    plt.ylabel('Cross-shore (m)')
    plt.title('Time = ' + str(int(t[i])) + ' yr')
    axes.set_xlim(0, (ny-1) * dy / 1000)
    # axes.set_ylim(-2000, 6000)  # KA, will need to update later - placeholder
    # ratio = 0.4   # aspect ratio
    # axes.set_aspect(0.05)
    # axes.margins(0.5)
    # axes.set_aspect(1/axes.get_data_ratio())
    axes.set_ylim(ymin, ymax)

    # Here I chose to only mark the inlets (all indices), and not where the barrier volume is less than 0...
    # need to go back and fix
    # if np.size(inlet_idx) != 0 :
    #    axes.plot(y[np.hstack(inlet_idx)]/1000, x_s[np.hstack(inlet_idx)], 'kD')

    if make_gif:
        filename = 'GIF/year_' + str(int(t[i])) + '.png'
        fig.savefig(filename, dpi=200)

        # use imageio to create gif
        frames.append(imageio.imread(filename))
        os.remove(filename)

    axes.clear()  # alternatively, plt.pause(0.0001)

# make gif from saved png files
if make_gif:
    imageio.mimsave(name + '.gif', frames, 'GIF-FI')
    os.rmdir(directory + '/GIF')

###############################################################################
# save data
###############################################################################
name = 'CASCADE_BRIE_COMPARISON'
filename = name+'.npz'

# loop through and create separate variable names

#np.savez(filename, DuneDomain = DuneDomain, DomainTS = DomainTS, x_s_TS = x_s_TS, x_b_TS = x_b_TS, x_t_TS = x_t_TS, s_sf_TS = s_sf_TS, InteriorWidth_AvgTS = InteriorWidth_AvgTS, QowTS = QowTS, QsfTS = QsfTS, Hd_AverageTS = Hd_AverageTS,
#         PercentCoverTS = PercentCoverTS, DeadPercentCoverTS = DeadPercentCoverTS, ShrubArea = ShrubArea, ShrubDomainAll = ShrubDomainAll, ShrubDeadTS = ShrubDeadTS, StormCount = StormCount, t = t, RunUpCount = RunUpCount, InundationCount = InundationCount,
#         ShorelineChange = ShorelineChange, Hd_Loss_TS = Hd_Loss_TS, SimDuration = SimDuration, StormSeries = StormSeries, Dmax = Dmax, SL = SL, MaxAvgSlope = MaxAvgSlope, fluxLimit = fluxLimit, SimParams = SimParams)

#output = np.empty( (np.size(param1), np.size(param2)), dtype=object)
# next to do: make the cross-sections of brie versus CASCADE to see if this was a good example run, make colorplots to check if the shoreline change makes sense

np.savez(filename, brie=brie, brieLTA=brieLTA, barrier3d=barrier3d)

# blah blah blah


# FIGURES TO DO
# - make movies for brie and CASCADE
# - stats summary, Qsf, Qow, (a) shoreline position, average width, dune height (for both brie and CASCADE)
# - mean barrier width (brie and CASCADE) [future, with and without inlets]
# - barrier evolution (block diagram?, like Jaap Figure 1a) for brie and CASCADE (maybe 3 snapshots?)
# - cross shore transects for an example dy brie vs CASCADE, first and last time steps?
# - [eventually with inlets, Qow vs Qinlet]

# 1) highlight different processes in models
# - 500 yr & 1000 yr run with alongshore homogenous dune line, Qow = 40, K = 30,000 (1 m wave height, 8 sec period), critical barrier height = 2 m, critical barrier width = 200
# - make a run in brie to see if we get the expected behavior from LTA14

# 2) what is the effect of dunes (from Ian's paper)
# 3) what is the effect of the alongshore variability of dunes (30 km)
#   - vary the growth parameter by varying rmin and rmax, but keep difference (range) constant
#        - [rmin = 0.35, raverage = 0.6, and rmax = 0.85 everywhere as control case] with diffusive wave parameters (look at Brie paper to see what conditions are considered diffusive, or high angle)
#        - 2 B3Ds at raverage = 0.45 (0.3) and 2 B3Ds at raverage=0.75 (0.9), all along the barrier, check that raverage is 0.6 across the barrier
#   - show A outcome, but any conclusion need to come from averaging of different storm scenarios
#   - hyopthesis is that it will prevent punctuated retreat
# 4) what is the effect of dunes on relative importance of tidal inlets (need to check brie discretization)
