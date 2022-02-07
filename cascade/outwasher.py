# For this class, we want to create bay to ocean flow and sediment transport
# Use Murray/Paola to carve channels through dune gaps and transport sediment?
import numpy as np

from barrier3d import Barrier3d

# --------------------------------------------------------------------------------------------------------------------
# 1. increase bay level by standard amount which will also change the interior domain in B3d
# this will only need to happen once, not every time step
# combine the dune domain with the interior domain and flip them so that the dune is at the back
# def combine_domain(interior_domain, dune_domain):
#     whole_domain = np.append(interior_domain, dune_domain, axis=0)
#     flipped_domain = np.flip(whole_domain, 0)
#     return flipped_domain
#
# # this and everything after will need to update every time step
# def bay_surge(flipped_domain, bay_level):
#     flipped_domain = flipped_domain - bay_level
#     return flipped_domain, bay_level

# --------------------------------------------------------------------------------------------------------------------
# 2. route water through interior: create channels in the interior domain to/through the dune domain using Murray/Paola
# seems like in the paper they track bed elevation and sediment transport (outputs of this function?)
# will we need to add some kind of force to get water to flow uphill, or will we need to submerge the entire interior?

# Currently happening:
# for each storm: erode dunes, route water from ocean to bay, use sed transport to route sand

#within our flow routing we are going to assume inundation overwash
# 1. set inundation rules
# 2. create water domain (970 in other code)- one in front of interior domain and one behind
# ------ sensitivity model for how many cells we need in front and behind
# ------ also assume no dune rn to see how flow routing works
# 3. create a wall of discharge water from the bay (constant mag)
# ------ later we will make this into a hydrograph

# set inundation regime rules
Rin = 0.1
Si = 0.001      # directional slope positive bc uphill
# how to find/calc this? what is the elevation of interior domain?
# what does the interior domain variable hold?
Slim = 0.25     # currently set in the yaml
n = 0.5         # currently set in the yaml

# creating the water domains
interior_domain = Barrier3d.InteriorDomain
bay_water = []
ocean_water = []
for i in range(len(interior_domain[0,:])):
    bay_water[0:20, i] = 10     # whatever magnitude we want this to be
    ocean_water[0, i] = 5       # whatever magnitude we want this to be
domain = ocean_water.append(interior_domain, 0)     # might need to flip interior domain first but I am not sure
domain = domain.append(bay_water, 0)

#  set discharge- currently set to teh Qdune value?
Qo = 0.5  # m3/hr

# water flow rules Murray
Qi = ((Qo - Rin)*abs(Si)**-n/sum(abs(Si)**-n))*(1-(abs(Si)/Slim))

# sediment transport Murray
Ki = 7.5e-6     # from paper
C = 0.72        # order of the average slope (Murray), supposed to be 10x average slope of barrier
m = 2           # >1 usually 2.5 (Murray), from yaml
Qsi = Ki*(Qi*(Si+C))**m

# --------------------------------------------------------------------------------------------------------------------
# 3. track sediment transport at the dune gaps using Nienhuis?
# only at the dune gap itself and does not alter the gap width

# --------------------------------------------------------------------------------------------------------------------
# 4. add a beach and deposit sediment (in different ways to see what happens?)

# other questions:
# will we allow any migration?
# will this be added into cascade so that it only occurs at a specific time interval?

# then we will make the class

# class Outwasher:
#     def __init__(self):  # add anything that will be an attribute to the class
#         self.____
#
#     def update(self, b3d):  # use the above functions to get outputs, will likely be using some b3d variables
#        self.___
