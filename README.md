🌀 🌊 🏄‍♀️ 🏚️ 🌀 🌊 🏄‍♀️ 🏚️ 🌀 🌊 🏄‍♀️ 🏚️ 🌀 🌊 🏄‍♀️ 🏚️
# cascade

The CoAStal Community-lAnDscape Evolution (*cascade*) model is a coupled landscape and human-dynamics modeling framework.
*cascade* combines elements of two exploratory morphodynamic models of barrier evolution -- *barrier3d* 
(Reeves et al., 2021) and the BarrierR Inlet Environment (*brie*) model (Nienhuis & Lorenzo-Trueba, 2019) -- into a 
single model framework (Figure 1). *barrier3d*, a spatially-explicit cellular exploratory model, is the core of *cascade*. 
It is used within the *cascade* framework to simulate the effects of individual storm events and SLR on shoreface 
evolution; dune dynamics, including dune growth, erosion, and migration; and overwash deposition by individual storms. 
BRIE is used to simulate large-scale coastline evolution arising from alongshore sediment transport processes; this is
accomplished by connecting individual *barrier3d* models through diffusive alongshore sediment transport. Human dynamics 
are incorporated in *cascade* in two separate modules. The first module simulates strategies for preventing roadway 
pavement damage during overwashing events, including rebuilding roadways at sufficiently low elevations to allow for 
burial by overwash, constructing large dunes, and relocating the road into the barrier interior. The second module 
incorporates management strategies for maintaining a coastal community, including beach nourishment, dune construction, 
and overwash removal. For a full description of model dynamics, please see the pre-print of "The Future of Developed 
Barrier Systems: Pathways Toward Uninhabitability, Drowning, and Rebound" by Anarde et al., (in review) at _____.

[PLACEHOLDER FOR FIGURE 1]

In development: *cascade* represents decisions about coastal land-use (e.g., housing markets) and community-level 
mitigation measures using an empirically-grounded agent-based real estate model – the Coastal Home Ownership Model (*chom*). 
*chom* receives information about the coastal environment and acts on that information to cause change to the environment, 
including decisions about beach nourishment and dune construction and maintenance.

## Model coupling
*cascade* can initialize a series of *barrier3d* models, each describing a barrier segment with different initial conditions 
or management strategies (detailed below). The *barrier3d* segments are then coupled alongshore through a 
diffusive wave-driven sediment transport model (with periodic boundary conditions; i.e., Ashton & Murray, 2006) 
housed within the *brie* model, which distributes sediment alongshore amongst the different barrier segments. 
This coupling is possible because both models describe shoreface and shoreline dynamics using the formulations of 
Lorenzo-Trueba and Ashton (2014). Functionally, this coupling of *barrier3d*’s cross-shore morphodynamics with *brie*’s 
alongshore transport model requires 1) initializing both models with equivalent barrier geometry and environmental 
parameters, 2) separating dune migration within *barrier3d* from the other model processes in the one year time step 
(Figure 1), and 3) turning off all other model processes within *brie* 
(i.e., cross-shore barrier model and tidal inlet model). While the version of *barrier3d* in the *cascade* framework
produces equivalent results to the version used in Reeves et al., (2021; version testing is automated in *cascade*, 
see the `tests` folder), the default parameters are modified to match the shoreface configuration in *brie*, which depends 
on local wave and sediment characteristics as well as the offshore wave climate (Hallermeier, 1980; Ferguson & Church, 2004; 
Lorenzo-Trueba & Ashton, 2014; Ortiz & Ashton, 2016). For ease of model coupling, *brie* and *chom* were rewritten in Python 
and all models (*barrier3d*, *brie*, *chom*) were appended with a basic-model interface with the help of the 
Community Surface Dynamics Modeling System. The repositories for the models coupled within *cascade* are noted here:
- *barrier3d*: [GitHub Python Repository - Version 2.0 (BMI)](https://github.com/UNC-CECL/Barrier3D)
- *brie*: [GitHub Python Repository - Version 1.0 (BMI)](https://github.com/UNC-CECL/brie)
- *chom*: [GitHub Python Repository - Version 0.0.1.dev0 (BMI)](https://github.com/UNC-CECL/CHOM)

## Installation 
This ReadMe corresponds to the development version of *cascade* (version = 0.0.1.dev0) used for the simulations detailed in 
"The Future of Developed Barrier Systems: Pathways Toward Uninhabitability, Drowning, and Rebound" by Anarde et al., 
(in review). Prior to publication, *cascade* will be made available for easy download on `conda` and the following
instructions will be updated. Reviewers can follow the instructions provided below for installation of *cascade* (version = 0.0.1.dev0).

[PLACEHOLDER FOR NEW INSTRUCTIONS FROM ERIC]

1. Fork the *cascade*, *brie*, *chom*, and *barrier3d* repos on GitHub.
2. Clone your fork locally. As an example, for *cascade*::

    $ git clone git@github.com:your_name_here/cascade.git

3. Install your local copy into a conda environment. Assuming you have conda
   installed, this is how you set up your fork for local development::

    $ conda create -n cascade python
    $ conda activate cascade
    $ cd cascade/
    $ conda install --file=requirements.txt

    $ pip install -e .

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the
   tests::

    $ flake8 cascade
    $ pytest

   To get flake8, just conda install it into your environment.

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

## Example simulations
For a more complete set of example model runs and description of module functionality, we direct the use to the examples
provided in `notebooks`. 

Example (default) data inputs for cascade are provided in the `data` directory:
```
from cascade.cascade import Cascade

datadir = "cascade/data/"  
```
To initialize an instance of *cascade* with no human dynamics, 3 barrier segments (each 500-m long), and 
default *barrier3d* and *brie* parameters:
```
cascade = Cascade(
    datadir,
    name="no_human_dynamics_3_barrier_segments",
    alongshore_section_count=3,
    roadway_management_module=False,
    alongshore_transport_module=True,
    beach_nourishment_module=False,
    community_economics_module=False, 
)
```
To initialize an instance of *cascade* with roadway barrier management on 1 barrier segment:
```
cascade = Cascade(
    datadir,
    name="roadway_mgmt_1_barrier_segments",
    alongshore_section_count=1,
    roadway_management_module=True,
    alongshore_transport_module=False,
    beach_nourishment_module=False,
    community_economics_module=False, 
)
```
To initialize *cascade* with community barrier management on 1 barrier segment:
```
cascade = Cascade(
    datadir,
    name="community_mgmt_1_barrier_segments",
    alongshore_section_count=1,
    roadway_management_module=False,
    alongshore_transport_module=False,
    beach_nourishment_module=True,
    community_economics_module=False, 
)
```
Once initialized, a *cascade* time loop can be completed as follows:
```
for time_step in range(cascade.time_step_count - 1):
    cascade.update()
    if cascade.b3d_break:
        break
```
