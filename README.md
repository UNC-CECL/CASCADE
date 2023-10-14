DOI: 10.5281/zenodo.10003562

🌀 🌊 🏄‍♀️ 🏚️ 🌀 🌊 🏄‍♀️ 🏚️ 🌀 🌊 🏄‍♀️ 🏚️ 🌀 🌊 🏄‍♀️ 🏚️
# cascade

The CoAStal Community-lAnDscape Evolution (*cascade*) model is a coupled landscape and human-dynamics modeling framework.
*cascade* combines elements of two exploratory morphodynamic models of barrier evolution -- *barrier3d*
(Reeves et al., 2021) and the BarrierR Inlet Environment (*brie*) model (Nienhuis & Lorenzo-Trueba, 2019) -- into a
single model framework (figure below). *barrier3d*, a spatially-explicit cellular exploratory model, is the core of *cascade*.
It is used within the *cascade* framework to simulate the effects of individual storm events and SLR on shoreface
evolution; dune dynamics, including dune growth, erosion, and migration; and overwash deposition by individual storms.
*brie* is used to simulate large-scale coastline evolution arising from alongshore sediment transport processes; this is
accomplished by connecting individual *barrier3d* models through diffusive alongshore sediment transport. Human dynamics
are incorporated in *cascade* in two separate modules. The first module simulates strategies for preventing roadway
pavement damage during overwashing events, including rebuilding roadways at sufficiently low elevations to allow for
burial by overwash, constructing large dunes, and relocating the road into the barrier interior. The second module
incorporates management strategies for maintaining a coastal community, including beach nourishment, dune construction,
and overwash removal. For a full description of model dynamics, please see the pre-print of "The Future of Developed
Barrier Systems: Pathways Toward Uninhabitability, Drowning, and Rebound" by Anarde et al., (in review, [Earth ArXiv preprint](https://doi.org/10.31223/X5P947)).

![ModelTimeLoop-01](https://user-images.githubusercontent.com/57640439/226623608-d0c58437-d44f-4dca-8f43-0b92623fcda6.png)

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

This ReadMe corresponds to the development version of *cascade* used for the
simulations detailed in *"The Future of Developed Barrier Systems: Pathways Toward
Uninhabitability, Drowning, and Rebound" by Anarde et al., (in review, [Earth ArXiv preprint](https://doi.org/10.31223/X5P947))*. Prior
to publication, *cascade* will be made available for easy installation using either
`pip` or `conda`. Reviewers can follow the instructions provided below for installation
of *cascade*.

To install the latest release of *cascade* using *pip*, simply run the following
in your terminal of choice:

      pip install coastal-cascade

You can also use `conda`:

      conda install coastal-cascade

### From Source

*cascade* is actively being developed on GitHub, where the code is freely available.
If you would like to modifying code or contributing new code to *cascade*, you will first
need to get *cascade*'s source code, and then install *cascade* from that code.

To get the source code you can either clone the repository with *git*:

      git clone git@github.com/UNC-CECL/cascade

or download a [zip file](https://github.com/UNC-CECL/CASCADE/archive/refs/heads/main.zip):

      curl -OL https://github.com/UNC-CECL/CASCADE/archive/refs/heads/main.zip

Once you have a copy of the source code, you can install it into your current
environment,

      pip install -e .

We use [nox] to automate routine maintenance tasks like running the tests,
removing lint, etc. Install [nox] with *pip*::

      pip install nox

When you're done making changes, you can now run [nox] to check that the tests
pass and that there isn't any lint:

      nox -s test  # run the unit tests
      nox -s test-notebooks  # test that the notebooks run successfully
      nox -s lint  # find and, where possible, remove lint (black, flake8, etc.)

To run all of the above in a single command:

      nox

[nox]: https://nox.thea.codes/

## Example simulations
For a more complete set of example model runs and description of module functionality, we direct the use to the examples
provided in `notebooks`.

Example (default) data inputs for cascade are provided in the `data` directory:
```
from cascade.cascade import Cascade

datadir = "data/"
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
