🌀 🌊 🏄‍♀️ 🏚️ 
# CASCADE

The CoAStal Community-lAnDscape Evolution (CASCADE) model is a coupled landscape and human-dynamics modeling framework. CASCADE couples elements of several large-scale models -- each critical in understanding barrier response to a subset of the aforementioned processes -- into a single geomorphic model of barrier evolution. Barrier3D, a quasi-3D spatially-explicit cellular exploratory model, is the core of CASCADE. It is used within the CASCADE framework to simulate the effects of individual storm events and RSLR on shoreface evolution; dune dynamics, including dune growth, erosion, and migration; and overwash deposition by individual storms (Reeves et al., 2021). The BarrierR Inlet and Environment model (BRIE, Nienhuis and Lorenzo-Trueba 2019) is used to connect individual Barrier3D models with alongshore sediment transport. Human dynamics are incorporated in CASCADE in two separate modules. The first module simulates strategies for preventing roadway pavement damage during overwashing events, including rebuilding roadways at sufficiently low elevations to allow for burial by overwash, constructing large dunes, and relocating the road into the barrier interior. The second module incorporates management strategies for maintaining a coastal community, including beach nourishment, dune construction, and overwash removal. CASCADE represents decisions about coastal land-use (e.g., housing markets) and community-level mitigation measures using an empirically-grounded agent-based real estate model – the Coastal Home Ownership Model (CHOM). CHOM receives information about the coastal environment and acts on that information to cause change to the environment, including decisions about beach nourishment and dune construction and maintenance.


## Installation (work in progress, my rough notes for now)
- requires Python 3.8 (you can download directly here)
- we suggest using a virtual environment for the following steps:
- requires Barrier3D, brie, CHOM: each must be installed using "pip install -e ." in command line in the top of each directory
- then "pip install -e ." in CASCADE directory
- if you're not using something like PyCharm, extend your system paths to brie, CHOM, Barrier3D
- run model from CASCADE directory
