# Dynamic contact network

## Introduction
This repository is a tutorial to calculate **contact frequency data** from a MD simulation trajectory of a protein system and then use the frequency data to build a **dynamic contact network**. 

Contact frequency calculation are a routine calculation when analyzing a trajectory and tring to gain some basic statistically meaningful interactions, whereas dynamic contact network can help us use the network theory to understand how contacts are connected or grouped as a network based on the dynamics sampled in a trajectory. 

The dynamic network will only use contact frequencies as edges and C-alpha atoms of residues as nodes to build a network, and then construct a network for further betweenness centrality (importance of a node or edge, "how many times the shortest paths go through a specific node or edge"). Of note, this concept is different from the "traditional" [dynamic network analysis](https://luthey-schulten.chemistry.illinois.edu/tutorials/network/network_tutorial.pdf) using covariance matrix calculated based on the dynamic movement of residues in a simulation system.

## Tools and dependencies
1. [getcontacts](https://getcontacts.github.io/)
2. any dependencies required by getcontacts: numpy scipy expat matplotlib scikit-learn pytest pandas seaborn cython vmd-python

```python
conda install numpy scipy expat matplotlib scikit-learn pytest pandas seaborn cython

conda install -c conda-forge vmd-python

# networkx package: https://networkx.org/
pip install networkx
# networkx is supposed to require python > 3.8 version; but python 3.6 works fine
```
3. python==3.6 (only 3.6 was tested, other versions might also work just fine). 

4. VMD

## Contact frequency calculation

