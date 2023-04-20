# Dynamic contact network

## Introduction
This repository is a tutorial to calculate **contact frequency data** from a MD simulation trajectory of a protein system and the use the frequency data to build a **dynamic contact network**. 

Contact frequency calculations are kind of a routine calculation when analyzing a trajectory and tring to gain some basic statistically meaningful interactions, whereas dynamic contact network can help us use network theory to understand how contacts are connected or grouped as a network. 

The dynamic network will be only using contact frequencies as edges and C-alpha atoms of residues as nodes to build a network, and then construct a network for further betweenness centrality (importance of a node or edge, "how many times the shortest paths go through a specific node or edge")

## Tools and dependencies
1. [getcontacts](https://getcontacts.github.io/)
2. any dependencies required by getcontacts: numpy scipy expat matplotlib scikit-learn pytest pandas seaborn cython vmd-python

```python
conda install numpy scipy expat matplotlib scikit-learn pytest pandas seaborn cython

conda install -c conda-forge vmd-python
```
3. python==3.6 (only 3.6 was tested, other versions might also work just fine)



## Contact frequency calculation

