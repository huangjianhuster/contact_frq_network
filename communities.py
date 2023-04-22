import networkx as nx
from networkx.algorithms.community import girvan_newman
import pickle
import sys

# 
# Usage:
# python communities.py [frq.tsv] [number_of_communities] [output]
# example: python communities.py all_frq.tsv 6 6commm

# This function is from getcontacts
def create_graph(contact_frequency):
    """
    Create networkx graph from edge list with contact frequencies

    Parameters
    ----------
    contact_frequency: string
        Path to contact_frequency file containing weighted edge list in format of <res1> <res2> <frequency>

    Return
    ------
    graph: networkx object
        Residue interaction graph object

    """

    # Parse the contact_frequency graph file
    with open(contact_frequency, 'r') as f:
        nodes, edges = set(), set()
        for line in f:
            if not line[0] == "#":
                linfo = line.strip().split("\t")
                res1 = linfo[0]
                res2 = linfo[1]
                freq = float(linfo[2])

                nodes.add(res1)
                nodes.add(res2)
                edges.add((res1, res2, freq))

    # Construct networkx graph
    graph = nx.Graph()
    for res in nodes:
        graph.add_node(res)

    for res1, res2, freq in edges:
        graph.add_edge(res1, res2, weight=freq)

    return graph


contact_freq = sys.argv[1]
no_comm = int(sys.argv[2])
output = sys.argv[3]

# Create graph
G = create_graph(contact_freq)

# Use the Girvan-Newman algorithm to find the communities
communities = girvan_newman(G)

# Stop after no_comm communities have been detected
for i in range(no_comm-1):
    current_level_communities = next(communities)

# Print the number of communities and their members
n_communities = len(current_level_communities)
for i in range(n_communities):
    members = list(current_level_communities[i])
    print(f"Community {i+1} ({len(members)} nodes): {members}")

# write tuple to file, which could be further used for plotting
with open(f'{output}.pkl', 'wb') as f:
    pickle.dump(current_level_communities, f)

# pre-define the position of each node [avoid randomness in the latter plotting script]
# this is also for plotting
pos = nx.spring_layout(G, seed=50, dim=2, iterations=100)
with open(f'{output}_pos.pkl', 'wb') as f:
    pickle.dump(pos, f)

####### comments
# Other analysis might be interesting
# 1. Community evolution
# You can track the evolution of communities over time, or across different conditions or contexts, to understand how they change and adapt. For example, you could compare the communities identified in different time periods or in different experimental treatments to see how they differ.
# 2. Community visualization: You can visualize the communities in various ways to better understand their structure and organization. For example, you could use a force-directed layout or a circular layout to visualize the nodes within each community and their connections.


