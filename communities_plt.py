import networkx as nx
import matplotlib.pyplot as plt
import mpld3
from mpld3 import plugins
import pandas as pd
import mplcursors
import pickle
import sys

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

# color list
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors)

# create graph from my file
contact_freq = sys.argv[1]
G = create_graph(contact_freq)

# read communities tuple from file
communities = sys.argv[2]
with open(f'{sys.argv[2]}.pkl', 'rb') as f:
   current_level_communities = pickle.load(f)

# read pre-computed positions for nodes
with open(f'{sys.argv[2]}_pos.pkl', 'rb') as f:
    pos = pickle.load(f)


######
# ANALYSIS --  for each community

# for annotation
# this is from: https://stackoverflow.com/questions/70340499/networkx-and-matplotlib-how-to-access-node-attributes-and-show-them-as-annotati
def update_annot(sel):
    node_index = sel.target.index
    node_name = list(G.nodes)[node_index]
    node_attr = G.nodes[node_name]
    text = node_name + ' ' + '\n'.join(f'{k}: {v}' for k, v in node_attr.items())
    sel.annotation.set_text(text)

# write out communities
def community_to_txt(subgraph, file):
    with open(file, "w") as f:
        for u,v,d in subgraph.edges(data=True):
            f.write(f"{u}\t{v}\t{d['weight']}\n")

# Get the size of each community
community_sizes = [len(community) for community in current_level_communities]
print(f"Community sizes: {community_sizes}")

# Draw graphs of each community
subgraphs = [G.subgraph(c) for c in current_level_communities]
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(20, 15))
ax = axes.flatten()

# pos = nx.spring_layout(G, seed=100)
# pos = nx.kamada_kawai_layout(G)
for i, community in enumerate(current_level_communities):
    # Extract the subgraph corresponding to the community
    subgraph = G.subgraph(community)

    # save communities
    community_to_txt(subgraph, f"{sys.argv[2]}-community-{i+1}.txt")
    # Calculate the degree for the nodes in the subgraph
    degree_dict = dict(subgraph.degree())
    pd.DataFrame.from_dict(degree_dict, orient="index").to_csv(f"{sys.argv[2]}_degree-community{i+1}.csv", header=False, sep="\t")
    # Calculate the betweenness for the nodes in the subgraph
    betweenness_dict = nx.betweenness_centrality(subgraph)
    pd.DataFrame.from_dict(betweenness_dict, orient="index").to_csv(f"{sys.argv[2]}_betweenness-community{i+1}.csv", header=False, sep="\t")

    # plot
    node_sizes = [1000 * betweenness_dict[node] for node in subgraph.nodes()]
    nodes = nx.draw_networkx_nodes(subgraph, pos, ax=ax[i], node_size=node_sizes, node_color='C'+str(i))
    # scale weights by 2
    nx.draw_networkx_edges(subgraph, pos, ax=ax[i],  width=[d['weight'] * 2 for (u, v, d) in G.edges(data=True)])
    
    # lables are too messy
    # labels = {node: node for node in subgraph.nodes()}
    # nx.draw_networkx_labels(subgraph, pos, labels, ax=ax[i], font_size=8, font_family='sans-serif')
    ax[i].set_title('Community {}'.format(i+1))
    cursor = mplcursors.cursor(nodes, hover=True)
    cursor.connect('add', update_annot)
plt.tight_layout()
plt.savefig(f"{sys.argv[2]}_sep.svg")
# plt.show()

# exit()

###  VISUALIZATION -- the whole graph
edge_betweenness = nx.edge_betweenness_centrality(G)
sorted_betweenness = sorted(edge_betweenness.items(), key=lambda x: x[1], reverse=True)
with open(f"{sys.argv[2]}_edge_betweenness.txt", 'w') as f:
    for edge, centrality in sorted_betweenness:
        node1, node2 = edge
        f.write(f"{node1}\t{node2}\t{centrality}\n")

# calculate cross community nodes
# create a new subgraph with edges connecting nodes from different communities
cross_community_edges = nx.Graph()
for edge in G.edges():
    node1, node2 = edge
    node1_community = None
    node2_community = None
    for i, community in enumerate(current_level_communities):
        if node1 in community:
            node1_community = i
        if node2 in community:
            node2_community = i
    if node1_community != node2_community:
        cross_community_edges.add_edge(node1, node2, weight=G.get_edge_data(node1, node2)['weight'])
with open(f"{sys.argv[2]}_crossing_nodes.dat", "w") as file:
    for (u, v, d) in cross_community_edges.edges(data=True):
        file.write(f"{u}    {v}    {d['weight']}\n")


# Set up the figure
fig, ax = plt.subplots(figsize=(10, 10))
# Set up the layout
# pos = nx.spring_layout(G, seed=42)
# pos = nx.kamada_kawai_layout(G)
# draw crossing nodes and edges first
nx.draw_networkx_edges(cross_community_edges, pos, width=[d['weight'] * 2 for (u, v, d) in cross_community_edges.edges(data=True)], edge_color='orange')
# Create a dummy scatter plot with fixed size for all communities
for i, community in enumerate(current_level_communities):
    ax.scatter([], [], s=50, c=f"C{i}", label=f"Community {i+1}")
# Loop over each community
for i, community in enumerate(current_level_communities):
    # Extract the subgraph corresponding to the community
    subgraph = G.subgraph(community)
    # Draw the nodes with different colors for each community
    betweenness_dict =  nx.betweenness_centrality(subgraph)
    node_sizes = [1000 * betweenness_dict[node] for node in subgraph.nodes()]
    nodes = nx.draw_networkx_nodes(subgraph, pos, ax=ax, node_size=node_sizes, node_color='C'+str(i) ) # , label=f"Community {i+1}")
    # Draw the edges with width proportional to the edge weights
    widths = [d['weight'] * 2 for (u, v, d) in subgraph.edges(data=True)]
    nx.draw_networkx_edges(subgraph, pos, ax=ax, width=widths)
    # Add labels for the nodes
    labels = {node: node for node in subgraph.nodes()}
    # nx.draw_networkx_labels(subgraph, pos, labels, font_size=10, font_family='sans-serif', hover=True)
    # nx.draw_networkx_labels(G, pos, labels=labels, font_size=10, font_color='black', alpha=0.7, ax=ax, bbox=dict(facecolor='white', edgecolor='none', alpha=0.5), horizontalalignment='center', verticalalignment='center', hover=True)
    cursor = mplcursors.cursor(nodes, hover=True)
    cursor.connect('add', update_annot)
# Set the title and show the plot
ax.set_title('Communities')
plt.legend()
plt.tight_layout()
plt.savefig(f"{sys.argv[2]}_toget.svg")
# plt.show()



####### comments
# Other analysis might be interesting
# 1. Community evolution
# You can track the evolution of communities over time, or across different conditions or contexts, to understand how they change and adapt. For example, you could compare the communities identified in different time periods or in different experimental treatments to see how they differ.
# 2. Community visualization: You can visualize the communities in various ways to better understand their structure and organization. For example, you could use a force-directed layout or a circular layout to visualize the nodes within each community and their connections.


