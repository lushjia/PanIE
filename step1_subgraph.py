"""
filter redundancy by bubble puncture method from **all** edges 

wkdir = ""

Input:
# 1. all edge depth file raw count 
    ./{chr}.hprc-v2.0-mc-grch38.depth_per_edge.added_chroms.txt
# 2. chr gfa 
    ./$chr.gfa

Output: 
    independent edge file: (not in GRCh38 edges, and independent)
        {chr}.independent_edge_list.txt

Usage:
python step1_subgraph.py chr21

# split at articualtion points and save subgraph into files for next step 
"""

import sys 
import gzip
from itertools import chain
import re 
import networkx as nx
import numpy as np
import pandas as pd 

chr = sys.argv[1]

# construct graph from all edges read depth raw count file 
edge_raw_file = gzip.open(f'./{chr}.hprc-v2.0-mc-grch38.depth_per_edge.added_chroms.txt.gz', 'rt')

# create graph from edge raw count file 
# edge_dict = {} # edge: depth in samples
G = nx.Graph()
node_edge_dict = {} # node1, node2: edge

edge_ave_depth = {} # create a dictionary to store edge and average depth of the edge

for line in edge_raw_file:
    line = line.strip().split('\t')
    edge = line[1]
    depth = line[2:]
    # edge_dict[edge] = depth
    edge_ave_depth[edge] = sum([int(i) for i in depth])/430 # '>69345828<69345829': 21.89069
    nodes = re.split(">|<", edge)[1:]
    # smaller node first 
    if nodes[0] > nodes[1]:
        nodes[0], nodes[1] = nodes[1], nodes[0]
    G.add_edge(nodes[0], nodes[1])
    node_edge_dict[tuple([nodes[0], nodes[1]])] = edge 

# # create a dictionary to store edge and average depth of the edge
# edge_ave_depth = {}
# for edge in edge_dict:
#     edge_ave_depth[edge] = sum([int(i) for i in edge_dict[edge]])/430 # '>69345828<69345829': 21.89069

# read in GRCh38 edges from gfa file 
grch38_edge_file = open(f'./{chr}.gfa', 'r')
grch38_edge_set = set()  # to store GRCh38 edges # 102643231-102643232 ... 
for line in grch38_edge_file:
    if line.startswith('W'):
        line_split = line.strip().split('\t')
        if line_split[1] == "GRCh38":
            path = line_split[6]
            # path is like: >102643231>102643232>102643233>102643234>102643235
            path_edges = re.split(">|<", path)[1:]
            # recording all edges in the path
            for i in range(len(path_edges) - 1):
                # smaller node first 
                node1 = path_edges[i]
                node2 = path_edges[i+1]
                if node1 > node2:
                    node1, node2 = node2, node1
                edge = f"{node1}-{node2}"
                grch38_edge_set.add(edge)

# get all brdige edges from origianl graph 
ori_bridge_edges = set(nx.bridges(G))  # set of tuples (node1, node2) # len 11886

# independent edge file output 

# ori bridge edges - independent expect reference edges 
with open(f"{chr}.original_bridge_edge.txt", 'w') as original_bridge_edge_file:
    for edge in ori_bridge_edges:
        # pop edge from graph 
        G.remove_edge(*edge)
        # write into independent edge file if not in grch38_edge_set
        node1, node2 = edge
        # smaller node first
        if node1 > node2:
            node1, node2 = node2, node1
        edge_mod = f"{node1}-{node2}"
        if edge_mod not in grch38_edge_set:
            # write to file 
            edge_ori = node_edge_dict[(node1, node2)]
            original_bridge_edge_file.write(f"{edge_ori}\n")  # write original edge format >node1>node2
            original_bridge_edge_file.flush()

# decompose the remaining graph into biconnected component and save the subgraph 

# copy graph to avoid modifying the original graph
G_copy = G.copy()  # create a copy of the graph to avoid modifying the original graph

n = 1
for bicomponent in nx.biconnected_components(G_copy):
    subgraph = G.subgraph(bicomponent)
    subgraph = nx.Graph(subgraph)  # unfreeze the graph
    # if nodes number <50, then keep in graph and save as a whole later 
    if len(subgraph.nodes) < 100:
        continue  # skip small subgraphs
    else:
        # save into a gml file 
        nx.write_gml(subgraph, f"./{chr}/subgraph_{n}.graph.gml")
        ebunch = list(subgraph.edges())
        G.remove_edges_from(ebunch)
        n += 1
        # read it back 
        # G_loaded = nx.read_gml("graph.gml")
    
# remove the isolated nodes 
G.remove_nodes_from(list(nx.isolates(G)))
# save the remain graph as a whole
nx.write_gml(G, f"./{chr}/subgraph_remaining.graph.gml")
    




