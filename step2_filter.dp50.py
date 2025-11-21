"""
speed up filter process by only check subgraph within 50 steps 

filter redundancy by bubble puncture method from **all** edges 

wkdir = ""

Input:
# 1. all edge depth file raw count 
    ./{chr}.hprc-v2.0-mc-grch38.depth_per_edge.added_chroms.txt
# 2. chr gfa -> GRCh38 edges
    ./$chr.gfa
# 3. subgraph file generated from step1 
    ./subgraph_100.graph.gml

Output: 
    independent edge file: (not in GRCh38 edges, and independent)
        {chr}.independent_edge_list.txt

python step2_filter.dp50.py chr21

# read in subgraph from step1 and apply filter to each subgraph
# step1: split at articualtion points and save subgraph into files for next step 
# step2: find independent/ representitive edge in each subgraph
# step2 compromise: only check dependent edge within 50 steps to speed up process
"""


import sys 
import gzip
from itertools import chain
import re 
import networkx as nx
import numpy as np
import pandas as pd 

chr = sys.argv[1]
graph_file = sys.argv[2]  # ./chr21/subgraph_100.graph.gml

# find prefic of graph file: subgraph_100
prefix = graph_file.split('/')[-1].split('.')[0]  # subgraph_100

# load graph 
G = nx.read_gml(graph_file)



#all edges read depth raw count file 
edge_raw_file = gzip.open(f'./{chr}.hprc-v2.0-mc-grch38.depth_per_edge.added_chroms.txt.gz', 'rt')

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
                edge = (node1, node2) # f"{node1}-{node2}"
                grch38_edge_set.add(edge) # (102643231,102643232)


# find edges in subgraph and find their depth 
# read in edges that in subgraph 
graph_edge_list = set(G.edges())
edge_ave_depth = {} # create a dictionary to store edge and average depth of the edge
node_edge_dict = {} # node1, node2: edge
# for edge in edge_raw_file, if in graph_edge_list, then save the edge and its depth
for line in edge_raw_file:
    line = line.strip().split('\t')
    edge = line[1]
    depth = line[2:]
    # if edge in graph_edge_list, then save the edge and its depth
    nodes = re.split(">|<", edge)[1:]
    # smaller node first 
    if nodes[0] > nodes[1]:
        nodes[0], nodes[1] = nodes[1], nodes[0]
    if (nodes[0], nodes[1]) in graph_edge_list or (nodes[1], nodes[0]) in graph_edge_list:
        edge_ave_depth[edge] = sum([int(i) for i in depth])/430 # '>69345828<69345829': 21.89069
        node_edge_dict[(nodes[0], nodes[1])] = edge # ('69345828','69345829'): '>69345828<69345829'
    else:
        continue
# find edges to consider (not in GRCh38 edges) 
edge_to_consider = set(node_edge_dict.keys()) - grch38_edge_set # ('69345828','69345829')
# order edges by average depth, descending
edge_to_consider = sorted(edge_to_consider, key=lambda x: edge_ave_depth[node_edge_dict[x]], reverse=True)
# loop through edges in edge_to_consider, and find independent edges
with open(f'./{chr}/{chr}.{prefix}.independent_edge_list.dp50.txt', 'w') as out_file:
    while edge_to_consider:
        # pop the first edge 
        edge = edge_to_consider.pop(0)  # ('69345828','69345829')
        # Get the local neighborhood (BFS) around both endpoints: depth < 50
        start_node = edge[0]
        local_nodes = set(nx.single_source_shortest_path_length(G, start_node, cutoff=50).keys())
        # get local subgraph 
        local_G = G.subgraph(local_nodes)
        local_G = nx.Graph(local_G) # unfreeze the graph 
        # remove bridge edges in the local subgraph
        local_bridge_edges = set(nx.bridges(local_G))
        # delete edge if edge is in local bridge edges
        local_bridge_edges.discard(edge)  # remove the edge itself if it is a bridge edge
        local_bridge_edges.discard((edge[1], edge[0]))  # remove the reverse edge if it exists
        local_G.remove_edges_from(local_bridge_edges) 
        # remove the edge from the local graph 
        local_G.remove_edge(*edge)  # remove the edge itself
        G.remove_edge(*edge)  
        # write edge to output file
        out_file.write(f"{node_edge_dict[edge]}\n")  # '>69345828<69345829'
        out_file.flush()
        # find bridge edges in the remaining subgraph (dependent edges)
        graph_bridge_edges = list(nx.bridges(local_G)) 
        # remove bridge edges from edge_to_consider and graph 
        for bri_edge in graph_bridge_edges:
            G.remove_edge(*bri_edge)  # remove from graph
            if bri_edge in edge_to_consider:
                edge_to_consider.remove(bri_edge)
            elif (bri_edge[1], bri_edge[0]) in edge_to_consider:
                edge_to_consider.remove((bri_edge[1], bri_edge[0]))
            else:
                continue


            
