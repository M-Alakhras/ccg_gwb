# network_tools.py
"""
functions to calculate network measures.
"""
from collections import Counter

import networkx as nx
import numpy as np
from networkx.algorithms.approximation import average_clustering, clique, connectivity
from networkx.classes.graph import Graph
from networkx.classes.multigraph import MultiGraph

all_correlation_network_measures = [
    "density",
    "connectivity",
    "large_clique_size",
    "average_clustering",
    #    'degree_assortativity_coefficient',
    "maximum_degree_centrality",
    "maximum_closeness_centrality",
    "maximum_betweenness_centrality",
    "maximum_eigenvector_centrality",
]
all_multilayer_network_measures = ["average_edge_overlap"]


def _create_network(matrix, labels):
    G = nx.from_numpy_array(matrix)
    G = nx.relabel_nodes(G, lambda x: labels[x])
    return G


def _create_multilayer_network(psrs):
    G1 = nx.visibility_graph(psrs[0].residuals)
    G = nx.MultiGraph()
    G.add_nodes_from(G1.nodes())
    for psr in psrs:
        G1 = nx.visibility_graph(psr.residuals)
        _ = G.add_edges_from(G1.edges, layer=psr.name)
    return G


def _thresholding_network(G, th):
    H = G.copy()
    for psr1, psr2, weight in G.edges(data=True):
        if weight["weight"] < th:
            H.remove_edge(psr1, psr2)
    return H


def _calculate_network_measures(G):
    measures = {}
    if isinstance(G, Graph):
        # network denisty
        measures.update({"density": nx.density(G)})
        # network connectivity
        measures.update({"connectivity": connectivity.node_connectivity(G)})
        # large clique size
        measures.update({"large_clique_size": clique.large_clique_size(G)})
        # average clustering
        measures.update({"average_clustering": average_clustering(G)})
        # maximum degree_centrality
        try:
            deg = nx.degree_centrality(G)
            max_deg = np.max(list(deg.values()))
        except:
            max_deg = np.nan
        measures.update({"maximum_degree_centrality": max_deg})
        # maximum closeness_centrality
        try:
            closeness = nx.closeness_centrality(G)
            max_closeness = np.max(list(closeness.values()))
        except:
            max_closeness = np.nan
        measures.update({"maximum_closeness_centrality": max_closeness})
        # maximum betweenness_centrality
        try:
            betweenness = nx.betweenness_centrality(G)
            max_betweenness = np.max(list(betweenness.values()))
        except:
            max_betweenness = np.nan
        measures.update({"maximum_betweenness_centrality": max_betweenness})
        # maximum eigenvector_centrality
        try:
            eig = nx.eigenvector_centrality(G)
            max_eig = np.max(list(eig.values()))
        except:
            max_eig = np.nan
        measures.update({"maximum_eigenvector_centrality": max_eig})
    elif isinstance(G, MultiGraph):
        # average edge overlap
        nlayers = len(G.nodes())
        top = np.sum(list(Counter(G.edges()).values()))
        bot = 0
        for i in range(len(G.nodes())):
            for j in range(i + 1, len(G.nodes())):
                if (i, j) in G.edges():
                    bot += 1
        w = top / (nlayers * bot)
        measures.update({"average_edge_overlap": w})

    return measures
