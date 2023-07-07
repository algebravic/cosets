"""
Generation of graphs, and writing them.
"""
from typing import Tuple, Iterable, Hashable, FrozenSet, List
from itertools import product, chain
from sympy import binomial
import numpy as np
import networkx as nx

def remove_node_and_neighbors(ogph: nx.Graph, node: Hashable) -> nx.Graph:
    """
    Remove a node and its neighbors.

    If we want a maximum independent set containing node,
    this will be the same as a maximum indpendent set in
    the graph obtained by removing that node and its neighbors
    along with the original node.
    """
    gph = ogph.copy()
    gph.remove_nodes_from(list(gph.neighbors(node)) + [node])
    return gph

def truncate(gph: nx.Graph) -> nx.Graph:
    """
    Remove minimal vertex and all neighbors
    """
    return remove_node_and_neighbors(gph, min(gph.nodes))

def independent(num: int) -> nx.Graph:
    """ Independent graph on n nodes"""
    gph = nx.Graph()
    for elt in range(num):
        gph.add_node(elt)
    return gph

def heuristic_partition(gph: nx.Graph) -> List[FrozenSet[Hashable]]:
    """
    Apply the heuristic partition algorithm.
    """
    parts = {} # A dict key: clique vertices
    # value: all vertices not in clique
    # which are adjancent to all vertices in clique

    ngph = gph.copy()

    while ngph.nodes:
        # counts[node] will have the values (ncliques, deg(node))
        # where ncliques is the number of cliques to which it's
        # adjacent
        # Find the node with the minimal value
        choice = None
        cnode = None
        for node in ngph.nodes:
            val = (sum((int(node in adj)
                        for adj in parts.values())),
                   ngph.degree(node))
            if choice is None or val < choice:
                choice = val
                cnode = node
        neighbors = set(gph.neighbors(cnode))
        # in original graph, so as to contain all clique nodes
        ngph.remove_node(cnode)

        # Look to see if can be inserted into a clique
        fclique = None
        for clique in parts.keys():
            if neighbors.issuperset(clique):
                fclique = clique
                break
        new_neighbors = neighbors.intersection(ngph.nodes)
        if fclique is not None:
            newclique = fclique.union([cnode])
            newnbrs = parts[fclique].union(
                new_neighbors).difference(newclique)
            del parts[fclique]
            parts[newclique] = newnbrs
        else:
            parts[frozenset([cnode])] = new_neighbors

    return list(parts.keys())
