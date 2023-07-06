"""
Generation of graphs, and writing them.
"""
from typing import Tuple, Iterable, Hashable
from itertools import product, chain
from sympy import binomial
import numpy as np
import networkx as nx

def remove_node_and_neighbors(ogph: nx.Graph, node: Hashable):
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

def independent(num: int) -> nx.Graph:
    """ Independent graph on n nodes"""
    gph = nx.Graph()
    for elt in range(num):
        gph.add_node(elt)
    return gph
