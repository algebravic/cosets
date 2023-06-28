"""
Generation of graphs, and writing them.
"""
from typing import Tuple, Iterable, Any
from itertools import product, chain
from sympy import binomial
import numpy as np
import networkx as nx
from .connection import dn_neighbors

def dn_graph(num: int, removal: int = 0) -> nx.Graph:
    """
    The Dn graph.
    Vertices are (z,b) where z in {-1,0,1,2}
    and b is in {0,1}^(n-1).
    Use the Gray embedding 0 -> 00, 1 -> 01, 2 -> 11, -1 -> 10
    """
    gph = nx.Graph()
    for eltx in dn_neighbors(num):
        for elt in product(range(2), repeat=num+1):
            nelt = np.array(elt, dtype=np.int8)
            gph.add_edge(elt, tuple((nelt ^ eltx).tolist()))
    if removal == 1:
        # Remove 0 and all neighbors
        zero = (num + 1) * (0,)
        nbrs = list(gph.neighbors(zero))
        gph.remove_node(zero)
        for nbr in nbrs:
            gph.remove_node(nbr)
    return gph

def independent(num: int) -> nx.Graph:
    """ Independent graph on n nodes"""
    gph = nx.Graph()
    for elt in range(num):
        gph.add_node(elt)
    return gph
