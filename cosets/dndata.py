"""
Graphs and Groups specific to the Dn problem.
"""
from typing import Tuple, Iterable, Callable
from itertools import product, chain, combinations
from functools import partial
import networkx as nx
import numpy as np
from lazytree import LazyTree
from sympy.combinatorics import PermutationGroup, Permutation
from sympy import binomial
from .schreier import make_tree, transposition
from .maxsat import maxsat_mis_tree

VEC = Tuple[int,...]

def small_weight(num: int, wgt: int) -> Iterable[VEC]:
    """
    All binary vectors of weight wgt of length n.
    """
    if wgt == 0:
        yield num * (0,)
    elif wgt == num:
        yield num * (1,)
    elif 0 < wgt < num:
        for elt in range(0,2):
            for rest in small_weight(num - 1, wgt - elt):
                yield (elt,) + rest

def dn_counting(num: int) -> int:
    """
    Count the number of edges per Veit Elser.
    """
    cnt = 2
    cnt += 4 * binomial(num-1,1)
    cnt += 4 * binomial(num-1,2)
    cnt += 2 * binomial(num-1,3)
    return (2 ** num) * cnt
                
def dn_neighbors(num: int) -> Iterable[VEC]:
    """
    Neighbors of 0 in the Dn graph    
    """

    yield from (delta + (num - 1) * (0,)
                for delta in [(0,1), (1,0)])
    yield from (pref + delta
                for pref in product(range(2), repeat=2)
                for delta in chain(small_weight(num - 1, 1),
                                   small_weight(num - 1, 2)))
    yield from (ndelta + delta
                for ndelta in [(0,0), (1,1)]
                for delta in small_weight(num - 1, 3))

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
    return nx.convert_node_labels_to_integers(gph,
                                              ordering='sorted')

def dn_group(num: int) -> PermutationGroup:
    """
    Construct the above permutation group.
    """
    trans = {elt: ind for ind, elt in enumerate(
        product(range(2), repeat=num + 1))}
    transpos = [(0, 1)] + [(_, _+1) for _ in range(2,num)]
    perms = list(map(partial(transposition, num+1),
                     transpos))
    bperms = []
    for perm in perms:
        # Find action on individual elements
        # The brute force way.  We can be more intelligent
        nperm = []
        # Action on n+1 tuples: which are subsets
        for elt in product(range(2), repeat=num+1):
            nperm.append(
                trans[tuple([elt[_] for _ in perm])])
        bperms.append(Permutation(nperm))

    return PermutationGroup(bperms)

def small_distance(num: int, dist: int) -> nx.Graph():
    """
    Graph: nodes - binary n-tuples
    edges: points of distance <= dist
    """
    gph = nx.Graph()
    for elt1, elt2 in combinations(product(range(2), repeat = num),2):
        if sum((_[0] ^ _[1] for _ in zip(elt1, elt2))) <= dist:
            gph.add_edge(elt1, elt2)
    return gph

def make_dn_tree(num: int) -> LazyTree:
    """
    Make the tree for Dn graph.
    """
    return make_tree(dn_graph(num), dn_group(num))

def dn_mis_tree(num: int,
                depth: int = 1,
                test: Callable[[LazyTree], bool] = lambda _: True,
                trace: int = 0,
                **kwds) -> Iterable[int]:
    """
    Solve the dn_graph MIS problem with the symmetry tree.
    """
    return maxsat_mis_tree(dn_graph(num),
                           dn_group(num),
                           depth,
                           test,
                           trace = trace,
                           **kwds)
