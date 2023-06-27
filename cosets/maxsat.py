"""
Use Max Sat for independent set.
"""
from typing import List, Tuple, Iterable, Any
from itertools import product
import networkx as nx
import numpy as np
from pysat.formula import WCNF, IDPool
from pysat.examples.rc2 import RC2
from .connection import dn_neighbors

def maxsat_mis_model(gph: nx.Graph) -> Tuple[WCNF, IDPool]:
    """
    Simple maxsat formulation.
    """
    cnf = WCNF()
    pool = IDPool()
    for node1, node2 in gph.edges:
        cnf.append([-pool.id(('x', node1)),
                    -pool.id(('x', node2))])
    for node in gph.nodes:
        cnf.append([pool.id(('x', node))], weight=1)
    return cnf, pool

def solve_maxsat(cnf: WCNF, pool: IDPool,
                 stem: str = 'x', **kwds) -> Iterable[Any]:
    """
    Solve maxsat
    """
    solver = RC2(cnf, **kwds)
    soln = solver.compute()
    pos = [pool.obj(_) for _ in soln if _ > 0]
    answer = [_[1] for _ in pos if _ is not None and _[0] == stem]
    print("Time = {}".format(solver.oracle_time()))
    return answer

def maxsat_mis(gph: nx.Graph, **kwds) -> Iterable[Tuple[int, ...]]:
    """
    Independent sets in a graph via Max Sat
    """
    cnf, pool = maxsat_mis_model(gph)
    return solve_maxsat(cnf, pool, stem = 'x', **kwds)

def breaking_model(num: int) -> Tuple[WCNF, IDPool]:
    """
    Model for independent set in Dn graph with some symmetry
    breaking constraints.
    """
    cnf = WCNF()
    pool = IDPool()
    for delta in dn_neighbors(num+1):
        for elt in product(range(2), repeat=num + 1):
            nelt = np.array(elt, dtype=np.int8)
            other = tuple((nelt ^ delta).tolist())
            if other > elt:
                cnf.append([-pool.id(('x', elt)),
                            -pool.id(('x', other))])
    for elt in product(range(2), repeat=num + 1):
        cnf.append([pool.id(('x', elt))], weight=1)
    # Now put in the symmetry breaking constraints
    cnf.append([pool.id(('x', (num + 1) * (0,)))])
    # The last n-1 bits can be permuted in any way
    # When we have chosen the first k possibilities
    # the bit pattern for each is right shifted in
    # sectors.
    candidates = [ _ * (0,) + (num - 1 - _) * (1,)
                  for _ in range(num - 1)]
    front = [(0,0), (0,1), (1,1)]
    possible = [_[0] + _[1] for _ in product(front, candidates)]
    cnf.append([pool.id(('x', _)) for _ in possible])

    return cnf, pool

def breaking_mis(num: int, **kwds) -> Iterable[Tuple[int, ...]]:
    """
    Use some symmetry breaking.
    """
    cnf, pool = breaking_model(num)
    return solve_maxsat(cnf, pool, stem = 'x', **kwds)
