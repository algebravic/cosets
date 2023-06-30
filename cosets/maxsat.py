"""
Use Max Sat for independent set.
"""
from typing import List, Tuple, Iterable, Any
from itertools import product
import networkx as nx
import numpy as np
from pysat.formula import WCNF, IDPool
from pysat.examples.rc2 import RC2
from .connection import dn_neighbors, candidates

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
    for delta in dn_neighbors(num):
        for elt in product(range(2), repeat = num + 1):
            nelt = np.array(elt, dtype=np.int8)
            other = tuple((nelt ^ delta).tolist())
            if other > elt:
                cnf.append([-pool.id(('x', elt)),
                            -pool.id(('x', other))])
    for elt in product(range(2), repeat = num + 1):
        cnf.append([pool.id(('x', elt))], weight=1)

    cnf.append([pool.id(('x', (num + 1) * (0,)))])
    # One of the nodes from the remaining classes must be there
    cnf.append([pool.id(('x', _)) for _ in candidates(num)])

    return cnf, pool

def breaking_mis(num: int, **kwds) -> Iterable[Tuple[int, ...]]:
    """
    Use some symmetry breaking.
    """
    cnf, pool = breaking_model(num)
    return solve_maxsat(cnf, pool, stem = 'x', **kwds)
