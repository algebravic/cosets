"""
Use Max Sat for independent set.
"""
from typing import List, Tuple, Iterable, Any, Callable
from itertools import product
from functools import partial
from time import time
import networkx as nx
import numpy as np
from lazytree import LazyTree
from sympy.combinatorics import PermutationGroup
from pysat.formula import WCNF, IDPool
from pysat.examples.rc2 import RC2
from .schreier import make_tree, tree_clauses

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
    if soln is None:
        print("Formula is UNSAT!")
        return None
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

def trace_iterable(count: int, data: Iterable[Any]) -> Iterable[Any]:
    """
    Progress indication for a stream.
    """
    ind = 1
    for elt in data:
        yield elt
        ind += 1
        if count > 0 and ind % count == 0:
            print('X', end='', flush=True)
    print('')
    
def maxsat_mis_tree(gph: nx.Graph,
                    grp: PermutationGroup,
                    depth: int,
                    test: Callable[[LazyTree], bool] = lambda _: True,
                    trace: int = 0,
                    **kwds) -> Iterable[Any]:
    """
    Use the tree model for symmetry breaking.
    """
    start = time()
    cnf, pool = maxsat_mis_model(gph)
    cnf.append([pool.id(('x', min(gph.nodes)))])
    trace_it = partial(trace_iterable, trace)
    cnf.extend(trace_it(tree_clauses(pool, test, depth,
                                     make_tree(gph, grp))))
    end = time()
    print(f"model time = {end - start}")
    return solve_maxsat(cnf, pool, stem = 'x', **kwds)
