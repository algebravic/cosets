"""
Use the greedy method (and refinements) from
'On Reducing Maximum Independent set to minimum satisfiability'
by Ignatiev, Morgado and Marques-Silva
"""
from typing import List, Tuple, Iterable, FrozenSet, Hashable
from itertools import combinations, chain, product
from collections import defaultdict
import networkx as nx
from pysat.formula import CNF, WCNF, IDPool
from .maxsat import solve_maxsat
from .graphs import heuristic_partition

CLAUSE = List[int]

def greedy(gph: nx.Graph) -> Tuple[WCNF, IDPool]:
    """
    Simple greedy algorithm.
    Produce a Minsat instance.
    Convert it to MaxSat separately.
    """
    cnf = WCNF()
    pool = IDPool()
    ngph = gph.copy()
    formula = {node: [] for node in gph.nodes}

    while len(ngph.edges) > 0:
        maxdeg, maxnode = max(((ngph.degree(_), _) for _ in ngph.nodes))
        formula[maxnode].append(pool.id(('x', maxnode)))
        nbrs = set(list(ngph.neighbors(maxnode)))
        for nbr in nbrs:
            ngph.remove_edge(maxnode,nbr)
            formula[nbr].append(-pool.id(('x', maxnode)))
        ngph.remove_node(maxnode)
    for clause in formula.values():
        cnf.append(clause, weight=1)
    return cnf, pool

def ne_encoding(cnf: WCNF) -> WCNF:
    """
    Convert according to natural encoding.
    See: 'A New Encoding from MinSAT into MaxSAT'
    by Zhu, Li, Manya and Argelich
    """
    ncnf = WCNF()
    if cnf.hard:
        ncnf.extend(ncnf.hard)
    for cls, wgt in zip(cnf.soft, cnf.wght):
        front = []
        for lit in cls:
            ncnf.append(front + [-lit], weight=wgt)
            front.append(lit)
    return ncnf

                

def maxflow_conversion(soft: Tuple[CLAUSE, int]):
    """
    Use the maxflow rendering of MinSat to MaxSat.
    First construct a directed weight graph.
    
    """
    pass

def aux_graph(clauses: List[CLAUSE]) -> nx.Graph:
    """
    Construct the auxilliary graph.
    See 'Exact MinSAT Solving'

    The nodes of the graph are the clauses.
    Their is an edge between two nodes if and only if
    the elementwise complement of one has a non empty
    intersection with the other.
    """
    gph = nx.Graph()
    for node1, node2 in combinations(clauses, 2):
        complement = set([-_ for _ in node1])
        if complement.intersection(node2):
            gph.add_edge(tuple(node1), tuple(node2))
    return gph

def compatible_sets(gph: nx.Graph) -> List[Tuple[Hashable, ...]]:
    """
    See 'On Reducing MIS to MinSat
    """
    cnf, pool = greedy(gph)
    incompatible = defaultdict(lambda : set())
    igph = nx.Graph()
    for clause in cnf.soft:
        for lit1, lit2 in combinations(clause, 2):
            if lit1 * lit2 < 0:
                igph.add_edge(lit1, lit2)
    cgph = nx.complement(aux_graph(cnf.soft))
    for clause1, clause2 in cgph.edges:
        for lit1, lit2 in product(clause1, clause2):
            igph.add_edge(lit1, -lit2)
            igph.add_edge(lit2, -lit1)
    agph = nx.complement(igph)
    # positive literals
    positive = set(map(abs,chain(*cnf.soft)))
    classes = set()
    for node in positive:
        classes.add(frozenset(list(agph.neighbors(_)) + [_]))
    return classes

def clique_encoding(gph: nx.Graph) -> Tuple[CNF, IDPool]:
    """
    Clique partition encoding form 'Exact MinSat Solving'
    """
    parts = heuristic_partition(gph)
    cnf = WCNF()
    pool = IDPool()
    for node1, node2 in gph.edges:
        cnf.append([-pool.id(('c', node1)),
                    -pool.id(('c', node2))])
    for part in parts:
        cnf.append([pool.id(('c', _)) for _ in part], weight=1)
    return cnf, pool

def new_solve(gph: nx.Graph, **kwds) -> Tuple[CNF, IDPool]:
    """
    Independent set via MinSat
    """
    cnf, pool = greedy(gph)
    xcnf, xpool = clique_encoding(aux_graph(cnf.soft))
    asoln = solve_maxsat(xcnf, xpool, stem='c', **kwds)
    # Find the soft clauses
    psoln = set(chain(*asoln))
    # psoln is the set of positive literal in the graph cnf encoding
    npos = [pool.obj(_) for _ in psoln if _ > 0]
    return [_[1] for _ in npos if _ is not None and _[0] == 'x']
