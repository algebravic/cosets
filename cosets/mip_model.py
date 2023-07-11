"""
Use the mip package to model independent set
"""
from typing import Iterable, Hashable, Tuple, Dict, List
from itertools import chain
import networkx as nx
from pysat.formula import CNF, WCNF, IDPool
from mip import Model, INTEGER, CONTINUOUS, BINARY, xsum, MAXIMIZE
from mip.entities import LinExpr

def mip_model(gph: nx.Graph) -> Tuple[Model, Dict[Hashable, int]]:
    """
    Use the standard mip_model.
    """
    ngph = nx.convert_node_labels_to_integers(gph,
                                              ordering='sorted')
    dct = dict(enumerate(sorted(gph.nodes)))
    mvars = {}
    model = Model(sense=MAXIMIZE)
    mvars = [model.add_var(name = f'x{node}', var_type=BINARY)
             for node in ngph.nodes]

    for node1, node2 in ngph.edges:
        model += mvars[node1] + mvars[node2] <= 1
    model.objective = xsum(mvars)

    return model, dct

def _get_clause(clause: List[int], variables: List[int]) -> LinExpr:
    """
    get clause inequality
    """
    rhs = 1 - len([_ for _ in clause if _ < 0])
    lhs = xsum([(2 * int(_ > 0) - 1) * variables[abs(_) - 1]
                for _ in clause])
    return lhs >= rhs

def cnf_model(cnf: WCNF) -> Tuple[Model, Dict[int, int]]:
    """
    Render a weight CNF as a MIP
    """
    model = Model(sense=MAXIMIZE)
    top = cnf.nv
    variables = [model.add_var(name = f'x{ind}', var_type=BINARY)
                 for ind in range(1, top+1)]

    # First the hard clauses
    for clause in cnf.hard:
        model += _get_clause(clause, variables)

    # Now the soft clauses
    objective = []
    for ind, (clause, wgt) in enumerate(zip(cnf.soft, cnf.wght), start=1):
        svar = model.add_var(name='s{ind}', var_type=BINARY)
        variables.append(svar)
        try:
            model += _get_clause([-(top+ind)] + clause, variables)
        except IndexError:
            print(f"top_ind = {top+ind}, len(variables) = {len(variables)}")
            raise
        objective.append(wgt * svar)

    model.objective = xsum(objective)
    return model
